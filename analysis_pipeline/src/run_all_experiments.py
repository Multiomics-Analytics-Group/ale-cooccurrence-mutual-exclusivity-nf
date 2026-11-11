"""
run_all_experiments.py

Runs the workflow independently for every experiment with at least the
configured minimum number of endpoint samples, then concatenates the results
into combined tables:
  pairwise_all_methods_<unit_tag>.csv
  per_experiment_summary_<unit_tag>.csv
  significant_pairs_<unit_tag>.csv
  combined_summary_<unit_tag>.csv
  skipped_experiments_<unit_tag>.csv (experiments below the minimum number of endpoints)

For each experiment: build & save the shared matrix -> DISCOVER + Fisher
significance + Pearson/phi (descriptive) on that matrix, plus a per-experiment
Pearson correlation matrix.

Every output of the run goes to output/analysis/<run_name>/. While a run is
going, the folder <run_name>.lock beside it keeps other runs with the same
run_name out. A finished run (one that wrote run_finished.txt) is never
overwritten. A folder left by an unfinished run of this pipeline (one that holds
run_started.txt) is removed and the run starts again. Any other folder that is
not empty is left alone and the run stops. An error in any experiment stops the
run, and the log names the experiment and where the error happened.

The run folder records how the run was made: run_started.txt gives the start
time, the command and the code version (git commit), and run_config.yaml is a
copy of the config the run read.
"""
import argparse
import importlib.util
import os
import shlex
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
REPO_ROOT = SRC.parent.parent


def _load(mod_name, file_name):
    spec = importlib.util.spec_from_file_location(mod_name, SRC / file_name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


cfgmod = _load("analysis_config_mod", "config.py")
iomod = _load("analysis_io_mod", "io.py")
mb = _load("analysis_matrix_mod", "matrix_builder.py")
pt = _load("analysis_pairwise_mod", "pairwise_tests.py")
sm = _load("analysis_summary_mod", "summarize_results.py")

# Written first and last, so a run folder shows whether this pipeline made it and whether the
# run finished.
STARTED = "run_started.txt"
FINISHED = "run_finished.txt"


def run_one(df, cfg, exp_name, logger):
    """The full workflow for one experiment: tables A (pairs), B (per-experiment summary),
    C (significant pairs) and the correlation matrix."""
    df_scope = df[df[cfg.exp_col] == exp_name]
    n_input_rows = len(df_scope)
    if not n_input_rows:
        raise ValueError(f"experiment {exp_name!r} has no rows in {cfg.input_table}")

    mb.build_and_save(df_scope, cfg, exp_name, logger=logger)
    art = mb.load_matrix(cfg, exp_name)
    M, freq, meta = art["matrix"], art["frequencies"], art["meta"]
    tested = freq.loc[freq["passed_filter"], "unit"].tolist()
    corr = pt.correlation_matrix(M, tested)

    logger.info("validation: input_rows=%d samples=%d units_total=%d units_tested=%d pairs=%d density=%.4f "
                "burden(mean/median/max)=%.2f/%.1f/%d",
                n_input_rows, meta["n_samples"], meta["n_units_total"], meta["n_units_tested"],
                meta["n_pairs_tested"], meta["density"], meta["mean_mutation_burden"],
                meta["median_mutation_burden"], meta["max_mutation_burden"])

    if len(tested) < 2:
        note = "too small" if meta["n_samples"] < 2 or meta["n_units_total"] < 2 else "<2 testable units"
        summary = sm.build_experiment_summary(pd.DataFrame(), meta, cfg, exp_name, n_input_rows, note=note)
        return {"pairs": pd.DataFrame(), "summary": summary, "significant": pd.DataFrame(),
                "correlation_matrix": corr, "meta": meta}

    contingency = pt.pair_contingency(M, tested)
    dco = pt.run_discover(M, tested, significance_threshold=cfg.significance_threshold, fdr_method=cfg.discover_fdr_method)
    fis = pt.run_fisher(M, tested, significance_threshold=cfg.significance_threshold)
    pea = pt.run_pearson(M, tested)

    pairs = sm.build_pair_table(contingency, dco, fis, pea, meta, freq, cfg, exp_name)
    summary = sm.build_experiment_summary(pairs, meta, cfg, exp_name, n_input_rows, note="ok")
    sig = sm.build_significant_table(pairs)
    logger.info("DISCOVER sig: co=%d mutex=%d | Fisher: co=%d mutex=%d | Pearson r reported (descriptive)",
                int(pairs["discover_sig_co"].sum()), int(pairs["discover_sig_mutex"].sum()),
                int(pairs["fisher_sig_co"].sum()), int(pairs["fisher_sig_mutex"].sum()))
    return {"pairs": pairs, "summary": summary, "significant": sig, "correlation_matrix": corr,
            "meta": meta}


def _utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _code_version():
    """The git commit of the code, and whether analysis_pipeline differs from it. Without git, or
    outside a git checkout, it is recorded as unknown. A checkout that git cannot read stops the
    run."""
    if shutil.which("git") is None:
        return "unknown (git is not installed)"
    if not (REPO_ROOT / ".git").exists():
        return f"unknown ({REPO_ROOT} is not a git checkout)"
    head = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
                          capture_output=True, text=True)
    if head.returncode:
        raise RuntimeError(f"git could not read the commit of {REPO_ROOT}: {head.stderr.strip()}")
    changes = subprocess.run(["git", "-C", str(REPO_ROOT), "status", "--porcelain", "--",
                              "analysis_pipeline"], capture_output=True, text=True)
    if changes.returncode:
        raise RuntimeError(f"git could not list the changes in {REPO_ROOT / 'analysis_pipeline'}: "
                           f"{changes.stderr.strip()}")
    commit = head.stdout.strip()
    return f"{commit} with uncommitted changes in analysis_pipeline" if changes.stdout.strip() else commit


def _running(pid):
    """True when a process with this id exists on this machine. Windows cannot be asked without
    side effects, so there a locked run is taken as still going."""
    if os.name == "nt":
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:     # it exists, under another user
        return True
    return True


def _take_lock(lock_dir, run_dir):
    """Take the lock: a folder beside the run folder, made with mkdir, which only one run can do.
    A lock left by a run that no longer exists is replaced."""
    try:
        lock_dir.mkdir()
    except FileExistsError:
        owner = lock_dir / "owner.txt"
        parts = owner.read_text().split() if owner.exists() else []
        pid = int(parts[0]) if parts and parts[0].isdigit() else None
        host = parts[1] if len(parts) > 1 else None
        if pid is None or host is None:
            raise RuntimeError(f"{run_dir} is locked by a run that cannot be identified ({lock_dir}). "
                               f"If no run is going, delete {lock_dir} and start again") from None
        if host != socket.gethostname():
            raise RuntimeError(f"{run_dir} is in use by a run on {host}. If no run is going, delete "
                               f"{lock_dir} and start again") from None
        if _running(pid):
            raise RuntimeError(f"{run_dir} is in use by a run that is still going (process {pid}). If "
                               f"no run is going, delete {lock_dir} and start again") from None
        shutil.rmtree(lock_dir)
        try:
            lock_dir.mkdir()
        except FileExistsError:
            raise RuntimeError(f"{run_dir} was just locked by another run. Nothing was run") from None
    iomod.write_new(lock_dir / "owner.txt",
                    lambda tmp: tmp.write_text(f"{os.getpid()} {socket.gethostname()} {_utc_now()}\n"))


def start_run_folder(run_dir):
    """Take the lock and make the run folder ready. A finished run is never overwritten. A folder
    left by an unfinished run of this pipeline (it holds run_started.txt) is removed. Any other
    folder that is not empty is left alone and the run stops. Returns the lock folder and whether
    an unfinished run was removed."""
    lock_dir = run_dir.with_name(run_dir.name + ".lock")
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    _take_lock(lock_dir, run_dir)
    try:
        cleared = False
        if run_dir.exists():
            finished = run_dir / FINISHED
            if finished.exists():
                lines = finished.read_text().splitlines()
                when = lines[0] if lines else f"{FINISHED} is empty"
                raise FileExistsError(f"{run_dir} holds a finished run ({when}). A finished run is never "
                                      "overwritten: use a new run_name, or delete the folder to run it again")
            if any(run_dir.iterdir()):
                if not (run_dir / STARTED).exists():
                    raise FileExistsError(f"{run_dir} exists but was not made by a run of this pipeline "
                                          f"(it has no {STARTED}). Move it away, or choose another run_name")
                cleared = True
            shutil.rmtree(run_dir)
        run_dir.mkdir()
    except BaseException:
        shutil.rmtree(lock_dir)
        raise
    return lock_dir, cleared


def run(cfg, logger):
    df = iomod.load_table(cfg, logger)

    ec = df.groupby(cfg.exp_col)[cfg.sample_col].nunique()
    exps = sorted(ec[ec >= cfg.min_endpoints].index.tolist())
    logger.info("experiments with >= %d endpoints: %d", cfg.min_endpoints, len(exps))
    below = ec[ec < cfg.min_endpoints].sort_index()
    skipped = pd.DataFrame({"exp": below.index, "n_samples": below.values,
                            "n_input_rows": df.groupby(cfg.exp_col).size().loc[below.index].values})
    if not exps:
        raise ValueError(f"no experiment in {cfg.input_table} has at least {cfg.min_endpoints} "
                         "endpoints (min_endpoints)")
    # Each experiment gets folders and files named after it. Two names that become the same folder
    # name (also when upper and lower case are taken as equal) would mix their files.
    folders = {}
    for exp in exps:
        key = mb._safe(exp).casefold()
        if key in folders:
            raise ValueError(f"experiments {folders[key]!r} and {exp!r} would share the folder name "
                             f"{mb._safe(exp)!r}. Rename one of them in the input table")
        folders[key] = exp
    iomod.console(f"Analyzing {len(exps)} experiments (>= {cfg.min_endpoints} endpoints)...")

    out_root = cfg.output_dir / "all_experiments"
    skipped_path = out_root / sm.out_name(cfg, "skipped_experiments")
    iomod.write_table(skipped, skipped_path)
    logger.info("skipped %d experiments with fewer than %d endpoints (%s)",
                len(skipped), cfg.min_endpoints, skipped_path)
    all_pairs, all_summary, all_sig = [], [], []
    corr_dir = out_root / "correlation_matrices"
    sig_dir_root = cfg.output_dir / "significant_experiments"
    n_sig_exps = 0
    bar = iomod.progress(exps, desc="experiments", unit="exp")
    for i, exp in enumerate(bar, 1):
        try:
            res = run_one(df, cfg, exp, logger)
            # Descriptive Pearson/phi correlation matrix for this experiment. Every analyzed
            # experiment gets one, empty when no unit is tested.
            iomod.write_matrix(res["correlation_matrix"], corr_dir / f"{mb._safe(exp)}_{cfg.unit_tag}.csv")
            all_summary.append(res["summary"])
            if len(res["pairs"]):
                all_pairs.append(res["pairs"])
                # Experiments with at least one significant pair get a dedicated
                # folder with a table of those pairs and all their metrics (Pearson,
                # Fisher, DISCOVER column order) and the experiment summary.
                # significant_in says which test and direction.
                sig_rows = res["pairs"][res["pairs"]["significant_in"] != "none"]
                if len(sig_rows):
                    n_sig_exps += 1
                    exp_dir = sig_dir_root / mb._safe(exp)
                    iomod.write_table(
                        sig_rows.sort_values(["discover_sig_any", "significant_in_both"],
                                             ascending=False),
                        exp_dir / sm.out_name(cfg, "significant_pairs_detail"))
                    iomod.write_table(res["summary"], exp_dir / sm.out_name(cfg, "experiment_summary"))
            if len(res["significant"]):
                all_sig.append(res["significant"])
        except Exception as e:
            logger.error("experiment '%s' failed, run stopped: %s: %s", exp, type(e).__name__, e)
            raise
        if hasattr(bar, "set_postfix"):
            bar.set_postfix(significant=n_sig_exps)
        if i % 20 == 0:
            logger.info("progress: %d / %d experiments", i, len(exps))

    if not all_pairs:
        raise ValueError(f"no experiment has at least 2 tested units (fixed_min_count = "
                         f"{cfg.fixed_min_count}), so no pair was tested")
    pairs = pd.concat(all_pairs, ignore_index=True)
    summary = pd.concat(all_summary, ignore_index=True)
    # The significant pairs are a subset of the pair rows, so an empty table keeps their columns.
    sig = pd.concat(all_sig, ignore_index=True) if all_sig else pairs.iloc[0:0]

    iomod.write_table(pairs, out_root / sm.out_name(cfg, "pairwise_all_methods"))
    iomod.write_table(summary, out_root / sm.out_name(cfg, "per_experiment_summary"))
    iomod.write_table(sig, out_root / sm.out_name(cfg, "significant_pairs"))

    # Overview: one row per experiment, with the Fisher and DISCOVER significant
    # pair counts.
    combined = sm.build_combined_summary(all_summary)
    iomod.write_table(combined, out_root / sm.out_name(cfg, "combined_summary"))

    n_disc = int(summary["discover_mutex_sig_local"].fillna(0).sum() + summary["discover_co_sig_local"].fillna(0).sum())
    logger.info("done. experiments analyzed=%d, experiments with significant pairs=%d, total DISCOVER significant pairs=%d",
                len(exps), n_sig_exps, n_disc)
    logger.info("combined_summary (top by signal):\n%s",
                combined.head(10)[["exp", "endpoints", "n_units_tested",
                                   "fisher_co_sig", "fisher_mutex_sig",
                                   "discover_co_sig", "discover_mutex_sig"]].to_string(index=False))
    logger.info("wrote combined outputs to %s, per-experiment significant folders under %s",
                out_root, sig_dir_root)

    iomod.console("")
    iomod.console(f"Done. {len(exps)} experiments analyzed | "
                  f"{n_sig_exps} with significant pairs | {n_disc} DISCOVER significant pairs.")
    iomod.console(f"Outputs: {out_root}")


def main():
    ap = argparse.ArgumentParser(allow_abbrev=False)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = cfgmod.load_config(args.config)
    command = f"cd {shlex.quote(os.getcwd())} && {shlex.join([sys.executable] + sys.argv)}"
    code = _code_version()
    lock_dir, cleared = start_run_folder(cfg.output_dir)
    try:
        iomod.write_new(cfg.output_dir / STARTED, lambda tmp: tmp.write_text(
            f"started {_utc_now()}\nconfig {cfg.path}\ncommand {command}\ncode {code}\n"))
        iomod.write_new(cfg.output_dir / "run_config.yaml", lambda tmp: tmp.write_bytes(cfg.raw))
        logger = iomod.setup_logging("run_all_experiments", cfg.output_dir / "logs")
        if cleared:
            msg = f"{cfg.output_dir} held an unfinished run, which was removed before this run."
            logger.info(msg)
            iomod.console(msg)
        run(cfg, logger)
        iomod.write_new(cfg.output_dir / FINISHED,
                        lambda tmp: tmp.write_text(f"finished {_utc_now()}\nconfig {cfg.path}\n"))
    finally:
        shutil.rmtree(lock_dir)


if __name__ == "__main__":
    main()
