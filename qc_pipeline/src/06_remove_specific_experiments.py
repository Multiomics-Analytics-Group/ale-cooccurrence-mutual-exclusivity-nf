"""
Step 06: Remove manually specified problematic experiments.

Experiment names and the reason for each removal live in
config/manual_experiment_exclusion_table.tsv (columns: exp, match_type, remove, reason,
notes). Nothing is hard-coded here. match_type is startswith, equals (also exact or ==)
or contains. With contains, exp is read as a regular expression. remove is yes or no,
and notes is not read.
Deposits that copy another one are listed in that table too, with the reason.
A rule whose experiments left the run in an earlier step (for example a deposit that step 01
skipped, or an experiment of another species) is only reported, as a warning. A rule that
matches no experiment that step 01 recorded stops the run, since it is most likely a typo.
Step 01 records a deposit whose mutation file has no rows among the skipped deposits, so its
rule is only reported. An experiment that ALEdb no longer holds is not recorded, so its rule
needs remove set to no.

Its result is the final table of the QC pipeline, qc_mutations.tsv, which the analysis
reads. The final checks run just before it is written, so a table that fails them is
never written.
"""
import argparse
import importlib.util
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "06"
STEP_NAME = "remove_specific_experiments"


def read_step01_names():
    """Experiment names that step 01 recorded: its table, and its records of skipped deposits and
    of rows removed while loading. Only the exp column is read, with pd.read_csv's default
    low_memory parsing (common.read_table turns it off), so the large step 01 table needs
    little memory."""
    paths = [common.out_dir() / "01_aledb_mutation_df.tsv"] + [
        common.removed_dir() / name for name in ("01_skipped_duplicate_deposits.tsv",
                                                 "01_removed_rows_without_gene.tsv",
                                                 "01_removed_duplication_rows.tsv")]
    return pd.concat([pd.read_csv(p, sep="\t", usecols=["exp"], dtype=str)["exp"]
                      for p in paths]).drop_duplicates()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    in_path = common.out_dir() / "05_gene_names_normalized.tsv"
    # The final table of the QC pipeline, read by the analysis. A step added after this one
    # takes over the name.
    out_path = common.out_dir() / "qc_mutations.tsv"
    excl_path = common.resolve(cfg["config_tables"]["manual_experiment_exclusion"])

    df = common.read_table(in_path)
    n_before = len(df)
    n_exp_before = df["exp"].nunique()
    logger.info("loaded %d rows from %s", n_before, in_path)

    excl = pd.read_csv(excl_path, sep="\t", dtype=str, keep_default_na=False)
    removed_parts = []
    removal_reasons = {}
    removed_rows_by_experiment = {}
    removed_earlier = []
    step01_names = None

    # A rule that matches no experiment of this step's input is checked against the names step 01
    # recorded, read only then. If it matches one, its experiments left the run in an earlier
    # step, which is reported. If it matches none, the run stops. The check reads this step's
    # input, not the table left by the rules before it, so it does not depend on the rule order.
    step_input = df
    for _, r in excl.iterrows():
        remove = str(r["remove"]).strip().lower()
        if remove not in ("yes", "true", "1", "no", "false", "0"):
            raise ValueError(f"{excl_path}: remove must be yes or no, found {r['remove']!r} "
                             f"for {r['exp']!r}")
        if remove in ("no", "false", "0"):
            continue
        pat, mtype = r["exp"], r["match_type"].strip().lower()
        if not pat.strip():
            raise ValueError(f"{excl_path}: a rule has an empty exp, which would match every row")
        if mtype == "startswith":
            match = lambda exp: exp.str.startswith(pat)
        elif mtype in ("equals", "exact", "=="):
            match = lambda exp: exp == pat
        elif mtype == "contains":
            match = lambda exp: exp.str.contains(pat, na=False)
        else:
            raise ValueError(f"{excl_path}: unknown match_type {mtype!r} for {pat!r}. Use "
                             "startswith, equals, exact, == or contains")
        if not match(step_input["exp"]).any():
            if step01_names is None:
                step01_names = read_step01_names()
            if not match(step01_names).any():
                raise ValueError(f"{excl_path}: the rule [{mtype} {pat}] matches no experiment that "
                                 "step 01 loaded or recorded in removed/01_*.tsv. Check the spelling. "
                                 "If the experiment is no longer in ALEdb, set remove to no.")
            logger.warning("rule [%s %s]: nothing to remove, its experiments left the run before "
                           "step 06 (%s)", mtype, pat, r["reason"])
            removed_earlier.append(f"{mtype}:{pat}")
            continue
        m = match(df["exp"])
        n = int(m.sum())
        logger.info("rule [%s %s]: removes %d rows (%s)", mtype, pat, n, r["reason"])
        if n:
            part = df[m].copy()
            part["removal_reason"] = f"{r['reason']} [{mtype}:{pat}]"
            removed_parts.append(part)
            removal_reasons[pat] = r["reason"]
            for e, c in part["exp"].value_counts().items():
                removed_rows_by_experiment[e] = int(c)
            df = df[~m]

    # Final checks of the table the analysis reads. Each repeats a guarantee of an earlier step
    # (step 02 keeps only the configured species, step 05 removes rows left without a gene name),
    # so they only stop the run if a later change breaks one. Condition_group is checked with the
    # required columns below.
    species = cfg["organism"]["species"]
    problems = []
    other = df[df["species_norm"] != species]
    if len(other):
        found = sorted(other["species_norm"].astype(str).unique())
        problems.append(f"{len(other)} rows are not {species}: {found}")
    no_gene = int(df["Gene_name_normalized"].isna().sum())
    if no_gene:
        problems.append(f"{no_gene} rows have no gene name (Gene_name_normalized is empty)")
    if problems:
        raise ValueError("final checks failed: " + " | ".join(problems))

    removed_all = pd.concat(removed_parts, ignore_index=True) if removed_parts else df.iloc[0:0].copy()
    common.save_removed(removed_all, common.removed_dir() / "06_removed_specific_experiments.tsv",
                        "manually excluded experiment")

    n_after = len(df)
    common.check_conservation(n_before, n_after, len(removed_all), logger)
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_file": str(in_path), "output_file": str(out_path),
        "n_rows_before": n_before, "n_rows_after": n_after, "n_removed": len(removed_all),
        "n_experiments_before": int(n_exp_before),
        "n_experiments_after": int(df["exp"].nunique()),
        "removed_experiments": sorted(removed_all["exp"].unique().tolist()) if len(removed_all) else [],
        "removed_rows_by_experiment": removed_rows_by_experiment,
        "removal_reasons": removal_reasons,
        "rules_with_experiments_removed_earlier": removed_earlier,
    }
    common.validate_common(df, logger, required_cols=["exp", "Position", "Gene_name_normalized",
                                                      "Condition_group"])
    logger.info("final checks passed: every row is %s and has a gene name, and the Condition_group "
                "column is present", species)
    common.write_table(df, out_path)
    common.write_checkpoint(metrics, common.checkpoints_dir() / "06_remove_specific_exps_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
