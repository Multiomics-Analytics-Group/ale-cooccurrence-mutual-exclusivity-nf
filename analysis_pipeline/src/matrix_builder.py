"""
matrix_builder.py  (UNIT-AGNOSTIC, used by run_all_experiments.py)

Builds the generic binary matrix and saves it, together with unit-frequency and
mutation-burden summaries, so the three analysis methods all load the exact same
matrix instead of rebuilding it.

Matrix orientation (fixed and documented everywhere):
    rows    = test units   (whatever unit_col points at)
    columns = samples       (sample_col)
    entry   = 1 if the unit is mutated at least once in that sample, else 0

The unit is defined only by unit_col. There is no gene/operon switch.
"""
import json
import re
from itertools import combinations
from pathlib import Path

import importlib.util

import numpy as np
import pandas as pd


def _load_local(mod_name, file_name):
    """Load a sibling module by path so io.py does not collide with stdlib io."""
    spec = importlib.util.spec_from_file_location(mod_name, Path(__file__).with_name(file_name))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


iomod = _load_local("analysis_io_mod", "io.py")


def _safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")[:80] or "scope"


def scope_dir(cfg, scope_name):
    return cfg.matrices_dir / cfg.unit_tag / _safe(scope_name)


def build_matrix(df_scope, cfg, scope_name, logger=None):
    """Build the binary matrix and summaries for one analysis scope."""
    gl = df_scope.groupby([cfg.unit_col, cfg.sample_col]).size().reset_index(name="count")
    gl["mutated"] = 1
    gm = gl.pivot_table(index=cfg.unit_col, columns=cfg.sample_col,
                        values="mutated", fill_value=0, aggfunc="max")
    gm = (gm > 0).astype(int)
    if gm.shape[0]:
        gm = gm.loc[gm.sum(1) >= 1, gm.sum(0) > 0]

    n_units_total, n_samples = gm.shape
    counts = gm.sum(1) if n_units_total else pd.Series(dtype=int)
    burden = gm.sum(0) if n_samples else pd.Series(dtype=int)        # per-sample mutated-unit count

    # A unit is tested when it is mutated in at least fixed_min_count samples.
    count_min = int(cfg.fixed_min_count) if n_units_total else None
    passed = (counts >= count_min) if n_units_total else pd.Series(dtype=bool)
    n_tested = int(passed.sum()) if n_units_total else 0
    n_pairs = int(n_tested * (n_tested - 1) // 2)

    # unit members: composite units expand via members_col and atomic units map to themselves
    if cfg.members_col:
        mem = (df_scope.groupby(cfg.unit_col)[cfg.members_col]
               .apply(lambda x: ";".join(sorted({str(v) for v in x.dropna()}))).to_dict())
    else:
        mem = {}

    def members_of(u):
        return mem.get(u, str(u))

    # mean sample burden when a unit is mutated
    def mean_burden(u):
        if counts.get(u, 0) == 0:
            return float("nan")
        return float((gm.loc[u] * burden).sum() / counts[u])

    freq_rows = []
    for u in gm.index:
        c = int(counts[u])
        f = c / n_samples if n_samples else float("nan")
        freq_rows.append({
            "unit": u,
            "unit_members": members_of(u),
            "count_mutated_samples": c,
            "frequency": f,
            "passed_filter": bool(passed[u]),
            "n_possible_partners": (n_tested - 1) if passed[u] else 0,
            "mean_sample_burden_when_mutated": mean_burden(u),
        })
    freq_df = pd.DataFrame(freq_rows)
    burden_df = pd.DataFrame({"sample": burden.index, "mutation_burden": burden.values})

    meta = {
        "scope": scope_name,
        "unit_tag": cfg.unit_tag,
        "sample_col": cfg.sample_col,
        "unit_col": cfg.unit_col,
        "orientation": "rows=units, columns=samples",
        "n_samples": int(n_samples),
        "n_units_total": int(n_units_total),
        "n_units_tested": int(n_tested),
        "n_pairs_tested": int(n_pairs),
        "count_min": None if count_min is None else int(count_min),
        "significance_threshold": cfg.significance_threshold,
        "density": float(gm.values.mean()) if gm.size else float("nan"),
        "mean_mutation_burden": float(burden.mean()) if n_samples else float("nan"),
        "median_mutation_burden": float(burden.median()) if n_samples else float("nan"),
        "max_mutation_burden": int(burden.max()) if n_samples else 0,
    }
    if logger:
        logger.info("scope '%s': units_total=%d samples=%d count_min=%s units_tested=%d pairs=%d density=%.4f",
                    scope_name, n_units_total, n_samples, meta["count_min"], n_tested, n_pairs, meta["density"])
    return {"matrix": gm, "frequencies": freq_df, "burden": burden_df, "meta": meta}


def save_artifacts(artifacts, cfg, scope_name):
    d = scope_dir(cfg, scope_name)
    iomod.write_matrix(artifacts["matrix"], d / "matrix.csv")      # index = unit, columns = samples
    iomod.write_table(artifacts["frequencies"], d / "unit_frequencies.csv")
    iomod.write_table(artifacts["burden"], d / "sample_burden.csv")
    iomod.write_json(artifacts["meta"], d / "meta.json")
    return d


def load_matrix(cfg, scope_name):
    d = scope_dir(cfg, scope_name)
    gm = pd.read_csv(d / "matrix.csv", index_col=0)
    gm.columns = gm.columns.astype(str)
    freq = pd.read_csv(d / "unit_frequencies.csv")
    burden = pd.read_csv(d / "sample_burden.csv")
    with open(d / "meta.json") as f:
        meta = json.load(f)
    return {"matrix": gm, "frequencies": freq, "burden": burden, "meta": meta}


def build_and_save(df_scope, cfg, scope_name, logger=None):
    artifacts = build_matrix(df_scope, cfg, scope_name, logger=logger)
    out = save_artifacts(artifacts, cfg, scope_name)
    if logger:
        logger.info("saved matrix artifacts to %s", out)
    return artifacts


if __name__ == "__main__":
    raise SystemExit("matrix_builder.py has no command of its own. Put the experiments you want in an "
                     "input table and run analysis_pipeline/scripts/run_all.sh on it")
