"""
Step 04: Remove hypermutator ALE lineages and apply frequency filtering.

Uses the mutil helpers.

Filter sequence (order matters):
  1. drop rows with empty Details (mostly multi-gene deletions), recorded in
     removed/04_removed_empty_details_rows.tsv
  2. remove hypermutator ALEs. An ALE is removed only if it is BOTH a mutation-burden
     outlier (a matplotlib boxplot flier, more than 1.5 * IQR above Q3 or below Q1) AND
     carries a coding mutation in a known hypermutator gene (listed in
     config_tables.hypermutator_genes). Both use the calls whose presence contains
     hypermutator_frequency_token, in the ALE's last flask with such calls
  3. remove high-burden flasks (>= high_mutation_count_threshold RAW mutation rows, counted
     over all samples of the flask), keyed on exp_ale_flask, so every flask of an ALE is
     checked, not only its endpoint. This runs BEFORE the frequency filter on purpose: it
     is the only check that sees the sub-fixation load, which is where hypermutation
     shows up
  4. drop non-standard gene annotations
  5. frequency filter (keep a row when one of its presence values is at/above
     population_frequency_threshold)
  6. extract ALE-unique endpoint mutations and build sample_endpoint = exp_ale_flask
     (rows not kept are recorded in removed/04_removed_not_endpoint_rows.tsv)

Two REPORT-ONLY diagnostics are written afterwards (no rows removed):
  - 04_endpoint_burden_report.csv        : per-endpoint burden after extraction, flagging
                                           Tukey outliers (> Q3 + 1.5 * IQR) and endpoints
                                           at/above high_mutation_count_threshold
  - 04_starting_strain_candidates_*.csv  : mutations near-fixed within an experiment
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.cbook import boxplot_stats

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "04"
STEP_NAME = "remove_hypermutators"
GENE_COL = "Gene (Scrollable)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    common.add_mutil_to_path(cfg)
    from mutil.mut import get_clean_mut_gene_list, is_coding_mut
    from mutil.alemutdf import get_multi_exp_max_freq_mut_df

    s4 = cfg["step04"]
    in_path = common.out_dir() / "03_starting_strain_removed.tsv"
    out_path = common.out_dir() / "04_hypermutators_removed.tsv"
    hyp_path = common.resolve(cfg["config_tables"]["hypermutator_genes"])

    HYPERMUTATOR_GENES = set()
    for line in Path(hyp_path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            HYPERMUTATOR_GENES.add(line)
    if not HYPERMUTATOR_GENES:
        raise ValueError(f"{hyp_path} lists no gene (config_tables.hypermutator_genes)")
    logger.info("loaded %d hypermutator genes", len(HYPERMUTATOR_GENES))

    muts_df = common.read_table(in_path)
    muts_df["exp ALE"] = muts_df.apply(lambda m: str(m["exp"]) + " " + str(m["ale"]), axis=1)
    n_before = len(muts_df)
    logger.info("loaded %d rows from %s", n_before, in_path)

    # --- drop NaN Details ---
    # Rows with an empty Details field (mostly multi-gene deletions) are not used.
    # They are written to removed/04_removed_empty_details_rows.tsv below.
    n0 = len(muts_df)
    nan_details = muts_df[muts_df["Details"].isna()]
    muts_df = muts_df[muts_df["Details"].notna()].copy()
    logger.info("removed %d rows with NaN Details", n0 - len(muts_df))

    # --- coding flag and noncoding gene-name harmonization ---
    muts_df["coding"] = muts_df.Details.apply(is_coding_mut)
    muts_df[GENE_COL] = muts_df.apply(
        lambda row: row[GENE_COL].replace(", ", "/") if not row["coding"] else row[GENE_COL], axis=1)
    # Genes as they read here, before the removals below, for the drop_gene_names check.
    genes_in_step = set(muts_df[GENE_COL])

    # --- clonal set, used for hypermutator detection only ---
    # The calls whose presence contains hypermutator_frequency_token.
    token = str(s4["hypermutator_frequency_token"])
    clone_pop_norm = muts_df[muts_df["presence"].str.contains(token, regex=False)]
    if clone_pop_norm.empty:
        raise ValueError(f"step04.hypermutator_frequency_token {token!r} is in no presence value")
    ep_max_freq = get_multi_exp_max_freq_mut_df(clone_pop_norm, endpoint_flask_only=True)

    def mutated_hyp_genes(m):
        if m["coding"]:
            genes = set(get_clean_mut_gene_list(m[GENE_COL]))
            return genes & HYPERMUTATOR_GENES
        return set()
    ep_max_freq["hypermutator gene hit"] = ep_max_freq.apply(
        lambda m: ", ".join(mutated_hyp_genes(m)), axis=1)

    # --- per-ALE burden and outlier detection ---
    data = []
    for exp, exp_mut_df in ep_max_freq.groupby("exp"):
        for exp_ale, ale_df in exp_mut_df.groupby("exp ALE"):
            data.append([exp, exp_ale, len(ale_df)])
    exp_ale_mut_cnt_df = pd.DataFrame(data, columns=["experiment", "exp ALE", "mutation count"])

    stats = boxplot_stats(exp_ale_mut_cnt_df["mutation count"]).pop(0)
    outlier_mut_cnts = stats["fliers"]
    whisker_hi = stats["whishi"]
    outlier_mut_count_ALEs = set(
        exp_ale_mut_cnt_df[exp_ale_mut_cnt_df["mutation count"].isin(outlier_mut_cnts)]["exp ALE"])
    ALEs_with_hyp_mut_gene = set(
        ep_max_freq[ep_max_freq["hypermutator gene hit"] != ""]["exp ALE"].unique())
    pred_hyp_mut_ALEs = ALEs_with_hyp_mut_gene & outlier_mut_count_ALEs
    # The upper whisker is the largest count at or below Q3 + 1.5 * IQR, not the cut itself.
    logger.info("total ALEs=%d upper whisker=%.0f hypermutator ALEs removed=%d",
                len(exp_ale_mut_cnt_df), whisker_hi, len(pred_hyp_mut_ALEs))

    # --- remove hypermutator ALEs ---
    m = muts_df["exp ALE"].isin(pred_hyp_mut_ALEs)
    removed_hyp = muts_df[m].copy()
    removed_hyp["removal_reason"] = "hypermutator ALE (burden outlier and hypermutator-gene hit)"
    n_removed_hyp_ales = removed_hyp["exp ALE"].nunique()
    muts_df = muts_df[~m]
    logger.info("removed %d hypermutator rows from %d ALEs", len(removed_hyp), n_removed_hyp_ales)

    # --- remove high-burden flasks ---
    # Keyed on exp_ale_flask, built like sample_endpoint below but for every flask of an ALE,
    # not only its endpoint. The sample column (ale flask isolate tech_rep) holds no exp, so it
    # is not used as the key. A flask's count is its number of rows, over all its samples.
    # This filter runs BEFORE the frequency filter on purpose: it is the only check that
    # sees the raw (sub-fixation) mutation load, which is where hypermutation shows up.
    thr = int(s4["high_mutation_count_threshold"])
    burden_key = (muts_df["exp"].astype(str) + "_"
                  + muts_df["ale"].astype(str) + "_" + muts_df["flask"].astype(str))
    flask_counts = burden_key.value_counts()
    high = flask_counts[flask_counts >= thr].index
    hb_mask = burden_key.isin(high)
    removed_hb = muts_df[hb_mask].copy()
    removed_hb["removal_reason"] = f"flask with >= {thr} mutations (raw, pre-frequency-filter)"
    muts_df = muts_df[~hb_mask].copy()
    logger.info("removed %d rows from %d high-burden flasks (>= %d raw mutations)",
                len(removed_hb), len(high), thr)

    # --- drop strange gene annotations ---
    drop_genes = s4["drop_gene_names"]
    absent = [g for g in drop_genes if g not in genes_in_step]
    if absent:
        raise ValueError(f"step04.drop_gene_names: no row has the gene {absent}")
    removed_genes = muts_df[muts_df[GENE_COL].isin(drop_genes)].copy()
    if len(removed_genes):
        removed_genes["removal_reason"] = "non-standard gene annotation " + str(drop_genes)
    muts_df = muts_df[~muts_df[GENE_COL].isin(drop_genes)]

    # --- low-frequency filter ---
    thr_f = float(s4["population_frequency_threshold"])

    def presence_filter(val):
        numbers = [float(x) for x in str(val).split("/")]
        return any(num >= thr_f for num in numbers)
    keep_freq = muts_df["presence"].apply(presence_filter)
    removed_freq = muts_df[~keep_freq].copy()
    if len(removed_freq):
        removed_freq["removal_reason"] = f"presence below frequency cutoff (>= {thr_f})"
    muts_df = muts_df[keep_freq]
    logger.info("frequency filter removed %d rows", len(removed_freq))

    # --- ALE-unique endpoint mutations ---
    if muts_df.empty:
        raise ValueError("no row is left for the endpoint pick after the removals above. Check "
                         "step04.high_mutation_count_threshold and population_frequency_threshold")
    pre_pick = muts_df.assign(_row=range(len(muts_df)))
    muts_df = get_multi_exp_max_freq_mut_df(
        pre_pick, endpoint_flask_only=True, representative_isolate_only=False)
    # Rows not kept as an ALE's endpoint row (earlier flasks, and further samples of the same
    # mutation in the endpoint flask). Recorded so that every input row is accounted for.
    not_endpoint = pre_pick[~pre_pick["_row"].isin(muts_df["_row"])].drop(columns="_row")
    muts_df = muts_df.drop(columns="_row")
    muts_df["sample_endpoint"] = (muts_df["exp"].astype(str) + "_"
                                  + muts_df["ale"].astype(str) + "_" + muts_df["flask"].astype(str))

    # --- post-endpoint burden diagnostic (REPORT ONLY, no rows are removed here) ---
    # The high-burden filter above acts on the RAW mutation load, before the frequency
    # filter. This report re-checks the burden on the analysis unit (sample_endpoint) after
    # extraction, so residual high-burden endpoints can be reviewed. Two flags are emitted:
    #   - exceeds_high_mutation_count_threshold : same absolute cutoff as the filter above
    #   - is_burden_outlier                     : Tukey rule, burden > Q3 + 1.5 * IQR
    ep_burden = muts_df.groupby("sample_endpoint").size().rename("n_mutations")
    ep_q1, ep_q3 = float(ep_burden.quantile(0.25)), float(ep_burden.quantile(0.75))
    ep_outlier_thr = ep_q3 + 1.5 * (ep_q3 - ep_q1)
    burden_report = ep_burden.reset_index()
    burden_report["is_burden_outlier"] = burden_report["n_mutations"] > ep_outlier_thr
    burden_report["exceeds_high_mutation_count_threshold"] = burden_report["n_mutations"] >= thr
    meta_cols = [c for c in ["exp", "strain", "strain_details", "clonal_or_population"]
                 if c in muts_df.columns]
    if meta_cols:
        ep_meta = muts_df.groupby("sample_endpoint")[meta_cols].first().reset_index()
        burden_report = burden_report.merge(ep_meta, on="sample_endpoint", how="left")
    # Equal counts are ordered by name, so the row order is fixed.
    burden_report = burden_report.sort_values(["n_mutations", "sample_endpoint"],
                                              ascending=[False, True])
    common.write_table(burden_report, common.qc_dir() / "04_endpoint_burden_report.csv")

    n_ep_outliers = int(burden_report["is_burden_outlier"].sum())
    n_ep_over_thr = int(burden_report["exceeds_high_mutation_count_threshold"].sum())
    logger.info("endpoint burden after extraction: median=%.0f Q3=%.0f outlier_threshold=%.0f max=%d",
                float(ep_burden.median()), ep_q3, ep_outlier_thr, int(ep_burden.max()))
    logger.info("endpoint burden report (REPORT ONLY, nothing removed): %d outliers (> Q3+1.5*IQR), "
                "%d endpoints >= %d mutations. See qc/04_endpoint_burden_report.csv",
                n_ep_outliers, n_ep_over_thr, thr)

    # --- starting-strain candidate report: mutations near-fixed within an experiment ---
    # (the same exact mutation in >= threshold of an experiment's endpoint samples is a
    #  starting-strain signature. It is reported for review, not removed here)
    thr_prev = float(s4["high_prevalence_report_threshold"])
    min_eps = int(s4["high_prevalence_min_endpoints"])
    # Reference Seq is not part of the key, because some deposits leave it empty and the user
    # starting-strain table does not match on it either. dropna=False counts a row with an empty
    # field, which groupby would otherwise leave out without a word.
    mut_key = ["exp", "Position", "Mutation Type", "Sequence Change", "Details", GENE_COL]
    mut_key = [c for c in mut_key if c in muts_df.columns]
    eps_per_exp = muts_df.groupby("exp")["sample_endpoint"].nunique().rename("n_endpoint_samples_in_exp")
    occ = (muts_df.groupby(mut_key, dropna=False)["sample_endpoint"].nunique()
           .rename("n_samples_with_mutation").reset_index())
    occ = occ.merge(eps_per_exp, on="exp", how="left")
    occ["prevalence"] = occ["n_samples_with_mutation"] / occ["n_endpoint_samples_in_exp"]
    cand = occ[(occ["prevalence"] >= thr_prev) & (occ["n_endpoint_samples_in_exp"] >= min_eps)] \
        .sort_values(["prevalence", "n_endpoint_samples_in_exp"], ascending=False)
    common.write_table(cand, common.qc_dir() / "04_starting_strain_candidates_high_prevalence.csv")
    logger.info("starting-strain candidate report: %d mutations at prevalence >= %.0f%% in experiments with >= %d endpoints (see qc/)",
                len(cand), 100 * thr_prev, min_eps)

    # --- removed-row tables (one file per removal reason) ---
    common.save_removed(nan_details, common.removed_dir() / "04_removed_empty_details_rows.tsv",
                        "empty Details")
    common.save_removed(removed_hyp, common.removed_dir() / "04_removed_hypermutator_rows.tsv",
                        "hypermutator ALE")
    common.save_removed(removed_hb, common.removed_dir() / "04_removed_high_burden_flasks.tsv",
                        f"flask with >= {thr} raw mutations")
    common.save_removed(removed_freq, common.removed_dir() / "04_removed_low_frequency_rows.tsv",
                        "presence below frequency cutoff")
    common.save_removed(removed_genes, common.removed_dir() / "04_removed_nonstandard_genes.tsv",
                        "non-standard gene annotation")
    common.save_removed(not_endpoint, common.removed_dir() / "04_removed_not_endpoint_rows.tsv",
                        "not the endpoint row (earlier flask, or another sample of the same mutation)")

    n_after = len(muts_df)
    common.check_conservation(n_before, n_after,
                              len(nan_details) + len(removed_hyp) + len(removed_hb) + len(removed_genes)
                              + len(removed_freq) + len(not_endpoint), logger)
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_file": str(in_path), "output_file": str(out_path),
        "n_rows_before": n_before, "n_rows_after": n_after,
        "n_removed_hypermutator_rows": int(len(removed_hyp)),
        "n_removed_hypermutator_ALEs": int(n_removed_hyp_ales),
        "n_removed_low_frequency_rows": int(len(removed_freq)),
        "n_removed_high_burden_rows": int(len(removed_hb)),
        "n_removed_nonstandard_gene_rows": int(len(removed_genes)),
        "n_removed_not_endpoint_rows": int(len(not_endpoint)),
        "n_removed_nan_details": int(len(nan_details)),
        "mutation_burden_q1": float(stats["q1"]),
        "mutation_burden_q3": float(stats["q3"]),
        "mutation_burden_iqr": float(stats["q3"] - stats["q1"]),
        "mutation_burden_upper_whisker": float(whisker_hi),
        "high_mutation_count_threshold": int(thr),
        "n_removed_high_burden_flasks": int(len(high)),
        "endpoint_burden_median": float(ep_burden.median()),
        "endpoint_burden_q3": float(ep_q3),
        "endpoint_burden_outlier_threshold": float(ep_outlier_thr),
        "endpoint_burden_max": int(ep_burden.max()),
        "n_endpoint_burden_outliers_reported": int(n_ep_outliers),
        "n_endpoints_over_high_mutation_threshold_reported": int(n_ep_over_thr),
        "hypermutator_genes_used": sorted(HYPERMUTATOR_GENES),
        "hypermutator_rows_by_experiment": removed_hyp["exp"].value_counts().to_dict(),
        "high_burden_key": "exp_ale_flask (every flask of an ALE, not only its endpoint)",
        "population_frequency_threshold": thr_f,
    }
    common.validate_common(muts_df, logger, required_cols=["exp", "Position", "sample_endpoint", "coding"])
    common.write_table(muts_df, out_path)
    common.write_checkpoint(metrics, common.checkpoints_dir() / "04_hypermutator_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
