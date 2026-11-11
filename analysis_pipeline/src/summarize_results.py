"""
summarize_results.py

Assembles the three output tables from the shared matrix and the three method
blocks:
  A. pairwise_all_methods_<unit_tag>.csv   (one row per tested pair)
  B. per_experiment_summary_<unit_tag>.csv (one row per experiment)
  C. significant_pairs_<unit_tag>.csv      (pairs significant in Fisher or DISCOVER)
and condenses B into combined_summary_<unit_tag>.csv, sorted by the DISCOVER and then
the Fisher significant pair counts, largest first.

DISCOVER is the primary, burden-aware test. Fisher is the exact test. Pearson/phi
is reported as a descriptive effect size (pearson_r) only, not a significance
test, so it never promotes a pair into the significant set.
"""
import numpy as np
import pandas as pd


def build_pair_table(contingency, discover_df, fisher_df, pearson_df, meta, freq_df, cfg,
                     scope_name):
    df = contingency.merge(discover_df, on=["unit1", "unit2"], how="left")
    df = df.merge(fisher_df, on=["unit1", "unit2"], how="left")
    df = df.merge(pearson_df, on=["unit1", "unit2"], how="left")

    # Composite units only: emit the member-gene expansion. For atomic (gene-level)
    # units the members equal the unit name, so these columns are omitted.
    if cfg.members_col:
        members = dict(zip(freq_df["unit"].astype(str), freq_df["unit_members"].astype(str)))
        df["unit1_members"] = df["unit1"].astype(str).map(members).fillna(df["unit1"].astype(str))
        df["unit2_members"] = df["unit2"].astype(str).map(members).fillna(df["unit2"].astype(str))

    analysis_id = f"{scope_name}__{cfg.unit_tag}"
    df.insert(0, "analysis_id", analysis_id)
    df.insert(1, "unit_tag", cfg.unit_tag)
    df.insert(2, "exp", scope_name)
    df["n_samples"] = meta["n_samples"]
    df["n_units_total"] = meta["n_units_total"]
    df["n_units_tested"] = meta["n_units_tested"]
    df["n_pairs_tested"] = meta["n_pairs_tested"]
    df["significance_threshold"] = meta["significance_threshold"]

    def b(x):
        return bool(x) if pd.notna(x) else False
    dco = df["discover_sig_co"].map(b); dme = df["discover_sig_mutex"].map(b)
    fco = df["fisher_sig_co"].map(b); fme = df["fisher_sig_mutex"].map(b)

    # Only Fisher and DISCOVER count as significance tests. Pearson/phi is a
    # descriptive effect size (pearson_r) and never contributes here.
    df["n_methods_support_co"] = fco.astype(int) + dco.astype(int)
    df["n_methods_support_mutex"] = fme.astype(int) + dme.astype(int)

    # Which test(s) is the pair significant in, and in which direction.
    def sig_in(i):
        parts = []
        if fco[i]: parts.append("Fisher(co)")
        if fme[i]: parts.append("Fisher(mutex)")
        if dco[i]: parts.append("DISCOVER(co)")
        if dme[i]: parts.append("DISCOVER(mutex)")
        return "; ".join(parts) if parts else "none"

    def direction(i):
        co_any = fco[i] or dco[i]
        me_any = fme[i] or dme[i]
        if co_any and me_any: return "mixed"
        if co_any: return "co"
        if me_any: return "mutex"
        return "none"

    df["significant_in"] = [sig_in(i) for i in df.index]
    df["significant_in_both"] = [bool((fco[i] or fme[i]) and (dco[i] or dme[i])) for i in df.index]
    df["direction"] = [direction(i) for i in df.index]
    return order_pair_table(df)


# Column order: metadata, pair identity, shared contingency, then the descriptive
# Pearson metric, then the two significance tests (Fisher, DISCOVER), then summary.
PAIR_COLUMN_ORDER = (
    ["analysis_id", "unit_tag", "exp", "n_samples", "n_units_total",
     "n_units_tested", "n_pairs_tested", "significance_threshold"]
    + ["unit1", "unit2", "unit1_members", "unit2_members"]
    + ["count_unit1", "count_unit2", "only_unit1", "only_unit2", "neither",
       "observed_overlap", "freq_unit1", "freq_unit2", "jaccard_index"]
    + ["pearson_r", "pearson_direction"]
    + ["fisher_expected_overlap", "fisher_observed_minus_expected", "fisher_odds_ratio",
       "fisher_log2_odds_ratio", "fisher_p_two_sided",
       "fisher_p_co", "fisher_q_co_bh_local", "fisher_sig_co",
       "fisher_p_mutex", "fisher_q_mutex_bh_local", "fisher_sig_mutex",
       "fisher_direction", "fisher_sig_any"]
    + ["discover_expected_overlap", "discover_observed_minus_expected", "discover_enrichment_ratio",
       "discover_z", "discover_p_co", "discover_q_co_local", "discover_sig_co",
       "discover_p_mutex", "discover_q_mutex_local", "discover_sig_mutex",
       "discover_direction", "discover_sig_any"]
    + ["significant_in", "significant_in_both", "n_methods_support_co", "n_methods_support_mutex",
       "direction"]
)


def order_pair_table(df):
    cols = [c for c in PAIR_COLUMN_ORDER if c in df.columns]
    extra = [c for c in df.columns if c not in cols]
    return df[cols + extra]


# Combined per-experiment summary, sorted by the DISCOVER and then the Fisher significant
# pair counts, largest first. Significance is Fisher and DISCOVER only, so Pearson/phi
# (descriptive) has no column here.
COMBINED_RENAME = {
    "n_samples": "endpoints",
    "fisher_co_sig_local": "fisher_co_sig",
    "fisher_mutex_sig_local": "fisher_mutex_sig",
    "discover_co_sig_local": "discover_co_sig",
    "discover_mutex_sig_local": "discover_mutex_sig",
}
COMBINED_ORDER = ["exp", "endpoints", "n_units_total", "n_units_tested", "n_pairs_tested",
                  "mean_mutation_burden", "median_mutation_burden", "max_mutation_burden",
                  "count_min", "fisher_co_sig", "fisher_mutex_sig",
                  "discover_co_sig", "discover_mutex_sig", "note"]


def build_combined_summary(summary_frames):
    if not summary_frames:
        return pd.DataFrame()
    df = pd.concat(summary_frames, ignore_index=True).rename(columns=COMBINED_RENAME)
    for c in COMBINED_ORDER:
        if c not in df.columns:
            df[c] = np.nan
    df = df[COMBINED_ORDER]
    sort_keys = ["discover_mutex_sig", "discover_co_sig", "fisher_mutex_sig", "fisher_co_sig"]
    return df.sort_values(sort_keys, ascending=False).reset_index(drop=True)


def build_experiment_summary(pair_df, meta, cfg, scope_name, n_input_rows, note="ok"):
    def s(col):
        return int(pair_df[col].map(lambda x: bool(x) if pd.notna(x) else False).sum()) if (col in pair_df and len(pair_df)) else 0
    row = {
        "exp": scope_name, "unit_tag": cfg.unit_tag,
        "n_input_rows": int(n_input_rows), "n_samples": meta["n_samples"],
        "n_units_total": meta["n_units_total"], "n_units_tested": meta["n_units_tested"],
        "n_pairs_tested": meta["n_pairs_tested"],
        "mean_mutation_burden": meta["mean_mutation_burden"],
        "median_mutation_burden": meta["median_mutation_burden"],
        "max_mutation_burden": meta["max_mutation_burden"],
        "count_min": meta["count_min"],
        "discover_co_sig_local": s("discover_sig_co"),
        "discover_mutex_sig_local": s("discover_sig_mutex"),
        "fisher_co_sig_local": s("fisher_sig_co"),
        "fisher_mutex_sig_local": s("fisher_sig_mutex"),
        "note": note,
    }
    return pd.DataFrame([row])


def build_significant_table(pair_df):
    if not len(pair_df):
        return pair_df
    def b(col):
        return pair_df[col].map(lambda x: bool(x) if pd.notna(x) else False) if col in pair_df else pd.Series(False, index=pair_df.index)
    keep = b("discover_sig_any") | b("fisher_sig_any")
    return pair_df[keep].copy()


def out_name(cfg, base):
    return f"{base}_{cfg.unit_tag}.csv"
