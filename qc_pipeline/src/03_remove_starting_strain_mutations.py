"""
Step 03: Remove starting-strain and known-artifact mutations.

First the rule tables are applied: config_tables.starting_strain_mutation_table
(one rule per row, all non-empty fields ANDed) and, when it is set,
config_tables.starting_strain_mutation_table_user (one mutation of one experiment
per row). Then the sample column is rewritten as "ale flask isolate tech_rep"
(without the experiment name), and in each experiment the sample equal to
step03.starting_strain_sample_label is removed together with the mutations
inherited from it. Last, rows whose Reference Seq is in step03.drop_reference_seqs
are removed, step03.gene_name_replacements rewrites the gene annotations it lists,
and sample_id (exp, "_" and sample) is added.

A later row is removed as inherited when the starting-strain sample of its
experiment has the same call: the columns of step03.strict_key_cols, compared
together (experiment, reference sequence, position and sequence change).
"""
import argparse
import importlib.util
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "03"
STEP_NAME = "remove_starting_strain_mutations"

GENE_COL = "Gene (Scrollable)"


def _split(val, cast=str):
    return [cast(x) for x in str(val).split(";") if x != ""]


# Columns every user rule must fill, and the only ones it is matched on. Their values are
# the same in the step 03 input and in every later table. 'Gene' is not matched: step 04
# writes the gene of a non-coding mutation with "/" in place of ", ", so a gene copied
# from step 04's starting-strain candidate report or a later table can differ from the
# step 03 input.
USER_REQUIRED = ["exp", "Position", "Sequence Change", "Details"]

# Columns of the curated rule table. rule_mask reads all of them except rule_id and reason.
RULE_COLUMNS = ["rule_id", "exp", "exp_set", "exp_not", "Reference_Seq", "Mutation_Type",
                "project", "Condition_normalized", "sub_species", "Position", "Position_set",
                "Sequence_Change", "SeqChange_set", "Details", "Gene", "Gene_set",
                "Gene_contains", "ale_set", "reason"]


def load_user_rules(path_cfg):
    """Load the user starting strain mutation table.

    Matches on exp + Position + Sequence Change + Details (all stable across the
    pipeline). 'Gene' (and Reference Seq / Mutation Type) are kept only for human
    readability and are NOT used for matching. A row missing a required field stops
    the run.
    """
    if not path_cfg:
        return pd.DataFrame()
    p = common.resolve(path_cfg)
    u = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for i, ur in u.iterrows():
        if not any(str(x).strip() for x in ur.values):
            continue        # an empty row, as spreadsheet exports leave, is a blank line
        missing = [c for c in USER_REQUIRED if not str(ur.get(c, "")).strip()]
        if missing:
            # Named by its content: a line number would be off after a blank line, which
            # pandas skips.
            written = ", ".join(f"{c}={ur[c]}" for c in u.columns if str(ur[c]).strip())
            raise ValueError(f"{p}: a rule has no {', '.join(missing)} ({written}). Each rule "
                             "needs exp + Position + Sequence Change + Details")
        gene_info = str(ur.get("Gene", "")).strip()
        reason = ur.get("reason", "") or "user-supplied starting strain mutation"
        rows.append({
            "rule_id": f"user_{i + 1}",
            "exp": ur["exp"].strip(),
            "Position": ur["Position"].strip(),
            "Sequence_Change": ur["Sequence Change"].strip(),
            "Details": ur["Details"].strip(),
            "reason": reason + (f" [gene: {gene_info}]" if gene_info else ""),
        })
    return pd.DataFrame(rows)


def rule_mask(df, r):
    """Build the AND mask for one rule row. Empty fields are ignored.

    Works for both the curated table (flexible columns) and the strict user
    table (exp, Position, Sequence Change and Details). Missing columns are
    treated as empty.
    """
    def g(k):
        try:
            v = r[k]
        except (KeyError, IndexError):
            return ""
        return "" if v is None else v

    # Some deposits write text fields with non-breaking spaces (\xa0). Compare them as normal
    # spaces, otherwise a line that names the right mutation silently matches nothing.
    def text(col):
        return df[col].astype(str).str.replace("\xa0", " ").str.strip()

    def plain(v):
        return str(v).replace("\xa0", " ").strip()

    mask = pd.Series(True, index=df.index)
    if g("exp"):
        mask &= df["exp"] == g("exp")
    if g("exp_set"):
        mask &= df["exp"].isin(_split(g("exp_set")))
    if g("exp_not"):
        mask &= df["exp"] != g("exp_not")
    if g("Reference_Seq"):
        mask &= df["Reference Seq"] == g("Reference_Seq")
    if g("Mutation_Type"):
        mask &= df["Mutation Type"] == g("Mutation_Type")
    if g("project"):
        mask &= df["project"] == g("project")
    if g("Condition_normalized"):
        mask &= df["Condition_normalized"] == g("Condition_normalized")
    if g("sub_species"):
        mask &= df["sub_species"] == g("sub_species")
    if g("Position"):
        mask &= df["Position"] == int(g("Position"))
    if g("Position_set"):
        mask &= df["Position"].isin(_split(g("Position_set"), int))
    if g("Sequence_Change"):
        mask &= text("Sequence Change") == plain(g("Sequence_Change"))
    if g("SeqChange_set"):
        mask &= text("Sequence Change").isin([plain(x) for x in _split(g("SeqChange_set"))])
    if g("Details"):
        mask &= text("Details") == plain(g("Details"))
    if g("Gene"):
        mask &= df[GENE_COL] == g("Gene")
    if g("Gene_set"):
        mask &= df[GENE_COL].isin(_split(g("Gene_set")))
    if g("Gene_contains"):
        mask &= df[GENE_COL].str.contains(g("Gene_contains"), na=False)
    if g("ale_set"):
        mask &= df["ale"].isin(_split(g("ale_set"), int))
    return mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    in_path = common.out_dir() / "02_species_selected.tsv"
    out_path = common.out_dir() / "03_starting_strain_removed.tsv"
    curated_path = common.resolve(cfg["config_tables"]["starting_strain_mutation_table"])
    user_path_cfg = cfg["config_tables"]["starting_strain_mutation_table_user"]

    muts_df = common.read_table(in_path)
    muts_df["exp ALE"] = muts_df.apply(lambda m: str(m["exp"]) + " " + str(m["ale"]), axis=1)
    muts_df = muts_df.reset_index(drop=True)
    n_before = len(muts_df)
    logger.info("loaded %d rows from %s", n_before, in_path)
    # The "matches nothing" checks below test the step input, so a value whose rows an earlier
    # removal already took is not mistaken for a typo.
    step_input = muts_df

    # --- assemble rules: curated table (flexible) + user table (strict, validated) ---
    rules = pd.read_csv(curated_path, sep="\t", dtype=str, keep_default_na=False)
    # A column rule_mask does not know would be ignored, which widens every rule that uses it,
    # and a rule with no field filled in would match every row.
    unknown = [c for c in rules.columns if c not in RULE_COLUMNS]
    if unknown:
        raise ValueError(f"{curated_path}: unknown columns {unknown}. The rule columns are "
                         + ", ".join(RULE_COLUMNS))
    fields = [c for c in rules.columns if c not in ("rule_id", "reason")]
    empty = rules.loc[(rules[fields].apply(lambda s: s.str.strip()) == "").all(axis=1), "rule_id"]
    if len(empty):
        raise ValueError(f"{curated_path}: rules {list(empty)} have no field filled in, so they "
                         "would remove every row")

    user_rules = load_user_rules(user_path_cfg)
    n_user = len(user_rules)
    user_ids = set(user_rules["rule_id"]) if n_user else set()
    if n_user:
        rules = pd.concat([rules, user_rules], ignore_index=True).fillna("")
    logger.info("loaded %d curated rules + %d valid user rules", len(rules) - n_user, n_user)

    removed_parts = []
    removed_by_rule = {}
    rules_no_match = []

    # --- apply rules ---
    for _, r in rules.iterrows():
        m = rule_mask(muts_df, r)
        n = int(m.sum())
        if n:
            part = muts_df[m].copy()
            part["removal_reason"] = f"rule {r['rule_id']}: {r['reason']}"
            removed_parts.append(part)
            removed_by_rule[r["rule_id"]] = n
            muts_df = muts_df[~m]
        else:
            rules_no_match.append(r["rule_id"])
    logger.info("rule table removed %d rows across %d active rules",
                sum(removed_by_rule.values()), len(removed_by_rule))
    # A user rule names one mutation of one experiment, so a user rule whose mutation is not in
    # the step input is written wrongly. One whose rows an earlier rule took is only redundant.
    user_no_match = rules[rules["rule_id"].isin(set(rules_no_match) & user_ids)]
    absent = user_no_match[[not rule_mask(step_input, r).any() for _, r in user_no_match.iterrows()]]
    if len(absent):
        raise ValueError(
            f"{len(absent)} rule(s) of the user starting-strain table "
            f"({common.resolve(user_path_cfg)}, config_tables.starting_strain_mutation_table_user) "
            "name a mutation that is not in the data, for example:\n"
            + absent[["exp", "Position", "Sequence_Change", "Details"]].head(5).to_string(index=False))
    if rules_no_match:
        # A rule can match nothing because an earlier rule already removed its rows, and a
        # curated rule also because its mutation is absent from the data or it is written wrongly.
        logger.warning("%d rule(s) matched nothing: %s", len(rules_no_match), rules_no_match)

    # --- generic per-experiment starting-strain removal ---
    label = cfg["step03"]["starting_strain_sample_label"]
    muts_df = muts_df.reset_index(drop=True)
    muts_df["sample"] = muts_df.apply(
        lambda r: f"{int(r.ale)} {int(r.flask)} {int(r.isolate)} {int(r.tech_rep)}", axis=1)

    starting_strains = muts_df[muts_df["sample"] == label]
    if starting_strains.empty:
        raise ValueError(f"step03.starting_strain_sample_label {label!r} matches no sample")
    ss_rows = starting_strains.copy()
    ss_rows["removal_reason"] = f"starting-strain sample ({label})"
    removed_parts.append(ss_rows)
    muts_df = muts_df[muts_df["sample"] != label]

    n_loop_removed = 0
    key_cols = cfg["step03"]["strict_key_cols"]
    ss_keys = set(map(tuple, starting_strains[key_cols].astype(str).itertuples(index=False, name=None)))
    tup = muts_df[key_cols].astype(str).itertuples(index=False, name=None)
    m = pd.Series([t in ss_keys for t in tup], index=muts_df.index)
    if m.any():
        part = muts_df[m].copy()
        part["removal_reason"] = "inherited starting-strain mutation"
        removed_parts.append(part)
        n_loop_removed = int(m.sum())
        muts_df = muts_df[~m]
    logger.info("starting-strain loop (key %s): removed %d sample rows + %d inherited rows",
                key_cols, len(starting_strains), n_loop_removed)

    # --- drop excluded reference sequences ---
    drop_refs = cfg["step03"]["drop_reference_seqs"]
    input_refs = set(step_input["Reference Seq"])
    absent = [r for r in drop_refs if r not in input_refs]
    if absent:
        raise ValueError(f"step03.drop_reference_seqs: no row has the Reference Seq {absent}")
    if drop_refs:
        m = muts_df["Reference Seq"].isin(drop_refs)
        if m.any():
            part = muts_df[m].copy()
            part["removal_reason"] = "Reference Seq in drop list " + str(drop_refs)
            removed_parts.append(part)
            logger.info("dropped %d rows with Reference Seq in %s", int(m.sum()), drop_refs)
            muts_df = muts_df[~m]

    # --- gene-name string replacements (annotation harmonization, no row change) ---
    replacements = cfg["step03"]["gene_name_replacements"]
    input_genes = set(step_input[GENE_COL])
    absent = [old for old in replacements if old not in input_genes]
    if absent:
        raise ValueError(f"step03.gene_name_replacements: no row has the gene {absent}")
    for old, new in replacements.items():
        muts_df[GENE_COL] = muts_df[GENE_COL].apply(lambda a: new if a == old else a)

    # --- globally unique sample id ---
    muts_df = muts_df.reset_index(drop=True)
    muts_df["sample_id"] = muts_df["exp"].astype(str) + "_" + muts_df["sample"].astype(str)

    removed_all = pd.concat(removed_parts, ignore_index=True) if removed_parts else muts_df.iloc[0:0].copy()
    n_removed = len(removed_all)
    common.check_conservation(n_before, len(muts_df), n_removed, logger)
    common.validate_common(muts_df, logger, required_cols=["exp", "Position", "sample_id"])

    common.write_table(muts_df, out_path)
    common.save_removed(removed_all, common.removed_dir() / "03_removed_starting_strain_mutations.tsv",
                        "starting-strain / artifact removal")
    logger.info("wrote %s with shape %s, removed %d rows total", out_path, muts_df.shape, n_removed)

    n_removed_by_exp = removed_all["exp"].value_counts().to_dict() if n_removed else {}
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_file": str(in_path), "output_file": str(out_path),
        "n_rows_before": n_before, "n_rows_after": len(muts_df), "n_removed": n_removed,
        "n_removed_by_experiment": n_removed_by_exp,
        "n_removed_unique_mutations": int(
            removed_all[key_cols].drop_duplicates().shape[0]) if n_removed else 0,
        "matching_key_cols": key_cols,
        "n_rules_active": len(removed_by_rule),
        "n_removed_by_rule_table": int(sum(removed_by_rule.values())),
    }
    common.write_checkpoint(metrics, common.checkpoints_dir() / "03_remove_starting_strain_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
