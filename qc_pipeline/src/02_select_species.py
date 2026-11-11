"""
Step 02: Select the configured species and attach normalized conditions.

Keeps the rows of organism.species (and organism.subspecies when it is set),
attaches Condition_normalized / Condition_group via a two-step merge (first on the
columns of step02.match_cols, then on exp alone for rows that did not match), KEEPS
rows that cannot be condition-mapped (blank condition), and removes only the
condition groups listed in step02.remove_condition_groups.
The condition table must give every full key one label, since a key with two
labels would copy each matching mutation row, so the step stops if it does not.
"""
import argparse
import importlib.util
import re
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "02"
STEP_NAME = "select_species"
PLACEHOLDER = "__MISSING__"


def clean_column(series):
    """Strip whitespace, drop trailing .0, fill empties with the merge placeholder."""
    def clean_val(x):
        if pd.isna(x):
            return PLACEHOLDER
        if isinstance(x, str):
            x = x.strip()
            if x == "":
                return PLACEHOLDER
            x = re.sub(r"\.0$", "", x)
            return x
        if isinstance(x, float) and x == int(x):
            return str(int(x))
        return x
    return series.apply(clean_val)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    in_path = common.out_dir() / "01_aledb_mutation_df.tsv"
    out_path = common.out_dir() / "02_species_selected.tsv"
    cond_path = common.resolve(cfg["inputs"]["condition_normalization"])

    logger.info("reading input mutations: %s", in_path)
    muts_df = common.read_table(in_path)
    supplement_norm = common.read_table(cond_path)

    muts_df["exp ALE"] = muts_df.apply(lambda m: str(m["exp"]) + " " + str(m["ale"]), axis=1)

    n_before = len(muts_df)
    n_exp_before = muts_df[["species_norm", "exp"]].drop_duplicates().shape[0]
    org = cfg["organism"]
    species_value = org["species"]
    subspecies_sel = org["subspecies"]

    # --- species selection ---
    all_ecoli = muts_df[muts_df["species_norm"] == species_value].copy()
    if all_ecoli.empty:
        raise ValueError(f"no row has species_norm {species_value!r} (organism.species)")
    # A subspecies the data does not have would select no row, so stop before any row is removed.
    absent = [s for s in subspecies_sel if s not in set(all_ecoli["sub_species"])]
    if absent:
        raise ValueError(f"organism.subspecies: no {species_value} row has the sub_species {absent}")
    non_ecoli = muts_df[muts_df["species_norm"] != species_value]
    n_removed_non_ecoli = len(non_ecoli)
    common.save_removed(non_ecoli, common.removed_dir() / "02_removed_non_species.tsv",
                        "species_norm is not " + str(species_value))
    logger.info("species filter: kept %d of %d rows (removed %d not '%s')",
                len(all_ecoli), n_before, n_removed_non_ecoli, species_value)

    # --- optional subspecies selection (organism.subspecies) ---
    n_removed_sub = 0
    if subspecies_sel:
        before_sub = len(all_ecoli)
        in_sub = all_ecoli["sub_species"].isin(subspecies_sel)
        removed_sub = all_ecoli[~in_sub]
        n_removed_sub = len(removed_sub)
        common.save_removed(removed_sub, common.removed_dir() / "02_removed_other_subspecies.tsv",
                            "sub_species not in " + str(subspecies_sel))
        all_ecoli = all_ecoli[in_sub].copy()
        logger.info("subspecies filter %s: kept %d of %d rows", subspecies_sel, len(all_ecoli), before_sub)
    else:
        # No subspecies filter. The table is still written, empty, so a run always has it.
        common.save_removed(all_ecoli.iloc[0:0], common.removed_dir() / "02_removed_other_subspecies.tsv",
                            "no subspecies filter")

    # --- two-step condition merge ---
    match_cols = cfg["step02"]["match_cols"]
    for col in match_cols:
        all_ecoli[col] = clean_column(all_ecoli[col])
        supplement_norm[col] = clean_column(supplement_norm[col])
    supplement_norm["Condition_normalized"] = (
        supplement_norm["Condition_normalized"].astype(str).str.replace(r"\.0$", "", regex=True)
        .str.strip())
    supplement_cols = match_cols + ["Condition_normalized", "Condition_group"]

    # A key with more than one label would copy every mutation row that matches it.
    condition_keys = supplement_norm[supplement_cols].drop_duplicates()
    conflicts = condition_keys[condition_keys.duplicated(subset=match_cols, keep=False)]
    if len(conflicts):
        raise ValueError(
            f"condition table {cond_path} gives {conflicts.groupby(match_cols).ngroups} key(s) more "
            f"than one label (key = {match_cols}). Keep one row per key:\n"
            + conflicts.sort_values(match_cols).to_string(index=False))

    condition_group_mut_df = all_ecoli.merge(
        condition_keys, on=match_cols, how="left", indicator=True)

    left_only_mask = condition_group_mut_df["_merge"] == "left_only"
    logger.info("merge 1: matched=%d unmatched=%d",
                (~left_only_mask).sum(), left_only_mask.sum())
    unmatched = condition_group_mut_df[left_only_mask].drop(
        columns=["Condition_normalized", "Condition_group", "_merge"])
    matched = condition_group_mut_df[~left_only_mask].drop(columns=["_merge"])

    unmatched_remapped = unmatched.merge(
        supplement_norm[["exp", "Condition_normalized", "Condition_group"]].drop_duplicates(
            subset="exp", keep="first"),
        on=["exp"], how="left", indicator=True)
    logger.info("merge 2 (exp-only fallback): matched=%d unmatched=%d",
                (unmatched_remapped["_merge"] != "left_only").sum(),
                (unmatched_remapped["_merge"] == "left_only").sum())

    combined = pd.concat([matched, unmatched_remapped], ignore_index=True)
    # Rows that match no condition are kept, with a blank Condition_normalized and
    # Condition_group, and their number is logged. Only the groups in
    # step02.remove_condition_groups are removed below.
    if "_merge" in combined.columns:
        n_unmapped = int((combined["_merge"] == "left_only").sum())
        if n_unmapped:
            logger.info("condition-unmapped rows kept (blank condition): %d", n_unmapped)
    combined = combined.drop(columns=["_merge"], errors="ignore")
    # The placeholder only lets a blank cell match a blank cell of the condition table. The
    # blanks are put back, so it does not reach the output.
    combined[match_cols] = combined[match_cols].mask(combined[match_cols] == PLACEHOLDER)

    # --- remove configured condition groups ---
    remove_groups = cfg["step02"]["remove_condition_groups"]
    absent = [g for g in remove_groups if g not in set(combined["Condition_group"])]
    if absent:
        raise ValueError(f"step02.remove_condition_groups: no selected {species_value} row has the "
                         f"Condition_group {absent}")
    drop_mask = combined["Condition_group"].isin(remove_groups)
    n_removed_cond = int(drop_mask.sum())
    common.save_removed(combined[drop_mask],
                        common.removed_dir() / "02_removed_condition_groups.tsv",
                        "Condition_group in " + str(remove_groups))
    selected = combined[~drop_mask].copy()
    common.check_conservation(n_before, len(selected),
                              n_removed_non_ecoli + n_removed_sub + n_removed_cond, logger)
    common.validate_common(selected, logger, required_cols=["exp", "Position", "Condition_group"])

    common.write_table(selected, out_path)
    logger.info("wrote %s with shape %s", out_path, selected.shape)

    # --- checkpoint metrics ---
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_file": str(in_path), "output_file": str(out_path),
        "n_rows_before": n_before, "n_rows_after": len(selected),
        "n_removed_non_ecoli": n_removed_non_ecoli,
        "n_removed_other_subspecies": n_removed_sub,
        "n_removed_condition_groups": n_removed_cond,
        "n_unique_experiments_before": int(n_exp_before),
        "n_unique_experiments_after": int(selected["exp"].nunique()),
        "species_counts": muts_df["species_norm"].value_counts(dropna=False).to_dict(),
        "subspecies_counts": selected["sub_species"].value_counts(dropna=False).to_dict(),
        "condition_group_counts": selected["Condition_group"].value_counts(dropna=False).to_dict(),
        "unclear_condition_count": int((selected["Condition_group"] == "Unclear").sum()),
        "gene_deletion_condition_count": int((selected["Condition_group"] == "Gene_deletion").sum()),
    }
    common.write_checkpoint(metrics, common.checkpoints_dir() / "02_select_species_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
