"""
Step 05: Normalize gene names into test units.

Unit construction, in order:

  1. INTERGENIC calls are assigned to the nearer flanking gene, using the
     "intergenic (d1/d2)" distances in Details, so every test unit is a gene. A call
     whose distances cannot be read is removed as unresolved.
  2. Identifiers are CLEANED: a trailing " CDS" is stripped, a duplicated "X|X" is
     collapsed to "X", and a leading "[" and a trailing "]" are dropped. None of these
     forms names a distinct gene: some references write gene names with a " CDS"
     suffix, "X|X" is what breseq emits when a position falls inside two overlapping
     feature records, and breseq puts a gene in brackets when a mutation covers only
     part of it.
  3. Calls are mapped to standard gene names via the STRING alias table.
  4. An identifier STRING cannot name is KEPT under its own annotation, not discarded,
     when it holds no "," or "/" (so a pair such as "X|Y" stays one unit). Gene lists
     joined by "," or "/" cannot be named as one gene and are removed with the unmapped
     rows.

Steps 2 and 4 are config-gated under step05 (clean_gene_strings and
keep_unmapped_as_annotated). With clean_gene_strings false, identifiers are used as
annotated. With keep_unmapped_as_annotated false, every identifier STRING cannot name
is removed.

Rows of the projects listed in step05.remove_projects are removed after the mapping.
A unique locus tag (the b-number for E. coli, in the column step05.locus_tag_column) is
attached only to STRING-mapped units.
"""
import argparse
import gzip
import importlib.util
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "05"
STEP_NAME = "normalize_gene_names"
GENE_COL = "Gene (Scrollable)"


def clean_token(s):
    """Clean ONE gene identifier: drop a trailing ' CDS', collapse 'X|X', trim brackets."""
    s = str(s).strip()
    s = re.sub(r"\s+CDS$", "", s, flags=re.IGNORECASE)
    m = re.match(r"^([^|]+)\|\1$", s)          # X|X -> X (one feature annotated twice)
    if m:
        s = m.group(1)
    return re.sub(r"^\[|\]$", "", s).strip()


def assign_closest_gene(row):
    """Assign an intergenic mutation to the closer flanking gene."""
    details = str(row["Details"])
    gene_col = str(row[GENE_COL])
    if "intergenic" not in details.lower():
        return gene_col
    genes = gene_col.split("/")
    if len(genes) != 2:
        return gene_col
    match = re.search(r"intergenic\s*\(\s*([^/]+)/([^)]+)\)", details, re.IGNORECASE)
    if not match:
        return gene_col
    try:
        d1 = match.group(1).strip().replace("‑", "-").replace("–", "-").replace("−", "-")
        d2 = match.group(2).strip().replace("‑", "-").replace("–", "-").replace("−", "-")
        dist1, dist2 = abs(float(d1)), abs(float(d2))
    except ValueError:
        return gene_col
    return genes[0].strip() if dist1 <= dist2 else genes[1].strip()


def build_string_maps(string_path, valid_sources, primary_sources):
    string_mapping, string_mapping_priority, string_id_to_info = {}, {}, {}
    seen_sources = set()
    with gzip.open(string_path, "rt") as f:
        next(f)
        for line in f:
            data = line.rstrip().split("\t")
            string_id, alias = data[0], data[1]
            sources = set(data[2:])
            seen_sources |= sources
            if len(valid_sources & sources) == 0:
                continue
            is_primary = len(primary_sources & sources) > 0
            if alias not in string_mapping:
                string_mapping[alias] = string_id
                string_mapping_priority[alias] = is_primary
            elif is_primary and not string_mapping_priority.get(alias, False):
                string_mapping[alias] = string_id
                string_mapping_priority[alias] = True
            if string_id not in string_id_to_info:
                string_id_to_info[string_id] = {"b_tag": None, "gene_name": None}
            # gene_name = first primary alias that is not a b-number-looking token
            if not re.match(r"^b[0-9]+$", alias):
                if is_primary and string_id_to_info[string_id]["gene_name"] is None:
                    string_id_to_info[string_id]["gene_name"] = alias

    # Locus tag read directly from the STRING protein id (<taxid>.<locus_tag>), for any species
    # (e.g. 511145.b0259 -> b0259, 208964.PA0001 -> PA0001).
    for sid in string_id_to_info:
        string_id_to_info[sid]["b_tag"] = sid.split(".", 1)[1] if "." in sid else sid

    valid_gene_names = {info["gene_name"] for info in string_id_to_info.values() if info["gene_name"]}
    alias_to_standard = {}
    for sid, info in string_id_to_info.items():
        std = info["gene_name"]
        if std is None:
            continue
        alias_to_standard[std] = std
        if info["b_tag"]:
            alias_to_standard[info["b_tag"]] = std
    for alias, sid in string_mapping.items():
        std = string_id_to_info.get(sid, {}).get("gene_name")
        if std and alias not in alias_to_standard:
            alias_to_standard[alias] = std
    absent = sorted((valid_sources | primary_sources) - seen_sources)
    if absent:
        raise ValueError(f"{string_path} has no alias from the source(s) {absent} "
                         "(step05.string_valid_sources, step05.string_primary_sources)")
    return string_mapping, string_id_to_info, valid_gene_names, alias_to_standard


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    s5 = cfg["step05"]
    in_path = common.out_dir() / "04_hypermutators_removed.tsv"
    out_path = common.out_dir() / "05_gene_names_normalized.tsv"
    string_path = common.resolve(cfg["inputs"]["string_aliases"])

    df = common.read_table(in_path)
    n_before = len(df)
    absent = [p for p in s5["remove_projects"] if p not in set(df["project"])]
    if absent:
        raise ValueError(f"step05.remove_projects: no row has the project {absent}")
    n_genes_before = df[GENE_COL].nunique()
    logger.info("loaded %d rows from %s", n_before, in_path)

    removed_parts = []

    clean_strings = s5["clean_gene_strings"]
    keep_unmapped = s5["keep_unmapped_as_annotated"]

    # --- assign each intergenic mutation to its closer flanking gene ---
    df["Intergenic_process_gene_name"] = df.apply(assign_closest_gene, axis=1)
    if clean_strings:
        df["Intergenic_process_gene_name"] = df["Intergenic_process_gene_name"].map(clean_token)
    df["is_intergenic"] = df["Details"].str.contains("intergenic", case=False, na=False)
    mask_unresolved = (
        df["Details"].str.contains("intergenic", case=False, na=False)
        & df["Intergenic_process_gene_name"].str.contains("/", na=False))
    if mask_unresolved.any():
        part = df[mask_unresolved].copy()
        part["removal_reason"] = "unresolved intergenic (could not split to a single gene)"
        removed_parts.append(part)
    df = df[~mask_unresolved].reset_index(drop=True)
    logger.info("removed %d unresolved intergenic rows", int(mask_unresolved.sum()))

    # --- build STRING alias -> standard map ---
    valid_sources = set(s5["string_valid_sources"])
    primary_sources = set(s5["string_primary_sources"])
    if not primary_sources <= valid_sources:
        raise ValueError(f"step05.string_primary_sources {sorted(primary_sources - valid_sources)} "
                         "must also be in step05.string_valid_sources, lines of other sources are "
                         "skipped before the primary test")
    string_mapping, string_id_to_info, valid_gene_names, alias_to_standard = build_string_maps(
        string_path, valid_sources, primary_sources)
    logger.info("STRING reference: %d valid gene names, %d alias mappings",
                len(valid_gene_names), len(alias_to_standard))

    def normalize_gene_name(gene):
        if pd.isna(gene):
            return None
        g = gene.strip()
        if g in valid_gene_names:
            return g
        return alias_to_standard.get(g, None)

    df["Gene_name_normalized"] = df["Intergenic_process_gene_name"].apply(normalize_gene_name)
    df["mapped_by_string"] = df["Gene_name_normalized"].notna()
    unmapped_mask = df["Gene_name_normalized"].isna()

    # Keep a single-gene identifier that STRING cannot name, rather than dropping it.
    n_kept_raw = 0
    if keep_unmapped:
        gid = df["Intergenic_process_gene_name"].astype(str)
        single = ~gid.str.contains(",", regex=False) & ~gid.str.contains("/", regex=False)
        rec = unmapped_mask & single
        df.loc[rec, "Gene_name_normalized"] = gid[rec]
        n_kept_raw = int(rec.sum())
        logger.info("kept %d unmapped single-gene rows under their own identifier", n_kept_raw)
        unmapped_mask = df["Gene_name_normalized"].isna()

    n_unmapped = int(unmapped_mask.sum())

    unmapped_report = df[unmapped_mask]["Intergenic_process_gene_name"].value_counts().reset_index()
    unmapped_report.columns = ["Gene_name", "Count"]
    # Equal counts are ordered by name, so the row order is fixed.
    unmapped_report = unmapped_report.sort_values(["Count", "Gene_name"], ascending=[False, True])
    common.write_table(unmapped_report, common.qc_dir() / "05_unmapped_gene_names.csv")
    if unmapped_mask.any():
        part = df[unmapped_mask].copy()
        part["removal_reason"] = "gene name could not be mapped to a standard name"
        removed_parts.append(part)
    df = df[~unmapped_mask].copy().reset_index(drop=True)
    logger.info("removed %d unmapped-gene rows, %d rows remain", n_unmapped, len(df))

    # --- remove the projects listed in step05.remove_projects ---
    for proj in s5["remove_projects"]:
        pm = df["project"] == proj
        part = df[pm].copy()
        part["removal_reason"] = f"project removed during gene normalization: {proj}"
        removed_parts.append(part)
        logger.info("removed %d rows from project %s", int(pm.sum()), proj)
        df = df[~pm].copy()

    # --- locus-tag assignment (unique only: multi-copy genes, such as IS genes, get none) ---
    gn_to_bnums = defaultdict(set)
    for info in string_id_to_info.values():
        if info["gene_name"] and info["b_tag"]:
            gn_to_bnums[info["gene_name"]].add(info["b_tag"])
    multicopy = {g for g, b in gn_to_bnums.items() if len(b) > 1}
    gene_name_to_bnum = {g: next(iter(b)) for g, b in gn_to_bnums.items() if len(b) == 1}
    locus_col = s5["locus_tag_column"]
    df = df.copy()
    # Only STRING-mapped units get a locus tag. Kept-raw identifiers have none, and inheriting
    # one would imply a resolution that was never made.
    df[locus_col] = df["Gene_name_normalized"].where(df["mapped_by_string"]).map(gene_name_to_bnum)

    # mapping diagnostics
    mapping_tbl = (df[["Intergenic_process_gene_name", "Gene_name_normalized", locus_col]]
                   .drop_duplicates().sort_values("Gene_name_normalized"))
    common.write_table(mapping_tbl, common.qc_dir() / "05_gene_name_mapping_table.tsv")
    removed_all = pd.concat(removed_parts, ignore_index=True) if removed_parts else df.iloc[0:0].copy()
    common.save_removed(removed_all, common.removed_dir() / "05_removed_multigene_or_unmapped_rows.tsv",
                        "unmapped / unresolved / project-removed")

    n_after = len(df)
    common.check_conservation(n_before, n_after, len(removed_all), logger)
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_file": str(in_path), "output_file": str(out_path),
        "n_rows_before": n_before, "n_rows_after": n_after, "n_removed": n_before - n_after,
        "n_gene_names_before": int(n_genes_before),
        "n_gene_names_after": int(df["Gene_name_normalized"].nunique()),
        "n_mapped_by_alias": int((df["Intergenic_process_gene_name"] != df["Gene_name_normalized"]).sum()),
        "n_unmapped": n_unmapped,
        "n_unresolved_intergenic_removed": int(mask_unresolved.sum()),
        "clean_gene_strings": bool(clean_strings),
        "keep_unmapped_as_annotated": bool(keep_unmapped),
        "n_rows_intergenic": int(df["is_intergenic"].sum()),
        "n_intergenic_units": int(df.loc[df["is_intergenic"], "Gene_name_normalized"].nunique()),
        "n_mapped_by_string": int(df["mapped_by_string"].sum()),
        "n_kept_raw_unmapped": n_kept_raw,
        "unmapped_rows_by_gene_string": unmapped_report.set_index("Gene_name")["Count"].to_dict(),
        "n_multicopy_genes_without_bnumber": len(multicopy),
        "n_rows_with_bnumber": int(df[locus_col].notna().sum()),
        "locus_tag_column": locus_col,
    }
    common.validate_common(df, logger, required_cols=["exp", "Position", "Gene_name_normalized"])
    common.write_table(df, out_path)
    common.write_checkpoint(metrics, common.checkpoints_dir() / "05_gene_name_normalization_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
