"""
Step 01: Load the ALEdb mutation dataframe.

Loads raw per-experiment breseq mutation CSVs (via mutil), builds basic sample
identifiers, reads positions as whole numbers, merges experiment metadata and
species normalization, and attaches subspecies / gene-deletion annotations.
Rows of every species are kept: step 02 selects the species. While loading, rows
without a gene name and rows whose Details is "Duplication" are removed, and sample
columns whose name repeats another column (name.1) are dropped. The log counts the
calls in those columns that are not loaded.
Each ALEdb deposit (ExpID) has its own species: the one the species table gives the
ExpID of its metadata file (source_file). A mutation row takes the metadata and the
species of its own deposit, so deposits that share an experiment name can be of
different species.
A deposit whose mutation file has no rows is skipped, and of the other deposits only one
per experiment name and species is loaded: the one with a DOI, then the smallest ExpID,
unless step01.deposit_overrides names one. Names are compared exactly.
The files of each folder are read in sorted order of their names (by character code), so
the row order of the tables does not depend on the order in which the file system lists
them.
The step stops, naming what to fix, when a deposit with mutation rows has no species in
the species table, when a row of the configured species has no metadata row or no
subspecies, or when the species or subspecies table gives one key two different values.
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

STEP_ID = "01"
STEP_NAME = "load_aledb_mutations"

REQUIRED_COLUMNS = ["exp", "ale", "flask", "isolate", "tech_rep", "presence",
                    "Position", "Mutation Type", "Sequence Change", "Details",
                    "Gene (Scrollable)", "Reference Seq"]

META_FILE_RE = re.compile(r"^md_Proj_(.*?)_Exp_(.*)_ExpID(\d+)\.csv$")
MUTS_FILE_RE = re.compile(r"_ExpID(\d+)_mut\.csv$")


def get_exp_sample_name(exp_name, ale, flask, isolate, tech_rep):
    return f"{exp_name} {int(ale)} {int(flask)} {int(isolate)} {int(tech_rep)}"


def concat_csv_files(folder_path, file_names):
    folder = Path(folder_path)
    all_dfs = []
    for csv_file in sorted(folder.glob("*.csv")):
        if csv_file.name not in file_names:
            continue
        d = pd.read_csv(csv_file)
        d["source_file"] = csv_file.name
        all_dfs.append(d)
    if not all_dfs:
        raise ValueError(f"No CSV files found in {folder_path} among the selected file names")
    return pd.concat(all_dfs, ignore_index=True)


def count_rows(path):
    """The mutation rows of a deposit: the lines of its mutation file after the header."""
    with open(path, "rb") as fh:
        return max(sum(1 for line in fh if line.strip()) - 1, 0)


def expid_of(meta_file, where):
    """The ExpID in the name of a metadata file. where is the table or folder that holds the name,
    shown in the error message."""
    m = META_FILE_RE.match(str(meta_file))
    if m is None:
        raise ValueError(f"{where}: cannot parse the ExpID of metadata file name {meta_file!r}")
    return int(m.group(3))


def select_deposits(meta_dir, muts_dir, overrides, species):
    """Pick one ALEdb deposit (ExpID) per experiment name and species, names compared exactly.

    The later steps tell experiments apart by name, so deposits of one species that
    share a name would be read as one experiment. species gives the species of each
    ExpID. A deposit whose mutation file has no rows is skipped before the deposits are
    grouped, so its kept_expid is empty. A deposit with mutation rows needs a species.
    Of each name and species, keep the deposit with a DOI, then the smallest ExpID, or the
    ExpID given in overrides. Returns every deposit, with its species in species_norm and
    the ExpID kept for its name and species in kept_expid.
    """
    rows = []
    for f in sorted(Path(meta_dir).glob("*.csv")):
        m = META_FILE_RE.match(f.name)
        if m is None:
            raise ValueError(f"cannot parse the ExpID of metadata file {f.name}")
        # Read as text, so an empty cell stays empty instead of becoming the word "nan".
        md = pd.read_csv(f, dtype=str, keep_default_na=False)
        names = [v for v in md["experiment"].unique() if v != ""]
        # a header-only metadata file has no experiment value, so take the name from the file name
        exp = names[0].replace(" ", "_") if len(names) else m.group(2)
        doi = sorted({v.strip() for v in md["doi"] if v.strip()})
        rows.append({"exp": exp, "expid": int(m.group(3)), "meta_file": f.name,
                     "doi": "; ".join(doi)})
    if not rows:
        raise ValueError(f"No CSV files found in {meta_dir}")
    deposits = pd.DataFrame(rows)
    # Every mutation file needs its metadata file and the reverse, otherwise a deposit
    # would be left out without a word.
    muts = {int(MUTS_FILE_RE.search(f.name).group(1)): f for f in Path(muts_dir).glob("*.csv")}
    unpaired = sorted(set(muts) ^ set(deposits["expid"]))
    if unpaired:
        raise ValueError(f"ExpIDs {unpaired[:10]} have a mutation file in {muts_dir} or a "
                         f"metadata file in {meta_dir}, but not both")
    deposits["mutation_rows"] = deposits["expid"].map(lambda e: count_rows(muts[e]))
    deposits["species_norm"] = deposits["expid"].map(species)
    has_rows = deposits["mutation_rows"] > 0
    no_species = deposits[deposits["species_norm"].isna() & has_rows]
    if len(no_species):
        raise ValueError(f"{len(no_species)} deposit(s) with mutation rows have no species in the species "
                         "table (no row for the ExpID of their metadata file, or an empty species). Add "
                         "one row per deposit, with the metadata file as source_file:\n"
                         + no_species.sort_values(["exp", "expid"])[["exp", "expid", "meta_file"]]
                         .to_string(index=False))
    # A deposit with no mutation rows is in no group, so its kept_expid is empty.
    ranked = deposits[has_rows].assign(has_doi=deposits["doi"] != "")
    ranked = ranked.sort_values(["has_doi", "expid"], ascending=[False, True])
    deposits["kept_expid"] = (ranked.groupby(["exp", "species_norm"])["expid"].transform("first")
                              .astype("Int64"))
    deposits.loc[has_rows, "kept_by"] = "DOI, then smallest ExpID"
    for exp, expid in overrides.items():
        with_rows = deposits.loc[(deposits["exp"] == exp) & has_rows, "expid"]
        if int(expid) not in set(with_rows):
            raise ValueError(f"deposit_overrides: experiment {exp} has no ExpID {expid} of that name with "
                             "mutation rows")
        # It replaces the deposit picked for its name and species, so it applies to the deposits
        # that kept the same ExpID as it.
        mask = deposits["kept_expid"] == deposits.loc[deposits["expid"] == int(expid), "kept_expid"].iloc[0]
        deposits.loc[mask, "kept_expid"] = int(expid)
        deposits.loc[mask, "kept_by"] = "deposit_overrides"
    return deposits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    common.add_mutil_to_path(cfg)
    from mutil.alemutdf import get_all_sample_mut_df

    inputs = cfg["inputs"]
    muts_dir = common.resolve(inputs["muts_all_dir"])
    meta_dir = common.resolve(inputs["metadata_dir"])
    species_path = common.resolve(inputs["species_normalization"])
    out_path = common.out_dir() / "01_aledb_mutation_df.tsv"

    org = cfg["organism"]
    subspecies_cfg = org["subspecies_table"]
    common.organism_input_warning(cfg, logger)

    # --- species of each deposit ---
    species_norm = common.read_table(species_path)
    species_norm["species_norm"] = species_norm["species_norm"].str.rstrip()
    species_norm["sub_species"] = species_norm["sub_species"].str.rstrip()
    # A cell that holds only spaces counts as empty.
    for col in ("species_norm", "sub_species"):
        species_norm[col] = species_norm[col].mask(species_norm[col] == "")
    # A deposit is found by the ExpID of its metadata file (source_file).
    species_norm["_expid"] = species_norm["source_file"].map(lambda f: expid_of(f, species_path))
    # One species and subspecies per deposit. With two, every mutation row of that deposit would be copied.
    deposit_species = species_norm[["_expid", "species_norm", "sub_species"]].drop_duplicates()
    conflicts = species_norm[species_norm["_expid"].isin(
        deposit_species.loc[deposit_species["_expid"].duplicated(), "_expid"])]
    if len(conflicts):
        shown = (conflicts[["source_file", "species_norm", "sub_species"]].drop_duplicates()
                 .sort_values("source_file").assign(species_norm=lambda d: d["species_norm"].map(repr),
                                                    sub_species=lambda d: d["sub_species"].map(repr)))
        raise ValueError(f"{species_path} gives {conflicts['_expid'].nunique()} deposit(s) (the ExpID of "
                         "source_file) more than one species or subspecies (values in quotes, so spaces "
                         "show). Keep one value per deposit:\n" + shown.to_string(index=False))
    species_of = deposit_species.set_index("_expid")["species_norm"]

    # --- one deposit per experiment name and species ---
    step_cfg = cfg["step01"]
    deposits = select_deposits(meta_dir, muts_dir, step_cfg["deposit_overrides"], species_of)
    # A deposit with no mutation rows has an empty kept_expid, so it is not kept.
    is_kept = (deposits["expid"] == deposits["kept_expid"]).fillna(False)
    kept_ids = set(deposits.loc[is_kept, "expid"])
    meta_files = set(deposits.loc[is_kept, "meta_file"])
    muts_files = {f.name for f in muts_dir.glob("*.csv")
                  if int(MUTS_FILE_RE.search(f.name).group(1)) in kept_ids}
    skipped = deposits[~is_kept].sort_values(["exp", "expid"])
    skipped["removal_reason"] = "same experiment name and species as the kept ExpID"
    skipped.loc[skipped["mutation_rows"] == 0, "removal_reason"] = "no mutation rows"
    for r in skipped.itertuples():
        if r.mutation_rows == 0:
            logger.info("skipped deposit ExpID%d (%s): %s", r.expid, r.exp, r.removal_reason)
        else:
            logger.info("skipped deposit ExpID%d (%s, %s): kept ExpID%d", r.expid, r.exp, r.species_norm,
                        r.kept_expid)
    common.save_removed(skipped, common.removed_dir() / "01_skipped_duplicate_deposits.tsv",
                        "same experiment name and species as the kept ExpID")
    # A name is kept once for each of its species, and the later steps tell experiments apart by
    # name, so each name kept for more than one species is listed.
    for exp, picks in deposits[is_kept].sort_values("expid").groupby("exp"):
        if picks["species_norm"].nunique() > 1:
            shown = [f"ExpID{r.expid} ({r.species_norm})" for r in picks.itertuples()]
            logger.info("%s: kept %s and %s", exp, ", ".join(shown[:-1]), shown[-1])
    logger.info("kept %d of %d deposits, one per experiment name and species", len(kept_ids), len(deposits))

    # --- load raw mutations ---
    logger.info("loading raw mutations from %s", muts_dir)
    muts_df, loader_removed, loader_counts = get_all_sample_mut_df(str(muts_dir), logger.info,
                                                                   file_names=muts_files)
    logger.info("loaded %d raw mutation rows from %d files (%d files without mutation rows)",
                len(muts_df), loader_counts["files_loaded"], loader_counts["files_without_rows"])
    logger.info("removed %d sample columns whose names repeat another column, in %d files: %d of them "
                "differ from the column they repeat, with %d calls that are not loaded",
                loader_counts["duplicated_columns_removed"], loader_counts["files_with_duplicated_columns"],
                loader_counts["duplicated_columns_differing"],
                loader_counts["calls_in_differing_duplicated_columns"])
    # Rows the loader cannot use are recorded like every other removal.
    removed_while_loading = {}
    for reason, fname in (("no gene name", "01_removed_rows_without_gene.tsv"),
                          ("Details is Duplication", "01_removed_duplication_rows.tsv")):
        part = loader_removed[loader_removed["removal_reason"] == reason]
        common.save_removed(part, common.removed_dir() / fname, reason)
        removed_while_loading[reason] = int(len(part))
        logger.info("removed %d rows while loading: %s (removed/%s)", len(part), reason, fname)
    # The deposit of each row: the ExpID in the name of the file it was loaded from.
    muts_df["_expid"] = muts_df.pop("mut_file").str.extract(MUTS_FILE_RE)[0].astype(int)

    muts_df["sample"] = muts_df.apply(
        lambda r: get_exp_sample_name(r.exp, r.ale, r.flask, r.isolate, r.tech_rep), axis=1)
    muts_df["sample type"] = muts_df.apply(
        lambda r: "population" if r["isolate"] == 0 else "clone", axis=1)

    muts_df.Position = muts_df.Position.apply(lambda x: int(str(x).replace(",", "")))
    muts_df.Position = muts_df.Position.astype(int)
    muts_df = muts_df.drop(columns=["Mut ID"])

    n_raw_rows = len(muts_df)
    # Samples and ALEs are counted within species, since one experiment name can be kept for two
    # species.
    raw_species = muts_df["_expid"].map(species_of)
    n_raw_samples = muts_df.groupby([raw_species, "sample"], dropna=False).ngroups
    n_raw_ales = muts_df.groupby([raw_species, "exp", "ale"], dropna=False).ngroups
    n_raw_genes = muts_df["Gene (Scrollable)"].nunique()

    # --- merge metadata ---
    combined_metadata = concat_csv_files(meta_dir, meta_files)
    combined_metadata["exp"] = combined_metadata["experiment"].str.replace(" ", "_")
    combined_metadata["_expid"] = combined_metadata["source_file"].map(lambda f: expid_of(f, meta_dir))
    parts = combined_metadata["sample_name"].str.extract(r"A(\d+)\s*F(\d+)\s*I(\d+)\s*R(\d+)")
    combined_metadata[["ale", "flask", "isolate", "tech_rep"]] = parts.astype(int)
    combined_metadata["sample"] = combined_metadata.apply(
        lambda r: get_exp_sample_name(r.exp, r.ale, r.flask, r.isolate, r.tech_rep), axis=1)
    # Deposits can share a name, so a row takes the metadata of its sample in its own deposit.
    combined_metadata = combined_metadata.drop_duplicates(subset=["_expid", "sample"], keep="first")

    redundant_cols = ["exp", "ale", "flask", "isolate", "tech_rep"]
    muts_metadata_df = muts_df.merge(
        combined_metadata.drop(columns=redundant_cols), on=["_expid", "sample"], how="left", indicator=True)
    n_unmatched_meta = int((muts_metadata_df["_merge"] != "both").sum())
    # Kept until the checks below, which stop on a row of the configured species without metadata.
    muts_metadata_df["_has_metadata"] = muts_metadata_df["_merge"] == "both"
    muts_metadata_df = muts_metadata_df.drop(columns=["_merge"])

    # --- merge species normalization: each row takes the species of its deposit ---
    muts_metadata_df = muts_metadata_df.merge(deposit_species, on=["_expid"], how="left")

    # Spaces around a strain name would make one strain look like two.
    muts_metadata_df["strain_details"] = muts_metadata_df["strain_details"].map(
        lambda v: v.strip() if isinstance(v, str) else v)

    # --- annotate subspecies / gene-deletion (NO species filtering here: selection is step 02) ---
    annotated = muts_metadata_df
    if subspecies_cfg:
        subspecies_path = common.resolve(subspecies_cfg)
        subspecies = common.read_table(subspecies_path)
        # One subspecies per (exp, source_file). With two, every matching mutation row would be copied.
        sub_rows = subspecies[["exp", "source_file", "sub_species", "Gene_deletion_record"]].drop_duplicates()
        conflicts = sub_rows[sub_rows.duplicated(subset=["exp", "source_file"], keep=False)]
        if len(conflicts):
            n_keys = conflicts.groupby(["exp", "source_file"], dropna=False).ngroups
            shown = conflicts.sort_values(["exp", "source_file"]).assign(
                sub_species=lambda d: d["sub_species"].map(repr),
                Gene_deletion_record=lambda d: d["Gene_deletion_record"].map(repr))
            raise ValueError(f"{subspecies_path} gives {n_keys} key(s) (exp, source_file) more than one "
                             "value (values in quotes, so spaces show). Keep one row per key:\n"
                             + shown.to_string(index=False))
        annotated = annotated.drop(columns=["_merge", "sub_species"], errors="ignore")
        annotated = annotated.merge(
            subspecies[["exp", "source_file", "sub_species", "Gene_deletion_record"]].drop_duplicates(),
            on=["exp", "source_file"], how="left", indicator=True)
        annotated = annotated.drop(columns=["_merge"])
        logger.info("annotated subspecies from %s", subspecies_cfg)
    else:
        logger.info("no subspecies_table configured, keeping species-level annotation only")

    # Every row has a species, since select_deposits stops on a deposit with mutation rows and no
    # species. Rows of the configured species also need their metadata row and a subspecies. Rows of
    # other species are removed at step 02, so their gaps do not stop this run.
    species = org["species"]
    sub_source = subspecies_path if subspecies_cfg else species_path

    def listed(mask):
        counts = annotated.loc[mask, "exp"].value_counts().sort_index()
        shown = ", ".join(f"{e} ({n} rows)" for e, n in counts.head(20).items())
        more = f", and {len(counts) - 20} more" if len(counts) > 20 else ""
        return f"{int(mask.sum())} rows in {len(counts)} experiment(s): {shown}{more}"

    problems = []
    in_run = annotated["species_norm"] == species
    no_meta = in_run & ~annotated["_has_metadata"]
    if no_meta.any():
        problems.append(f"{species} rows whose sample has no row in its metadata file in {meta_dir}, "
                        f"{listed(no_meta)}")
    # With a subspecies table, the subspecies is looked up by the source file that the metadata row
    # gives, so a row without metadata is reported once, as missing metadata, and not again here.
    needs_meta = annotated["_has_metadata"] if subspecies_cfg else True
    sub = annotated["sub_species"]
    no_sub = in_run & needs_meta & (sub.isna() | (sub.astype(str).str.strip() == ""))
    if no_sub.any():
        problems.append(f"{species} rows with no subspecies in {sub_source}, {listed(no_sub)}")
    if problems:
        raise ValueError("lookup checks failed. Add the missing entries and run again:\n"
                         + "\n".join(problems))
    annotated = annotated.drop(columns=["_has_metadata", "_expid"])
    tables = "the species and subspecies tables give" if subspecies_cfg else "the species table gives"
    logger.info("lookup checks passed: every loaded deposit has a species, every %s row has its metadata "
                "row and a subspecies, and %s one value per key", species, tables)

    muts_ecoli = annotated  # all species, annotated. Step 02 does the species/subspecies selection
    common.validate_common(muts_ecoli, logger, required_cols=REQUIRED_COLUMNS)
    common.write_table(muts_ecoli, out_path)
    logger.info("wrote %s with shape %s (all species, annotated)", out_path, muts_ecoli.shape)

    missing_required = [c for c in REQUIRED_COLUMNS if c not in muts_ecoli.columns]
    metrics = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "input_dir": str(muts_dir), "output_file": str(out_path),
        "n_rows": int(len(muts_ecoli)), "n_columns": int(muts_ecoli.shape[1]),
        # Experiments, samples and ALEs are counted within species, as are raw_samples and raw_ales.
        "n_unique_experiments": int(muts_ecoli[["species_norm", "exp"]].drop_duplicates().shape[0]),
        "n_unique_samples": int(muts_ecoli[["species_norm", "sample"]].drop_duplicates().shape[0]),
        "n_unique_ales": int(muts_ecoli[["species_norm", "exp", "ale"]].drop_duplicates().shape[0]),
        "n_unique_genes_raw": int(n_raw_genes),
        "raw_rows_loaded": int(n_raw_rows), "raw_samples": int(n_raw_samples),
        "raw_ales": int(n_raw_ales), "metadata_unmatched_rows": n_unmatched_meta,
        "deposits_total": len(deposits), "deposits_skipped": len(skipped),
        **loader_counts, "rows_removed_while_loading": removed_while_loading,
        "required_columns_present": [c for c in REQUIRED_COLUMNS if c in muts_ecoli.columns],
        "missing_required_columns": missing_required,
    }
    common.write_checkpoint(metrics, common.checkpoints_dir() / "01_load_aledb_mutations_summary.json")
    common.finish_step(STEP_ID, STEP_NAME, logger)


if __name__ == "__main__":
    main()
