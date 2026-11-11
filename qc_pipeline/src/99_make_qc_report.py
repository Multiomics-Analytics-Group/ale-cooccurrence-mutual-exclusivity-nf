"""
Step 99: Write the QC pipeline report.

Reads the six step tables from the run folder, output/qc/<run_name>/, and writes
qc_pipeline_report.{tsv,json} there, with one row per step: the numbers of rows,
experiments, samples and genes before and after it. Samples are counted on sample_endpoint
(one per ALE endpoint) from step 04 on, and on exp and sample together before that, since
step 03 rewrites sample as ALE, flask, isolate and tech_rep without exp. Genes are counted on
Gene_name_normalized from step 05 on, and on Gene (Scrollable) before that. Step 01 holds
every species, so its experiments and samples are counted within species: an experiment name
kept for two species counts as two experiments. The final table of the QC pipeline,
qc_mutations.tsv, is written by step 06.
"""
import argparse
import importlib.util
from pathlib import Path

import pandas as pd

SRC = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("qc_common", SRC / "00_common.py")
common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(common)

STEP_ID = "99"
STEP_NAME = "make_qc_report"

STEPS = [
    ("01", "load_aledb_mutations", "01_aledb_mutation_df.tsv"),
    ("02", "select_species", "02_species_selected.tsv"),
    ("03", "remove_starting_strain_mutations", "03_starting_strain_removed.tsv"),
    ("04", "remove_hypermutators", "04_hypermutators_removed.tsv"),
    ("05", "normalize_gene_names", "05_gene_names_normalized.tsv"),
    ("06", "remove_specific_experiments", "qc_mutations.tsv"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    common.load_config(args.config)
    logger = common.setup_logging(STEP_ID, STEP_NAME)

    # --- per-step report rows ---
    report = []
    prev = None
    for step_id, name, fname in STEPS:
        out_file = common.out_dir() / fname
        df = common.read_table(out_file)
        m = common.basic_metrics(df)
        # Step 01 holds every species, so its experiments and samples are counted within species.
        if step_id == "01":
            m["n_unique_experiments"] = int(df[["species_norm", "exp"]].drop_duplicates().shape[0])
            m["n_unique_samples"] = int(df[["species_norm", "exp", "sample"]].drop_duplicates().shape[0])
        row = {
            "step_id": step_id, "step_name": name,
            "input_file": prev["output_file"] if prev else "raw ALEdb mutation export",
            "output_file": str(out_file),
            "n_rows_before": prev["n_rows_after"] if prev else None,
            "n_rows_after": m["n_rows"],
            "n_rows_removed": (prev["n_rows_after"] - m["n_rows"]) if prev else None,
            "n_experiments_before": prev["n_experiments_after"] if prev else None,
            "n_experiments_after": m["n_unique_experiments"],
            "n_samples_before": prev["n_samples_after"] if prev else None,
            "n_samples_after": m["n_unique_samples"],
            "n_genes_before": prev["n_genes_after"] if prev else None,
            "n_genes_after": m["n_unique_genes"],
            "status": "ok",
        }
        report.append(row)
        prev = {"output_file": str(out_file), "n_rows_after": m["n_rows"],
                "n_experiments_after": m["n_unique_experiments"],
                "n_samples_after": m["n_unique_samples"], "n_genes_after": m["n_unique_genes"]}

    report_df = pd.DataFrame(report)
    common.write_table(report_df, common.out_dir() / "qc_pipeline_report.tsv")
    common.write_json(report, common.out_dir() / "qc_pipeline_report.json")
    logger.info("wrote qc_pipeline_report.tsv / .json (%d steps)", len(report))

    summary = {
        "step_id": STEP_ID, "step_name": STEP_NAME,
        "report_files": [str(common.out_dir() / "qc_pipeline_report.tsv"),
                         str(common.out_dir() / "qc_pipeline_report.json")],
        "n_steps_reported": len(report),
    }
    common.write_checkpoint(summary, common.checkpoints_dir() / "99_make_qc_report_summary.json")
    logger.info("step 99 (make_qc_report) finished")


if __name__ == "__main__":
    main()
