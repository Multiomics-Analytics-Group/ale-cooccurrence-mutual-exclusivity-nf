/*
 * The pairwise analysis of analysis_pipeline/ on the final QC table.
 *
 * run_name is the run_name of the analysis config, and qc_run_name that of the QC config. The
 * analysis reads input_table from the folder it is started in and writes its run folder
 * there, as output/analysis/<run_name>/, with run_started.txt, run_config.yaml, logs/ and
 * run_finished.txt. The task places the QC table in qc_run/ and links output/qc/<qc_run_name> to
 * it, so input_table, output/qc/<qc_run_name>/qc_mutations.tsv, finds it. The run folder is
 * published to <outdir>/analysis/<run_name>/.
 * folder_state is a checksum of the path, size and date of every file in analysis_pipeline/ and
 * docker/ and of nextflow.config, and of the image ID. Nextflow checks only the folders
 * themselves and the image name, so this input is what runs a cached task again after a file
 * inside them changes or the image is rebuilt.
 */

process ANALYSIS {
    tag { run_name }
    publishDir path: { params.outdir }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^output\//, '') }

    input:
    val run_name
    val qc_run_name
    path code, stageAs: 'analysis_pipeline'
    val folder_state
    path config, stageAs: 'analysis_config.yaml'
    path qc_table, stageAs: 'qc_run/qc_mutations.tsv'

    output:
    path 'output/analysis/**'

    script:
    """
    mkdir -p output/qc
    ln -s ../../qc_run "output/qc/${qc_run_name}"
    python3 analysis_pipeline/src/run_all_experiments.py --config analysis_config.yaml
    """
}
