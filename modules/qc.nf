/*
 * The QC steps of qc_pipeline/, one process each.
 *
 * A step reads every path of its config from the folder it is started in, and finds its run
 * folder there as output/qc/<run_name>/. Each task is started in its task folder, which holds:
 *   qc_pipeline/   the step code, the mutil package and the config tables
 *   data/          the input data
 *   qc_run/        the run folder: run_config.yaml, the config every step reads, and the files
 *                  of the earlier steps that this step reads
 * The task links output/qc/<run_name> to qc_run, because Nextflow places input files only at
 * names fixed in the process. The files the step writes into qc_run/ are the task outputs, and
 * they are published to <outdir>/qc/<run_name>/ in the same layout. Everything the step prints
 * (stdout and stderr) is also kept as logs/<step>.console.log. The task scripts run
 * with pipefail (see nextflow.config), so a step that fails also fails its task.
 * folder_state is a checksum of the path, size and date of every file in qc_pipeline/, data/ and
 * docker/ and of nextflow.config, and of the image ID. Nextflow checks only the folders
 * themselves and the image name, so this input is what runs a cached task again after a file
 * inside them changes or the image is rebuilt.
 */

process QC_01_LOAD {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_config.yaml'

    output:
    path 'qc_run/run_config.yaml', emit: config
    path 'qc_run/01_aledb_mutation_df.tsv', emit: table
    path 'qc_run/removed/*', emit: removed
    path 'qc_run/{checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    cp qc_config.yaml qc_run/run_config.yaml
    python3 qc_pipeline/src/01_load_aledb_mutations.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/01_load_aledb_mutations.console.log
    """
}

process QC_02_SPECIES {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path table_01, stageAs: 'qc_run/*'

    output:
    path 'qc_run/02_species_selected.tsv', emit: table
    path 'qc_run/{removed,checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/02_select_species.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/02_select_species.console.log
    """
}

process QC_03_STARTING_STRAIN {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path table_02, stageAs: 'qc_run/*'

    output:
    path 'qc_run/03_starting_strain_removed.tsv', emit: table
    path 'qc_run/{removed,checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/03_remove_starting_strain_mutations.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/03_remove_starting_strain_mutations.console.log
    """
}

process QC_04_HYPERMUTATORS {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path table_03, stageAs: 'qc_run/*'

    output:
    path 'qc_run/04_hypermutators_removed.tsv', emit: table
    path 'qc_run/{removed,qc,checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/04_remove_hypermutators.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/04_remove_hypermutators.console.log
    """
}

process QC_05_GENE_NAMES {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path table_04, stageAs: 'qc_run/*'

    output:
    path 'qc_run/05_gene_names_normalized.tsv', emit: table
    path 'qc_run/{removed,qc,checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/05_normalize_gene_names.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/05_normalize_gene_names.console.log
    """
}

// When an exclusion rule matches no experiment of its input, step 06 also reads the step 01
// table and removed/01_*.tsv, to tell a rule whose experiments left the run in an earlier step
// from a rule that matches no experiment at all.
process QC_06_EXCLUSIONS {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path table_05, stageAs: 'qc_run/*'
    path table_01, stageAs: 'qc_run/*'
    path removed_01, stageAs: 'qc_run/removed/*'

    output:
    path 'qc_run/qc_mutations.tsv', emit: table
    path 'qc_run/{removed,checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/06_remove_specific_experiments.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/06_remove_specific_experiments.console.log
    """
}

// The report reads the six step tables.
process QC_99_REPORT {
    tag { run_name }
    publishDir path: { "${params.outdir}/qc/${run_name}" }, mode: 'copy', saveAs: { f -> f.replaceFirst(/^qc_run\//, '') }

    input:
    val run_name
    path code, stageAs: 'qc_pipeline'
    path data, stageAs: 'data'
    val folder_state
    path config, stageAs: 'qc_run/run_config.yaml'
    path tables, stageAs: 'qc_run/*'

    output:
    path 'qc_run/qc_pipeline_report.{tsv,json}'
    path 'qc_run/{checkpoints,logs}/*'

    script:
    """
    mkdir -p qc_run/logs output/qc
    ln -s ../../qc_run "output/qc/${run_name}"
    python3 qc_pipeline/src/99_make_qc_report.py --config "output/qc/${run_name}/run_config.yaml" 2>&1 \\
        | tee qc_run/logs/99_make_qc_report.console.log
    """
}
