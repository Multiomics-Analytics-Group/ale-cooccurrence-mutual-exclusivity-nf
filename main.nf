/*
 * ALEdb mutation QC and pairwise co-occurrence analysis.
 *
 * Runs each QC step of qc_pipeline/ as one task, then the analysis of analysis_pipeline/ on
 * the final QC table, qc_mutations.tsv. Every file the steps write is published to
 * <outdir>/qc/<QC run_name>/ and <outdir>/analysis/<analysis run_name>/, in the layout the plain
 * runners use under output/, and so is run_config.yaml, the copy of the QC config the QC steps
 * read. Each run_name is read from its config.
 *
 * A run never writes into an existing run folder: when one of the two exists, the run stops
 * before any task starts. The QC run folder also receives run_started.txt, written when the run
 * starts, and run_finished.txt, written only when the whole run has succeeded.
 *
 * Run it from the repository folder, since relative paths are read from the folder Nextflow
 * is started in:
 *   nextflow run main.nf -params-file params/ecoli.yaml -resume
 */

include { QC_01_LOAD; QC_02_SPECIES; QC_03_STARTING_STRAIN; QC_04_HYPERMUTATORS;
          QC_05_GENE_NAMES; QC_06_EXCLUSIONS; QC_99_REPORT } from './modules/qc'
include { ANALYSIS } from './modules/analysis'

// Every parameter is required. A run without one of them stops before any task starts.
params {
    // The QC config, for example qc_pipeline/config/qc_config.yaml.
    qc_config: Path
    // The analysis config. Its input_table must be output/qc/<QC run_name>/qc_mutations.tsv,
    // where its task finds the final QC table.
    analysis_config: Path
    // The folder with the input data. Every QC task sees it as data/, where the input paths of
    // the QC config point.
    data_dir: Path
    // The folder that receives qc/<QC run_name>/ and analysis/<analysis run_name>/.
    outdir: String
}

workflow {
    main:
    checkValues()
    def names = runNames(params.qc_config, params.analysis_config)
    def qc_dir = file("${params.outdir}/qc/${names.qc}")
    def analysis_dir = file("${params.outdir}/analysis/${names.analysis}")
    def existing = [qc_dir, analysis_dir].findAll { dir -> dir.exists() }
    if( existing )
        error "The run folder(s) ${existing.join(' and ')} already exist, and a run never writes " +
              "into an existing run folder. Delete the folder(s), or choose another outdir or run_name"
    def image = imageRecord()
    def code = codeRecord()

    def qc_code = file("${projectDir}/qc_pipeline")
    def analysis_code = file("${projectDir}/analysis_pipeline")
    // The files the image is built from and the settings of the tasks.
    def environment = [file("${projectDir}/docker"), file("${projectDir}/nextflow.config")]
    def qc_state = folderState([qc_code, params.data_dir] + environment, image)
    def analysis_state = folderState([analysis_code] + environment, image)

    // A preview run (-preview) starts no task and writes no record.
    if( !workflow.preview ) {
        qc_dir.mkdirs()
        qc_dir.resolve('run_started.txt').text =
            "started ${utc(workflow.start)}\n" +
            "qc_config ${params.qc_config.toAbsolutePath()}\n" +
            "analysis_config ${params.analysis_config.toAbsolutePath()}\n" +
            "data_dir ${params.data_dir.toAbsolutePath()}\n" +
            "outdir ${file(params.outdir)}\n" +
            "command cd ${workflow.launchDir} && ${workflow.commandLine}\n" +
            "container ${image}\n" +
            "code ${code}\n"
    }

    QC_01_LOAD(names.qc, qc_code, params.data_dir, qc_state, params.qc_config)
    // Every later step reads the copy of the config that step 01 wrote into the run folder.
    def run_config = QC_01_LOAD.out.config
    QC_02_SPECIES(names.qc, qc_code, params.data_dir, qc_state, run_config, QC_01_LOAD.out.table)
    QC_03_STARTING_STRAIN(names.qc, qc_code, params.data_dir, qc_state, run_config,
                          QC_02_SPECIES.out.table)
    QC_04_HYPERMUTATORS(names.qc, qc_code, params.data_dir, qc_state, run_config,
                        QC_03_STARTING_STRAIN.out.table)
    QC_05_GENE_NAMES(names.qc, qc_code, params.data_dir, qc_state, run_config,
                     QC_04_HYPERMUTATORS.out.table)
    QC_06_EXCLUSIONS(names.qc, qc_code, params.data_dir, qc_state, run_config,
                     QC_05_GENE_NAMES.out.table, QC_01_LOAD.out.table, QC_01_LOAD.out.removed)
    def tables = QC_01_LOAD.out.table.mix(QC_02_SPECIES.out.table, QC_03_STARTING_STRAIN.out.table,
                                          QC_04_HYPERMUTATORS.out.table, QC_05_GENE_NAMES.out.table,
                                          QC_06_EXCLUSIONS.out.table).collect()
    QC_99_REPORT(names.qc, qc_code, params.data_dir, qc_state, run_config, tables)

    ANALYSIS(names.analysis, names.qc, analysis_code, analysis_state, params.analysis_config,
             QC_06_EXCLUSIONS.out.table)

    onComplete:
    // Nextflow runs this after every output file has been published.
    if( workflow.success && !workflow.preview )
        qc_dir.resolve('run_finished.txt').text =
            "finished ${utc(workflow.complete)}\nconfig ${params.qc_config.toAbsolutePath()}\n"
}

// Nextflow stops a run that lacks a parameter, and a run whose qc_config, analysis_config or
// data_dir is blank text or names no existing file. This check stops the other runs without a
// value: any of the four parameters written without a value in a params file arrives as null,
// and outdir, which is text, can arrive as blank text, or as the text true when it is given on
// the command line without a value.
def checkValues() {
    def empty = ['qc_config', 'analysis_config', 'data_dir', 'outdir'].findAll { name ->
        params[name] == null || params[name].toString().trim() in ['', 'true']
    }
    if( empty )
        error "No value for the parameter(s) ${empty.join(', ')}"
}

// The run_name of each config. The QC steps write to output/qc/<QC run_name>/, and the analysis
// reads its input_table there, so input_table must be the final table of that folder.
def runNames(qc_config, analysis_config) {
    def qc = readYaml(qc_config)
    def analysis = readYaml(analysis_config)
    def qc_name = qc instanceof Map ? qc.run_name : null
    if( !qc_name )
        error "The QC config ${qc_config} has no run_name"
    def analysis_name = analysis instanceof Map ? analysis.run_name : null
    if( !analysis_name )
        error "The analysis config ${analysis_config} has no run_name"
    def table = "output/qc/${qc_name}/qc_mutations.tsv".toString()
    if( analysis.input_table != table )
        error "The analysis config ${analysis_config} has input_table ${analysis.input_table}. " +
              "Set it to ${table}, the final table of the QC run ${qc_name}"
    return [qc: qc_name, analysis: analysis_name]
}

// A YAML file, read whole. A file that is not valid YAML stops the run with its name.
def readYaml(path) {
    try {
        return new groovy.yaml.YamlSlurper().parse(path)
    }
    catch( Exception e ) {
        error "Cannot read ${path} as YAML: ${e.message}"
    }
}

// The image the tasks run in, by name and image ID, or none when Docker is off.
def imageRecord() {
    if( !workflow.containerEngine )
        return 'none (Docker is off)'
    def id = hostCommand(['docker', 'image', 'inspect', '--format', '{{.Id}}', workflow.container])
    return "${workflow.container} ${id}".toString()
}

// The code version: the git commit of the repository folder, and whether any of its files
// differs from it. Without git, or outside a git checkout, the version is recorded as unknown.
def codeRecord() {
    if( !file("${projectDir}/.git").exists() )
        return "unknown (${projectDir} is not a git checkout)".toString()
    if( !gitInstalled() )
        return 'unknown (git is not installed)'
    def git = ['git', '-C', projectDir.toString()]
    def commit = hostCommand(git + ['rev-parse', 'HEAD'])
    if( hostCommand(git + ['status', '--porcelain']) != '' )
        return "${commit} with uncommitted changes".toString()
    return commit
}

// Whether git can be started on the host.
def gitInstalled() {
    try {
        ['git', '--version'].execute().waitFor()
        return true
    }
    catch( IOException _e ) {
        return false
    }
}

// Runs a command on the host and returns what it printed, without the final line break. A
// command that cannot run or that fails stops the run with its message.
def hostCommand(command) {
    def out = new StringBuilder()
    def err = new StringBuilder()
    def proc = command.execute()
    proc.waitForProcessOutput(out, err)
    if( proc.exitValue() != 0 )
        error "${command.join(' ')} failed: ${err.toString().trim()}"
    return out.toString().trim()
}

// A checksum of the path, size and date of every file in the folders and of each other file,
// and of the image record. Nextflow checks only the size and date of a folder input itself,
// which an edit to a file inside it leaves unchanged, and it knows the image by its name only.
// The tasks take this checksum as an input, so after such an edit or a new image -resume runs
// them again.
def folderState(paths, image) {
    def lines = [image]
    paths.each { p ->
        if( p.isDirectory() )
            p.eachFileRecurse(groovy.io.FileType.FILES) { f -> lines << "${f} ${f.size()} ${f.lastModified()}".toString() }
        else
            lines << "${p} ${p.size()} ${p.lastModified()}".toString()
    }
    return lines.sort().join('\n').md5()
}

// A time in UTC, to the second, as in 2026-01-15T09:30:12Z.
def utc(time) {
    return time.toInstant().truncatedTo(java.time.temporal.ChronoUnit.SECONDS).toString()
}
