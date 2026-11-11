# ALEdb mutation QC and gene co-occurrence analysis

A Nextflow pipeline in two parts, run in one Docker image:

- `qc_pipeline/` cleans the mutation calls of adaptive laboratory evolution (ALE) experiments
  exported from ALEdb into one table of endpoint mutations, `qc_mutations.tsv`.
- `analysis_pipeline/` tests, within each experiment, pairs of mutated genes for
  co-occurrence and mutual exclusivity across the endpoints of the experiment.

## What it does

The QC steps run in order. Each one writes its table and the rows it removed, with the reason.

1. Load the mutation calls and the metadata of the ALEdb deposits in the input data, one deposit
   per experiment name and species, with species and subspecies.
2. Keep the configured species, attach normalized growth conditions and remove the condition
   groups listed in the config.
3. Remove starting-strain mutations: those in the two starting-strain tables and, in each
   experiment, those inherited from its starting-strain sample.
4. Remove hypermutator lineages and high-burden flasks, keep the calls that reach the frequency
   threshold and keep the mutations of each ALE endpoint.
5. Map gene names to standard names with STRING aliases. An intergenic call goes to the nearer
   flanking gene.
6. Remove the experiments listed in the exclusion table. The result is `qc_mutations.tsv`.

A report step then counts rows, experiments, samples and genes before and after each step.

The analysis takes each experiment with at least `min_endpoints` endpoints and builds a binary
matrix of genes by endpoints. Every pair of genes that are each mutated in at least
`fixed_min_count` endpoints is tested for co-occurrence and for mutual exclusivity, with
DISCOVER and with Fisher's exact test. Within each experiment and direction, q-values are
computed over the pairs that can reach significance. The Pearson (phi) correlation of each pair
is reported without a test.

## Requirements

- Nextflow 26.04 or later, which needs Java 17 or later
- Docker, usable by the user who starts Nextflow
- An x86-64 machine, or a Mac with Apple silicon, where Docker Desktop runs the x86-64 image
  through emulation. The machine needs more than 8 GB of total memory: the first QC step asks for
  8 GB, and Nextflow stops a run with a task that asks for more memory than the machine has. The
  second QC step has the highest peak memory use of any task on the E. coli data, about 0.9 GB.
  On a Mac the tasks run in Docker Desktop's virtual machine, so the memory limit set in Docker
  Desktop must also be above that peak.
- git (optional): when the repository folder is a git checkout, the run records the commit of
  the code. Without git the code version is recorded as unknown.
- Free disk: about 1.2 GB for the E. coli run (0.6 GB in `work/` and 0.6 GB in `results/`), about
  0.7 GB for the P. putida run, and about 5 GB in Docker's storage for the image and its build

## Inputs

`data.tar.gz` holds the input data. `tar -xzf data.tar.gz` unpacks it into `data/`:

- `muts_all/` and `metadata_all/`: the mutation calls and the metadata of 410 ALEdb deposits,
  one CSV file each (see Data sources and citations)
- `Species_normalization_latest.tsv`: the species and subspecies of each deposit
- `E_coli_sub_species_check.tsv`: the subspecies and gene deletions of each E. coli deposit that
  step 1 loads
- `Condition_normalization_table.tsv`: normalized growth conditions
- `nonpublic_gene_by_endpoint.tsv`: the mutated genes of each endpoint of four experiments whose
  deposits are not included (see Four experiments that are not public)
- `String/511145.protein.aliases.v12.0.txt.gz`: STRING v12.0 gene aliases of E. coli K-12 MG1655
- `String/160488.protein.aliases.v12.0.txt.gz`: STRING v12.0 gene aliases of P. putida KT2440

`data.md5` lists the md5 checksum of each file in `data/`. To check the unpacked files, run
`md5sum -c --quiet data.md5`: it prints only the files that differ or are missing, and nothing
when all match. It does not report files that it does not list, so unpack `data.tar.gz` where no
earlier `data/` folder exists.

The E. coli settings are in `qc_pipeline/config/qc_config.yaml`, with its tables in
`qc_pipeline/config/` (starting-strain mutations, hypermutator genes, excluded experiments), and
in `analysis_pipeline/config/analysis_config.yaml`. `params/ecoli.yaml` names the two configs, the
data folder and the output folder. The P. putida run reads
`qc_pipeline/config/qc_config_pputida.yaml` and
`analysis_pipeline/config/analysis_config_pputida.yaml`, named in `params/pputida.yaml`. Its
thresholds are those of E. coli. Only the settings that name species-specific things differ: the
species, the subspecies table (none for P. putida), the gene alias file and the alias sources read
from it, the name of the locus-tag column, the deposits to load, the lists of names to remove or
rename, and the starting-strain and exclusion tables (the files ending in `_pputida.tsv`). The
hypermutator gene list is the same file for both. The code has no default values: a missing or
wrong setting or parameter stops the run with a message that names it. Each config has its own
`run_name`, which names its run folder, and the `input_table` of the analysis config must be
`output/qc/<run_name of the QC config>/qc_mutations.tsv`. A QC task finds only `qc_pipeline/` and
the data folder (as `data/`) in its task folder, so every path in the QC config must point into
one of them.

## Build the image

From the repository folder:

```bash
docker build -t ale-cooccurrence:dev docker
```

The image holds the Python environment of both pipelines, with every package pinned, and no
pipeline code or data. `docker/README.md` lists its contents and what the build downloads.

The build downloads 46 conda packages from Anaconda's default channel (repo.anaconda.com) and 3
from conda-forge, and so do the conda commands under Rebuild the environment. Anaconda states that
use of repo.anaconda.com is subject to its Terms of Service
(https://www.anaconda.com/legal/terms/terms-of-service). The terms require a paid plan for use on
behalf of a for-profit organization with more than 200 employees or contractors, including
affiliates, unless the use is one of the free uses they list, which `docker/README.md` summarizes.
Check the terms before you build the image or the environment.

## Run E. coli and P. putida

From the repository folder:

```bash
tar -xzf data.tar.gz
nextflow run main.nf -params-file params/ecoli.yaml -resume
nextflow run main.nf -params-file params/pputida.yaml -resume
```

The two runs are independent and can run in either order, one at a time. They write different run
folders under `results/`. Each QC step and the analysis run as one task in the image, as the user
who started Nextflow. Each task uses one CPU. The QC steps run one after another, and only the
report step and the analysis run at the same time, so more than two CPUs do not shorten a run. The
E. coli run takes about 2 minutes on a Linux server and longer on an Apple M1 Pro
laptop, where Docker runs the x86-64 image through emulation. The P. putida run takes less than a
minute on the same server.

Start Nextflow in the repository folder, since relative paths in the params file and in the
configs are read from the folder it is started in.

Add `-resume` to every run. Nextflow then reuses each task whose inputs have not changed: the
code, the data, the configs, the files in `docker/`, `nextflow.config` and the image. Runs started
with `-resume` in this folder share one cache, so the E. coli and P. putida runs can take turns
and each still reuses its own tasks. Every file in `qc_pipeline/` and `analysis_pipeline/` counts,
so an edit to a config or table of one species also makes tasks of the other species run again. On
the first run there is nothing to reuse, and Nextflow says so in a warning. A run started without
`-resume` starts a new, empty cache, even when it stops before any task, so the next `-resume`
reuses only the tasks run since then.

A run never writes into an existing run folder. When one of its own run folders exists
(`results/qc/ecoli/` or `results/analysis/ecoli/` for E. coli, `results/qc/pputida/` or
`results/analysis/pputida/` for P. putida), the run names the folder and stops before any task
starts. The run folders of the other species do not stop it. To run again, delete the run folders,
or choose another `outdir` in the params file or another `run_name` in the config of the folder (a
new QC `run_name` also needs the matching `input_table` in the analysis config). A new `outdir`
inside the repository folder should be under `results/`, for example `results/run2`. Git ignores
`results/`, so outputs there never count as uncommitted changes in `run_started.txt`. After the
run folders are deleted, `-resume` reuses every task and only copies the outputs again.

Every run writes both run folders, so a run with another analysis config also needs another
`outdir`, or the QC run folder deleted. With `-resume`, it reuses the QC tasks and runs only the
analysis.

A run that stops before the end leaves `run_started.txt` in the QC run folder but no
`run_finished.txt`, next to the outputs of the tasks that finished. Delete its run folders before
the next run, which with `-resume` reuses the tasks that finished. When a step fails, Nextflow
shows what the step printed and names its task folder under `work/`.

## Outputs

Nextflow copies the outputs to `results/` (`outdir` in the params file), into `qc/<run_name>/` and
`analysis/<run_name>/`, with the `run_name` of the QC config and of the analysis config. For the
E. coli run these are the two folders below. A P. putida run writes the same files into
`qc/pputida/` and `analysis/pputida/`. The paths in each table are relative to the run folder in
its heading.

### QC run folder: `results/qc/ecoli/`

| File or folder | Contents |
| --- | --- |
| `qc_mutations.tsv` | The final table: the endpoint mutations that passed every step. The analysis reads it. |
| `qc_pipeline_report.tsv` | The number of rows, experiments, samples and genes before and after each step (also as `.json`) |
| `removed/` | The rows each step removed, in files whose names start with the step number, with the reason for each row in the column `removal_reason`. `01_skipped_duplicate_deposits.tsv` instead lists the ALEdb deposits that step 1 did not load. |
| `qc/` | Tables for review. Step 4 writes the number of mutations of every endpoint, with unusually high numbers flagged, and the mutations found in at least a set fraction of the endpoints of an experiment with enough endpoints (possible starting-strain mutations, see `high_prevalence_report_threshold` and `high_prevalence_min_endpoints` in the QC config). Step 5 writes the gene name mapping and the gene names it could not map, whose rows it removed. |
| `01_aledb_mutation_df.tsv` to `05_gene_names_normalized.tsv` | The table after each of steps 1 to 5 |
| `checkpoints/` | A summary of each step (JSON) |
| `run_config.yaml` | The config the steps read |
| `logs/` | The log of each step and `<step>.console.log`, which holds everything the step printed (stdout and stderr) |
| `run_started.txt` | The start time, the four parameters of the params file, the command, the image by name and ID (none when Docker is off), and the code version (the git commit of the repository folder, with whether any of its files differs from it, or unknown without git or outside a git checkout) |
| `run_finished.txt` | The end time and the path of the QC config. It is written only when the whole run, QC and analysis, has succeeded. |

### Analysis run folder: `results/analysis/ecoli/`

In the file and column names, `unit` means a gene, `sample` means an endpoint, `co` means
co-occurrence, `mutex` means mutual exclusivity and `local` means within one experiment. The
settings named below are in the analysis config, `analysis_pipeline/config/analysis_config.yaml`
(for P. putida `analysis_pipeline/config/analysis_config_pputida.yaml`, with the same values).

| File or folder | Contents |
| --- | --- |
| `all_experiments/pairwise_all_methods_gene.csv` | Every tested gene pair, one row per pair and experiment, with the results of each method. A q-value is empty in a direction in which the pair cannot reach significance. |
| `all_experiments/significant_pairs_gene.csv` | The pairs with a q-value below `significance_threshold` for Fisher's exact test or DISCOVER, in either direction, with the same columns. The column `significant_in` names the tests and directions, such as `Fisher(co)`. |
| `significant_experiments/<experiment>/` | For each experiment with significant pairs: its summary and its significant pairs. A run without significant pairs has no `significant_experiments/` folder. |
| `all_experiments/per_experiment_summary_gene.csv` | One row per analyzed experiment: endpoints, genes, tested genes and pairs, mutated genes per endpoint (mean, median and maximum) and the number of significant pairs of each test and direction. The column `note` says why an experiment has no tested pair. |
| `all_experiments/combined_summary_gene.csv` | The same summary in fewer columns, with `endpoints` for `n_samples` and without `_local` in the names of the significant pair counts, sorted largest first by the number of significant pairs of DISCOVER mutual exclusivity, then DISCOVER co-occurrence, Fisher mutual exclusivity and Fisher co-occurrence |
| `all_experiments/skipped_experiments_gene.csv` | The experiments in `qc_mutations.tsv` with fewer endpoints than `min_endpoints`, which are not analyzed. An experiment that QC removed is in no analysis file: look for its name in `removed/` of the QC run folder. |
| `all_experiments/correlation_matrices/` | One file per analyzed experiment, `<experiment>_gene.csv`: the Pearson (phi) correlation matrix of its tested genes. A gene mutated in every endpoint does not vary, so its row and column are empty. |
| `matrices/gene/<experiment>/` | For each analyzed experiment: the binary matrix of genes by endpoints (`matrix.csv`), the number and fraction of endpoints in which each gene is mutated, and whether it is tested, that is, mutated in at least `fixed_min_count` endpoints (`unit_frequencies.csv`), the number of mutated genes of each endpoint (`sample_burden.csv`), and the counts and settings of the matrix (`meta.json`) |
| `run_config.yaml` | The config the analysis read |
| `logs/` | The log of the analysis |
| `run_started.txt` | The start time, the path of the config, the command and the code version (the git commit, with whether any file in `analysis_pipeline/` differs from it). The image has no git, so a run in the image gives the code version as unknown. |
| `run_finished.txt` | The end time and the path of the config. It is written only when the analysis has succeeded. |

In folder and file names, `<experiment>` is the experiment name with each sequence of characters
other than A to Z, a to z and 0 to 9 replaced by `_`, any `_` at the start or end removed, and the
result cut to 80 characters. `meta.json` gives the full name as `scope`.

### Both run folders

Tasks reused with `-resume` keep the logs and records of the run that made them, so only
`run_finished.txt` in the QC run folder shows that the whole run succeeded. The outputs of a step
that fails are not published. They stay in its task folder, where `qc_run/logs/<step>.console.log`
holds everything the step printed.

Each task runs in its own folder under `work/`. The logs, `checkpoints/`, the QC report and the
analysis `run_started.txt` give paths in those folders. `work/` can be deleted after a run, and
`-resume` then has nothing to reuse.

Results can differ slightly between computers, in the later digits of some numbers.

## Rebuild the environment

The image is built from pinned lists only: `docker/conda-explicit.txt` gives each conda package
by URL and md5, `docker/requirements.txt` each PyPI wheel by sha256, and discover comes from
`docker/vendor/`. The two base images are pinned by digest. To build the image again from
scratch:

```bash
docker build --no-cache -t ale-cooccurrence:dev docker
```

Nextflow finds the image by its name. The run records its image ID, and `-resume` runs every task
again when the image ID changes.

The same environment can be made without Docker, with conda and pip, on Linux x86-64 only:

```bash
sed "s|file:///tmp/build/vendor/|file://$PWD/docker/vendor/|" docker/conda-explicit.txt \
  > /tmp/ale-env.txt
conda create --yes --prefix "$HOME/envs/ale-cooccurrence" --file /tmp/ale-env.txt
"$HOME/envs/ale-cooccurrence/bin/python" -m pip install --no-deps --only-binary=:all: \
  --require-hashes --requirement docker/requirements.txt
```

To run the pipeline in it, put it first on `PATH` and turn Docker off:

```bash
PATH="$HOME/envs/ale-cooccurrence/bin:$PATH" \
  nextflow run main.nf -params-file params/ecoli.yaml -without-docker -resume
```

## Without Nextflow

`qc_pipeline/scripts/run_qc_pipeline.sh` and `analysis_pipeline/scripts/run_all.sh` run the same
steps in one Python environment and write to `output/qc/<run_name>/` and
`output/analysis/<run_name>/`. `docker/README.md` shows how to run them in the image.

## Four experiments that are not public

The article that describes this pipeline also analyzed four E. coli experiments whose ALEdb
deposits are not public: STE_BOP (ExpID 215), GLV (ExpID 68), stALE (ExpID 1024) and
MS56_ALE_2022 (ExpID 2170). Their mutation and metadata CSV files are not in `data.tar.gz`, so
the runs above leave these experiments out. For every other experiment they give the article's
results. Numbers that QC computes from the mutation calls of all experiments differ from those in
the article, such as the row counts of the QC report and the hypermutator cutoff of step 4
(`mutation_burden_upper_whisker` in `checkpoints/04_hypermutator_summary.json`).
`Species_normalization_latest.tsv`, `E_coli_sub_species_check.tsv` and
`Condition_normalization_table.tsv` keep the rows of the four deposits, so a run with their CSV
files added to `data/muts_all/` and `data/metadata_all/` repeats the article's run. The article's
data availability statement says how to request these files.

`data/nonpublic_gene_by_endpoint.tsv` lists the genes that are mutated in each endpoint of the
four experiments, one row per gene and endpoint, as the QC steps of this pipeline left them in the
article's run. The analysis uses only which genes are mutated in which endpoints, so these commands
give the article's analysis results for the four experiments, in
`output/analysis/ecoli_nonpublic/`. Run them from the repository folder, with the image built (see
Build the image) and `data.tar.gz` unpacked:

```bash
mkdir -p output
sed -e 's|^input_table:.*|input_table: data/nonpublic_gene_by_endpoint.tsv|' \
  -e 's|^run_name:.*|run_name: ecoli_nonpublic|' \
  analysis_pipeline/config/analysis_config.yaml > output/analysis_config_nonpublic.yaml
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/analysis_pipeline:/work/analysis_pipeline:ro" \
  -v "$PWD/data:/work/data:ro" \
  -v "$PWD/output:/work/output" \
  -w /work \
  -e KMP_DUPLICATE_LIB_OK=TRUE -e MKL_THREADING_LAYER=GNU -e OMP_NUM_THREADS=1 \
  ale-cooccurrence:dev \
  bash analysis_pipeline/scripts/run_all.sh --config output/analysis_config_nonpublic.yaml
```

Without Docker, in the conda environment made under Rebuild the environment, the last command is:

```bash
OMP_NUM_THREADS=1 MKL_THREADING_LAYER=GNU PATH="$HOME/envs/ale-cooccurrence/bin:$PATH" \
  bash analysis_pipeline/scripts/run_all.sh --config output/analysis_config_nonpublic.yaml
```

In this run, `n_input_rows` of the experiment summaries counts the rows of that table, one per
gene and endpoint, not mutations.

## Data sources and citations

The mutation calls and their metadata come from ALEdb (https://aledb.org). Please cite: Phaneuf
PV, Gosting D, Palsson BO, Feist AM. ALEdb 1.0: a database of mutations from adaptive laboratory
evolution experimentation. Nucleic Acids Research 47(D1), D1164-D1171 (2019).
https://doi.org/10.1093/nar/gky983

The 410 ALEdb deposits in `data.tar.gz` are all public in ALEdb. The mutation file of GYD
(ExpID 682) is taken from a previous publication of these data and differs from the current
ALEdb file.

The gene aliases are the unmodified STRING v12.0 alias files of E. coli K-12 MG1655,
`511145.protein.aliases.v12.0.txt.gz`, and of P. putida KT2440,
`160488.protein.aliases.v12.0.txt.gz` (https://string-db.org), under the Creative Commons
Attribution 4.0 license (CC BY 4.0). Please cite: Szklarczyk D, Kirsch R, Koutrouli M, et al.
The STRING database in 2023: protein-protein association networks and functional enrichment
analyses for any sequenced genome of interest. Nucleic Acids Research 51(D1), D638-D646 (2023).
https://doi.org/10.1093/nar/gkac1000

The included code of mutil is under the MIT License, and discover under the Apache License 2.0
(see License).

## License

This project is licensed under the [MIT License](LICENSE).

The included code keeps its own license:

- `qc_pipeline/mutil/alemutdf.py` and `qc_pipeline/mutil/mut.py` are adapted from mutil
  (https://github.com/Aletechdev/mutil), MIT License, Copyright (c) 2021 Patrick Phaneuf. See
  `qc_pipeline/mutil/LICENSE`.
- `docker/vendor/discover-0.9.6-py312h3866932_0.conda` is the unmodified conda package of
  DISCOVER 0.9.6 by Sander Canisius, Netherlands Cancer Institute
  (http://ccb.nki.nl/software/discover), Apache License 2.0. See `docker/vendor/SOURCE.txt` and
  `docker/vendor/LICENSE-Apache-2.0.txt`, which the image also holds in
  `/opt/env/share/licenses/discover/`.

The STRING alias files in `data.tar.gz` are under the Creative Commons Attribution 4.0 license
(CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), as given under Data sources and
citations.
