# Container image

One image holds the Python environment for both pipelines, `qc_pipeline` and `analysis_pipeline`.
The image has no pipeline code and no data: mount them when you run it.

## Build

From the repository folder:

```bash
docker build -t ale-cooccurrence:dev docker
```

The image is built for linux/amd64 (x86-64). With BuildKit, Docker warns that the platform is
fixed (FromPlatformFlagConstDisallowed): this is expected. Add `--no-cache` to build again from
scratch. The build downloads:

- the base images from Docker Hub, both pinned by digest: `mambaorg/micromamba:2.9.0-ubuntu24.04`
  (build stage only) and `ubuntu:24.04`
- the conda packages in `conda-explicit.txt`: 46 from Anaconda's default channel
  (repo.anaconda.com) and 3 from conda-forge (conda.anaconda.org)
- the wheels in `requirements.txt` from PyPI

discover is installed from `vendor/`. A package whose md5 or sha256 does not match stops the build.
Both lists are installed as written, with no dependency solving, so a change to either must keep
it complete. `pip check` at the end of the build stops it when a Python package lacks a dependency.

Anaconda states that use of repo.anaconda.com is subject to its Terms of Service
(https://www.anaconda.com/legal/terms/terms-of-service). The terms dated July 15, 2025 make use
free for an individual's personal, non-commercial purposes, for use on behalf of or in association
with an eligible academic institution or an eligible non-profit or research organization (as
defined in, and subject to the conditions of, Anaconda's Academic Policy and Non-Profit and
Research Policy), and for use on behalf of a for-profit organization with 200 or fewer employees or
contractors, including affiliates. Other use on behalf of a for-profit organization with more than
200 employees or contractors, including affiliates, needs a paid Business plan. Check the current
terms before you build the image.

## Contents

- `/opt/env`, first on `PATH`: Python 3.12.12 with numpy 1.26.4, pandas 2.3.3, scipy 1.16.3,
  statsmodels 0.14.5, PyYAML 6.0.3, matplotlib 3.10.7 and discover 0.9.6, which uses Intel MKL 2023.1.0.
- `/opt/env/share/licenses/discover/`: the license text of discover (`LICENSE-Apache-2.0.txt`)
  and the note on where its package comes from (`SOURCE.txt`).
- The environment variables `PYTHONDONTWRITEBYTECODE=1`, `LANG=C.UTF-8`, `HOME=/tmp` and
  `MPLCONFIGDIR=/tmp/matplotlib` (a writable home folder and matplotlib settings folder for any
  user).
- The user 65534 (nobody) and group 65534 (nogroup), unless another user is given with `--user`.
- No git, so a runner started in the image records the code version as unknown. The Nextflow
  pipeline records the git commit of the repository folder from the host, in `run_started.txt` of
  the QC run folder.

## Files

- `Dockerfile`: micromamba creates `/opt/env` in a build stage, and the final image copies it onto Ubuntu 24.04.
- `conda-explicit.txt`: the conda packages, each by URL and md5.
- `requirements.txt`: the PyPI wheels, each by sha256.
- `vendor/`: the discover conda package, `SOURCE.txt` (where the package comes from, its checksums
  and its license) and the Apache License 2.0 text.

## Run

The Nextflow pipeline (`main.nf`, described in the README of the repository folder) runs each
step in this image. The two runners can also run in it without Nextflow. From the repository
folder, with the input data unpacked in `data/`:

```bash
mkdir -p output
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/qc_pipeline:/work/qc_pipeline:ro" \
  -v "$PWD/analysis_pipeline:/work/analysis_pipeline:ro" \
  -v "$PWD/data:/work/data:ro" \
  -v "$PWD/output:/work/output" \
  -w /work \
  -e KMP_DUPLICATE_LIB_OK=TRUE -e MKL_THREADING_LAYER=GNU -e OMP_NUM_THREADS=1 \
  ale-cooccurrence:dev bash -c '
    bash qc_pipeline/scripts/run_qc_pipeline.sh --config qc_pipeline/config/qc_config.yaml &&
    bash analysis_pipeline/scripts/run_all.sh --config analysis_pipeline/config/analysis_config.yaml'
```

`--user` makes the output files yours (the image's own user, 65534, cannot write to `output`), and
`output` is made first because Docker would make a missing folder as root. The `-e` settings run
MKL, OpenBLAS and numexpr on one thread, with MKL on the GNU OpenMP runtime.
