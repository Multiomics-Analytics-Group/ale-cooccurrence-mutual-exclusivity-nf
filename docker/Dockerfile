# The Python environment that runs both pipelines, qc_pipeline and analysis_pipeline.
# The pipeline code and the input data are not in the image: mount them at run time.
#
# Build from the repository folder:
#   docker build -t ale-cooccurrence:dev docker
#
# Every package is pinned: conda-explicit.txt lists the conda packages by URL and md5 and
# requirements.txt the PyPI wheels by sha256. Nothing is solved at build time, and a download
# that does not match its checksum stops the build. Both stages are built for linux/amd64, the
# platform of the pinned packages.

# Stage 1: micromamba creates the environment in /opt/env.
FROM --platform=linux/amd64 mambaorg/micromamba:2.9.0-ubuntu24.04@sha256:d586ee0cd12ffd035e36ffb28f245b0288de088f44a77a2e4e4d6feb670299e5 AS build
USER root
COPY conda-explicit.txt requirements.txt /tmp/build/
COPY vendor/discover-0.9.6-py312h3866932_0.conda /tmp/build/vendor/
RUN micromamba create --yes --prefix /opt/env --file /tmp/build/conda-explicit.txt \
 && /opt/env/bin/python -m pip install --no-cache-dir --disable-pip-version-check \
      --root-user-action=ignore --no-deps --only-binary=:all: --require-hashes \
      --requirement /tmp/build/requirements.txt \
 && /opt/env/bin/python -m pip check

# Stage 2: the environment on Ubuntu 24.04.
FROM --platform=linux/amd64 ubuntu:24.04@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55
COPY --from=build /opt/env /opt/env
# The license text of discover and the note on where its package comes from.
COPY vendor/LICENSE-Apache-2.0.txt vendor/SOURCE.txt /opt/env/share/licenses/discover/
# HOME and MPLCONFIGDIR: a user given with --user may have no home folder, and matplotlib needs a
# writable settings folder.
ENV PATH=/opt/env/bin:$PATH \
    LANG=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1 \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib
# Stop the build if a package the pipelines import does not load on this base image.
RUN python -c "import numpy, scipy, pandas, statsmodels, yaml, matplotlib.cbook, discover, _discover" \
 && rm -rf /tmp/matplotlib
# The image runs as user and group 65534 (nobody) unless the user is given, for example with
# docker run --user "$(id -u):$(id -g)".
USER 65534:65534
