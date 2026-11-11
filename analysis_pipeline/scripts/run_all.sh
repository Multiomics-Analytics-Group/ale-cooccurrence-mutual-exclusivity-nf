#!/usr/bin/env bash
# Run the pairwise analysis for all experiments with >= min_endpoints samples.
#   bash analysis_pipeline/scripts/run_all.sh --config analysis_pipeline/config/analysis_config.yaml
# --config is required. The script moves to the repository folder itself, and every output
# goes to output/analysis/<run_name>/ there, with run_name taken from the config.
# A finished run (one with run_finished.txt) is never overwritten. A folder left by an
# unfinished run is removed and the run starts again.
# The interpreter must have the discover, numpy, pandas, scipy, statsmodels and pyyaml
# packages (tqdm is optional and shows a progress bar). Activate your env first, or
# override with PYTHON=/path/to/python.
set -euo pipefail

# DISCOVER ships a recent Intel MKL. On macOS and older (non-AVX) Intel CPUs,
# mixing MKL's OpenMP with another OpenMP runtime, or dispatching AVX code on a
# CPU that lacks it, can crash the process with "Segmentation fault: 11". These
# settings make MKL safe on such machines and are harmless elsewhere (the impact
# is negligible for this pipeline's small matrices). Override any by pre-setting.
export KMP_DUPLICATE_LIB_OK="${KMP_DUPLICATE_LIB_OK:-TRUE}"
export MKL_THREADING_LAYER="${MKL_THREADING_LAYER:-SEQUENTIAL}"
export MKL_ENABLE_INSTRUCTIONS="${MKL_ENABLE_INSTRUCTIONS:-SSE4_2}"
# Python's UTF-8 mode: text is read and written as UTF-8 whatever the system locale.
# Override by pre-setting PYTHONUTF8.
export PYTHONUTF8="${PYTHONUTF8:-1}"
# Python would otherwise write compiled __pycache__ files into the code folders.
export PYTHONDONTWRITEBYTECODE=1

PYTHON="${PYTHON:-python3}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SRC="${REPO_ROOT}/analysis_pipeline/src"

# A relative --config or interpreter path is read from the folder the script was called from.
# A bare interpreter name such as python3 is still looked up on PATH.
case "${PYTHON}" in
  /*) ;;
  */*) PYTHON="${PWD}/${PYTHON}";;
esac
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 ]]; then ARGS+=("$1"); shift; continue; fi
      v="$2"; [[ -z "${v}" || "${v}" = /* ]] || v="${PWD}/${v}"
      ARGS+=(--config "${v}"); shift 2;;
    --config=*)
      v="${1#--config=}"; [[ -z "${v}" || "${v}" = /* ]] || v="${PWD}/${v}"
      ARGS+=(--config "${v}"); shift;;
    *) ARGS+=("$1"); shift;;
  esac
done
cd "${REPO_ROOT}"
"${PYTHON}" "${SRC}/run_all_experiments.py" ${ARGS[@]+"${ARGS[@]}"}
