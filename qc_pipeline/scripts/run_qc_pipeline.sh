#!/usr/bin/env bash
# =====================================================================
# ALEdb data-QC pipeline runner.
#
# Usage (from any folder, the runner moves to the repository folder itself):
#   bash qc_pipeline/scripts/run_qc_pipeline.sh --config qc_pipeline/config/qc_config.yaml
#
# --config is required. Every output of the run goes to output/qc/<run_name>/ in the
# repository folder, with run_name taken from the config.
# A finished run (one with run_finished.txt) is never overwritten. A folder left by an
# unfinished run of this pipeline (one with run_started.txt) is removed and the run starts
# again from step 01. While a run is going, the folder <run_name>.lock beside the run folder
# keeps other runs with the same run_name out.
# The run folder records how the run was made: run_started.txt gives the start time, the command
# and the code version (git commit), run_config.yaml is the config the steps read, and
# logs/run_qc_pipeline.log holds everything the steps showed on the terminal.
#
# The Python interpreter must have pandas, numpy, matplotlib and pyyaml.
# Override it with:  PYTHON=/path/to/python bash qc_pipeline/scripts/run_qc_pipeline.sh ...
# =====================================================================
set -euo pipefail

# Force Python UTF-8 mode so an unset or non-UTF-8 terminal locale cannot crash
# the interpreter with "nl_langinfo(CODESET) failed" or mis-decode input tables.
# Override by pre-setting PYTHONUTF8.
export PYTHONUTF8="${PYTHONUTF8:-1}"
# Python would otherwise write compiled __pycache__ files into the code folders.
export PYTHONDONTWRITEBYTECODE=1

PYTHON="${PYTHON:-python3}"

# Resolve the pipeline root (parent of this scripts/ directory) and the repo root.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${PIPELINE_DIR}/.." && pwd -P)"
SRC="${PIPELINE_DIR}/src"

USAGE="usage: run_qc_pipeline.sh --config <config file>"
CALLED_FROM="${PWD}"
CALLED_AS=("$0" "$@")
CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 || "$2" == --* ]]; then
        echo "${USAGE}" >&2
        echo "$1 needs a value" >&2
        exit 2
      fi
      CONFIG="$2"
      shift 2;;
    --config=*) CONFIG="${1#--config=}"; shift;;
    *) echo "${USAGE}" >&2; echo "Unknown argument: $1" >&2; exit 2;;
  esac
done
if [[ -z "${CONFIG}" ]]; then
  echo "${USAGE}" >&2
  echo "--config is required and must not be empty" >&2
  exit 2
fi
# A relative config or interpreter path is read from the folder the runner was called from.
# A bare interpreter name such as python3 is still looked up on PATH.
case "${CONFIG}" in
  /*) ;;
  *) CONFIG="${PWD}/${CONFIG}";;
esac
case "${PYTHON}" in
  /*) ;;
  */*) PYTHON="${PWD}/${PYTHON}";;
esac

# Every path in the config is read from the repository folder, and the outputs go there.
cd "${REPO_ROOT}"

# Load and check the config before anything is written, and get the run folder from it.
if ! RUN_DIR="$("${PYTHON}" - "${SRC}/00_common.py" "${CONFIG}" <<'PY'
import importlib.util
import sys

spec = importlib.util.spec_from_file_location("qc_common", sys.argv[1])
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)
common.load_config(sys.argv[2])
print(common.out_dir())
PY
)"; then
  echo "Checking the config ${CONFIG} with ${PYTHON} failed (see the message above). Nothing was run." >&2
  exit 1
fi
if ! PYTHON_PATH="$(command -v -- "${PYTHON}")"; then
  echo "Cannot find the interpreter ${PYTHON}. Nothing was run." >&2
  exit 1
fi

# The code version: the git commit, and whether qc_pipeline differs from it. Without git, or
# outside a git checkout, it is recorded as unknown. A checkout that git cannot read stops the run.
code_version() {
  local commit changes
  if ! command -v git > /dev/null 2>&1; then
    echo "unknown (git is not installed)"
    return 0
  fi
  if [[ ! -e "${REPO_ROOT}/.git" ]]; then
    echo "unknown (${REPO_ROOT} is not a git checkout)"
    return 0
  fi
  if ! commit="$(git -C "${REPO_ROOT}" rev-parse HEAD)"; then
    echo "git could not read the commit of ${REPO_ROOT} (see the message above)." >&2
    return 1
  fi
  if ! changes="$(git -C "${REPO_ROOT}" status --porcelain -- qc_pipeline)"; then
    echo "git could not list the changes in ${REPO_ROOT}/qc_pipeline (see the message above)." >&2
    return 1
  fi
  if [[ -n "${changes}" ]]; then
    echo "${commit} with uncommitted changes in qc_pipeline"
  else
    echo "${commit}"
  fi
}
if ! CODE_VERSION="$(code_version)"; then
  echo "Nothing was run." >&2
  exit 1
fi
COMMAND="$(printf 'cd %q && PYTHON=%q bash' "${CALLED_FROM}" "${PYTHON_PATH}"; printf ' %q' "${CALLED_AS[@]}")"

# 1. The run folder and its lock. The lock is a folder beside the run folder, <run folder>.lock,
# made with mkdir, which only one run can do, so two runs never use the same run folder. A
# finished run is never overwritten. A folder left by an unfinished run of this pipeline (it holds
# run_started.txt) is removed and the run starts again. Any other folder that is not empty
# is left alone and the run stops.
LOCK_DIR="${RUN_DIR}.lock"
OWNER="${LOCK_DIR}/owner.txt"
STARTED="${RUN_DIR}/run_started.txt"
FINISHED="${RUN_DIR}/run_finished.txt"
STARTED_AT="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
COPY_PIPE="${LOCK_DIR}/step_output"
STEP_PID=""
COPY_PID=""
PIPELINE_LOG=""
LOCK_NOTE=""

# True when a process with this id exists on this machine, under any user.
is_alive() {
  local out
  [[ "${1:-}" =~ ^[0-9]+$ ]] || return 1
  out="$(kill -0 "$1" 2>&1)" && return 0
  [[ "${out}" == *"not permitted"* ]]
}

# The lock records the runner's process id, the machine, the start time and the running step.
write_owner() {
  printf '%s %s %s\n%s\n' "$$" "$(hostname)" "${STARTED_AT}" "${STEP_PID}" > "${OWNER}.tmp"
  mv -f -- "${OWNER}.tmp" "${OWNER}"
}

mkdir -p "$(dirname "${RUN_DIR}")"
if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
  OWNER_PID=""; OWNER_HOST=""; OWNER_STEP=""
  if [[ -f "${OWNER}" ]]; then
    read -r OWNER_PID OWNER_HOST _ < "${OWNER}" || true
    OWNER_STEP="$(sed -n 2p "${OWNER}")"
  fi
  if [[ -z "${OWNER_HOST}" ]]; then
    echo "${RUN_DIR} is locked by a run that cannot be identified (${LOCK_DIR}). If no run is" \
         "going, delete ${LOCK_DIR} and start again." >&2
    exit 1
  fi
  if [[ "${OWNER_HOST}" != "$(hostname)" ]]; then
    echo "${RUN_DIR} is in use by a run on ${OWNER_HOST}. If no run is going, delete" \
         "${LOCK_DIR} and start again." >&2
    exit 1
  fi
  if is_alive "${OWNER_PID}" || is_alive "${OWNER_STEP}"; then
    WHO="process ${OWNER_PID}"
    [[ -z "${OWNER_STEP}" ]] || WHO="${WHO}, step process ${OWNER_STEP}"
    echo "${RUN_DIR} is in use by a run that is still going (${WHO})." >&2
    exit 1
  fi
  LOCK_NOTE="Removed the lock of a run that no longer exists: ${LOCK_DIR}"
  rm -rf -- "${LOCK_DIR}"
  if ! mkdir "${LOCK_DIR}" 2>/dev/null; then
    echo "${RUN_DIR} was just locked by another run. Nothing was run." >&2
    exit 1
  fi
fi
write_owner

# A run that ends for any reason releases its lock. A signal also stops the step that is running,
# so no step keeps writing into a folder that a later run may remove.
release_lock() { rm -rf -- "${LOCK_DIR}"; }
on_signal() {
  local waited=0
  if [[ -n "${STEP_PID}" ]] && is_alive "${STEP_PID}"; then
    kill -TERM "${STEP_PID}" 2>/dev/null || true
    wait "${STEP_PID}" 2>/dev/null || true
  fi
  # The copy of the step's output ends when the step's end of the pipe closes. It still waits only
  # when the step was stopped before it opened the pipe, so after 2 seconds it is stopped.
  if [[ -n "${COPY_PID}" ]]; then
    while is_alive "${COPY_PID}" && [[ ${waited} -lt 20 ]]; do
      sleep 0.1
      waited=$((waited + 1))
    done
    kill -TERM "${COPY_PID}" 2>/dev/null || true
    wait "${COPY_PID}" 2>/dev/null || true
  fi
  local msg="Stopped by a signal. The run is unfinished: the same command removes it and starts again."
  if [[ -n "${PIPELINE_LOG}" && -f "${PIPELINE_LOG}" ]]; then
    echo "${msg}" | tee -a "${PIPELINE_LOG}" >&2
  else
    echo "${msg}" >&2
  fi
  exit "$1"
}
trap release_lock EXIT
trap 'on_signal 130' INT
trap 'on_signal 143' TERM
trap 'on_signal 129' HUP

# The config is copied before an old run folder is removed, since it may sit inside that folder.
cp -- "${CONFIG}" "${LOCK_DIR}/run_config.yaml"

CLEARED=""
if [[ -d "${RUN_DIR}" ]]; then
  if [[ -f "${FINISHED}" ]]; then
    WHEN="$(head -n 1 "${FINISHED}")"
    echo "${RUN_DIR} holds a finished run (${WHEN:-run_finished.txt is empty}). A finished run is" \
         "never overwritten: use a new run_name, or delete the folder to run it again." >&2
    exit 1
  fi
  if [[ -n "$(ls -A "${RUN_DIR}")" ]]; then
    if [[ ! -f "${STARTED}" ]]; then
      echo "${RUN_DIR} exists but was not made by a run of this pipeline (it has no" \
           "run_started.txt). Move it away, or choose another run_name." >&2
      exit 1
    fi
    CLEARED="${RUN_DIR} held an unfinished run, which was removed before this run."
  fi
  case "${RUN_DIR}" in
    "${REPO_ROOT}"/output/qc/?*) ;;
    *) echo "Not removing ${RUN_DIR}: it is not a run folder under ${REPO_ROOT}/output/qc." >&2
       exit 1;;
  esac
  rm -rf -- "${RUN_DIR}"
fi
mkdir -p "${RUN_DIR}/logs"
printf 'started %s\nconfig %s\ncommand %s\ncode %s\n' "${STARTED_AT}" "${CONFIG}" "${COMMAND}" \
  "${CODE_VERSION}" > "${STARTED}.tmp"
mv -- "${STARTED}.tmp" "${STARTED}"
# The steps read this copy, so it is the config the run used even if the original is edited
# while the run is going.
RUN_CONFIG="${RUN_DIR}/run_config.yaml"
mv -- "${LOCK_DIR}/run_config.yaml" "${RUN_CONFIG}"
COMMON_ARGS=(--config "${RUN_CONFIG}")
PIPELINE_LOG="${RUN_DIR}/logs/run_qc_pipeline.log"
echo "QC pipeline run started: ${STARTED_AT}" | tee "${PIPELINE_LOG}"
for note in "${LOCK_NOTE}" "${CLEARED}"; do
  if [[ -n "${note}" ]]; then
    echo "${note}" | tee -a "${PIPELINE_LOG}"
  fi
done

# Each step runs in the background and the runner waits for it, so a signal to the runner
# stops the step as well. Everything the step shows on the terminal, errors included, goes
# through a named pipe to tee, which shows it and copies it into the run log.
run_step() {
  local status=0 copy_status=0
  rm -f -- "${COPY_PIPE}"
  mkfifo -- "${COPY_PIPE}"
  tee -i -a "${PIPELINE_LOG}" < "${COPY_PIPE}" &
  COPY_PID=$!
  "${PYTHON}" "$@" > "${COPY_PIPE}" 2>&1 &
  STEP_PID=$!
  write_owner
  wait "${STEP_PID}" || status=$?
  STEP_PID=""
  write_owner
  wait "${COPY_PID}" || copy_status=$?
  COPY_PID=""
  rm -f -- "${COPY_PIPE}"
  if [[ "${copy_status}" -ne 0 ]]; then
    echo "Copying the output of $1 into ${PIPELINE_LOG} failed." >&2
    return 1
  fi
  return "${status}"
}

STEPS=(
  "01_load_aledb_mutations.py"
  "02_select_species.py"
  "03_remove_starting_strain_mutations.py"
  "04_remove_hypermutators.py"
  "05_normalize_gene_names.py"
  "06_remove_specific_experiments.py"
)

# 2. Run each QC step in order and stop at the first that fails. Per-step logs and
# checkpoint summaries are written by the step scripts themselves.
for step in "${STEPS[@]}"; do
  echo "----------------------------------------------------------------" | tee -a "${PIPELINE_LOG}"
  echo "Running ${step}" | tee -a "${PIPELINE_LOG}"
  if ! run_step "${SRC}/${step}" "${COMMON_ARGS[@]}"; then
    echo "STEP FAILED: ${step}. Stopping." | tee -a "${PIPELINE_LOG}" >&2
    exit 1
  fi
  echo "Completed ${step}" | tee -a "${PIPELINE_LOG}"
done

# 3. Final QC report.
echo "----------------------------------------------------------------" | tee -a "${PIPELINE_LOG}"
echo "Running 99_make_qc_report.py" | tee -a "${PIPELINE_LOG}"
if ! run_step "${SRC}/99_make_qc_report.py" "${COMMON_ARGS[@]}"; then
  echo "REPORT STEP FAILED. Stopping." | tee -a "${PIPELINE_LOG}" >&2
  exit 1
fi

echo "QC pipeline run finished: $(date -u +'%Y-%m-%dT%H:%M:%SZ')" | tee -a "${PIPELINE_LOG}"
echo "Outputs under: ${RUN_DIR}" | tee -a "${PIPELINE_LOG}"
# Written last, so it exists only when every step finished.
printf 'finished %s\nconfig %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "${CONFIG}" > "${FINISHED}.tmp"
mv -- "${FINISHED}.tmp" "${FINISHED}"
