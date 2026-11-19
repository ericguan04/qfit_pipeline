#!/bin/bash
# Runs the pipeline on the current machine, with the same stage scripts and
# status files as submit_slurm.sh. Structures run one after another, or
# several at once with -p.
set -uo pipefail
REPO=$(cd "$(dirname "$0")" && pwd)

usage() {
  cat >&2 <<EOF
Usage: $0 -r RUN_DIR [-c CONFIG] [-p N] IDS_FILE

  -r RUN_DIR  run directory. New: created with CONFIG. Existing: continues the
              run with its own config, running only unfinished structures.
  -c CONFIG   config for a new run (default: ${REPO}/config.sh)
  -p N        structures to run at the same time (default 1). Each stage
              still uses the cores set in the config.
  IDS_FILE    one structure ID per line
EOF
  exit 1
}

RUN_DIR=""; CONFIG=""; PARALLEL=1
while getopts "r:c:p:h" opt; do
  case "${opt}" in
    r) RUN_DIR="${OPTARG}" ;;
    c) CONFIG="${OPTARG}" ;;
    p) PARALLEL="${OPTARG}" ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))
[ -n "${RUN_DIR}" ] && [ $# -eq 1 ] || usage

# shellcheck source=lib/setup_run.sh
source "${REPO}/lib/setup_run.sh"
setup_run "${RUN_DIR}" "${CONFIG}" "$1"
echo "Run:        ${RUN_DIR}"
echo "Structures: ${N_TODO} to run"
[ "${N_TODO}" -gt 0 ] || { echo "Nothing to run."; exit 0; }

# All stages of one structure; stops at the first failure.
export REPO RUN_DIR
run_structure() {
  local s
  for s in pre_qfit qfit post_qfit; do
    bash "${REPO}/lib/run_stage.sh" "${RUN_DIR}" "${s}" "$1" || return 1
  done
}
export -f run_structure

xargs -P "${PARALLEL}" -I{} bash -c 'run_structure "$1"' _ {} < "${TODO_FILE}" \
  2>&1 | tee "${SUB_DIR}/run_local.log"
echo
bash "${REPO}/status.sh" -r "${RUN_DIR}"
