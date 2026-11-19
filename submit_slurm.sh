#!/bin/bash
# Submits the pipeline to SLURM as three job arrays (pre_qfit -> qfit -> post_qfit).
# Task i of each array handles the i-th structure, and only starts once task i
# of the previous array succeeded (--dependency=aftercorr).
set -uo pipefail
REPO=$(cd "$(dirname "$0")" && pwd)

usage() {
  cat >&2 <<EOF
Usage: $0 -r RUN_DIR [-c CONFIG] [-n] IDS_FILE

  -r RUN_DIR  run directory. New: created with CONFIG. Existing: continues the
              run with its own config, submitting only unfinished structures.
  -c CONFIG   config for a new run (default: ${REPO}/config.sh)
  -n          dry run: set up and print the sbatch commands without submitting
  IDS_FILE    one structure ID per line (blank lines and # comments ignored)

Example:
  $0 -r /data/\$USER/qfit_runs/batch1 -c my_config.sh ids.txt
To retry failed or interrupted structures, run the same command again
(without -c): finished structures and finished stages are skipped.
EOF
  exit 1
}

RUN_DIR=""; CONFIG=""; DRY_RUN=false
while getopts "r:c:nh" opt; do
  case "${opt}" in
    r) RUN_DIR="${OPTARG}" ;;
    c) CONFIG="${OPTARG}" ;;
    n) DRY_RUN=true ;;
    *) usage ;;
  esac
done
shift $((OPTIND - 1))
[ -n "${RUN_DIR}" ] && [ $# -eq 1 ] || usage
${DRY_RUN} || command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found." >&2; exit 1; }

# shellcheck source=lib/setup_run.sh
source "${REPO}/lib/setup_run.sh"
setup_run "${RUN_DIR}" "${CONFIG}" "$1"

echo "Run:        ${RUN_DIR}"
echo "Launch:     ${SUB_DIR}"
echo "Structures: ${N_TODO} to run (of $(wc -l < "${SUB_DIR}/requested_ids.txt" | tr -d ' ') requested)"
if [ "${N_TODO}" -eq 0 ]; then
  echo "Nothing to submit."
  exit 0
fi

common=(--parsable --export=ALL,QFIT_PIPELINE_REPO="${REPO}")
[ -n "${SLURM_ACCOUNT}" ] && common+=(--account="${SLURM_ACCOUNT}")
[ -n "${SLURM_PARTITION}" ] && common+=(--partition="${SLURM_PARTITION}")
# shellcheck disable=SC2206
[ -n "${SLURM_EXTRA_ARGS}" ] && common+=(${SLURM_EXTRA_ARGS})

submit() {  # submit STAGE CPUS MEM TIME LIMIT DEPENDENCY_JOB
  local stage="$1" cpus="$2" mem="$3" time="$4" limit="$5" dep="$6" args
  args=("${common[@]}"
        --job-name="qfit_${stage}"
        --array="1-${N_TODO}%${limit}"
        --cpus-per-task="${cpus}" --mem="${mem}" --time="${time}"
        --output="${SUB_DIR}/slurm/${stage}_%a.out")
  [ -n "${dep}" ] && args+=(--dependency="aftercorr:${dep}" --kill-on-invalid-dep=yes)
  args+=("${REPO}/lib/slurm_task.sh" "${stage}" "${RUN_DIR}" "${TODO_FILE}")
  echo "sbatch ${args[*]}" >> "${SUB_DIR}/sbatch_commands.txt"
  if ${DRY_RUN}; then
    echo "DRYRUN_${stage}"
  else
    local out
    out=$(sbatch "${args[@]}") || return 1
    echo "${out%%;*}"   # --parsable prints "jobid" or "jobid;cluster"
  fi
}

mkdir -p "${SUB_DIR}/slurm"
job_pre=$(submit pre_qfit "${PRE_CPUS}" "${PRE_MEM}" "${PRE_TIME}" "${LIMIT_PRE}" "") \
  || { echo "sbatch failed for pre_qfit." >&2; exit 1; }
job_qfit=$(submit qfit "${QFIT_CPUS}" "${QFIT_MEM}" "${QFIT_TIME}" "${LIMIT_QFIT}" "${job_pre}") \
  || { echo "sbatch failed for qfit; cancel job ${job_pre}." >&2; exit 1; }
job_post=$(submit post_qfit "${POST_CPUS}" "${POST_MEM}" "${POST_TIME}" "${LIMIT_POST}" "${job_qfit}") \
  || { echo "sbatch failed for post_qfit; cancel jobs ${job_pre} ${job_qfit}." >&2; exit 1; }

printf 'pre_qfit=%s\nqfit=%s\npost_qfit=%s\n' "${job_pre}" "${job_qfit}" "${job_post}" > "${SUB_DIR}/jobs.txt"
if ${DRY_RUN}; then
  echo "Dry run. sbatch commands:"
  sed 's/^/  /' "${SUB_DIR}/sbatch_commands.txt"
else
  echo "Submitted: pre_qfit ${job_pre}, qfit ${job_qfit}, post_qfit ${job_post}"
  echo "Status:    ${REPO}/status.sh -r ${RUN_DIR}"
fi
