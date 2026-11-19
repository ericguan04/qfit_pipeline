#!/bin/bash
# Runs one stage for one structure, on the host (a SLURM task or run_local.sh).
# Usage: run_stage.sh RUN_DIR STAGE ID
#
# State lives in RUN_DIR/structures/<ID>/status/:
#   <stage>.done     written only when the stage exited 0 and its outputs exist
#   <stage>.failed   exit code, reason and last log lines
#   <stage>.running  host, pid and SLURM job while the stage runs
# A stage with a .done marker is skipped. Running a stage removes the .done
# markers of later stages, since their inputs are about to change.

set -uo pipefail
RUN_DIR="$1"; STAGE="$2"; ID="$3"
REPO=$(cd "$(dirname "$0")/.." && pwd)
CONFIG="${RUN_DIR}/run_info/config.sh"
# shellcheck source=/dev/null
source "${CONFIG}"

STAGES="pre_qfit qfit post_qfit"
case "${STAGE}" in
  pre_qfit)  prev="";         nproc="${PRE_CPUS}";  later="qfit post_qfit" ;;
  qfit)      prev="pre_qfit"; nproc="${QFIT_CPUS}"; later="post_qfit" ;;
  post_qfit) prev="qfit";     nproc="${POST_CPUS}"; later="" ;;
  *) echo "Unknown stage '${STAGE}' (expected one of: ${STAGES})" >&2; exit 2 ;;
esac
nproc="${SLURM_CPUS_PER_TASK:-${nproc}}"

dir="${RUN_DIR}/structures/${ID}"
status="${dir}/status"
log="${dir}/logs/${STAGE}.log"
mkdir -p "${status}" "${dir}/logs/previous"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }
job_ref() {
  if [ -n "${SLURM_ARRAY_JOB_ID:-}" ]; then
    echo "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  else
    echo "${SLURM_JOB_ID:-}"
  fi
}

if [ -f "${status}/${STAGE}.done" ]; then
  echo "[$(stamp)] ${ID} ${STAGE}: already done, skipping."
  exit 0
fi

fail() {  # fail EXIT_CODE REASON
  {
    echo "time=$(stamp)"
    echo "exit_code=$1"
    echo "reason=$2"
    echo "log=${log}"
    echo "job=$(job_ref)"
    echo "--- last lines of log ---"
    [ -f "${log}" ] && tail -n 15 "${log}"
  } > "${status}/${STAGE}.failed"
  rm -f "${status}/${STAGE}.running"
  echo "[$(stamp)] ${ID} ${STAGE}: FAILED ($2). Log: ${log}" >&2
  exit "$1"
}

if [ -n "${prev}" ] && [ ! -f "${status}/${prev}.done" ]; then
  fail 1 "previous stage ${prev} is not done"
fi

# Start this attempt: clear old state, keep the previous log.
rm -f "${status}/${STAGE}.failed"
for s in ${later}; do rm -f "${status}/${s}.done"; done
if [ -f "${log}" ]; then
  mv "${log}" "${dir}/logs/previous/${STAGE}.$(date -r "${log}" '+%Y%m%d_%H%M%S' 2>/dev/null || date '+%Y%m%d_%H%M%S').log"
fi
printf 'host=%s\npid=%s\njob=%s\nstart=%s\n' "$(hostname)" "$$" "$(job_ref)" "$(stamp)" \
  > "${status}/${STAGE}.running"

# SLURM sends SIGTERM at the time limit and on scancel.
trap 'fail 143 "terminated by a signal (time limit or scancel?)"' TERM INT

cmd=(bash "${REPO}/stages/${STAGE}.sh" "${CONFIG}" "${ID}" "${nproc}")
if [ -n "${CONTAINER}" ]; then
  cmd=(apptainer exec --bind "${APPTAINER_BINDS}" "${CONTAINER}" "${cmd[@]}")
fi

start=$(date +%s)
echo "[$(stamp)] ${ID} ${STAGE}: started on $(hostname) with ${nproc} cores. Log: ${log}"
{
  echo "# ${STAGE} for ${ID}, started $(stamp) on $(hostname), job $(job_ref)"
  echo "# command: ${cmd[*]}"
} > "${log}"
# Foreground on purpose: SLURM signals every process of the job and Ctrl-C
# reaches the whole foreground group, so the tools stop and the trap above then
# records the failure. (Background children of a script would ignore Ctrl-C.)
( cd "${dir}" && "${cmd[@]}" ) >> "${log}" 2>&1
rc=$?
elapsed=$(( $(date +%s) - start ))

if [ ${rc} -ne 0 ]; then
  reason="exit code ${rc}"
  [ ${rc} -eq 137 ] && reason="killed (exit 137): out of memory?"
  fail "${rc}" "${reason}"
fi

{
  echo "time=$(stamp)"
  echo "seconds=${elapsed}"
  echo "host=$(hostname)"
  echo "job=$(job_ref)"
  echo "cores=${nproc}"
} > "${status}/${STAGE}.done"
rm -f "${status}/${STAGE}.running"
echo "[$(stamp)] ${ID} ${STAGE}: done in ${elapsed}s."
