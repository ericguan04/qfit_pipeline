#!/bin/bash
# One SLURM array task: picks its ID from the launch's ID list and runs a stage.
# Usage (by submit_slurm.sh): sbatch --array=... slurm_task.sh STAGE RUN_DIR IDS_FILE
set -uo pipefail
STAGE="$1"; RUN_DIR="$2"; IDS_FILE="$3"
ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${IDS_FILE}")
if [ -z "${ID}" ]; then
  echo "No ID on line ${SLURM_ARRAY_TASK_ID} of ${IDS_FILE}" >&2
  exit 1
fi
# submit_slurm.sh passes the repository path, because SLURM runs a copy of
# this script from its spool directory.
exec bash "${QFIT_PIPELINE_REPO}/lib/run_stage.sh" "${RUN_DIR}" "${STAGE}" "${ID}"
