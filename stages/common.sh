# Helpers shared by the stage scripts. Sourced, not executed.

die() {
  echo "ERROR: $*" >&2
  exit 1
}

# Fails unless every argument is a non-empty file.
require_files() {
  local f
  for f in "$@"; do
    [ -s "${f}" ] || die "expected output ${f} is missing or empty."
  done
}

# Runs a command, echoing it first, and exits with the command's own exit code
# on failure (137 = killed, e.g. out of memory, which run_stage.sh reports).
run() {
  echo "+ $*"
  "$@"
  local rc=$?
  if [ ${rc} -ne 0 ]; then
    echo "ERROR: command failed with exit code ${rc}: $1" >&2
    exit ${rc}
  fi
}

# Sets XRAY_LABEL and RFREE_LABEL from the column table of an .mtz file.
# The data-label logic is the one in the original pre_qfit.sh and in
# qfit_final_refine_xray.sh (first matching column in file order). The R-free
# label is matched on the whole column name, as in the fixed
# qfit_final_refine_xray.sh, so both refinements use the same column
# (R-free-flags must not match R-free-flags-1).
detect_labels() {
  local meta field ampfields
  local obstypes=("FP" "FOBS" "F-obs" "I" "IOBS" "I-obs" "F(+)" "I(+)" "F-obs-filtered(+)")
  XRAY_LABEL=""
  RFREE_LABEL=""
  meta=$(phenix.mtz.dump "$1") || return 1
  ampfields=$(grep -E "amplitude|intensity|F\(\+\)|I\(\+\)" <<< "${meta}" | awk '{$1=$1};1' | cut -d " " -f 1)
  for field in ${ampfields}; do
    if [[ " ${obstypes[*]} " =~ " ${field} " ]] && grep -F -q -w "SIG${field}" <<< "${meta}"; then
      XRAY_LABEL="${field},SIG${field}"
      break
    fi
  done
  for field in FREE R-free-flags FreeR_flag; do
    if awk -v label="${field}" '$1 == label {found=1} END {exit !found}' <<< "${meta}"; then
      RFREE_LABEL="${field}"
      break
    fi
  done
}

# Free-set percentage phenix.refine reports in a log, e.g. "4.82".
free_percent() {
  [ -f "$1" ] || return 0
  grep -o '[0-9.]* *% free' "$1" | tail -n 1 | awk '{print $1}'
}
