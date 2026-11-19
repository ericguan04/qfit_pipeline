#!/bin/bash
# Pipeline stage 5: final refinement with qfit_final_refine_xray.sh from the
# qFit fork (fixes for occupancy constraint groups, the cull loop and R-free
# flags). Input: multiconformer_model2.pdb with exact copies already merged.
# Usage: post_qfit.sh CONFIG ID NPROC   (run inside the structure directory)
# Outputs: <id>_qFit.{pdb,mtz,log}
#
# No R-free filter. The only extra check is that the final refinement used the
# same free set as the pre-qFit refinement (the fork fix for
# xray_data.r_free_flags.generate); a mismatch is an error.

set -uo pipefail
CONFIG="$1"; ID="$2"; NPROC="$3"
# shellcheck source=/dev/null
source "${CONFIG}"
source "$(dirname "$0")/common.sh"
id=$(echo "${ID}" | tr '[:upper:]' '[:lower:]')

require_files "${id}.mtz" multiconformer_model2.pdb
rm -f "${id}_qFit.pdb" "${id}_qFit.mtz" "${id}_qFit.log"

# qfit_final_refine_xray.sh does not stop at failed steps, so its exit code
# alone is not enough: the output checks below catch what it misses.
run qfit_final_refine_xray.sh "${id}.mtz" multiconformer_model2.pdb
require_files "${id}_qFit.pdb" "${id}_qFit.mtz" "${id}_qFit.log"

pre_free=$(free_percent "${id}.updated_refine_001.log")
final_free=$(free_percent "${id}_qFit.log")
echo "free set: pre-qFit ${pre_free:-?} %, final ${final_free:-?} %"
if [ -z "${pre_free}" ] || [ -z "${final_free}" ]; then
  die "could not read the free-set percentage from the refinement logs."
fi
[ "${pre_free}" = "${final_free}" ] || \
  die "final refinement used a different free set (${final_free} %) from the pre-qFit refinement (${pre_free} %)."

# Record R-free against the pre-qFit baseline (same free set, checked above),
# not against the start of the script's last refinement round.
rfree_of() { grep "^Final R-work" "$1" | tail -n 1 | sed -n 's/.*R-free = *\([0-9.]*\).*/\1/p'; }
pre_rfree=$(rfree_of "${id}.updated_refine_001.log")
final_rfree=$(rfree_of "${id}_qFit.log")
[ -n "${pre_rfree}" ] && [ -n "${final_rfree}" ] || \
  die "could not read the final R-free from the refinement logs."
echo "R-free: pre-qFit ${pre_rfree}, final ${final_rfree}, change $(awk -v a="${pre_rfree}" -v b="${final_rfree}" 'BEGIN {printf "%+.4f", b - a}')"
