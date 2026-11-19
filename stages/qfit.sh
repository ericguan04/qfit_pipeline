#!/bin/bash
# Pipeline stages 3 and 4: qFit, then removal of exact-copy conformers.
# Usage: qfit.sh CONFIG ID NPROC   (run inside the structure directory)
# Outputs: multiconformer_model2.qfit_raw.pdb   qFit's own output, unchanged
#          multiconformer_model2.pdb            the same model with exact copies
#                                               merged; input to final refinement

set -uo pipefail
CONFIG="$1"; ID="$2"; NPROC="$3"
# shellcheck source=/dev/null
source "${CONFIG}"
source "$(dirname "$0")/common.sh"
id=$(echo "${ID}" | tr '[:upper:]' '[:lower:]')

refined="${id}.updated_refine_001.pdb"
require_files composite_omit_map.mtz "${refined}"
rm -f multiconformer_model2.pdb multiconformer_model2.qfit_raw.pdb

# Stage 3. --scale-rmask is always passed, so the log records the value used.
run qfit_protein composite_omit_map.mtz -l 2FOFCWT,PH2FOFCWT "${refined}" -p "${NPROC}" \
  --scale-rmask "${QFIT_SCALE_RMASK:-1.0}"
require_files multiconformer_model2.pdb
run mv multiconformer_model2.pdb multiconformer_model2.qfit_raw.pdb

# Stage 4. qFit's segment step pads residues with exact copies of their
# conformers, which drift apart in refinement. Merge them (occupancies added)
# before final refinement. The model's calculated structure factors are
# unchanged. The cleaned model keeps qFit's file name because
# qfit_final_refine_xray.sh relies on it.
run phenix.python "$(dirname "$0")/remove_copies.py" \
  multiconformer_model2.qfit_raw.pdb multiconformer_model2.pdb
require_files multiconformer_model2.pdb
