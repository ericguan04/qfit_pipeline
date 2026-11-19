#!/bin/bash
# Pipeline stages 1 and 2: inputs, space-group standardisation, ready_set,
# pre-qFit refinement, composite omit map.
# Usage: pre_qfit.sh CONFIG ID NPROC   (run inside the structure directory)
#
# Same steps and phenix.refine settings as the original pre_qfit.sh, without
# the R-free filter. Outputs: <id>.updated_refine_001.{pdb,log},
# composite_omit_map.mtz.

set -uo pipefail
CONFIG="$1"; ID="$2"; NPROC="$3"
# shellcheck source=/dev/null
source "${CONFIG}"
source "$(dirname "$0")/common.sh"
id=$(echo "${ID}" | tr '[:upper:]' '[:lower:]')

refined="${id}.updated_refine_001.pdb"
refine_log="${id}.updated_refine_001.log"
# Remove this stage's outputs from an earlier attempt, so they cannot be
# mistaken for new ones.
rm -f "${id}.updated.pdb" "${refined}" "${refine_log}" "${id}.updated_refine_001.mtz" \
      composite_omit_map.mtz

#______________________________ 1. INPUTS ______________________________
# Raw inputs never change, so they are only fetched or copied once.
if [ ! -s "${id}.pdb" ] || [ ! -s "${id}.mtz" ]; then
  if [ -n "${INPUT_DIR:-}" ]; then
    src=""
    for d in "${INPUT_DIR}/${ID}" "${INPUT_DIR}/${id}"; do
      [ -s "${d}/${id}.mtz" ] && src="${d}" && break
    done
    [ -n "${src}" ] || die "no ${id}.mtz under ${INPUT_DIR}/${ID} or ${INPUT_DIR}/${id}."
    for f in "${id}.pdb" "${id}.mtz" "${id}-sf.cif" "${id}.fa"; do
      [ -f "${src}/${f}" ] && run cp "${src}/${f}" .
    done
  else
    run phenix.fetch_pdb "${id}" action=all convert_to_mtz=True
  fi
fi
require_files "${id}.pdb" "${id}.mtz"

# qFit loses non-standard space-group settings (UnitCell.copy rebuilds the
# group from its number), which misplaces symmetry mates in its map scaling.
# Re-express model and data in the standard setting: a change of axes only.
# No-op for structures already in a standard setting; originals are kept as
# <id>.deposited_setting.{pdb,mtz}.
run phenix.python "$(dirname "$0")/standardize_setting.py" "${id}.pdb" "${id}.mtz"
require_files "${id}.pdb" "${id}.mtz"

#______________________________ 2. LABELS ______________________________
detect_labels "${id}.mtz" || die "phenix.mtz.dump failed on ${id}.mtz."
[ -n "${XRAY_LABEL}" ] || die "no X-ray data label (with SIG column) found in ${id}.mtz."
[ -n "${RFREE_LABEL}" ] || die "no R-free flags (FREE, R-free-flags, FreeR_flag) found in ${id}.mtz."
echo "X-ray data label: ${XRAY_LABEL}"
echo "R-free label:     ${RFREE_LABEL}"

#______________________________ 3. READY_SET ______________________________
run phenix.ready_set "${id}.pdb"
# ready_set writes no updated file when there is nothing to change. Copy the
# input so that every later file name is the same either way.
if [ ! -s "${id}.updated.pdb" ]; then
  echo "ready_set wrote no ${id}.updated.pdb; using ${id}.pdb unchanged."
  run cp "${id}.pdb" "${id}.updated.pdb"
fi
ligand_cif=$(ls ./*.ligands.cif 2>/dev/null | head -n 1)

#______________________________ 4. PRE-QFIT REFINEMENT ______________________________
cmd=(
  phenix.refine
  "${id}.updated.pdb"
  "${id}.mtz"
  "refinement.refine.strategy=individual_sites+individual_adp+occupancies"
  "refinement.main.nqh_flips=true"
  "refinement.main.nproc=${NPROC}"
  "hydrogens.refine=riding"
  "refinement.main.number_of_macro_cycles=${PRE_MACRO_CYCLES}"
  "refinement.input.xray_data.labels=${XRAY_LABEL}"
  "refinement.input.xray_data.r_free_flags.label=${RFREE_LABEL}"
  "data_manager.miller_array.labels.name=${XRAY_LABEL}"
  "data_manager.miller_array.labels.name=${RFREE_LABEL}"
  "refinement.output.write_def_file=false"
  "refinement.output.write_eff_file=false"
  "refinement.output.write_geo_file=false"
  "overwrite=true"
)
if [ "${OPTIMIZE_WEIGHTS}" = true ]; then
  cmd+=("refinement.target_weights.optimize_xyz_weight=true"
        "refinement.target_weights.optimize_adp_weight=true")
fi
[ -n "${ligand_cif}" ] && cmd+=("${ligand_cif}")
run "${cmd[@]}"
require_files "${refined}" "${refine_log}"
grep -E "^(Start|Final) R-work" "${refine_log}" | tail -n 2
echo "free set: $(free_percent "${refine_log}") %"

#______________________________ 5. COMPOSITE OMIT MAP ______________________________
run phenix.composite_omit_map "${refined}" "${id}.mtz" \
  "r_free_flags.label=${RFREE_LABEL}" "input.xray_data.labels=${XRAY_LABEL}" \
  "use_experimental_phases=False" "nproc=${NPROC}"
require_files composite_omit_map.mtz
