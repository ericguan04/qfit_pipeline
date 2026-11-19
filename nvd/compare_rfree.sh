#!/bin/bash
# NVD analysis helper (nvd-experiments branch only; not part of the pipeline).
#
# Tabulates R-free and free-set size per structure from the refinement logs:
#   deposited  start R-free of the pre-qFit refinement (prepared, unrefined model)
#   pre        final R-free of the pre-qFit refinement
#   final      final R-free of qfit_final_refine_xray.sh
# All three use the same free set when the % free columns agree.
#
# Usage: compare_rfree.sh RUN_DIR [OLD_DIR]
#   OLD_DIR: earlier results laid out as OLD_DIR/<ID>/ (e.g. test_data/), added
#            as old_free% / old_final columns.
set -uo pipefail
RUN_DIR="${1:?Usage: $0 RUN_DIR [OLD_DIR]}"
OLD_DIR="${2:-}"

rfree() {  # rfree LOG Start|Final
  [ -f "$1" ] && grep "$2 R-work" "$1" | tail -n 1 | awk -F'R-free = ' '{print $2}' | awk '{print $1}'
}
pct() {
  [ -f "$1" ] && grep -o '[0-9.]* *% free' "$1" | tail -n 1 | awk '{print $1}'
}

{
  printf 'ID\tpre_free%%\tfinal_free%%\tdeposited\tpre\tfinal\tfinal-pre\tfinal-deposited'
  [ -n "${OLD_DIR}" ] && printf '\told_free%%\told_final'
  printf '\n'
  while IFS= read -r ID; do
    id=$(echo "${ID}" | tr '[:upper:]' '[:lower:]')
    d="${RUN_DIR}/structures/${ID}"
    pre_log="${d}/${id}.updated_refine_001.log"
    fin_log="${d}/${id}_qFit.log"
    dep=$(rfree "${pre_log}" Start); pre=$(rfree "${pre_log}" Final); fin=$(rfree "${fin_log}" Final)
    d1=$(awk -v a="${fin}" -v b="${pre}" 'BEGIN { if (a != "" && b != "") printf "%+.4f", a - b }')
    d2=$(awk -v a="${fin}" -v b="${dep}" 'BEGIN { if (a != "" && b != "") printf "%+.4f", a - b }')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' "${ID}" "$(pct "${pre_log}")" "$(pct "${fin_log}")" \
      "${dep}" "${pre}" "${fin}" "${d1}" "${d2}"
    if [ -n "${OLD_DIR}" ]; then
      o="${OLD_DIR}/${ID}"; [ -d "${o}" ] || o="${OLD_DIR}/${id}"
      printf '\t%s\t%s' "$(pct "${o}/${id}_qFit.log")" "$(rfree "${o}/${id}_qFit.log" Final)"
    fi
    printf '\n'
  done < "${RUN_DIR}/run_info/ids.txt"
} | awk -F'\t' '{ for (i = 1; i <= NF; i++) { v = ($i == "" ? "-" : $i); printf "%-17s", v } print "" }'
