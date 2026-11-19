#!/bin/bash
# Shows where each structure of a run is: finished, running, or which stage
# failed and why. Reads only the status files, so it works at any time.
set -uo pipefail
REPO=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=lib/state.sh
source "${REPO}/lib/state.sh"

usage() {
  cat >&2 <<EOF
Usage: $0 -r RUN_DIR [-l STATE]

  -r RUN_DIR  run directory
  -l STATE    only print the IDs in STATE (complete, running, failed,
              interrupted, pending), one per line
EOF
  exit 1
}

RUN_DIR=""; LIST=""
while getopts "r:l:h" opt; do
  case "${opt}" in
    r) RUN_DIR="${OPTARG}" ;;
    l) LIST="${OPTARG}" ;;
    *) usage ;;
  esac
done
[ -n "${RUN_DIR}" ] || usage
ids="${RUN_DIR}/run_info/ids.txt"
[ -f "${ids}" ] || { echo "Not a run directory (no run_info/ids.txt): ${RUN_DIR}" >&2; exit 1; }

# One line per structure: ID, overall state, stage, detail. (A function, not
# inline in $(...): bash 3.2 cannot parse `case` inside command substitution.)
build_rows() {
  local id d overall where detail s st
  while IFS= read -r id; do
    d="${RUN_DIR}/structures/${id}"
    overall=complete; where="-"; detail=""
    for s in ${STAGES}; do
      st=$(stage_state "${d}" "${s}")
      [ "${st}" = done ] && continue
      overall="${st}"; where="${s}"
      case "${st}" in
        failed)
          detail="$(marker_value "${d}/status/${s}.failed" reason); $(grep -m1 '^ERROR:' "${d}/logs/${s}.log" 2>/dev/null | cut -c1-100)"
          detail="${detail%; } | ${d}/logs/${s}.log" ;;
        running)
          detail="since $(marker_value "${d}/status/${s}.running" start) (job $(marker_value "${d}/status/${s}.running" job))" ;;
        interrupted)
          detail="job $(marker_value "${d}/status/${s}.running" job) ended without a status (killed, timed out or out of memory?) | ${d}/logs/${s}.log" ;;
      esac
      break
    done
    printf '%s\t%s\t%s\t%s\n' "${id}" "${overall}" "${where}" "${detail}"
  done < "${ids}"
}
rows=$(build_rows)

if [ -n "${LIST}" ]; then
  awk -F'\t' -v s="${LIST}" '$2 == s {print $1}' <<< "${rows}"
  exit 0
fi

{
  printf 'ID\tSTATE\tSTAGE\tDETAIL\n'
  echo "${rows}"
} | awk -F'\t' '{ printf "%-10s %-12s %-10s %s\n", $1, $2, $3, $4 }'
echo
awk -F'\t' '{ n[$2]++ } END { for (s in n) printf "%s: %d  ", s, n[s]; print "" }' <<< "${rows}"
