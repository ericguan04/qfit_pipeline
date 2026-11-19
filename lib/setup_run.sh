# Run setup shared by submit_slurm.sh and run_local.sh. Sourced, not executed.
#
# setup_run RUN_DIR CONFIG IDS_FILE (CONFIG may be empty) sets:
#   RUN_DIR    absolute run directory
#   SUB_DIR    run_info/submissions/<time>/ for this launch
#   TODO_FILE  IDs in this launch: those without a finished post_qfit stage
#   N_TODO     number of lines in TODO_FILE

die() { echo "ERROR: $*" >&2; exit 1; }

# shellcheck source=state.sh
source "${REPO}/lib/state.sh"

# Upper-case IDs, one per line, without blanks, comments, CR characters or
# duplicates. awk also keeps a last line that has no newline.
normalize_ids() {
  awk '{ sub(/\r$/, ""); sub(/#.*/, ""); gsub(/[ \t]+/, "");
         if ($0 == "") next; u = toupper($0); if (!(u in seen)) { seen[u] = 1; print u } }' "$1"
}

# Runs a bash command where the stages run: in the container, or natively.
in_env() {
  if [ -n "${CONTAINER}" ]; then
    apptainer exec --bind "${APPTAINER_BINDS}" "${CONTAINER}" bash -c "$1"
  else
    bash -c "$1"
  fi
}

# Checks that the tools work and writes $1 with everything needed to say
# later which code produced a run.
check_and_record_env() {
  local out="$1" err qfit_version qfit_commit
  if [ -n "${CONTAINER}" ]; then
    command -v apptainer >/dev/null 2>&1 || die "apptainer not found (CONTAINER is set in the config)."
    [ -f "${CONTAINER}" ] || die "container not found: ${CONTAINER}"
  fi
  # Judged by exit code, not by output: apptainer and Python may print
  # harmless warnings on stderr.
  local rc=0
  err=$(in_env 'bad=0
                for c in phenix.refine phenix.mtz.dump phenix.ready_set phenix.composite_omit_map \
                         qfit_protein qfit_final_refine_xray.sh create_restraints_file.py \
                         redistribute_cull_low_occupancies; do
                  command -v "$c" >/dev/null || { echo "not found: $c"; bad=1; }
                done
                [ -n "${PHENIX:-}" ] || { echo "PHENIX is not set"; bad=1; }
                out=$(python3 -c "import qfit, qfit.solvers" 2>&1) || { echo "cannot import qfit:"; echo "$out" | tail -n 3; bad=1; }
                exit $bad' 2>&1) || rc=$?
  [ ${rc} -eq 0 ] || die "the $( [ -n "${CONTAINER}" ] && echo container || echo environment ) is not usable:
${err}"

  qfit_version=$(in_env 'python3 -m pip show qfit 2>/dev/null | awk "/^Version:/ {print \$2}"')
  qfit_commit=$(in_env 'cat /opt/qfit-3.0.commit 2>/dev/null' || true)
  {
    echo "date:            $(date)"
    echo "user@host:       $(whoami)@$(hostname)"
    echo "pipeline commit: $(git -C "${REPO}" log -1 --format='%h %s' 2>/dev/null || echo unknown)"
    echo "pipeline dirty:  $(git -C "${REPO}" status --porcelain 2>/dev/null | wc -l | tr -d ' ') changed files"
    if [ -n "${CONTAINER}" ]; then
      echo "container:       $(ls -lL "${CONTAINER}")"
      apptainer inspect --labels "${CONTAINER}" 2>/dev/null | sed 's/^/  label: /'
    else
      echo "container:       none (native environment)"
    fi
    echo "qfit version:    ${qfit_version}"
    echo "qfit commit:     ${qfit_commit:-unknown}"
    in_env 'echo "PHENIX:          ${PHENIX}"
            echo "qfit_protein:    $(command -v qfit_protein)"
            phenix.version 2>/dev/null | head -n 3 | sed "s/^/phenix.version:  /"' 2>/dev/null
  } > "${out}"

  if [ -n "${QFIT_EXPECTED_COMMIT}" ]; then
    local want="${QFIT_EXPECTED_COMMIT:0:7}"
    case "${qfit_commit}:${qfit_version}" in
      "${want}"*|*"+g${want}"*) ;;
      *) die "qFit is ${qfit_version} (commit ${qfit_commit:-unknown}), not the expected ${QFIT_EXPECTED_COMMIT}." ;;
    esac
  fi
}

setup_run() {
  local run_dir="$1" config="$2" ids_file="$3" bad id
  [ -f "${ids_file}" ] || die "ID list not found: ${ids_file}"
  [ -z "${config}" ] || [ -f "${config}" ] || die "config not found: ${config}"
  mkdir -p "${run_dir}/run_info" "${run_dir}/structures" || die "cannot create ${run_dir}"
  RUN_DIR=$(cd "${run_dir}" && pwd)

  # The run keeps the config it was created with.
  if [ -f "${RUN_DIR}/run_info/config.sh" ]; then
    if [ -n "${config}" ] && ! cmp -s "${config}" "${RUN_DIR}/run_info/config.sh"; then
      die "${RUN_DIR} already uses the config in run_info/config.sh, which differs from ${config}.
Leave out -c to continue this run, or use a new run directory for new settings."
    fi
  else
    cp "${config:-${REPO}/config.sh}" "${RUN_DIR}/run_info/config.sh"
  fi
  # shellcheck source=/dev/null
  source "${RUN_DIR}/run_info/config.sh"
  [ -z "${INPUT_DIR}" ] || [ -d "${INPUT_DIR}" ] || die "INPUT_DIR not found: ${INPUT_DIR}"

  SUB_DIR="${RUN_DIR}/run_info/submissions/$(date '+%Y%m%d_%H%M%S')"
  mkdir -p "${SUB_DIR}"
  normalize_ids "${ids_file}" > "${SUB_DIR}/requested_ids.txt"
  bad=$(grep -v -E '^[A-Z0-9][A-Z0-9_.-]*$' "${SUB_DIR}/requested_ids.txt")
  [ -z "${bad}" ] || die "invalid IDs in ${ids_file}: $(echo ${bad})"
  [ -s "${SUB_DIR}/requested_ids.txt" ] || die "no IDs in ${ids_file}."

  # All IDs ever requested in this run, in first-seen order (status.sh reads it).
  touch "${RUN_DIR}/run_info/ids.txt"
  cat "${RUN_DIR}/run_info/ids.txt" "${SUB_DIR}/requested_ids.txt" | awk '!seen[$0]++' \
    > "${RUN_DIR}/run_info/ids.txt.tmp" && mv "${RUN_DIR}/run_info/ids.txt.tmp" "${RUN_DIR}/run_info/ids.txt"

  # Skip finished structures and ones a live job is still working on.
  TODO_FILE="${SUB_DIR}/ids.txt"
  : > "${TODO_FILE}"
  local s state busy
  while IFS= read -r id; do
    [ -f "${RUN_DIR}/structures/${id}/status/post_qfit.done" ] && continue
    busy=""
    for s in ${STAGES}; do
      state=$(stage_state "${RUN_DIR}/structures/${id}" "${s}")
      [ "${state}" = running ] && busy="${s}"
    done
    if [ -n "${busy}" ]; then
      echo "Skipping ${id}: ${busy} is still running." >&2
      continue
    fi
    echo "${id}" >> "${TODO_FILE}"
  done < "${SUB_DIR}/requested_ids.txt"
  N_TODO=$(wc -l < "${TODO_FILE}" | tr -d ' ')

  check_and_record_env "${SUB_DIR}/versions.txt"
}
