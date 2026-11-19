#!/bin/bash
# Tests the pipeline logic with stub versions of Phenix, qFit, apptainer,
# sbatch and squeue, so it runs anywhere in seconds. Nothing scientific is
# tested here: the stubs only write the files the real tools would.
# Usage: tests/run_tests.sh
set -uo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
T=$(mktemp -d "${TMPDIR:-/tmp}/qfit_pipeline_test.XXXXXX")
# KEEP=1 tests/run_tests.sh keeps the test directory for debugging.
trap 'if [ -n "${KEEP:-}" ]; then echo "Kept ${T}"; else rm -rf "${T}"; fi' EXIT
PASS=0; FAIL=0

check() {  # check DESCRIPTION COMMAND...
  if "${@:2}" >/dev/null 2>&1; then PASS=$((PASS + 1)); echo "  ok   $1"
  else FAIL=$((FAIL + 1)); echo "  FAIL $1"; fi
}
contains() { grep -q -- "$2" "$1"; }

#______________________________ STUBS ______________________________
B="${T}/bin"; mkdir -p "${B}" "${T}/py/qfit"
touch "${T}/py/qfit/__init__.py" "${T}/py/qfit/solvers.py"
stub() { printf '#!/bin/bash\n%s\n' "$2" > "${B}/$1"; chmod +x "${B}/$1"; }
here='id=$(basename "$PWD" | tr A-Z a-z)'

stub phenix.mtz.dump 'printf "%s\n" "    label  #valid  %valid  min  max type" \
  "    FOBS     100 100.00%  1 2 F: amplitude" "    SIGFOBS  100 100.00%  1 2 Q: standard deviation" \
  "    R-free-flags-1 100 100.00% 0 19 I: integer" "    R-free-flags 100 100.00% 0 1 I: integer"'
stub phenix.fetch_pdb 'echo "fetch should not be called" >&2; exit 1'
# No updated file for 2XYZ, to test the ready_set fallback.
stub phenix.ready_set "${here}"'; [ "$id" = 2xyz ] || cp "$1" "${1%.pdb}.updated.pdb"'
stub phenix.refine 'printf "%s\n" "$@" > refine_args.txt
b=$(basename "$1" .pdb)_refine_001; cp "$1" $b.pdb
printf "|--(resolution: 1.5 - 30 A, n_refl.=1000 (all), 5.00  %% free)--|\nStart R-work = 0.2000, R-free = 0.2300\nFinal R-work = 0.1900, R-free = 0.2200\n" > $b.log'
stub phenix.composite_omit_map 'echo map > composite_omit_map.mtz'
# phenix.python: remove_copies.py runs for real (standard library only);
# standardize_setting.py (needs cctbx) is only recorded.
stub phenix.python 'case "$1" in *remove_copies.py) exec python3 "$@" ;; esac
printf "%s\n" "$@" >> std_args.txt; echo "standardize_setting: stub"'
stub qfit_protein "${here}"'; printf "%s\n" "$@" > qfit_args.txt
[ "$id" = "${STUB_FAIL_QFIT:-}" ] && { echo "qfit crashed" >&2; exit 3; }
# Residue 1: conformer B is an exact copy of C.
printf "%s\n" \
"ATOM      1  CA ALEU A   1       1.000   2.000   3.000  0.34 20.00           C" \
"ATOM      2  CA BLEU A   1       1.500   2.000   3.000  0.21 20.00           C" \
"ATOM      3  CA CLEU A   1       1.500   2.000   3.000  0.45 20.00           C" END > multiconformer_model2.pdb'
stub qfit_final_refine_xray.sh 'n=${1%.mtz}; echo final > ${n}_qFit.pdb; echo mtz > ${n}_qFit.mtz
printf "|--(resolution: 1.5 - 30 A, n_refl.=1000 (all), ${STUB_FINAL_FREE:-5.00}  %% free)--|\nStart R-work = 0.21, R-free = 0.2400\nFinal R-work = 0.18, R-free = 0.2250\n" > ${n}_qFit.log'
for c in create_restraints_file.py redistribute_cull_low_occupancies remove_duplicates; do stub "$c" 'exit 0'; done
stub sbatch 'n=$(cat '"${T}"'/sbatch_n 2>/dev/null || echo 100); n=$((n + 1)); echo $n > '"${T}"'/sbatch_n
printf "%s\n" "$@" > '"${T}"'/sbatch_$n.args; echo "$n;cluster1"'
stub squeue 'case "$*" in *"${STUB_ALIVE_JOB:-none}"*) echo "  $STUB_ALIVE_JOB R" ;; esac'
# apptainer exec --bind BINDS IMAGE CMD...  ->  run CMD natively
stub apptainer 'case "$1" in
  exec) echo "$@" >> '"${T}"'/apptainer_calls; shift 4; exec "$@" ;;
  inspect) echo "qfit.commit abc1234" ;; esac'
export PATH="${B}:${PATH}" PYTHONPATH="${T}/py" PHENIX=/stub/phenix

# Inputs, laid out like test_data/
for id in 1abc 2xyz; do
  ID=$(echo ${id} | tr a-z A-Z); mkdir -p "${T}/inputs/${ID}"
  echo "ATOM" > "${T}/inputs/${ID}/${id}.pdb"; echo "mtz" > "${T}/inputs/${ID}/${id}.mtz"
done
cfg() {  # cfg FILE [CONTAINER]
  sed -e "s|^CONTAINER=.*|CONTAINER=\"${2:-}\"|" -e "s|^INPUT_DIR=.*|INPUT_DIR=\"${T}/inputs\"|" \
      -e 's|^QFIT_CPUS=.*|QFIT_CPUS=3|' \
      -e 's|^APPTAINER_BINDS=.*|APPTAINER_BINDS="/data,/scratch"|' "${REPO}/config.sh" > "$1"
}
cfg "${T}/config.sh"
# IDs with blank lines, CRLF, lower case, a duplicate, a comment, no final newline
printf '1abc\r\n\n  # comment\n2XYZ\n1ABC\n2xyz' > "${T}/ids.txt"

#______________________________ TESTS ______________________________
echo "1. Local run, one structure's qFit fails"
R="${T}/run1"
STUB_FAIL_QFIT=2xyz bash "${REPO}/run_local.sh" -r "${R}" -c "${T}/config.sh" "${T}/ids.txt" > "${T}/out1" 2>&1
check "IDs normalized (1ABC, 2XYZ)" test "$(cat "${R}/run_info/ids.txt" | tr '\n' ,)" = "1ABC,2XYZ,"
check "config frozen in run" cmp -s "${T}/config.sh" "${R}/run_info/config.sh"
check "versions recorded" contains "${R}"/run_info/submissions/*/versions.txt "pipeline commit"
for s in pre_qfit qfit post_qfit; do check "1ABC ${s} done" test -f "${R}/structures/1ABC/status/${s}.done"; done
check "1ABC final model" test -s "${R}/structures/1ABC/1abc_qFit.pdb"
check "pre-qFit refine got R-free-flags (not -1)" contains "${R}/structures/1ABC/refine_args.txt" "r_free_flags.label=R-free-flags$"
check "weight optimization off by default" bash -c "! grep -q optimize_xyz '${R}/structures/1ABC/refine_args.txt'"
check "qFit used QFIT_CPUS=3" contains "${R}/structures/1ABC/qfit_args.txt" "^3$"
check "post_qfit log compares final R-free with pre-qFit" contains "${R}/structures/1ABC/logs/post_qfit.log" "R-free: pre-qFit 0.2200, final 0.2250, change +0.0050"
check "qFit got --scale-rmask from config" contains "${R}/structures/1ABC/qfit_args.txt" "^--scale-rmask$"
check "qFit raw output kept" contains "${R}/structures/1ABC/multiconformer_model2.qfit_raw.pdb" " CA BLEU"
check "exact copy merged before final refinement" bash -c "! grep -q ' CA BLEU' '${R}/structures/1ABC/multiconformer_model2.pdb' && grep -q ' CA CLEU.* 0.66' '${R}/structures/1ABC/multiconformer_model2.pdb'"
check "qfit log reports merged copies" contains "${R}/structures/1ABC/logs/qfit.log" "merged 1 exact-copy conformers in 1 residues"
check "space-group standardisation ran on raw inputs" contains "${R}/structures/1ABC/std_args.txt" "standardize_setting.py"
check "2XYZ ready_set fallback copied the model" test -s "${R}/structures/2XYZ/2xyz.updated_refine_001.pdb"
check "2XYZ qfit failed marker" test -f "${R}/structures/2XYZ/status/qfit.failed"
check "2XYZ post_qfit not run" test ! -e "${R}/structures/2XYZ/status/post_qfit.done"
check "failed marker has exit code 3" contains "${R}/structures/2XYZ/status/qfit.failed" "exit_code=3"
bash "${REPO}/status.sh" -r "${R}" > "${T}/status1"
check "status: 1ABC complete" contains "${T}/status1" "^1ABC  *complete"
check "status: 2XYZ failed at qfit with log" contains "${T}/status1" "^2XYZ  *failed  *qfit .*command failed with exit code 3.*logs/qfit.log"
check "status -l failed lists 2XYZ" test "$(bash "${REPO}/status.sh" -r "${R}" -l failed)" = "2XYZ"

echo "2. Rerun the same command: only 2XYZ, pre_qfit skipped"
bash "${REPO}/run_local.sh" -r "${R}" "${T}/ids.txt" > "${T}/out2" 2>&1
check "only 2XYZ to run" contains "${T}/out2" "Structures: 1 to run"
check "2XYZ pre_qfit skipped" contains "${T}/out2" "2XYZ pre_qfit: already done"
check "2XYZ now complete" test -f "${R}/structures/2XYZ/status/post_qfit.done"
check "previous qfit log kept" bash -c "ls '${R}'/structures/2XYZ/logs/previous/qfit.*.log"
check "old failed marker removed" test ! -e "${R}/structures/2XYZ/status/qfit.failed"
check "status: all complete" contains <(bash "${REPO}/status.sh" -r "${R}") "complete: 2"
check "different -c on existing run refused" bash -c "! bash '${REPO}/run_local.sh' -r '${R}' -c '${REPO}/config.sh' '${T}/ids.txt'"

echo "3. Final refinement on a different free set"
R="${T}/run3"
STUB_FINAL_FREE=9.07 bash "${REPO}/run_local.sh" -r "${R}" -c "${T}/config.sh" "${T}/ids.txt" > "${T}/out3" 2>&1
check "post_qfit failed" test -f "${R}/structures/1ABC/status/post_qfit.failed"
check "reason in log" contains "${R}/structures/1ABC/logs/post_qfit.log" "different free set (9.07 %) from the pre-qFit refinement (5.00 %)"

echo "4. Rerunning a stage invalidates later stages"
R="${T}/run1"
rm "${R}/structures/1ABC/status/qfit.done"
bash "${REPO}/lib/run_stage.sh" "${R}" qfit 1ABC > /dev/null 2>&1
check "qfit done again" test -f "${R}/structures/1ABC/status/qfit.done"
check "post_qfit.done removed" test ! -e "${R}/structures/1ABC/status/post_qfit.done"
check "status: 1ABC pending at post_qfit" contains <(bash "${REPO}/status.sh" -r "${R}") "^1ABC  *pending  *post_qfit"

echo "5. Interrupted and running stages"
printf 'host=%s\npid=999999\njob=\nstart=x\n' "$(hostname)" > "${R}/structures/1ABC/status/post_qfit.running"
check "dead local pid -> interrupted" contains <(bash "${REPO}/status.sh" -r "${R}") "^1ABC  *interrupted"
printf 'host=node1\npid=1\njob=555_1\nstart=x\n' > "${R}/structures/1ABC/status/post_qfit.running"
check "job gone from squeue -> interrupted" contains <(bash "${REPO}/status.sh" -r "${R}") "^1ABC  *interrupted"
check "job in squeue -> running" contains <(STUB_ALIVE_JOB=555_1 bash "${REPO}/status.sh" -r "${R}") "^1ABC  *running"
STUB_ALIVE_JOB=555_1 bash "${REPO}/run_local.sh" -r "${R}" "${T}/ids.txt" > "${T}/out5" 2>&1
check "running structure not relaunched" contains "${T}/out5" "Skipping 1ABC: post_qfit is still running"
bash "${REPO}/run_local.sh" -r "${R}" "${T}/ids.txt" > "${T}/out5b" 2>&1
check "interrupted structure relaunched and finished" test -f "${R}/structures/1ABC/status/post_qfit.done"

echo "6. SLURM submission (stub sbatch) and array tasks"
R="${T}/run6"; rm -f "${T}"/sbatch_*
bash "${REPO}/submit_slurm.sh" -r "${R}" -c "${T}/config.sh" "${T}/ids.txt" > "${T}/out6" 2>&1
check "three arrays submitted" test "$(ls "${T}"/sbatch_*.args | wc -l | tr -d ' ')" = 3
check "pre array 1-2%100, 1 cpu" bash -c "grep -q -- '--array=1-2%100' '${T}/sbatch_101.args' && grep -q -- '--cpus-per-task=1' '${T}/sbatch_101.args'"
check "qfit after pre (aftercorr:101)" contains "${T}/sbatch_102.args" "^--dependency=aftercorr:101$"
check "post after qfit (aftercorr:102)" contains "${T}/sbatch_103.args" "^--dependency=aftercorr:102$"
check "qfit array 3 cpus, 18G" bash -c "grep -q -- '--cpus-per-task=3' '${T}/sbatch_102.args' && grep -q -- '--mem=18G' '${T}/sbatch_102.args'"
check "job ids recorded without cluster suffix" contains "${R}"/run_info/submissions/*/jobs.txt "^qfit=102$"
ids_file=$(ls "${R}"/run_info/submissions/*/ids.txt)
for s in pre_qfit qfit post_qfit; do
  for i in 1 2; do
    QFIT_PIPELINE_REPO="${REPO}" SLURM_ARRAY_JOB_ID=7 SLURM_ARRAY_TASK_ID=$i SLURM_CPUS_PER_TASK=2 \
      bash "${REPO}/lib/slurm_task.sh" "${s}" "${R}" "${ids_file}" > /dev/null 2>&1
  done
done
check "array tasks completed both structures" contains <(bash "${REPO}/status.sh" -r "${R}") "complete: 2"
check "SLURM cores override config (2)" contains "${R}/structures/2XYZ/qfit_args.txt" "^2$"
check "job recorded in done marker" contains "${R}/structures/2XYZ/status/qfit.done" "job=7_2"
bash "${REPO}/submit_slurm.sh" -r "${R}" "${T}/ids.txt" > "${T}/out6b" 2>&1
check "resubmitting a finished run submits nothing" contains "${T}/out6b" "Nothing to submit"

echo "7. Container mode and environment checks"
R="${T}/run7"; touch "${T}/image.sif"; cfg "${T}/config_c.sh" "${T}/image.sif"
bash "${REPO}/run_local.sh" -r "${R}" -c "${T}/config_c.sh" "${T}/ids.txt" > "${T}/out7" 2>&1
check "stages ran through apptainer exec --bind" contains "${T}/apptainer_calls" "exec --bind /data,/scratch ${T}/image.sif bash ${REPO}/stages/pre_qfit.sh"
check "container labels recorded" contains "${R}"/run_info/submissions/*/versions.txt "label: qfit.commit abc1234"
check "complete in container mode" contains <(bash "${REPO}/status.sh" -r "${R}") "complete: 2"
sed 's|^QFIT_EXPECTED_COMMIT=.*|QFIT_EXPECTED_COMMIT="deadbeef"|' "${T}/config.sh" > "${T}/config_bad.sh"
check "wrong qFit commit refused" bash -c "! bash '${REPO}/run_local.sh' -r '${T}/run7b' -c '${T}/config_bad.sh' '${T}/ids.txt'"
echo 'raise ImportError("numpy.lib.array_utils")' > "${T}/py/qfit/solvers.py"
bash "${REPO}/run_local.sh" -r "${T}/run7c" -c "${T}/config.sh" "${T}/ids.txt" > "${T}/out7c" 2>&1
check "broken qFit import stops before running" contains "${T}/out7c" "array_utils"
check "... and no structure was started" test ! -d "${T}/run7c/structures/1ABC"

echo "8. Space-group standardisation (real cctbx, not stubs)"
# Needs a Python with cctbx; skipped unless CCTBX_PYTHON is set, e.g.
#   CCTBX_PYTHON=phenix.python tests/run_tests.sh
# (PATH still holds the stubs here, so give a full path for phenix.python.)
if [ -n "${CCTBX_PYTHON:-}" ]; then
  "${CCTBX_PYTHON}" "${REPO}/tests/test_standardize_setting.py" > "${T}/out8" 2>&1
  check "standardize_setting.py tests (see tests/test_standardize_setting.py)" contains "${T}/out8" "all standardize_setting tests passed"
  grep FAIL "${T}/out8" || true
else
  echo "  skip (set CCTBX_PYTHON to run)"
fi

echo "9. Exact-copy removal (real script, standard library only)"
python3 "${REPO}/tests/test_remove_copies.py" > "${T}/out9" 2>&1
check "remove_copies.py tests (see tests/test_remove_copies.py)" contains "${T}/out9" "all remove_copies tests passed"
grep FAIL "${T}/out9" || true

echo
echo "${PASS} passed, ${FAIL} failed"
[ ${FAIL} -eq 0 ]
