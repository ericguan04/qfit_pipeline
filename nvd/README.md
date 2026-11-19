# NVD runs (branch `nvd-experiments`)

Project-specific files for the NVD natural vs designed protein comparison.
The pipeline itself is unchanged from `main`; this folder only adds configs,
ID lists and analysis helpers. To pick up later pipeline fixes:
`git merge main`.

| File | What it is |
|---|---|
| `config_nvd.sh` | config for the HPC runs (container and bind paths filled in), `--scale-rmask 1.0` |
| `config_nvd_rmask2.sh` | identical, except `QFIT_SCALE_RMASK=2.0` |
| `scale_rmask_subset.txt` | 15 structures for the `--scale-rmask` test |
| `compare_rfree.sh` | tabulates R-free (deposited, pre-qFit, final) and free-set % per structure |
| `REVISED_PIPELINE.md` | the five-stage pipeline and the evidence for each change |

## 1. Main run (`--scale-rmask 1.0`)

```bash
./submit_slurm.sh -r /data/sbi/eguan/qfit_runs/nvd_rmask1 -c nvd/config_nvd.sh nvd_ids.txt
./status.sh -r /data/sbi/eguan/qfit_runs/nvd_rmask1
```

Check that 5GYC and 5GYB were standardised (other structures should report
"already the standard setting"). A structure whose final refinement used a
different free set fails in `post_qfit` and shows up in `status.sh`.

```bash
cd /data/sbi/eguan/qfit_runs/nvd_rmask1/structures
grep -H "standardize_setting:" 5GYC/logs/pre_qfit.log 5GYB/logs/pre_qfit.log
```

## 2. `--scale-rmask` test

qFit's map scaling leaves the solvent region too high in many structures,
which adds spurious alternates. `--scale-rmask 2.0` corrected this in sampling
tests (multi-conformer residues 20.6% → 6.5% on 28 structures) but has not
been checked by refinement. This test does that.

**Decide before looking at natural vs designed results.** Run the 2.0 arm if
the main run still shows, for example, multi-conformer rates tracking
resolution or the scaled solvent level, or the same outliers (7UR7, 7L1Y,
1N0R, 6W70) far above the median.

The subset spans 1.15–2.0 Å, 4 synthetic and 11 natural structures. It covers:

- high-background outliers: 7UR7, 7L1Y, 6W70, 1N0R
- structures that lost deposited alternates at 2.0: 2IC6, 2VM1, 4LG6, 7MCC
- low-background structures with deposited alternates (should not change): 4OO5, 6S18
- low-background controls: 1CPQ, 1FNJ, 3NGJ, 5L33
- 5GYC (space-group standardisation check)

```bash
./submit_slurm.sh -r /data/sbi/eguan/qfit_runs/rmask2_subset -c nvd/config_nvd_rmask2.sh nvd/scale_rmask_subset.txt
```

The 1.0 arm is the main run, so nothing else needs to be submitted. Compare
the two runs on the subset:

```bash
nvd/compare_rfree.sh /data/sbi/eguan/qfit_runs/nvd_rmask1
nvd/compare_rfree.sh /data/sbi/eguan/qfit_runs/rmask2_subset
```

Also compare multi-conformer counts, and look for difference density
(mFo-DFc) where 2.0 removed alternates that 1.0 kept, especially the deposited
alternates 2IC6 B24, 2VM1 C102, 4LG6 A284, 7L1Y A54 and 7MCC A111.

**Proposed criterion:** adopt 2.0 if final R-free does not worsen (within
about 0.002) and no significant difference density appears where alternates
were removed. If adopted, rerun all structures with `config_nvd_rmask2.sh` in a
new run directory, and report the 1.0 run as the comparison.
