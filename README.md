# qfit_pipeline

Batch qFit runs on a SLURM cluster (or one machine), using Phenix 2.0 and a
fork of qFit with fixes to the final refinement step.

For each structure, five stages run in order:

| # | Pipeline stage | Job | What it does | Main outputs |
|---|---|---|---|---|
| 1 | Standardise setting | `pre_qfit` | fetch inputs (or copy them from `INPUT_DIR`), move model and data to the standard space-group setting | `<id>.pdb`, `<id>.mtz` (originals: `<id>.deposited_setting.*`, if changed) |
| 2 | Pre-qFit refinement | `pre_qfit` | `phenix.ready_set`, `phenix.refine` (8 macro cycles, sites + ADP + occupancies, riding H), `phenix.composite_omit_map` | `<id>.updated_refine_001.pdb`, `composite_omit_map.mtz` |
| 3 | qFit | `qfit` | `qfit_protein` on the composite omit map, default settings (`--scale-rmask` from the config, 1.0) | `multiconformer_model2.qfit_raw.pdb` |
| 4 | Remove exact copies | `qfit` | merge alternate conformers that are exact copies of another (`stages/remove_copies.py`) | `multiconformer_model2.pdb` |
| 5 | Final refinement | `post_qfit` | `qfit_final_refine_xray.sh` from the fork | `<id>_qFit.pdb`, `<id>_qFit.mtz`, `<id>_qFit.log` |

The five stages run as three jobs per structure (`pre_qfit`, `qfit`,
`post_qfit`).

There are no R-free filters. A stage fails only when a tool fails or an
expected output is missing. One consistency check: `post_qfit` fails if the
final refinement used a different free set than the pre-qFit refinement. Its
log records the final R-free next to the pre-qFit R-free.

**Space-group setting.** qFit rebuilds space groups from their number, so a
structure deposited in a non-standard setting (e.g. P 2 21 21) is scaled with
the wrong symmetry operators. `pre_qfit` therefore re-expresses the model and
data in the standard setting (`stages/standardize_setting.py`; a change of axes
only, R-free flags included). Structures already in a standard setting are not
touched. The originals are kept as `<id>.deposited_setting.{pdb,mtz}`, and all
outputs are in the standard setting.

**Exact copies.** qFit's segment step pads residues with exact copies of their
conformers (its identical-conformer removal is commented out), and refinement
then moves the copies apart. Before final refinement, `remove_copies.py`
merges any altloc whose heavy atoms lie within 0.01 Å RMSD (qFit's own
`--rmsd-cutoff`) of a higher-occupancy altloc, adding its occupancy to that
one. The calculated structure factors do not change. The `qfit` log reports
how many conformers were merged; qFit's own output is kept as
`multiconformer_model2.qfit_raw.pdb`.

## 1. Build the container (once per version)

The image holds Phenix and qFit, so the cluster needs nothing else installed.
Phenix cannot be downloaded automatically or shared, so you supply its
installer.

```bash
git clone https://github.com/ericguan04/qfit_pipeline.git
cd qfit_pipeline/container
ln -s /path/to/Phenix-2.0-5936-Linux-x86_64.sh .     # your Phenix installer
apptainer build qfit_pipeline.sif qfit.def          # add --fakeroot if required
apptainer inspect --labels qfit_pipeline.sif        # qFit commit, Phenix version
mv qfit_pipeline.sif /shared/lab/path/              # anywhere jobs can read
```

The build ends with a self-test (qFit imports, all tools on `PATH`); run it again
with `apptainer test qfit_pipeline.sif`. The qFit commit is set by `QFIT_REF`
in `container/qfit.def`.

**qFit version.** The container uses a fork of qFit
([ericguan04/qfit-3.0](https://github.com/ericguan04/qfit-3.0), branch `dev`),
pinned to one commit. It is upstream qFit commit `33fd3e9` with three fixes to
the final refinement scripts (occupancy constraint groups, cull-loop
convergence, R-free flag preservation); qFit's sampling code is unchanged.
Once the fixes are merged upstream, `QFIT_REF` will point to upstream qFit.

## 2. Configure

Copy `config.sh` and edit at least `CONTAINER` and `APPTAINER_BINDS`. The binds
must cover this repository, your run directories and `INPUT_DIR`, if used.
Optionally set `SLURM_ACCOUNT`, `SLURM_PARTITION`, resources and `INPUT_DIR`.

```bash
cp config.sh ~/my_config.sh
```

## 3. Run

Write the structure IDs into a file, one per line:

```bash
./submit_slurm.sh -r /data/$USER/qfit_runs/batch1 -c ~/my_config.sh ids.txt
```

This submits three job arrays. Task *i* of each stage starts when task *i* of
the previous stage succeeded. Use `-n` to set up the run and print the `sbatch`
commands without submitting.

On a single machine (no SLURM), use the same options; `-p N` runs N structures
at a time:

```bash
./run_local.sh -r ~/qfit_runs/test -c ~/my_config.sh -p 2 ids.txt
```

Leave `CONTAINER` empty to use Phenix and qFit from the current environment
instead (`source phenix_env.sh` and activate the qFit conda env first).

## 4. Check progress and retry

```bash
./status.sh -r /data/$USER/qfit_runs/batch1
```

```
ID         STATE        STAGE      DETAIL
1ABC       complete     -
2XYZ       failed       qfit       exit code 3; ERROR: ... | .../structures/2XYZ/logs/qfit.log
3DEF       running      post_qfit  since 2026-11-18 09:12:40 (job 4411_3)
4GHI       interrupted  qfit       job 4410_4 ended without a status (killed, timed out or out of memory?) | ...
```

`-l failed` (or `interrupted`, `pending`, ...) prints just those IDs.

**To retry, run the same submit command again** (without `-c`). Finished
structures are skipped, and so are finished stages: a structure that failed in
`qfit` restarts at `qfit`.

## Results

Each structure's files are in `<run_dir>/structures/<ID>/`:

| File | What it is |
|---|---|
| `<id>_qFit.pdb`, `<id>_qFit.mtz` | final multiconformer model and map coefficients |
| `multiconformer_model2.qfit_raw.pdb` | qFit's own output, before exact copies were merged |
| `multiconformer_model2.pdb` | qFit output with exact copies merged (input to final refinement) |
| `<id>.updated_refine_001.pdb` | pre-qFit refined model (qFit's input) |
| `<id>.deposited_setting.{pdb,mtz}` | original files, only if the setting was changed |

Per-structure numbers worth checking are in the logs:

```bash
cd <run_dir>/structures
grep -H "remove_copies:" */logs/qfit.log     # exact-copy conformers merged
grep -H "^R-free:" */logs/post_qfit.log      # final R-free vs the pre-qFit model
grep -H "free set:" */logs/post_qfit.log     # free-set % in both refinements
grep -H "standardize_setting:" */logs/pre_qfit.log  # space-group setting changed or not
```

## Run directory

```
<run_dir>/
├── run_info/
│   ├── config.sh              # the config this run uses, frozen at creation
│   ├── ids.txt                # every ID ever requested in this run
│   └── submissions/<time>/    # one per submit: IDs, sbatch commands, job IDs,
│                              #   versions.txt (pipeline/qFit/Phenix/container), SLURM .out files
└── structures/<ID>/
    ├── <all stage inputs and outputs>
    ├── logs/<stage>.log        # full tool output; earlier attempts in logs/previous/
    └── status/<stage>.done | .failed | .running
```

Rules that keep runs trustworthy:

- **A stage counts as finished only if it wrote `status/<stage>.done`.** The
  marker is written after the stage exits 0 and its outputs exist and are not
  empty. Output files alone never count, and each stage deletes its old outputs
  before starting.
- **Rerunning a stage invalidates later stages.** Their `.done` markers are
  removed.
- **A run keeps the config it was created with.** To try other settings, start
  a new run directory. `submit_slurm.sh` refuses a different `-c` on an
  existing run.
- **Don't edit or delete files in a run by hand.** To redo a stage for a
  structure, delete its `status/<stage>.done` and resubmit.

## Tests

```bash
tests/run_tests.sh
```

These tests use stub versions of Phenix, qFit, `sbatch`, `squeue` and
`apptainer`, so they run anywhere in seconds. They check the pipeline logic
(markers, resume, failure reporting, SLURM arrays and dependencies, container
calls), not the science.

The space-group standardisation is also tested for real (structure factors,
R-free flags, refusals) when a Python with cctbx is given:

```bash
CCTBX_PYTHON=/path/to/phenix.python tests/run_tests.sh
```

The exact-copy removal is always tested for real (`tests/test_remove_copies.py`).
