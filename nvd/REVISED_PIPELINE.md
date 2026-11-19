# Revised qFit pipeline

The revised pipeline for multiconformer modelling of the NVD dataset (154 X-ray
structures, 1.0–2.0 Å), and the evidence for each change.

qFit's sampling code is used unmodified. Every structure gets exactly the same
steps. Each change either leaves the science unchanged (equivalent input,
correct bookkeeping, removing exact copies) or fixes a script that did not do
what it claimed.

Full investigation, with all tests and per-structure data:
[qFit modelling issues: NVD 154-structure run](https://claude.ai/artifact/P6dc1GVtTqgZWPE2tjcecx).

## Pipeline

| # | Stage | What it does | Changed? |
|---|---|---|---|
| 1 | Standardise setting | Re-express the deposited model and data in the standard space-group setting | **New** (Issue 1) |
| 2 | Pre-qFit refinement | `phenix.ready_set`, `phenix.refine`, composite omit map | Unchanged |
| 3 | qFit | `qfit_protein`, default settings | Unchanged |
| 4 | Remove exact copies | Merge alternate conformers that are exact copies of another conformer | **New** (Issue 5) |
| 5 | Final refinement | `qfit_final_refine_xray.sh` from the fork's `dev` branch | **Bug fixes** (Issues 2–4) |

qFit and the final-refinement scripts come from the fork
`ericguan04/qfit-3.0`, branch `dev`, commit `f9edb86`. That commit is upstream
`33fd3e9`, the version used in the original runs, plus three fixes to the
post-refinement scripts. The sampling code is byte-for-byte identical to
upstream.

> **Order change from the earlier plan.** Removing copies comes **before** the
> final refinement, not after it. The copies are exact in qFit's raw output, but
> refinement moves them apart, after which they can no longer be told apart
> from genuine close alternates (see Issue 5).

---

## Stage 1: Standardise the space-group setting

### Issue 1: qFit loses non-standard space-group settings

`UnitCell.copy()` rebuilds the space group from its number alone. A structure
deposited in a non-standard setting (for example P 2 21 21, No. 18) is
therefore treated as the standard setting (P 21 21 2) during map scaling.
That places three of the four symmetry mates wrongly, the map is scaled badly,
and qFit fills the mis-scaled density with spurious conformers.

**Affected:** 5GYC and 5GYB in this dataset. The bug is still present in
upstream qFit.

**Fix:** `standardize_setting.py` changes the axes of the model and the data
(R-free flags included) to the standard setting before any refinement. It is
a change of labels only: the structure and the measurements are identical.
Structures already in a standard setting are left untouched, so it runs on
every structure.

### Evidence

Map scaling for 5GYC (S = scale, k = offset; a typical structure has S ≈ 0.24,
k ≈ 0.5):

| 5GYC | S | k | Solvent level after scaling |
|---|---|---|---|
| As run | 0.09 | 1.15 | 1.11 (dataset median 0.41) |
| Bug patched in qFit's memory | 0.39 | 0.46 | – |
| **Standardised input, unmodified qFit** | **0.43** | **0.36** | – |

Full 5GYC run with the in-memory patch:

| 5GYC | As run | Patched |
|---|---|---|
| Multi-conformer residues | 94% | 43% |
| Residues with genuinely distinct alternates (> 0.5 Å) | 86% | 13% |
| Mean conformers per residue | 3.35 | 1.70 |

Standardising gave better scaling than the patch did, because the patch missed
a second place where qFit drops the setting. The full standardised run is part
of the HPC validation.

**Correctness checks:** structure factors computed from the model before and
after standardisation agree (exactly for axis swaps; within coordinate
rounding, about 0.1%, for monoclinic cases). R-free flags stay with their
reflections, and the free-set size is unchanged. Unit tests cover these
checks and the refusal cases.

---

## Stage 4: Remove exact-copy conformers

### Issue 5: qFit pads residues with exact copies

In qFit's segment step, the call that removes identical conformers is
commented out, although the log still says "after removal of identical
conformers". Each residue in a segment is then given the same set of
alternate labels. A residue with one or two real conformers is padded with
exact copies to fill the set.

The copies are not harmless. In final refinement they drift apart, creating
apparent heterogeneity from nothing.

**Fix:** after qFit and before final refinement, merge any alternate
conformer that is an exact copy of a higher-occupancy one. Its occupancy is
added to the conformer it copies. Copies are identified by heavy-atom RMSD
below 0.01 Å, which is qFit's own duplicate threshold (`--rmsd-cutoff`). This
restores the step qFit was meant to run.

### Evidence

qFit raw output, all 154 structures:

| Measure | Count |
|---|---|
| Multi-conformer residues | 6,087 of 52,107 (11.7%) |
| … with at least one exact copy | 1,951 (32%), in 68 structures |
| Distinct conformers added by the segment step | 0 |
| Alternate conformers that are exact copies | 2,977 of 15,094 (20%) |

The copies really are exact: cutoffs of 0.01 Å and 0.1 Å remove the same
conformers. The choice of cutoff therefore does not matter.

Effect on the top outliers (raw qFit output):

| PDB | Multi-conformer residues: before → after | Conformers per multi-conformer residue: before → after |
|---|---|---|
| 5GYC | 94.1% → 87.9% | 3.49 → 2.14 |
| 7UR7 | 82.5% → 79.4% | 3.27 → 2.08 |
| 5GYB | 80.5% → 74.6% | 3.23 → 2.07 |
| 7L1Y | 75.0% → 67.6% | 3.54 → 2.19 |
| 1N0R | 62.7% → 57.1% | 3.24 → 2.11 |
| 6W70 | 60.3% → 54.4% | 3.01 → 2.13 |
| Dataset | 11.7% → 11.1% | 2.48 → 2.04 |

Removing copies mainly removes phantom conformers. It changes the share of
multi-conformer residues only a little, because most padded residues also
have a real alternate. The excess of multi-conformer residues comes from the
map background (see the last section) and, in 5GYC and 5GYB, from Issue 1.

### Why before final refinement

Residues with a conformer pair within 0.1 Å or 0.3 Å, among multi-conformer
residues:

| Stage | < 0.1 Å | < 0.3 Å |
|---|---|---|
| qFit raw output | 32% | 38% |
| After stock final refinement | 7% | 31% |

In 7UR7, refined with the fixed scripts, residues with a pair within 0.1 Å
fell from 39 in the raw output to 4 after refinement.

After refinement, the former copies sit 0.1–0.3 Å apart, where no cutoff can
separate them from genuine close alternates. Before refinement they are exact,
so merging them is unambiguous. As a bonus, the final R-factors then describe
the final model directly.

---

## Stage 5: Final refinement with bug fixes

The stock final-refinement step has three bugs. Each distorts occupancies or
R-free, the quantities a multiconformer comparison relies on. The fork's
`dev` branch fixes them; nothing else in the step changes (the 0.09 occupancy
cull, refinement strategy, ADP choice and output names are all the same).

| Issue | Bug | Fix (commit) |
|---|---|---|
| 2 | Each altloc gets its own occupancy group, so altloc occupancies are not constrained to sum to 1 | One constrained group per multi-conformer segment, one selection per altloc (`763f52d`) |
| 3 | The cull-loop counter counts every residue, so the "refine until occupancies converge" loop always runs its maximum 50 rounds | Count only residues that lost a conformer (`11233c8`) |
| 4 | The first refinement call generates a new R-free set, discarding the one used before qFit | Reuse the existing R-free flags; generate new ones only if none exist, with a warning (`0e97a61`) |

### Evidence

From the original 154-structure run (stock scripts):

| Issue | Evidence |
|---|---|
| 2 | 5GYC: 236 of 251 multi-conformer residues have occupancy sums above 1.05 (mean 1.16, max 1.20). Inflated occupancies also keep conformers above the 0.09 cull. |
| 3 | Every completed run hit the 50-round limit. R-free rose over the loop in 137/154 structures (overfitting). |
| 4 | The free-set percentage changed in 141/154 structures (7UR7: 4.82% → 9.07%). R-free was below R-work at the start of final refinement in 81/154, the signature of a free set that includes reflections already fitted. |

Stock vs fixed scripts on the same qFit output:

| | 5GYC stock | 5GYC fixed | 7UR7 stock | 7UR7 fixed |
|---|---|---|---|---|
| Free reflections (pre-qFit / final) | 5.09% / 4.06% | 5.09% / 5.09% | 4.82% / 9.07% | 4.82% / 4.82% |
| Residues with occupancy sum > 1.05 | 236 | 0 | 0 | 0 |
| Mean occupancy sum | 1.16 | 1.00 | 0.99 | 1.00 |
| Cull-loop rounds | 50 | 4 | 50 | 2 |
| Final R-free (pre-qFit R-free, same free set) | not comparable | 0.1997 (0.1830) | not comparable | 0.2069 (0.2222) |

With the fixes, final R-free can be compared with the pre-qFit model for the
first time: 5GYC gets worse after qFit (its spurious conformers, Issue 1), and
7UR7 improves. The pipeline checks that both refinements used the same free set
and records the change in R-free for every structure.

**Why the fixes are defensible:** each restores what the script already
claimed to do, and each matches standard practice (occupancies of alternates
sum to 1; one free set throughout). All three have unit tests and will be
submitted upstream as separate pull requests.

---

## Not fixed: map background offset (`--scale-rmask`)

qFit scales the input map to the model only inside a mask around the atoms.
When the model–map fit is modest, this leaves the solvent region of the scaled
map too high (dataset median 0.41, against about 0.33 for a well-scaled map).
qFit then fills the excess density with extra conformers. This affects the
whole dataset, is worst at high resolution, and is the main source of excess
alternates outside 5GYC and 5GYB.

The pipeline runs qFit with default settings for now. The documented qFit
option `--scale-rmask` (enlarges the scaling mask) corrects the offset
structure by structure, and is worth testing.

Dataset-wide test (122 structures, 10 residues each), solvent level moved to
0.33:

| | As run | Corrected |
|---|---|---|
| Multi-conformer residues | 7.9% | 1.8% |
| Correlation with solvent level (Spearman) | +0.47 | +0.10 |
| Correlation with resolution (Spearman) | −0.50 | −0.30 |
| Deposited single-conformer residues given alternates | 7% | 1% |

qFit option test (28 structures weighted toward problem cases, 10 residues
each):

| Setting | Median solvent level | Multi-conformer residues | Deposited alternates kept |
|---|---|---|---|
| Default (`--scale-rmask 1.0`) | 0.48 | 20.6% | 9/18 |
| `--scale-rmask 1.5` | 0.36 | 9.0% | 6/18 |
| **`--scale-rmask 2.0`** | **0.33** | **6.5%** | **4/18** |
| `--scale-rmask 3.0` | 0.28 | 5.1% | 3/18 |

At 2.0, 7UR7 went from 7/10 multi-conformer residues to 0/10. Low-background
structures were unaffected. 4OO5 and 6S18 kept all their deposited alternates.

**Open question:** at 2.0, five deposited alternates were lost, all in
high-background structures. Whether they are real needs Phenix refinement.

**Proposed HPC test:** about 15 structures, full runs at `--scale-rmask` 1.0
and 2.0. Compare R-free, multi-conformer counts, and difference density where
alternates were removed. Adopt 2.0 if R-free does not worsen (within about
0.002) and no significant difference density appears.

**Caveat for the analysis:** without this correction, multi-conformer rates
partly track map background and resolution, not only biology. Control for
resolution in any natural vs designed comparison.
