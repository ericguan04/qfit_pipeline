"""Tests for stages/standardize_setting.py on small synthetic structures.

Needs a Python with cctbx (phenix.python, or a qFit conda env):
    CCTBX_PYTHON=phenix.python tests/run_tests.sh
or directly:
    phenix.python tests/test_standardize_setting.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

import iotbx.mtz
import iotbx.pdb
from cctbx import miller, sgtbx
from scitbx.array_family import flex

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "stages", "standardize_setting.py")
FAILS = []

# Two residues with an ANISOU record, a water and an SSBOND to itself (1555).
MODEL = """\
HEADER    TEST
SSBOND   1 CYS A    1    CYS A    2                          1555   1555  2.05
CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 {sg:<11s} 4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.033333  0.000000  0.000000        0.00000
SCALE2      0.000000  0.025000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.020000        0.00000
ATOM      1  N   CYS A   1       3.000   4.000   5.000  1.00 20.00           N
ANISOU    1  N   CYS A   1     2000   2500   3000    300   -200    100       N
ATOM      2  CA  CYS A   1       4.200   4.500   5.800  1.00 21.00           C
ATOM      3  SG  CYS A   1       5.100   6.000   7.200  1.00 22.00           S
ATOM      4  N   CYS A   2       7.000   9.000  11.000  0.60 23.00           N
ATOM      5  CA  CYS A   2       8.100   9.700  12.100  0.60 24.00           C
HETATM    6  O   HOH A 101      12.000  14.000  16.000  1.00 30.00           O
END
"""


def check(name, cond):
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")
    if not cond:
        FAILS.append(name)


def make_inputs(d, sg_symbol, with_phases=False):
    """Model in setting sg_symbol plus an mtz of FOBS/SIGFOBS/R-free-flags computed from it."""
    pdb = os.path.join(d, "x.pdb")
    mtz = os.path.join(d, "x.mtz")
    open(pdb, "w").write(MODEL.format(sg=sg_symbol))
    xrs = iotbx.pdb.input(pdb).xray_structure_simple()
    fc = xrs.structure_factors(d_min=2.0, algorithm="direct").f_calc()
    fobs = abs(fc).customized_copy(sigmas=flex.double(fc.size(), 1.0))
    flags = fobs.generate_r_free_flags(fraction=0.1)
    ds = fobs.as_mtz_dataset(column_root_label="FOBS")
    ds.add_miller_array(flags.customized_copy(data=flags.data().as_int()), column_root_label="R-free-flags")
    if with_phases:
        ds.add_miller_array(fc, column_root_label="2FOFCWT")
    ds.mtz_object().write(mtz)
    return pdb, mtz, xrs, fobs, flags


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


def arrays(mtz):
    return {a.info().label_string(): a for a in iotbx.mtz.object(mtz).as_miller_arrays()}


def test_nonstandard():
    print("non-standard setting (P 2 21 21 -> P 21 21 2)")
    d = tempfile.mkdtemp()
    pdb, mtz, xrs_old, fobs_old, flags_old = make_inputs(d, "P 2 21 21")
    r = run(pdb, mtz)
    check("exit 0", r.returncode == 0)
    check("originals kept", os.path.exists(os.path.join(d, "x.deposited_setting.pdb"))
          and os.path.exists(os.path.join(d, "x.deposited_setting.mtz")))
    check("deposited copy unchanged", open(os.path.join(d, "x.deposited_setting.pdb")).read()
          == MODEL.format(sg="P 2 21 21"))

    inp = iotbx.pdb.input(pdb)
    cb = sgtbx.space_group_info("P 2 21 21").change_of_basis_op_to_reference_setting()
    check("model now P 21 21 2", str(inp.crystal_symmetry().space_group_info()) == "P 21 21 2")
    check("cell permuted", tuple(round(x, 3) for x in inp.crystal_symmetry().unit_cell().parameters()[:3])
          == (40.0, 50.0, 30.0))
    text = open(pdb).read()
    check("SCALE/ORIGX dropped", "SCALE1" not in text and "ORIGX1" not in text)
    check("other records kept", text.startswith("HEADER    TEST") and "SSBOND   1 CYS A" in text)
    check("atom count kept", inp.atoms().size() == 6)
    check("occupancies kept", list(inp.atoms().extract_occ()) == [1, 1, 1, 0.6, 0.6, 1])

    # Same structure: complex F of the new model equals the old F, re-indexed.
    xrs_new = inp.xray_structure_simple()
    f_old = xrs_old.structure_factors(d_min=2.0, algorithm="direct").f_calc().change_basis(cb)
    f_old = f_old.customized_copy(crystal_symmetry=xrs_new.crystal_symmetry())
    f_new = xrs_new.structure_factors(d_min=2.0, algorithm="direct").f_calc()
    a, b = f_old.common_sets(f_new)
    rel = flex.max(flex.abs(a.data() - b.data())) / flex.mean(flex.abs(a.data()))
    check(f"model gives the same structure factors (max rel. diff {rel:.1e})", a.size() >= 0.99 * f_new.size() and rel < 1e-3)

    # Data: every column re-indexed, R-free flags follow their reflections.
    arr = arrays(mtz)
    check("same columns", iotbx.mtz.object(mtz).column_labels() == ["H", "K", "L", "FOBS", "SIGFOBS", "R-free-flags"])
    check("data now P 21 21 2", str(arr["FOBS,SIGFOBS"].space_group_info()) == "P 21 21 2")
    exp_f = fobs_old.change_basis(cb).customized_copy(crystal_symmetry=arr["FOBS,SIGFOBS"].crystal_symmetry())
    m1, m2 = exp_f.common_sets(arr["FOBS,SIGFOBS"])
    check("amplitudes follow their reflections", m1.size() == exp_f.size() and flex.max(flex.abs(m1.data() - m2.data())) < 1e-3)
    exp_flags = flags_old.change_basis(cb).customized_copy(crystal_symmetry=arr["R-free-flags"].crystal_symmetry())
    m1, m2 = exp_flags.common_sets(arr["R-free-flags"])
    check("R-free flags follow their reflections", m1.size() == exp_flags.size()
          and (m1.data() == (m2.data() != 0)).all_eq(True))
    check("free-set size unchanged", (arr["R-free-flags"].data() != 0).count(True) == flags_old.data().count(True))

    # Model and data agree in the new setting.
    fm = abs(f_new).customized_copy(crystal_symmetry=arr["FOBS,SIGFOBS"].crystal_symmetry())
    m1, m2 = fm.common_sets(arr["FOBS,SIGFOBS"])
    check("model matches re-indexed data", flex.max(flex.abs(m1.data() - m2.data())) / flex.mean(m2.data()) < 1e-3)

    before = (open(pdb).read(), open(mtz, "rb").read())
    r = run(pdb, mtz)
    check("second run is a no-op", r.returncode == 0 and "already the standard setting" in r.stdout
          and (open(pdb).read(), open(mtz, "rb").read()) == before)
    shutil.rmtree(d)


def test_standard_untouched():
    print("standard setting (P 21 21 21)")
    d = tempfile.mkdtemp()
    pdb, mtz, *_ = make_inputs(d, "P 21 21 21")
    before = (open(pdb).read(), open(mtz, "rb").read())
    r = run(pdb, mtz)
    check("exit 0, reported standard", r.returncode == 0 and "already the standard setting" in r.stdout)
    check("files unchanged", (open(pdb).read(), open(mtz, "rb").read()) == before)
    check("no copies made", not os.path.exists(os.path.join(d, "x.deposited_setting.pdb")))
    shutil.rmtree(d)


def test_refusals():
    print("refusals")
    d = tempfile.mkdtemp()
    pdb, mtz, *_ = make_inputs(d, "P 2 21 21", with_phases=True)
    before = (open(pdb).read(), open(mtz, "rb").read())
    r = run(pdb, mtz)
    check("phase columns refused", r.returncode != 0 and "phase columns" in r.stderr)
    check("... and nothing written", (open(pdb).read(), open(mtz, "rb").read()) == before)
    shutil.rmtree(d)

    d = tempfile.mkdtemp()
    pdb, mtz, *_ = make_inputs(d, "P 2 21 21")
    open(pdb, "w").write(MODEL.format(sg="P 2 21 21").replace("1555   1555", "1555   2555"))
    r = run(pdb, mtz)
    check("SSBOND to a symmetry mate refused", r.returncode != 0 and "symmetry mate" in r.stderr)
    shutil.rmtree(d)

    d = tempfile.mkdtemp()
    pdb, mtz, *_ = make_inputs(d, "P 21 21 21")
    open(pdb, "w").write(MODEL.format(sg="P 2 21 21"))
    r = run(pdb, mtz)
    check("model/data setting mismatch refused", r.returncode != 0 and "different space-group settings" in r.stderr)
    shutil.rmtree(d)


if __name__ == "__main__":
    test_nonstandard()
    test_standard_untouched()
    test_refusals()
    print(f"{'FAILED: ' + ', '.join(FAILS) if FAILS else 'all standardize_setting tests passed'}")
    sys.exit(1 if FAILS else 0)
