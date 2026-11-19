"""Move a model and its data to the reference (standard) space-group setting.

Usage: phenix.python standardize_setting.py MODEL.pdb DATA.mtz

qFit rebuilds space groups from their number when it copies a unit cell
(qfit.xtal.unitcell.UnitCell.copy), so a structure deposited in a non-standard
setting (e.g. P 2 21 21, No. 18) gets the operators of the standard setting
(P 21 21 2) during map scaling, and symmetry mates end up in the wrong place.
Re-expressing model and data in the standard setting avoids this without
changing qFit. It is a change of axes only: the structure, the reflections and
their R-free flags are unchanged.

If the setting is already standard, nothing is written. Otherwise the original
files are kept as <name>.deposited_setting.pdb / .mtz and MODEL / DATA are
overwritten with the transformed versions. Running it again is a no-op.

Model: CRYST1, ATOM/HETATM coordinates and ANISOU records are transformed;
SCALEn and ORIGXn records are dropped (they would be stale; Phenix does not
need them); every other record is kept as is. Data: all columns, including
the R-free flags, are re-indexed with iotbx.mtz change_basis_in_place.
"""
import os
import shutil
import sys

import iotbx.mtz
import iotbx.pdb
from scitbx import matrix


def die(msg):
    sys.exit(f"ERROR: standardize_setting: {msg}")


def cartesian_transform(cb_op, uc_old, uc_new):
    """Rotation M and translation t with x_new = M x_old + t (Cartesian, Angstrom)."""
    c = cb_op.c()
    r = matrix.sqr(c.r().as_double())
    t = matrix.col(c.t().as_double())
    o_new = matrix.sqr(uc_new.orthogonalization_matrix())
    f_old = matrix.sqr(uc_old.fractionalization_matrix())
    return o_new * r * f_old, o_new * t


def transform_pdb_lines(lines, m, t, cryst1):
    out = []
    for line in lines:
        rec = line[:6]
        if rec.startswith(("SCALE", "ORIGX")):
            continue
        if rec == "CRYST1":
            out.append(cryst1)
        elif rec in ("ATOM  ", "HETATM"):
            x = matrix.col([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            xn = m * x + t
            out.append(f"{line[:30]}{xn[0]:8.3f}{xn[1]:8.3f}{xn[2]:8.3f}{line[54:]}")
        elif rec == "ANISOU":
            u = [int(line[28 + 7 * i: 35 + 7 * i]) for i in range(6)]  # U11 U22 U33 U12 U13 U23 (x1e4)
            u_mat = matrix.sqr([u[0], u[3], u[4], u[3], u[1], u[5], u[4], u[5], u[2]])
            un = m * u_mat * m.transpose()
            vals = [un[0], un[4], un[8], un[1], un[2], un[5]]
            out.append(line[:28] + "".join(f"{int(round(v)):7d}" for v in vals) + line[70:])
        else:
            out.append(line)
    return out


def main():
    if len(sys.argv) != 3:
        die("usage: standardize_setting.py MODEL.pdb DATA.mtz")
    pdb_path, mtz_path = sys.argv[1:3]

    mtz = iotbx.mtz.object(mtz_path)
    sg_mtz = mtz.space_group_info()
    phases = [c.label() for c in mtz.columns() if c.type() == "P"]
    cs_pdb = iotbx.pdb.input(pdb_path).crystal_symmetry()
    if cs_pdb is None or cs_pdb.space_group_info() is None:
        die(f"{pdb_path} has no usable CRYST1 record.")
    sg_pdb = cs_pdb.space_group_info()
    if sg_pdb.group() != sg_mtz.group():
        die(f"model ({sg_pdb}) and data ({sg_mtz}) are in different space-group settings.")

    cb_op = sg_mtz.change_of_basis_op_to_reference_setting()
    if cb_op.is_identity_op():
        print(f"standardize_setting: {sg_mtz} is already the standard setting; files unchanged.")
        return
    if phases:
        die(f"{mtz_path} has phase columns ({', '.join(phases)}), which cannot be re-indexed "
            f"in place. Run this on the raw data, before any maps are made.")
    sg_new = sg_mtz.change_basis(cb_op)
    if not sg_new.change_of_basis_op_to_reference_setting().is_identity_op():
        die(f"change of basis {cb_op} did not give a standard setting ({sg_new}).")

    lines = open(pdb_path).read().splitlines()
    for line in lines:
        if line.startswith(("LINK  ", "SSBOND")) and any(
                op.strip() not in ("", "1555") for op in (line[59:65], line[66:72])):
            die(f"{pdb_path} has a LINK/SSBOND to a symmetry mate, which this script "
                f"cannot re-express: {line.strip()}")
        if line.startswith("MTRIX") and not line.rstrip().endswith("1"):
            die(f"{pdb_path} has non-given MTRIX records, which this script cannot re-express.")

    uc_old = cs_pdb.unit_cell()
    uc_new = uc_old.change_basis(cb_op)
    m, t = cartesian_transform(cb_op, uc_old, uc_new)
    if abs(m.determinant() - 1) > 1e-6 or not m.is_r3_rotation_matrix(1e-4):
        die(f"change of basis {cb_op} is not a pure rotation in Cartesian space.")
    old_cryst1 = next(l for l in lines if l.startswith("CRYST1"))
    z = old_cryst1[66:70].ljust(4)
    a, b, c, al, be, ga = uc_new.parameters()
    cryst1 = f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{al:7.2f}{be:7.2f}{ga:7.2f} {sg_new.type().lookup_symbol():<11s}{z}"

    new_lines = transform_pdb_lines(lines, m, t, cryst1)
    mtz.change_basis_in_place(cb_op=cb_op, new_space_group_info=sg_new)

    for path in (pdb_path, mtz_path):
        root, ext = os.path.splitext(path)
        kept = f"{root}.deposited_setting{ext}"
        if not os.path.exists(kept):
            shutil.copy2(path, kept)
    with open(pdb_path, "w") as f:
        f.write("\n".join(new_lines) + "\n")
    mtz.write(mtz_path)
    print(f"standardize_setting: {sg_mtz} -> {sg_new} (change of basis {cb_op}); "
          f"unit cell {uc_old.parameters()} -> {uc_new.parameters()}. "
          f"Originals kept as *.deposited_setting.pdb/.mtz.")


if __name__ == "__main__":
    main()
