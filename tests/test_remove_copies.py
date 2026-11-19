"""Tests for stages/remove_copies.py on small synthetic models.

Standard library only:  python3 tests/test_remove_copies.py
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "stages", "remove_copies.py")
FAILS = []


def atom(serial, name, alt, resn, chain, resi, xyz, occ, element, record="ATOM"):
    x, y, z = xyz
    return (f"{record:<6s}{serial:5d} {name:<4s}{alt}{resn:>3s} {chain}{resi:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{20.0:6.2f}          {element:>2s}")


def anisou(atom_line):
    return "ANISOU" + atom_line[6:27] + "  2000   2500   3000    300   -200    100" + atom_line[70:]


def residue(serial, resn, chain, resi, confs):
    """confs: list of (altloc, occupancy, shift); atoms N, CA, CB, HA."""
    lines = []
    for alt, occ, dx in confs:
        for name, el, base in (("N", "N", (0, 0, 0)), ("CA", "C", (1.5, 0, 0)),
                               ("CB", "C", (2.0, 1.4, 0)), ("HA", "H", (1.9, -0.9, 0.5))):
            line = atom(serial, name, alt, resn, chain, resi, (base[0] + dx, base[1] + resi, base[2]), occ, el)
            lines += [line, anisou(line)]
            serial += 1
    return lines, serial


def build():
    lines, s = ["CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1"], 1
    blocks = [
        ("ALA", "A", 1, [(" ", 1.00, 0.0)]),                                      # single
        ("SER", "A", 2, [("A", 0.34, 0.0), ("B", 0.21, 0.8), ("C", 0.45, 0.8)]),  # B copies C
        ("VAL", "A", 3, [("A", 0.30, 0.0), ("B", 0.30, 0.0), ("C", 0.40, 0.0)]),  # all copies
        ("LEU", "A", 4, [("A", 0.50, 0.0), ("B", 0.50, 0.005)]),                 # 0.005 A: copy
        ("ILE", "A", 5, [("A", 0.50, 0.0), ("B", 0.50, 0.05)]),                  # 0.05 A: distinct
        ("SER", "B", 2, [("A", 0.60, 0.0), ("B", 0.40, 1.0)]),                   # distinct
    ]
    for resn, chain, resi, confs in blocks:
        block, s = residue(s, resn, chain, resi, confs)
        lines += block
    lines += ["TER", "END"]
    return lines


def atoms(path):
    out = {}
    for line in open(path):
        if line.startswith("ATOM"):
            out[(line[21], int(line[22:26]), line[12:16].strip(), line[16])] = (line[30:54], float(line[54:60]))
    return out


def check(name, cond):
    print(f"  {'ok  ' if cond else 'FAIL'} {name}")
    if not cond:
        FAILS.append(name)


def main():
    d = tempfile.mkdtemp()
    src, dst = os.path.join(d, "in.pdb"), os.path.join(d, "out.pdb")
    open(src, "w").write("\n".join(build()) + "\n")
    r = subprocess.run([sys.executable, SCRIPT, src, dst], capture_output=True, text=True)
    check("exit 0", r.returncode == 0)
    check("summary line", "merged 4 exact-copy conformers in 3 residues (2 now single-conformer)" in r.stdout)
    a, b = atoms(src), atoms(dst)

    check("single-conformer residue unchanged", all(b[k] == a[k] for k in a if k[1] == 1))
    check("copy B of A2 removed, hydrogens included", not any(k[:2] == ("A", 2) and k[3] == "B" for k in b))
    check("its occupancy added to C (0.45 + 0.21)", b[("A", 2, "CA", "C")][1] == 0.66 and b[("A", 2, "HA", "C")][1] == 0.66)
    check("A2 conformer A untouched", b[("A", 2, "CA", "A")] == a[("A", 2, "CA", "A")])
    check("altloc letters not relabelled (A, C)", sorted({k[3] for k in b if k[:2] == ("A", 2)}) == ["A", "C"])
    check("all-copy residue becomes single, blank altloc, occupancy 1",
          sorted(k[3] for k in b if k[:2] == ("A", 3)) == [" "] * 4 and b[("A", 3, "CA", " ")][1] == 1.0)
    check("... keeps the highest-occupancy conformer's coordinates", b[("A", 3, "CA", " ")][0] == a[("A", 3, "CA", "C")][0])
    check("0.005 A pair merged", ("A", 4, "CA", "B") not in b and b[("A", 4, "CA", " ")][1] == 1.0)
    check("0.05 A pair kept (above 0.01 A cutoff)", b[("A", 5, "CA", "B")] == a[("A", 5, "CA", "B")])
    check("distinct conformers in chain B kept", all(b[k] == a[k] for k in a if k[0] == "B"))

    text = open(dst).read().splitlines()
    check("one ANISOU per remaining atom, right after it",
          all(text[i + 1].startswith("ANISOU") and text[i + 1][6:27] == text[i][6:27]
              for i in range(len(text)) if text[i].startswith("ATOM")))
    check("no orphan ANISOU", sum(l.startswith("ANISOU") for l in text) == sum(l.startswith("ATOM") for l in text))
    check("other records kept", text[0].startswith("CRYST1") and text[-2:] == ["TER", "END"])

    r = subprocess.run([sys.executable, SCRIPT, dst, os.path.join(d, "again.pdb")], capture_output=True, text=True)
    check("second pass finds nothing", "merged 0 exact-copy conformers" in r.stdout
          and open(os.path.join(d, "again.pdb")).read() == open(dst).read())

    r = subprocess.run([sys.executable, SCRIPT, src, dst, "--cutoff", "0.1"], capture_output=True, text=True)
    check("--cutoff 0.1 also merges the 0.05 A pair", "merged 5 exact-copy conformers in 4 residues" in r.stdout)

    print("all remove_copies tests passed" if not FAILS else "FAILED: " + ", ".join(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
