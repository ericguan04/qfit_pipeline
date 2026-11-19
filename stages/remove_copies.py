"""Merge alternate conformers that are exact copies of another conformer.

Usage: python remove_copies.py IN.pdb OUT.pdb [--cutoff 0.01]

qFit's segment step gives every residue in a segment the same set of altlocs,
padding residues that have fewer real conformers with exact copies (its
identical-conformer removal is commented out in QFitSegment). This script
undoes that padding before final refinement.

For each residue (chain, number, insertion code, name), altlocs are visited
from highest to lowest occupancy (ties: altloc letter). An altloc is a copy of
an earlier, kept altloc if both have the same heavy atoms and their heavy-atom
RMSD is below the cutoff (default 0.01 A, qFit's own --rmsd-cutoff). A copy is
removed (with its hydrogens and ANISOU records) and its occupancy is added to
the conformer it copies. A residue left with one conformer becomes a normal
single-conformer residue (blank altloc, occupancy = the summed occupancy).

Atoms without an altloc, other altlocs, every other record and the order of
the file are unchanged; remaining altloc letters are not relabelled. Prints
one summary line per run. Standard library only.
"""
import argparse
import math
import sys
from collections import OrderedDict, defaultdict


def is_atom(line):
    return line.startswith(("ATOM  ", "HETATM"))


def residue_key(line):
    return (line[21], line[22:26], line[26], line[17:20])


def is_hydrogen(line):
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element in ("H", "D")
    return line[12:16].strip().lstrip("0123456789").startswith(("H", "D"))


def rmsd(a, b):
    s = sum((p - q) ** 2 for n in a for p, q in zip(a[n], b[n]))
    return math.sqrt(s / len(a))


def find_copies(lines, cutoff):
    """Returns {line index: new occupancy}, the set of line indices to drop,
    the indices to make single-conformer (blank altloc), and counts."""
    # residue -> altloc -> {"idx": [line indices], "heavy": {name: xyz}, "occ": max occupancy}
    residues = defaultdict(OrderedDict)
    for i, line in enumerate(lines):
        if not is_atom(line) or line[16] == " ":
            continue
        alt = residues[residue_key(line)].setdefault(line[16], {"idx": [], "heavy": {}, "occ": 0.0})
        alt["idx"].append(i)
        alt["occ"] = max(alt["occ"], float(line[54:60]))
        if not is_hydrogen(line):
            alt["heavy"][line[12:16]] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    new_occ, drop, make_single = {}, set(), set()
    n_residues = n_removed = n_single = 0
    for alts in residues.values():
        if len(alts) < 2:
            continue
        order = sorted(alts, key=lambda a: (-alts[a]["occ"], a))
        kept, total = [], {}
        for a in order:
            heavy = alts[a]["heavy"]
            target = next((k for k in kept if heavy and set(heavy) == set(alts[k]["heavy"])
                           and rmsd(heavy, alts[k]["heavy"]) < cutoff), None)
            if target is None:
                kept.append(a)
                total[a] = alts[a]["occ"]
            else:
                total[target] += alts[a]["occ"]
                drop.update(alts[a]["idx"])
                n_removed += 1
        if len(kept) == len(alts):
            continue
        n_residues += 1
        for k in kept:
            for i in alts[k]["idx"]:
                new_occ[i] = total[k]
        if len(kept) == 1:
            make_single.update(alts[kept[0]]["idx"])
            n_single += 1
    return new_occ, drop, make_single, (n_residues, n_removed, n_single)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("pdb_in")
    p.add_argument("pdb_out")
    p.add_argument("--cutoff", type=float, default=0.01, help="heavy-atom RMSD (A) below which two altlocs are copies")
    args = p.parse_args()

    lines = open(args.pdb_in).read().splitlines()
    new_occ, drop, make_single, (n_res, n_removed, n_single) = find_copies(lines, args.cutoff)

    out, dropped_serials = [], set()
    last_dropped = last_single = False
    for i, line in enumerate(lines):
        if line.startswith("ANISOU"):
            # An ANISOU record belongs to the ATOM/HETATM record just before it.
            if not last_dropped:
                out.append(line[:16] + " " + line[17:] if last_single else line)
            continue
        last_dropped, last_single = i in drop, i in make_single
        if last_dropped:
            dropped_serials.add(line[6:11].strip())
            continue
        if i in new_occ:
            line = line[:54] + f"{min(new_occ[i], 1.0):6.2f}" + line[60:]
        if last_single:
            line = line[:16] + " " + line[17:]
        if line.startswith("CONECT") and dropped_serials & set(line[6:].split()):
            continue
        out.append(line)

    with open(args.pdb_out, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"remove_copies: merged {n_removed} exact-copy conformers in {n_res} residues "
          f"({n_single} now single-conformer); cutoff {args.cutoff} A")


if __name__ == "__main__":
    sys.exit(main())
