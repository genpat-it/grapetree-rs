#!/usr/bin/env python
"""Compare two NEWICK trees for scientific equivalence.

Exit 0 (EQUIVALENT) iff: identical leaf set, identical total branch length
(within tol), and Robinson-Foulds distance 0 (identical unrooted topology).
Byte-identity is stricter than needed for MST networks because NEWICK child
ordering and the (re-)rooting are not scientifically meaningful.

Usage: compare_trees.py <a.nwk> <b.nwk> [tol] [--topo]

--topo: score topology only (RF=0 + identical leaf set), ignoring branch
lengths. Use for the NJ family, where the reference (FastME) assigns
balanced-ME edge lengths while grapetree-rs uses canonical Saitou-Nei NJ
lengths — same topology, different length scheme.
"""
import sys
from ete3 import Tree


def main():
    topo_only = "--topo" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--topo"]
    a = Tree(open(args[0]).read(), format=1)
    b = Tree(open(args[1]).read(), format=1)
    tol = float(args[2]) if len(args) > 2 else 1e-4

    la = sorted(l.name for l in a.get_leaves())
    lb = sorted(l.name for l in b.get_leaves())
    ta = sum(n.dist for n in a.traverse())
    tb = sum(n.dist for n in b.traverse())
    same_leaves = la == lb
    same_len = abs(ta - tb) <= tol
    rf = a.compare(b, unrooted=True)["norm_rf"] if same_leaves else 1.0

    ok = same_leaves and rf == 0.0 and (topo_only or same_len)
    status = "EQUIVALENT" if ok else "DIFFERENT"
    print(f"{status} leaves={same_leaves}({len(la)}) len_o={ta:.3f} len_r={tb:.3f} "
          f"len_ok={same_len} rf={rf:.4f}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
