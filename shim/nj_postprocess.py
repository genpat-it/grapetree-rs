#!/usr/bin/env python
"""ete3 post-processing for the NJ family, byte-identical to GrapeTree.

Reproduces methods.NJ/RapidNJ/ninja's ete3 tail + backend() post-processing:
midpoint-root, unroot, rename leaves to representative names, collapse tiny
branches when the tree spans >3, expand `embeded` duplicate groups, and emit
NEWICK format=1. Reads a tree file + names/embeded JSON; prints the final tree.

Args: <tree.nwk> <names.json> <embeded.json> [scale]
  scale: multiply every branch length by this first (ninja divides dist by
  n_loci then multiplies back; pass n_loci). Default 1.
"""
import sys, json
from ete3 import Tree

def main():
    nwk, names_f, emb_f = sys.argv[1], sys.argv[2], sys.argv[3]
    scale = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    names = json.load(open(names_f))
    embeded = json.load(open(emb_f))
    tree = Tree(nwk)
    if scale != 1.0:
        for node in tree.traverse():
            node.dist *= scale
    try:
        tree.set_outgroup(tree.get_midpoint_outgroup())
        tree.unroot()
    except Exception:
        pass
    for leaf in tree.get_leaves():
        leaf.name = names[int(leaf.name.strip("'"))]
    # backend() post-processing
    maxDist = 0.0
    for node in tree.iter_descendants():
        if node.dist > maxDist:
            maxDist = node.dist
    if maxDist > 3:
        for node in tree.iter_descendants('postorder'):
            if 0 < node.dist < 0.1:
                for s in node.get_sisters():
                    s.dist += node.dist
                node.dist = 0
    for leaf in tree.get_leaves():
        group = embeded[leaf.name]
        if len(group) > 1:
            leaf.name = ''
            for n in group:
                leaf.add_child(name=n, dist=0.0)
    sys.stdout.write(tree.write(format=1).replace("'", ""))

if __name__ == "__main__":
    main()
