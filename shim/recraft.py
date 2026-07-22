#!/usr/bin/env python
"""branch_recraft via NumPy — the second bit grapetree-rs cannot reproduce.

GrapeTree's `contemporary` likelihood test uses `np.log`, and NumPy's float64
`log` is its own polynomial (not correctly rounded, ~1 ULP off ~21% of the time,
and CPU/SIMD-dispatched). Rust's `f64::ln` uses the platform libm and differs, so
near-tied `p1 >= p2` decisions flip. "What is external stays external": for the
bit-identical (default) mode we run branch_recraft here in NumPy — a verbatim port
of the reference `methods._branch_recraft` / `contemporary` — so `np.log` is the
same on the same environment. Requires the same NumPy the reference used.

Usage: recraft.py <net.tsv> <dist.f32.bin> <N> <weight.f32.bin> <n_loci> <out.tsv>
  net.tsv  : network edges, one "src\\ttgt\\tbrlen" per line (ints, float)
  dist.bin : row-major N*N float32 asymmetric distance matrix
  weight.bin: N float32 harmonic weights
  out.tsv  : recrafted edges written here, "src\\ttgt\\tbrlen" per line
"""
import sys
import numpy as np


def contemporary(a, b, c, n_loci):
    a[0], a[1] = max(min(a[0], n_loci - 0.5), 0.5), max(min(a[1], n_loci - 0.5), 0.5)
    b, c = max(min(b, n_loci - 0.5), 0.5), max(min(c, n_loci - 0.5), 0.5)
    if b >= a[0] + c and b >= a[1] + c:
        return False
    elif b == c:
        return True
    s11, s12 = np.sqrt(1 - a[0] / n_loci), (2 * n_loci - b - c) / 2 / np.sqrt(n_loci * (n_loci - a[0]))
    v = 1 - ((n_loci - a[1]) * (n_loci - c) / n_loci + (n_loci - b)) / 2 / n_loci
    s21, s22 = 1 + a[1] * v / (b - 2 * n_loci * v), 1 + c * v / (b - 2 * n_loci * v)
    p1 = a[0] * np.log(1 - s11 * s11) + (n_loci - a[0]) * np.log(s11 * s11) + (b + c) * np.log(1 - s11 * s12) + (2 * n_loci - b - c) * np.log(s11 * s12)
    p2 = a[1] * np.log(1 - s21) + (n_loci - a[1]) * np.log(s21) + b * np.log(1 - s21 * s22) + (n_loci - b) * np.log(s21 * s22) + c * np.log(1 - s22) + (n_loci - c) * np.log(s22)
    return p1 >= p2


def branch_recraft(branches, dist, weights, n_loci):
    group_id = {b: b for br in branches for b in br[:2]}
    groups = {b: [b] for br in branches for b in br[:2]}
    childrens = {b: [] for br in branches for b in br[:2]}
    branches = sorted(branches, key=lambda br: [dist[br[0], br[1]]] + sorted([weights[br[0]], weights[br[1]]]))
    i = 0
    while i < len(branches):
        src, tgt, brlen = branches[i]
        sources, targets = groups[group_id[src]], groups[group_id[tgt]]
        tried = {}
        if len(sources) > 1:
            for w, d, s in sorted(zip(weights[sources], dist[sources, tgt], sources))[:3]:
                if s == src:
                    break
                if d < 1.5 * dist[src, tgt]:
                    if contemporary([dist[s, src], dist[src, s]], d, dist[src, tgt], n_loci):
                        tried[src], src = s, s
                        break
            while src not in tried:
                tried[src] = src
                mid_nodes = sorted([[weights[s], dist[s, tgt], s] for s in childrens[src] if s not in tried and dist[s, tgt] < 2 * dist[src, tgt]])
                for w, d, s in mid_nodes:
                    if d < dist[src, tgt]:
                        if not contemporary([dist[src, s], dist[s, src]], dist[src, tgt], d, n_loci):
                            tried[src], src = s, s
                            break
                    elif w < weights[src]:
                        if contemporary([dist[s, src], dist[src, s]], d, dist[src, tgt], n_loci):
                            tried[src], src = s, s
                            break
                    tried[s] = src
        if len(targets) > 1:
            for w, d, t in sorted(zip(weights[targets], dist[src, targets], targets))[:3]:
                if t == tgt:
                    break
                if d < 1.5 * dist[src, tgt]:
                    if contemporary([dist[t, tgt], dist[tgt, t]], d, dist[src, tgt], n_loci):
                        tried[tgt], tgt = t, t
                        break
            while tgt not in tried:
                tried[tgt] = tgt
                mid_nodes = sorted([[weights[t], dist[src, t], t] for t in childrens[tgt] if t not in tried and dist[src, t] < 2 * dist[src, tgt]])
                for w, d, s in mid_nodes:
                    if d < dist[src, tgt]:
                        if not contemporary([dist[tgt, t], dist[t, tgt]], dist[src, tgt], d, n_loci):
                            tried[tgt], tgt = t, t
                            break
                    elif w < weights[tgt]:
                        if contemporary([dist[t, tgt], dist[tgt, t]], d, dist[src, tgt], n_loci):
                            tried[tgt], tgt = t, t
                            break
                    tried[t] = tgt
        brlen = dist[src, tgt]
        branches[i] = [src, tgt, brlen]
        if i >= len(branches) - 1 or branches[i + 1][2] >= brlen:
            tid = group_id[tgt]
            for t in targets:
                group_id[t] = group_id[src]
            groups[group_id[src]].extend(groups.pop(tid, []))
            childrens[src].append(tgt)
            childrens[tgt].append(src)
            i += 1
        else:
            branches[i:] = sorted(branches[i:], key=lambda br: br[2])
    return branches


def main():
    net_p, dist_p, n_s, w_p, nloci_s, out_p = sys.argv[1:7]
    n = int(n_s)
    n_loci = int(nloci_s)
    # memmap (not fromfile): branch_recraft touches only a subset of cells (the
    # local recraft neighbourhood), so the OS pages in far less than the full
    # 13 GB — same f32 values, so the result is bit-identical.
    dist = np.memmap(dist_p, dtype="<f4", mode="r", shape=(n, n))
    weights = np.fromfile(w_p, dtype="<f4")
    branches = []
    with open(net_p) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 3:
                branches.append([int(p[0]), int(p[1]), float(p[2])])
    out = branch_recraft(branches, dist, weights, n_loci)
    with open(out_p, "w") as fh:
        for s, t, d in out:
            fh.write(f"{int(s)}\t{int(t)}\t{float(d)}\n")


if __name__ == "__main__":
    main()
