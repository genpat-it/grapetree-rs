#!/usr/bin/env python
"""Harmonic weights via NumPy — the one bit grapetree-rs cannot reproduce.

GrapeTree ranks nodes by harmonic centrality `N / sum(1/(dist+0.1))` computed in
float32. NumPy's float32 `sum` uses a SIMD (AVX) reduction whose addition order
depends on the CPU/build; floating-point addition is not associative, so the
last bit is not portably reproducible in Rust. On large inputs this flips ~1 ULP
on many nodes, which flips near-tied harmonic *ranks*, which perturbs the reduced
arborescence matrix and the final tree (~0.8% of leaves at 63k).

"What is external stays external": for the bit-identical (default) mode we hand
the distance matrix to NumPy and let it compute the exact reference weights —
just like we delegate the arborescence to the `edmonds` binary. Requires the
same NumPy the reference used (set GT_PYTHON to that env).

Usage: harmonic_weights.py <dist.f32.bin> <N> <n_str.txt> <out.f32.bin>
  dist.f32.bin : row-major N*N float32 asymmetric distance matrix
  n_str.txt    : N integers (one per line), embeded group sizes per node
  out.f32.bin  : N float32 weights written here
"""
import sys
import numpy as np


def main():
    dist_path, n_s, nstr_path, out_path = sys.argv[1:5]
    n = int(n_s)
    dist = np.fromfile(dist_path, dtype="<f4").reshape(n, n)
    n_str = np.loadtxt(nstr_path, dtype=int).reshape(-1)
    assert n_str.size == n, f"n_str size {n_str.size} != N {n}"

    # EXACT reference `distance_matrix.harmonic` (module/MSTrees.py):
    weights = dist.shape[0] / np.sum(1.0 / (dist + 0.1), 1)
    cw = np.vstack([-np.asarray(n_str), weights])
    weights[np.lexsort(cw)] = np.arange(dist.shape[0], dtype=float) / dist.shape[0]

    weights.astype("<f4").tofile(out_path)


if __name__ == "__main__":
    main()
