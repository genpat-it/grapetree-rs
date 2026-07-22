# Benchmarks — grapetree-rs vs Python GrapeTree

Host: 80-core Linux (AlmaLinux). Both at `n_proc=5`. Synthetic clonal datasets
(`gen_synth.py`), 1500 loci, 2% mutation, 5% duplicates. Wall-clock seconds
(single run), lower is better.

## distance method (symmetric PHYLIP)

| N     | Python | Rust  | speedup |
|-------|--------|-------|---------|
| 100   | 0.40   | 0.03  | 13.3×   |
| 300   | 0.59   | 0.12  | 4.9×    |
| 600   | 0.96   | 0.29  | 3.3×    |
| 1000  | 1.64   | 0.62  | 2.6×    |

Rust wins across the board; the rayon-parallel bit-loop dominates. (Python also
parallelises the distance step, so the gap narrows at large N but stays >2×.)

## MSTreeV2 (default method)

With the efficient O(E log V) Chu-Liu/Edmonds (lazy skew heap + rollback
union-find, KACTL-style), which replaced a slow recursive/HashMap version:

| N     | Python | Rust  | speedup |
|-------|--------|-------|---------|
| 100   | 0.78   | 0.08  | 9.8×    |
| 300   | 1.11   | 0.16  | 6.9×    |
| 600   | 1.74   | 0.59  | 2.9×    |
| 1000  | 3.14   | 1.43  | 2.2×    |
| 2000  | 9.71   | 7.54  | 1.3×    |

Rust wins at every size (largest at the small/medium N of interactive use).
At N=1000 the efficient Edmonds cut MSTreeV2 from 47.2s → 1.43s (**33×**). Both
implementations still trend together at large N because Python offloads the
arborescence to a compiled C binary and the distance step to multiple processes.

### Edmonds rewrite

`src/edmonds.rs` now uses the Tarjan/Gabow formulation: a lazy-add skew heap of
incoming edges per node and a rollback union-find, O(E log V) with reconstruction
via cycle expansion. It stays fuzz-verified 40/40 against the reference `edmonds`
binary, and MSTreeV2 output remains byte-identical to the oracle. No external
binary, no temp-file round trip.

## Scaling to large datasets (10k–80k)

Two things dominate at scale, and they are the *same* for grapetree-rs and the
reference:

1. **The distance matrix is dense O(n²·L).** Both build the full `n × n` matrix.
   The reference's own `-c` estimate is ~`103·n²` bytes (≈660 GB at n=80k) — this
   is server-scale work by design (GrapeTree runs inside EnteroBase).
2. **`get_shortcut` decides how big the Edmonds step is.** Before the minimum
   arborescence, near-identical "leaf" isolates within a distance cutoff
   (2/5/10/20, growing with N) are absorbed. On **clonal** data (real cgMLST,
   isolates a few alleles apart) this collapses the vast majority of nodes, so
   the arborescence runs on a small backbone. On data with no close pairs the
   shortcut does nothing and Edmonds runs on all n nodes — slow for *both* tools.

### `distance` scaling (grapetree-rs, 40 threads, 1500 loci)

| N     | wall     | peak RAM |
|-------|----------|----------|
| 10000 | 16 s     | 1.7 GB   |
| 20000 | 57 s     | 6.9 GB   |
| 40000 | 3 m 30 s | 24 GB    |

### MSTreeV2 on realistic clonal data — head to head (N=10000)

Single founder, isolates 1–4 alleles apart (a clonal expansion, like surveillance
cgMLST); 40 threads:

| tool                | wall    | peak RAM | tree            |
|---------------------|---------|----------|-----------------|
| GrapeTree (Python)  | 30.8 s  | 2.1 GB   | —               |
| **grapetree-rs**    | **9.4 s** | **1.7 GB** | RF≈0.0002 vs oracle, total length 24807 vs 24808 |

≈ **3.3× faster**, less memory, topologically equivalent. On the same-size but
adversarial "several maximally-distant clusters" input, neither tool finishes
quickly (shortcuts never fire) — so scalability comes from the data being clonal,
which is the regime GrapeTree targets.
