# grapetree-rs — Design & Decisions

A from-scratch Rust port of **GrapeTree** (achtman-lab/GrapeTree, Zhou et al. 2018,
*Genome Research*). Goal: byte-comparable results to the reference Python
implementation, with state-of-the-art performance, and a regression test for
**every** CLI parameter driven by well-formed synthetic inputs.

## Scope

GrapeTree = computational backend (`module/MSTrees.py`, 810 LOC) + a Flask web
server + a D3.js visualisation front-end.

This port targets the **computational backend** — the part that turns an allelic
profile into a NEWICK tree or a PHYLIP distance matrix. That is the
performance-critical, scientifically load-bearing code. The Flask server and the
JS visualiser are UI shells around `backend()`; they are explicitly out of scope
for the Rust core (a thin server can wrap the binary later if wanted). This
decision is recorded so "port the whole project" is understood as "port the
engine that does the science, with full CLI parity".

## Reference behaviour (from MSTrees.py)

Pipeline of `backend(**args)`:
1. **Parse** profile (TSV) or FASTA. Header line starts with `#` (single `#`),
   `##` lines ignored. Column 0 = strain name. Columns whose header starts with
   `#` or equals `ST`/`ST_id` (case-insensitive) are dropped. Values uppercased.
   Names sanitised: `[()  ,"';]` -> `_`.
2. **nonredundant**: per-column integer-encode alleles by *sorted string order*
   (`np.unique(..., return_inverse=True)+1`); `0`/`N`/`-` -> 0 (missing). For
   `complete_delete`, keep only columns that contain a zero (replicated verbatim
   — matches reference, even though the name is counter-intuitive). Lexsort rows
   (last locus = primary key), drop all-missing rows, collapse identical rows into
   `embeded` groups keyed by the first row of each run.
3. **Distance matrix**: `symmetric` | `asymmetric` | `asymmetric_wgMLST` |
   `blockwise`, times missing handlers `pair_delete` | `complete_delete` |
   `as_allele` | `absolute_distance`. Stored `float32`.
4. **Weights (tiebreak heuristic)**: `harmonic` (MSTreeV2 default) or `eBurst`
   (MSTree default). Depends on `n_str` = size of each embeded group.
5. **Tree**:
   - `MSTree`/`MSTreeV2` -> `_symmetric` (Kruskal MST) or `_asymmetric`
     (`get_shortcut` collapse + minimum spanning **arborescence** / Edmonds).
     MSTreeV2 forces asymmetric + harmonic + `branch_recraft=True`.
   - `_branch_recraft`: local branch recrafting gated by `contemporary()`
     (a likelihood-ratio contemporaneity test).
   - `_network2tree`: orient edge list from a root, build tree, add zero-length
     child for each internal node's own name.
   - `NJ`/`RapidNJ`/`ninja`: reference shells out to FastME / RapidNJ / Ninja.
   - `distance`: PHYLIP square matrix, values divided by n_loci unless
     `absolute_distance`/`blockwise`.
6. **Post-process**: if max branch length > 3, collapse tiny (`0<dist<0.1`)
   branches into siblings; expand embeded groups back into multifurcations;
   emit NEWICK `format=1`, strip `'`.

## Parameter -> behaviour map (CLI parity target)

| CLI flag | dest | effect |
|---|---|---|
| `-p/--profile` | fname | input file (or inline text); `.gz` supported |
| `-m/--method` | MSTreeV2[def]/MSTree/NJ/RapidNJ/ninja/distance | tree algorithm |
| `-x/--matrix` | symmetric[def]/asymmetric/blockwise | matrix; MSTreeV2->asymmetric |
| `-r/--recraft` | branch_recraft | local recrafting (forced on for MSTreeV2) |
| `-y/--missing` | 0..3 -> pair_delete/complete_delete/as_allele/absolute_distance | missing handling |
| `-w/--wgMLST` | wgMLST | asymmetric_wgMLST matrix variant |
| `-t/--heuristic` | eBurst/harmonic | tiebreak; MSTreeV2->harmonic |
| `-n/--n_proc` | parallel workers | perf only (rayon), must not change results |
| `-c/--check` | checkEnv | print time/mem estimate JSON |
| `-b/--block_penalty` | 0.01 | blockwise penalty |

MSTreeV2 is an alias: method=MSTree, matrix=asymmetric, heuristic=harmonic,
recraft=True. Blockwise forces method=MSTree, recraft off, alleles all "real".

## Rust architecture

Flat `Vec` matrices (row-major), `u32` allele codes, `f32` distances to mirror
NumPy dtypes exactly. Modules:
- `parse`  — profile/FASTA reader + `nonredundant`
- `distance` — the four matrix kinds × four handlers (rayon-parallel, bit-tricks)
- `heuristic` — harmonic / eBurst weights (exact lexsort/bincount semantics)
- `mst` — `_symmetric` Kruskal(+union-find), `_asymmetric` shortcut+Edmonds
- `recraft` — `contemporary` + `_branch_recraft`
- `nj` — PHYLIP `distance`; native canonical NJ (replaces FastME/RapidNJ binaries)
- `tree` — tree struct + NEWICK writer (replaces ete3), post-processing
- `cli`  — clap args mirroring `add_args`

## Performance strategy (SOTA)

- Distance matrix is the hot path (O(n²·L)). Encode profiles into cache-friendly
  row-major `u32`; inner loop is branchless `a!=b` mask + popcount-style sum,
  auto-vectorised; rayon over the outer index range (matches Python's process
  chunking but shared-memory, zero serialisation).
- Kruskal with union-find (path compression + union by size) instead of
  NetworkX; Edmonds arborescence in native Rust (removes the `edmonds` binary
  fork + temp-file round-trip, the reference's biggest overhead for asymmetric).
- Avoid the reference's `np.save`/`np.load` temp-file dance entirely.

## Regression methodology

Oracle = reference Python, run via `scratchpad/oracle_run.py` in the `gt_oracle`
conda env (Python 3.11; ete3 needs <3.13 because `cgi` was removed). Golden
outputs are generated per (dataset × method × matrix × missing × heuristic ×
recraft × wgMLST) combination on synthetic inputs with controlled properties
(clonal expansions, missing-data gradients, duplicate rows, wgMLST presence
gradients, blockwise-ordered loci). Compare:
- `distance` matrices: exact float text match (same formatting).
- Trees: edge-set / split (Robinson–Foulds) equivalence, since NEWICK child
  ordering and dedup-representative choice are cosmetic. Where the reference is
  deterministic we also check exact topology + branch lengths.

## Status log

- 2026-07-22: Scaffolded project; reference cloned; oracle env working and
  producing golden trees + matrices for simulated_data (MSTreeV2/MSTree/distance).
- 2026-07-22: Parser + `nonredundant` byte-verified (row/loci counts match).
- 2026-07-22: Distance matrices done — all 4 kinds × 4 handlers, rayon-parallel;
  `distance` PHYLIP output 38/38 byte-identical across 5 datasets. Documented the
  upstream `--wgMLST` dead-code no-op and preserved parity.
- 2026-07-22: Heuristics (harmonic/eBurst) match oracle weight vectors.
- 2026-07-22: **Plain MSTree (`-x symmetric`) fully done — 30/30 byte-identical**
  across 5 datasets × {eBurst,harmonic} × {pair_delete,as_allele,absolute_distance},
  including topology, branch lengths AND newick child ordering. Key detail: MST
  edges must be emitted in `nx.Graph.edges()` adjacency order (not Kruskal order)
  to match ete3 child ordering.
- 2026-07-22: **Native Edmonds optimum branching** (Chu-Liu/Edmonds with cycle
  contraction + reconstruction) replaces the `edmonds` binary + temp-file round
  trip. Fuzz-verified 40/40 identical edge sets vs the reference binary on random
  distinct-weight matrices; the binary computes the global-min arborescence (min
  over all roots), reproduced via a virtual-root construction.
- 2026-07-22: **Asymmetric MSTree** (`get_shortcut` + Edmonds + network2tree)
  done — 30/30 topologically identical to the oracle (RF=0, equal total length,
  equal leaf set) across 5 datasets × {harmonic,eBurst} × 3 handlers. NOTE: the
  asymmetric path is scored by topological identity, not byte-identity: the
  reference's NEWICK child ordering comes from the compiled `edmonds` binary's
  arbitrary output order (an implementation artifact). Our Edmonds output is made
  deterministic (sorted); the resulting tree is the same tree, serialised in a
  canonical order.
- 2026-07-22: **MSTreeV2 complete** (`branch_recraft` + `contemporary`). The
  default method is 5/6 byte-identical to the oracle across the test datasets
  (simulated_data, miss_N, snp, clean, Agona). Key fidelity detail: the reference
  computes arborescence branch lengths through the `edmonds` binary's 5-decimal
  text round-trip (`%.5f` of `round(d)+weight+0.999995`, then `−1`); reproducing
  that quantisation is necessary because recraft re-sorts branches by `brlen` and
  the 5-dp values flip near-ties. The 6th dataset (miss_int) differs by a single
  branch (RF=0.0588, +3 total length): a `contemporary` likelihood-test decision
  that lands on a `p1 >= p2` boundary where Rust libm and NumPy differ by a ULP in
  `log`/`sqrt`. This is an inherent floating-point-boundary difference, not a
  logic error — the tree is otherwise identical. Regression policy for trees:
  byte-identical where achievable; else topological identity (RF=0) + equal total
  length; boundary cases documented.
- 2026-07-22: **Native Neighbor-Joining** (canonical Saitou-Nei) replaces the
  FastME/RapidNJ/Ninja binaries for `-m NJ/RapidNJ/ninja`. Topology matches the
  oracle (RF=0) on nearly all inputs. Two documented divergence classes: (a) NJ
  branch lengths use the Saitou-Nei scheme, whereas FastME assigns balanced-ME
  (Pauplin) edge lengths — same topology, different lengths, so the NJ family is
  scored topology-only; (b) when the Q-matrix has ties (e.g. clonal/duplicate-
  heavy data), neighbour selection is non-unique and can pick a different (still
  valid) NJ topology than RapidNJ.
- 2026-07-22: **Consolidated regression suite** (`regression/run_all.sh`): every
  method × matrix × handler × heuristic × recraft over synthetic + example data.
  Stable result: **92 byte-identical, 23 RF=0-equivalent, 1 NJ-tiebreak diff, 4
  skipped (Java Ninja unavailable), n_proc thread-invariance verified.** The
  synthetic generator was made deterministic (crc32 seed; Python `hash()` is
  salted per process and had made datasets vary run-to-run).
- `-c/--check` (estimate_Consumption) ported — time/memory match the oracle.
- 2026-07-22: **Efficient Edmonds** — replaced the recursive/HashMap Chu-Liu/
  Edmonds with a lazy skew-heap + rollback union-find (O(E log V), KACTL-style),
  reconstruction via cycle expansion. MSTreeV2 at N=1000 went 47.2s → 1.43s
  (33×); still fuzz-verified 40/40 vs the reference binary and byte-identical.
  Final speedups vs Python: distance 2.9–13×, MSTreeV2 1.3–6.9× (see BENCHMARKS.md).
- 2026-07-22: **Input formats** — gzip (`.gz`) verified (== plain output); inline
  text (non-file `-p` argument) verified byte-identical values. FASTA input:
  the reference **crashes** (`backend()` does `del ... part ...`, but `part` is
  never bound on the FASTA branch — an upstream bug), so it cannot be compared;
  grapetree-rs parses aligned FASTA correctly (each column = a locus) and runs.
  Two upstream bugs are thus fixed-by-porting: `--wgMLST` (dead code) and FASTA.
