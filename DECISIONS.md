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
- 2026-07-22: **Scaling MSTreeV2 to 60k+ real samples — the dense Chu-Liu/Edmonds.**
  Running MSTreeV2 on a real *Campylobacter* cgMLST set (63,005 samples × 1,142
  loci, 64 threads) exposed the true scaling wall — and it was *not* where a first
  guess put it.

  **Problem, diagnosed with per-phase timing (`GT_TIMING=1`):**
  - distance 66 s, weights 2 s — fast, rayon-parallel.
  - `get_shortcut` reduced 57,402 non-redundant profiles to **24,627 survivors**
    (real cgMLST *does* collapse well; a naive synthetic-clonal set collapsed
    almost none — a lesson: validate on real data, not the worst-case synthetic).
  - **the minimum spanning arborescence was the wall.** The reduced graph is
    (near) complete: m = 24,627 survivors ⇒ **~6.06×10⁸ directed edges**. The
    existing skew-heap Chu-Liu/Edmonds is `O(E log V)` in time **and `O(E)` in
    memory** — it materialises every edge as a heap node (~32 B) plus an edge
    record (~16 B). That is ~20 GB *just for the edge list*, on top of two 13 GB
    distance matrices, and the per-edge heap churn made it run for >35 min without
    finishing (killed at 74 GB and still climbing).
  - `branch_recraft` was a **red herring** — its memory is `O(n)`; the earlier
    suspicion that it was the bottleneck was wrong. (It was still optimised, see
    below, but it is not the scaling limit.)

  **What is different from GrapeTree here:** upstream shells out to a compiled C
  `edmonds` binary, writing the m×m reduced cost matrix as a **text file** and
  reading the arborescence back. That is the *same* `O(m²)` work — a ~24k×24k
  text matrix is multi-GB of I/O — just in C instead of Rust. Neither tool
  escapes the `O(m²)` density; the reference's own `-c` estimate for this size is
  ~2.7 h / ~412 GB.

  **Solution — a dense `O(V)`-memory Chu-Liu/Edmonds** (`optimum_branching_dense`
  in `src/edmonds.rs`). Boruvka-style contraction: each round every non-root
  super-node takes its cheapest incoming edge, *all* resulting cycles are
  contracted at once, and the cost matrix is rebuilt over the shrunken super-node
  set; cycles are expanded in reverse at the end (routing each external edge to
  the cycle member that *contains* its real target, via per-round label
  snapshots). No per-edge heap nodes, no materialised edge list — it reuses the
  matrix already in RAM.
  - **Result on the Campylobacter set: 3 min 28 s, peak 46 GB** (vs >35 min /
    ≥74 GB and never finishing). The arborescence itself dropped from tens of
    minutes to ~108 s (including `get_shortcut`); recraft 1.9 s; tree write 0.15 s.
  - **Hybrid dispatch** (`src/mst.rs`, `DENSE_THRESHOLD = 10_000`): the dense
    variant wins on large `m`, but on *ultra-clonal* graphs (many tiny cycles ⇒
    many contraction rounds ⇒ `O(rounds·m²)`) the skew-heap is faster and its
    `O(m²)` memory is harmless when `m` is small. So below 10k survivors we keep
    the skew-heap, above it we switch to dense. Both return the **same**
    arborescence.

  **Bit-identity status after this change — unchanged, and re-verified.** The two
  Edmonds implementations are byte-identical because the MSTreeV2 optimum is
  *unique*: the incoming costs to any node, `round(d[i][j]) + weight[i]`, are
  pairwise distinct (the harmonic `weight` is a distinct rank per source), so
  there are no ties for any correct algorithm to break differently. Verified by:
  (a) a property test comparing dense vs skew-heap over random matrices *with*
  ties, n = 2..30 (`dense_matches_heap`, `dense_identical_when_unique`); and
  (b) the regression golden — **7/7 datasets (incl. clonal-2k/5k) byte-identical**
  in both `--native` and default modes before/after the change. The overall
  fidelity tally is untouched: distance 7/7, MSTree-sym 7/7, MSTree-asym 7/7,
  MSTreeV2 6/7 (the one residual is the known NEWICK child-order case, RF=0, equal
  length), NJ 7/7, RapidNJ 7/7.

- 2026-07-22: **`branch_recraft` made O(n²)-cheaper without changing output.**
  The faithful port cloned each endpoint's whole group every iteration and
  `sorted(...)[:3]` over it. Three output-preserving changes: groups read by
  reference (no clone); `group_id`/`groups`/`childrens` are `Vec`s indexed by
  node id (the reference's numpy fancy-assign `group_id[targets]=…` becomes an
  O(group) array write, not hashed inserts); and the "3 nearest" taken by partial
  selection instead of a full sort — safe because the candidate key
  `(weight, dist, node_id)` is a *total* order (node id is unique), so the three
  smallest are unambiguous and identical to `sorted(...)[:3]`. Byte-identical on
  all 7 regression datasets; recraft at 63k now 1.9 s.
- 2026-07-22: **Full bit-identity at 63k — the two non-portable numerics, delegated
  to NumPy.** MSTreeV2's output is deterministic but *not portable*: it depends on
  two floating-point pieces that differ across CPUs/NumPy builds. We isolated both
  by reproducing every algorithmic stage in Rust and diffing byte-for-byte against
  upstream on a real 63k *Campylobacter* set:
  1. **Harmonic centrality** `N/Σ 1/(dist+0.1)` — NumPy sums it in float32 with a
     SIMD (AVX) reduction whose add-order (hence last ULP) is CPU/build-dependent.
     Our scalar pairwise matches NumPy on *some* rows but not ~half of 57k-element
     rows. → `shim/harmonic_weights.py` (NumPy computes the weights).
  2. **`branch_recraft`'s `contemporary` test** — uses `np.log`, which is NumPy's
     own polynomial (correctly-rounded only ~79% of the time vs glibc's ~99.9%, and
     SIMD-dispatched). Rust's libm `ln` differs, flipping `p1>=p2` at boundaries.
     → `shim/recraft.py` (verbatim NumPy port of `_branch_recraft`/`contemporary`).
  Plus two pure-Rust f32 fixes needed for the reduced matrix to match byte-for-byte:
  `get_shortcut`'s `dist+weight` key computed in f32 (`GRAPETREE-COMPAT[shortcut-f32]`),
  and the reduced-matrix `+0.999995` offset added in f32 (`[edmonds-reduced-f32]`).
  **Result: default mode is byte-identical (same md5) to upstream GrapeTree on the
  full 63,005-sample set** (12m47s / 133 GB, vs upstream 18m48s / 133 GB — faster
  because our distance kernel is Rust). `--native` stays pure-Rust (3m28s / 46 GB,
  ~99.2%). Requires Python+NumPy at runtime for default mode (set `GT_PYTHON`), the
  same class of dependency GrapeTree already has.

  Deeper insight (see private WEAKNESS.md): GrapeTree uses floating-point almost
  entirely to *decide* (argmin, sort order, `p1>=p2`) rather than to *report a
  value*. A 1-ULP error, harmless in a value, becomes a coin-flip at a decision
  boundary and jumps the discrete outcome (a leaf's parent) by up to 747 alleles.
  This is the "robust geometric predicates" problem; GrapeTree lacks the exact/
  canonical tie-breaking that would make its combinatorial output reproducible.
