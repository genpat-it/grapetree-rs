<h1 align="center">🍇🌳 grapetree-rs</h1>

<p align="center">
  <em>A fast, self-contained Rust port of GrapeTree's tree &amp; distance engine.</em>
</p>

<p align="center">
  <a href="https://github.com/genpat-it/grapetree-rs/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/genpat-it/grapetree-rs/actions/workflows/ci.yml/badge.svg"></a>
  <a href="LICENSE"><img alt="License: GPL v3" src="https://img.shields.io/badge/license-GPLv3-blue.svg"></a>
  <img alt="Rust 2021" src="https://img.shields.io/badge/rust-2021-orange.svg?logo=rust">
  <img alt="Status: experimental" src="https://img.shields.io/badge/status-experimental-yellow.svg">
  <img alt="Fidelity" src="https://img.shields.io/badge/vs%20reference-byte--identical%20(most%20params)-brightgreen.svg">
</p>

> ⚠️ **Experimental.** This is an independent, work-in-progress reimplementation
> for research and evaluation. Output matches the reference on the vast majority
> of inputs (see *Fidelity* below), but it is not yet a validated drop-in
> replacement — verify results against upstream GrapeTree for production use.

Turns an allelic profile (cgMLST / wgMLST / SNP) into a **NEWICK tree** or a
**PHYLIP distance matrix** — a Rust port of
[GrapeTree](https://github.com/achtman-lab/GrapeTree)'s computational backend.

### Highlights

- 🧬 **All methods** — `MSTreeV2`, `MSTree`, `NJ`, `RapidNJ`, `ninja`, `distance`, full CLI parity
- 🎯 **Faithful** — byte-identical to the reference on most parameter combinations
- 🦀 **No `edmonds` binary** — MSTreeV2's arborescence is a byte-identical pure-Rust port of it (default); `--native` also replaces the NJ binaries
- ⚡ **Fast** — rayon-parallel distance kernel; `distance` 2.9–13×, `MSTreeV2` 1.3–6.9× vs Python
- ✅ **Tested** — per-parameter regression harness against the Python reference

> Scope: this ports the scientific engine (`module/MSTrees.py`). The Flask web
> server and the D3 visualiser that ship with GrapeTree are UI shells around it
> and are out of scope.

## Build

```bash
cd grapetree-rs
cargo build --release      # binary at target/release/grapetree
cargo test --release       # unit + property tests
```

## Docker

A prebuilt image is published to GHCR on every release. Pull `:latest` or a
version tag:

```bash
# if the GHCR package is private, authenticate first (PAT with read:packages):
#   echo "$PAT" | docker login ghcr.io -u <username> --password-stdin
docker pull ghcr.io/genpat-it/grapetree-rs:latest      # or :0.1.1
```

The image runs the **default byte-identical mode**: the minimum spanning
arborescence is the **pure-Rust `edmonds` port** (no C binary), and the two
non-portable numerics use the bundled NumPy shims — so the output is
**byte-identical** to upstream GrapeTree.

Mount the directory with your profile and pass grapetree's usual flags after the
image name:

```bash
# MSTreeV2 (default) — byte-identical to upstream, edmonds in Rust, no C binary:
docker run --rm -v "$PWD":/data ghcr.io/genpat-it/grapetree-rs:latest \
    -p /data/profiles.tsv -m MSTreeV2 > tree.nwk

# PHYLIP distance matrix:
docker run --rm -v "$PWD":/data ghcr.io/genpat-it/grapetree-rs:latest \
    -p /data/profiles.tsv -m distance > dist.tsv

# symmetric MST:
docker run --rm -v "$PWD":/data ghcr.io/genpat-it/grapetree-rs:latest \
    -p /data/profiles.tsv -m MSTree -x symmetric > sym.nwk
```

`distance` / `MSTree` / `MSTreeV2` work out of the box; the NJ family needs the
optional `ete3` + FastME/RapidNJ/Ninja layer (see the `Dockerfile`). Pass
`--native` for a fully self-contained (topology-only, not byte-identical) run.

## Usage

```bash
grapetree -p <profile> [-m METHOD] [options] > tree.nwk
```

Input may be a MLST/SNP profile TSV, an aligned FASTA (each column = a locus),
a `.gz` of either, or the file contents passed inline.

| Flag | Meaning | Default |
|------|---------|---------|
| `-p, --profile` | input file / `.gz` / inline text | (required) |
| `-m, --method` | `MSTreeV2`, `MSTree`, `NJ`, `RapidNJ`, `ninja`, `distance` | `MSTreeV2` |
| `-x, --matrix` | `symmetric`, `asymmetric`, `blockwise` | `symmetric` (MSTreeV2 forces `asymmetric`) |
| `-y, --missing` | `0` pair-delete, `1` complete-delete, `2` as-allele, `3` absolute | `0` |
| `-t, --heuristic` | `eBurst`, `harmonic` (tie-break for MST) | `eBurst` (MSTreeV2 forces `harmonic`) |
| `-r, --recraft` | local branch recrafting | off (MSTreeV2 forces on) |
| `-w, --wgMLST` | wgMLST support (see note) | off |
| `-n, --n_proc` | worker threads (performance only) | 5 |
| `-b, --block_penalty` | penalty for `blockwise` | 0.01 |
| `-c, --check` | print estimated time/memory JSON | off |

`MSTreeV2` is an alias for `MSTree -x asymmetric -t harmonic -r`.

### All options (`grapetree --help`)

```text
Rust port of GrapeTree: NEWICK tree / distance matrix from allelic profiles

Usage: grapetree [OPTIONS] --profile <PROFILE>

Options:
  -p, --profile <PROFILE>              Input file (MLST/SNP profile TSV, or aligned FASTA). `.gz` supported
  -m, --method <METHOD>                MSTreeV2 [default], MSTree, NJ, RapidNJ, ninja, distance [default: MSTreeV2]
  -x, --matrix <MATRIX_TYPE>           symmetric [default for MSTree/NJ], asymmetric, blockwise [default: symmetric]
  -r, --recraft                        Trigger local branch recrafting (forced on for MSTreeV2)
  -y, --missing <HANDLER>              Missing-data handler: 0 pair_delete [default], 1 complete_delete, 2 as_allele, 3 absolute_distance. (symmetric distance matrix only) [default: 0]
  -w, --wgMLST                         [experimental] better support for wgMLST schemes
  -t, --heuristic <HEURISTIC>          Tiebreak heuristic: eBurst [default MSTree] or harmonic [default MSTreeV2] [default: eBurst]
  -n, --n_proc <N_PROC>                Number of parallel worker threads (performance only; results unchanged) [default: 5]
  -c, --check                          Only report estimated time/memory requirements
  -b, --block_penalty <BLOCK_PENALTY>  Penalty for a different locus led by another difference (blockwise only) [default: 0.01]
      --native                         Use the pure-Rust native arborescence/NJ instead of the bundled reference binaries. Faster and self-contained, but only topologically equivalent (RF=0), not bit-identical to upstream GrapeTree. Default: bit-identical
  -h, --help                           Print help
```

Diagnostics: set `GT_TIMING=1` to print per-phase timings (distance / weights /
arborescence / recraft / tree) to stderr — useful on large inputs.

## Fidelity vs the reference

Verified against the Python reference (run in an isolated env) over synthetic
datasets (clonal structure, duplicates, missing-data gradients) and the bundled
examples. See `regression/run_all.sh`.

- **`distance`** — byte-identical across all matrix types × missing handlers.
- **`MSTree`** — byte-identical (symmetric); topologically identical (asymmetric).
- **`MSTreeV2`** — **byte-identical by default**, verified on a real 63,005-sample
  *Campylobacter* cgMLST set (same md5 as upstream). The minimum spanning
  arborescence is now a **faithful pure-Rust port of the `edmonds` C binary**
  (`src/edmonds_tofigh.rs`, Ali Tofigh's `edmonds_optimum_branching` — same
  critical-edge tie-break, cycle contraction and edge *emission order*), so **no
  external binary is needed** and the arborescence is byte-identical (verified: the
  full regression suite 116/116 byte-identical, and campy 63k same md5). Only the
  two numerical steps that are *not* portably reproducible in Rust stay in NumPy —
  the harmonic weights (`shim/harmonic_weights.py`; NumPy's float32 SIMD `sum`) and
  `branch_recraft` (`shim/recraft.py`; NumPy's `np.log` in the `contemporary`
  test). Needs Python 3 + NumPy at run time (set `GT_PYTHON`). `GT_EDMONDS_BINARY=1`
  forces the bundled C binary instead of the port (for cross-validation).
  `--native` forces a fully pure-Rust path (native weights + recraft too):
  self-contained and fast, topologically ~99% but not byte-identical (see below).
- **`NJ` / `RapidNJ`** — **bit-identical by default** by delegating to the bundled
  FastME/RapidNJ binaries + a small ete3 post-processing shim (`shim/`), exactly
  the reference toolchain. Requires the binaries + Python 3 with `ete3` at run
  time (set `GT_PYTHON`), just like GrapeTree. `--native` uses a pure-Rust
  canonical NJ instead (self-contained, topology identical, lengths differ).
- **`ninja`** — same delegation path (Java + `Ninja.jar`); note upstream `ninja`
  is broken on Java ≥ 9 (`java -d64` was removed), so it can't be validated
  against the reference on modern JVMs.

Aggregate (default mode, Rust `edmonds` port): **116 byte-identical, 0 failed**,
thread-count invariant, gzip==plain, FASTA runs (`regression/run_all.sh`).

Two upstream bugs are fixed simply by porting cleanly:
- `--wgMLST` is dead code in the reference (a local variable never propagates),
  so it silently does nothing; we preserve that behaviour for parity but also
  implement the intended `asymmetric_wgMLST` matrix.
- FASTA input crashes the reference CLI (`del part` on an unbound variable);
  grapetree-rs parses aligned FASTA correctly.

## Performance

Native Chu-Liu/Edmonds and a rayon-parallel distance kernel. Faster than the
Python + compiled-binary pipeline at every size tested (see `BENCHMARKS.md`):

| N | `distance` | `MSTreeV2` |
|---|-----------|-----------|
| 300 | 4.9× | 6.9× |
| 1000 | 2.6× | 2.2× (47s → 1.4s vs the old naive Edmonds) |
| 2000 | 2.9× | 1.3× |

**At 60k+ samples** the minimum spanning arborescence (not the distance step)
dominates, because the reduced graph is near-complete. The default bit-identical
path runs the arborescence through the **pure-Rust port of the `edmonds` binary**,
which is both far faster and far lighter than shelling out to the C binary. On a
real *Campylobacter* set (63,005 samples × 1,142 loci, 64 threads):

| MSTreeV2, byte-identical (same md5) | wall | peak RAM |
|---|---|---|
| upstream GrapeTree (Python + `edmonds` C binary) | 18 min 48 s | 139 GB |
| grapetree-rs via the bundled `edmonds` C binary | 11 min 35 s | 139 GB |
| **grapetree-rs default (Rust `edmonds` port)** | **3 min 48 s** | **39 GB** |

The Rust port does the arborescence in ~19 s vs the C binary's ~480 s (no
5 GB text round-trip, no 606M-edge Boost graph). `--native` (topological, not
byte-identical) is comparable in speed and RAM. See `BENCHMARKS.md`.

## Regression suite

```bash
# needs a Python 3.11 env with numpy/networkx/numba/ete3 (ete3 needs <3.13)
GT_ORIG=/path/to/GrapeTree bash regression/run_all.sh
```

Generates deterministic synthetic inputs (`gen_synth.py`), runs every method ×
parameter through both implementations, and scores byte-identity or topological
equivalence (`compare_trees.py`).

## Layout

```
src/parse.rs      profile/FASTA reader + nonredundant encoding
src/distance.rs   symmetric/asymmetric/wgMLST/blockwise matrices (rayon)
src/heuristic.rs  harmonic / eBurst weights
src/mst.rs        Kruskal MST + get_shortcut + symmetric_link
src/edmonds.rs    optimum branching (Chu-Liu/Edmonds): skew-heap + dense O(V)-mem
src/edmonds_tofigh.rs  pure-Rust port of the `edmonds` binary (byte-identical, default)
src/recraft.rs    branch_recraft + contemporary (native, pure Rust)
src/nj.rs         canonical Neighbor-Joining
src/tree.rs       tree assembly + NEWICK writer
src/cli.rs        argument parsing
shim/harmonic_weights.py  NumPy harmonic weights   (default-mode bit-identity)
shim/recraft.py           NumPy branch_recraft     (default-mode bit-identity)
shim/nj_postprocess.py    ete3 NJ post-processing  (NJ-family bit-identity)
```

See `DECISIONS.md` for the full design rationale and a dated build log.

## Credits

All credit for GrapeTree — the algorithms (MSTreeV2, branch recrafting, the
contemporaneity test), the science, and the original implementation — belongs to
its authors. This project is only a Rust reimplementation of their work.

- **GrapeTree** by **Zhemin Zhou**, Nabil-Fareed Alikhan, Martin J. Sergeant,
  Nina Luhmann, Cátia Vaz, Alexandre P. Francisco, João André Carriço and
  **Mark Achtman** — EnteroBase / achtman-lab.
- Paper: Zhou Z. *et al.* "GrapeTree: visualization of core genomic relationships
  among 100,000 bacterial pathogens." *Genome Research* 28:1395–1404 (2018).
  doi:10.1101/gr.232397.117
- Upstream: <https://github.com/achtman-lab/GrapeTree>

This port (`grapetree-rs`) is developed and maintained by the **GenPat Team**
(IZSAM). It is not affiliated with or endorsed by the original GrapeTree authors.

## AI disclosure

This Rust port was carried out with the assistance of **Claude** (Anthropic,
Claude Code). The AI helped with translating the reference algorithms, the
optimisation work, and the byte-for-byte regression testing against upstream
GrapeTree. All output was reviewed and validated by the maintainers, and the
scientific credit for the original methods remains entirely with the GrapeTree
authors (see *Credits*).

## Licence

GPL-3.0-or-later, matching upstream GrapeTree.
