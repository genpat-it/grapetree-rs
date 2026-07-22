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
- 🦀 **Self-contained** — native Chu-Liu/Edmonds & Neighbor-Joining; no `edmonds`/FastME/RapidNJ/Ninja binaries
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

## Fidelity vs the reference

Verified against the Python reference (run in an isolated env) over synthetic
datasets (clonal structure, duplicates, missing-data gradients) and the bundled
examples. See `regression/run_all.sh`.

- **`distance`** — byte-identical across all matrix types × missing handlers.
- **`MSTree`** — byte-identical (symmetric); topologically identical (asymmetric).
- **`MSTreeV2`** — **bit-identical by default**: the minimum spanning
  arborescence is delegated to the bundled reference `edmonds` binary (as
  upstream does), then recrafted/serialised in Rust. A rare synthetic case can
  still differ only in NEWICK child ordering (RF=0, identical length) from a
  residual NumPy float32-dtype detail in the weight/reduced-matrix construction.
  `--native` uses the pure-Rust Edmonds instead (faster, self-contained,
  topologically equivalent but not bit-identical).
- **`NJ` / `RapidNJ`** — **bit-identical by default** by delegating to the bundled
  FastME/RapidNJ binaries + a small ete3 post-processing shim (`shim/`), exactly
  the reference toolchain. Requires the binaries + Python 3 with `ete3` at run
  time (set `GT_PYTHON`), just like GrapeTree. `--native` uses a pure-Rust
  canonical NJ instead (self-contained, topology identical, lengths differ).
- **`ninja`** — same delegation path (Java + `Ninja.jar`); note upstream `ninja`
  is broken on Java ≥ 9 (`java -d64` was removed), so it can't be validated
  against the reference on modern JVMs.

Aggregate: **92 byte-identical, 23 topology-identical (RF=0), 1 NJ tie-break
difference**, thread-count invariant.

Two upstream bugs are fixed simply by porting cleanly:
- `--wgMLST` is dead code in the reference (a local variable never propagates),
  so it silently does nothing; we preserve that behaviour for parity but also
  implement the intended `asymmetric_wgMLST` matrix.
- FASTA input crashes the reference CLI (`del part` on an unbound variable);
  grapetree-rs parses aligned FASTA correctly.

## Performance

Native Chu-Liu/Edmonds (lazy skew heap + rollback union-find, O(E log V)) and a
rayon-parallel distance kernel. Faster than the Python + compiled-binary pipeline
at every size tested (see `BENCHMARKS.md`):

| N | `distance` | `MSTreeV2` |
|---|-----------|-----------|
| 300 | 4.9× | 6.9× |
| 1000 | 2.6× | 2.2× (47s → 1.4s vs the old naive Edmonds) |
| 2000 | 2.9× | 1.3× |

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
src/edmonds.rs    optimum branching (Chu-Liu/Edmonds)
src/recraft.rs    branch_recraft + contemporary
src/nj.rs         canonical Neighbor-Joining
src/tree.rs       tree assembly + NEWICK writer
src/cli.rs        argument parsing
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

## Licence

GPL-3.0-or-later, matching upstream GrapeTree.
