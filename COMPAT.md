# Compatibility anchors — deliberate reference-matching in grapetree-rs

To reproduce GrapeTree bit-for-bit we sometimes copy behaviour that is
counter-intuitive, or is an outright upstream bug (see `BUGS.md`). Every such
place is marked in the source with an anchor comment so it is never mistaken for
our own logic:

```
// GRAPETREE-COMPAT[<id>]: <one-line why>  (see COMPAT.md / BUGS.md#<id>)
```

Grep the tree with `rg "GRAPETREE-COMPAT"` to list them all. Registry:

| id | kind | where | what we replicate |
|----|------|-------|-------------------|
| `wgmlst-noop` | bug (BUGS.md#1) | `src/distance.rs` `MatrixKind::resolve` | `--wgMLST` does nothing (dead-code in upstream); we keep the no-op for parity and also implement the intended matrix |
| `fasta-ok` | bug (BUGS.md#2) | `src/parse.rs` | upstream FASTA CLI crashes; we parse it correctly (divergence *in our favour*, documented) |
| `complete-delete-cols` | quirk | `src/parse.rs` `nonredundant` | `complete_delete` keeps only columns that *contain* a missing marker — verbatim from the reference, counter-intuitive but required |
| `missing-encode-shift` | quirk | `src/parse.rs` `nonredundant` | missing markers (`0`/`N`/`-`) participate in the per-column sorted-unique set, shifting allele codes — matches `np.unique(...)+1` exactly |
| `nx-edge-order` | quirk | `src/mst.rs` `symmetric_mst` | MST edges emitted in `nx.Graph.edges()` adjacency order (not Kruskal order) so NEWICK child ordering matches ete3 |
| `edmonds-brlen-roundtrip` | bug (BUGS.md#3) | `src/mst.rs` `asymmetric_network` | recraft branch length = the `edmonds` binary's value through a `%.5f` text round-trip + `astype(int)` truncation + `-1`; only fully bit-exact via the binary (see BUGS.md#3) |
| `distance-trailing-nl` | quirk | `src/main.rs` | PHYLIP output newline handling matches the reference CLI |

## Status toward full bit-identity

- `distance`, `MSTree` (symmetric): **bit-identical**.
- `MSTreeV2`: bit-identical on most inputs; the residual is `edmonds-brlen-roundtrip`
  (BUGS.md#3). Guaranteed bit-identity requires the reference `edmonds` binary.
- `NJ`/`RapidNJ`/`ninja`: topologically identical (RF=0); bit-identity requires
  the FastME/RapidNJ/Ninja binaries (different branch-length scheme).

A `--exact` mode that shells out to the bundled reference binaries would make all
methods bit-identical; the native default trades that for pure-Rust speed and
self-containment with scientific equivalence. See DECISIONS.md.
