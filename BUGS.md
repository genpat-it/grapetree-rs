# Upstream GrapeTree findings

While building this port we characterised three discrepancies in the reference
implementation (`achtman-lab/GrapeTree`, `module/MSTrees.py`). Reported here
constructively — GrapeTree is an excellent, widely‑used tool; these are edge
cases we had to understand to reason about fidelity. Line references are to the
`master` `module/MSTrees.py` at the time of writing; upstream may differ.

| # | Kind | Symptom | Scientific impact |
|---|------|---------|-------------------|
| 1 | Functional bug | `--wgMLST` silently does nothing | flag is a no‑op |
| 2 | Functional bug | FASTA input crashes the CLI | feature unusable standalone |
| 3 | Latent / non‑determinism‑by‑artifact | MSTreeV2 tree depends on implementation accidents | none (valid tree), but not reproducible by an independent reimplementation |

---

## Bug 1 — `--wgMLST` is dead code (silent no‑op)

**Where.** `backend()`:

```python
if params['wgMLST'] and params['matrix_type'] == 'asymmetric' :
    matrix_type = 'asymmetric_wgMLST'
```

**What's wrong.** This assigns a **local** variable `matrix_type`. Everything
downstream reads `params['matrix_type']` (or receives it via `**params`), which
is still `'asymmetric'`. So the `asymmetric_wgMLST` matrix — a whole implemented
code path (`distance_matrix.asymmetric_wgMLST`) — is **never reached**, and
`--wgMLST` has no effect on the output.

**Reproduce.**
```bash
grapetree.py -p profile.tsv -m distance -x asymmetric            > a.txt
grapetree.py -p profile.tsv -m distance -x asymmetric --wgMLST   > b.txt
diff a.txt b.txt      # identical
```
We verified: with missing data present, `--wgMLST` output is byte‑identical to
plain asymmetric.

**Fix.** `params['matrix_type'] = 'asymmetric_wgMLST'` (assign into `params`).

**In grapetree‑rs.** We reproduce the upstream behaviour for parity (`--wgMLST`
is a no‑op) **and** implement the intended `asymmetric_wgMLST` matrix, reachable
internally, so the corrected behaviour is available. See
`src/distance.rs::MatrixKind::resolve`.

---

## Bug 2 — FASTA input crashes the CLI

**Where.** `backend()`, right after the input‑parsing block:

```python
del fin, line, line_id, part, header
```

**What's wrong.** `part` (and `header`) are only bound in the **profile‑TSV**
branch. On a **FASTA** input those names are never assigned, so this `del`
raises `UnboundLocalError: cannot access local variable 'part'`. FASTA input is
advertised in `--help` ("*OR a fasta file containing aligned sequences*") but
cannot be processed by the standalone CLI.

**Reproduce.**
```bash
grapetree.py -p aligned.fasta -m distance
# UnboundLocalError: cannot access local variable 'part' where it is not associated with a value
```

**Fix.** Guard the `del` (e.g. `del fin, line, line_id`), or only delete names
that are bound on both code paths.

**In grapetree‑rs.** FASTA is parsed correctly (each aligned column = a locus),
so `-m MSTreeV2/distance/...` work on FASTA where the reference CLI aborts. See
`src/parse.rs`.

---

## Bug 3 — MSTreeV2 depends on implementation artifacts (latent)

**Not a correctness bug.** The tree produced is always a valid minimum spanning
tree with branch recrafting applied. The issue is that *which* of several
equally‑valid trees you get depends on incidental implementation details, not on
the algorithm — so the result is deterministic **for a given build** but is not
reproducible by an independent reimplementation.

**The mechanism.** In `methods._asymmetric`, the minimum spanning arborescence is
delegated to the compiled `edmonds` C binary via a text round‑trip:

```python
# reduced cost matrix written at 5 decimals, offset by ~1:
fout.write('{0:.5f}'.format(round(dmod[i][j]) + weight[i] + 0.999995))   # e.g. "4.00000"
mstree = Popen([edmonds_binary, dist_file]).communicate()[0]
mstree = np.array([...], dtype=float).astype(int)   # (a) TRUNCATES to int
mstree.T[2] -= 1                                     # then subtracts 1
```

Two accidents combine to set each branch's `brlen`:

- **(a)** `.astype(int)` **truncates** the branch length before `-= 1`
  (rounding would give a different integer at the `x.99999x` boundary).
- **(b)** the `edmonds` binary reports, for edges absorbed into a **contracted
  cycle** during Chu‑Liu/Edmonds, the **contraction‑reduced** weight (original
  cost minus the cycle's entering‑edge discount), not the original cost.

So an edge with reduced‑matrix cost 3 (`round(dmod)+weight`) can come back with
`brlen = 2` instead of `3`.

**Why it matters.** `methods._branch_recraft` is **order‑sensitive**: it sorts
branches by `brlen` and its commit condition is `branches[i+1][2] >= brlen`.
With integer allelic distances, equal‑`brlen` ties are everywhere. A `brlen`
that is off by the contraction discount flips a commit‑vs‑resort decision at a
tie, and because each commit mutates the group/adjacency state, the difference
**cascades** into a handful of different re‑attachments.

Crucially, the `brlen` from this step is **thrown away** later — final branch
lengths are recomputed as raw allelic (hamming) differences by
`distance_matrix.symmetric_link`. So `brlen` is only ever a *tie‑break key*, yet
its value is defined by (a) an int‑truncation and (b) a C binary's internal
weight‑reduction reporting. That is the smell: a tie‑break driven by a quantity
that was never meant to be semantically meaningful.

**Observed size.** On a clonal 10k dataset, reference vs an independent
reimplementation differ on ≈0.02–0.5% of splits (Robinson‑Foulds), with nearly
identical total tree length. Scientifically immaterial; both are valid MSTs.

**In grapetree‑rs.** We match the deterministic core exactly (distance matrix,
heuristic weights, and the Edmonds **edge set** are bit‑identical). MSTreeV2 is
byte‑identical on most inputs; on the rest it is topologically equivalent
(RF ≈ 0). To reach *bit*‑identity we would have to reproduce the C binary's
contraction‑reduced `brlen` for every edge (our skew‑heap Edmonds already tracks
those reductions internally via lazy deltas — a targeted future change). See
`src/edmonds.rs`, `src/recraft.rs`, and `DECISIONS.md`.

---

*This document is part of an independent, experimental reimplementation and is
not affiliated with or endorsed by the GrapeTree authors.*
