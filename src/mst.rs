//! Minimum spanning tree construction.
//!
//! `symmetric` MST (Kruskal + union-find) reproducing NetworkX's edge ordering
//! so tie-breaking matches the reference. Edge weights for the final tree are
//! recomputed as raw allelic (hamming) differences by [`symmetric_link`],
//! mirroring `distance_matrix.symmetric_link`.

use crate::distance::DistMatrix;
use crate::edmonds;
use crate::parse::Parsed;
use crate::HandleMissing;

/// Union-find with union-by-size + path compression.
struct DisjointSet {
    parent: Vec<usize>,
    size: Vec<usize>,
}
impl DisjointSet {
    fn new(n: usize) -> Self {
        DisjointSet {
            parent: (0..n).collect(),
            size: vec![1; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        // path compression
        let mut c = x;
        while self.parent[c] != r {
            let next = self.parent[c];
            self.parent[c] = r;
            c = next;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        let (big, small) = if self.size[ra] >= self.size[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[small] = big;
        self.size[big] += self.size[small];
        true
    }
}

/// `_symmetric`: MST edge topology over `w(i,j) = round(d[i,j]) + min(w_i, w_j)`.
///
/// Edges are enumerated in row-major upper-triangle order and stably sorted by
/// weight — identical to `nx.Graph(matrix)` insertion order + `nx.kruskal`'s
/// stable sort — so tie-breaks match the reference. Zero-weight cells are
/// treated as absent edges, exactly as `nx.Graph` does.
pub fn symmetric_mst(dist: &DistMatrix, weight: &[f64]) -> Vec<(usize, usize)> {
    let n = dist.n;
    // (weight_edge, insertion_index, i, j)
    let mut edges: Vec<(f64, usize, usize, usize)> = Vec::with_capacity(n * (n - 1) / 2);
    let mut ins = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let rd = (dist.get(i, j) as f64).round();
            let we = rd + weight[i].min(weight[j]);
            if we != 0.0 {
                edges.push((we, ins, i, j));
            }
            ins += 1;
        }
    }
    // stable sort by weight; ties preserve insertion order (== networkx behaviour)
    edges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut ds = DisjointSet::new(n);
    // adjacency in Kruskal-add order, to reproduce NetworkX graph insertion order
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (_, _, i, j) in edges {
        if ds.union(i, j) {
            adj[i].push(j);
            adj[j].push(i);
        }
    }
    // Emit edges in `nx.Graph.edges()` order: iterate nodes 0..n and, for each,
    // its neighbours in adjacency-insertion order, yielding each undirected edge
    // from whichever endpoint is visited first (NetworkX marks a node "seen"
    // only after iterating it). This is the order the reference passes to
    // `symmetric_link` / `_network2tree`, so NEWICK child ordering matches.
    let mut done = vec![false; n];
    let mut mst = Vec::with_capacity(n.saturating_sub(1));
    for u in 0..n {
        for &v in &adj[u] {
            if !done[v] {
                mst.push((u, v));
            }
        }
        done[u] = true;
    }
    mst
}

/// `get_shortcut`: for each target `t`, the lowest-cost more-central source `s`
/// (`weight[s] < weight[t]`, `dist[s][t] < cutoff+1`) minimising
/// `dist[s][t] + weight[s]`. Returned ordered by target ascending (np.unique).
fn get_shortcut(dist: &DistMatrix, weight: &[f64]) -> Vec<(usize, usize)> {
    let n = dist.n;
    let cutoff = if n < 3000 {
        2
    } else if n < 10000 {
        5
    } else if n < 30000 {
        10
    } else {
        20
    } as f64;

    let mut shortcuts = Vec::new();
    for t in 0..n {
        let mut best: Option<(f64, usize)> = None; // (dist+weight[s], s)
        for s in 0..n {
            if s == t {
                continue;
            }
            if weight[s] < weight[t] && (dist.get(s, t) as f64) < cutoff + 1.0 {
                let key = dist.get(s, t) as f64 + weight[s];
                match best {
                    Some((bk, _)) if bk <= key => {}
                    _ => best = Some((key, s)),
                }
            }
        }
        if let Some((_, s)) = best {
            shortcuts.push((s, t));
        }
    }
    shortcuts
}

/// `_asymmetric`: shortcut collapse + minimum spanning arborescence over the
/// reduced node set. Returns network edges `(from, to, brlen)` over the
/// deduplicated profile indices, where `brlen` reproduces the reference's
/// initial branch lengths (arborescence: `round(d_mod)+weight[src]-0.000005`
/// after the binary's `+0.999995`/`-1` round trip; shortcut: `int(d+weight)`).
/// These feed `branch_recraft`; the final weights are set by [`symmetric_link`].
pub fn asymmetric_network(dist: &DistMatrix, weight: &[f64]) -> Vec<(usize, usize, f64)> {
    let n = dist.n;
    let shortcuts = get_shortcut(dist, weight);

    // working copy: collapse each target into its source (row-wise min).
    let mut d: Vec<f32> = dist.data.clone();
    for &(s, t) in &shortcuts {
        for k in 0..n {
            let (a, b) = (d[s * n + k], d[t * n + k]);
            if a > b {
                d[s * n + k] = b;
            }
        }
    }

    // surviving nodes (shortcut targets removed), ascending original order.
    let mut removed = vec![false; n];
    for &(_, t) in &shortcuts {
        removed[t] = true;
    }
    let kept: Vec<usize> = (0..n).filter(|&i| !removed[i]).collect();
    let m = kept.len();
    if std::env::var("GT_TIMING").is_ok() {
        eprintln!(
            "[timing] n={n} shortcuts={} survivors(m)={m} edmonds_edges={}",
            shortcuts.len(),
            m * m
        );
    }

    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    if m >= 2 {
        let cost = |i: usize, j: usize| -> f64 {
            let (oi, oj) = (kept[i], kept[j]);
            (d[oi * n + oj] as f64).round() + weight[oi]
        };
        // Two byte-identical Chu-Liu/Edmonds implementations, chosen by size:
        //   - skew-heap (O(E log V) time, O(E)=O(m²) memory): fast for small m,
        //     and robust when the reduced graph is ultra-clonal (many tiny cycles);
        //   - dense matrix (O(V) memory beyond the matrix): the skew-heap's m²
        //     edge list blows up past ~10k survivors (~20 GB at m≈25k), so above
        //     the threshold we reuse the matrix already in RAM instead.
        // Both return the same arborescence (the optimum is unique for MSTreeV2).
        const DENSE_THRESHOLD: usize = 10_000;
        let arb = if m >= DENSE_THRESHOLD {
            edmonds::optimum_branching_dense(m, cost)
        } else {
            edmonds::optimum_branching(m, |i, j| Some(cost(i, j)))
        };
        for (u, v) in arb {
            let (os, ot) = (kept[u], kept[v]);
            // Reproduce the reference's edmonds text round-trip exactly: the
            // cost is written as `%.5f` of `round(d)+weight+0.999995`, the binary
            // echoes it, then `-= 1`. The 5-dp quantisation matters because it
            // can flip near-ties in the recraft brlen re-sort.
            let raw = (d[os * n + ot] as f64).round() + weight[os] + 0.999995;
            let brlen = (raw * 1e5).round() / 1e5 - 1.0;
            edges.push((os, ot, brlen));
        }
    }
    // shortcut edges: brlen = int(dist_orig[s][t] + weight[s])
    for &(s, t) in &shortcuts {
        let brlen = (dist.get(s, t) as f64 + weight[s]).trunc();
        edges.push((s, t, brlen));
    }
    edges
}

/// Round half to even (NumPy's `np.round` convention), so the reduced-matrix
/// text we hand to the `edmonds` binary matches the reference byte-for-byte.
fn np_round(x: f64) -> f64 {
    let f = x.floor();
    let diff = x - f;
    if diff < 0.5 {
        f
    } else if diff > 0.5 {
        f + 1.0
    } else if (f as i64) % 2 == 0 {
        f
    } else {
        f + 1.0
    }
}

/// Exact `_asymmetric`: identical to [`asymmetric_network`] but the minimum
/// spanning arborescence is delegated to the reference `edmonds` binary, so the
/// per-edge branch lengths (and hence the order-sensitive branch recrafting)
/// are **bit-identical** to upstream GrapeTree — see BUGS.md#3 / COMPAT.md
/// (`edmonds-brlen-roundtrip`). Falls back to `None` if the binary can't run.
pub fn asymmetric_network_exact(
    dist: &DistMatrix,
    weight: &[f64],
    edmonds_path: &std::path::Path,
) -> Option<Vec<(usize, usize, f64)>> {
    use std::process::Command;
    let n = dist.n;
    let shortcuts = get_shortcut(dist, weight);

    // row-min collapse of shortcut targets into their sources (same as native).
    let mut d: Vec<f32> = dist.data.clone();
    for &(s, t) in &shortcuts {
        for k in 0..n {
            if d[s * n + k] > d[t * n + k] {
                d[s * n + k] = d[t * n + k];
            }
        }
    }
    let mut removed = vec![false; n];
    for &(_, t) in &shortcuts {
        removed[t] = true;
    }
    let kept: Vec<usize> = (0..n).filter(|&i| !removed[i]).collect();
    let m = kept.len();

    let mut edges: Vec<(usize, usize, f64)> = Vec::new();
    if m >= 2 {
        // reduced cost matrix = round(dmod)+weight[src], diagonal 0, written as
        // "%.5f" of (value + 0.999995) — the exact reference file format.
        let mut buf = String::with_capacity(m * m * 8);
        for i in 0..m {
            for j in 0..m {
                if j > 0 {
                    buf.push('\t');
                }
                // GRAPETREE-COMPAT[edmonds-reduced-f32]: the reference builds the
                // reduced matrix as `np.round(dmod) + weight2` in float32 (both
                // operands float32), then writes `value + (1 - 0.000005)` at
                // %.5f. We must do the round+weight add in f32 too, else a ~1e-8
                // difference flips the 5th decimal and the binary's tie-breaks.
                let v = if i == j {
                    0.0f64
                } else {
                    ((np_round(d[kept[i] * n + kept[j]] as f64) as f32) + (weight[kept[i]] as f32))
                        as f64
                };
                buf.push_str(&format!("{:.5}", v + (1.0 - 0.000005)));
            }
            buf.push('\n');
        }
        let path = std::env::temp_dir().join(format!("gt_edmonds_{}.list", std::process::id()));
        std::fs::write(&path, &buf).ok()?;
        let out = Command::new(edmonds_path).arg(&path).output();
        let _ = std::fs::remove_file(&path);
        let out = out.ok()?;
        if !out.status.success() {
            return None;
        }
        let stdout = String::from_utf8_lossy(&out.stdout);
        let mut any = false;
        for line in stdout.lines() {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() < 3 {
                continue;
            }
            // reference parses the row as float then `.astype(int)` (truncates),
            // then subtracts 1 from the weight column.
            let s = p[0].parse::<f64>().ok()?.trunc() as usize;
            let t = p[1].parse::<f64>().ok()?.trunc() as usize;
            let w = p[2].parse::<f64>().ok()?.trunc() - 1.0;
            edges.push((kept[s], kept[t], w));
            any = true;
        }
        if !any {
            return None;
        }
    }
    // shortcut edges: brlen = int(dist_orig[s][t] + weight[s])
    for &(s, t) in &shortcuts {
        edges.push((s, t, (dist.get(s, t) as f64 + weight[s]).trunc()));
    }
    Some(edges)
}

/// Per-row presence used by `symmetric_link`, matching the reference.
fn presence_matrix(p: &Parsed, hm: HandleMissing) -> PresenceMask {
    match hm {
        HandleMissing::AsAllele => PresenceMask::All,
        HandleMissing::PairDelete | HandleMissing::AbsoluteDistance => PresenceMask::PerCell,
        HandleMissing::CompleteDelete => {
            let l = p.n_cols;
            let mask: Vec<bool> = (0..l)
                .map(|c| (0..p.n_rows).all(|r| p.codes[r * l + c] > 0))
                .collect();
            PresenceMask::PerColumn(mask)
        }
    }
}

enum PresenceMask {
    All,
    PerCell,
    PerColumn(Vec<bool>),
}
impl PresenceMask {
    #[inline]
    fn present(&self, code: u32, col: usize) -> bool {
        match self {
            PresenceMask::All => true,
            PresenceMask::PerCell => code > 0,
            PresenceMask::PerColumn(m) => m[col],
        }
    }
}

/// `symmetric_link`: recompute each edge's weight as the raw allelic difference
/// count (hamming over loci present in both, per the handler).
pub fn symmetric_link(
    p: &Parsed,
    edges: &[(usize, usize)],
    hm: HandleMissing,
) -> Vec<(usize, usize, u32)> {
    let l = p.n_cols;
    let mask = presence_matrix(p, hm);
    edges
        .iter()
        .map(|&(s, t)| {
            let ps = p.row(s);
            let pt = p.row(t);
            let mut d = 0u32;
            for k in 0..l {
                if mask.present(ps[k], k) && mask.present(pt[k], k) && ps[k] != pt[k] {
                    d += 1;
                }
            }
            (s, t, d)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::DistMatrix;

    #[test]
    fn mst_is_a_spanning_tree() {
        // path-graph-ish distances -> MST has n-1 edges connecting all nodes
        let n = 5;
        let mut data = vec![0f32; n * n];
        for i in 0..n {
            for j in 0..n {
                data[i * n + j] = ((i as i32 - j as i32).abs()) as f32;
            }
        }
        let m = DistMatrix { n, data };
        let w = vec![0.0, 0.2, 0.4, 0.6, 0.8];
        let mst = symmetric_mst(&m, &w);
        assert_eq!(mst.len(), n - 1, "spanning tree must have n-1 edges");
    }
}
