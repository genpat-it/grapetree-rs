//! Faithful Rust port of GrapeTree's bundled `edmonds` binary — Ali Tofigh's
//! `edmonds_optimum_branching<false, true, false>` (minimum, spanning, sparse),
//! built on Boost Graph. Ported so the bit-identical MSTreeV2 path can drop the
//! external C binary (~480 s / ~130 GB) while producing the **same** arborescence
//! edges **in the same emission order** the binary prints — which matters because
//! the downstream `branch_recraft` is order-sensitive.
//!
//! The port mirrors `edmonds_optimum_branching_impl.hpp` structure-for-structure:
//! per-target in-edge lists sorted by source ascending, a LIFO `roots` stack over
//! `0..n`, critical-edge = min weight with ties broken by lowest source, cycle
//! contraction breaking the maximum-weight cycle edge, and the F/λ/children
//! expansion that fixes the output order. Edge weights are held as `f32`-widened
//! `f64` because the binary reads them via `float atof(...)`.
//!
//! The graph is the reduced cost matrix as a *complete* digraph on `n` super-nodes
//! (every `i -> j`, `i != j`); no virtual root — spanning falls out of contraction.

/// Boost `disjoint_sets_with_storage`-style union-find: union by rank + full path
/// compression. `link(x, y)` assumes `x`, `y` are set representatives.
struct Dsu {
    parent: Vec<u32>,
    rank: Vec<u8>,
}
impl Dsu {
    fn new(n: usize) -> Self {
        Dsu {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
        }
    }
    fn find(&mut self, x: usize) -> usize {
        // iterative find + full path compression
        let mut r = x;
        while self.parent[r] as usize != r {
            r = self.parent[r] as usize;
        }
        let mut c = x;
        while self.parent[c] as usize != r {
            let nxt = self.parent[c] as usize;
            self.parent[c] = r as u32;
            c = nxt;
        }
        r
    }
    /// Link two representatives (Boost semantics: higher rank stays root; on tie,
    /// `y` stays root and its rank grows).
    fn link(&mut self, x: usize, y: usize) {
        if x == y {
            return;
        }
        if self.rank[x] > self.rank[y] {
            self.parent[y] = x as u32;
        } else {
            self.parent[x] = y as u32;
            if self.rank[x] == self.rank[y] {
                self.rank[y] += 1;
            }
        }
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        self.link(ra, rb);
    }
}

const NIL: i64 = -1;

/// Minimum spanning arborescence over a complete dense cost function, matching
/// the bundled `edmonds` binary edge-for-edge and in the same emission order.
/// `cost(i, j)` is the reduced-matrix value for `i -> j` (already the `%.5f`
/// quantised value); it is stored as `f32`-widened `f64`, as the binary does.
/// Returns arborescence edges `(source, target)` in the binary's print order.
pub fn optimum_branching_tofigh(
    n: usize,
    cost: impl Fn(usize, usize) -> f64,
) -> Vec<(usize, usize)> {
    use std::collections::{HashMap, HashSet};
    if n <= 1 {
        return Vec::new();
    }
    // Per-edge arrays for every i->j (i != j). Only source/target/weight need to
    // be dense; the F-forest links (parent/children) and the removed flag touch
    // only the O(V) critical edges, so they live in sparse maps — saving ~20 GB at
    // 63k (a dense children `Vec<Vec>` alone would be ~606M near-empty headers).
    // Edge ids are u32 (n_edges < 2^32); in_edges hold u32 too. in_edges[v] is
    // built source-ascending (== the radix-sort-by-source result).
    let n_edges = n * (n - 1);
    let mut e_src: Vec<u32> = Vec::with_capacity(n_edges);
    let mut e_tgt: Vec<u32> = Vec::with_capacity(n_edges);
    let mut e_w: Vec<f64> = Vec::with_capacity(n_edges);
    let mut in_edges: Vec<Vec<u32>> = vec![Vec::new(); n];
    for v in 0..n {
        in_edges[v].reserve(n - 1);
        for u in 0..n {
            if u == v {
                continue;
            }
            let id = e_src.len() as u32;
            e_src.push(u as u32);
            e_tgt.push(v as u32);
            e_w.push(cost(u, v) as f32 as f64);
            in_edges[v].push(id);
        }
    }
    // sparse F-forest bookkeeping, keyed by edge id (only critical edges appear).
    let mut parent: HashMap<u32, u32> = HashMap::new();
    let mut children: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut removed_f: HashSet<u32> = HashSet::new();

    let mut s = Dsu::new(n);
    let mut w_dsu = Dsu::new(n);
    let mut min_: Vec<usize> = (0..n).collect();
    let mut enter: Vec<i64> = vec![NIL; n];
    let mut lambda: Vec<i64> = vec![NIL; n];
    let mut cycle: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut edge_weight_change: Vec<f64> = vec![0.0; n];
    let mut f_list: Vec<u32> = Vec::new();
    let mut roots: Vec<usize> = (0..n).collect(); // stack: pop from back
    let mut final_roots: Vec<usize> = Vec::new();

    while let Some(cur_root) = roots.pop() {
        if in_edges[cur_root].is_empty() {
            final_roots.push(min_[cur_root]);
            continue;
        }
        // critical_edge = min-weight in-edge; ties keep the earliest (lowest source).
        let mut critical = in_edges[cur_root][0] as usize;
        for &en in &in_edges[cur_root] {
            let en = en as usize;
            if e_w[en] < e_w[critical] {
                critical = en;
            }
        }
        // TAttemptToSpan = true → always insert the critical edge.
        f_list.push(critical as u32);
        // children of critical_edge = cycle[cur_root]
        let cyc_children = std::mem::take(&mut cycle[cur_root]);
        for &en in &cyc_children {
            parent.insert(en, critical as u32);
            children.entry(critical as u32).or_default().push(en);
        }
        if cyc_children.is_empty() {
            lambda[cur_root] = critical as i64;
        }

        let cs = e_src[critical] as usize;
        let ct = e_tgt[critical] as usize;
        if w_dsu.find(cs) != w_dsu.find(ct) {
            enter[cur_root] = critical as i64;
            w_dsu.union(cs, ct);
        } else {
            // critical_edge closes a cycle; contract it.
            let mut cycle_edges: Vec<u32> = Vec::new();
            let mut cycle_repr: Vec<usize> = Vec::new();
            let mut least_costly = critical;
            enter[cur_root] = NIL;
            cycle_edges.push(critical as u32);
            cycle_repr.push(s.find(ct));
            let mut v = s.find(cs);
            while enter[v] != NIL {
                let ev = enter[v] as usize;
                cycle_edges.push(ev as u32);
                cycle_repr.push(v);
                // minimum-branching: least_costly = the MAX-weight edge in the cycle.
                if e_w[ev] > e_w[least_costly] {
                    least_costly = ev;
                }
                v = s.find(e_src[ev] as usize);
            }
            // reweight offset for edges entering each cycle vertex.
            for &en in &cycle_edges {
                let key = s.find(e_tgt[en as usize] as usize);
                edge_weight_change[key] = e_w[least_costly] - e_w[en as usize];
            }
            let cycle_root = min_[s.find(e_tgt[least_costly] as usize)];
            // union all cycle components into one representative.
            let mut new_repr = cycle_repr[0];
            for &vv in &cycle_repr {
                let (rv, rn) = (s.find(vv), s.find(new_repr));
                s.link(rv, rn);
                new_repr = s.find(new_repr);
            }
            min_[new_repr] = cycle_root;
            roots.push(new_repr);
            cycle[new_repr] = cycle_edges;
            // apply the reweighting to each cycle component's in-edges.
            for &vv in &cycle_repr {
                let delta = edge_weight_change[vv];
                if delta != 0.0 {
                    let list = std::mem::take(&mut in_edges[vv]);
                    for &en in &list {
                        e_w[en as usize] += delta;
                    }
                    in_edges[vv] = list;
                }
            }
            // merge the cycle components' (source-sorted) in-edge lists, dropping
            // now-internal edges (source in new_repr) and keeping the better weight
            // on equal source. Running pairwise merge into cycle_repr[i].
            for i in 1..cycle_repr.len() {
                let l1 = std::mem::take(&mut in_edges[cycle_repr[i]]);
                let l2 = std::mem::take(&mut in_edges[cycle_repr[i - 1]]);
                let mut merged: Vec<u32> = Vec::with_capacity(l1.len() + l2.len());
                let (mut i1, mut i2) = (0usize, 0usize);
                loop {
                    while i1 < l1.len() && s.find(e_src[l1[i1] as usize] as usize) == new_repr {
                        i1 += 1;
                    }
                    while i2 < l2.len() && s.find(e_src[l2[i2] as usize] as usize) == new_repr {
                        i2 += 1;
                    }
                    if i1 == l1.len() && i2 == l2.len() {
                        break;
                    }
                    if i1 == l1.len() {
                        merged.push(l2[i2]);
                        i2 += 1;
                    } else if i2 == l2.len() {
                        merged.push(l1[i1]);
                        i1 += 1;
                    } else {
                        let s1 = e_src[l1[i1] as usize];
                        let s2 = e_src[l2[i2] as usize];
                        if s1 < s2 {
                            merged.push(l1[i1]);
                            i1 += 1;
                        } else if s1 > s2 {
                            merged.push(l2[i2]);
                            i2 += 1;
                        } else {
                            // same source: keep the better (minimum) weight.
                            if e_w[l1[i1] as usize] < e_w[l2[i2] as usize] {
                                merged.push(l1[i1]);
                            } else {
                                merged.push(l2[i2]);
                            }
                            i1 += 1;
                            i2 += 1;
                        }
                    }
                }
                in_edges[cycle_repr[i]] = merged;
            }
            in_edges[new_repr] = std::mem::take(&mut in_edges[*cycle_repr.last().unwrap()]);
            edge_weight_change[new_repr] = 0.0;
        }
    }

    // ---- extract the optimum branching, matching the binary's print order ----
    let mut f_roots: Vec<usize> = Vec::new();
    for &en in &f_list {
        if !parent.contains_key(&en) {
            f_roots.push(en as usize);
        }
    }
    // remove edges entering the final roots.
    for &vtx in &final_roots {
        if lambda[vtx] != NIL {
            remove_from_f(
                lambda[vtx] as usize,
                &mut f_roots,
                &mut removed_f,
                &mut parent,
                &mut children,
            );
        }
    }
    let mut out: Vec<(usize, usize)> = Vec::new();
    while let Some(en) = f_roots.pop() {
        if removed_f.contains(&(en as u32)) {
            continue;
        }
        out.push((e_src[en] as usize, e_tgt[en] as usize));
        let t = e_tgt[en] as usize;
        if lambda[t] != NIL {
            remove_from_f(
                lambda[t] as usize,
                &mut f_roots,
                &mut removed_f,
                &mut parent,
                &mut children,
            );
        }
    }
    out
}

/// `remove_from_F`: mark `en` and all its ancestors removed; newly orphaned
/// children become F-roots. Sparse maps keyed by edge id.
fn remove_from_f(
    mut en: usize,
    f_roots: &mut Vec<usize>,
    removed: &mut std::collections::HashSet<u32>,
    parent: &mut std::collections::HashMap<u32, u32>,
    children: &mut std::collections::HashMap<u32, Vec<u32>>,
) {
    loop {
        removed.insert(en as u32);
        if let Some(kids) = children.remove(&(en as u32)) {
            for c in kids {
                f_roots.push(c as usize);
                parent.remove(&c);
            }
        }
        match parent.get(&(en as u32)) {
            Some(&p) => en = p as usize,
            None => break,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::process::Command;

    /// Format a cost matrix the way `asymmetric_network_exact` writes it and run
    /// the bundled binary, returning its (source, target) output in print order.
    fn run_binary(n: usize, mat: &[f64]) -> Option<Vec<(usize, usize)>> {
        let bin = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("binaries/edmonds-linux");
        if !bin.exists() {
            return None;
        }
        let mut buf = String::new();
        for i in 0..n {
            for j in 0..n {
                if j > 0 {
                    buf.push('\t');
                }
                buf.push_str(&format!("{:.5}", mat[i * n + j]));
            }
            buf.push('\n');
        }
        let path = std::env::temp_dir().join(format!("gt_tofigh_test_{}.list", std::process::id()));
        std::fs::File::create(&path)
            .unwrap()
            .write_all(buf.as_bytes())
            .unwrap();
        let out = Command::new(&bin).arg(&path).output().unwrap();
        let _ = std::fs::remove_file(&path);
        let text = String::from_utf8_lossy(&out.stdout);
        let mut edges = Vec::new();
        for line in text.lines() {
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() >= 3 {
                edges.push((
                    p[0].parse::<f64>().unwrap().trunc() as usize,
                    p[1].parse::<f64>().unwrap().trunc() as usize,
                ));
            }
        }
        Some(edges)
    }

    #[test]
    fn tofigh_matches_binary_random() {
        let mut seed: u64 = 0xC0FF_EE12_3456_789A;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut ran = false;
        for n in 2..=40usize {
            for _ in 0..30 {
                // integer-ish costs with the +0.999995 offset, quantised to %.5f,
                // so both sides see identical values (and plenty of ties).
                let mut mat = vec![0f64; n * n];
                for i in 0..n {
                    for j in 0..n {
                        let base = if i == j { 0.0 } else { (next() % 6) as f64 };
                        let dd = (base as f32 + (1.0_f64 - 0.000005_f64) as f32) as f64;
                        mat[i * n + j] = format!("{:.5}", dd).parse::<f64>().unwrap();
                    }
                }
                let bin = match run_binary(n, &mat) {
                    Some(b) => b,
                    None => return, // binary not present in this environment
                };
                ran = true;
                let ours = optimum_branching_tofigh(n, |i, j| mat[i * n + j]);
                assert_eq!(
                    ours, bin,
                    "tofigh port != binary at n={n} (order-sensitive compare)"
                );
            }
        }
        assert!(ran, "binary never ran");
    }
}
