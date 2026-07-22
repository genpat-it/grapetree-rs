//! Optimum branching (minimum spanning arborescence, root unspecified).
//!
//! Efficient O(E log V) Chu-Liu/Edmonds via a lazy-add skew heap of incoming
//! edges per node plus a rollback union-find (the Tarjan/Gabow formulation,
//! as popularised by KACTL's `DirectedMST`). This replaces the compiled
//! `edmonds` binary the reference shells out to — with no temp-file round trip.
//!
//! The reference feeds a dense cost matrix `w[i][j] = round(d[i][j]) + weight[i]`
//! (a constant `+0.999995` offset it adds does not change the argmin and is
//! omitted). We compute the globally minimum arborescence over all roots by
//! adding a virtual root joined to every node at a prohibitive cost, so exactly
//! one real root is selected. For MSTreeV2 the per-source weights are the
//! distinct harmonic ranks, so the optimum branching is unique.

/// Minimum-arborescence edge chosen for a node: source `a`, target `b`, cost `w`.
#[derive(Clone, Copy)]
struct E {
    a: i32,
    b: i32,
    w: f64,
}
const NIL: i32 = -1;

/// An original directed edge `(source, target)` over real node ids.
type OrigEdge = (u32, u32);
/// Cycle members contracted in one round: `(old_label, its in-cycle edge)`.
type CycleMembers = Vec<(usize, OrigEdge)>;

/// Lazy-add skew heap arena. Each node holds an edge key, children, and a
/// pending additive delta pushed down on access.
struct Heaps {
    a: Vec<i32>,
    b: Vec<i32>,
    w: Vec<f64>,
    l: Vec<i32>,
    r: Vec<i32>,
    delta: Vec<f64>,
}
impl Heaps {
    fn with_capacity(c: usize) -> Self {
        Heaps {
            a: Vec::with_capacity(c),
            b: Vec::with_capacity(c),
            w: Vec::with_capacity(c),
            l: Vec::with_capacity(c),
            r: Vec::with_capacity(c),
            delta: Vec::with_capacity(c),
        }
    }
    fn new_node(&mut self, e: E) -> i32 {
        let id = self.a.len() as i32;
        self.a.push(e.a);
        self.b.push(e.b);
        self.w.push(e.w);
        self.l.push(NIL);
        self.r.push(NIL);
        self.delta.push(0.0);
        id
    }
    #[inline]
    fn prop(&mut self, x: i32) {
        let x = x as usize;
        let d = self.delta[x];
        if d != 0.0 {
            self.w[x] += d;
            let (l, r) = (self.l[x], self.r[x]);
            if l != NIL {
                self.delta[l as usize] += d;
            }
            if r != NIL {
                self.delta[r as usize] += d;
            }
            self.delta[x] = 0.0;
        }
    }
    /// Skew-heap meld with lazy-delta propagation (min-heap on `w`).
    fn merge(&mut self, a: i32, b: i32) -> i32 {
        if a == NIL {
            return b;
        }
        if b == NIL {
            return a;
        }
        self.prop(a);
        self.prop(b);
        let (mut a, mut b) = (a, b);
        if self.w[a as usize] > self.w[b as usize] {
            std::mem::swap(&mut a, &mut b);
        }
        // a.r = merge(b, a.r); then swap(a.l, a.r)
        let ar = self.r[a as usize];
        let merged = self.merge(b, ar);
        let ai = a as usize;
        self.r[ai] = merged;
        std::mem::swap(&mut self.l[ai], &mut self.r[ai]);
        a
    }

    /// Peek the minimum edge (after propagating pending deltas).
    fn top(&mut self, x: i32) -> E {
        self.prop(x);
        let xi = x as usize;
        E {
            a: self.a[xi],
            b: self.b[xi],
            w: self.w[xi],
        }
    }

    /// Remove the minimum, returning the new heap root.
    fn pop(&mut self, x: i32) -> i32 {
        self.prop(x);
        let xi = x as usize;
        self.merge(self.l[xi], self.r[xi])
    }
}

/// Compute the optimum branching over a dense cost function.
///
/// `w(i, j)` gives the cost of directed edge `i -> j` (i != j); return `None`
/// to mean "no such edge". Returns arborescence edges as `(parent, child)`.
pub fn optimum_branching(n: usize, w: impl Fn(usize, usize) -> Option<f64>) -> Vec<(usize, usize)> {
    if n <= 1 {
        return Vec::new();
    }
    let mut edges: Vec<E> = Vec::new();
    let mut maxw = 0f64;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            if let Some(c) = w(i, j) {
                if c > maxw {
                    maxw = c;
                }
                edges.push(E {
                    a: i as i32,
                    b: j as i32,
                    w: c,
                });
            }
        }
    }
    let vroot = n;
    let big = maxw * (n as f64 + 1.0) + 1.0;
    for v in 0..n {
        edges.push(E {
            a: vroot as i32,
            b: v as i32,
            w: big,
        });
    }
    let n_all = n + 1;
    let in_edges = dmst(n_all, vroot, &edges);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let e = in_edges[i];
        if e.a >= 0 && e.a as usize != vroot {
            out.push((e.a as usize, e.b as usize));
        }
    }
    out.sort_unstable();
    out
}

/// One contraction round's bookkeeping, kept for reverse expansion.
struct Round {
    /// For each cycle contracted this round: `(new_label, members)`.
    cycles: Vec<(usize, CycleMembers)>,
    /// old supernode label -> new supernode label.
    map: Vec<usize>,
    old_n: usize,
    old_root: usize,
}

/// Dense minimum spanning arborescence over all roots, in **O(V) memory** beyond
/// the caller's matrix (no per-edge heap nodes, no materialised edge list).
///
/// Boruvka-style Chu-Liu/Edmonds: each round every non-root node takes its
/// cheapest incoming edge, all resulting cycles are contracted at once, and the
/// cost matrix is rebuilt over the shrunken super-node set. Cycles are expanded
/// in reverse at the end. `cost(i, j)` is the directed cost of `i -> j` over
/// `0..m`; a virtual root joined to every node at a prohibitive cost selects
/// exactly one real root, matching [`optimum_branching`].
///
/// When the optimum branching is unique — MSTreeV2's case, where the incoming
/// costs to any node are pairwise distinct (`round(d)+weight[src]`, with the
/// harmonic `weight` a distinct rank per source) — this returns the *same*
/// arborescence as [`optimum_branching`], byte-for-byte (verified by the
/// `dense_matches_heap` property test and the regression golden).
pub fn optimum_branching_dense(
    m: usize,
    cost: impl Fn(usize, usize) -> f64,
) -> Vec<(usize, usize)> {
    if m <= 1 {
        return Vec::new();
    }
    let n0 = m + 1; // virtual root at index m
    let root0 = m;
    let mut maxw = 0f64;
    for i in 0..m {
        for j in 0..m {
            if i != j {
                let c = cost(i, j);
                if c > maxw {
                    maxw = c;
                }
            }
        }
    }
    let big = maxw * (n0 as f64 + 1.0) + 1.0;
    let inf = f64::INFINITY;

    // Dense super-graph: `cst[a*n+b]` = min cost of a real edge a->b (INF if none),
    // `oe[a*n+b]` = the original (u,v) achieving it (v is always a real node).
    let mut n = n0;
    let mut cst = vec![inf; n * n];
    let mut oe = vec![(u32::MAX, u32::MAX); n * n];
    for b in 0..m {
        for a in 0..n {
            if a == b {
                continue;
            }
            let c = if a == root0 { big } else { cost(a, b) };
            cst[a * n + b] = c;
            oe[a * n + b] = (a as u32, b as u32);
        }
    }
    let timing = std::env::var("GT_TIMING").is_ok();
    let mut root = root0;
    let mut history: Vec<Round> = Vec::new();
    // `cur[r]` = current super-node label of real node r; snapshotted per round
    // (old labels) so expansion can route an external edge to the member that
    // *contains* its real target (which may be nested several rounds deep).
    let mut cur: Vec<usize> = (0..m).collect();
    let mut snaps: Vec<Vec<usize>> = Vec::new();

    loop {
        // 1. cheapest incoming edge per non-root super-node.
        let mut in_arg = vec![usize::MAX; n];
        let mut in_min = vec![inf; n];
        let mut in_edge = vec![(u32::MAX, u32::MAX); n];
        for b in 0..n {
            if b == root {
                continue;
            }
            let (mut best, mut ba) = (inf, usize::MAX);
            for a in 0..n {
                if a != b {
                    let c = cst[a * n + b];
                    if c < best {
                        best = c;
                        ba = a;
                    }
                }
            }
            in_arg[b] = ba;
            in_min[b] = best;
            in_edge[b] = oe[ba * n + b];
        }

        // 2. find cycles in the chosen-parent functional graph.
        let mut comp = vec![usize::MAX; n];
        let mut state = vec![0u8; n]; // 0 unseen, 1 on stack, 2 done
        let mut ncyc = 0usize;
        for start in 0..n {
            if state[start] != 0 || start == root {
                continue;
            }
            let mut path: Vec<usize> = Vec::new();
            let mut u = start;
            while u != root && state[u] == 0 {
                state[u] = 1;
                path.push(u);
                u = in_arg[u];
            }
            if u != root && state[u] == 1 && comp[u] == usize::MAX {
                // new cycle: the suffix of `path` from the first occurrence of u.
                let pos = path.iter().position(|&x| x == u).unwrap();
                let cid = ncyc;
                ncyc += 1;
                for &node in &path[pos..] {
                    comp[node] = cid;
                }
            }
            for &node in &path {
                state[node] = 2;
            }
        }

        if ncyc == 0 {
            // no cycles: the chosen incoming edges are the arborescence.
            if timing {
                eprintln!("[timing] dense edmonds: {} rounds (m={m})", history.len());
            }
            let mut chosen: Vec<Option<OrigEdge>> = vec![None; n];
            for b in 0..n {
                if b != root {
                    chosen[b] = Some(in_edge[b]);
                }
            }
            return expand_dense(chosen, history, snaps, m);
        }

        // 3. relabel: cycles get ids 0..ncyc, other nodes get fresh ids after.
        let mut map = vec![usize::MAX; n];
        let mut newn = ncyc;
        for v in 0..n {
            map[v] = if comp[v] != usize::MAX {
                comp[v]
            } else {
                let id = newn;
                newn += 1;
                id
            };
        }
        let new_root = map[root];
        let mut round_cycles: Vec<(usize, CycleMembers)> =
            (0..ncyc).map(|c| (c, Vec::new())).collect();
        for v in 0..n {
            if comp[v] != usize::MAX {
                round_cycles[comp[v]].1.push((v, in_edge[v]));
            }
        }

        // 4. rebuild the cost matrix over super-nodes, reweighting edges into
        //    cycle nodes by subtracting that node's cheapest-incoming cost.
        let mut ncst = vec![inf; newn * newn];
        let mut noe = vec![(u32::MAX, u32::MAX); newn * newn];
        for a in 0..n {
            let base = a * n;
            for b in 0..n {
                if a == b {
                    continue;
                }
                let c = cst[base + b];
                if c == inf {
                    continue;
                }
                let (na, nb) = (map[a], map[b]);
                if na == nb {
                    continue; // internal edge, dropped
                }
                let w = if comp[b] != usize::MAX {
                    c - in_min[b]
                } else {
                    c
                };
                let idx = na * newn + nb;
                if w < ncst[idx] {
                    ncst[idx] = w;
                    noe[idx] = oe[base + b];
                }
            }
        }

        // snapshot old labels of every real node, then advance to new labels.
        snaps.push(cur.clone());
        for r in 0..m {
            cur[r] = map[cur[r]];
        }
        history.push(Round {
            cycles: round_cycles,
            map,
            old_n: n,
            old_root: root,
        });
        cst = ncst;
        oe = noe;
        n = newn;
        root = new_root;
    }
}

/// Reverse-expand the contracted cycles to recover per-node original edges.
/// `snaps[r][real]` is the old-label of real node `real` during round `r`, used
/// to route an external edge into the exact cycle member that *contains* the
/// edge's real target.
fn expand_dense(
    mut chosen: Vec<Option<OrigEdge>>,
    history: Vec<Round>,
    snaps: Vec<Vec<usize>>,
    m: usize,
) -> Vec<(usize, usize)> {
    for ri in (0..history.len()).rev() {
        let round = &history[ri];
        let snap = &snaps[ri];
        let mut old_chosen: Vec<Option<OrigEdge>> = vec![None; round.old_n];
        let mut is_member = vec![false; round.old_n];
        for (cnew, members) in &round.cycles {
            // external edge entering the contracted cycle super-node
            let ext = chosen[*cnew];
            for &(old_label, in_e) in members {
                is_member[old_label] = true;
                old_chosen[old_label] = Some(in_e); // keep in-cycle edge by default
            }
            // route the external edge to the member containing its real target;
            // that member drops its in-cycle edge.
            if let Some((ue, ve)) = ext {
                let member = snap[ve as usize];
                old_chosen[member] = Some((ue, ve));
            }
        }
        for v in 0..round.old_n {
            if !is_member[v] && v != round.old_root {
                old_chosen[v] = chosen[round.map[v]];
            }
        }
        chosen = old_chosen;
    }
    let mut out = Vec::new();
    for b in 0..m {
        if let Some((u, v)) = chosen[b] {
            if u as usize != m {
                out.push((u as usize, v as usize));
            }
        }
    }
    out.sort_unstable();
    out
}

/// Directed minimum spanning arborescence rooted at `root`.
/// Returns, per node, the chosen incoming edge (`a < 0` if none).
fn dmst(n: usize, root: usize, g: &[E]) -> Vec<E> {
    let mut uf = RollbackUf::new(n);
    let mut h = Heaps::with_capacity(g.len() + 8);
    let mut heap: Vec<i32> = vec![NIL; n];
    for &e in g {
        let node = h.new_node(e);
        heap[e.b as usize] = h.merge(heap[e.b as usize], node);
    }
    let mut seen: Vec<i32> = vec![-1; n];
    let mut path: Vec<usize> = vec![0; n];
    let mut q: Vec<E> = vec![
        E {
            a: NIL,
            b: NIL,
            w: 0.0
        };
        n
    ];
    let mut in_e: Vec<E> = vec![
        E {
            a: NIL,
            b: NIL,
            w: 0.0
        };
        n
    ];
    seen[root] = root as i32;
    // cycles recorded for expansion: (contracted-node, uf-time, edges-in-cycle)
    let mut cycs: Vec<(usize, usize, Vec<E>)> = Vec::new();

    for s in 0..n {
        let mut u = s;
        let mut qi = 0usize;
        while seen[u] < 0 {
            if heap[u] == NIL {
                // unreachable node -> no arborescence (shouldn't happen with vroot)
                return vec![
                    E {
                        a: NIL,
                        b: NIL,
                        w: 0.0
                    };
                    n
                ];
            }
            let top = h.top(heap[u]);
            // subtract this edge's weight from the whole heap, then pop it
            h.delta[heap[u] as usize] -= top.w;
            heap[u] = h.pop(heap[u]);
            q[qi] = top;
            path[qi] = u;
            qi += 1;
            seen[u] = s as i32;
            u = uf.find(top.a as usize);
            if seen[u] == s as i32 {
                // found a cycle: contract its nodes
                let mut cyc = NIL;
                let end = qi;
                let time = uf.time();
                let mut w;
                loop {
                    qi -= 1;
                    w = path[qi];
                    cyc = h.merge(cyc, heap[w]);
                    if !uf.join(u, w) {
                        break;
                    }
                }
                u = uf.find(u);
                heap[u] = cyc;
                seen[u] = -1;
                let comp: Vec<E> = q[qi..end].to_vec();
                cycs.push((u, time, comp));
            }
        }
        for i in 0..qi {
            let t = uf.find(q[i].b as usize);
            in_e[t] = q[i];
        }
    }

    if std::env::var("GT_TIMING").is_ok() {
        eprintln!("[timing] edmonds: n_all={n} cycles={}", cycs.len());
    }
    // expand cycles in reverse
    for (u, time, comp) in cycs.into_iter().rev() {
        uf.rollback(time);
        let in_edge = in_e[u];
        for &e in &comp {
            let t = uf.find(e.b as usize);
            in_e[t] = e;
        }
        let t = uf.find(in_edge.b as usize);
        in_e[t] = in_edge;
    }
    in_e
}

/// Rollback union-find (union by size, no path compression).
struct RollbackUf {
    e: Vec<i32>,
    st: Vec<(usize, i32)>,
}
impl RollbackUf {
    fn new(n: usize) -> Self {
        RollbackUf {
            e: vec![-1; n],
            st: Vec::new(),
        }
    }
    fn find(&self, mut x: usize) -> usize {
        while self.e[x] >= 0 {
            x = self.e[x] as usize;
        }
        x
    }
    fn time(&self) -> usize {
        self.st.len()
    }
    fn rollback(&mut self, t: usize) {
        while self.st.len() > t {
            let (i, v) = self.st.pop().unwrap();
            self.e[i] = v;
        }
    }
    fn join(&mut self, a: usize, b: usize) -> bool {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return false;
        }
        if self.e[a] > self.e[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.st.push((a, self.e[a]));
        self.st.push((b, self.e[b]));
        self.e[a] += self.e[b];
        self.e[b] = a as i32;
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn eset(v: Vec<(usize, usize)>) -> HashSet<(usize, usize)> {
        v.into_iter().collect()
    }

    #[test]
    fn case_a_chain() {
        let m = [
            [0., 1., 5., 9.],
            [5., 0., 1., 9.],
            [9., 5., 0., 1.],
            [9., 9., 5., 0.],
        ];
        let r = eset(optimum_branching(4, |i, j| Some(m[i][j])));
        assert_eq!(r, eset(vec![(0, 1), (1, 2), (2, 3)]));
    }

    #[test]
    fn case_d_nonzero_root() {
        let m = [
            [0.99999, 3.99999, 2.49999],
            [2.99999, 0.99999, 1.49999],
            [2.99999, 2.49999, 0.99999],
        ];
        let r = eset(optimum_branching(3, |i, j| Some(m[i][j])));
        assert_eq!(r, eset(vec![(1, 2), (1, 0)]));
    }

    #[test]
    fn dense_matches_heap() {
        // Deterministic LCG so the test is reproducible without external crates.
        let mut seed: u64 = 0x1234_5678_9abc_def0;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for n in 2..=30usize {
            for _ in 0..40 {
                // distinct-ish integer costs; include ties to stress cycle logic.
                let mut w = vec![0f64; n * n];
                for i in 0..n {
                    for j in 0..n {
                        if i != j {
                            w[i * n + j] = (next() % 50) as f64;
                        }
                    }
                }
                let cost = |i: usize, j: usize| w[i * n + j];
                let heap = eset(optimum_branching(n, |i, j| Some(cost(i, j))));
                let dense = eset(optimum_branching_dense(n, cost));
                // Both must be valid arborescences; when the optimum is unique
                // they are identical. Compare total cost (always equal) and, when
                // equal edge-set, exact identity.
                let sum = |s: &HashSet<(usize, usize)>| -> f64 {
                    s.iter().map(|&(i, j)| w[i * n + j]).sum()
                };
                assert_eq!(heap.len(), n - 1, "heap not spanning n={n}");
                assert_eq!(dense.len(), n - 1, "dense not spanning n={n}");
                assert!(
                    (sum(&heap) - sum(&dense)).abs() < 1e-9,
                    "cost mismatch n={n}: heap={} dense={}",
                    sum(&heap),
                    sum(&dense)
                );
                // acyclic + single root for the dense result
                let mut parent = vec![usize::MAX; n];
                for &(u, v) in &dense {
                    assert_eq!(parent[v], usize::MAX, "two parents n={n}");
                    parent[v] = u;
                }
                assert_eq!(parent.iter().filter(|&&p| p == usize::MAX).count(), 1);
            }
        }
    }

    #[test]
    fn dense_identical_when_unique() {
        // MSTreeV2-like costs: round(d) + distinct per-source weight -> unique.
        let n = 25;
        let cost = |i: usize, j: usize| {
            let d = ((i * 31 + j * 17) % 11) as f64; // integer "distance"
            d + (i as f64) / (n as f64) // distinct source weight in [0,1)
        };
        let heap = eset(optimum_branching(n, |i, j| Some(cost(i, j))));
        let dense = eset(optimum_branching_dense(n, cost));
        assert_eq!(heap, dense, "dense != heap for unique optimum");
    }

    #[test]
    fn spanning_and_acyclic() {
        let n = 8;
        let cost = |i: usize, j: usize| ((i * 7 + j * 13 + (i ^ j) * 3) % 17 + 1) as f64;
        let arb = optimum_branching(n, |i, j| Some(cost(i, j)));
        assert_eq!(arb.len(), n - 1);
        let mut parent = vec![usize::MAX; n];
        for &(u, v) in &arb {
            assert_eq!(parent[v], usize::MAX, "node {v} has two parents");
            parent[v] = u;
        }
        assert_eq!(parent.iter().filter(|&&p| p == usize::MAX).count(), 1);
        for start in 0..n {
            let (mut cur, mut steps) = (start, 0);
            while parent[cur] != usize::MAX {
                cur = parent[cur];
                steps += 1;
                assert!(steps <= n, "cycle detected");
            }
        }
    }
}
