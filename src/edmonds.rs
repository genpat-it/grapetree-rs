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
