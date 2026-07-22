//! Tree assembly + NEWICK output (replacing ete3).
//!
//! `network2tree` orients an MST/arborescence edge list into a rooted tree and
//! attaches, for every internal node, a zero-length pendant leaf carrying that
//! node's own strain (the reference's self-child trick). Post-processing then
//! collapses tiny branches and expands `embeded` duplicate groups, and the tree
//! is serialised as NEWICK `format=1` (internal node names + branch lengths,
//! `%g`-style numbers) exactly like ete3's writer.

use crate::parse::Parsed;
use std::collections::HashMap;

#[derive(Debug, Clone)]
struct TNode {
    /// Profile index this node represents (before naming); `usize::MAX` once
    /// the node has been renamed to a literal string with no profile identity.
    pid: usize,
    name: String,
    dist: f64,
    children: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<TNode>,
    root: usize,
}

impl Tree {
    fn add(&mut self, pid: usize, dist: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(TNode {
            pid,
            name: String::new(),
            dist,
            children: Vec::new(),
        });
        idx
    }

    /// Build a tree from a *directed-agnostic* weighted edge list.
    ///
    /// Port of `_network2tree`: sort edges by weight descending (stable), orient
    /// them outward from the first edge's first endpoint, then materialise.
    pub fn network2tree(edges: &[(usize, usize, u32)], names: &[String]) -> Tree {
        if edges.is_empty() {
            // Single node (or none): emit a lone leaf if there is one name.
            let mut t = Tree {
                nodes: Vec::new(),
                root: 0,
            };
            let r = t.add(0, 0.0);
            t.nodes[r].name = names.first().cloned().unwrap_or_default();
            t.root = r;
            return t;
        }

        // sort by weight desc, ties keep original order (Python reverse=True stable)
        let mut sorted: Vec<(usize, usize, u32)> = edges.to_vec();
        sorted.sort_by_key(|b| std::cmp::Reverse(b.2));

        // orient outward
        let root_pid = sorted[0].0;
        let mut in_use: HashMap<usize, bool> = HashMap::new();
        in_use.insert(root_pid, true);
        let mut oriented: Vec<(usize, usize, u32)> = Vec::with_capacity(sorted.len());
        let mut pending = sorted;
        while !pending.is_empty() {
            let mut remain = Vec::new();
            for (a, b, w) in pending.into_iter() {
                if in_use.contains_key(&a) {
                    oriented.push((a, b, w));
                    in_use.insert(b, true);
                } else if in_use.contains_key(&b) {
                    oriented.push((b, a, w));
                    in_use.insert(a, true);
                } else {
                    remain.push((a, b, w));
                }
            }
            pending = remain;
        }

        // materialise
        let mut t = Tree {
            nodes: Vec::new(),
            root: 0,
        };
        let root = t.add(root_pid, 0.0);
        t.root = root;
        let mut node_of: HashMap<usize, usize> = HashMap::new();
        node_of.insert(root_pid, root);
        for (src, tgt, w) in &oriented {
            let parent = node_of[src];
            let child = t.add(*tgt, *w as f64);
            t.nodes[parent].children.push(child);
            node_of.insert(*tgt, child);
        }

        // pendant self-child for internal nodes; name leaves.
        let count = t.nodes.len();
        for i in 0..count {
            if t.nodes[i].children.is_empty() {
                let pid = t.nodes[i].pid;
                t.nodes[i].name = names[pid].clone();
            } else {
                let pid = t.nodes[i].pid;
                t.nodes[i].name = String::new();
                let pendant = t.add(pid, 0.0);
                t.nodes[pendant].name = names[pid].clone();
                t.nodes[pendant].pid = usize::MAX;
                t.nodes[i].children.push(pendant);
            }
        }
        t
    }

    /// Post-process like `backend()`: collapse tiny fractional branches when the
    /// tree spans more than 3 units, then expand `embeded` duplicate groups.
    pub fn post_process(&mut self, p: &Parsed) {
        // max branch length among non-root nodes
        let max_dist = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != self.root)
            .map(|(_, n)| n.dist)
            .fold(0.0f64, f64::max);

        if max_dist > 3.0 {
            self.collapse_tiny_branches();
        }

        // expand embeded groups: leaf whose name is a representative with a
        // group size > 1 becomes an empty node with one child per member.
        let group_of: HashMap<&str, &Vec<String>> =
            p.embeded.iter().map(|g| (g[0].as_str(), g)).collect();

        let leaves: Vec<usize> = (0..self.nodes.len())
            .filter(|&i| self.nodes[i].children.is_empty())
            .collect();
        for leaf in leaves {
            let name = self.nodes[leaf].name.clone();
            if let Some(group) = group_of.get(name.as_str()) {
                if group.len() > 1 {
                    self.nodes[leaf].name = String::new();
                    for member in group.iter() {
                        let c = self.add(usize::MAX, 0.0);
                        self.nodes[c].name = member.clone();
                        self.nodes[leaf].children.push(c);
                    }
                }
            }
        }
    }

    /// Collapse branches with `0 < dist < 0.1` into their siblings (reference
    /// behaviour when `maxDist > 3`). No-op for integer branch lengths.
    fn collapse_tiny_branches(&mut self) {
        // parent lookup
        let mut parent = vec![usize::MAX; self.nodes.len()];
        for i in 0..self.nodes.len() {
            for &c in &self.nodes[i].children.clone() {
                parent[c] = i;
            }
        }
        // postorder over descendants
        let order = self.postorder();
        for node in order {
            if node == self.root {
                continue;
            }
            let d = self.nodes[node].dist;
            if d > 0.0 && d < 0.1 {
                let par = parent[node];
                if par != usize::MAX {
                    let sibs: Vec<usize> = self.nodes[par]
                        .children
                        .iter()
                        .copied()
                        .filter(|&c| c != node)
                        .collect();
                    for s in sibs {
                        self.nodes[s].dist += d;
                    }
                    self.nodes[node].dist = 0.0;
                }
            }
        }
    }

    fn postorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![(self.root, false)];
        while let Some((n, processed)) = stack.pop() {
            if processed {
                out.push(n);
            } else {
                stack.push((n, true));
                for &c in self.nodes[n].children.iter().rev() {
                    stack.push((c, false));
                }
            }
        }
        out
    }

    /// Serialise to NEWICK `format=1` (ete3-compatible), trailing `;`.
    pub fn to_newick(&self) -> String {
        let mut s = String::new();
        self.write_node(self.root, true, &mut s);
        s.push(';');
        s
    }

    fn write_node(&self, node: usize, is_root: bool, out: &mut String) {
        let n = &self.nodes[node];
        if n.children.is_empty() {
            out.push_str(&n.name);
            out.push(':');
            out.push_str(&fmt_dist(n.dist));
        } else {
            out.push('(');
            for (k, &c) in n.children.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                self.write_node(c, false, out);
            }
            out.push(')');
            out.push_str(&n.name);
            if !is_root {
                out.push(':');
                out.push_str(&fmt_dist(n.dist));
            }
        }
    }
}

/// Format a branch length like C `printf("%0.6g")` (ete3's default), so
/// integer-valued distances print without a decimal point.
fn fmt_dist(x: f64) -> String {
    if x.is_finite() && x == x.trunc() && x.abs() < 1e15 {
        return format!("{}", x as i64);
    }
    // 6 significant digits, strip trailing zeros (approx %g for our value range)
    let mut s = format!("{:.6}", x);
    while s.contains('.') && (s.ends_with('0') || s.ends_with('.')) {
        s.pop();
    }
    s
}
