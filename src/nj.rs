//! Neighbor-Joining (Saitou & Nei 1987) — a native replacement for the
//! FastME / RapidNJ / Ninja binaries the reference shells out to for the
//! `NJ` / `RapidNJ` / `ninja` methods.
//!
//! Byte-identity with those third-party tools is not a goal (different
//! implementations, tie-breaking and rooting). We build the canonical NJ tree
//! on GrapeTree's symmetric distance matrix and emit an unrooted NEWICK; the
//! reference's subsequent midpoint-rooting/unrooting does not change the
//! unrooted topology or branch lengths, so regression is scored on the
//! unrooted tree (RF + branch lengths).

use crate::distance::DistMatrix;
use crate::parse::Parsed;
use std::collections::HashMap;

struct NjNode {
    // NEWICK subtree string built so far for this active cluster
    repr: String,
    is_leaf: bool,
}

/// Build the canonical NJ tree over the (symmetric) distance matrix and return
/// an unrooted NEWICK string. Leaf labels are the representative names; caller
/// applies embeded expansion.
pub fn neighbor_joining(p: &Parsed, dist: &DistMatrix, names: &[String]) -> String {
    let n = dist.n;
    if n == 0 {
        return ";".to_string();
    }
    if n == 1 {
        return format!("{};", expand_leaf(&names[0], p));
    }

    // active distance matrix as f64, indexed by active-slot
    let mut d: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| dist.get(i, j) as f64).collect())
        .collect();
    let mut active: Vec<usize> = (0..n).collect(); // slot -> node id
    let mut nodes: HashMap<usize, NjNode> = (0..n)
        .map(|i| {
            (
                i,
                NjNode {
                    repr: expand_leaf(&names[i], p),
                    is_leaf: true,
                },
            )
        })
        .collect();
    let mut next_id = n;
    // d is indexed by slot; map node-id -> slot.
    let mut slot_of: HashMap<usize, usize> =
        active.iter().enumerate().map(|(s, &id)| (id, s)).collect();

    while active.len() > 2 {
        let r = active.len();
        // row sums S[id]
        let s: HashMap<usize, f64> = active
            .iter()
            .map(|&id| {
                let si = active
                    .iter()
                    .map(|&j| d[slot_of[&id]][slot_of[&j]])
                    .sum::<f64>();
                (id, si)
            })
            .collect();
        // find min Q(i,j) = (r-2) d(i,j) - S[i] - S[j]
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..active.len() {
            for b in (a + 1)..active.len() {
                let (i, j) = (active[a], active[b]);
                let q = (r as f64 - 2.0) * d[slot_of[&i]][slot_of[&j]] - s[&i] - s[&j];
                match best {
                    Some((bq, _, _)) if bq <= q => {}
                    _ => best = Some((q, i, j)),
                }
            }
        }
        let (_, i, j) = best.unwrap();
        let dij = d[slot_of[&i]][slot_of[&j]];
        // branch lengths from the new node u to i and j
        let li = 0.5 * dij + (s[&i] - s[&j]) / (2.0 * (r as f64 - 2.0));
        let lj = dij - li;

        let u = next_id;
        next_id += 1;
        let ni = nodes.remove(&i).unwrap();
        let nj = nodes.remove(&j).unwrap();
        let repr = format!("({}:{},{}:{})", ni.repr, fmt(li), nj.repr, fmt(lj));
        nodes.insert(
            u,
            NjNode {
                repr,
                is_leaf: false,
            },
        );

        // new distances d(u,k) = 0.5 (d(i,k)+d(j,k)-d(i,j))
        // reuse slot of i for u; drop slot of j
        let si = slot_of[&i];
        let sj = slot_of[&j];
        for &k in &active {
            if k == i || k == j {
                continue;
            }
            let sk = slot_of[&k];
            let duk = 0.5 * (d[si][sk] + d[sj][sk] - dij);
            d[si][sk] = duk;
            d[sk][si] = duk;
        }
        // update active list + slot map
        active.retain(|&x| x != i && x != j);
        active.push(u);
        slot_of.remove(&i);
        slot_of.remove(&j);
        slot_of.insert(u, si);
        let _ = sj;
    }

    // join the final two clusters
    match active.as_slice() {
        [a, b] => {
            let dab = d[slot_of[a]][slot_of[b]];
            let na = nodes.remove(a).unwrap();
            let nb = nodes.remove(b).unwrap();
            // unrooted: connect a and b with a single branch; emit as a
            // top-level clade so the tree is unrooted (b's subtree hangs off a).
            if na.is_leaf && nb.is_leaf {
                format!("({}:{},{}:0);", na.repr, fmt(dab), nb.repr)
            } else {
                format!("({}:{},{});", na.repr, fmt(dab), nb.repr)
            }
        }
        [a] => format!("({});", nodes.remove(a).unwrap().repr),
        _ => unreachable!("NJ terminates with 1 or 2 active clusters"),
    }
}

/// Expand an embeded group into a leaf or multifurcation, mirroring backend().
fn expand_leaf(name: &str, p: &Parsed) -> String {
    if let Some(group) = p.embeded.iter().find(|g| g[0] == name) {
        if group.len() > 1 {
            let kids: Vec<String> = group.iter().map(|m| format!("{m}:0")).collect();
            return format!("({})", kids.join(","));
        }
    }
    name.to_string()
}

/// Format a branch length like ete3's `%0.6g`.
fn fmt(x: f64) -> String {
    if x.is_finite() && x == x.trunc() && x.abs() < 1e15 {
        return format!("{}", x as i64);
    }
    let mut s = format!("{:.6}", x);
    while s.contains('.') && (s.ends_with('0') || s.ends_with('.')) {
        s.pop();
    }
    s
}
