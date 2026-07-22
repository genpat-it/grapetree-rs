//! Branch recrafting — port of `methods._branch_recraft` + `contemporary`.
//!
//! MSTreeV2's local optimisation: process branches in increasing (distance,
//! endpoint-weights) order, and try to re-attach each branch's source/target to
//! a better-fitting nearby node, gated by a contemporaneity likelihood test.
//! Groups are merged union-find-style as branches commit.

use crate::distance::DistMatrix;
use std::collections::HashMap;

/// Contemporaneity test (`contemporary`): returns true if the two candidate
/// nodes are plausibly contemporaneous given the pairwise allelic distances and
/// number of loci. Faithful numeric port (natural logs, same clamping).
fn contemporary(a0: f64, a1: f64, b: f64, c: f64, n_loci: f64) -> bool {
    let clamp = |x: f64| x.min(n_loci - 0.5).max(0.5);
    let (a0, a1, b, c) = (clamp(a0), clamp(a1), clamp(b), clamp(c));
    if b >= a0 + c && b >= a1 + c {
        return false;
    }
    if b == c {
        return true;
    }
    let s11 = (1.0 - a0 / n_loci).sqrt();
    let s12 = (2.0 * n_loci - b - c) / 2.0 / (n_loci * (n_loci - a0)).sqrt();
    let v = 1.0 - ((n_loci - a1) * (n_loci - c) / n_loci + (n_loci - b)) / 2.0 / n_loci;
    let s21 = 1.0 + a1 * v / (b - 2.0 * n_loci * v);
    let s22 = 1.0 + c * v / (b - 2.0 * n_loci * v);
    let p1 = a0 * (1.0 - s11 * s11).ln()
        + (n_loci - a0) * (s11 * s11).ln()
        + (b + c) * (1.0 - s11 * s12).ln()
        + (2.0 * n_loci - b - c) * (s11 * s12).ln();
    let p2 = a1 * (1.0 - s21).ln()
        + (n_loci - a1) * s21.ln()
        + b * (1.0 - s21 * s22).ln()
        + (n_loci - b) * (s21 * s22).ln()
        + c * (1.0 - s22).ln()
        + (n_loci - c) * s22.ln();
    p1 >= p2
}

#[inline]
fn d(dist: &DistMatrix, i: usize, j: usize) -> f64 {
    dist.get(i, j) as f64
}

/// Port of `_branch_recraft`. `branches` are `(src, tgt, brlen)`; `dist` is the
/// (asymmetric) distance matrix; `weights` the node weights; `n_loci` the loci
/// count. Returns the recrafted branch list.
pub fn branch_recraft(
    mut branches: Vec<(usize, usize, f64)>,
    dist: &DistMatrix,
    weights: &[f64],
    n_loci: f64,
) -> Vec<(usize, usize, f64)> {
    // group/adjacency bookkeeping over all endpoints
    let mut group_id: HashMap<usize, usize> = HashMap::new();
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut childrens: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(s, t, _) in &branches {
        for b in [s, t] {
            group_id.entry(b).or_insert(b);
            groups.entry(b).or_insert_with(|| vec![b]);
            childrens.entry(b).or_default();
        }
    }

    // initial sort by [dist[src][tgt], sorted(weights[src], weights[tgt])]
    let sort_key = |br: &(usize, usize, f64)| -> (f64, f64, f64) {
        let (s, t) = (br.0, br.1);
        let (w0, w1) = (weights[s], weights[t]);
        (d(dist, s, t), w0.min(w1), w0.max(w1))
    };
    branches.sort_by(|a, b| sort_key(a).partial_cmp(&sort_key(b)).unwrap());

    let mut i = 0usize;
    while i < branches.len() {
        let (mut src, mut tgt, _) = branches[i];
        let sources = groups[&group_id[&src]].clone();
        let targets = groups[&group_id[&tgt]].clone();
        let mut tried: HashMap<usize, usize> = HashMap::new();

        // --- try to re-attach the source end ---
        if sources.len() > 1 {
            let mut cand: Vec<(f64, f64, usize)> = sources
                .iter()
                .map(|&s| (weights[s], d(dist, s, tgt), s))
                .collect();
            cand.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for &(_, dd, s) in cand.iter().take(3) {
                if s == src {
                    break;
                }
                if dd < 1.5 * d(dist, src, tgt)
                    && contemporary(
                        d(dist, s, src),
                        d(dist, src, s),
                        dd,
                        d(dist, src, tgt),
                        n_loci,
                    )
                {
                    tried.insert(src, s);
                    src = s;
                    break;
                }
            }
            while !tried.contains_key(&src) {
                tried.insert(src, src);
                let dsrctgt = d(dist, src, tgt);
                let mut mid: Vec<(f64, f64, usize)> = childrens[&src]
                    .iter()
                    .filter(|&&s| !tried.contains_key(&s) && d(dist, s, tgt) < 2.0 * dsrctgt)
                    .map(|&s| (weights[s], d(dist, s, tgt), s))
                    .collect();
                mid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                for &(w, dd, s) in &mid {
                    if dd < d(dist, src, tgt) {
                        if !contemporary(
                            d(dist, src, s),
                            d(dist, s, src),
                            d(dist, src, tgt),
                            dd,
                            n_loci,
                        ) {
                            tried.insert(src, s);
                            src = s;
                            break;
                        }
                    } else if w < weights[src]
                        && contemporary(
                            d(dist, s, src),
                            d(dist, src, s),
                            dd,
                            d(dist, src, tgt),
                            n_loci,
                        )
                    {
                        tried.insert(src, s);
                        src = s;
                        break;
                    }
                    tried.insert(s, src);
                }
            }
        }

        // --- try to re-attach the target end (mirror of the source block) ---
        if targets.len() > 1 {
            let mut cand: Vec<(f64, f64, usize)> = targets
                .iter()
                .map(|&t| (weights[t], d(dist, src, t), t))
                .collect();
            cand.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for &(_, dd, t) in cand.iter().take(3) {
                if t == tgt {
                    break;
                }
                if dd < 1.5 * d(dist, src, tgt)
                    && contemporary(
                        d(dist, t, tgt),
                        d(dist, tgt, t),
                        dd,
                        d(dist, src, tgt),
                        n_loci,
                    )
                {
                    tried.insert(tgt, t);
                    tgt = t;
                    break;
                }
            }
            while !tried.contains_key(&tgt) {
                tried.insert(tgt, tgt);
                let dsrctgt = d(dist, src, tgt);
                let mut mid: Vec<(f64, f64, usize)> = childrens[&tgt]
                    .iter()
                    .filter(|&&t| !tried.contains_key(&t) && d(dist, src, t) < 2.0 * dsrctgt)
                    .map(|&t| (weights[t], d(dist, src, t), t))
                    .collect();
                mid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                for &(w, dd, t) in &mid {
                    if dd < d(dist, src, tgt) {
                        if !contemporary(
                            d(dist, tgt, t),
                            d(dist, t, tgt),
                            d(dist, src, tgt),
                            dd,
                            n_loci,
                        ) {
                            tried.insert(tgt, t);
                            tgt = t;
                            break;
                        }
                    } else if w < weights[tgt]
                        && contemporary(
                            d(dist, t, tgt),
                            d(dist, tgt, t),
                            dd,
                            d(dist, src, tgt),
                            n_loci,
                        )
                    {
                        tried.insert(tgt, t);
                        tgt = t;
                        break;
                    }
                    tried.insert(t, tgt);
                }
            }
        }

        let brlen = d(dist, src, tgt);
        branches[i] = (src, tgt, brlen);

        if i >= branches.len() - 1 || branches[i + 1].2 >= brlen {
            // commit: merge tgt's group into src's; record adjacency
            let gsrc = group_id[&src];
            let tid = group_id[&tgt];
            for &t in &targets {
                group_id.insert(t, gsrc);
            }
            if tid != gsrc {
                if let Some(mut popped) = groups.remove(&tid) {
                    groups.get_mut(&gsrc).unwrap().append(&mut popped);
                }
            }
            childrens.get_mut(&src).unwrap().push(tgt);
            childrens.get_mut(&tgt).unwrap().push(src);
            i += 1;
        } else {
            // not yet the minimum remaining branch: re-sort the tail by brlen
            branches[i..].sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        }
    }
    branches
}
