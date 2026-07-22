//! Branch recrafting — port of `methods._branch_recraft` + `contemporary`.
//!
//! MSTreeV2's local optimisation: process branches in increasing (distance,
//! endpoint-weights) order, and try to re-attach each branch's source/target to
//! a better-fitting nearby node, gated by a contemporaneity likelihood test.
//! Groups are merged union-find-style as branches commit.
//!
//! Performance note (GRAPETREE-COMPAT[recraft-perf]): the algorithm's greedy
//! spine is inherently sequential (each branch's decision depends on the merged
//! group state left by all prior commits, and the processing order itself is
//! re-derived as branch lengths change). What is *not* inherent is the original
//! port's per-iteration O(group) overhead. Three output-preserving speedups vs a
//! naive transcription of the reference:
//!   1. the endpoint groups are read by reference, never cloned;
//!   2. `group_id` / `groups` / `childrens` are `Vec`s indexed by node id (the
//!      endpoints are dense indices `0..n`), so the reference's numpy fancy
//!      assignment `group_id[targets] = group_id[src]` is an O(group) array write
//!      rather than O(group) hashed inserts;
//!   3. the "3 nearest candidates" are taken with a partial selection instead of
//!      a full sort. The candidate order is a *total* order (the tuple's last key
//!      is the unique node id), so selecting the 3 smallest is unambiguous and
//!      byte-identical to `sorted(...)[:3]`.
//!
//! None of these change which branches are chosen or in what order — the output
//! NEWICK is bit-identical to the reference (verified by the regression golden).

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

/// Reorder `cand` so its first `min(3, len)` entries are the three smallest, in
/// exact ascending order. The comparison is a *total* order (the final tuple
/// element is a unique node id), so this is byte-identical to a full stable sort
/// followed by `take(3)` — only cheaper (O(n) selection vs O(n log n) sort).
/// Entries past index 3 are left unspecified; the callers only read `take(3)`.
#[inline]
fn top3(cand: &mut [(f64, f64, usize)]) {
    let k = cand.len().min(3);
    if cand.len() > k {
        cand.select_nth_unstable_by(k - 1, |a, b| a.partial_cmp(b).unwrap());
    }
    cand[..k].sort_by(|a, b| a.partial_cmp(b).unwrap());
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
    // group/adjacency bookkeeping over all endpoints. Endpoints are dense node
    // ids in `0..n` (a spanning tree over `n` nodes touches them all), so `Vec`s
    // indexed by id replace the reference's dict/numpy-array bookkeeping without
    // changing behaviour — just dropping the hashing overhead.
    let n = weights.len();
    let mut group_id: Vec<usize> = (0..n).collect();
    let mut groups: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut childrens: Vec<Vec<usize>> = vec![Vec::new(); n];

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
        // group ids of the *original* endpoints (captured before re-attachment,
        // exactly like the reference snapshots `sources`/`targets` up front).
        let src_gid = group_id[src];
        let tgt_gid = group_id[tgt];
        let mut tried: HashMap<usize, usize> = HashMap::new();

        // --- try to re-attach the source end ---
        if groups[src_gid].len() > 1 {
            let mut cand: Vec<(f64, f64, usize)> = groups[src_gid]
                .iter()
                .map(|&s| (weights[s], d(dist, s, tgt), s))
                .collect();
            top3(&mut cand);
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
                let mut mid: Vec<(f64, f64, usize)> = childrens[src]
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
        if groups[tgt_gid].len() > 1 {
            let mut cand: Vec<(f64, f64, usize)> = groups[tgt_gid]
                .iter()
                .map(|&t| (weights[t], d(dist, src, t), t))
                .collect();
            top3(&mut cand);
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
                let mut mid: Vec<(f64, f64, usize)> = childrens[tgt]
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
            // commit: merge tgt's group into src's; record adjacency.
            let gsrc = group_id[src];
            let tid = group_id[tgt];
            // Re-point every member of the *original* target group to gsrc
            // (reference: `group_id[targets] = group_id[src]`). `group_id` and
            // `groups` are distinct `Vec`s, so we can read the member list while
            // writing ids. `groups[tgt_gid]` is untouched until the take below,
            // so it still holds the pre-commit snapshot the reference used.
            let m = groups[tgt_gid].len();
            for idx in 0..m {
                let t = groups[tgt_gid][idx];
                group_id[t] = gsrc;
            }
            if tid != gsrc {
                let mut popped = std::mem::take(&mut groups[tid]);
                groups[gsrc].append(&mut popped);
            }
            childrens[src].push(tgt);
            childrens[tgt].push(src);
            i += 1;
        } else {
            // not yet the minimum remaining branch: re-sort the tail by brlen
            branches[i..].sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        }
    }
    branches
}
