//! Tiebreak weight heuristics — port of `distance_matrix.harmonic` / `.eBurst`.
//!
//! Both assign each node a rank in `[0, 1)` (`rank / N`); lower weight = higher
//! priority as an MST/arborescence parent-or-root. The subtlety is the exact
//! `np.lexsort` ordering, reproduced here with a stable sort over an explicit
//! priority key (highest-priority component first).

use crate::distance::DistMatrix;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Heuristic {
    Harmonic,
    EBurst,
}

impl Heuristic {
    pub fn parse(s: &str) -> Heuristic {
        match s {
            "harmonic" => Heuristic::Harmonic,
            "eBurst" => Heuristic::EBurst,
            other => panic!("unknown heuristic {other:?}"),
        }
    }
}

/// Stable argsort of `0..n` by the key `key(i)` (ascending), ties by index.
fn stable_argsort<K: PartialOrd + Copy>(n: usize, key: impl Fn(usize) -> Vec<K>) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        let (ka, kb) = (key(a), key(b));
        for (x, y) in ka.iter().zip(kb.iter()) {
            match x.partial_cmp(y).unwrap() {
                std::cmp::Ordering::Equal => continue,
                o => return o,
            }
        }
        std::cmp::Ordering::Equal // stable: preserves index order
    });
    idx
}

/// Compute node weights for the given heuristic.
///
/// `n_str[i]` is the number of original strains collapsed into row `i`.
pub fn weights(m: &DistMatrix, n_str: &[usize], heuristic: Heuristic) -> Vec<f64> {
    match heuristic {
        Heuristic::Harmonic => harmonic(m, n_str),
        Heuristic::EBurst => eburst(m, n_str),
    }
}

/// `harmonic`: sort by (harmonic centrality asc, n_str desc); weight = rank/N.
fn harmonic(m: &DistMatrix, n_str: &[usize]) -> Vec<f64> {
    let n = m.n;
    let nf = n as f64;
    // h[i] = N / sum_j 1/(d[i,j] + 0.1)
    let h: Vec<f64> = (0..n)
        .map(|i| {
            let s: f64 = (0..n).map(|j| 1.0 / (m.get(i, j) as f64 + 0.1)).sum();
            nf / s
        })
        .collect();
    // lexsort keys = [-n_str, h] -> primary h asc, secondary -n_str asc.
    let order = stable_argsort(n, |i| vec![h[i], -(n_str[i] as f64)]);
    let mut w = vec![0f64; n];
    for (rank, &node) in order.iter().enumerate() {
        // GRAPETREE-COMPAT[harmonic-f32]: the reference stores `rank/N` as
        // float32 (numpy keeps the harmonic weights in float32), so we round the
        // rank/N through f32. The ~3e-8 difference vs f64 perturbs the reduced
        // arborescence matrix and flips edmonds tie-breaks — required for
        // bit-identity. See COMPAT.md.
        w[node] = (rank as f32 / nf as f32) as f64;
    }
    w
}

/// `eBurst`: goeBURST-style ordering by neighbour-distance histogram.
///
/// For each node, bin the (truncated-integer) row distances. `count[0]` is
/// augmented by `n_str`. Priority (high→low): count[1] desc, count[2] desc, …,
/// count[maxv] desc, then count[0] desc. weight = rank/N.
fn eburst(m: &DistMatrix, n_str: &[usize]) -> Vec<f64> {
    let n = m.n;
    let nf = n as f64;
    // truncate-to-int distances; find global max
    let maxv: usize = (0..n)
        .flat_map(|i| (0..n).map(move |j| m.get(i, j) as i64))
        .max()
        .unwrap_or(0)
        .max(0) as usize;

    // per-node histogram of counts over 0..=maxv
    let mut hist: Vec<Vec<i64>> = vec![vec![0i64; maxv + 1]; n];
    for i in 0..n {
        for j in 0..n {
            let v = m.get(i, j) as i64;
            if v >= 0 {
                hist[i][v as usize] += 1;
            }
        }
        hist[i][0] += n_str[i] as i64;
    }

    // priority key (ascending): [-count1, -count2, ..., -count_maxv, -count0]
    let order = stable_argsort(n, |i| {
        let mut k: Vec<f64> = Vec::with_capacity(maxv + 1);
        for v in 1..=maxv {
            k.push(-(hist[i][v] as f64));
        }
        k.push(-(hist[i][0] as f64));
        k
    });

    let mut w = vec![0f64; n];
    for (rank, &node) in order.iter().enumerate() {
        w[node] = rank as f64 / nf;
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::DistMatrix;

    fn mat(rows: &[&[f32]]) -> DistMatrix {
        let n = rows.len();
        let mut data = vec![0f32; n * n];
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                data[i * n + j] = v;
            }
        }
        DistMatrix { n, data }
    }

    // Oracle values captured from reference for prof [[1,1,1,1,1],[1,2,1,1,1],
    // [1,2,2,1,1],[3,2,2,1,1],[1,1,1,2,2]], n_str=[1,3,1,2,1].
    #[test]
    fn harmonic_matches_oracle() {
        let m = mat(&[
            &[0.0, 1.008, 2.006, 3.004, 2.006],
            &[1.008, 0.0, 1.008, 2.006, 3.004],
            &[2.006, 1.008, 0.0, 1.008, 4.002],
            &[3.004, 2.006, 1.008, 0.0, 5.0],
            &[2.006, 3.004, 4.002, 5.0, 0.0],
        ]);
        let w = weights(&m, &[1, 3, 1, 2, 1], Heuristic::Harmonic);
        // weights are stored via f32 (reference parity), so compare approximately
        let expect = [0.4, 0.0, 0.2, 0.6, 0.8];
        for (a, b) in w.iter().zip(expect.iter()) {
            assert!((a - b).abs() < 1e-6, "got {a}, want {b}");
        }
    }

    #[test]
    fn eburst_matches_oracle() {
        let m = mat(&[
            &[0.0, 1.0, 2.0, 3.0, 2.0],
            &[1.0, 0.0, 1.0, 2.0, 3.0],
            &[2.0, 1.0, 0.0, 1.0, 4.0],
            &[3.0, 2.0, 1.0, 0.0, 5.0],
            &[2.0, 3.0, 4.0, 5.0, 0.0],
        ]);
        let w = weights(&m, &[1, 3, 1, 2, 1], Heuristic::EBurst);
        assert_eq!(w, vec![0.4, 0.0, 0.2, 0.6, 0.8]);
    }
}
