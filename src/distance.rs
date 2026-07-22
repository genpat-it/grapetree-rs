//! Distance matrices — faithful port of `distance_matrix.*` (MSTrees.py).
//!
//! Four matrix kinds (`symmetric`, `asymmetric`, `asymmetric_wgMLST`,
//! `blockwise`) × four missing handlers (`pair_delete`, `complete_delete`,
//! `as_allele`, `absolute_distance`). Values are stored `f32` to mirror NumPy's
//! `float32` matrices exactly. Per-cell locus sums run in a fixed order so
//! results are identical regardless of thread count (rayon parallelises over
//! the outer index only).

use crate::parse::Parsed;
use crate::HandleMissing;
use rayon::prelude::*;

/// Which matrix formulation to build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatrixKind {
    Symmetric,
    Asymmetric,
    AsymmetricWgMlst,
    Blockwise,
}

impl MatrixKind {
    /// Resolve from the CLI `matrix_type` string plus the `wgMLST` flag.
    ///
    /// UPSTREAM PARITY NOTE: in reference `backend()` the line
    /// `if wgMLST and matrix=='asymmetric': matrix_type = 'asymmetric_wgMLST'`
    /// assigns a *local* variable that is never propagated to `params`, so
    /// `asymmetric_wgMLST` is unreachable and `--wgMLST` is effectively a no-op
    /// (verified: its `distance` output is byte-identical to plain asymmetric).
    /// We reproduce that behaviour so regression stays exact. The corrected
    /// wgMLST matrix is implemented (see [`MatrixKind::AsymmetricWgMlst`] /
    /// [`asymmetric_wgmlst_into`]) and can be selected explicitly, but is not
    /// produced by this default resolver.
    pub fn resolve(matrix_type: &str, _wg_mlst: bool) -> MatrixKind {
        match matrix_type {
            "symmetric" => MatrixKind::Symmetric,
            "asymmetric" => MatrixKind::Asymmetric,
            "blockwise" => MatrixKind::Blockwise,
            other => panic!("unknown matrix type {other:?}"),
        }
    }
}

/// A dense row-major `n × n` `f32` distance matrix over the (deduplicated) rows.
#[derive(Debug, Clone)]
pub struct DistMatrix {
    pub n: usize,
    pub data: Vec<f32>,
}

impl DistMatrix {
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f32 {
        self.data[i * self.n + j]
    }
    #[inline]
    fn set(&mut self, i: usize, j: usize, v: f32) {
        self.data[i * self.n + j] = v;
    }
}

/// Compute the distance matrix for the deduplicated profiles in `p`.
///
/// `block_penalty` is only consulted for `Blockwise`.
pub fn compute(
    p: &Parsed,
    kind: MatrixKind,
    handle_missing: HandleMissing,
    block_penalty: f64,
) -> DistMatrix {
    let n = p.n_rows;
    let l = p.n_cols;
    let codes = &p.codes;

    // Per-row presence (allele code > 0). For `complete_delete`+symmetric the
    // reference instead uses a per-column "present in all rows" mask.
    let mut data = vec![0f32; n * n];

    match kind {
        MatrixKind::Symmetric => {
            symmetric_into(&mut data, codes, n, l, handle_missing);
        }
        MatrixKind::Asymmetric => {
            asymmetric_into(
                &mut data,
                codes,
                n,
                l,
                handle_missing == HandleMissing::AbsoluteDistance,
            );
        }
        MatrixKind::AsymmetricWgMlst => {
            asymmetric_wgmlst_into(
                &mut data,
                codes,
                n,
                l,
                handle_missing == HandleMissing::AbsoluteDistance,
            );
        }
        MatrixKind::Blockwise => {
            blockwise_into(&mut data, codes, n, l, block_penalty);
        }
    }

    let mut m = DistMatrix { n, data };
    if kind == MatrixKind::Symmetric {
        // get_distance symmetrises symmetric matrices with an elementwise max.
        // Our symmetric builder is already symmetric, but we apply it defensively
        // so any floating asymmetry collapses exactly as the reference would.
        symmetrise_max(&mut m);
    }
    m
}

#[inline]
fn present(code: u32) -> bool {
    code > 0
}

/// `symmetric` for all four handlers. Only the lower triangle is meaningful;
/// we fill the full matrix symmetrically (diagonal stays 0).
fn symmetric_into(data: &mut [f32], codes: &[u32], n: usize, l: usize, hm: HandleMissing) {
    // complete_delete uses a per-column "present in every row" mask.
    let col_all_present: Option<Vec<bool>> = if hm == HandleMissing::CompleteDelete {
        Some(
            (0..l)
                .map(|c| (0..n).all(|r| present(codes[r * l + c])))
                .collect(),
        )
    } else {
        None
    };

    // Compute row-major; parallelise over the query row `i`, fill [i][j] for j<i.
    let rows: Vec<(usize, Vec<f32>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let pi = &codes[i * l..(i + 1) * l];
            let mut out = vec![0f32; i];
            for j in 0..i {
                let pj = &codes[j * l..(j + 1) * l];
                let v = match hm {
                    HandleMissing::PairDelete => {
                        let mut comparable = 0u64;
                        let mut diff = 0u64;
                        for k in 0..l {
                            if present(pi[k]) && present(pj[k]) {
                                comparable += 1;
                                if pi[k] != pj[k] {
                                    diff += 1;
                                }
                            }
                        }
                        ((diff as f64 + 0.01) * l as f64 / (comparable as f64 + 0.01)) as f32
                    }
                    HandleMissing::AsAllele => {
                        let mut diff = 0u64;
                        for k in 0..l {
                            if pi[k] != pj[k] {
                                diff += 1;
                            }
                        }
                        diff as f32
                    }
                    HandleMissing::AbsoluteDistance => {
                        let mut diff = 0u64;
                        for k in 0..l {
                            if present(pi[k]) && present(pj[k]) && pi[k] != pj[k] {
                                diff += 1;
                            }
                        }
                        diff as f32
                    }
                    HandleMissing::CompleteDelete => {
                        let mask = col_all_present.as_ref().unwrap();
                        let mut diff = 0u64;
                        for k in 0..l {
                            if mask[k] && pi[k] != pj[k] {
                                diff += 1;
                            }
                        }
                        diff as f32
                    }
                };
                out[j] = v;
            }
            (i, out)
        })
        .collect();

    for (i, out) in rows {
        for (j, v) in out.into_iter().enumerate() {
            data[i * n + j] = v;
            data[j * n + i] = v;
        }
    }
}

/// `asymmetric`: `d[a][b]` restricted to loci present in the reference `b`.
fn asymmetric_into(data: &mut [f32], codes: &[u32], n: usize, l: usize, absolute: bool) {
    // present count per reference column b
    let present_count: Vec<u64> = (0..n)
        .map(|b| (0..l).filter(|&k| present(codes[b * l + k])).count() as u64)
        .collect();

    let cols: Vec<(usize, Vec<f32>)> = (0..n)
        .into_par_iter()
        .map(|b| {
            let pb = &codes[b * l..(b + 1) * l];
            let pc = present_count[b];
            let mut col = vec![0f32; n];
            for a in 0..n {
                if a == b {
                    continue;
                }
                let pa = &codes[a * l..(a + 1) * l];
                let mut diff = 0u64;
                for k in 0..l {
                    if present(pb[k]) && pa[k] != pb[k] {
                        diff += 1;
                    }
                }
                col[a] = if absolute {
                    diff as f32
                } else {
                    (diff as f64 * l as f64 / pc as f64) as f32
                };
            }
            (b, col)
        })
        .collect();

    for (b, col) in cols {
        for (a, v) in col.into_iter().enumerate() {
            data[a * n + b] = v;
        }
    }
}

/// `asymmetric_wgMLST`: adds a fractional penalty for loci present in the
/// reference but absent in the compared strain.
fn asymmetric_wgmlst_into(data: &mut [f32], codes: &[u32], n: usize, l: usize, absolute: bool) {
    // pp[k] = P*(P-1)/(N*(N-1)) with P = #rows present at locus k.
    let nn = n as f64;
    let pp: Vec<f64> = (0..l)
        .map(|k| {
            let p_cnt = (0..n).filter(|&r| present(codes[r * l + k])).count() as f64;
            if nn > 1.0 {
                p_cnt * (p_cnt - 1.0) / (nn * (nn - 1.0))
            } else {
                0.0
            }
        })
        .collect();
    let present_count: Vec<u64> = (0..n)
        .map(|b| (0..l).filter(|&k| present(codes[b * l + k])).count() as u64)
        .collect();

    let cols: Vec<(usize, Vec<f32>)> = (0..n)
        .into_par_iter()
        .map(|b| {
            let pb = &codes[b * l..(b + 1) * l];
            let pc = present_count[b];
            let mut col = vec![0f32; n];
            for a in 0..n {
                let pa = &codes[a * l..(a + 1) * l];
                if absolute {
                    // reference's absolute branch ignores wgMLST extra term
                    let mut diff = 0u64;
                    for k in 0..l {
                        if present(pb[k]) && pa[k] != pb[k] {
                            diff += 1;
                        }
                    }
                    col[a] = diff as f32;
                } else {
                    let mut acc = 0f64;
                    for k in 0..l {
                        let (ka, kb) = (present(pa[k]), present(pb[k]));
                        if ka && kb && pa[k] != pb[k] {
                            acc += 1.0;
                        }
                        if !ka && kb {
                            acc += pp[k];
                        }
                    }
                    col[a] = (acc * l as f64 / pc as f64) as f32;
                }
            }
            (b, col)
        })
        .collect();

    for (b, col) in cols {
        for (a, v) in col.into_iter().enumerate() {
            data[a * n + b] = v;
        }
    }
}

/// `blockwise`: consecutive differing loci after the first pay only `penalty`.
fn blockwise_into(data: &mut [f32], codes: &[u32], n: usize, l: usize, penalty: f64) {
    let cols: Vec<(usize, Vec<f32>)> = (0..n)
        .into_par_iter()
        .map(|b| {
            let pb = &codes[b * l..(b + 1) * l];
            let mut col = vec![0f32; n];
            for a in 0..n {
                if a == b {
                    continue;
                }
                let pa = &codes[a * l..(a + 1) * l];
                // diff sequence with sentinel 0 on both ends: prev starts at 0.
                let mut d1 = 0u64; // block starts (transition into a nonzero diff)
                let mut total_nonzero = 0u64;
                let mut prev: i64 = 0;
                for k in 0..l {
                    let cur: i64 = pa[k] as i64 - pb[k] as i64;
                    if cur != 0 {
                        total_nonzero += 1;
                        if cur != prev {
                            d1 += 1;
                        }
                    }
                    prev = cur;
                }
                // trailing sentinel 0 vs last cur: (0 != cur) & (0 != 0) -> false, no d1.
                let d2 = total_nonzero - d1;
                col[a] = (d1 as f64 + d2 as f64 * penalty) as f32;
            }
            (b, col)
        })
        .collect();

    for (b, col) in cols {
        for (a, v) in col.into_iter().enumerate() {
            data[a * n + b] = v;
        }
    }
}

/// Elementwise `m[i][j] = max(m[i][j], m[j][i])` (reference symmetrisation).
fn symmetrise_max(m: &mut DistMatrix) {
    let n = m.n;
    for i in 0..n {
        for j in 0..i {
            let a = m.get(i, j);
            let b = m.get(j, i);
            let v = a.max(b);
            m.set(i, j, v);
            m.set(j, i, v);
        }
    }
}

/// Render the PHYLIP square matrix for the `distance` method, expanding the
/// `embeded` groups back to every original strain.
///
/// Mirrors `methods.distance`: rows/cols are all original names ordered by
/// (representative kept-index, name); values divided by `n_loci` unless the
/// matrix is `absolute_distance`/`blockwise`.
pub fn phylip(
    p: &Parsed,
    m: &DistMatrix,
    kind: MatrixKind,
    handle_missing: HandleMissing,
) -> String {
    // original name -> kept representative index
    let mut pairs: Vec<(String, usize)> = Vec::new();
    for (kept_idx, group) in p.embeded.iter().enumerate() {
        for name in group {
            pairs.push((name.clone(), kept_idx));
        }
    }
    // sort by (representative index, name)
    pairs.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let names: Vec<&str> = pairs.iter().map(|(n, _)| n.as_str()).collect();
    let idx: Vec<usize> = pairs.iter().map(|(_, i)| *i).collect();

    let normalise =
        !(handle_missing == HandleMissing::AbsoluteDistance || kind == MatrixKind::Blockwise);
    let l = p.n_cols as f32;

    let mut out = String::new();
    out.push_str(&format!("    {}\n", names.len()));
    for (row, &i2) in names.iter().zip(idx.iter()) {
        // name left-justified to width 10 then a space (Python "{0!s:10} ")
        out.push_str(&format!("{:<10} ", row));
        let cells: Vec<String> = idx
            .iter()
            .map(|&j2| {
                let mut v = m.get(i2, j2);
                if normalise {
                    v /= l;
                }
                format!("{:.6}", v as f64)
            })
            .collect();
        out.push_str(&cells.join(" "));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::parse_and_reduce;

    // Corrected wgMLST (the behaviour upstream *intended* but never reached).
    // Hand-computed on the tiny 4×5 case: output[a][b] = acc(a,b)/present(b).
    #[test]
    fn wgmlst_corrected_values() {
        let text = "#Strain\tL1\tL2\tL3\tL4\tL5\nA\t1\t1\t1\t1\t1\nB\t1\t2\t0\t1\t1\nC\t1\t1\t1\t2\t0\nD\t2\t0\t0\t1\t1\n";
        let p = parse_and_reduce(text, HandleMissing::PairDelete);
        let m = compute(
            &p,
            MatrixKind::AsymmetricWgMlst,
            HandleMissing::PairDelete,
            0.01,
        );
        let out = phylip(
            &p,
            &m,
            MatrixKind::AsymmetricWgMlst,
            HandleMissing::PairDelete,
        );
        // corrected wgMLST differs from the (buggy) upstream no-op; row A is all-present
        // so its distances equal plain asymmetric and are a stable anchor.
        assert!(out.contains("A         "), "expected strain A row present");
        // acc(C,D)=2.5, present(D)=3 -> 0.833333 (vs upstream-asymmetric 1.0)
        assert!(
            out.contains("0.833333"),
            "corrected wgMLST C->D cell missing:\n{out}"
        );
    }

    #[test]
    fn symmetric_pair_delete_is_symmetric() {
        let text = "#Strain\tL1\tL2\tL3\nA\t1\t1\t1\nB\t1\t2\t1\nC\t2\t2\t2\n";
        let p = parse_and_reduce(text, HandleMissing::PairDelete);
        let m = compute(&p, MatrixKind::Symmetric, HandleMissing::PairDelete, 0.01);
        for i in 0..m.n {
            for j in 0..m.n {
                assert_eq!(
                    m.get(i, j),
                    m.get(j, i),
                    "symmetric matrix must be symmetric"
                );
            }
        }
    }
}
