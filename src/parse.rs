//! Input parsing and `nonredundant` allele encoding.
//!
//! Faithful port of the reference `backend()` reader and `nonredundant()`
//! (module/MSTrees.py). We reproduce the exact column-selection, uppercasing,
//! name-sanitisation, integer-encoding (sorted-string order), lexsort ordering
//! and duplicate collapsing so that downstream node identities match.

use crate::HandleMissing;
use std::collections::BTreeMap;

/// A parsed, deduplicated, integer-encoded profile matrix.
#[derive(Debug, Clone)]
pub struct Parsed {
    /// Representative name of each kept (unique) row, in lexsorted order.
    pub names: Vec<String>,
    /// Row-major `n_rows * n_cols` allele codes (0 = missing).
    pub codes: Vec<u32>,
    pub n_rows: usize,
    pub n_cols: usize,
    /// Members collapsed into each kept row (`embeded` groups). `group[i]` are
    /// the original names identical to `names[i]`, including `names[i]` itself.
    pub embeded: Vec<Vec<String>>,
}

impl Parsed {
    #[inline]
    pub fn row(&self, i: usize) -> &[u32] {
        &self.codes[i * self.n_cols..(i + 1) * self.n_cols]
    }
    /// `n_str`: number of original strains collapsed into each kept row.
    pub fn n_str(&self) -> Vec<usize> {
        self.embeded.iter().map(|g| g.len()).collect()
    }
}

/// Sanitise a strain name: `[()  ,"';]` -> `_` (regex `[\(\)\ \,\"\';]`).
fn sanitize_name(n: &str) -> String {
    n.chars()
        .map(|c| match c {
            '(' | ')' | ' ' | ',' | '"' | '\'' | ';' => '_',
            other => other,
        })
        .collect()
}

/// A column header should be dropped if it starts with `#` or equals `ST`/`ST_id`.
fn is_dropped_col(col: &str) -> bool {
    let low = col.to_lowercase();
    col.starts_with('#') || low == "st_id" || low == "st"
}

/// Raw parse result before encoding.
struct Raw {
    names: Vec<String>,
    /// Row-major string matrix, already uppercased.
    values: Vec<String>,
    n_rows: usize,
    n_cols: usize,
}

/// Parse profile-TSV or FASTA text into a raw uppercased string matrix.
///
/// `text` is the full file contents (already decompressed if it was `.gz`).
fn parse_raw(text: &str) -> Raw {
    let lines: Vec<&str> = text.split('\n').collect();

    // ----- header / format detection (mirrors reference loop) -----
    let mut allele_cols: Option<Vec<usize>> = None;
    let mut is_fasta = false;
    let mut data_start = 0usize;

    for (line_id, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            if !line.starts_with("##") {
                let header: Vec<&str> = line.trim_end_matches(['\r']).trim().split('\t').collect();
                allele_cols = Some(
                    header
                        .iter()
                        .enumerate()
                        .filter(|(id, col)| *id > 0 && !is_dropped_col(col))
                        .map(|(id, _)| id)
                        .collect(),
                );
            }
            continue;
        }
        if line.starts_with('>') {
            is_fasta = true;
            data_start = line_id;
        } else {
            // profile row
            if allele_cols.is_none() {
                let header: Vec<&str> = line.trim().split('\t').collect();
                allele_cols = Some(
                    header
                        .iter()
                        .enumerate()
                        .filter(|(id, col)| *id > 0 && !is_dropped_col(col))
                        .map(|(id, _)| id)
                        .collect(),
                );
                data_start = line_id + 1; // this line was the header
            } else {
                data_start = line_id;
            }
        }
        break;
    }

    let mut names: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<String>> = Vec::new();

    if is_fasta {
        for line in &lines[data_start..] {
            if let Some(rest) = line.strip_prefix('>') {
                let name = rest.split_whitespace().next().unwrap_or("").to_string();
                names.push(name);
                rows.push(Vec::new());
            } else if !rows.is_empty() {
                // Each character of the (concatenated) sequence is a locus.
                for ch in line.split_whitespace().flat_map(|s| s.chars()) {
                    rows.last_mut()
                        .unwrap()
                        .push(ch.to_ascii_uppercase().to_string());
                }
            }
        }
    } else {
        let cols = allele_cols.as_ref();
        for line in &lines[data_start..] {
            let part: Vec<&str> = line.trim_end_matches(['\r']).trim().split('\t').collect();
            if part.is_empty() || part[0].is_empty() {
                continue;
            }
            names.push(part[0].to_string());
            let row: Vec<String> = match cols {
                Some(c) => c
                    .iter()
                    .map(|&idx| part.get(idx).copied().unwrap_or("").to_uppercase())
                    .collect(),
                None => part[1..].iter().map(|s| s.to_uppercase()).collect(),
            };
            rows.push(row);
        }
    }

    let names: Vec<String> = names.iter().map(|n| sanitize_name(n)).collect();
    let n_rows = rows.len();
    let n_cols = if n_rows > 0 { rows[0].len() } else { 0 };
    let mut values = Vec::with_capacity(n_rows * n_cols);
    for r in &rows {
        assert_eq!(
            r.len(),
            n_cols,
            "ragged profile row: expected {n_cols} cols"
        );
        values.extend_from_slice(r);
    }
    Raw {
        names,
        values,
        n_rows,
        n_cols,
    }
}

/// True if a raw allele string denotes missing data (`0`, `N`, `-`).
#[inline]
fn is_missing(v: &str) -> bool {
    v == "0" || v == "N" || v == "-"
}

/// Full parse + `nonredundant`. Returns the deduplicated encoded matrix.
pub fn parse_and_reduce(text: &str, handle_missing: HandleMissing) -> Parsed {
    let raw = parse_raw(text);
    nonredundant(raw, handle_missing)
}

/// Port of `nonredundant`: integer-encode by sorted-string order, mark missing,
/// (optionally) drop columns for complete_delete, lexsort rows, drop all-missing
/// rows, and collapse identical rows into `embeded` groups.
fn nonredundant(raw: Raw, handle_missing: HandleMissing) -> Parsed {
    let Raw {
        names,
        values,
        n_rows,
        n_cols,
    } = raw;

    // --- per-column integer encoding by sorted (uppercased) string order ---
    // np.unique sorts distinct values; code = sorted-index + 1; missing -> 0.
    let mut encoded = vec![0u32; n_rows * n_cols];
    for c in 0..n_cols {
        // distinct values in this column (BTreeMap keeps them sorted)
        let mut order: BTreeMap<&str, u32> = BTreeMap::new();
        for r in 0..n_rows {
            let v = values[r * n_cols + c].as_str();
            order.entry(v).or_insert(0);
        }
        for (idx, code) in (1u32..).zip(order.values_mut()) {
            *code = idx;
        }
        for r in 0..n_rows {
            let v = values[r * n_cols + c].as_str();
            encoded[r * n_cols + c] = if is_missing(v) { 0 } else { order[v] };
        }
    }

    // --- complete_delete: keep only columns that contain a zero (reference verbatim) ---
    let kept_cols: Vec<usize> = if handle_missing == HandleMissing::CompleteDelete {
        (0..n_cols)
            .filter(|&c| (0..n_rows).any(|r| encoded[r * n_cols + c] == 0))
            .collect()
    } else {
        (0..n_cols).collect()
    };
    let m_cols = kept_cols.len();
    let mut mat = vec![0u32; n_rows * m_cols];
    for r in 0..n_rows {
        for (j, &c) in kept_cols.iter().enumerate() {
            mat[r * m_cols + j] = encoded[r * n_cols + c];
        }
    }

    // --- lexsort rows: last column is the primary key (np.lexsort semantics) ---
    let mut order: Vec<usize> = (0..n_rows).collect();
    order.sort_by(|&a, &b| {
        for j in (0..m_cols).rev() {
            let (x, y) = (mat[a * m_cols + j], mat[b * m_cols + j]);
            if x != y {
                return x.cmp(&y);
            }
        }
        std::cmp::Ordering::Equal
    });

    // --- drop all-missing rows, then collapse identical consecutive rows ---
    let mut names_out: Vec<String> = Vec::new();
    let mut codes_out: Vec<u32> = Vec::new();
    let mut embeded: Vec<Vec<String>> = Vec::new();

    let mut prev: Option<usize> = None;
    for &r in &order {
        let has_present = (0..m_cols).any(|j| mat[r * m_cols + j] > 0);
        if !has_present {
            continue;
        }
        let same_as_prev = match prev {
            Some(p) => (0..m_cols).all(|j| mat[r * m_cols + j] == mat[p * m_cols + j]),
            None => false,
        };
        if same_as_prev {
            embeded.last_mut().unwrap().push(names[r].clone());
        } else {
            names_out.push(names[r].clone());
            codes_out.extend_from_slice(&mat[r * m_cols..(r + 1) * m_cols]);
            embeded.push(vec![names[r].clone()]);
            prev = Some(r);
        }
    }

    Parsed {
        n_rows: names_out.len(),
        n_cols: m_cols,
        names: names_out,
        codes: codes_out,
        embeded,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_profile() {
        let text = "#Strain\tL1\tL2\tL3\nA\t1\t1\t1\nB\t1\t2\t1\nC\t1\t1\t1\n";
        let p = parse_and_reduce(text, HandleMissing::PairDelete);
        // A and C are identical -> collapsed; B distinct.
        assert_eq!(p.n_rows, 2);
        assert_eq!(p.n_cols, 3);
        // the AC group has size 2 somewhere
        let sizes: Vec<usize> = p.n_str();
        assert!(sizes.contains(&2) && sizes.contains(&1));
    }

    #[test]
    fn drops_st_and_hash_columns() {
        let text = "#Strain\tST\tL1\t#note\tL2\nA\t7\t1\tx\t1\nB\t7\t1\ty\t2\n";
        let p = parse_and_reduce(text, HandleMissing::PairDelete);
        assert_eq!(p.n_cols, 2, "ST and #note columns must be dropped");
    }

    #[test]
    fn missing_markers_encode_to_zero() {
        let text = "#Strain\tL1\tL2\nA\t0\tN\nB\t1\t-\n";
        let p = parse_and_reduce(text, HandleMissing::PairDelete);
        // A: [0,0] all-missing -> dropped; B kept. L1 distinct sorted = {"0","1"}
        // so "1" -> code 2 (matches numpy unique+1), then "-"/"N" -> 0.
        assert_eq!(p.n_rows, 1);
        assert_eq!(p.row(0), &[2, 0]);
    }
}
