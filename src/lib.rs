//! grapetree-rs — a Rust port of GrapeTree's computational backend.
//!
//! See `DECISIONS.md` for the design rationale and the parameter→behaviour map.
//!
//! Row-major matrix code intentionally uses index loops for clarity/parity with
//! the reference implementation.
#![allow(clippy::needless_range_loop)]

pub mod cli;
pub mod distance;
pub mod edmonds;
pub mod edmonds_tofigh;
pub mod heuristic;
pub mod mst;
pub mod nj;
pub mod nj_exact;
pub mod parse;
pub mod recraft;
pub mod tree;

/// How missing alleles (`0`, `N`, `-`) are handled in pairwise comparisons.
///
/// Mirrors the reference's `handle_missing` strings and the `-y/--missing`
/// integer codes (0..3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HandleMissing {
    /// `-y 0` (default): ignore missing in each pairwise comparison.
    PairDelete,
    /// `-y 1`: keep only loci genotyped in every strain.
    CompleteDelete,
    /// `-y 2`: treat the missing marker as a real allele.
    AsAllele,
    /// `-y 3`: absolute number of allelic differences.
    AbsoluteDistance,
}

impl HandleMissing {
    pub fn from_code(code: i64) -> Self {
        match code {
            0 => HandleMissing::PairDelete,
            1 => HandleMissing::CompleteDelete,
            2 => HandleMissing::AsAllele,
            3 => HandleMissing::AbsoluteDistance,
            other => panic!("invalid --missing code {other}, expected 0..3"),
        }
    }
    pub fn as_str(self) -> &'static str {
        match self {
            HandleMissing::PairDelete => "pair_delete",
            HandleMissing::CompleteDelete => "complete_delete",
            HandleMissing::AsAllele => "as_allele",
            HandleMissing::AbsoluteDistance => "absolute_distance",
        }
    }
}
