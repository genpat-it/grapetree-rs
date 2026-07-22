//! Command-line interface mirroring the reference `add_args()`.

use crate::HandleMissing;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "grapetree",
    about = "Rust port of GrapeTree: NEWICK tree / distance matrix from allelic profiles",
    long_about = None,
)]
pub struct Args {
    /// Input file (MLST/SNP profile TSV, or aligned FASTA). `.gz` supported.
    #[arg(short = 'p', long = "profile")]
    pub profile: String,

    /// MSTreeV2 [default], MSTree, NJ, RapidNJ, ninja, distance.
    #[arg(short = 'm', long = "method", default_value = "MSTreeV2")]
    pub method: String,

    /// symmetric [default for MSTree/NJ], asymmetric, blockwise.
    #[arg(short = 'x', long = "matrix", default_value = "symmetric")]
    pub matrix_type: String,

    /// Trigger local branch recrafting (forced on for MSTreeV2).
    #[arg(short = 'r', long = "recraft", default_value_t = false)]
    pub branch_recraft: bool,

    /// Missing-data handler: 0 pair_delete [default], 1 complete_delete,
    /// 2 as_allele, 3 absolute_distance. (symmetric distance matrix only)
    #[arg(short = 'y', long = "missing", default_value_t = 0)]
    pub handler: i64,

    /// [experimental] better support for wgMLST schemes.
    #[arg(short = 'w', long = "wgMLST", default_value_t = false)]
    pub wg_mlst: bool,

    /// Tiebreak heuristic: eBurst [default MSTree] or harmonic [default MSTreeV2].
    #[arg(short = 't', long = "heuristic", default_value = "eBurst")]
    pub heuristic: String,

    /// Number of parallel worker threads (performance only; results unchanged).
    #[arg(short = 'n', long = "n_proc", default_value_t = 5)]
    pub n_proc: usize,

    /// Only report estimated time/memory requirements.
    #[arg(short = 'c', long = "check", default_value_t = false)]
    pub check_env: bool,

    /// Penalty for a different locus led by another difference (blockwise only).
    #[arg(short = 'b', long = "block_penalty", default_value_t = 0.01)]
    pub block_penalty: f64,
}

/// Fully resolved run parameters after applying the MSTreeV2/blockwise aliases.
#[derive(Debug, Clone)]
pub struct Params {
    pub profile: String,
    pub method: String,
    pub matrix_type: String,
    pub heuristic: String,
    pub branch_recraft: bool,
    pub handle_missing: HandleMissing,
    pub wg_mlst: bool,
    pub n_proc: usize,
    pub check_env: bool,
    pub block_penalty: f64,
}

impl Args {
    /// Resolve aliases exactly like `add_args()` + `backend()`:
    /// - MSTreeV2 -> MSTree + asymmetric + harmonic + recraft
    /// - blockwise -> MSTree, recraft off, alleles all real (handle via block_penalty)
    pub fn resolve(self) -> Params {
        let mut method = self.method.clone();
        let mut matrix_type = self.matrix_type.clone();
        let mut heuristic = self.heuristic.clone();
        let mut branch_recraft = self.branch_recraft;
        let mut handle_missing = HandleMissing::from_code(self.handler);

        if matrix_type == "blockwise" {
            if method == "MSTreeV2" {
                method = "MSTree".to_string();
            }
            branch_recraft = false;
            // reference stashes block_penalty into handle_missing here; we keep
            // handle_missing typed and pass block_penalty separately.
        }
        if method == "MSTreeV2" {
            method = "MSTree".to_string();
            matrix_type = "asymmetric".to_string();
            heuristic = "harmonic".to_string();
            branch_recraft = true;
            // MSTreeV2 ignores the missing handler in the same way the reference does
            let _ = &mut handle_missing;
        }

        Params {
            profile: self.profile,
            method,
            matrix_type,
            heuristic,
            branch_recraft,
            handle_missing,
            wg_mlst: self.wg_mlst,
            n_proc: self.n_proc,
            check_env: self.check_env,
            block_penalty: self.block_penalty,
        }
    }
}
