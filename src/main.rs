use anyhow::{Context, Result};
use clap::Parser;
use grapetree::cli::Args;
use grapetree::distance::{self, MatrixKind};
use grapetree::heuristic::{self, Heuristic};
use grapetree::mst;
use grapetree::parse;
use grapetree::recraft;
use grapetree::tree::Tree;
use std::io::Read;

/// Port of `estimate_Consumption` (Linux coefficients). `method`/`matrix` are
/// the already-resolved values. Returns `(seconds, bytes)`.
fn estimate_consumption(
    method: &str,
    matrix: &str,
    n_proc: usize,
    n_loci: usize,
    n_profile: usize,
) -> (f64, f64) {
    let (p, l, np) = (n_profile as f64, n_loci as f64, n_proc as f64);
    let (time, memory) = match method {
        "MSTree" | "RapidNJ" | "ninja" => {
            if matrix == "asymmetric" {
                (
                    2.431284e-6 * p * p + 2.701426667e-9 * l * p * p / np + 33.753,
                    103.77 * p * p + 516_625_000.0,
                )
            } else {
                (
                    2.272428e-6 * p * p + 32.625 + 2.52492e-9 * l * p * p / np,
                    66.297 * p * p + 429_570_000.0,
                )
            }
        }
        "NJ" => (
            1.1042e-8 * p * p * p,
            (0.058292 * p * p * p).max(1.39e6 * p - 9.86e8),
        ),
        _ => (5.0, 50.0 * 1024.0 * 1024.0),
    };
    (time.max(5.0), memory.max(50.0 * 1024.0 * 1024.0))
}

/// Available RAM in bytes (Linux `/proc/meminfo` MemAvailable), if readable.
fn available_memory_bytes() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in s.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kb: u64 = rest.trim().trim_end_matches("kB").trim().parse().ok()?;
            return Some(kb * 1024);
        }
    }
    None
}

fn read_input(path: &str) -> Result<String> {
    // Reference accepts a filename or the literal file contents. If the path
    // exists on disk, read it (decompressing `.gz`); otherwise treat as inline.
    if std::path::Path::new(path).is_file() {
        if path.to_lowercase().ends_with(".gz") {
            let f = std::fs::File::open(path).with_context(|| format!("open {path}"))?;
            let mut gz = flate2::read::MultiGzDecoder::new(f);
            let mut s = String::new();
            gz.read_to_string(&mut s)?;
            Ok(s)
        } else {
            std::fs::read_to_string(path).with_context(|| format!("read {path}"))
        }
    } else {
        Ok(path.to_string())
    }
}

fn main() -> Result<()> {
    let params = Args::parse().resolve();

    // Respect -n/--n_proc for the rayon thread pool (performance only).
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(params.n_proc.max(1))
        .build_global();

    let text = read_input(&params.profile)?;
    let parsed = parse::parse_and_reduce(&text, params.handle_missing);

    let kind = MatrixKind::resolve(&params.matrix_type, params.wg_mlst);

    if params.check_env {
        // -c/--check: report the estimated time/memory requirement, like the
        // reference's estimate_Consumption + backend() JSON.
        let (time, memory) = estimate_consumption(
            &params.method,
            &params.matrix_type,
            params.n_proc.max(1),
            parsed.n_cols,
            parsed.n_rows,
        );
        // `affordable` is machine-dependent (free RAM); report against this host.
        let affordable = available_memory_bytes()
            .map(|m| (m as f64) >= memory)
            .unwrap_or(true);
        println!(
            "{{\"time\": {}, \"memory\": {}, \"affordable\": {}}}",
            time, memory, affordable
        );
        return Ok(());
    }

    match params.method.as_str() {
        "distance" => {
            let m = distance::compute(&parsed, kind, params.handle_missing, params.block_penalty);
            let out = distance::phylip(&parsed, &m, kind, params.handle_missing);
            print!("{out}");
        }
        "MSTree" => {
            let dm = distance::compute(&parsed, kind, params.handle_missing, params.block_penalty);
            let n_str = parsed.n_str();
            let heur = Heuristic::parse(&params.heuristic);
            let weight = heuristic::weights(&dm, &n_str, heur);
            let edges: Vec<(usize, usize)> = match kind {
                MatrixKind::Symmetric | MatrixKind::Blockwise => mst::symmetric_mst(&dm, &weight),
                MatrixKind::Asymmetric => {
                    let net = mst::asymmetric_network(&dm, &weight);
                    let net = if params.branch_recraft {
                        recraft::branch_recraft(net, &dm, &weight, parsed.n_cols as f64)
                    } else {
                        net
                    };
                    net.into_iter().map(|(s, t, _)| (s, t)).collect()
                }
                MatrixKind::AsymmetricWgMlst => unreachable!("wgMLST is not produced by resolve()"),
            };
            let linked = mst::symmetric_link(&parsed, &edges, params.handle_missing);
            let mut tree = Tree::network2tree(&linked, &parsed.names);
            tree.post_process(&parsed);
            println!("{}", tree.to_newick());
        }
        "NJ" | "RapidNJ" | "ninja" => {
            // The reference builds all NJ variants on the symmetric matrix via
            // FastME/RapidNJ/Ninja; we use a native canonical NJ instead.
            let dm = distance::compute(
                &parsed,
                MatrixKind::Symmetric,
                params.handle_missing,
                params.block_penalty,
            );
            let nwk = grapetree::nj::neighbor_joining(&parsed, &dm, &parsed.names);
            println!("{nwk}");
        }
        other => {
            eprintln!(
                "[grapetree-rs] method={other} matrix={:?} heuristic={} recraft={} missing={} :: {} unique rows × {} loci (this method lands in a later iteration)",
                kind,
                params.heuristic,
                params.branch_recraft,
                params.handle_missing.as_str(),
                parsed.n_rows,
                parsed.n_cols,
            );
        }
    }
    Ok(())
}
