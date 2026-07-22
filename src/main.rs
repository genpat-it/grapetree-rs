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

/// Locate a bundled file under `binaries/` or `shim/` (next to the executable,
/// up to a few levels up, or in the source tree).
fn resolve_bundled(subdir: &str, name: &str) -> Option<std::path::PathBuf> {
    use std::path::PathBuf;
    if let Ok(exe) = std::env::current_exe() {
        let mut dir = exe.parent().map(|p| p.to_path_buf());
        for _ in 0..4 {
            if let Some(d) = &dir {
                let pb = d.join(subdir).join(name);
                if pb.is_file() {
                    return Some(pb);
                }
                dir = d.parent().map(|p| p.to_path_buf());
            }
        }
    }
    let pb = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(subdir)
        .join(name);
    pb.is_file().then_some(pb)
}

fn resolve_edmonds() -> Option<std::path::PathBuf> {
    if let Ok(p) = std::env::var("GT_EDMONDS") {
        let pb = std::path::PathBuf::from(p);
        if pb.is_file() {
            return Some(pb);
        }
    }
    resolve_bundled("binaries", "edmonds-linux")
}

/// Full bit-identical NJ pipeline via the bundled binaries + ete3 shim.
fn nj_exact_pipeline(
    parsed: &parse::Parsed,
    dm: &distance::DistMatrix,
    method: &str,
) -> Option<String> {
    use grapetree::nj_exact;
    let shim = resolve_bundled("shim", "nj_postprocess.py")?;
    let python = std::env::var("GT_PYTHON").unwrap_or_else(|_| "python3".to_string());
    let n_loci = parsed.n_cols as f64;
    let (raw, scale) = match method {
        "NJ" => (
            nj_exact::run_fastme(dm, &resolve_bundled("binaries", "fastme-2.1.5-linux64")?)?,
            1.0,
        ),
        "RapidNJ" => (
            nj_exact::run_rapidnj(dm, &resolve_bundled("binaries", "rapidnj")?)?,
            1.0,
        ),
        "ninja" => {
            let java = std::env::var("GT_JAVA").unwrap_or_else(|_| "java".to_string());
            (
                nj_exact::run_ninja(
                    dm,
                    &resolve_bundled("binaries", "Ninja.jar")?,
                    &java,
                    n_loci,
                )?,
                n_loci,
            )
        }
        _ => return None,
    };
    nj_exact::neighbor_joining_exact(parsed, &raw, &shim, &python, scale)
}

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
            let timing = std::env::var("GT_TIMING").is_ok();
            let t0 = std::time::Instant::now();
            let dm = distance::compute(&parsed, kind, params.handle_missing, params.block_penalty);
            if timing {
                eprintln!("[timing] distance: {:.2}s", t0.elapsed().as_secs_f64());
            }
            let t1 = std::time::Instant::now();
            let n_str = parsed.n_str();
            let heur = Heuristic::parse(&params.heuristic);
            let weight = heuristic::weights(&dm, &n_str, heur);
            if timing {
                eprintln!("[timing] weights: {:.2}s", t1.elapsed().as_secs_f64());
            }
            let t2 = std::time::Instant::now();
            let edges: Vec<(usize, usize)> = match kind {
                MatrixKind::Symmetric | MatrixKind::Blockwise => mst::symmetric_mst(&dm, &weight),
                MatrixKind::Asymmetric => {
                    // Bit-identical by default: delegate the arborescence to the
                    // bundled `edmonds` binary (as upstream does). `--native`
                    // uses the pure-Rust Edmonds (equivalent, not bit-identical).
                    let net = if params.native {
                        mst::asymmetric_network(&dm, &weight)
                    } else {
                        match resolve_edmonds() {
                            Some(path) => mst::asymmetric_network_exact(&dm, &weight, &path)
                                .unwrap_or_else(|| {
                                    eprintln!("[grapetree-rs] edmonds binary failed; falling back to native (equivalent, not bit-identical)");
                                    mst::asymmetric_network(&dm, &weight)
                                }),
                            None => {
                                eprintln!("[grapetree-rs] edmonds binary not found; using native (equivalent, not bit-identical). Set GT_EDMONDS or pass --native to silence.");
                                mst::asymmetric_network(&dm, &weight)
                            }
                        }
                    };
                    if timing {
                        eprintln!("[timing] arborescence: {:.2}s", t2.elapsed().as_secs_f64());
                    }
                    let tr = std::time::Instant::now();
                    let net = if params.branch_recraft {
                        recraft::branch_recraft(net, &dm, &weight, parsed.n_cols as f64)
                    } else {
                        net
                    };
                    if timing {
                        eprintln!("[timing] recraft: {:.2}s", tr.elapsed().as_secs_f64());
                    }
                    net.into_iter().map(|(s, t, _)| (s, t)).collect()
                }
                MatrixKind::AsymmetricWgMlst => unreachable!("wgMLST is not produced by resolve()"),
            };
            let t3 = std::time::Instant::now();
            let linked = mst::symmetric_link(&parsed, &edges, params.handle_missing);
            let mut tree = Tree::network2tree(&linked, &parsed.names);
            tree.post_process(&parsed);
            let nwk = tree.to_newick();
            if timing {
                eprintln!("[timing] link+tree: {:.2}s", t3.elapsed().as_secs_f64());
            }
            println!("{}", nwk);
        }
        "NJ" | "RapidNJ" | "ninja" => {
            let dm = distance::compute(
                &parsed,
                MatrixKind::Symmetric,
                params.handle_missing,
                params.block_penalty,
            );
            // Bit-identical by default: delegate the tree to the bundled
            // FastME/RapidNJ/Ninja binary + the ete3 post-processing shim (as
            // upstream does). `--native` uses the pure-Rust canonical NJ.
            let exact = if params.native {
                None
            } else {
                nj_exact_pipeline(&parsed, &dm, &params.method)
            };
            match exact {
                Some(nwk) => println!("{}", nwk.trim_end()),
                None => {
                    if !params.native {
                        eprintln!("[grapetree-rs] NJ external toolchain unavailable; using native canonical NJ (topologically equivalent, not bit-identical). Pass --native to silence.");
                    }
                    let nwk = grapetree::nj::neighbor_joining(&parsed, &dm, &parsed.names);
                    println!("{nwk}");
                }
            }
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
