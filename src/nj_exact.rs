//! Bit-identical NJ family (`NJ` / `RapidNJ` / `ninja`).
//!
//! The reference builds these trees with external tools and manipulates them
//! with ete3. We do the same — "what is external stays external": the tree is
//! built by the bundled FastME / RapidNJ / Ninja binary, and the ete3 tail
//! (midpoint-root, unroot, rename, collapse, `embeded` expansion, NEWICK
//! `format=1`) is done by the bundled `shim/nj_postprocess.py`, verbatim from
//! GrapeTree. Result: byte-identical output, with the same external runtime
//! dependencies GrapeTree itself has (the binaries + Python/ete3).
//!
//! grapetree-rs only contributes the fast Rust distance matrix here; everything
//! that could diverge is delegated to the reference's own toolchain.

use crate::distance::DistMatrix;
use crate::parse::Parsed;
use std::path::Path;
use std::process::Command;

/// Distance file (PHYLIP square) with integer indices as taxon names — exactly
/// the reference NJ writer (`'{0!s:10} {1}'`).
fn write_phylip(dist: &DistMatrix, scale: f64) -> String {
    let n = dist.n;
    let mut s = format!("    {}\n", n);
    for i in 0..n {
        s.push_str(&format!("{:<10} ", i));
        let row: Vec<String> = (0..n)
            .map(|j| format!("{:.6}", dist.get(i, j) as f64 / scale))
            .collect();
        s.push_str(&row.join(" "));
        s.push('\n');
    }
    s
}

/// Run FastME (`-m N`) on the distance matrix and return the raw NEWICK.
pub fn run_fastme(dist: &DistMatrix, fastme: &Path) -> Option<String> {
    let base = std::env::temp_dir().join(format!("gt_nj_{}.list", std::process::id()));
    std::fs::write(&base, write_phylip(dist, 1.0)).ok()?;
    let out_name = format!("{}_fastme_tree.nwk", base.file_name()?.to_str()?);
    let outp = base.with_file_name(&out_name);
    Command::new(fastme)
        .args(["-i", base.to_str()?, "-m", "N"])
        .output()
        .ok()?;
    let nwk = std::fs::read_to_string(&outp).ok();
    let stat = base.with_file_name(format!("{}_fastme_stat.txt", base.file_name()?.to_str()?));
    for f in [&base, &outp, &stat] {
        let _ = std::fs::remove_file(f);
    }
    nwk
}

/// Run RapidNJ (`-i pd -n -x out`) on the distance matrix; return raw NEWICK.
pub fn run_rapidnj(dist: &DistMatrix, rapidnj: &Path) -> Option<String> {
    let base = std::env::temp_dir().join(format!("gt_rnj_{}.list", std::process::id()));
    std::fs::write(&base, write_phylip(dist, 1.0)).ok()?;
    let outp = base.with_file_name(format!("{}_rapidnj.nwk", base.file_name()?.to_str()?));
    Command::new(rapidnj)
        .args(["-n", "-x", outp.to_str()?, "-i", "pd", base.to_str()?])
        .output()
        .ok()?;
    let nwk = std::fs::read_to_string(&outp).ok();
    let _ = std::fs::remove_file(&base);
    let _ = std::fs::remove_file(&outp);
    nwk
}

/// Run Ninja (`java -jar Ninja.jar --in_type d`). The reference divides the
/// matrix by `n_loci` first and multiplies branch lengths back afterwards
/// (the shim `scale` does the multiply). Returns raw NEWICK.
pub fn run_ninja(dist: &DistMatrix, jar: &Path, java: &str, n_loci: f64) -> Option<String> {
    let base = std::env::temp_dir().join(format!("gt_ninja_{}.list", std::process::id()));
    std::fs::write(&base, write_phylip(dist, n_loci)).ok()?;
    let out = Command::new(java)
        .args(["-jar", jar.to_str()?, "--in_type", "d", base.to_str()?])
        .output();
    let _ = std::fs::remove_file(&base);
    let out = out.ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// Full bit-identical NJ post-processing: hand the raw tree + names + embeded to
/// the bundled ete3 shim and return the final NEWICK. `scale` multiplies branch
/// lengths first (n_loci for `ninja`, else 1). Returns `None` on any failure.
pub fn neighbor_joining_exact(
    p: &Parsed,
    tree_nwk: &str,
    shim: &Path,
    python: &str,
    scale: f64,
) -> Option<String> {
    let names_json: String = {
        let items: Vec<String> = p.names.iter().map(|n| json_str(n)).collect();
        format!("[{}]", items.join(","))
    };
    let emb_json: String = {
        let items: Vec<String> = p
            .embeded
            .iter()
            .map(|g| {
                let members: Vec<String> = g.iter().map(|m| json_str(m)).collect();
                format!("{}:[{}]", json_str(&g[0]), members.join(","))
            })
            .collect();
        format!("{{{}}}", items.join(","))
    };
    let base = std::env::temp_dir().join(format!("gt_njpp_{}", std::process::id()));
    let nwk_p = base.with_extension("nwk");
    let names_p = base.with_extension("names.json");
    let emb_p = base.with_extension("emb.json");
    std::fs::write(&nwk_p, tree_nwk).ok()?;
    std::fs::write(&names_p, names_json).ok()?;
    std::fs::write(&emb_p, emb_json).ok()?;
    let out = Command::new(python)
        .arg(shim)
        .arg(&nwk_p)
        .arg(&names_p)
        .arg(&emb_p)
        .arg(format!("{scale}"))
        .output();
    for f in [&nwk_p, &names_p, &emb_p] {
        let _ = std::fs::remove_file(f);
    }
    let out = out.ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// Minimal JSON string escaping (names are already sanitised by `sanitize_name`,
/// but escape defensively).
fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            _ => out.push(c),
        }
    }
    out.push('"');
    out
}
