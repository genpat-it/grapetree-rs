// Reads a dense cost matrix (tab-separated rows) from a file arg; prints the
// optimum-branching edges as "from\tto", one per line. Diagonal ignored.
use std::io::Read;
fn main() {
    let path = std::env::args().nth(1).unwrap();
    let mut s = String::new();
    std::fs::File::open(&path)
        .unwrap()
        .read_to_string(&mut s)
        .unwrap();
    let rows: Vec<Vec<f64>> = s
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            l.split('\t')
                .map(|x| x.trim().parse::<f64>().unwrap())
                .collect()
        })
        .collect();
    let n = rows.len();
    let arb = grapetree::edmonds::optimum_branching(n, |i, j| Some(rows[i][j]));
    let mut e: Vec<(usize, usize)> = arb;
    e.sort();
    for (u, v) in e {
        println!("{}\t{}", u, v);
    }
}
