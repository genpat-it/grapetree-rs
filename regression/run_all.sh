#!/usr/bin/env bash
# Comprehensive per-parameter regression: grapetree-rs vs the Python oracle.
#
# Covers every CLI method × matrix × missing handler × heuristic × recraft, on
# synthetic datasets with controlled properties (clonal structure, duplicates,
# missing-data gradients) plus the bundled examples. Scoring:
#   - distance / MSTree / MSTreeV2 : byte-identical NEWICK/PHYLIP
#   - MSTree asymmetric            : topological identity (RF=0) + equal length
#   - NJ family                    : topology only (RF=0); FastME uses balanced-ME
#                                    edge lengths, we use canonical Saitou-Nei NJ
#   - -n/--n_proc                  : output invariant to thread count
#
# Requires the `gt_oracle` conda env (Python 3.11 + ete3/numpy/numba/networkx).
#
# Configure via env vars:
#   GT_ORIG   path to a GrapeTree checkout (default: ./grapetree-orig)
#   GT_BIN    path to the grapetree-rs binary (default: target/release/grapetree)
#   GT_ENV    conda env with the Python oracle deps (default: gt_oracle)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
source ~/miniconda3/etc/profile.d/conda.sh && conda activate "${GT_ENV:-gt_oracle}"
ORIG="${GT_ORIG:-$ROOT/grapetree-orig}"
if [ ! -f "$ORIG/module/MSTrees.py" ]; then
  echo "GrapeTree reference not found at '$ORIG'. Set GT_ORIG to a checkout of"
  echo "https://github.com/achtman-lab/GrapeTree (needs module/MSTrees.py)." >&2
  exit 2
fi
MP=$ORIG/module/MSTrees.py
G=${GT_BIN:-$ROOT/target/release/grapetree}
CMP=$HERE/compare_trees.py
ORUN=$HERE/oracle_run.py
GEN=$HERE/gen_synth.py
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

byte=0; equiv=0; fail=0; skip=0; ninv=0; nfail=0
declare -a FAILED

# --- build synthetic datasets ---
DSDIR=$TMP/synth; mkdir -p $DSDIR
python $GEN miss   30 120 0.10 0.15 int 0 > $DSDIR/miss_int.profile
python $GEN missN  30 120 0.10 0.15 int N > $DSDIR/miss_N.profile
python $GEN snp    25 200 0.05 0.10 nuc - > $DSDIR/snp.profile
python $GEN clean  40 150 0.00 0.20 int 0 > $DSDIR/clean.profile
python $GEN big    80 300 0.03 0.10 int 0 > $DSDIR/big.profile
DATASETS=("$ORIG/examples/simulated_data.profile" $DSDIR/*.profile)

bytecmp() { diff -q "$1" "$2" >/dev/null 2>&1; }
treecmp() { python "$CMP" "$1" "$2" >/dev/null 2>&1; }
topocmp() { python "$CMP" "$1" "$2" --topo >/dev/null 2>&1; }

score() { # kind desc oracle rust  (kind: byte|tree|topo)
  local kind="$1" desc="$2" o="$3" r="$4"
  if [ ! -s "$o" ]; then skip=$((skip+1)); return; fi
  if bytecmp "$o" "$r"; then byte=$((byte+1)); return; fi
  if [ "$kind" = tree ] && treecmp "$o" "$r"; then equiv=$((equiv+1)); return; fi
  if [ "$kind" = topo ] && topocmp "$o" "$r"; then equiv=$((equiv+1)); return; fi
  fail=$((fail+1)); FAILED+=("$desc"); echo "  ❌ $desc"
}

for pf in "${DATASETS[@]}"; do
  b=$(basename "$pf" .profile)
  # ---- distance: matrix × handler ----
  for x in symmetric asymmetric blockwise; do
    for y in 0 1 2 3; do
      hm=(pair_delete complete_delete as_allele absolute_distance); hmn=${hm[$y]}
      [ "$x" = blockwise ] && { [ $y -ne 0 ] && continue; oa="matrix_type=blockwise handle_missing=0.01"; ra="-x blockwise"; } || { oa="matrix_type=$x handle_missing=$hmn"; ra="-x $x -y $y"; }
      python $ORUN $MP "$pf" distance $oa n_proc=1 2>/dev/null > $TMP/o.txt
      $G -p "$pf" -m distance $ra > $TMP/r.txt 2>/dev/null
      score byte "distance $b $x $hmn" $TMP/o.txt $TMP/r.txt
    done
  done
  # ---- MSTree symmetric: heuristic × handler (byte) ----
  for t in eBurst harmonic; do for y in 0 2 3; do
    hm=(pair_delete _ as_allele absolute_distance); hmn=${hm[$y]}
    python $ORUN $MP "$pf" MSTree matrix_type=symmetric heuristic=$t handle_missing=$hmn branch_recraft=False n_proc=1 2>/dev/null > $TMP/o.txt
    $G -p "$pf" -m MSTree -x symmetric -t $t -y $y > $TMP/r.txt 2>/dev/null
    score byte "MSTree-sym $b $t $hmn" $TMP/o.txt $TMP/r.txt
  done; done
  # ---- MSTree asymmetric (no recraft): tree ----
  for t in harmonic eBurst; do
    python $ORUN $MP "$pf" MSTree matrix_type=asymmetric heuristic=$t handle_missing=pair_delete branch_recraft=False n_proc=1 2>/dev/null > $TMP/o.txt
    $G -p "$pf" -m MSTree -x asymmetric -t $t -y 0 > $TMP/r.txt 2>/dev/null
    score tree "MSTree-asym $b $t" $TMP/o.txt $TMP/r.txt
  done
  # ---- MSTreeV2 (default): byte, else tree ----
  python $ORUN $MP "$pf" MSTreeV2 n_proc=1 2>/dev/null > $TMP/o.txt
  $G -p "$pf" -m MSTreeV2 > $TMP/r.txt 2>/dev/null
  score tree "MSTreeV2 $b" $TMP/o.txt $TMP/r.txt
  # ---- NJ family: tree (RF=0) ----
  for meth in NJ RapidNJ; do
    python $ORUN $MP "$pf" $meth n_proc=1 2>/dev/null > $TMP/o.txt
    $G -p "$pf" -m $meth > $TMP/r.txt 2>/dev/null
    score topo "$meth $b" $TMP/o.txt $TMP/r.txt
  done
  # ---- -n/--n_proc invariance (rust output must not depend on threads) ----
  $G -p "$pf" -m MSTreeV2 -n 1 > $TMP/n1.txt 2>/dev/null
  $G -p "$pf" -m MSTreeV2 -n 8 > $TMP/n8.txt 2>/dev/null
  if bytecmp $TMP/n1.txt $TMP/n8.txt; then ninv=$((ninv+1)); else nfail=$((nfail+1)); echo "  ❌ n_proc invariance $b"; fi
done

# --- input-format checks (oracle can't do FASTA: upstream `del part` crash) ---
fmtok=0; fmtfail=0
ex=$ORIG/examples/simulated_data.profile
gzip -c "$ex" > $TMP/ex.profile.gz
$G -p $TMP/ex.profile.gz -m MSTreeV2 > $TMP/gz.txt 2>/dev/null
$G -p "$ex" -m MSTreeV2 > $TMP/pl.txt 2>/dev/null
if bytecmp $TMP/gz.txt $TMP/pl.txt; then fmtok=$((fmtok+1)); else fmtfail=$((fmtfail+1)); echo "  ❌ gzip != plain"; fi
# FASTA: rust must produce a non-empty NEWICK where the reference crashes
python - > $TMP/aln.fasta <<'PY'
import random; random.seed(7)
L,N=100,15; base=[random.choice("ACGT") for _ in range(L)]
for i in range(N):
    s=[c if random.random()>=0.08 else random.choice("ACGT") for c in base]
    print(f">s{i}"); print("".join(s))
PY
$G -p $TMP/aln.fasta -m MSTreeV2 > $TMP/fa.txt 2>/dev/null
if [ -s $TMP/fa.txt ] && grep -q ';' $TMP/fa.txt; then fmtok=$((fmtok+1)); else fmtfail=$((fmtfail+1)); echo "  ❌ FASTA produced no tree"; fi

echo "======================================================"
echo "REGRESSION SUMMARY"
echo "  byte-identical : $byte"
echo "  RF=0 equivalent: $equiv"
echo "  failed         : $fail"
echo "  skipped(oracle): $skip"
echo "  n_proc invariant: $ninv (fail $nfail)"
echo "  input formats   : $fmtok ok (fail $fmtfail)  [gzip==plain, FASTA runs]"
[ -n "${FAILED[*]:-}" ] && printf '  FAIL: %s\n' "${FAILED[@]}"
echo "======================================================"
exit $((fail + nfail + fmtfail))
