#!/usr/bin/env bash
# Created: 2026-09-26T16:38+02:00 | by: CC audit (manager) | purpose: capture the non-timing CSV columns of bls_analyze over a fixed matrix of ORIGINAL inputs, so a later binary can be checked bit for bit against a pre-fix one ("PBC none reproduces pre-audit outputs")
#
# Usage: capture_reference.sh /path/to/bls_analyze OUTDIR
# Writes OUTDIR/<case>.csv (columns 1-14 of the CSV: frame..refined_voxels; elapsed_ms, replicate and
# anything appended after them are dropped) and OUTDIR/cases.txt (one line per case: name + command).
# Compare two captures with: diff -r OUT_A OUT_B
set -uo pipefail
BIN="${1:?usage: $0 BIN OUTDIR}"
OUT="${2:?usage: $0 BIN OUTDIR}"
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"   # outer tree
SYS="$ROOT/systems"
W="$ROOT/results/E7_optimal/decks/opt_optimal.in"   # working deck (4.0 A, ALPHA 2.45, cubic F)
mkdir -p "$OUT"
: > "$OUT/cases.txt"
fail=0

args_for() {  # algo -> extra args used by the campaign runners (skip_dfs1 = unit stride)
  case "$1" in
    kmeans) echo "--algo kmeans --algo-k ${2:-1000}" ;;
    hierarchical) echo "--algo hierarchical --algo-threshold 4.0" ;;
    dbscan) echo "--algo dbscan --algo-eps 3.0 --algo-minpts 10" ;;
    hdbscan) echo "--algo hdbscan --algo-minclustersize 5 --algo-minsamples 5" ;;
    skip_dfs1) echo "--algo skip_dfs --algo-skip 1" ;;
    skip_dfs3) echo "--algo skip_dfs --algo-skip 3" ;;
    cc3d|cc3d_optimized) echo "--algo $1 --algo-connectivity 6" ;;
    *) echo "--algo $1" ;;
  esac
}
run_case() {  # name pdb deck algo [k]
  local name="$1" pdb="$2" deck="$3" algo="$4" k="${5:-}"
  local a; a="$(args_for "$algo" "$k")"
  echo "$name bls_analyze --system $pdb --conf $deck $a" >> "$OUT/cases.txt"
  # shellcheck disable=SC2086
  if ! "$BIN" --system "$pdb" --conf "$deck" $a --quiet --out "$OUT/$name.full.csv" 2>>"$OUT/stderr.log"; then
    echo "FAILED: $name" | tee -a "$OUT/stderr.log"; fail=1; return
  fi
  cut -d, -f1-14 "$OUT/$name.full.csv" > "$OUT/$name.csv" && rm -f "$OUT/$name.full.csv"
}

ALL=(bls traditional_dfs cc3d_optimized rle_ccl_optimized gcbd skip_dfs1 skip_dfs3 vccs_optimized hierarchical dbscan cc3d rle_ccl vccs)
for s in Ic Ih ld-asw vii; do
  for a in "${ALL[@]}" kmeans; do run_case "E1W_${s}_${a}" "$SYS/E1-lattice/$s.pdb" "$W" "$a" 1000; done
  run_case "E1P_${s}_bls" "$SYS/E1-lattice/$s.pdb" "$SYS/E1-lattice/opt_${s}.in" bls
  run_case "E1D35_${s}_bls" "$SYS/E1-lattice/$s.pdb" "$SYS/E1-lattice/opt.in" bls
  run_case "E1D35_${s}_dfs" "$SYS/E1-lattice/$s.pdb" "$SYS/E1-lattice/opt.in" traditional_dfs
done
for n in 1 100 1000; do
  for a in bls traditional_dfs cc3d_optimized; do run_case "E3_${n}_${a}" "$SYS/E3-Scaling/${n}Ic.pdb" "$SYS/E3-Scaling/opt.in" "$a"; done
done
for n in 100 250 750; do
  for a in bls traditional_dfs gcbd; do run_case "E4_${n}_${a}" "$SYS/E4-Density/den_nc${n}.pdb" "$SYS/E4-Density/opt.in" "$a"; done
done
for f in res_gap7 res_gap8 res_multi_gap5; do
  for a in "${ALL[@]}" hdbscan kmeans; do run_case "E5_${f}_${a}" "$SYS/E5-Resolution/$f.pdb" "$SYS/E5-Resolution/opt.in" "$a" 2; done
  run_case "E5W_${f}_bls" "$SYS/E5-Resolution/$f.pdb" "$W" bls
done
echo "cases: $(wc -l < "$OUT/cases.txt"), failures: $fail"
exit $fail
