#!/usr/bin/env bash
# Created: 2026-09-26T12:53+02:00 | by: CC audit (manager) | purpose: multi-frame (trajectory) and same-grid tests through the real bls_analyze binary, for every reported algorithm
#
# Guards claims-ledger IDs:
#   TRJ-ID    N identical frames give identical per-frame results, equal to the single-frame run
#   TRJ-ORD   frame order A,B vs B,A gives identical per-frame results (no state leaks between frames)
#   TRJ-CNT   processed frame count == input frame count
#   TRJ-BOX   per-frame box honoured: (a) BOX AUTO without CRYST1 -> each frame's own bounding box;
#             (b) per-MODEL CRYST1 records with different cells -> each frame's own cell
#   S-5       "all methods label the same grid" (§3.1 l.426-428): the grid fingerprint -- NX,NY,NZ, per-axis
#             voxel edge h_x,h_y,h_z, origin, periodic flags and SHA-256 of the occupancy bits -- identical
#             across all methods on every frame (fingerprint columns added by the audit, 2026-09-26)
#   TRJ-ORIG  frames with equal grid dimensions but different origins (a translated copy) are each gridded
#             at their own origin (pre-audit BLS re-used frame 1's box/origin whenever the dimensions matched)
#   D12       BOX MANUAL is honoured by every method (Defect 12), BOX CELL takes each MODEL's CRYST1
#
# Usage: test_audit_trajectory.sh /path/to/bls_analyze
set -uo pipefail
BIN="${1:?usage: $0 /path/to/bls_analyze}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/bls_audit_traj.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
PY=python3

# ---- inputs -----------------------------------------------------------------
# Frame A: 4 compact 8-atom blobs; frame B: 6 blobs over a larger extent; both deterministic.
$PY - "$WORK" <<'EOF'
import sys, os, math
w = sys.argv[1]
def blobs(centres, spacing=1.0):
    atoms = []
    for (cx, cy, cz) in centres:
        for dx in (0, spacing):
            for dy in (0, spacing):
                for dz in (0, spacing):
                    atoms.append((cx + dx, cy + dy, cz + dz))
    return atoms
A = blobs([(2, 2, 2), (12, 3, 2), (3, 13, 4), (14, 14, 12)])
B = blobs([(1, 1, 1), (18, 2, 3), (4, 20, 5), (21, 19, 17), (10, 10, 22), (25, 6, 9)])
def atom_line(i, x, y, z):
    return "ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           O\n" % (i, i // 8 + 1, x, y, z)
def cryst(a, b, c):
    return "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n" % (a, b, c, 90.0, 90.0, 90.0, "P 1", 1)
def write(name, frames, cells=None):
    with open(os.path.join(w, name), "w") as f:
        f.write("REMARK   1 CREATED by test_audit_trajectory.sh (CC audit)\n")
        for m, fr in enumerate(frames, 1):
            if cells: f.write(cryst(*cells[m - 1]))
            f.write("MODEL     %4d\n" % m)
            for i, (x, y, z) in enumerate(fr, 1):
                f.write(atom_line(i, x, y, z))
            f.write("ENDMDL\n")
        f.write("END\n")
write("A.pdb", [A]); write("B.pdb", [B])
write("AAA.pdb", [A, A, A]); write("AB.pdb", [A, B]); write("BA.pdb", [B, A])
# per-MODEL CRYST1: same atoms (A), two different cells -> grids must differ per frame
write("A_cellsmall.pdb", [A], [(30.0, 30.0, 30.0)])
write("A_cellbig.pdb", [A], [(44.0, 40.0, 36.0)])
write("A_twocells.pdb", [A, A], [(30.0, 30.0, 30.0), (44.0, 40.0, 36.0)])
# C = A translated by +1.0 A in x (half a voxel): same extents -> same grid dims, different origin
C = [(x + 1.0, y + 0.3, z + 0.6) for (x, y, z) in A]
write("C.pdb", [C]); write("AC.pdb", [A, C])
EOF

cat > "$WORK/deck.in" <<'EOF'
# Created by test_audit_trajectory.sh (CC audit): small deck, dNN = 0.5*2.4/2.0 = 0.6 voxel so BLS probes every voxel
BLS ...
  GROUP ATOMS=all
  BOX AUTO
  GRID_SPACING 2.0
  CONNECTIVITY 6
  SKIP 1
  ALPHA 0.5
  DNN 0
  RADII 1.5, 1.2
  CUTOFF 0.25
  OCCUPANCY ANY
  LATTICE cubic
  CENTERING F
... BLS
EOF

sed -e 's/^  BOX AUTO/  BOX MANUAL XLO -5 XHI 45 YLO -5 YHI 45 ZLO -5 ZHI 45/' "$WORK/deck.in" > "$WORK/deck_manual.in"
sed -e 's/^  BOX AUTO/  BOX CELL/' "$WORK/deck.in" > "$WORK/deck_cell.in"

ALGOS=(bls traditional_dfs cc3d_optimized rle_ccl_optimized gcbd skip_dfs dbscan hierarchical kmeans vccs_optimized hdbscan)
extra_args() {
  case "$1" in
    kmeans) echo "--algo-k 4" ;;
    skip_dfs) echo "--algo-skip 1" ;;
    *) echo "" ;;
  esac
}
run() {  # algo pdb out [deck]
  # shellcheck disable=SC2046
  "$BIN" --system "$WORK/$2" --format pdb --conf "$WORK/${4:-deck.in}" --algo "$1" $(extra_args "$1") --quiet --out "$WORK/$3" 2>>"$WORK/stderr.log"
}
# grid fingerprint of frame n: NX,NY,NZ + h_x,h_y,h_z,origin_x,origin_y,origin_z,pbc,occ_sha256 (cols 4-6, 18-25).
# Empty when the binary predates the fingerprint columns (then the S-5/TRJ-ORIG checks fail, as they should).
fp() { awk -F, -v n="$2" 'NR==n+1 {if (NF>=25) print $4","$5","$6","$18","$19","$20","$21","$22","$23","$24","$25; else print ""}' "$WORK/$1"; }

FAILS=0; CHECKS=0
check() {  # claim ok message
  CHECKS=$((CHECKS + 1))
  if [[ "$2" == "1" ]]; then return 0; fi
  FAILS=$((FAILS + 1)); echo "  FAIL [$1] $3"
  eval "FAILED_$(echo "$1" | tr -c 'A-Za-z0-9' '_')=1"
}
# comparable columns: NX,NY,NZ (4-6), seeds (10), seed_hits (11), nclusters (12), max_cluster (13), refined (14)
row() { awk -F, -v n="$2" 'NR==n+1 {print $4","$5","$6","$10","$11","$12","$13","$14}' "$WORK/$1"; }
nrows() { echo $(( $(wc -l < "$WORK/$1") - 1 )); }

declare -A DIMS_A DIMS_B FP_A FP_B FP_MAN1 FP_MAN2 FP_CELL1 FP_CELL2
for a in "${ALGOS[@]}"; do
  run "$a" A.pdb "${a}_A.csv"; run "$a" B.pdb "${a}_B.csv"
  run "$a" AAA.pdb "${a}_AAA.csv"; run "$a" AB.pdb "${a}_AB.csv"; run "$a" BA.pdb "${a}_BA.csv"
  run "$a" A_cellsmall.pdb "${a}_cs.csv"; run "$a" A_cellbig.pdb "${a}_cb.csv"; run "$a" A_twocells.pdb "${a}_tc.csv"

  check TRJ-CNT "$([[ $(nrows "${a}_AAA.csv") == 3 && $(nrows "${a}_AB.csv") == 2 ]] && echo 1)" "$a: frame count"
  rA=$(row "${a}_A.csv" 1); rB=$(row "${a}_B.csv" 1)
  for f in 1 2 3; do
    check TRJ-ID "$([[ "$(row "${a}_AAA.csv" $f)" == "$rA" ]] && echo 1)" "$a: AAA frame $f ($(row "${a}_AAA.csv" $f)) != single A ($rA)"
  done
  check TRJ-ORD "$([[ "$(row "${a}_AB.csv" 1)" == "$rA" && "$(row "${a}_AB.csv" 2)" == "$rB" ]] && echo 1)" \
        "$a: A,B per-frame != singles (got $(row "${a}_AB.csv" 1) | $(row "${a}_AB.csv" 2); want $rA | $rB)"
  check TRJ-ORD "$([[ "$(row "${a}_BA.csv" 1)" == "$rB" && "$(row "${a}_BA.csv" 2)" == "$rA" ]] && echo 1)" \
        "$a: B,A per-frame != singles (got $(row "${a}_BA.csv" 1) | $(row "${a}_BA.csv" 2); want $rB | $rA)"
  # (a) auto box follows each frame: A and B have different extents -> different dims
  check TRJ-BOX "$([[ "${rA%%,*}" != "${rB%%,*}" || "$(echo "$rA" | cut -d, -f2)" != "$(echo "$rB" | cut -d, -f2)" ]] && echo 1)" \
        "$a: frames with different extents got the same grid dims"
  # (b) per-MODEL CRYST1: frame 2 must use its own cell (dims of A_cellbig)
  cs=$(row "${a}_cs.csv" 1 | cut -d, -f1-3); cb=$(row "${a}_cb.csv" 1 | cut -d, -f1-3)
  t1=$(row "${a}_tc.csv" 1 | cut -d, -f1-3); t2=$(row "${a}_tc.csv" 2 | cut -d, -f1-3)
  check TRJ-BOX "$([[ "$t1" == "$cs" && "$t2" == "$cb" ]] && echo 1)" \
        "$a: per-MODEL CRYST1 not honoured (frame dims $t1 | $t2; single-cell runs $cs | $cb)"
  DIMS_A[$a]=$(echo "$rA" | cut -d, -f1-3); DIMS_B[$a]=$(echo "$rB" | cut -d, -f1-3)
  # TRJ-ORIG: A then C (same dims, other origin) -> each frame equals its single-frame run, fingerprint included
  run "$a" C.pdb "${a}_C.csv"; run "$a" AC.pdb "${a}_AC.csv"
  check TRJ-ORIG "$([[ "$(row "${a}_AC.csv" 2)" == "$(row "${a}_C.csv" 1)" && -n "$(fp "${a}_C.csv" 1)" && \
                       "$(fp "${a}_AC.csv" 2)" == "$(fp "${a}_C.csv" 1)" && "$(fp "${a}_AC.csv" 1)" == "$(fp "${a}_A.csv" 1)" ]] && echo 1)" \
        "$a: frame 2 of A,C differs from single C (row $(row "${a}_AC.csv" 2) vs $(row "${a}_C.csv" 1); fp '$(fp "${a}_AC.csv" 2)' vs '$(fp "${a}_C.csv" 1)')"
  # D12: BOX MANUAL and BOX CELL (per MODEL) through every method
  run "$a" AB.pdb "${a}_man.csv" deck_manual.in; run "$a" A_twocells.pdb "${a}_cell.csv" deck_cell.in
  FP_MAN1[$a]=$(fp "${a}_man.csv" 1); FP_MAN2[$a]=$(fp "${a}_man.csv" 2)
  FP_CELL1[$a]=$(fp "${a}_cell.csv" 1); FP_CELL2[$a]=$(fp "${a}_cell.csv" 2)
  FP_A[$a]=$(fp "${a}_AB.csv" 1); FP_B[$a]=$(fp "${a}_AB.csv" 2)
  check D12 "$([[ "$(row "${a}_man.csv" 1 | cut -d, -f1-3)" == "25,25,25" && "$(row "${a}_man.csv" 2 | cut -d, -f1-3)" == "25,25,25" ]] && echo 1)" \
        "$a: BOX MANUAL (50 A, h 2.0) not honoured: dims $(row "${a}_man.csv" 1 | cut -d, -f1-3)"
  check D12 "$([[ "$(row "${a}_cell.csv" 1 | cut -d, -f1-3)" == "15,15,15" && "$(row "${a}_cell.csv" 2 | cut -d, -f1-3)" == "22,20,18" ]] && echo 1)" \
        "$a: BOX CELL per MODEL not honoured: dims $(row "${a}_cell.csv" 1 | cut -d, -f1-3) | $(row "${a}_cell.csv" 2 | cut -d, -f1-3)"
done

# same grid across methods: full fingerprint (dims, h_i, origin, pbc, occupancy SHA-256), every frame, every deck
for a in "${ALGOS[@]}"; do
  check S-5 "$([[ "${DIMS_A[$a]}" == "${DIMS_A[bls]}" && "${DIMS_B[$a]}" == "${DIMS_B[bls]}" ]] && echo 1)" \
        "$a: grid dims differ from BLS (${DIMS_A[$a]} vs ${DIMS_A[bls]})"
  for v in FP_A FP_B FP_MAN1 FP_MAN2 FP_CELL1 FP_CELL2; do
    declare -n arr="$v"
    check S-5 "$([[ -n "${arr[$a]}" && "${arr[$a]}" == "${arr[bls]}" ]] && echo 1)" \
          "$a: grid fingerprint ($v) differs from BLS: '${arr[$a]}' vs '${arr[bls]}'"
    unset -n arr
  done
done

# exact methods: identical nclusters and max_cluster on each frame (BLS at dNN 0.6 voxel included)
for a in traditional_dfs cc3d_optimized rle_ccl_optimized gcbd skip_dfs; do
  for f in A B; do
    x=$(row "${a}_${f}.csv" 1 | cut -d, -f6-7); y=$(row "bls_${f}.csv" 1 | cut -d, -f6-7)
    check TRJ-ID "$([[ "$x" == "$y" ]] && echo 1)" "$a frame $f: nclusters,max ($x) != BLS ($y)"
  done
done

echo
echo "== test_audit_trajectory: per-claim summary =="
for c in TRJ-CNT TRJ-ID TRJ-ORD TRJ-BOX TRJ-ORIG D12 S-5; do
  v="FAILED_$(echo "$c" | tr -c 'A-Za-z0-9' '_')"
  if [[ -n "${!v:-}" ]]; then echo "  FAIL $c"; else echo "  PASS $c"; fi
done
if [[ -s "$WORK/stderr.log" ]]; then echo "  (stderr of bls_analyze, first lines:)"; head -n 5 "$WORK/stderr.log"; fi
echo "== test_audit_trajectory: $CHECKS checks, $FAILS failed =="
[[ $FAILS -eq 0 ]]
