#!/usr/bin/env bash
# Created: 2026-09-26T13:40+02:00 | by: CC audit (manager) | purpose: multi-frame (trajectory) and same-grid tests through the real bls_analyze binary, for every reported algorithm
#
# Guards claims-ledger IDs:
#   TRJ-ID    N identical frames give identical per-frame results, equal to the single-frame run
#   TRJ-ORD   frame order A,B vs B,A gives identical per-frame results (no state leaks between frames)
#   TRJ-CNT   processed frame count == input frame count
#   TRJ-BOX   per-frame box honoured: (a) BOX AUTO without CRYST1 -> each frame's own bounding box;
#             (b) per-MODEL CRYST1 records with different cells -> each frame's own cell
#   S-5       "all methods label the same grid" (§3.1 l.426-428): NX,NY,NZ identical across methods per frame
#             (proxy for the grid fingerprint; an occupancy hash needs a code change — see STOP A fix plan)
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

ALGOS=(bls traditional_dfs cc3d_optimized rle_ccl_optimized gcbd skip_dfs dbscan hierarchical kmeans vccs_optimized hdbscan)
extra_args() {
  case "$1" in
    kmeans) echo "--algo-k 4" ;;
    skip_dfs) echo "--algo-skip 1" ;;
    *) echo "" ;;
  esac
}
run() {  # algo pdb out
  # shellcheck disable=SC2046
  "$BIN" --system "$WORK/$2" --format pdb --conf "$WORK/deck.in" --algo "$1" $(extra_args "$1") --quiet --out "$WORK/$3" 2>>"$WORK/stderr.log"
}

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

declare -A DIMS_A DIMS_B
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
done

# same grid across methods (dims proxy)
for a in "${ALGOS[@]}"; do
  check S-5 "$([[ "${DIMS_A[$a]}" == "${DIMS_A[bls]}" && "${DIMS_B[$a]}" == "${DIMS_B[bls]}" ]] && echo 1)" \
        "$a: grid dims differ from BLS (${DIMS_A[$a]} vs ${DIMS_A[bls]})"
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
for c in TRJ-CNT TRJ-ID TRJ-ORD TRJ-BOX S-5; do
  v="FAILED_$(echo "$c" | tr -c 'A-Za-z0-9' '_')"
  if [[ -n "${!v:-}" ]]; then echo "  FAIL $c"; else echo "  PASS $c"; fi
done
if [[ -s "$WORK/stderr.log" ]]; then echo "  (stderr of bls_analyze, first lines:)"; head -n 5 "$WORK/stderr.log"; fi
echo "== test_audit_trajectory: $CHECKS checks, $FAILS failed =="
[[ $FAILS -eq 0 ]]
