#!/usr/bin/env bash
# Created: 2026-09-26T17:48+02:00 | by: CC audit (manager) | purpose: PBC end to end through the real bls_analyze binary (WP2: deck PBC xyz + BOX CELL, commensurate probe lattice, periodic refinement, LATTICE_ORIGIN)
#
# Guards decision / claim IDs:
#   E2E-PBC    a chain of atoms crossing the x face of a 40 A CRYST1 cell is ONE BLS cluster with PBC xyz
#              and TWO with PBC none (same cell, BOX CELL); a separate blob stays separate in both, and the
#              joined chain is the largest cluster only under PBC (sizes are voxels, not atoms)
#   E2E-META   under PBC xyz the CSV reports pbc=xyz, an effective dNN <= the requested one (D10) and probes > 0
#   E2E-ORIG   LATTICE_ORIGIN is accepted and applied: at full coverage the answer does not change; on a sparse
#              lattice (dNN 2.4 voxels) the seed set moves while the periodic probe count n_x n_y n_z |offsets|
#              stays fixed
#   E2E-EXACT  the exact comparison methods with periodic adjacency (WP3: traditional_dfs, skip_dfs at unit
#              stride, gcbd, cc3d_optimized, rle_ccl_optimized) agree with BLS under PBC xyz and PBC none
#   E2E-REFUSE methods without periodic adjacency refuse PBC xyz (non-zero exit, message): withdrawn cc3d,
#              rle_ccl, vccs, dbscan, hdbscan; k-means, hierarchical; a hexagonal probe lattice under
#              PBC xyz is refused (D10); vccs_optimized (WP4) runs under PBC and covers the occupied set
#
# Usage: test_audit_pbc_e2e.sh /path/to/bls_analyze
set -uo pipefail
BIN="${1:?usage: $0 /path/to/bls_analyze}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/bls_audit_pbc.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

python3 - "$WORK" <<'EOF'
import sys, os
w = sys.argv[1]
# Chain along x across the face at x = 0 / 40 (2 A steps = one voxel at GRID_SPACING 2.0), plus a
# compact blob in the middle of the cell. Atoms outside [0, 40) are written as they are; the voxeliser
# folds them into the cell when periodic.
chain = [(x, 20.0, 20.0) for x in (29.0, 31.0, 33.0, 35.0, 37.0, 39.0, 1.0, 3.0, 5.0, 7.0, 9.0, 11.0)]
blob = [(18.0 + dx, 8.0 + dy, 30.0 + dz) for dx in (0, 2) for dy in (0, 2) for dz in (0, 2)]
with open(os.path.join(w, "cross.pdb"), "w") as f:
    f.write("REMARK   1 CREATED by test_audit_pbc_e2e.sh (CC audit)\n")
    f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n" % (40.0, 40.0, 40.0, 90.0, 90.0, 90.0, "P 1", 1))
    for i, (x, y, z) in enumerate(chain + blob, 1):
        f.write("ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00           O\n" % (i, i, x, y, z))
    f.write("END\n")
EOF

deck() {  # name pbc lattice [extra line]
  cat > "$WORK/$1" <<EOF
# Created by test_audit_pbc_e2e.sh (CC audit): dNN = 0.5*2.4/2.0 = 0.6 voxel, so BLS probes every voxel
BLS ...
  GROUP ATOMS=all
  BOX CELL
  PBC $2
  GRID_SPACING 2.0
  CONNECTIVITY 6
  SKIP 1
  ALPHA 0.5
  DNN 0
  RADII 1.5, 1.2
  CUTOFF 0.25
  OCCUPANCY ANY
  LATTICE $3
  CENTERING F
${4:-}
... BLS
EOF
}
deck pbc.in xyz cubic
deck none.in none cubic
deck pbc_orig.in xyz cubic "  LATTICE_ORIGIN 0.5 0.25 0.75"
deck pbc_hex.in xyz hexagonal
sed -e 's/ALPHA 0.5/ALPHA 2.0/' "$WORK/pbc.in" > "$WORK/sparse.in"            # dNN 2.4 voxels
sed -e 's/ALPHA 0.5/ALPHA 2.0/' "$WORK/pbc_orig.in" > "$WORK/sparse_orig.in"

FAILS=0; CHECKS=0
check() {  # claim ok message
  CHECKS=$((CHECKS + 1))
  if [[ "$2" == "1" ]]; then return 0; fi
  FAILS=$((FAILS + 1)); echo "  FAIL [$1] $3"
  eval "FAILED_$(echo "$1" | tr -c 'A-Za-z0-9' '_')=1"
}
run() {  # deck out algo
  "$BIN" --system "$WORK/cross.pdb" --format pdb --conf "$WORK/$1" --algo "$3" --quiet --out "$WORK/$2" 2>>"$WORK/stderr_$2.log"
}
col() { awk -F, -v c="$2" 'NR==2 {print $c}' "$WORK/$1"; }   # frame 1, column c

run pbc.in pbc.csv bls;       rc_pbc=$?
run none.in none.csv bls;     rc_none=$?
run pbc_orig.in orig.csv bls; rc_orig=$?
check E2E-PBC "$([[ $rc_pbc == 0 && $(col pbc.csv 12) == 2 ]] && echo 1)" \
      "PBC xyz: nclusters = $(col pbc.csv 12) (want 2: chain joined across the face + blob) rc=$rc_pbc $(cat "$WORK/stderr_pbc.csv.log")"
check E2E-PBC "$([[ $rc_none == 0 && $(col none.csv 12) == 3 ]] && echo 1)" \
      "PBC none: nclusters = $(col none.csv 12) (want 3: chain split at the face + blob) rc=$rc_none"
check E2E-PBC "$([[ $rc_pbc == 0 && $rc_none == 0 && $(col pbc.csv 13) -gt $(col none.csv 13) ]] && echo 1)" \
      "max_cluster PBC $(col pbc.csv 13) not above PBC none $(col none.csv 13) (the joined chain should dominate)"
check E2E-META "$([[ $(col pbc.csv 24) == xyz && $(col none.csv 24) == none ]] && echo 1)" \
      "pbc column: '$(col pbc.csv 24)' / '$(col none.csv 24)'"
check E2E-META "$(awk -v d="$(col pbc.csv 7)" 'BEGIN{print (d > 0 && d <= 0.6 + 1e-9) ? 1 : 0}')" \
      "PBC dNN_vox $(col pbc.csv 7) not in (0, 0.6]"
check E2E-META "$([[ $(col pbc.csv 26) -gt 0 && $(col none.csv 26) -gt 0 ]] && echo 1)" \
      "probes column: PBC $(col pbc.csv 26), none $(col none.csv 26)"
check E2E-ORIG "$([[ $rc_orig == 0 && $(col orig.csv 12) == 2 && $(col orig.csv 13) == $(col pbc.csv 13) ]] && echo 1)" \
      "LATTICE_ORIGIN under PBC: nclusters,max = $(col orig.csv 12),$(col orig.csv 13) (want 2,$(col pbc.csv 13)) rc=$rc_orig"
run sparse.in sparse.csv bls; run sparse_orig.in sparse_orig.csv bls
check E2E-ORIG "$([[ $(col sparse.csv 10) != "$(col sparse_orig.csv 10)" && $(col sparse.csv 26) == "$(col sparse_orig.csv 26)" && \
                     $(col sparse.csv 12) == 2 && $(col sparse_orig.csv 12) == 2 ]] && echo 1)" \
      "sparse PBC: seeds $(col sparse.csv 10) vs $(col sparse_orig.csv 10) with LATTICE_ORIGIN (want different), probes $(col sparse.csv 26) vs $(col sparse_orig.csv 26) (want equal), nclusters $(col sparse.csv 12)/$(col sparse_orig.csv 12) (want 2)"

# Exact methods with periodic adjacency: same clusters as BLS (full coverage) under both settings.
run_algo() {  # deck out algo extra...
  local d="$1" o="$2" a="$3"; shift 3
  "$BIN" --system "$WORK/cross.pdb" --format pdb --conf "$WORK/$d" --algo "$a" "$@" --quiet --out "$WORK/$o" 2>>"$WORK/stderr_$o.log"
}
for a in traditional_dfs skip_dfs gcbd cc3d_optimized rle_ccl_optimized; do
  extra=(); [[ $a == skip_dfs ]] && extra=(--algo-skip 1); [[ $a == cc3d_optimized ]] && extra=(--algo-connectivity 6)
  run_algo pbc.in "x_$a.csv" "$a" "${extra[@]}"; rc1=$?
  run_algo none.in "n_$a.csv" "$a" "${extra[@]}"; rc2=$?
  check E2E-EXACT "$([[ $rc1 == 0 && $rc2 == 0 && $(col "x_$a.csv" 12),$(col "x_$a.csv" 13) == $(col pbc.csv 12),$(col pbc.csv 13) && \
                        $(col "n_$a.csv" 12),$(col "n_$a.csv" 13) == $(col none.csv 12),$(col none.csv 13) && $(col "x_$a.csv" 24) == xyz ]] && echo 1)" \
        "$a: PBC xyz $(col "x_$a.csv" 12),$(col "x_$a.csv" 13) vs BLS $(col pbc.csv 12),$(col pbc.csv 13); none $(col "n_$a.csv" 12),$(col "n_$a.csv" 13) vs BLS $(col none.csv 12),$(col none.csv 13) (rc $rc1/$rc2) $(cat "$WORK/stderr_x_$a.csv.log")"
done

# Refusals: methods without periodic adjacency, hexagonal probe lattice under PBC (D10).
for a in cc3d rle_ccl vccs dbscan hdbscan kmeans hierarchical; do
  extra=(); [[ $a == kmeans ]] && extra=(--algo-k 2)
  run_algo pbc.in "refuse_$a.csv" "$a" "${extra[@]}"; rc=$?
  check E2E-REFUSE "$([[ $rc != 0 ]] && grep -q "PBC is not implemented" "$WORK/stderr_refuse_$a.csv.log" && echo 1)" \
        "$a ran under PBC xyz (rc=$rc)"
done
# vccs_optimized (WP4) runs under PBC and, as a segmentation, covers every occupied voxel.
run_algo pbc.in vccs.csv vccs_optimized --algo-connectivity 26; rc=$?
check E2E-REFUSE "$([[ $rc == 0 && $(col vccs.csv 24) == xyz && $(col vccs.csv 14) == $(col pbc.csv 14) ]] && echo 1)" \
      "vccs_optimized under PBC xyz: rc=$rc pbc='$(col vccs.csv 24)' covered $(col vccs.csv 14) of $(col pbc.csv 14) $(cat "$WORK/stderr_vccs.csv.log")"
run pbc_hex.in hex.csv bls; rc=$?
check E2E-REFUSE "$([[ $rc != 0 ]] && grep -qi "cubic" "$WORK/stderr_hex.csv.log" && echo 1)" \
      "hexagonal probe lattice ran under PBC xyz (rc=$rc): $(cat "$WORK/stderr_hex.csv.log")"

echo
echo "== test_audit_pbc_e2e: per-claim summary =="
for c in E2E-PBC E2E-META E2E-ORIG E2E-EXACT E2E-REFUSE; do
  v="FAILED_$(echo "$c" | tr -c 'A-Za-z0-9' '_')"
  if [[ -n "${!v:-}" ]]; then echo "  FAIL $c"; else echo "  PASS $c"; fi
done
echo "== test_audit_pbc_e2e: $CHECKS checks, $FAILS failed =="
[[ $FAILS == 0 ]]
