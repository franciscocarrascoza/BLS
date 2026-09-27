#!/usr/bin/env bash
# Created: 2026-09-27T18:55+02:00 | by: CC audit (manager) | purpose: the origin-offset tools bls_sizefloor and bls_e7 after WP6b (finding S2-3: full 3-D offset set; HARN-10: E7 timing convention; PBC path)
#
# Guards finding / decision IDs:
#   OFF-SET    the offset set must be named (--offsets halton64|diagonal8); indices outside the set are refused
#   OFF-GT     the DFS ground truth (grid, occupied voxels, components, size histogram) is identical at all 64
#              halton64 offsets: offsets move the probe lattice, never the atoms
#   OFF-3D     halton64 offset 0 is (0,0,0); the 64 origins are distinct, lie in one conventional cell of the
#              probe lattice (cubic: 0 <= o_i < a) and are not confined to the diagonal
#   OFF-DIAG   diagonal8 offset k (lattice at (-k,-k,-k)) gives the counts of the pre-audit mechanism: the atoms
#              of the fixed cell (--dump-pdb) translated by (k,k,k) voxels and run through bls_analyze
#   OFF-PBC    PBC xyz: the deck's cell is used as it is (nx = ceil(L/h)), origins lie in the commensurate cell,
#              a hexagonal probe lattice is refused (D10), --dump-pdb is refused
#   E7-TIME    bls_e7: elapsed <= total for BLS and DFS (labelling-only timers, D9), recall = detected/present,
#              and the same BLS/DFS counts as bls_sizefloor at the same deck and offset
#
# Usage: test_audit_offset_tools.sh /path/to/bls_analyze   (bls_sizefloor and bls_e7 are taken from its directory)
set -uo pipefail
BIN="${1:?usage: $0 /path/to/bls_analyze}"
SF="$(dirname "$BIN")/bls_sizefloor"
E7="$(dirname "$BIN")/bls_e7"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/bls_audit_off.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

python3 - "$WORK" <<'PY'
import os, random, sys
w = sys.argv[1]
rnd = random.Random(2027)
# 40 small, well separated "molecules" (1-4 atoms, O/H) in a 60 A cell: small enough that a sparse probe
# lattice misses some of them at some offsets, so the offsets matter.
atoms = []
centres = []
while len(centres) < 40:
    c = [rnd.uniform(4.0, 56.0) for _ in range(3)]
    if all(sum((a - b) ** 2 for a, b in zip(c, d)) > 7.0 ** 2 for d in centres):
        centres.append(c)
for c in centres:
    for k in range(rnd.randint(1, 4)):
        atoms.append(("O" if k == 0 else "H", [c[0] + 1.0 * k, c[1] + 0.5 * (k % 2), c[2]]))
with open(os.path.join(w, "mol.pdb"), "w") as f:
    f.write("REMARK   1 CREATED by test_audit_offset_tools.sh (CC audit)\n")
    f.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n" % (60.0, 60.0, 60.0, 90.0, 90.0, 90.0, "P 1", 1))
    for i, (el, (x, y, z)) in enumerate(atoms, 1):
        f.write("ATOM  %5d  %-2s  HOH A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (i, el, i, x, y, z, el))
    f.write("END\n")
PY

deck() {  # name box pbc
  cat > "$WORK/$1" <<DECK
# Created by test_audit_offset_tools.sh (CC audit): dNN = 2.0*2.4/2.0 = 2.4 voxels (sparse)
BLS ...
  GROUP ATOMS=all
  BOX $2
  PBC $3
  GRID_SPACING 2.0
  CONNECTIVITY 6
  SKIP 1
  ALPHA 2.0
  DNN 0
  RADII 1.5, 1.2
  CUTOFF 0.25
  OCCUPANCY ANY
  LATTICE cubic
  CENTERING F
... BLS
DECK
}
deck none.in AUTO none
deck pbc.in CELL xyz

FAILS=0; CHECKS=0
check() {  # claim ok message
  CHECKS=$((CHECKS + 1))
  if [[ "$2" == "1" ]]; then return 0; fi
  FAILS=$((FAILS + 1)); echo "  FAIL [$1] $3"
  eval "FAILED_$(echo "$1" | tr -c 'A-Za-z0-9' '_')=1"
}
field() { tr ' ' '\n' <<<"$1" | awk -F= -v k="$2" '$1==k {print $2}'; }   # CELL line, key

# OFF-SET
"$SF" "$WORK/mol.pdb" "$WORK/none.in" cubic F 0 >/dev/null 2>&1; rc=$?
check OFF-SET "$([[ $rc == 2 ]] && echo 1)" "bls_sizefloor without --offsets: rc=$rc (want 2)"
"$SF" "$WORK/mol.pdb" "$WORK/none.in" cubic F 64 --offsets halton64 >/dev/null 2>&1; rc=$?
check OFF-SET "$([[ $rc == 2 ]] && echo 1)" "halton64 offset 64 accepted: rc=$rc"
"$SF" "$WORK/mol.pdb" "$WORK/none.in" cubic F 8 --offsets diagonal8 >/dev/null 2>&1; rc=$?
check OFF-SET "$([[ $rc == 2 ]] && echo 1)" "diagonal8 offset 8 accepted: rc=$rc"
"$E7" "$WORK/mol.pdb" "$WORK/none.in" cubic F 2.0 2.0 0 >/dev/null 2>&1; rc=$?
check OFF-SET "$([[ $rc == 2 ]] && echo 1)" "bls_e7 without --offsets: rc=$rc (want 2)"

# OFF-GT, OFF-3D (PBC none) and OFF-PBC (PBC xyz) over the 64 halton offsets
for mode in none pbc; do
  : > "$WORK/sf_$mode.txt"
  for i in $(seq 0 63); do
    "$SF" "$WORK/mol.pdb" "$WORK/$mode.in" cubic F "$i" --offsets halton64 >> "$WORK/sf_$mode.txt" 2>>"$WORK/err_$mode.log" \
      || echo "FAILED offset $i" >> "$WORK/sf_$mode.txt"
  done
done
for mode in none pbc; do
  out="$(python3 - "$WORK/sf_$mode.txt" "$mode" <<'PY'
import sys
path, mode = sys.argv[1], sys.argv[2]
cells, cur = [], None
for line in open(path):
    if line.startswith("FAILED"):
        print("FAILED"); sys.exit()
    if line.startswith("CELL "):
        cur = dict(t.split("=", 1) for t in line.split()[1:]); cur["hist"] = []; cells.append(cur)
    elif line.startswith("HIST "):
        cur["hist"].append(tuple(line.split()[1:3]))
gt = {(c["nx"], c["ny"], c["nz"], c["occupied"], c["dfs_nclusters"], c["dfs_max"], tuple(c["hist"])) for c in cells}
orig = [tuple(float(v) for v in c["origin"].split(",")) for c in cells]
dnn = float(cells[0]["dnn_vox"])
a = dnn * 2 ** 0.5                 # cubic F: conventional edge = sqrt(2) dNN
inside = all(-1e-9 <= o[k] < a + 1e-9 for o in orig for k in range(3))
offdiag = any(abs(o[0] - o[1]) > 1e-6 or abs(o[1] - o[2]) > 1e-6 for o in orig)
found = {c["bls_nclusters"] for c in cells}
print(len(cells), len(gt), int(orig[0] == (0.0, 0.0, 0.0)), len(set(orig)), int(inside), int(offdiag),
      len(found), cells[0]["nx"], cells[0]["pbc"], "%.6f" % dnn,
      int(all(c["grid_ok"] == "1" and c["merged"] == "0" and c["split"] == "0" and c["truncated"] == "0" for c in cells)))
PY
)"
  read -r n ngt zero ndist inside offdiag nfound nx pbc dnn clean <<<"$out"
  claim=OFF-GT; [[ $mode == pbc ]] && claim=OFF-PBC
  check "$claim" "$([[ $n == 64 && $ngt == 1 && $clean == 1 ]] && echo 1)" \
        "$mode: $n cells, $ngt ground-truth variants (want 64, 1), grid/refinement clean=$clean $(head -c 300 "$WORK/err_$mode.log")"
  check OFF-3D "$([[ $zero == 1 && $ndist == 64 && $inside == 1 && $offdiag == 1 ]] && echo 1)" \
        "$mode: offset0 at origin=$zero, distinct=$ndist, inside cell=$inside, off-diagonal=$offdiag"
  check OFF-3D "$([[ $nfound -gt 1 ]] && echo 1)" "$mode: BLS count identical at all 64 offsets ($nfound value): the offsets do not move the lattice"
  if [[ $mode == pbc ]]; then
    check OFF-PBC "$([[ $nx == 30 && $pbc == xyz ]] && echo 1)" "PBC: nx=$nx pbc=$pbc (want 30 = 60 A / 2.0 A, xyz)"
    check OFF-PBC "$(awk -v d="$dnn" 'BEGIN{print (d > 0 && d <= 2.4 + 1e-9) ? 1 : 0}')" "PBC effective dNN $dnn not in (0, 2.4]"
  fi
done
"$SF" "$WORK/mol.pdb" "$WORK/pbc.in" hexagonal P 0 --offsets halton64 >/dev/null 2>"$WORK/err_hex.log"; rc=$?
check OFF-PBC "$([[ $rc != 0 ]] && grep -qi cubic "$WORK/err_hex.log" && echo 1)" "hexagonal under PBC xyz: rc=$rc $(cat "$WORK/err_hex.log")"
"$SF" "$WORK/mol.pdb" "$WORK/pbc.in" cubic F 0 --offsets halton64 --dump-pdb "$WORK/x.pdb" >/dev/null 2>&1; rc=$?
check OFF-PBC "$([[ $rc == 2 ]] && echo 1)" "--dump-pdb under PBC xyz: rc=$rc (want 2)"

# OFF-DIAG: the pre-audit mechanism replayed through bls_analyze
"$SF" "$WORK/mol.pdb" "$WORK/none.in" cubic F 0 --offsets diagonal8 --dump-pdb "$WORK/fixed.pdb" >/dev/null 2>&1
diag_ok=1; diag_msg=""
for k in $(seq 0 7); do
  line="$("$SF" "$WORK/mol.pdb" "$WORK/none.in" cubic F "$k" --offsets diagonal8 2>/dev/null | head -1)"
  python3 - "$WORK/fixed.pdb" "$WORK/shift.pdb" "$k" <<'PY'
import sys
src, dst, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(src) as f, open(dst, "w") as g:
    for l in f:
        if l.startswith("ATOM"):
            x, y, z = (float(l[30:38]) + 2.0 * k, float(l[38:46]) + 2.0 * k, float(l[46:54]) + 2.0 * k)
            l = l[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + l[54:]
        g.write(l)
PY
  "$BIN" --system "$WORK/shift.pdb" --format pdb --conf "$WORK/none.in" --algo bls --quiet --out "$WORK/shift.csv" 2>/dev/null
  got="$(awk -F, 'NR==2 {print $12","$13","$10}' "$WORK/shift.csv")"
  want="$(field "$line" bls_nclusters),$(field "$line" bls_max),$(field "$line" bls_seeds)"
  [[ "$got" == "$want" ]] || { diag_ok=0; diag_msg+=" k=$k analyze=$got probe=$want"; }
done
check OFF-DIAG "$diag_ok" "diagonal8 differs from the atom-translation mechanism:$diag_msg"

# E7-TIME
for mode in none pbc; do
  e7="$("$E7" "$WORK/mol.pdb" "$WORK/$mode.in" cubic F 2.0 2.0 5 --offsets halton64 2>&1)"
  sf="$(sed -n '/^CELL/p' "$WORK/sf_$mode.txt" | sed -n 6p)"
  ok="$(python3 - "$e7" "$sf" <<'PY'
import sys
e = dict(t.split("=", 1) for t in sys.argv[1].split()[1:])
s = dict(t.split("=", 1) for t in sys.argv[2].split()[1:])
f = lambda k: float(e[k])
ok = (f("bls_elapsed_ms") <= f("bls_total_ms") and f("dfs_elapsed_ms") <= f("dfs_total_ms")
      and abs(f("recall") - int(e["detected"]) / int(e["present"])) < 1e-6
      and e["detected"] == s["bls_nclusters"] and e["present"] == s["dfs_nclusters"]
      and e["bls_seeds"] == s["bls_seeds"] and e["origin"] == s["origin"] and e["grid_ok"] == "1")
print(int(ok))
PY
)"
  check E7-TIME "$ok" "$mode: bls_e7 '$e7' vs bls_sizefloor '$sf'"
done

echo
echo "== test_audit_offset_tools: per-claim summary =="
for c in OFF-SET OFF-GT OFF-3D OFF-DIAG OFF-PBC E7-TIME; do
  v="FAILED_$(echo "$c" | tr -c 'A-Za-z0-9' '_')"
  if [[ -n "${!v:-}" ]]; then echo "  FAIL $c"; else echo "  PASS $c"; fi
done
echo "== test_audit_offset_tools: $CHECKS checks, $FAILS failed =="
[[ $FAILS == 0 ]]
