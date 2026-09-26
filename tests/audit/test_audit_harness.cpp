// Created: 2026-09-26T16:45+02:00 | by: CC audit (manager) | purpose: harness unit tests added with the audit fixes: SHA-256, PDB CRYST1 round-trip per MODEL, deck parser strictness (BOX, PBC), shared grid builder
//
// Guards claims-ledger / finding IDs:
//   FP-SHA    SHA-256 used by the grid fingerprint matches FIPS 180-4 test vectors
//   PDB-RT    a PDB written in the audit's exact format (REMARK line 1, CRYST1 line 2, per MODEL) is
//             read back with every MODEL's own cell; REMARK records are ignored (brief §6.1)
//   HARN-2    BOX MANUAL parsing is strict (leading MANUAL token, all six bounds, unknown keys rejected)
//   PBC-DECK  PBC xyz|none parsed; PBC with BOX AUTO rejected; per-axis PBC rejected
//   D12       the shared grid builder honours MANUAL, CELL and AUTO exactly as specified

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "audit_common.hpp"
#include "bls/Options.hpp"
#include "common/Box.hpp"
#include "config/Parser.hpp"
#include "grid/GridSpec.hpp"
#include "io/TrajectoryReader.hpp"
#include "util/Sha256.hpp"

namespace {

std::string tmpPath(const std::string& name) {
  const char* d = std::getenv("TMPDIR");
  return std::string(d ? d : "/tmp") + "/bls_audit_harness_" + name;
}

void writeFile(const std::string& path, const std::string& text) {
  std::ofstream f(path);
  f << text;
}

std::string sha(const std::string& s) {
  return bls::sha256Hex(reinterpret_cast<const std::uint8_t*>(s.data()), s.size());
}

void testSha() {
  AUDIT_CHECK("FP-SHA", sha("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "empty");
  AUDIT_CHECK("FP-SHA", sha("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "abc");
  AUDIT_CHECK("FP-SHA",
              sha("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
                  "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
              "448-bit message");
  AUDIT_CHECK("FP-SHA", sha(std::string(1000000, 'a')) ==
                            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
              "one million 'a'");
  // Occupancy fingerprint: bit-packed, so any non-zero byte is 1 and the order matters.
  std::vector<std::uint8_t> a{0, 1, 0, 0, 0, 0, 0, 0, 1}, b{0, 7, 0, 0, 0, 0, 0, 0, 1}, c{1, 0, 0, 0, 0, 0, 0, 0, 1};
  AUDIT_CHECK("FP-SHA", bls::occupancySha256(a) == bls::occupancySha256(b), "non-zero bytes are one bit");
  AUDIT_CHECK("FP-SHA", bls::occupancySha256(a) != bls::occupancySha256(c), "bit order is significant");
}

std::string cryst1(double a, double b, double c, double al, double be, double ga) {
  char buf[128];  // the exact format of brief §6.1
  std::snprintf(buf, sizeof buf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d", a, b, c, al, be, ga, "P 1", 1);
  return buf;
}

void testPdbRoundTrip() {
  const std::string path = tmpPath("rt.pdb");
  std::string t;
  t += "REMARK   1 CREATED 2026-09-26T16:45+02:00 CC AUDIT BOX_SOURCE=test CRYST1 99 99 99\n";
  t += cryst1(30.0, 31.5, 32.25, 90.0, 90.0, 90.0) + "\n";
  t += "MODEL        1\n";
  t += "ATOM      1  O   HOH A   1       1.000   2.000   3.000  1.00  0.00           O\n";
  t += "ENDMDL\n";
  t += cryst1(44.0, 40.0, 36.0, 90.0, 100.0, 110.0) + "\n";  // triclinic second frame
  t += "MODEL        2\n";
  t += "ATOM      1  O   HOH A   1       4.000   5.000   6.000  1.00  0.00           O\n";
  t += "ENDMDL\nEND\n";
  writeFile(path, t);

  std::string err;
  auto reader = bls::makeTrajectoryReader(path, "pdb", err);
  AUDIT_CHECK("PDB-RT", reader != nullptr, "reader: " + err);
  if (!reader) return;
  bls::Frame f1, f2, f3;
  const bool r1 = reader->read(f1, err), r2 = reader->read(f2, err), r3 = reader->read(f3, err);
  AUDIT_CHECK("PDB-RT", r1 && r2 && !r3, "expected exactly two frames");
  const bls::Mat3 e1 = bls::buildTriclinicBox(30.0, 31.5, 32.25, 90.0, 90.0, 90.0);
  const bls::Mat3 e2 = bls::buildTriclinicBox(44.0, 40.0, 36.0, 90.0, 100.0, 110.0);
  auto same = [](const bls::Mat3& x, const bls::Mat3& y) {
    for (int c = 0; c < 3; ++c)
      for (int r = 0; r < 3; ++r)
        if (std::fabs(x.cols[c][r] - y.cols[c][r]) > 1e-9) return false;
    return true;
  };
  AUDIT_CHECK("PDB-RT", same(f1.box, e1), "frame 1 cell not read back (REMARK must be ignored)");
  AUDIT_CHECK("PDB-RT", same(f2.box, e2), "frame 2 must carry its own MODEL's CRYST1");
  AUDIT_CHECK("PDB-RT", f1.natoms == 1 && std::fabs(f1.xyz[0].z - 3.0) < 1e-12, "frame 1 atoms");
  AUDIT_CHECK("PDB-RT", f2.natoms == 1 && std::fabs(f2.xyz[0].x - 4.0) < 1e-12, "frame 2 atoms");
  std::remove(path.c_str());
}

bool parseDeck(const std::string& boxLine, const std::string& extra, bls::BLSConfig& cfg,
               std::string& err) {
  const std::string path = tmpPath("deck.in");
  writeFile(path, "BLS ...\n  " + boxLine + "\n  GRID_SPACING 2.0\n  LATTICE cubic\n  CENTERING F\n" +
                      extra + "... BLS\n");
  cfg = bls::BLSConfig{};
  bls::Parser p;
  const bool ok = p.parseFile(path, cfg, err);
  std::remove(path.c_str());
  return ok;
}

void testParser() {
  bls::BLSConfig cfg;
  std::string err;
  AUDIT_CHECK("HARN-2", parseDeck("BOX MANUAL XLO -5 XHI 45 YLO 0 YHI 10 ZLO 1 ZHI 2", "", cfg, err) &&
                            cfg.boxMode == bls::BoxMode::Manual && cfg.manualBox.xlo == -5 &&
                            cfg.manualBox.xhi == 45 && cfg.manualBox.zhi == 2,
              "BOX MANUAL with leading MANUAL token: " + err);
  AUDIT_CHECK("HARN-2", parseDeck("BOX XLO 0 XHI 1 YLO 0 YHI 1 ZLO 0 ZHI 1", "", cfg, err) &&
                            cfg.boxMode == bls::BoxMode::Manual && cfg.manualBox.yhi == 1,
              "BOX without MANUAL token: " + err);
  AUDIT_CHECK("HARN-2", !parseDeck("BOX MANUAL XLO 0 XHI 1 YLO 0 YHI 1", "", cfg, err), "missing bounds accepted");
  AUDIT_CHECK("HARN-2", !parseDeck("BOX MANUAL XLO 0 XHI 1 YLO 0 YHI 1 ZLO 0 ZZZ 1", "", cfg, err),
              "unknown key accepted");
  AUDIT_CHECK("HARN-2", !parseDeck("BOX MANUAL XLO 0 XHI 1 XLO 0 YHI 1 ZLO 0 ZHI 1", "", cfg, err),
              "repeated key accepted");
  AUDIT_CHECK("PBC-DECK", parseDeck("BOX CELL", "  PBC xyz\n", cfg, err) && cfg.boxMode == bls::BoxMode::Cell &&
                              cfg.pbc.all(),
              "BOX CELL + PBC xyz: " + err);
  AUDIT_CHECK("PBC-DECK", parseDeck("BOX AUTO", "", cfg, err) && !cfg.pbc.any(), "default must be PBC none");
  AUDIT_CHECK("PBC-DECK", !parseDeck("BOX AUTO", "  PBC xyz\n", cfg, err), "PBC with BOX AUTO accepted");
  AUDIT_CHECK("PBC-DECK", !parseDeck("BOX CELL", "  PBC xy\n", cfg, err), "per-axis PBC accepted");
  AUDIT_CHECK("PBC-DECK", parseDeck("BOX CELL", "  LATTICE_ORIGIN 0.5 1 1.5\n", cfg, err) &&
                              cfg.latticeOrigin[0] == 0.5 && cfg.latticeOrigin[2] == 1.5,
              "LATTICE_ORIGIN: " + err);
}

void testGridSpec() {
  bls::Frame fr;
  fr.xyz = {bls::Vec3{1.0, 2.0, 3.0}, bls::Vec3{11.0, 12.0, 13.0}};
  fr.natoms = 2;
  bls::BLSConfig cfg;
  cfg.gridSpacing = 2.0;
  bls::GridSpec s;
  std::string err;

  // AUTO without a cell: bounding box + 2h padding, origin at the minimum corner.
  AUDIT_CHECK("D12", bls::deriveGridSpec(cfg, fr, nullptr, s, err) && s.nx == 7 && s.ny == 7 && s.nz == 7 &&
                         s.origin.x == 1.0 && s.origin.z == 3.0 && !s.pbc.any(),
              "AUTO bounding box: " + err);
  // AUTO with a usable cell: the cell, origin 0.
  fr.box = bls::buildTriclinicBox(30.0, 30.0, 30.0, 90, 90, 90);
  AUDIT_CHECK("D12", bls::deriveGridSpec(cfg, fr, nullptr, s, err) && s.nx == 15 && s.origin.x == 0.0,
              "AUTO with cell");
  // MANUAL wins over the cell.
  cfg.boxMode = bls::BoxMode::Manual;
  cfg.manualBox = bls::ManualBox{-5, 45, 0, 20, 0, 21};
  AUDIT_CHECK("D12", bls::deriveGridSpec(cfg, fr, nullptr, s, err) && s.nx == 25 && s.ny == 10 && s.nz == 11 &&
                         s.origin.x == -5.0,
              "MANUAL box");
  AUDIT_CHECK("D12", std::fabs(s.h(2) - 21.0 / 11.0) < 1e-12 && s.h(2) <= 2.0, "per-axis h <= nominal");
  // CELL: the frame's cell; error without one.
  cfg.boxMode = bls::BoxMode::Cell;
  cfg.pbc = bls::PeriodicAxes{true, true, true};
  AUDIT_CHECK("D12", bls::deriveGridSpec(cfg, fr, nullptr, s, err) && s.nx == 15 && s.pbc.all() &&
                         s.periodicity() == bls::BoxPeriodicity::Periodic,
              "CELL box");
  fr.box = bls::Mat3{};
  AUDIT_CHECK("D12", !bls::deriveGridSpec(cfg, fr, nullptr, s, err), "CELL without a cell must fail");
}

}  // namespace

int main() {
  testSha();
  testPdbRoundTrip();
  testParser();
  testGridSpec();
  return audit::finish("test_audit_harness");
}
