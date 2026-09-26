// Created: 2026-09-26T11:20+02:00 | by: CC audit (manager) | purpose: shared helpers for the pre-release audit test suite (tests/audit/test_audit_*.cpp)
//
// Every audit test names the claims-ledger ID it guards (audit/CLAIMS_LEDGER_*.md)
// through AUDIT_CHECK(claimId, condition, message). A test binary exits non-zero
// iff any check failed, and prints one PASS/FAIL line per claim ID at the end.
//
// Partitions are compared CANONICALLY, never by count alone: every occupied voxel
// is relabelled with the smallest linear index in its component ("min-voxel
// relabelling"), so two labellings describe the same partition iff their canonical
// arrays are byte-identical. Identical count + max size is not evidence of an
// identical partition.

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"
#include "common/Types.hpp"
#include "lattice/Basis.hpp"
#include "lattice/Enumerator.hpp"
#include "refine/SkipDFS.hpp"

namespace audit {

// ---------------------------------------------------------------------------
// Claim-tagged checks
// ---------------------------------------------------------------------------
struct ClaimTally {
  int pass{0};
  int fail{0};
};

inline std::map<std::string, ClaimTally>& claimTallies() {
  static std::map<std::string, ClaimTally> t;
  return t;
}

inline int& failureBudget() {  // print at most this many individual failure lines per claim
  static int b = 20;
  return b;
}

inline bool check(const std::string& claim, bool ok, const std::string& msg) {
  auto& t = claimTallies()[claim];
  if (ok) {
    ++t.pass;
  } else {
    ++t.fail;
    if (t.fail <= failureBudget()) std::printf("  FAIL [%s] %s\n", claim.c_str(), msg.c_str());
  }
  return ok;
}

#define AUDIT_CHECK(claim, cond, msg) ::audit::check((claim), static_cast<bool>(cond), (msg))

// Prints the per-claim summary and returns the process exit code.
inline int finish(const char* testName) {
  int failedClaims = 0;
  std::printf("\n== %s: per-claim summary ==\n", testName);
  for (const auto& kv : claimTallies()) {
    const bool ok = kv.second.fail == 0;
    if (!ok) ++failedClaims;
    std::printf("  %-4s %-28s checks=%d failed=%d\n", ok ? "PASS" : "FAIL", kv.first.c_str(),
                kv.second.pass + kv.second.fail, kv.second.fail);
  }
  std::printf("== %s: %s (%d claim(s) failing) ==\n", testName, failedClaims ? "FAILED" : "PASSED",
              failedClaims);
  return failedClaims ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Grids
// ---------------------------------------------------------------------------
struct Grid3 {
  int nx{0}, ny{0}, nz{0};
  std::vector<uint8_t> occ;
  std::size_t size() const { return static_cast<std::size_t>(nx) * ny * nz; }
  std::size_t index(int i, int j, int k) const {  // x-major, z fastest -- same as the library
    return static_cast<std::size_t>(i) * ny * nz + static_cast<std::size_t>(j) * nz + k;
  }
  bool at(int i, int j, int k) const { return occ[index(i, j, k)] != 0; }
  void set(int i, int j, int k, bool v = true) { occ[index(i, j, k)] = v ? 1 : 0; }
};

inline Grid3 emptyGrid(int nx, int ny, int nz) {
  Grid3 g;
  g.nx = nx; g.ny = ny; g.nz = nz;
  g.occ.assign(g.size(), 0);
  return g;
}

// Bernoulli(p) occupancy. Deterministic in `rng`.
inline Grid3 randomGrid(std::mt19937_64& rng, int nx, int ny, int nz, double p) {
  Grid3 g = emptyGrid(nx, ny, nz);
  std::bernoulli_distribution b(p);
  for (auto& v : g.occ) v = b(rng) ? 1 : 0;
  return g;
}

// ---------------------------------------------------------------------------
// Canonical partitions
// ---------------------------------------------------------------------------
// Input: labels with -1 for empty voxels and any non-negative id per component.
// Output: -1 for empty; for occupied voxels, the smallest linear voxel index of the
// voxel's component. Throws nothing; asserts label/occupancy consistency via the
// returned flag.
inline std::vector<long long> canonical(const std::vector<int>& labels) {
  std::map<int, long long> minIdx;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    if (labels[i] < 0) continue;
    auto it = minIdx.find(labels[i]);
    if (it == minIdx.end()) minIdx.emplace(labels[i], static_cast<long long>(i));
  }
  std::vector<long long> out(labels.size(), -1);
  for (std::size_t i = 0; i < labels.size(); ++i)
    if (labels[i] >= 0) out[i] = minIdx[labels[i]];
  return out;
}

// True iff every occupied voxel carries a label >= 0 and every empty voxel -1.
inline bool labelsCoverOccupancy(const Grid3& g, const std::vector<int>& labels) {
  if (labels.size() != g.size()) return false;
  for (std::size_t i = 0; i < g.size(); ++i)
    if ((g.occ[i] != 0) != (labels[i] >= 0)) return false;
  return true;
}

// Sorted (descending) component sizes from a labelling.
inline std::vector<int> sizesOf(const std::vector<int>& labels) {
  std::map<int, int> n;
  for (int l : labels) if (l >= 0) ++n[l];
  std::vector<int> s;
  for (const auto& kv : n) s.push_back(kv.second);
  std::sort(s.begin(), s.end(), std::greater<int>());
  return s;
}

// Subset relation for a PARTIAL labelling (BLS labels only the components it
// reaches; unreached occupied voxels stay -1): every component of `partial`
// must coincide exactly with one component of `full` (no truncation, split or
// merge). Returns an empty string on success, otherwise a description.
inline std::string subsetViolation(const std::vector<int>& partial, const std::vector<int>& full) {
  if (partial.size() != full.size()) return "size mismatch";
  std::map<int, int> partialToFull;          // partial label -> full label
  std::map<int, int> fullToPartial;          // full label -> partial label (if reached)
  std::map<int, long long> partialCount, fullCount;
  for (std::size_t i = 0; i < full.size(); ++i) {
    if (full[i] >= 0) ++fullCount[full[i]];
    if (partial[i] < 0) continue;
    if (full[i] < 0) return "partial labels an empty voxel";
    ++partialCount[partial[i]];
    auto a = partialToFull.emplace(partial[i], full[i]);
    if (a.first->second != full[i]) return "a partial component spans two full components (merge)";
    auto b = fullToPartial.emplace(full[i], partial[i]);
    if (b.first->second != partial[i]) return "a full component is split across partial components";
  }
  for (const auto& kv : partialToFull)
    if (partialCount[kv.first] != fullCount[kv.second]) return "a partial component is a truncated full component";
  return {};
}

// ---------------------------------------------------------------------------
// Running the library
// ---------------------------------------------------------------------------
inline bls::ClusterParams paramsFor(const Grid3& g) {
  bls::ClusterParams p;
  p.nx = g.nx; p.ny = g.ny; p.nz = g.nz;
  p.connectivity = 6;
  return p;
}

// Exact methods with label output (TraditionalDFS, GCBD, CC3D, CC3DOptimized,
// RLECCL, RLECCLOptimized). `visited` is supplied zeroed, as the harness does.
inline bls::ClusterResult runLabelled(bls::ClusterAlgorithm a, const Grid3& g, std::vector<int>& labels,
                                      bls::ClusterParams p) {
  std::vector<uint8_t> visited(g.size(), 0);
  return bls::runClusterAlgorithm(a, p, g.occ, visited, &labels);
}
inline bls::ClusterResult runLabelled(bls::ClusterAlgorithm a, const Grid3& g, std::vector<int>& labels) {
  return runLabelled(a, g, labels, paramsFor(g));
}

// Any method, no labels (count + sizes only).
inline bls::ClusterResult runUnlabelled(bls::ClusterAlgorithm a, const Grid3& g, bls::ClusterParams p) {
  std::vector<uint8_t> visited(g.size(), 0);
  bls::ClusterResult r = bls::runClusterAlgorithm(a, p, g.occ, visited, nullptr);
  std::sort(r.clusterSizes.begin(), r.clusterSizes.end(), std::greater<int>());
  return r;
}

// BLS on a raw occupancy grid, reproducing Analyzer::processFrame's seeding and
// refinement exactly (bls/BLS.cpp:227-267): Enumerator built from the lattice
// scaled so that its nearest-neighbour distance is `dnnVoxels`, occupancy-driven;
// SkipDFS refinement with CONNECTIVITY 6 and stride `stride`; labels dense in
// seed order, -1 for unreached voxels.
struct BlsRaw {
  int seeds{0};
  int nclusters{0};
  std::vector<int> labels;  // -1 = empty or not reached
  std::vector<int> sizes;   // descending
};

inline bls::LatticeDescriptor latticeFor(bls::LatticeType t, bls::CenteringType c) {
  bls::LatticeSettings s;
  s.lattice = t;
  s.centering = c;
  s.latticeSet = s.centeringSet = true;
  return bls::buildLattice(s);  // triclinic default = reference basis a=1,b=1.2,c=1.4, 90/100/110
}

inline BlsRaw runBlsRaw(const Grid3& g, const bls::LatticeDescriptor& d, double dnnVoxels, int stride = 1) {
  BlsRaw o;
  const bls::Mat3 sb = d.basis * (dnnVoxels / d.dmin);  // as BLS.cpp Impl
  bls::Enumerator en(sb, d.offsets, g.nx, g.ny, g.nz, g.occ);
  o.seeds = en.count();
  std::vector<uint8_t> visited(g.size(), 0);
  o.labels.assign(g.size(), -1);
  bls::SkipDFSConfig cfg{g.nx, g.ny, g.nz, 6, stride};
  bls::SkipDFS dfs(cfg, g.occ, visited);
  en.forEach([&](const bls::Enumerator::Seed& s) {
    const std::size_t idx = g.index(s.x, s.y, s.z);
    if (!g.occ[idx] || visited[idx]) return;
    int n = dfs.runFrom(s.x, s.y, s.z, &o.labels, o.nclusters);
    if (n > 0) { ++o.nclusters; o.sizes.push_back(n); }
  });
  std::sort(o.sizes.begin(), o.sizes.end(), std::greater<int>());
  return o;
}

}  // namespace audit
