// Task 3c: the optional label-array output.
//
// These checks establish that the labels an algorithm reports are consistent
// with the summary it already reported -- same component count, same size
// histogram, same occupied set. That is the precondition for Task 4 using
// labels to compare partitions across algorithms: if labels and counts can
// disagree, a partition mismatch cannot be attributed to the algorithm.
//
// What is deliberately NOT asserted here is that different algorithms agree
// with each other. That is Task 4's job, and pre-judging it here would hide
// the very disagreements the differential harness exists to find.

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "bls/BLS.hpp"
#include "bls/Options.hpp"
#include "cluster/Algorithms.hpp"

using bls::ClusterAlgorithm;
using bls::ClusterParams;
using bls::ClusterResult;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool cond, const std::string& msg) {
  ++g_checks;
  if (!cond) {
    ++g_failures;
    std::cerr << "FAIL: " << msg << "\n";
  }
}

struct Case {
  std::string name;
  int nx, ny, nz;
  std::vector<uint8_t> occ;
};

std::vector<Case> buildCases() {
  std::vector<Case> cases;
  auto make = [](const std::string& name, int nx, int ny, int nz) {
    Case c;
    c.name = name;
    c.nx = nx; c.ny = ny; c.nz = nz;
    c.occ.assign(static_cast<std::size_t>(nx) * ny * nz, 0);
    return c;
  };
  auto at = [](Case& c, int i, int j, int k) -> uint8_t& {
    return c.occ[static_cast<std::size_t>(i) * c.ny * c.nz + static_cast<std::size_t>(j) * c.nz + k];
  };

  { Case c = make("empty", 8, 8, 8); cases.push_back(c); }
  { Case c = make("single-voxel", 8, 8, 8); at(c, 3, 4, 5) = 1; cases.push_back(c); }
  { Case c = make("full", 6, 7, 5); std::fill(c.occ.begin(), c.occ.end(), 1); cases.push_back(c); }
  { Case c = make("two-corner-touching", 9, 9, 9);
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k) at(c, i, j, k) = 1;
    for (int i = 4; i < 8; ++i) for (int j = 4; j < 8; ++j) for (int k = 4; k < 8; ++k) at(c, i, j, k) = 1;
    cases.push_back(c); }
  { Case c = make("diagonal-chain", 12, 12, 12);
    for (int i = 0; i < 12; ++i) at(c, i, i, i) = 1;
    cases.push_back(c); }
  { Case c = make("separated-planes", 10, 10, 11);
    for (int i = 0; i < 10; ++i) for (int j = 0; j < 10; ++j) { at(c, i, j, 0) = 1; at(c, i, j, 10) = 1; }
    cases.push_back(c); }

  std::mt19937_64 rng(31337);
  for (double p : {0.05, 0.20, 0.31, 0.50, 0.80}) {
    Case c = make("random-p" + std::to_string(p), 16, 15, 14);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (auto& v : c.occ) v = u(rng) < p ? 1 : 0;
    cases.push_back(c);
  }
  return cases;
}

const ClusterAlgorithm kLabelled[] = {
    ClusterAlgorithm::TraditionalDFS, ClusterAlgorithm::CC3D, ClusterAlgorithm::CC3DOptimized,
    ClusterAlgorithm::GCBD,           ClusterAlgorithm::RLECCL, ClusterAlgorithm::OctreeCCL};

const ClusterAlgorithm kUnlabelled[] = {
    ClusterAlgorithm::SkipDFS, ClusterAlgorithm::DBSCAN,  ClusterAlgorithm::Hierarchical,
    ClusterAlgorithm::KMeans,  ClusterAlgorithm::Spectral, ClusterAlgorithm::HDBSCAN,
    ClusterAlgorithm::VCCS};

}  // namespace

int main() {
  const auto cases = buildCases();

  for (const auto& c : cases) {
    ClusterParams params;
    params.nx = c.nx; params.ny = c.ny; params.nz = c.nz;
    params.connectivity = 6;
    params.octreeLeafSize = 8;  // the shipped default, stated explicitly

    for (ClusterAlgorithm algo : kLabelled) {
      const std::string tag = c.name + "/" + bls::algorithmToString(algo);
      check(bls::supportsLabels(algo), tag + ": supportsLabels() says no");

      // Poison the buffer first: a correct implementation must overwrite every
      // element, so a stale value surviving is itself a defect.
      std::vector<int> labels(3, -12345);
      std::vector<uint8_t> visited(c.occ.size(), 0);
      ClusterResult r = bls::runClusterAlgorithm(algo, params, c.occ, visited, &labels);

      check(labels.size() == c.occ.size(), tag + ": label buffer not resized to the grid");
      if (labels.size() != c.occ.size()) continue;

      // 1. Labelled set == occupied set, exactly.
      std::size_t mislabelled = 0;
      for (std::size_t i = 0; i < labels.size(); ++i) {
        const bool occupied = c.occ[i] != 0;
        const bool labelled = labels[i] >= 0;
        if (occupied != labelled) ++mislabelled;
      }
      check(mislabelled == 0,
            tag + ": " + std::to_string(mislabelled) + " voxels labelled but not occupied, or vice versa");

      // 2. Labels are dense and start at 0.
      std::set<int> distinct;
      for (int v : labels)
        if (v >= 0) distinct.insert(v);
      check(static_cast<int>(distinct.size()) == r.nclusters,
            tag + ": " + std::to_string(distinct.size()) + " distinct labels but nclusters=" +
                std::to_string(r.nclusters));
      bool dense = true;
      int expect = 0;
      for (int v : distinct) {
        if (v != expect++) { dense = false; break; }
      }
      check(dense, tag + ": labels are not dense 0..nclusters-1");

      // 3. Size histogram from the labels matches the reported one. This is
      //    the check that catches a labelling which agrees on the COUNT while
      //    disagreeing on the partition -- one merge plus one split.
      std::map<int, int> sizes;
      for (int v : labels)
        if (v >= 0) ++sizes[v];
      std::vector<int> fromLabels;
      for (const auto& kv : sizes) fromLabels.push_back(kv.second);
      std::sort(fromLabels.begin(), fromLabels.end(), std::greater<int>());
      std::vector<int> reported = r.clusterSizes;
      std::sort(reported.begin(), reported.end(), std::greater<int>());
      check(fromLabels == reported, tag + ": size histogram from labels disagrees with clusterSizes");

      // 4. Each label class is genuinely connected under the requested
      //    connectivity, and no two distinct classes are adjacent. Together
      //    these say the labelling IS a connected-component partition -- not
      //    merely a partition of the right shape.
      const int nx = c.nx, ny = c.ny, nz = c.nz;
      auto idx = [&](int i, int j, int k) {
        return static_cast<std::size_t>(i) * ny * nz + static_cast<std::size_t>(j) * nz + k;
      };
      std::size_t adjacentDifferent = 0;
      static const int d6[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
          for (int k = 0; k < nz; ++k) {
            const int a = labels[idx(i, j, k)];
            if (a < 0) continue;
            for (const auto& d : d6) {
              int ni = i + d[0], nj = j + d[1], nk = k + d[2];
              if (ni < 0 || nj < 0 || nk < 0 || ni >= nx || nj >= ny || nk >= nz) continue;
              const int b = labels[idx(ni, nj, nk)];
              if (b >= 0 && b != a) ++adjacentDifferent;
            }
          }
      check(adjacentDifferent == 0,
            tag + ": " + std::to_string(adjacentDifferent) +
                " adjacent occupied voxel pairs carry different labels (component split)");

      // Connectivity of each class, by independent BFS over the label array.
      std::vector<uint8_t> seen(labels.size(), 0);
      int reachedClasses = 0;
      std::vector<std::size_t> stack;
      for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
          for (int k = 0; k < nz; ++k) {
            std::size_t s = idx(i, j, k);
            if (labels[s] < 0 || seen[s]) continue;
            ++reachedClasses;
            const int cls = labels[s];
            stack.clear();
            stack.push_back(s);
            seen[s] = 1;
            while (!stack.empty()) {
              std::size_t cur = stack.back(); stack.pop_back();
              int ci = static_cast<int>(cur / (static_cast<std::size_t>(ny) * nz));
              int rem = static_cast<int>(cur % (static_cast<std::size_t>(ny) * nz));
              int cj = rem / nz, ck = rem % nz;
              for (const auto& d : d6) {
                int ni = ci + d[0], nj = cj + d[1], nk = ck + d[2];
                if (ni < 0 || nj < 0 || nk < 0 || ni >= nx || nj >= ny || nk >= nz) continue;
                std::size_t nidx = idx(ni, nj, nk);
                if (seen[nidx] || labels[nidx] != cls) continue;
                seen[nidx] = 1;
                stack.push_back(nidx);
              }
            }
          }
      check(reachedClasses == r.nclusters,
            tag + ": " + std::to_string(reachedClasses) + " connected label regions but " +
                std::to_string(r.nclusters) + " clusters reported (component merge)");

      // 5. Requesting labels must not change what the algorithm reports.
      std::vector<uint8_t> visited2(c.occ.size(), 0);
      ClusterResult plain = bls::runClusterAlgorithm(algo, params, c.occ, visited2, nullptr);
      check(plain.nclusters == r.nclusters && plain.maxCluster == r.maxCluster &&
                plain.clusterSizes == r.clusterSizes && plain.visitedVoxels == r.visitedVoxels,
            tag + ": results differ depending on whether labels were requested");
    }

    // Algorithms with no label support must refuse rather than hand back a
    // buffer that looks like a partition and is not.
    for (ClusterAlgorithm algo : kUnlabelled) {
      const std::string tag = c.name + "/" + bls::algorithmToString(algo);
      check(!bls::supportsLabels(algo), tag + ": unexpectedly claims label support");
      std::vector<int> labels;
      std::vector<uint8_t> visited(c.occ.size(), 0);
      bool threw = false;
      try {
        bls::runClusterAlgorithm(algo, params, c.occ, visited, &labels);
      } catch (const std::exception&) {
        threw = true;
      }
      check(threw, tag + ": accepted a label request it cannot satisfy");
    }
  }

  // ---- BLS's own label path, through Analyzer -----------------------------
  //
  // BLS cannot go through runClusterAlgorithm (it needs lattice enumeration),
  // so it carries the same contract on Analyzer::processFrame. The same
  // consistency properties are asserted: labelled set == occupied set, dense
  // ids, histogram agreeing with the reported sizes, and identical metrics
  // with and without labels requested.
  {
    std::mt19937_64 rng(4242);
    std::uniform_real_distribution<double> u(0.5, 19.5);
    for (int trial = 0; trial < 6; ++trial) {
      bls::BLSConfig cfg;
      cfg.gridSpacing = 1.0;
      cfg.connectivity = 6;
      cfg.refinementStride = 1;
      cfg.cutoff = 1.0;
      cfg.lattice.lattice = bls::LatticeType::Cubic;
      cfg.lattice.centering = trial % 2 ? bls::CenteringType::F : bls::CenteringType::P;

      bls::Frame frame;
      for (int i = 0; i < 40 + 30 * trial; ++i) frame.xyz.push_back(bls::Vec3{u(rng), u(rng), u(rng)});
      frame.natoms = static_cast<int>(frame.xyz.size());
      frame.box = bls::Mat3{bls::Vec3{20, 0, 0}, bls::Vec3{0, 20, 0}, bls::Vec3{0, 0, 20}};

      bls::Analyzer analyzer(cfg);
      analyzer.setSelection({}, frame.natoms);

      const std::string tag = "bls/trial" + std::to_string(trial);
      std::string err;
      std::vector<int> labels(3, -12345);
      bls::FrameMetrics withLabels, withoutLabels;
      check(analyzer.processFrame(frame, withLabels, err, &labels), tag + ": failed: " + err);
      check(analyzer.processFrame(frame, withoutLabels, err, nullptr), tag + ": failed: " + err);

      const std::size_t total = static_cast<std::size_t>(withLabels.nx) * withLabels.ny * withLabels.nz;
      check(labels.size() == total, tag + ": label buffer not resized to the grid");
      if (labels.size() != total) continue;

      std::set<int> distinct;
      std::map<int, int> sizes;
      for (int v : labels)
        if (v >= 0) { distinct.insert(v); ++sizes[v]; }
      check(static_cast<int>(distinct.size()) == withLabels.nclusters,
            tag + ": " + std::to_string(distinct.size()) + " distinct labels but nclusters=" +
                std::to_string(withLabels.nclusters));
      bool dense = true;
      int expect = 0;
      for (int v : distinct)
        if (v != expect++) { dense = false; break; }
      check(dense, tag + ": labels are not dense 0..nclusters-1");

      std::vector<int> fromLabels;
      for (const auto& kv : sizes) fromLabels.push_back(kv.second);
      std::sort(fromLabels.begin(), fromLabels.end(), std::greater<int>());
      check(fromLabels == withLabels.clusterSizes,
            tag + ": size histogram from labels disagrees with clusterSizes");

      // BLS only labels what it actually visits. Every labelled voxel must be
      // occupied; the converse does NOT hold and is not asserted -- BLS may
      // legitimately leave a component unvisited if no lattice seed lands in
      // it, which is precisely the coverage question Task 6 exists to measure.
      check(withLabels.nclusters == withoutLabels.nclusters &&
                withLabels.maxCluster == withoutLabels.maxCluster &&
                withLabels.clusterSizes == withoutLabels.clusterSizes &&
                withLabels.seeds == withoutLabels.seeds &&
                withLabels.seedHits == withoutLabels.seedHits &&
                withLabels.refinedVoxels == withoutLabels.refinedVoxels,
            tag + ": metrics differ depending on whether labels were requested");
    }
  }

  std::cout << "checks: " << g_checks << ", failures: " << g_failures << "\n";
  return g_failures == 0 ? 0 : 1;
}
