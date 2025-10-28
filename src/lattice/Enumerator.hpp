#pragma once

#include <array>
#include <functional>
#include <vector>

#include "common/Types.hpp"

namespace bls {

class Enumerator {
 public:
  Enumerator(const Mat3& basis, const std::vector<Vec3>& offsets, int nx, int ny, int nz);

  struct Seed {
    int x;
    int y;
    int z;
  };

  template <typename F>
  void forEach(F&& func) const {
    for (const auto& seed : seeds_) {
      func(seed);
    }
  }

  int count() const { return static_cast<int>(seeds_.size()); }

 private:
  std::vector<Seed> seeds_;
};

}  // namespace bls

