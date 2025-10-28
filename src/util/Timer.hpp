#pragma once

#include <chrono>

namespace bls {

class ScopedTimer {
 public:
  ScopedTimer() : start_(Clock::now()) {}

  void reset() { start_ = Clock::now(); }

  double elapsedMilliseconds() const {
    return std::chrono::duration<double, std::milli>(Clock::now() - start_)
        .count();
  }

 private:
  using Clock = std::chrono::high_resolution_clock;
  std::chrono::time_point<Clock> start_;
};

}  // namespace bls

