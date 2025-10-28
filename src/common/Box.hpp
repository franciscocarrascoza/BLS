#pragma once

#include <cmath>

#include "common/Types.hpp"

namespace bls {

inline Mat3 buildTriclinicBox(double a, double b, double c, double alpha_deg, double beta_deg,
                              double gamma_deg) {
  constexpr double pi = 3.14159265358979323846;
  const double deg2rad = pi / 180.0;
  const double alpha = alpha_deg * deg2rad;
  const double beta = beta_deg * deg2rad;
  const double gamma = gamma_deg * deg2rad;

  const double cx = c * std::cos(beta);
  const double cy =
      c * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
  const double cz_sq = c * c - cx * cx - cy * cy;
  const double cz = cz_sq > 0.0 ? std::sqrt(cz_sq) : 0.0;

  Vec3 a1{a, 0.0, 0.0};
  Vec3 a2{b * std::cos(gamma), b * std::sin(gamma), 0.0};
  Vec3 a3{cx, cy, cz};
  return Mat3{a1, a2, a3};
}

}  // namespace bls
