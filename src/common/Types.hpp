#pragma once

#include <cmath>
#include <ostream>
#include <stdexcept>

namespace bls {

struct Vec3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  Vec3() = default;
  Vec3(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

  Vec3 operator+(const Vec3& other) const {
    return Vec3{x + other.x, y + other.y, z + other.z};
  }

  Vec3 operator-(const Vec3& other) const {
    return Vec3{x - other.x, y - other.y, z - other.z};
  }

  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Vec3& operator-=(const Vec3& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  Vec3 operator*(double s) const { return Vec3{x * s, y * s, z * s}; }

  Vec3 operator/(double s) const { return Vec3{x / s, y / s, z / s}; }

  Vec3& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  Vec3& operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }
  double& operator[](int i) {
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
  }

  const double& operator[](int i) const {
    if (i == 0) return x;
    if (i == 1) return y;
    return z;
  }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }

inline double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return Vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x};
}

inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

struct Mat3 {
  Vec3 cols[3];

  Mat3() = default;
  Mat3(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
    cols[0] = c0;
    cols[1] = c1;
    cols[2] = c2;
  }

  Vec3 column(int i) const { return cols[i]; }
  void setColumn(int i, const Vec3& v) { cols[i] = v; }

  Vec3 row(int i) const { return Vec3{cols[0][i], cols[1][i], cols[2][i]}; }

  Vec3& operator[](int i) { return cols[i]; }
  const Vec3& operator[](int i) const { return cols[i]; }

  Vec3 operator*(const Vec3& v) const {
    return cols[0] * v.x + cols[1] * v.y + cols[2] * v.z;
  }

  Mat3 operator*(double s) const {
    return Mat3{cols[0] * s, cols[1] * s, cols[2] * s};
  }

  Mat3 transpose() const {
    return Mat3{Vec3{cols[0].x, cols[1].x, cols[2].x},
                Vec3{cols[0].y, cols[1].y, cols[2].y},
                Vec3{cols[0].z, cols[1].z, cols[2].z}};
  }
};

inline Mat3 outer(const Vec3& a, const Vec3& b) {
  return Mat3{Vec3{a.x * b.x, a.x * b.y, a.x * b.z},
              Vec3{a.y * b.x, a.y * b.y, a.y * b.z},
              Vec3{a.z * b.x, a.z * b.y, a.z * b.z}};
}

inline double determinant(const Mat3& m) {
  return dot(m.cols[0], cross(m.cols[1], m.cols[2]));
}

inline Mat3 inverse(const Mat3& m) {
  Vec3 c0 = cross(m.cols[1], m.cols[2]);
  Vec3 c1 = cross(m.cols[2], m.cols[0]);
  Vec3 c2 = cross(m.cols[0], m.cols[1]);
  double det = dot(m.cols[0], c0);
  if (std::abs(det) < 1e-12) {
    throw std::runtime_error("Singular matrix inversion request");
  }
  return Mat3{c0 / det, c1 / det, c2 / det};
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
  os << "(" << v.x << "," << v.y << "," << v.z << ")";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const Mat3& m) {
  os << "[" << m.cols[0] << "," << m.cols[1] << "," << m.cols[2] << "]";
  return os;
}

}  // namespace bls
