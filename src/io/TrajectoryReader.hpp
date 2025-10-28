#pragma once

#include <memory>
#include <string>
#include <vector>

#include "common/Types.hpp"

namespace bls {

struct Frame {
  std::vector<Vec3> xyz;
  Mat3 box{};
  int natoms{0};
  double time{0.0};
};

class TrajectoryReader {
 public:
  virtual ~TrajectoryReader() = default;
  virtual bool open(const std::string& path, std::string& err) = 0;
  virtual bool read(Frame& f, std::string& err) = 0;
  virtual void close() = 0;
};

using TrajectoryReaderPtr = std::unique_ptr<TrajectoryReader>;

TrajectoryReaderPtr makeTrajectoryReader(const std::string& path, const std::string& format,
                                         std::string& err);

}  // namespace bls

