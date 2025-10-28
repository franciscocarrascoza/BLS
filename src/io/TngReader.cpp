#include "io/TrajectoryReader.hpp"

#include <string>
#include <vector>

#include "util/Logging.hpp"

#ifdef USE_TNG
#include <tng/tng_io.h>
#endif

namespace bls {

#ifdef USE_TNG

class TngReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    if (tng_trajectory_init(&traj_) != TNG_SUCCESS) {
      err = "Failed to initialise TNG trajectory.";
      return false;
    }
    if (tng_trajectory_open_file(traj_, path.c_str(), 'r') != TNG_SUCCESS) {
      err = "Unable to open TNG file: " + path;
      return false;
    }
    int64_t nParticles = 0;
    if (tng_trajectory_num_particles_get(traj_, &nParticles) != TNG_SUCCESS) {
      err = "Unable to query number of particles in TNG.";
      return false;
    }
    natoms_ = static_cast<int>(nParticles);
    coords_.resize(static_cast<std::size_t>(natoms_) * 3);
    path_ = path;
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!traj_) {
      err = "TNG reader not opened.";
      return false;
    }

    tng_function_status status = tng_trajectory_next_frame(traj_);
    if (status == TNG_CRITICAL) {
      err = "Error while reading TNG frame.";
      return false;
    }
    if (status == TNG_FAILURE) {
      return false;
    }

    double time = 0.0;
    tng_trajectory_time_get(traj_, &time);

    float boxData[9] = {0.0f};
    tng_trajectory_box_shape_get(traj_, boxData);

    int64_t nParticles = 0;
    float* positions = nullptr;
    if (tng_trajectory_particle_positions_get(traj_, &positions, &nParticles) != TNG_SUCCESS) {
      err = "Unable to fetch particle positions from TNG.";
      return false;
    }
    if (nParticles != natoms_) {
      err = "Inconsistent particle count in TNG frame.";
      return false;
    }

    frame.xyz.resize(static_cast<std::size_t>(natoms_));
    for (int i = 0; i < natoms_; ++i) {
      frame.xyz[i] =
          Vec3{static_cast<double>(positions[3 * i + 0]),
               static_cast<double>(positions[3 * i + 1]),
               static_cast<double>(positions[3 * i + 2])};
    }
    frame.natoms = natoms_;
    frame.time = time;
    frame.box =
        Mat3{Vec3{boxData[0], boxData[3], boxData[6]}, Vec3{boxData[1], boxData[4], boxData[7]},
             Vec3{boxData[2], boxData[5], boxData[8]}};
    return true;
  }

  void close() override {
    if (traj_) {
      tng_trajectory_destroy(&traj_);
      traj_ = nullptr;
    }
  }

 private:
  tng_trajectory_t traj_{nullptr};
  int natoms_{0};
  std::vector<float> coords_;
  std::string path_;
};

#endif  // USE_TNG

TrajectoryReaderPtr makeTngReader() {
#ifdef USE_TNG
  return TrajectoryReaderPtr(new TngReader());
#else
  return nullptr;
#endif
}

}  // namespace bls

