#include "io/TrajectoryReader.hpp"

#include <string>
#include <vector>

#include "util/Logging.hpp"

#ifdef USE_XDRFILE
#include <xdrfile.h>
#include <xdrfile_trr.h>
#endif

namespace bls {

#ifdef USE_XDRFILE

class TrrReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    int natoms = 0;
    if (read_trr_natoms(const_cast<char*>(path.c_str()), &natoms) != exdrOK ||
        natoms <= 0) {
      err = "Failed to read number of atoms from TRR: " + path;
      return false;
    }
    handle_ = xdrfile_open(path.c_str(), "r");
    if (!handle_) {
      err = "Unable to open TRR file: " + path;
      return false;
    }
    natoms_ = natoms;
    coords_.resize(static_cast<std::size_t>(natoms_) * 3);
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!handle_) {
      err = "TRR reader not opened.";
      return false;
    }

    matrix box;
    int step = 0;
    float time = 0.0f;
    float lambda = 0.0f;
    int hasProps = 0;
    int status = read_trr(handle_, natoms_, &step, &time, &lambda, box,
                          reinterpret_cast<rvec*>(coords_.data()), nullptr,
                          nullptr, &hasProps);

    if (status == exdrENDOFFILE) {
      return false;
    }
    if (status != exdrOK) {
      err = "Error while reading TRR frame.";
      return false;
    }

    frame.xyz.resize(static_cast<std::size_t>(natoms_));
    for (int i = 0; i < natoms_; ++i) {
      frame.xyz[i] = Vec3{coords_[3 * i + 0], coords_[3 * i + 1],
                          coords_[3 * i + 2]};
    }
    frame.natoms = natoms_;
    frame.time = static_cast<double>(time);
    frame.box = Mat3{Vec3{box[0][0], box[0][1], box[0][2]},
                     Vec3{box[1][0], box[1][1], box[1][2]},
                     Vec3{box[2][0], box[2][1], box[2][2]}};
    return true;
  }

  void close() override {
    if (handle_) {
      xdrfile_close(handle_);
      handle_ = nullptr;
    }
  }

 private:
  XDRFILE* handle_{nullptr};
  int natoms_{0};
  std::vector<float> coords_;
};

#endif  // USE_XDRFILE

TrajectoryReaderPtr makeTrrReader() {
#ifdef USE_XDRFILE
  return TrajectoryReaderPtr(new TrrReader());
#else
  return nullptr;
#endif
}

}  // namespace bls
