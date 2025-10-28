#include "io/TrajectoryReader.hpp"

#include <string>
#include <vector>

#include "util/Logging.hpp"

#ifdef USE_XDRFILE
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#endif

namespace bls {

#ifdef USE_XDRFILE

class XtcReader final : public TrajectoryReader {
 public:
  bool open(const std::string& path, std::string& err) override {
    close();
    int natoms = 0;
    if (read_xtc_natoms(const_cast<char*>(path.c_str()), &natoms) != exdrOK || natoms <= 0) {
      err = "Failed to read number of atoms from XTC: " + path;
      return false;
    }
    handle_ = xdrfile_open(path.c_str(), "r");
    if (!handle_) {
      err = "Unable to open XTC file: " + path;
      return false;
    }
    natoms_ = natoms;
    coords_.resize(static_cast<std::size_t>(natoms_) * 3);
    return true;
  }

  bool read(Frame& frame, std::string& err) override {
    if (!handle_) {
      err = "XTC reader not opened.";
      return false;
    }

    matrix box;
    int step = 0;
    float time = 0.0f;
    float prec = 1000.0f;
    auto status = read_xtc(handle_, natoms_, &step, &time, box,
                           reinterpret_cast<rvec*>(coords_.data()), &prec);
    if (status == exdrENDOFFILE) {
      return false;
    }
    if (status != exdrOK) {
      err = "Error while reading XTC frame.";
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

TrajectoryReaderPtr makeXtcReader() {
#ifdef USE_XDRFILE
  return TrajectoryReaderPtr(new XtcReader());
#else
  return nullptr;
#endif
}

}  // namespace bls
