#include "io/TrajectoryReader.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <string>

#include "util/Logging.hpp"

namespace bls {

namespace {

std::string toLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

std::string extensionOf(const std::string& path) {
  auto pos = path.find_last_of('.');
  if (pos == std::string::npos) return {};
  return toLower(path.substr(pos + 1));
}

}  // namespace

TrajectoryReaderPtr makeXtcReader();
TrajectoryReaderPtr makeTrrReader();
TrajectoryReaderPtr makeTngReader();
TrajectoryReaderPtr makeGroReader();
TrajectoryReaderPtr makePdbReader();
TrajectoryReaderPtr makeXyzReader();
TrajectoryReaderPtr makeMolReader();

TrajectoryReaderPtr makeTrajectoryReader(const std::string& path, const std::string& format,
                                         std::string& err) {
  std::string chosen = toLower(format);
  if (chosen == "auto") {
    chosen = extensionOf(path);
  }

  TrajectoryReaderPtr reader;
  if (chosen == "xtc") {
    reader = makeXtcReader();
    if (!reader) err = "This build lacks XTC support (recompile with USE_XDRFILE).";
  } else if (chosen == "trr") {
    reader = makeTrrReader();
    if (!reader) err = "This build lacks TRR support (recompile with USE_XDRFILE).";
  } else if (chosen == "tng") {
    reader = makeTngReader();
    if (!reader) err = "This build lacks TNG support (recompile with USE_TNG).";
  } else if (chosen == "gro") {
    reader = makeGroReader();
  } else if (chosen == "pdb") {
    reader = makePdbReader();
  } else if (chosen == "xyz") {
    reader = makeXyzReader();
  } else if (chosen == "mol" || chosen == "sdf") {
    reader = makeMolReader();
  } else {
    err = "Unsupported molecular system format: " + chosen;
    return nullptr;
  }

  if (!reader) return nullptr;

  if (!reader->open(path, err)) {
    return nullptr;
  }

  return reader;
}

}  // namespace bls

