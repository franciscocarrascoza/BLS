#pragma once

#include <string>

#include "bls/Options.hpp"

namespace bls {

class Parser {
 public:
  bool parseFile(const std::string& path, BLSConfig& config, std::string& err);
};

}  // namespace bls

