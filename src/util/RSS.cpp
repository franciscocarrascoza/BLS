#include "util/RSS.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <string>

#if defined(__linux__)
#include <sys/resource.h>
#include <sys/sysinfo.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#endif

namespace bls {

std::size_t currentRSSBytes() {
#if defined(__linux__)
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    return static_cast<std::size_t>(usage.ru_maxrss) * 1024;
  }
  return 0;
#elif defined(__APPLE__)
  mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS) {
    return static_cast<std::size_t>(info.resident_size);
  }
  return 0;
#else
  return 0;
#endif
}

std::size_t availableSystemRAMBytes() {
#if defined(__linux__)
  struct sysinfo info;
  if (sysinfo(&info) == 0) {
    // Return available RAM (free + buffers + cached)
    std::ifstream meminfo("/proc/meminfo");
    if (meminfo.is_open()) {
      std::size_t memAvailable = 0;
      std::string line;
      while (std::getline(meminfo, line)) {
        if (line.compare(0, 13, "MemAvailable:") == 0) {
          std::size_t pos = line.find_first_of("0123456789");
          if (pos != std::string::npos) {
            memAvailable = std::stoull(line.substr(pos)) * 1024; // Convert kB to bytes
            break;
          }
        }
      }
      if (memAvailable > 0) {
        return memAvailable;
      }
    }
    // Fallback to total free RAM
    return static_cast<std::size_t>(info.freeram) * info.mem_unit;
  }
  return 0;
#elif defined(__APPLE__)
  int mib[2] = {CTL_HW, HW_MEMSIZE};
  std::size_t length = sizeof(std::size_t);
  std::size_t memsize = 0;
  if (sysctl(mib, 2, &memsize, &length, nullptr, 0) == 0) {
    return memsize;
  }
  return 0;
#else
  return 0;
#endif
}

std::size_t estimateGridMemoryBytes(std::size_t nx, std::size_t ny, std::size_t nz) {
  std::size_t totalVoxels = nx * ny * nz;
  // Two uint8_t arrays: occupancy and visited
  return totalVoxels * 2;
}

double maxBoxDimensionForRAM(double gridSpacing, std::size_t maxRAMBytes) {
  // Grid memory = 2 * (L/spacing)^3 bytes for a cubic box
  // Solve: 2 * (L/spacing)^3 = maxRAMBytes
  // L = spacing * (maxRAMBytes / 2)^(1/3)
  double maxVoxelsPerDim = std::cbrt(static_cast<double>(maxRAMBytes) / 2.0);
  return gridSpacing * maxVoxelsPerDim;
}

}  // namespace bls

