#include "util/RSS.hpp"

#include <cstddef>

#if defined(__linux__)
#include <sys/resource.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
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

}  // namespace bls

