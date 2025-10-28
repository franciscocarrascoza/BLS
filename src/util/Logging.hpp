#pragma once

#include <iostream>
#include <string>

namespace bls {

enum class LogLevel { Error = 0, Warn = 1, Info = 2, Debug = 3 };

class Logger {
 public:
  static void setLevel(LogLevel lvl) { instance().level_ = lvl; }

  static LogLevel level() { return instance().level_; }

  template <typename... Args>
  static void info(Args&&... args) {
    if (instance().level_ >= LogLevel::Info) {
      instance().write("[INFO] ", std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static void warn(Args&&... args) {
    if (instance().level_ >= LogLevel::Warn) {
      instance().write("[WARN] ", std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static void error(Args&&... args) {
    if (instance().level_ >= LogLevel::Error) {
      instance().write("[ERROR] ", std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static void debug(Args&&... args) {
    if (instance().level_ >= LogLevel::Debug) {
      instance().write("[DEBUG] ", std::forward<Args>(args)...);
    }
  }

 private:
  Logger() = default;

  template <typename... Args>
  void write(const std::string& prefix, Args&&... args) {
    std::cerr << prefix;
    (std::cerr << ... << std::forward<Args>(args)) << '\n';
  }

  static Logger& instance() {
    static Logger inst;
    return inst;
  }

  LogLevel level_{LogLevel::Info};
};

}  // namespace bls

