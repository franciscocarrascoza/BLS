#!/usr/bin/env bash
# Guarded build entry point.
#
# WHY THIS EXISTS: a CMakeLists.txt-embedded guard cannot, by itself, catch
# the sibling-tree trap. If build/ is a copy of another checkout's build
# directory, the generated Makefile's internal rerun rule cd's into the
# *source* directory recorded when that build/ was originally generated --
# a real, different, still-existing path on disk -- before ever invoking
# cmake. That reruns cmake against the OTHER tree's CMakeLists.txt, not this
# one, so a guard living only in this file's CMakeLists.txt never executes.
# (Verified empirically: a build/ directory generated from a second checkout
# and dropped into this tree's build/, then rebuilt with plain `make`, prints
# "Build files have been written to: <the other tree>/build" and exits 0.)
#
# This script is the check that actually fires in that scenario, because it
# runs before cmake/make ever gets a chance to redirect anywhere: it inspects
# build/CMakeCache.txt textually and compares its own realpath -- computed
# independently, from this script's own location, not from anything read out
# of build/ -- to what's recorded there.
#
# Usage: BLS/scripts/build.sh [extra cmake args...]
# Env: BLS_BUILD_SUBDIR=build-asan BLS/scripts/build.sh -DBLS_SANITIZE=ON ...
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
SOURCE_DIR="$(cd -- "${SCRIPT_DIR}/.." >/dev/null 2>&1 && pwd)"
BUILD_DIR="${SOURCE_DIR}/${BLS_BUILD_SUBDIR:-build}"

if [[ -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
  cached_home="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "${BUILD_DIR}/CMakeCache.txt")"
  if [[ -n "${cached_home}" ]]; then
    cached_home_real="$(cd -- "${cached_home}" >/dev/null 2>&1 && pwd || echo "${cached_home}")"
    if [[ "${cached_home_real}" != "${SOURCE_DIR}" ]]; then
      echo "error: ${BUILD_DIR} was generated from a different source tree." >&2
      echo "  build/CMakeCache.txt CMAKE_HOME_DIRECTORY = ${cached_home_real}" >&2
      echo "  this checkout                             = ${SOURCE_DIR}" >&2
      echo "This is the sibling-tree trap: an incremental 'make' here would silently" >&2
      echo "reconfigure and rebuild against ${cached_home_real} instead, leaving this" >&2
      echo "checkout's binaries untouched with no error. Fix:" >&2
      echo "  rm -rf ${BUILD_DIR} && ${SCRIPT_DIR}/build.sh" >&2
      exit 1
    fi
  fi
fi

mkdir -p "${BUILD_DIR}"
cmake -S "${SOURCE_DIR}" -B "${BUILD_DIR}" "$@"
cmake --build "${BUILD_DIR}" -j"$(nproc)"
