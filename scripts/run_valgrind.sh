#!/usr/bin/env bash
# Run memcheck under valgrind against a Debug build.
# Usage: BLS/scripts/run_valgrind.sh <program> [args...]
# The Debug build (build-debug/) is created separately, e.g.:
#   BLS_BUILD_SUBDIR=build-debug BLS/scripts/build.sh -DCMAKE_BUILD_TYPE=Debug
set -euo pipefail
exec valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all \
  --track-origins=yes --error-exitcode=97 "$@"
