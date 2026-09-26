#!/usr/bin/env bash
# Created: 2026-09-26T11:17+02:00 | by: CC audit | purpose: ONE command that builds and runs the whole test suite (existing + audit) in Debug+ASan/UBSan+_GLIBCXX_ASSERTIONS and in Release.
#
# Usage:  BLS/tests/audit/run_all.sh [extra ctest args, e.g. -L audit | -R test_audit_exact]
# Logs:   $BLS_TEST_LOGDIR (default BLS/build-audit-logs/<DDMMYY_HHMM>/)
# Exit:   0 iff both configurations build and every selected test passes.
set -uo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
SRC="$(cd -- "${HERE}/../.." >/dev/null 2>&1 && pwd)"
STAMP="$(date +%d%m%y_%H%M)"
LOGDIR="${BLS_TEST_LOGDIR:-${SRC}/build-audit-logs/${STAMP}}"
mkdir -p "${LOGDIR}"

# UBSan only prints by default; make every sanitizer report a test failure.
export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=1:halt_on_error=1:abort_on_error=0}"
export UBSAN_OPTIONS="${UBSAN_OPTIONS:-halt_on_error=1:print_stacktrace=1}"

rc=0
echo "[run_all] building Debug+ASan/UBSan (build-audit-asan) ..."
if ! BLS_BUILD_SUBDIR=build-audit-asan "${SRC}/scripts/build.sh" -DCMAKE_BUILD_TYPE=Debug -DBLS_SANITIZE=ON \
     "-DCMAKE_CXX_FLAGS=-D_GLIBCXX_ASSERTIONS" > "${LOGDIR}/build_asan.log" 2>&1; then
  echo "[run_all] ASan build FAILED (see ${LOGDIR}/build_asan.log)"; exit 2
fi
echo "[run_all] building Release (build-audit) ..."
if ! BLS_BUILD_SUBDIR=build-audit "${SRC}/scripts/build.sh" -DCMAKE_BUILD_TYPE=Release \
     > "${LOGDIR}/build_release.log" 2>&1; then
  echo "[run_all] Release build FAILED (see ${LOGDIR}/build_release.log)"; exit 2
fi
grep -H . "${SRC}/build-audit/bls_build_info.txt" "${SRC}/build-audit-asan/bls_build_info.txt" > "${LOGDIR}/build_info.txt"

for cfg in build-audit-asan build-audit; do
  echo "[run_all] ctest in ${cfg} ..."
  ctest --test-dir "${SRC}/${cfg}" --output-on-failure "$@" > "${LOGDIR}/ctest_${cfg}.log" 2>&1
  r=$?
  tail -n 40 "${LOGDIR}/ctest_${cfg}.log" | grep -E "tests passed|tests failed|\*\*\*|Failed|Passed" | tail -n 30
  san=$(grep -cE "runtime error|ERROR: AddressSanitizer|ERROR: LeakSanitizer" "${LOGDIR}/ctest_${cfg}.log")
  echo "[run_all] ${cfg}: ctest exit ${r}, sanitizer reports ${san}"
  if [[ ${r} -ne 0 || ${san} -ne 0 ]]; then rc=1; fi
done
echo "[run_all] logs in ${LOGDIR}; overall exit ${rc}"
exit ${rc}
