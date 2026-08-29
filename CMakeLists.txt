cmake_minimum_required(VERSION 3.16)
project(bls_analyzer LANGUAGES CXX)

# Guard against a known trap: a build/ directory copied or symlinked in from a
# different checkout carries a CMakeCache.txt whose CMAKE_HOME_DIRECTORY still
# points at the OLD source tree. On an incremental `make` (as opposed to a
# fresh `cmake -S/-B`), CMake silently trusts that stale value instead of
# re-deriving it from where the build actually lives -- so the build quietly
# reconfigures against, and recompiles, the wrong source tree. This bit twice:
# once via a copied build/ directory, once via scripts/run_E0B.sh hardcoding a
# sibling path. Every documented workflow (Dockerfile, README) builds with the
# binary directory as an immediate subdirectory of the source tree, so the
# source directory must equal the binary directory's parent; if it doesn't,
# CMAKE_HOME_DIRECTORY (== CMAKE_SOURCE_DIR here, at the top level) does not
# resolve to this tree.
get_filename_component(_bls_binary_parent "${CMAKE_BINARY_DIR}" DIRECTORY)
get_filename_component(_bls_binary_parent "${_bls_binary_parent}" REALPATH)
get_filename_component(_bls_source_real "${CMAKE_SOURCE_DIR}" REALPATH)
if(NOT _bls_binary_parent STREQUAL _bls_source_real)
  message(FATAL_ERROR
    "CMAKE_HOME_DIRECTORY resolves to '${_bls_source_real}', but the build "
    "directory '${CMAKE_BINARY_DIR}' lives under '${_bls_binary_parent}' -- "
    "they must be the same tree. This is almost always a build/ directory "
    "copied in from a different checkout (its CMakeCache.txt still points at "
    "the old source tree) or an out-of-tree build outside the documented "
    "layout. Delete this build directory and reconfigure with:\n"
    "    cmake -S ${CMAKE_CURRENT_LIST_DIR} -B <this-tree>/build")
endif()
unset(_bls_binary_parent)
unset(_bls_source_real)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BLS_USE_XDRFILE "Enable XDRFILE support" OFF)
option(BLS_USE_TNG "Enable TNG support" OFF)
option(BLS_SANITIZE "Build with -fsanitize=address,undefined" OFF)

# Requires a build directory of its own (sanitizer flags apply to every
# object file, so they can't coexist with a plain build in the same tree):
#   BLS/scripts/build.sh -DBLS_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug \
#     -- with BLS_BUILD_SUBDIR=build-asan (see build.sh) or, directly:
#   cmake -S . -B build-asan -DBLS_SANITIZE=ON -DCMAKE_BUILD_TYPE=Debug
#   cmake --build build-asan
# For valgrind, no special flags are required -- any build works, though a
# Debug build gives more useful stack traces. See scripts/run_valgrind.sh.
if(BLS_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=address,undefined)
endif()

option(BLS_COVERAGE "Build with --coverage (gcov) instrumentation" OFF)
if(BLS_COVERAGE)
  add_compile_options(--coverage -O0 -g)
  add_link_options(--coverage)
endif()

set(core_sources
    src/config/Parser.cpp
    src/io/TrajectoryFactory.cpp
    src/io/XtcReader.cpp
    src/io/TrrReader.cpp
    src/io/TngReader.cpp
    src/io/GroReader.cpp
    src/io/PdbReader.cpp
    src/io/XyzReader.cpp
    src/io/MolReader.cpp
    src/io/Topology.cpp
    src/lattice/Basis.cpp
    src/lattice/Enumerator.cpp
    src/grid/Grid.cpp
    src/refine/SkipDFS.cpp
    src/bls/BLS.cpp
    src/bls/Comparison.cpp
    src/cluster/Algorithms.cpp
    src/cluster/UnionFindBasic.cpp
    src/util/RSS.cpp)

add_library(bls_core STATIC ${core_sources})
target_include_directories(bls_core PUBLIC src)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(bls_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(BLS_USE_XDRFILE)
  find_package(PkgConfig REQUIRED)
  set(_xdr_pkgconfig "${CMAKE_SOURCE_DIR}/libxdrfile/install/lib/pkgconfig")
  if(EXISTS "${_xdr_pkgconfig}/xdrfile.pc")
    if(DEFINED ENV{PKG_CONFIG_PATH} AND NOT "$ENV{PKG_CONFIG_PATH}" STREQUAL "")
      set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${_xdr_pkgconfig}")
    else()
      set(ENV{PKG_CONFIG_PATH} "${_xdr_pkgconfig}")
    endif()
  endif()
  pkg_check_modules(XDRFILE REQUIRED xdrfile)
  target_compile_definitions(bls_core PUBLIC USE_XDRFILE)
  target_include_directories(bls_core PUBLIC ${XDRFILE_INCLUDE_DIRS})
  target_link_directories(bls_core PUBLIC ${XDRFILE_LIBRARY_DIRS})
  target_link_libraries(bls_core PUBLIC ${XDRFILE_LIBRARIES})
  set(_bls_xdr_library_dirs "${XDRFILE_LIBRARY_DIRS}")
  string(REPLACE ";" ":" _bls_xdr_rpath "${_bls_xdr_library_dirs}")
endif()

if(BLS_USE_TNG)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(TNG REQUIRED tng_io)
  target_compile_definitions(bls_core PUBLIC USE_TNG)
  target_include_directories(bls_core PUBLIC ${TNG_INCLUDE_DIRS})
  target_link_libraries(bls_core PUBLIC ${TNG_LIBRARIES})
endif()

add_executable(bls_analyze src/main.cpp)
target_link_libraries(bls_analyze PRIVATE bls_core)
if(BLS_USE_XDRFILE)
  if(DEFINED _bls_xdr_rpath)
    set_target_properties(bls_core PROPERTIES BUILD_RPATH "${_bls_xdr_rpath}"
                                            INSTALL_RPATH "${_bls_xdr_rpath}")
    set_target_properties(bls_analyze PROPERTIES BUILD_RPATH "${_bls_xdr_rpath}"
                                                INSTALL_RPATH "${_bls_xdr_rpath}")
  endif()
endif()

enable_testing()

add_executable(test_basic tests/test_basic.cpp)
target_link_libraries(test_basic PRIVATE bls_core)
if(BLS_USE_XDRFILE AND DEFINED _bls_xdr_rpath)
  set_target_properties(test_basic PROPERTIES BUILD_RPATH "${_bls_xdr_rpath}"
                                             INSTALL_RPATH "${_bls_xdr_rpath}")
endif()
add_test(NAME bls_basic COMMAND test_basic)

add_executable(test_cc3d_fairness tests/test_cc3d_fairness.cpp)
target_link_libraries(test_cc3d_fairness PRIVATE bls_core)
if(BLS_USE_XDRFILE AND DEFINED _bls_xdr_rpath)
  set_target_properties(test_cc3d_fairness PROPERTIES BUILD_RPATH "${_bls_xdr_rpath}"
                                                      INSTALL_RPATH "${_bls_xdr_rpath}")
endif()
add_test(NAME cc3d_fairness COMMAND test_cc3d_fairness)

add_executable(test_geometry tests/test_geometry.cpp)
target_link_libraries(test_geometry PRIVATE bls_core)
if(BLS_USE_XDRFILE AND DEFINED _bls_xdr_rpath)
  set_target_properties(test_geometry PROPERTIES BUILD_RPATH "${_bls_xdr_rpath}"
                                                 INSTALL_RPATH "${_bls_xdr_rpath}")
endif()
add_test(NAME geometry COMMAND test_geometry)
