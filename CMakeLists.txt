cmake_minimum_required(VERSION 3.16)
project(bls_analyzer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BLS_USE_XDRFILE "Enable XDRFILE support" OFF)
option(BLS_USE_TNG "Enable TNG support" OFF)

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
