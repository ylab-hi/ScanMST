cmake_minimum_required(VERSION 3.18...3.31)
project(scanmst LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Policies
# ---------------------------------------------------------------------------
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW) # do not use the deprecated FindPythonInterp
endif()
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # honour <Pkg>_ROOT in find_package
endif()

option(SCANMST_EMBED_RPATH
       "Embed an RPATH to the htslib prefix in the installed extension" ON)

# ---------------------------------------------------------------------------
# Python + pybind11
# ---------------------------------------------------------------------------
set(PYBIND11_FINDPYTHON ON)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# scikit-build-core puts the build env's site-packages on CMAKE_PREFIX_PATH,
# but ask pybind11 directly so a bare `cmake -S . -B build` also works.
if(NOT DEFINED pybind11_DIR)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE _pybind11_cmakedir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(_pybind11_cmakedir)
    set(pybind11_DIR "${_pybind11_cmakedir}")
  endif()
endif()
find_package(pybind11 CONFIG REQUIRED)

# ---------------------------------------------------------------------------
# Toolchain / language standards
# ---------------------------------------------------------------------------
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT MSVC)
  add_compile_options(-fvisibility=hidden)
endif()

# ---------------------------------------------------------------------------
# htslib / zlib prefix resolution
#
#   priority:  -DHTSLIB_ROOT  >  $ENV{HTSLIB_ROOT}
#           >  -DCONDA_PREFIX >  $ENV{CONDA_PREFIX}
#
# Empty strings are treated as unset. This is the bug being fixed relative to
# the origin/feat/uv prototype: `if(DEFINED CONDA_PREFIX)` is TRUE for
# -DCONDA_PREFIX="", which skipped the intended FATAL_ERROR and instead failed
# later with an opaque find_path error. scikit-build-core always passes the
# define (with default = ""), so the empty case is the common one, not an edge.
#
# HTSLIB_ROOT is the escape hatch for environments with no conda at all --
# cibuildwheel/manylinux (Step 2) and bioconda (Step 4) both need it.
# ---------------------------------------------------------------------------
set(_HTS_PREFIX "")
foreach(_cand IN ITEMS "${HTSLIB_ROOT}" "$ENV{HTSLIB_ROOT}" "${CONDA_PREFIX}"
                       "$ENV{CONDA_PREFIX}")
  if(NOT _HTS_PREFIX AND _cand)
    set(_HTS_PREFIX "${_cand}")
  endif()
endforeach()

if(NOT _HTS_PREFIX)
  message(
    FATAL_ERROR
      "scanmst: could not locate an htslib prefix.\n"
      "Activate a conda environment providing htslib, or configure with\n"
      "  -DHTSLIB_ROOT=/path/to/prefix\n"
      "The prefix must contain include/htslib/hts.h and lib/libhts.*")
endif()

message(STATUS "scanmst: using htslib prefix ${_HTS_PREFIX}")

set(_HTS_INC "${_HTS_PREFIX}/include")
set(_HTS_LIB_DIR "${_HTS_PREFIX}/lib")

find_path(
  HTS_INCLUDE_DIR htslib/hts.h
  HINTS "${_HTS_INC}"
  REQUIRED)
find_library(
  HTS_LIBRARY
  NAMES hts
  HINTS "${_HTS_LIB_DIR}"
  REQUIRED)
find_library(
  ZLIB_LIBRARY
  NAMES z zlib
  HINTS "${_HTS_LIB_DIR}"
  REQUIRED)

# ---------------------------------------------------------------------------
# Sources -- exact parity with the retired build.py
#
# DELIBERATELY NOT COMPILED:
#   src/scanmst/cppext/src/binding.cpp  -- declares a SECOND
#       PYBIND11_MODULE(_cppext, m); linking it alongside bindings/_cppext.cpp
#       produces a duplicate PyInit__cppext symbol.
#   src/scanmst/cppext/src/backward.cpp -- optional stacktrace helper, never
#       built by build.py.
#   src/scanmst/cppext/src/main.cpp     -- standalone test harness, untracked.
# ---------------------------------------------------------------------------
set(SCANMST_SOURCES
    src/scanmst/cppext/src/bam.cpp src/scanmst/cppext/src/rescuer.cpp
    src/scanmst/cppext/src/ssw.c src/scanmst/cppext/src/ssw_cpp.cpp)

# build.py globbed this directory recursively; keep that behaviour so newly
# added bindings are picked up the same way.
file(
  GLOB_RECURSE
  SCANMST_BINDINGS
  CONFIGURE_DEPENDS
  "src/scanmst/cppext/bindings/*.c"
  "src/scanmst/cppext/bindings/*.cpp")

# ---------------------------------------------------------------------------
# The extension module
# ---------------------------------------------------------------------------
pybind11_add_module(_cppext ${SCANMST_SOURCES} ${SCANMST_BINDINGS})

target_include_directories(
  _cppext PRIVATE ${HTS_INCLUDE_DIR}
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/scanmst/cppext/include)

target_link_libraries(_cppext PRIVATE ${HTS_LIBRARY} ${ZLIB_LIBRARY})

# Parity with build.py's extra_compile_args=["-DNSCDEBUG"].
# NOTE: NSCDEBUG is referenced by no C/C++ source in this repo -- it is a no-op
# kept only for parity with the old build, and can be dropped later.
target_compile_definitions(_cppext PRIVATE NSCDEBUG=1)

# ---------------------------------------------------------------------------
# RPATH
#
# The old build only ever got an RPATH by accident: conda's compiler
# activation exports LDFLAGS containing -Wl,-rpath,$CONDA_PREFIX/lib for the
# *base* env, so the .so picked one up without ever asking for it. That is an
# artifact of the developer's shell, not of the build.
#
# BUILD_RPATH alone does not survive `cmake --install`, so without
# INSTALL_RPATH the installed .so would have no RPATH of its own and fail with
# "libhts.so.3: cannot open shared object file".
#
# Note the conda LDFLAGS entry is additive and still appears in local dev
# builds even with SCANMST_EMBED_RPATH=OFF. It is absent in the manylinux
# container, where the .so is left with no RPATH for auditwheel to rewrite.
#
# Step 2 (cibuildwheel) must configure with -DSCANMST_EMBED_RPATH=OFF and let
# auditwheel/delocate bundle libhts and rewrite the paths.
# ---------------------------------------------------------------------------
set_target_properties(_cppext PROPERTIES BUILD_RPATH "${_HTS_LIB_DIR}")

if(SCANMST_EMBED_RPATH)
  set_target_properties(
    _cppext
    PROPERTIES INSTALL_RPATH "${_HTS_LIB_DIR}" BUILD_WITH_INSTALL_RPATH OFF
               INSTALL_RPATH_USE_LINK_PATH OFF)
endif()

# ---------------------------------------------------------------------------
# Install layout -- LOAD-BEARING
#
# src/scanmst/cppext/__init__.py does
#     from scanmst._cppext.cppext import *
# so the .so must land at   <wheel-root>/scanmst/_cppext.<abi>.so
# i.e. a SIBLING of scanmst/__init__.py, NOT inside scanmst/cppext/.
#
# Getting this wrong yields a wheel that builds and installs cleanly and then
# raises ImportError on first use.
# ---------------------------------------------------------------------------
install(
  TARGETS _cppext
  LIBRARY DESTINATION scanmst
  RUNTIME DESTINATION scanmst)
