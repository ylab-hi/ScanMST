cmake_minimum_required(VERSION 3.18...3.30)
project(scanmst LANGUAGES C CXX)

# --- Modern Python discovery to avoid CMP0148 warning ---
if(POLICY CMP0148)
  cmake_policy(SET CMP0148 NEW)
endif()
set(PYBIND11_FINDPYTHON ON)
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

# ---- Toolchain / standards ----
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT MSVC)
  add_compile_options(-fvisibility=hidden)
endif()

# ---- pybind11 ----
find_package(pybind11 CONFIG REQUIRED)

# ---- Find HTSlib (conda or custom prefix) ----
if(DEFINED HTSLIB_ROOT)
  set(_HTS_INC "${HTSLIB_ROOT}/include")
  set(_HTS_LIB_DIR "${HTSLIB_ROOT}/lib")
elseif(DEFINED CONDA_PREFIX)
  set(_HTS_INC "${CONDA_PREFIX}/include")
  set(_HTS_LIB_DIR "${CONDA_PREFIX}/lib")
else()
  message(FATAL_ERROR "No CONDA_PREFIX and no HTSLIB_ROOT; activate the conda env or set HTSLIB_ROOT.")
endif()

find_path(HTS_INCLUDE_DIR htslib/hts.h HINTS "${_HTS_INC}" REQUIRED)
find_library(HTS_LIBRARY NAMES hts HINTS "${_HTS_LIB_DIR}" REQUIRED)
find_library(ZLIB_LIBRARY NAMES z zlib HINTS "${_HTS_LIB_DIR}" REQUIRED)

# ---- Sources ----
set(SCANMST_SOURCES
  src/scanmst/cppext/src/bam.cpp
  src/scanmst/cppext/src/rescuer.cpp
  src/scanmst/cppext/src/ssw.c
  src/scanmst/cppext/src/ssw_cpp.cpp
)
file(GLOB_RECURSE SCANMST_BINDINGS
  CONFIGURE_DEPENDS
  "src/scanmst/cppext/bindings/*.c"
  "src/scanmst/cppext/bindings/*.cpp"
)

# ---- Extension ----
pybind11_add_module(_cppext
  ${SCANMST_SOURCES}
  ${SCANMST_BINDINGS}
)

target_include_directories(_cppext
  PRIVATE
    ${HTS_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scanmst/cppext/include
)

target_link_libraries(_cppext PRIVATE ${HTS_LIBRARY} ${ZLIB_LIBRARY})
target_compile_definitions(_cppext PRIVATE NSCDEBUG=1)
set_target_properties(_cppext PROPERTIES BUILD_RPATH "${_HTS_LIB_DIR}")

install(TARGETS _cppext
  LIBRARY DESTINATION scanmst
  RUNTIME DESTINATION scanmst
)
