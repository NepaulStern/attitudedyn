cmake_minimum_required(VERSION 3.15...3.28)
project(attitudedyn_python LANGUAGES CXX)

# The extension is a binding. The likelihood, gradient, draws, and
# post-estimation all come from the shared C++17 library in ../src/core.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED SKBUILD_PLATLIB_DIR)
  message(FATAL_ERROR "This CMake project is meant to be configured by scikit-build-core (pip install ./python).")
endif()

# A wheel sdist carries a copy of src/core at ./core. A checkout install
# (`pip install ./python`) uses the core at ../src/core.
set(ATTDYN_CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/core")
if(NOT EXISTS "${ATTDYN_CORE_DIR}/include/attitudedyn/api.hpp")
  set(ATTDYN_CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src/core")
endif()
if(NOT EXISTS "${ATTDYN_CORE_DIR}/include/attitudedyn/api.hpp")
  message(FATAL_ERROR
    "The attitudedyn C++ core was not found. "
    "Install from a checkout of the repository: pip install ./python")
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE ATTDYN_NUMPY_INCLUDE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY
)

pybind11_add_module(_core MODULE
  src/bindings.cpp
  "${ATTDYN_CORE_DIR}/src/layout.cpp"
  "${ATTDYN_CORE_DIR}/src/draws.cpp"
  "${ATTDYN_CORE_DIR}/src/likelihood.cpp"
  "${ATTDYN_CORE_DIR}/src/postestimation.cpp"
)

target_include_directories(_core PRIVATE
  "${ATTDYN_CORE_DIR}/include"
  "${ATTDYN_CORE_DIR}/src"
  "${ATTDYN_NUMPY_INCLUDE}"
)

target_compile_definitions(_core PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)

if(MSVC)
  target_compile_options(_core PRIVATE /W3 /EHsc /bigobj)
else()
  target_compile_options(_core PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

find_package(Threads REQUIRED)
target_link_libraries(_core PRIVATE Threads::Threads)

install(TARGETS _core DESTINATION "${SKBUILD_PLATLIB_DIR}/attitudedyn")
