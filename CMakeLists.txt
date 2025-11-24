cmake_minimum_required(VERSION 3.10)
project(PtolemyCpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include directories
include_directories(include)

# Source files
file(GLOB_RECURSE SOURCES src/*.cpp)

# Executable
add_executable(Cpp_AI
    src/main.cpp
    src/dwba/dwba.cpp
    src/dwba/potential_eval.cpp
    src/dwba/math_utils.cpp
    src/dwba/av18_potential.cpp
    src/dwba/rcwfn.cpp
    src/input/InputGenerator.cpp
    src/input/Potentials.cpp
    src/input/Isotope.cpp
    src/input/PtolemyParser.cpp
)

# Copy data files to build directory
file(COPY data/mass20.txt DESTINATION ${CMAKE_BINARY_DIR}/data)
