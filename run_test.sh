#!/bin/bash
# Build and run the C++ DWBA code
# Usage: ./run_test.sh [input.in] [angle_min] [angle_max] [angle_step]

INPUT="${1:-}"
AMIN="${2:-0}"
AMAX="${3:-180}"
ASTEP="${4:-5}"

# Compile
echo "Compiling..."
g++ -std=c++17 -O2 -Iinclude -DHAVE_INELDC_FR \
    src/main.cpp \
    src/dwba/wavelj.cpp \
    src/dwba/potential_eval.cpp \
    src/dwba/rcwfn.cpp \
    src/dwba/math_utils.cpp \
    src/dwba/bound.cpp \
    src/dwba/ineldc.cpp \
    src/dwba/ineldc_zr.cpp \
    src/dwba/a12.cpp \
    src/dwba/setup.cpp \
    src/dwba/grdset.cpp \
    src/dwba/xsectn.cpp \
    src/dwba/av18_potential.cpp \
    src/input/InputGenerator.cpp \
    src/input/Isotope.cpp \
    src/input/Potentials.cpp \
    src/input/PtolemyParser.cpp \
    -o dwba

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi
echo "Compiled OK -> dwba"

# Run
if [ -n "$INPUT" ]; then
    echo "Running: ./dwba $INPUT $AMIN $AMAX $ASTEP"
    ./dwba "$INPUT" "$AMIN" "$AMAX" "$ASTEP"
else
    echo "No input file specified. Usage: ./run_test.sh <input.in> [amin] [amax] [astep]"
fi
