#!/bin/bash
# Build test_inelastic.cpp using DWBA WAVELJ for distorted waves
# Working as of 2026-04-15: mean DCS error 0.36%, max 2.73% vs Fortran
cd "$(dirname "$0")"
g++ -O2 -std=c++17 -Iinclude \
  test_inelastic.cpp \
  src/dwba/wavelj.cpp \
  src/dwba/bound.cpp \
  src/dwba/potential_eval.cpp \
  src/dwba/rcwfn.cpp \
  src/dwba/math_utils.cpp \
  src/dwba/setup.cpp \
  src/dwba/dwba.cpp \
  src/dwba/xsectn.cpp \
  src/dwba/spline.cpp \
  src/dwba/grdset_ineldc_faithful.cpp \
  src/dwba/ineldc_collective.cpp \
  src/dwba/a12.cpp \
  src/dwba/av18_potential.cpp \
  src/dwba/coulin.cpp \
  src/input/Isotope.cpp \
  -o test_inel -lm "$@"
