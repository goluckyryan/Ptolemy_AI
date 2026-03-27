// main.cpp — Generic Ptolemy-compatible DWBA driver
// Reads Ptolemy .in format from stdin or file argument.
//
// Usage:
//   ./ptolemy < input.in
//   ./ptolemy input.in
//
// Output: angle (deg) and dσ/dΩ (mb/sr) table to stdout.

#include "dwba.h"
#include "PtolemyParser.h"
#include <iostream>
#include <iomanip>
#include <cstdio>

int main(int argc, char* argv[]) {
    DWBA dwba;

    // Defaults (can be overridden by input)
    dwba.SetAngles(0.0, 180.0, 5.0);

    PtolemyParser parser;
    if (argc > 1) {
        parser.ParseFile(argv[1], dwba);
    } else {
        parser.ParseStdin(dwba);
    }

    // Print parsed parameters
    dwba.PrintParameters();

    // Run the calculation
    dwba.Calculate();

    return 0;
}
