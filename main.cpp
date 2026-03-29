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
#include <fstream>
#include <sstream>
#include <complex>
#include <map>
#include <regex>

// Load Fortran elastic S-matrix from file for override
static void loadFortranSmat(const char* fn, DWBA& dwba) {
    std::ifstream f(fn);
    if (!f.is_open()) return;
    std::map<std::pair<int,int>, std::complex<double>> smat_in, smat_out;
    std::string line;
    std::regex rx(R"(FTN_SMAT ch=\s*(\d+)\s+L=\s*(\d+)\s+JP=\s*(\d+)\s+([-.+E\d]+)\s+([-.+E\d]+))");
    while (std::getline(f, line)) {
        std::smatch m;
        if (std::regex_search(line, m, rx)) {
            int ch = std::stoi(m[1]), L = std::stoi(m[2]), JP = std::stoi(m[3]);
            double sr = std::stod(m[4]), si = std::stod(m[5]);
            auto key = std::make_pair(L, JP);
            if (ch == 0) smat_in[key] = {sr, si};
            else         smat_out[key] = {sr, si};
        }
    }
    fprintf(stderr, "Loaded Fortran S-matrix: %zu incoming, %zu outgoing entries\n",
            smat_in.size(), smat_out.size());
    dwba.SetElasticSMatrixOverride(smat_in, smat_out);
}

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

#ifdef OVERRIDE_SMAT
    loadFortranSmat("/tmp/ftn_hg_smat_unique.txt", dwba);
#endif

    // Print parsed parameters
    dwba.PrintParameters();

    // Run the calculation
    dwba.Calculate();

    return 0;
}
