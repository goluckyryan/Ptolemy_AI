// main.cpp — Ptolemy++ main entry point
// Supports both DWBA transfer and elastic scattering modes.
// Usage:
//   ptolemy++ input.in                      (parse .in file)
//   ptolemy++ < input.in                    (read from stdin)
//   ptolemy++ input.in 0 180 1              (override angles)

#include <iostream>
#include <string>
#include <vector>
#include "InputGenerator.h"
#include "PtolemyParser.h"
#include "dwba.h"
#include "elastic.h"

int main(int argc, char* argv[]) {
    // Default values
    std::string inputFile;
    double angleMin = 0.0;
    double angleMax = 180.0;
    double angleStep = 1.0;

    // Parse command line arguments
    if (argc > 1) inputFile = argv[1];
    if (argc > 2) angleMin = std::stod(argv[2]);
    if (argc > 3) angleMax = std::stod(argv[3]);
    if (argc > 4) angleStep = std::stod(argv[4]);

    DWBA dwba;
    PtolemyParser parser;
    
    // Override angles from command line
    dwba.SetAngles(angleMin, angleMax, angleStep);

    if (inputFile.empty()) {
        // Read from stdin
        parser.ParseStdin(dwba);
    } else if (inputFile.size() > 3 && inputFile.substr(inputFile.size()-3) == ".in") {
        parser.ParseFile(inputFile, dwba);
    } else {
        // Generate Ptolemy input file from DWBA spec
        GenerateInput(inputFile, angleMin, angleMax, angleStep);
        std::string ptolemyInputFile = inputFile + ".in";
        parser.ParseFile(ptolemyInputFile, dwba);
    }

    // Check if this is an elastic scattering calculation
    if (parser.IsElastic()) {
        std::cout << "=== Ptolemy++ Elastic Scattering ===" << std::endl;
        parser.RunElastic(dwba);
    } else {
        std::cout << "=== Ptolemy++ DWBA Transfer ===" << std::endl;
        dwba.Calculate();
    }

    return 0;
}
