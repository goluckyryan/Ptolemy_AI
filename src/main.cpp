#include <iostream>
#include <string>
#include <vector>
#include "InputGenerator.h"
#include "PtolemyParser.h"
#include "dwba.h"

int main(int argc, char* argv[]) {
    // Default values
    std::string inputFile = "digios/analysis/working/DWBA";
    double angleMin = 0.0;
    double angleMax = 180.0;
    double angleStep = 1.0;

    // Parse command line arguments
    if (argc > 1) inputFile = argv[1];
    if (argc > 2) angleMin = std::stod(argv[2]);
    if (argc > 3) angleMax = std::stod(argv[3]);
    if (argc > 4) angleStep = std::stod(argv[4]);

    std::cout << "Input File: " << inputFile << std::endl;
    std::cout << "Angle Range: " << angleMin << " - " << angleMax << " step " << angleStep << std::endl;

    // Generate Ptolemy input file
    GenerateInput(inputFile, angleMin, angleMax, angleStep);
    
    // Parse the generated input file and run DWBA
    std::string ptolemyInputFile = inputFile + ".in";
    
    DWBA dwba;
    PtolemyParser parser;
    
    std::cout << "Parsing generated input file: " << ptolemyInputFile << std::endl;
    parser.ParseFile(ptolemyInputFile, dwba);
    
    std::cout << "Running DWBA calculation..." << std::endl;
    dwba.Calculate();
    
    std::cout << "Calculation complete." << std::endl;

    return 0;
}
