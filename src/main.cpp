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

    // Determine the Ptolemy input file to parse
    std::string ptolemyInputFile;
    
    // If the file already ends with .in, use it directly
    if (inputFile.size() > 3 && inputFile.substr(inputFile.size()-3) == ".in") {
        ptolemyInputFile = inputFile;
    } else {
        // Generate Ptolemy input file
        GenerateInput(inputFile, angleMin, angleMax, angleStep);
        ptolemyInputFile = inputFile + ".in";
    }
    
    DWBA dwba;
    PtolemyParser parser;
    
    // Override angles from command line
    dwba.SetAngles(angleMin, angleMax, angleStep);
    
    std::cout << "Parsing input file: " << ptolemyInputFile << std::endl;
    parser.ParseFile(ptolemyInputFile, dwba);
    
    std::cout << "Running DWBA calculation..." << std::endl;
    dwba.Calculate();
    
    std::cout << "Calculation complete." << std::endl;

    return 0;
}
