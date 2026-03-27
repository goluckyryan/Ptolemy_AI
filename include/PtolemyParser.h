#ifndef PTOLEMY_PARSER_H
#define PTOLEMY_PARSER_H

#include <string>
#include <vector>
#include "dwba.h"

class PtolemyParser {
public:
    PtolemyParser();
    ~PtolemyParser();

    void ParseFile(const std::string& filename, DWBA& dwba);
    void ParseStdin(DWBA& dwba);

private:
    void ParseLines(const std::vector<std::string>& lines, DWBA& dwba);
    void ParseReactionLine(const std::string& line, DWBA& dwba);
    void ParseParameterSet(const std::string& line, DWBA& dwba);

    // Reaction isotope names (stored for auto-BE calculation)
    std::string rxn_target, rxn_proj, rxn_ejectile, rxn_residual;
    double rxn_excitation = 0;
};

#endif
