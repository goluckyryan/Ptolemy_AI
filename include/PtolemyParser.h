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

    // Elastic scattering support
    bool IsElastic() const { return isElastic_; }
    void RunElastic(const DWBA& dwba);

private:
    void ParseLines(const std::vector<std::string>& lines, DWBA& dwba);
    void ParseReactionLine(const std::string& line, DWBA& dwba);
    void ParseParameterSet(const std::string& line, DWBA& dwba);

    // Reaction isotope names (stored for auto-BE calculation)
    std::string rxn_target, rxn_proj, rxn_ejectile, rxn_residual;
    double rxn_excitation = 0;

    // Elastic scattering data
    bool isElastic_ = false;
    std::string elastic_target_, elastic_proj_;
    int elastic_Ap_ = 0, elastic_Zp_ = 0, elastic_At_ = 0, elastic_Zt_ = 0;
    double elastic_Elab_ = 0;
    ChannelPotential elastic_pot_ = {};
    double elastic_angleMin_ = 0, elastic_angleMax_ = 180, elastic_angleStep_ = 5;
    int elastic_Lmax_ = -1;          // -1 = auto
    double elastic_asymptopia_ = 30; // Rmax proxy
};

#endif
