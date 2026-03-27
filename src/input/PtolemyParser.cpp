// PtolemyParser.cpp — Parse Ptolemy .in input files
// Supports the standard Ptolemy input format with REACTION:, PARAMETERSET,
// PROJECTILE, TARGET, INCOMING, OUTGOING blocks terminated by ';'
//
// Usage: PtolemyParser parser; parser.ParseFile("input.in", dwba);

#include "PtolemyParser.h"
#include "Isotope.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

PtolemyParser::PtolemyParser() {}
PtolemyParser::~PtolemyParser() {}

// ── Helpers ──────────────────────────────────────────────────────────────────

static std::string ToUpper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

static std::string Trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Parse "KEY=VALUE" from a token, case-insensitive key.
// Returns true if found, sets key (uppercased) and val.
static bool ParseKeyVal(const std::string &token, std::string &key, double &val) {
    size_t eq = token.find('=');
    if (eq == std::string::npos || eq == 0) return false;
    key = ToUpper(token.substr(0, eq));
    try {
        val = std::stod(token.substr(eq + 1));
    } catch (...) {
        return false;
    }
    return true;
}

// Parse "KEY=STRING" from a token.
static bool ParseKeyStr(const std::string &token, std::string &key, std::string &val) {
    size_t eq = token.find('=');
    if (eq == std::string::npos || eq == 0) return false;
    key = ToUpper(token.substr(0, eq));
    val = token.substr(eq + 1);
    return true;
}

// Apply a KEY=VALUE pair to a ChannelPotential.
static void ApplyPotParam(ChannelPotential &pot, const std::string &key, double val) {
    if      (key == "V")     pot.V     = val;
    else if (key == "R0")    pot.R0    = val;
    else if (key == "A")     pot.A     = val;
    else if (key == "VI")    pot.VI    = val;
    else if (key == "W")     pot.VI    = val;   // W is alias for VI (volume imag)
    else if (key == "RI0")   pot.RI0   = val;
    else if (key == "AI")    pot.AI    = val;
    else if (key == "VSO")   pot.VSO   = val;
    else if (key == "RSO0")  pot.RSO0  = val;
    else if (key == "ASO")   pot.ASO   = val;
    else if (key == "VSOI")  pot.VSOI  = val;
    else if (key == "RSOI0") pot.RSOI0 = val;
    else if (key == "ASOI")  pot.ASOI  = val;
    else if (key == "VSI")   pot.VSI   = val;
    else if (key == "WD")    pot.VSI   = val;   // WD is alias for VSI (surface imag)
    else if (key == "RSI0")  pot.RSI0  = val;
    else if (key == "ASI")   pot.ASI   = val;
    else if (key == "RC0")   pot.RC0   = val;
}

// ── Main entry point ─────────────────────────────────────────────────────────

void PtolemyParser::ParseFile(const std::string &filename, DWBA &dwba) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "PtolemyParser: cannot open " << filename << "\n";
        return;
    }

    // Read all lines
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    ParseLines(lines, dwba);
}

void PtolemyParser::ParseStdin(DWBA &dwba) {
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(std::cin, line)) {
        lines.push_back(line);
    }
    ParseLines(lines, dwba);
}

void PtolemyParser::ParseLines(const std::vector<std::string> &lines, DWBA &dwba) {

    // State machine
    enum Section { NONE, PROJECTILE, TARGET, INCOMING, OUTGOING };
    Section section = NONE;

    // Temporaries for bound-state blocks
    ChannelPotential bsPot = {};
    int bs_nodes = 0, bs_l = 0;
    double bs_j = 0.5, bs_binding = 0.0, bs_v = 0.0;
    bool bs_v_set = false;
    std::string prj_wavefunction;  // "av18", "reid", etc.

    // Temporaries for optical potentials (accumulated across multi-line blocks)
    ChannelPotential inPot = {}, outPot = {};

    // Angle tracking (set angles at the end)
    double angleMin = 0, angleMax = 180, angleStep = 5;
    bool anglesSet = false;

    // Track target BS for potential reuse
    int targetBS_n = 0, targetBS_l = 0;
    double targetBS_j = 0.5;
    ChannelPotential targetBS_pot = {};



    for (size_t i = 0; i < lines.size(); i++) {
        std::string raw = Trim(lines[i]);
        if (raw.empty()) continue;

        // Comment lines (Ptolemy uses ' as comment)
        if (raw[0] == '\'') continue;

        std::string upper = ToUpper(raw);

        // ── Section terminators ──
        if (raw == ";" || upper == "END" || upper == "END;") {
            if (section == PROJECTILE) {
                // Apply projectile bound state
                if (prj_wavefunction == "AV18" || prj_wavefunction == "av18") {
                    // AV18 deuteron: the DWBA code handles this via LoadDeuteronWavefunction
                    // Just set the bound state params for the potential bracket
                    bsPot.V = bs_v_set ? bs_v : 72.0;
                    dwba.SetProjectileBoundState(bs_nodes, bs_l, bs_j, 2.224, bsPot);
                } else {
                    bsPot.V = bs_v_set ? bs_v : 72.0;
                    dwba.SetProjectileBoundState(bs_nodes, bs_l, bs_j, bs_binding > 0 ? bs_binding : 2.224, bsPot);
                }
            } else if (section == TARGET) {
                bsPot.V = bs_v_set ? bs_v : 50.0;
                // Auto-compute binding energy if not explicitly given
                if (bs_binding < 1e-6 && !rxn_target.empty() && !rxn_residual.empty()) {
                    // For (d,p): BE = M(target) + M(n) - M(residual) - Ex
                    Isotope core(rxn_target);
                    Isotope residual(rxn_residual);
                    const double mn = 939.565346;  // neutron mass in MeV
                    bs_binding = core.Mass + mn - residual.Mass - rxn_excitation;
                    std::cerr << "Auto-computed target binding energy: "
                              << bs_binding << " MeV (Sn of " << rxn_residual << ")\n";
                }
                targetBS_n = bs_nodes; targetBS_l = bs_l; targetBS_j = bs_j;
                targetBS_pot = bsPot;
                dwba.SetTargetBoundState(bs_nodes, bs_l, bs_j, bs_binding, bsPot);
            } else if (section == INCOMING) {
                dwba.SetIncomingPotential(inPot);
            } else if (section == OUTGOING) {
                dwba.SetOutgoingPotential(outPot);
            }
            section = NONE;
            continue;
        }

        // ── Section headers ──
        if (upper == "PROJECTILE") {
            section = PROJECTILE;
            bsPot = {};
            bs_nodes = 0; bs_l = 0; bs_j = 0.5; bs_binding = 0; bs_v = 0; bs_v_set = false;
            prj_wavefunction = "";
            continue;
        }
        if (upper == "TARGET") {
            section = TARGET;
            bsPot = {};
            bs_nodes = 0; bs_l = 0; bs_j = 0.5; bs_binding = 0; bs_v = 0; bs_v_set = false;
            continue;
        }
        if (upper == "INCOMING") {
            section = INCOMING;
            inPot = {};
            inPot.RC0 = 1.25;  // default
            continue;
        }
        if (upper == "OUTGOING") {
            section = OUTGOING;
            outPot = {};
            outPot.RC0 = 1.25;  // default
            continue;
        }

        // ── REACTION: line ──
        if (upper.substr(0, 9) == "REACTION:") {
            ParseReactionLine(raw, dwba);
            continue;
        }

        // ── PARAMETERSET line ──
        if (upper.substr(0, 12) == "PARAMETERSET") {
            ParseParameterSet(raw, dwba);
            continue;
        }

        // ── RESET ──
        if (upper == "RESET") {
            // Reset is a no-op for us (DWBA starts clean)
            continue;
        }

        // ── HEADER (comment) ──
        if (upper.substr(0, 6) == "HEADER") {
            continue;
        }

        // ── Standalone key=value lines (angles, print, etc.) ──
        // Tokenize the line and process each KEY=VALUE pair
        std::stringstream ss(raw);
        std::string token;
        while (ss >> token) {
            std::string key;
            double val;
            std::string sval;

            // First try string-value parsing for JP (handles "5/2" format)
            if (ParseKeyStr(token, key, sval)) {
                if ((key == "JP" || key == "J") && (section == PROJECTILE || section == TARGET)) {
                    size_t slash = sval.find('/');
                    if (slash != std::string::npos) {
                        bs_j = std::stod(sval.substr(0, slash)) / std::stod(sval.substr(slash + 1));
                    } else {
                        bs_j = std::stod(sval);
                    }
                    continue;  // handled
                } else if (key == "WAVEFUNCTION" && section == PROJECTILE) {
                    prj_wavefunction = ToUpper(sval);
                    continue;
                }
            }

            if (ParseKeyVal(token, key, val)) {
                if (section == INCOMING) {
                    ApplyPotParam(inPot, key, val);
                } else if (section == OUTGOING) {
                    ApplyPotParam(outPot, key, val);
                } else if (section == PROJECTILE || section == TARGET) {
                    // Bound state parameters
                    if (key == "NODES" || key == "N") {
                        bs_nodes = (int)val;
                    } else if (key == "L") {
                        bs_l = (int)val;
                    } else if (key == "BINDING" || key == "BE") {
                        bs_binding = val;
                    } else if (key == "V") {
                        bs_v = val; bs_v_set = true;
                    } else {
                        // Could be R0, A, RC0, VSO, etc. for bound-state potential
                        ApplyPotParam(bsPot, key, val);
                    }
                } else {
                    // Top-level key=value
                    if (key == "ANGLEMIN")       { angleMin = val; anglesSet = true; }
                    else if (key == "ANGLEMAX")  { angleMax = val; anglesSet = true; }
                    else if (key == "ANGLESTEP") { angleStep = val; anglesSet = true; }
                    else if (key == "ELAB")      dwba.SetEnergy(val);
                    else if (key == "LMIN")      dwba.SetLmin((int)val);
                    else if (key == "LMAX")      dwba.SetLmax((int)val);
                    else if (key == "PRINT")     {} // ignore for now
                }
            } else {
                // Bare keyword (no '=')
                std::string utoken = ToUpper(token);
                if (utoken == "WAVEFUNCTION" && section == PROJECTILE) {
                    // Next token is the wavefunction type
                    if (ss >> token) {
                        prj_wavefunction = ToUpper(token);
                    }
                }
                // Other bare keywords: r0target, labangles, dpsb, etc. — handled in PARAMETERSET
            }
        }
    }

    // Apply angles
    if (anglesSet) {
        dwba.SetAngles(angleMin, angleMax, angleStep);
    }

}

// ── Parse REACTION: line ─────────────────────────────────────────────────────
// Format: REACTION: 16O(d,p)17O(5/2+ 0.000) ELAB=20.0
// Also:   REACTION: 16O(d,p)17O  ELAB=20.0

void PtolemyParser::ParseReactionLine(const std::string &line, DWBA &dwba) {
    // Skip "REACTION:" prefix
    size_t colon = line.find(':');
    if (colon == std::string::npos) return;
    std::string rest = Trim(line.substr(colon + 1));

    std::stringstream ss(rest);
    std::string rxn_str;
    ss >> rxn_str;

    // Parse reaction string: Target(proj,eject)Residual or Target(proj,eject)Residual(JP Ex)
    // Examples: 16O(d,p)17O    16O(d,p)17O(5/2+ 0.000)
    size_t p1 = rxn_str.find('(');
    size_t p2 = rxn_str.find(',');
    size_t p3 = rxn_str.find(')');
    if (p1 == std::string::npos || p2 == std::string::npos || p3 == std::string::npos) return;

    std::string target    = rxn_str.substr(0, p1);
    std::string proj      = rxn_str.substr(p1 + 1, p2 - p1 - 1);
    std::string ejectile  = rxn_str.substr(p2 + 1, p3 - p2 - 1);

    // Residual: everything after first ')', possibly with (JP Ex) appended
    std::string after = rxn_str.substr(p3 + 1);
    std::string residual;
    double residualJ = -1.0;
    double excitation = 0.0;

    size_t p4 = after.find('(');
    if (p4 != std::string::npos) {
        residual = after.substr(0, p4);
        // Parse (5/2+ 0.000) — state info
        size_t p5 = after.find(')', p4);
        if (p5 != std::string::npos) {
            std::string stateStr = after.substr(p4 + 1, p5 - p4 - 1);
            std::stringstream ss2(stateStr);
            std::string jpStr;
            ss2 >> jpStr;
            // Parse JP string like "5/2+" or "0+"
            std::string jpNum = jpStr;
            // Remove parity suffix (+/-)
            if (!jpNum.empty() && (jpNum.back() == '+' || jpNum.back() == '-')) {
                jpNum.pop_back();
            }
            size_t slash = jpNum.find('/');
            if (slash != std::string::npos) {
                residualJ = std::stod(jpNum.substr(0, slash)) / std::stod(jpNum.substr(slash + 1));
            } else {
                residualJ = std::stod(jpNum);
            }
            if (ss2 >> excitation) {
                dwba.SetExcitation(excitation);
            }
        }
    } else {
        residual = after;
    }

    // Map common names: d→2H, p→1H, t→3H, 3He→3He, a/alpha→4He
    auto MapIsotope = [](const std::string &name) -> std::string {
        std::string u = name;
        std::transform(u.begin(), u.end(), u.begin(), ::tolower);
        if (u == "d" || u == "deuteron") return "2H";
        if (u == "p" || u == "proton")   return "1H";
        if (u == "t" || u == "triton")   return "3H";
        if (u == "a" || u == "alpha")    return "4He";
        if (u == "n" || u == "neutron")  return "n";
        return name;  // Already in AZZ format like "16O"
    };

    std::string tgt  = MapIsotope(target);
    std::string prj  = MapIsotope(proj);
    std::string ejt  = MapIsotope(ejectile);
    std::string res  = MapIsotope(residual);

    // SetReaction(target, projectile, ejectile, recoil)
    // In Ptolemy .in: Target(proj,light_out)Residual
    // SetReaction convention: ejectile=heavy_residual, recoil=light_outgoing
    dwba.SetReaction(tgt, prj, res, ejt);

    // Store for auto-BE calculation (member variables)
    this->rxn_target = tgt;
    this->rxn_proj = prj;
    this->rxn_ejectile = ejt;
    this->rxn_residual = res;
    this->rxn_excitation = excitation;

    if (residualJ >= 0) {
        dwba.SetResidualSpin(residualJ);
    }

    // Check for ELAB= on the same line
    std::string tok;
    while (ss >> tok) {
        std::string key;
        double val;
        if (ParseKeyVal(tok, key, val)) {
            if (key == "ELAB") dwba.SetEnergy(val);
        }
    }
}

// ── Parse PARAMETERSET line ──────────────────────────────────────────────────
// Format: PARAMETERSET dpsb labangles r0target lstep=1 lmin=0 lmax=40 maxlextrap=0 asymptopia=30
// Tokens can be bare keywords (flags) or KEY=VALUE pairs.

void PtolemyParser::ParseParameterSet(const std::string &line, DWBA &dwba) {
    std::stringstream ss(line);
    std::string first;
    ss >> first;  // skip "PARAMETERSET"

    std::string token;
    while (ss >> token) {
        std::string key;
        double val;
        std::string utoken = ToUpper(token);

        if (ParseKeyVal(token, key, val)) {
            if      (key == "LMIN")        dwba.SetLmin((int)val);
            else if (key == "LMAX")        dwba.SetLmax((int)val);
            else if (key == "LSTEP")       {} // default 1, ignore for now
            else if (key == "MAXLEXTRAP")  {} // ignore for now
            else if (key == "ASYMPTOPIA")  {} // ignore for now
            else if (key == "PRINT")       {} // ignore for now
        } else {
            // Bare keyword flags
            if      (utoken == "LABANGLES")  {} // TODO: implement lab angle conversion
            else if (utoken == "R0TARGET")   {} // already default behavior
            else if (utoken == "DPSB")       {} // deuteron projectile stripping/breakup — default
            else if (utoken == "NONLOCALITY") {} // ignore for now
        }
    }
}


