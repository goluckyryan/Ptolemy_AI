// PtolemyParser.cpp — Parse Ptolemy .in input files
// Supports the standard Ptolemy input format with REACTION:, PARAMETERSET,
// PROJECTILE, TARGET, INCOMING, OUTGOING blocks terminated by ';'
//
// Usage: PtolemyParser parser; parser.ParseFile("input.in", dwba);

#include "PtolemyParser.h"
#include "Isotope.h"
#include "ptolemy_mass_table.h"
#include "elastic.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
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
    bool outVI_set = false, outVSO_set = false;  // track explicit set vs defaulted

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

        // Comment lines (Ptolemy uses ' and $ as comment chars)
        if (raw[0] == '\'' || raw[0] == '$') continue;

        // Strip inline $ comments (everything after $ is a comment)
        {
            size_t dollar = raw.find('$');
            if (dollar != std::string::npos) {
                raw = Trim(raw.substr(0, dollar));
                if (raw.empty()) continue;
            }
        }

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
                    // Compute BE from Ptolemy mass excesses (AME2003), matching Fortran SETCHN:
                    // BE = Sn - Ex = (ME_core + ME_trans - ME_res) - Ex
                    Isotope core(rxn_target);
                    Isotope residual(rxn_residual);
                    // For stripping (d,p): core=target, residual=residual, A_trans>0
                    // For pickup   (p,d): core=target(heavy), residual=residual(lighter)
                    // Always: lighter + neutron = heavier, so use abs()
                    int A_trans = std::abs(residual.A - core.A);
                    int Z_trans = std::abs(residual.Z - core.Z);
                    // lighter nucleus (core for pickup, residual for stripping)
                    Isotope &lighter = (residual.A < core.A) ? residual : core;
                    Isotope &heavier = (residual.A < core.A) ? core : residual;
                    double ME_lighter = PtolemyMass::MassExcess_MeV(lighter.Z, lighter.A);
                    double ME_heavier = PtolemyMass::MassExcess_MeV(heavier.Z, heavier.A);
                    double ME_trans   = PtolemyMass::MassExcess_MeV(Z_trans, A_trans);

                    // BE = separation energy = ME_lighter + ME_trans - ME_heavier - Ex
                    bs_binding = ME_lighter + ME_trans - ME_heavier - rxn_excitation;
                }
                targetBS_n = bs_nodes; targetBS_l = bs_l; targetBS_j = bs_j;
                targetBS_pot = bsPot;
                dwba.SetTargetBoundState(bs_nodes, bs_l, bs_j, bs_binding, bsPot);
            } else if (section == INCOMING) {
                dwba.SetIncomingPotential(inPot);
                if (isElastic_) elastic_pot_ = inPot;
            } else if (section == OUTGOING) {
                // Ptolemy convention: inherit missing VI and VSO from incoming when not explicitly set
                if (!outVI_set && inPot.VI != 0.0) {
                    outPot.VI = inPot.VI; outPot.RI0 = inPot.RI0; outPot.AI = inPot.AI;
                }
                if (!outVSO_set && inPot.VSO != 0.0) {
                    outPot.VSO = inPot.VSO; outPot.RSO0 = inPot.RSO0; outPot.ASO = inPot.ASO;
                    outPot.VSOI = inPot.VSOI; outPot.RSOI0 = inPot.RSOI0; outPot.ASOI = inPot.ASOI;
                }
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
            outVI_set = false; outVSO_set = false;
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
        // Pre-process: normalize "key = value" to "key=value"
        // This handles Ptolemy's loose formatting (spaces around '=')
        {
            std::string normalized;
            std::vector<std::string> parts;
            std::istringstream splitter(raw);
            std::string p;
            while (splitter >> p) parts.push_back(p);
            for (size_t pi = 0; pi < parts.size(); pi++) {
                if (pi > 0) normalized += ' ';
                // Pattern: TOKEN = VALUE  →  TOKEN=VALUE
                if (pi + 2 < parts.size() && parts[pi+1] == "=") {
                    normalized += parts[pi] + "=" + parts[pi+2];
                    pi += 2;
                }
                // Pattern: TOKEN= VALUE  →  TOKEN=VALUE
                else if (parts[pi].back() == '=' && pi + 1 < parts.size() && 
                         parts[pi+1].find('=') == std::string::npos) {
                    normalized += parts[pi] + parts[pi+1];
                    pi += 1;
                }
                else {
                    normalized += parts[pi];
                }
            }
            raw = normalized;
        }
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
                    if (key=="VI"||key=="W") outVI_set=true;
                    if (key=="VSO"||key=="VSOI") outVSO_set=true;
                } else if (section == PROJECTILE || section == TARGET) {
                    // Bound state parameters
                    if (key == "NODES" || key == "N") {
                        bs_nodes = (int)val;
                    } else if (key == "JBIGA" || key == "JBIGB" || key == "JA" || key == "JB") {
                        // Target/projectile spin: JBIGA=2*J_target, etc.
                        // Currently handled by auto-detection; store if needed later
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
                    else if (key == "ELAB")      { dwba.SetEnergy(val); elastic_Elab_ = val; }
                    else if (key == "LMIN")      dwba.SetLmin((int)val);
                    else if (key == "LMAX")      { dwba.SetLmax((int)val); elastic_Lmax_ = (int)val; }
                    else if (key == "ASYMPTOPIA") { dwba.SetAsymptopia(val); elastic_asymptopia_ = val; }
                    else if (key == "WYNN")      { elastic_wynn_ = (val != 0); }
                    else if (key == "LSTEP")     {} // default 1
                    else if (key == "MAXLEXTRAP") {} // ignore
                    else if (key == "BELX")      { dwba.BELx = val; }
                    else if (key == "LX")        { dwba.Lx  = (int)val; }
                    else if (key == "JBIGA")     {} // spin of target — handled by SetTargetSpin
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
                // Standalone bare keywords (also accepted on PARAMETERSET line)
                else if (utoken == "TMATCH")    { dwba.SetUseTMATCH(true); }
                else if (utoken == "WRONSKIAN") { dwba.SetUseTMATCH(false); }
            }
        }
    }

    // Apply angles
    if (anglesSet) {
        dwba.SetAngles(angleMin, angleMax, angleStep);
        if (isElastic_) {
            elastic_angleMin_ = angleMin;
            elastic_angleMax_ = angleMax;
            elastic_angleStep_ = angleStep;
        }
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

    // Extract everything before ELAB as the reaction string
    // Need to handle spaces inside parens: "16O(d,p)17O(1/2+ 0.871) ELAB=10"
    std::string rxn_str;
    {
        // Find ELAB (case-insensitive) to split reaction part from energy part
        std::string rest_upper = rest;
        std::transform(rest_upper.begin(), rest_upper.end(), rest_upper.begin(), ::toupper);
        size_t elab_pos = rest_upper.find("ELAB");
        if (elab_pos != std::string::npos) {
            rxn_str = Trim(rest.substr(0, elab_pos));
        } else {
            rxn_str = Trim(rest);
        }
    }
    std::stringstream ss(rest);

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

    // Detect elastic scattering: ejectile == projectile AND no excitation (elastic) or Ex=0
    // Note: inelastic (d,d') also has ejt==prj but with Ex>0 or BELX>0
    if (ejt == prj && excitation == 0.0) {
        isElastic_ = true;
        elastic_target_ = tgt;
        elastic_proj_ = prj;
        Isotope projIso(prj);
        Isotope tgtIso(tgt);
        elastic_Ap_ = projIso.A;
        elastic_Zp_ = projIso.Z;
        elastic_At_ = tgtIso.A;
        elastic_Zt_ = tgtIso.Z;
          }

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
        // For collective inelastic (ejt==prj, BELx>0): Lx = int(residualJ)
        // e.g. "3-" → residualJ=3 → Lx=3 (E3 transition)
        // Only set if BELx already parsed OR if this is inelastic reaction
              if (ejt == prj && excitation > 0.0) {
            dwba.Lx = (int)(residualJ + 0.5);  // round to nearest integer
        }
    }

    // Check for ELAB= on the same line (handles "ELAB=20" and "ELAB= 20")
    std::string tok;
    while (ss >> tok) {
        std::string key;
        double val;
        // Handle "ELAB= VALUE" (space after =)
        std::string utok = ToUpper(tok);
        if (utok == "ELAB" || utok == "ELAB=") {
            // Next token should be '=' or the value
            std::string next;
            if (utok == "ELAB" && ss >> next) {
                if (next[0] == '=') {
                    // "ELAB = 20" or "ELAB =20"
                    std::string valStr = next.substr(1);
                    if (valStr.empty() && ss >> valStr) {}
                    try { double e = std::stod(valStr); dwba.SetEnergy(e); elastic_Elab_ = e; } catch (...) {}
                }
            } else if (utok == "ELAB=" && ss >> next) {
                // "ELAB= 20"
                try { double e = std::stod(next); dwba.SetEnergy(e); elastic_Elab_ = e; } catch (...) {}
            }
            continue;
        }
        if (ParseKeyVal(tok, key, val)) {
            if (key == "ELAB") { dwba.SetEnergy(val); elastic_Elab_ = val; }
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
            else if (key == "ASYMPTOPIA")  dwba.SetAsymptopia(val);
            else if (key == "WYNN")        { elastic_wynn_ = (val != 0); }
            else if (key == "PRINT")       {} // ignore for now
        } else {
            // Bare keyword flags
            if      (utoken == "LABANGLES")  {} // TODO: implement lab angle conversion
            else if (utoken == "R0TARGET")   {} // already default behavior
            // ── Fortran parameterset grid presets (RGRIDS/IGRIDS) ─────────────────
            // Order: DWCUT SUMPTS STPSPR GAMSUM GAMDIF BNDASY ALMNMT ALMXMT SCTASY PHIMID AMDMLT
            //        NPSUM NPDIF NPPHI
            // STPSPR (step size = min(2pi/k,1)/STPSPR) from RGRIDS(3,row)
            else if (utoken == "ALPHA1") {
                dwba.SetSctAsy(-15.0); dwba.GrdSTEPSPR=12;
                dwba.GrdDWCUT=1e-3; dwba.GrdSUMPTS=6; dwba.GrdGAMSUM=3; dwba.GrdGAMDIF=5; dwba.GrdAMDMLT=2.0;
                dwba.GrdNPSUM=40; dwba.GrdNPDIF=25; dwba.GrdNPPHI=12;
            }
            else if (utoken == "ALPHA2") {
                dwba.SetSctAsy(-20.0); dwba.GrdSTEPSPR=20;
                dwba.GrdDWCUT=1e-4; dwba.GrdSUMPTS=7; dwba.GrdGAMSUM=3; dwba.GrdGAMDIF=5; dwba.GrdAMDMLT=2.0;
                dwba.GrdNPSUM=50; dwba.GrdNPDIF=30; dwba.GrdNPPHI=14;
            }
            else if (utoken == "ALPHA3") {
                dwba.SetSctAsy(-24.0); dwba.GrdSTEPSPR=25;
                dwba.GrdDWCUT=1e-5; dwba.GrdSUMPTS=8; dwba.GrdGAMSUM=3; dwba.GrdGAMDIF=5; dwba.GrdAMDMLT=2.0;
                dwba.GrdNPSUM=70; dwba.GrdNPDIF=35; dwba.GrdNPPHI=16;
            }
            else if (utoken == "DPSA") {
                dwba.SetSctAsy(-20.0); dwba.GrdSTEPSPR=5;
                dwba.GrdDWCUT=2e-5; dwba.GrdSUMPTS=6; dwba.GrdGAMSUM=2; dwba.GrdGAMDIF=12; dwba.GrdAMDMLT=0.90;
                dwba.GrdNPSUM=20; dwba.GrdNPDIF=20; dwba.GrdNPPHI=10;
            }
            else if (utoken == "DPDA") {
                dwba.SetSctAsy(-20.0); dwba.GrdSTEPSPR=5;
                dwba.GrdDWCUT=2e-6; dwba.GrdSUMPTS=6; dwba.GrdGAMSUM=2; dwba.GrdGAMDIF=20; dwba.GrdAMDMLT=0.90;
                dwba.GrdNPSUM=20; dwba.GrdNPDIF=30; dwba.GrdNPPHI=15;
            }
            else if (utoken == "DPSB") {
                dwba.SetSctAsy(-20.0); dwba.GrdSTEPSPR=8;
                dwba.GrdDWCUT=2e-6; dwba.GrdSUMPTS=8; dwba.GrdGAMSUM=2; dwba.GrdGAMDIF=12; dwba.GrdAMDMLT=0.90;
                dwba.GrdNPSUM=40; dwba.GrdNPDIF=40; dwba.GrdNPPHI=20;
            }
            else if (utoken == "DPDB") {
                dwba.SetSctAsy(-20.0); dwba.GrdSTEPSPR=8;
                dwba.GrdDWCUT=5e-7; dwba.GrdSUMPTS=8; dwba.GrdGAMSUM=2; dwba.GrdGAMDIF=20; dwba.GrdAMDMLT=0.90;
                dwba.GrdNPSUM=40; dwba.GrdNPDIF=60; dwba.GrdNPPHI=30;
            }
            else if (utoken == "NONLOCALITY") {} // ignore for now
            else if (utoken == "TMATCH")     { dwba.SetUseTMATCH(true); }
            else if (utoken == "WRONSKIAN")  { dwba.SetUseTMATCH(false); }
            else if (utoken == "INELOCA1" || utoken == "INELOCA2" || utoken == "INELOCA3") {
                dwba.ParameterSet = utoken;
                // STEPSPER=15 for INELOCA (source.f:28049-28057)
            }
        }
    }
}

// ── RunElastic: perform elastic scattering calculation using ElasticSolver ──

void PtolemyParser::RunElastic(const DWBA& dwba) {
    if (!isElastic_) return;

    ElasticSolver solver;
    solver.SetSystem(elastic_Ap_, elastic_Zp_, elastic_At_, elastic_Zt_, elastic_Elab_);

    // Grid: use Raphael-like defaults (h=0.1 fm, Rmax=30 fm)
    double h = 0.1;
    int N = 300;
    if (elastic_asymptopia_ > 0) {
        N = (int)(elastic_asymptopia_ / h + 0.5);
        if (N < 200) N = 200;
    }
    solver.SetGrid(h, N);

    // Lmax
    if (elastic_Lmax_ > 0) {
        solver.SetLmax(elastic_Lmax_);
    } else {
        solver.SetLmax(-1);  // auto
    }

    // Add potentials from the INCOMING OM
    const auto& p = elastic_pot_;
    if (p.V != 0 || p.VI != 0)
        solver.AddVolumeWS({-std::abs(p.V), -std::abs(p.VI)}, p.R0, p.A);
    if (p.VSI != 0)
        solver.AddSurfaceWS({0.0, -std::abs(p.VSI)}, p.RSI0, p.ASI);
    if (p.VSO != 0 || p.VSOI != 0)
        solver.AddSpinOrbit({-std::abs(p.VSO), -std::abs(p.VSOI)}, p.RSO0, p.ASO);
    if (p.RC0 > 0)
        solver.AddCoulomb(p.RC0);

    // Wynn epsilon option
    solver.SetWynn(elastic_wynn_);

    // Compute
    solver.CalcKinematics();
    solver.CalcScatteringMatrix();

    // Print header
    std::cout << std::endl;
    std::cout << "REACTION: " << elastic_proj_ << " + " << elastic_target_
              << "  ELASTIC SCATTERING" << std::endl;
    std::cout << "ELAB = " << std::fixed << std::setprecision(3) << elastic_Elab_ << " MeV" << std::endl;
    solver.PrintKinematics();

    // Print optical potential
    std::cout << std::endl;
    std::cout << "OPTICAL MODEL POTENTIAL:" << std::endl;
    if (p.V != 0 || p.VI != 0)
        std::cout << "  Volume WS:   V=" << p.V << "  R0=" << p.R0 << "  A=" << p.A
                  << "  VI=" << p.VI << std::endl;
    if (p.VSI != 0)
        std::cout << "  Surface WS:  VSI=" << p.VSI << "  RSI0=" << p.RSI0 << "  ASI=" << p.ASI << std::endl;
    if (p.VSO != 0 || p.VSOI != 0)
        std::cout << "  Spin-Orbit:  VSO=" << p.VSO << "  RSO0=" << p.RSO0 << "  ASO=" << p.ASO
                  << "  VSOI=" << p.VSOI << std::endl;
    if (p.RC0 > 0)
        std::cout << "  Coulomb:     RC0=" << p.RC0 << std::endl;

    // Print S-matrix (magnitude and phase)
    std::cout << std::endl;
    solver.PrintSMatrix();

    // Print DCS table
    double thMin = elastic_angleMin_;
    double thMax = elastic_angleMax_;
    double thStep = elastic_angleStep_;
    if (thMin < 0.5) thMin = thStep;  // avoid θ=0

    double k = solver.k();
    double eta = solver.eta();

    std::cout << std::endl;
    std::cout << std::setw(8) << "ANGLE" << std::setw(14) << "SIGMA/RUTH"
              << std::setw(14) << "SIGMA_CM" << std::setw(14) << "RUTHERFORD" << std::endl;
    std::cout << std::setw(8) << "(deg)" << std::setw(14) << ""
              << std::setw(14) << "(mb/sr)" << std::setw(14) << "(mb/sr)" << std::endl;

    double totalReaction = 0;
    for (double theta = thMin; theta <= thMax + 0.001; theta += thStep) {
        double sigma = solver.DCSUnpolarized(theta);

        // Rutherford cross section
        double th_rad = theta * M_PI / 180.0;
        double sin2 = std::sin(th_rad / 2.0);
        double ruth = (eta * eta) / (4.0 * k * k * sin2 * sin2 * sin2 * sin2) * 10.0;
        // fm²/sr → mb/sr: multiply by 10

        double ratio = (ruth > 0) ? sigma / ruth : 0;

        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << theta
                  << std::setprecision(4) << std::setw(14) << ratio
                  << std::scientific << std::setprecision(4) << std::setw(14) << sigma
                  << std::setw(14) << ruth << std::endl;
    }

    // Total reaction cross section from S-matrix
    // σ_R = (π/k²) Σ_{L,J} (2J+1) (1 - |S_{LJ}|²)
    double sigmaR = 0;
    double spin = solver.S();
    int Lmax = 0;
    // Find Lmax from S-matrix
    for (int L = 0; L < 200; L++) {
        bool found = false;
        for (int dj2 = -(int)(2*spin); dj2 <= (int)(2*spin); dj2 += 2) {
            int twoJ = 2*L + dj2;
            if (twoJ < 0) continue;
            auto S = solver.GetSMatrix(L, twoJ);
            if (std::abs(S) > 1e-15 || L < 5) {
                found = true;
                double J = twoJ / 2.0;
                sigmaR += (2*J + 1) * (1.0 - std::norm(S));
            }
        }
        if (!found && L > 5) break;
        Lmax = L;
    }
    sigmaR *= M_PI / (k * k) / (2.0 * spin + 1.0) * 10.0;  // fm² → mb

    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "TOTAL REACTION CROSS SECTION = " << sigmaR << " mb" << std::endl;
    std::cout << std::endl;
}
