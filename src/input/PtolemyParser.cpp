#include "PtolemyParser.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

PtolemyParser::PtolemyParser() {}
PtolemyParser::~PtolemyParser() {}

void PtolemyParser::ParseFile(const std::string &filename, DWBA &dwba) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::string line;
  while (std::getline(file, line)) {
    ParseLine(line, dwba);
  }
}

void PtolemyParser::ParseLine(const std::string &line, DWBA &dwba) {
  if (line.empty() || line[0] == '\'')
    return;

  std::stringstream ss(line);
  std::string command;
  ss >> command;

  // Remove colon if present
  if (command.back() == ':')
    command.pop_back();

  if (command == "reset") {
    // Reset
  } else if (command == "reaction") {
    std::string reaction;
    ss >> reaction;
    // Parse reaction string if needed, e.g. 206Hg(d,p)207Hg
    // dwba.SetReaction(...) is called by main? No, main calls InputGenerator.
    // InputGenerator generates the .in file.
    // We need to set the reaction in DWBA.
    // Parse reaction string: Target(Proj,Recoil)Ejectile
    size_t p1 = reaction.find("(");
    size_t p2 = reaction.find(",");
    size_t p3 = reaction.find(")");
    if (p1 != std::string::npos && p2 != std::string::npos &&
        p3 != std::string::npos) {
      std::string target = reaction.substr(0, p1);
      std::string proj = reaction.substr(p1 + 1, p2 - p1 - 1);
      std::string recoil = reaction.substr(p2 + 1, p3 - p2 - 1);
      std::string ejectile = reaction.substr(p3 + 1);
      dwba.SetReaction(target, proj, ejectile, recoil);

      // Set default Projectile Bound State if (d,p)
      // Set default Projectile Bound State if (d,p)
      if ((proj == "d" || proj == "2H") && (recoil == "p" || recoil == "1H")) {
        ChannelPotential pot;
        pot.V = 72.15; // Reid soft core approx or similar?
        pot.R0 = 1.25;
        pot.A = 0.65; // Standard WS
        // Deuteron binding 2.224 MeV
        dwba.SetProjectileBoundState(0, 0, 0.5, 2.224, pot);
      }
    }
  } else if (command == "energy" || command == "elab") {
    double e;
    ss >> e;
    dwba.SetEnergy(e);
  } else if (command == "excitation") {
    double ex;
    ss >> ex;
    dwba.SetExcitation(ex);
  } else if (command == "angles") {
    std::string part;
    double min = 0, max = 180, step = 1;
    while (ss >> part) {
      if (part.find("min=") == 0)
        min = std::stod(part.substr(4));
      if (part.find("max=") == 0)
        max = std::stod(part.substr(4));
      if (part.find("step=") == 0)
        step = std::stod(part.substr(5));
    }
    dwba.SetAngles(min, max, step);
  } else if (command == "incoming" || command == "outgoing") {
    ChannelPotential pot;
    // Initialize
    pot.V = 0;
    pot.R0 = 0;
    pot.A = 0;
    pot.VI = 0;
    pot.RI0 = 0;
    pot.AI = 0;
    pot.VSO = 0;
    pot.RSO0 = 0;
    pot.ASO = 0;
    pot.VSOI = 0;
    pot.RSOI0 = 0;
    pot.ASOI = 0;
    pot.VSI = 0;
    pot.RSI0 = 0;
    pot.ASI = 0;
    pot.RC0 = 1.25;

    std::string part;
    while (ss >> part) {
      size_t eq = part.find('=');
      if (eq != std::string::npos) {
        std::string key = part.substr(0, eq);
        double val = std::stod(part.substr(eq + 1));
        // Convert key to lowercase for comparison if needed, but InputGenerator
        // uses lowercase InputGenerator: v=... r0=...
        if (key == "v")
          pot.V = val;
        else if (key == "r0")
          pot.R0 = val;
        else if (key == "a")
          pot.A = val;
        else if (key == "vi")
          pot.VI = val;
        else if (key == "ri0")
          pot.RI0 = val;
        else if (key == "ai")
          pot.AI = val;
        else if (key == "vso")
          pot.VSO = val;
        else if (key == "rso0")
          pot.RSO0 = val;
        else if (key == "aso")
          pot.ASO = val;
        else if (key == "vsoi")
          pot.VSOI = val;
        else if (key == "rsoi0")
          pot.RSOI0 = val;
        else if (key == "asoi")
          pot.ASOI = val;
        else if (key == "vsi")
          pot.VSI = val;
        else if (key == "rsi0")
          pot.RSI0 = val;
        else if (key == "asi")
          pot.ASI = val;
        else if (key == "rc0")
          pot.RC0 = val;
      }
    }
    if (command == "incoming")
      dwba.SetIncomingPotential(pot);
    else
      dwba.SetOutgoingPotential(pot);
  } else if (command == "qvalue") {
    // Stored via SetEnergy after energy is set; no-op here as kinematics
    // computed in SetEnergy using Q from reaction parser or explicit qvalue.
    // For now accept and ignore — Q is set by SetEnergy from reaction string.
    // TODO: pass Q explicitly to DWBA if needed.
  } else if (command == "boundstate" || command == "projectilews") {
    int n = 0, l = 0;
    double j = 0.5, binding = 0;
    // Bound state potential parameters (read from input)
    ChannelPotential pot;
    pot.V = 50.0;   // Initial guess — will be adjusted by bound state solver
    pot.R0 = 1.25;
    pot.A = 0.65;
    pot.VSO = 6.0;
    pot.RSO0 = 1.10;
    pot.ASO = 0.65;
    pot.RC0 = 1.30;
    pot.VI = 0; pot.RI0 = 0; pot.AI = 0;
    pot.VSI = 0; pot.RSI0 = 0; pot.ASI = 0;
    pot.VSOI = 0; pot.RSOI0 = 0; pot.ASOI = 0;

    std::string part;
    while (ss >> part) {
      size_t eq = part.find('=');
      if (eq != std::string::npos) {
        std::string key = part.substr(0, eq);
        std::string valstr = part.substr(eq + 1);
        if (key == "n")       n = std::stoi(valstr);
        else if (key == "l")  l = std::stoi(valstr);
        else if (key == "j")  j = std::stod(valstr);
        else if (key == "binding") binding = std::stod(valstr);
        else if (key == "v")   pot.V    = std::stod(valstr);
        else if (key == "r0")  pot.R0   = std::stod(valstr);
        else if (key == "a")   pot.A    = std::stod(valstr);
        else if (key == "vso") pot.VSO  = std::stod(valstr);
        else if (key == "rso0") pot.RSO0 = std::stod(valstr);
        else if (key == "aso") pot.ASO  = std::stod(valstr);
        else if (key == "rc0") pot.RC0  = std::stod(valstr);
      }
    }

    if (command == "boundstate")
      dwba.SetTargetBoundState(n, l, j, binding, pot);
    else  // projectilews
      dwba.SetProjectileBoundState(n, l, j, binding, pot);
  } else if (command == "target_spin") {
    double spin;
    ss >> spin;
    dwba.SetTargetSpin(spin);
  } else if (command == "residual_spin") {
    double spin;
    ss >> spin;
    dwba.SetResidualSpin(spin);
  } else if (command == "projectile_wf_file") {
    // Read external projectile bound state wavefunction
    // Format: projectile_wf_file: <path> h=<step> spam=<spectroscopic_amplitude>
    std::string path;
    ss >> path;
    double h_ext = 0.0625;
    double spam_ext = 1.0;
    std::string part;
    while (ss >> part) {
      if (part.find("h=") == 0) h_ext = std::stod(part.substr(2));
      if (part.find("spam=") == 0) spam_ext = std::stod(part.substr(5));
    }
    // Read the WF file: columns are (index, r, phi=u/r)
    std::ifstream wffile(path);
    if (!wffile.is_open()) {
      std::cerr << "Warning: Cannot open projectile_wf_file: " << path << std::endl;
    } else {
      std::vector<std::pair<double,double>> wf_data;
      std::string wfline;
      while (std::getline(wffile, wfline)) {
        std::stringstream wfss(wfline);
        int idx; double r, phi;
        if (wfss >> idx >> r >> phi) {
          wf_data.push_back({r, phi});
        }
      }
      if (!wf_data.empty()) {
        dwba.SetProjectileWFFromFile(wf_data, h_ext, spam_ext);
      }
    }
  }
}
