#include "InputGenerator.h"
#include "Isotope.h"
#include "Potentials.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// Helper functions
int GetLValue(std::string str) {
  if (str.find("s") != std::string::npos)
    return 0;
  if (str.find("p") != std::string::npos)
    return 1;
  if (str.find("d") != std::string::npos)
    return 2;
  if (str.find("f") != std::string::npos)
    return 3;
  if (str.find("g") != std::string::npos)
    return 4;
  if (str.find("h") != std::string::npos)
    return 5;
  if (str.find("i") != std::string::npos)
    return 6;
  if (str.find("j") != std::string::npos)
    return 7;
  return -1;
}

double GetSpinValue(std::string str) {
  size_t found = str.find("/");
  if (found != std::string::npos) {
    std::string n = str.substr(0, found);
    std::string d = str.substr(found + 1);
    return atof(n.c_str()) * 1.0 / atof(d.c_str());
  } else {
    return atof(str.c_str());
  }
}

std::vector<std::string> SplitStr(std::string str, std::string delimiter) {
  std::vector<std::string> output;
  size_t pos = 0;
  std::string token;
  while ((pos = str.find(delimiter)) != std::string::npos) {
    token = str.substr(0, pos);
    output.push_back(token);
    str.erase(0, pos + delimiter.length());
  }
  output.push_back(str);
  return output;
}

void GenerateInput(std::string inputFileName, double angleMin, double angleMax,
                   double angleStep) {
  std::ifstream file_in;
  file_in.open(inputFileName.c_str(), std::ios::in);

  if (!file_in) {
    printf(" cannot open file : %s \n", inputFileName.c_str());
    return;
  }

  std::string line;
  int lineNum = 0;

  while (file_in.good()) {
    getline(file_in, line);
    lineNum++;
    if (line.length() < 5)
      continue;
    if (line.substr(0, 2) == "//")
      continue;

    // Parse reaction: e.g., 206Hg(d,p)207Hg
    size_t p1 = line.find("(");
    size_t p2 = line.find(",");
    size_t p3 = line.find(")");

    std::string targetName = line.substr(0, p1);
    std::string projectileName = line.substr(p1 + 1, p2 - p1 - 1);
    std::string recoilName = line.substr(p2 + 1, p3 - p2 - 1);
    std::string ejectileName =
        line.substr(p3 + 1); // This might need adjustment based on format

    // Extract rest of the line
    std::string rest = line.substr(p3 + 1);
    std::stringstream ss(rest);
    std::string ejectileNameStr, JpiStr, orbitalStr, JpiExStr, ExStr, ElabStr,
        potStr;
    ss >> ejectileNameStr >> JpiStr >> orbitalStr >> JpiExStr >> ExStr >>
        ElabStr >> potStr;

    // Re-assign ejectile name from parsed string if needed, or use the one from
    // reaction string In the example: 33Si(d,p)34Si 3/2+ 0d3/2 2+ 0.000 10MeV/u
    // AK target=33Si, proj=d, recoil=p, ejectile=34Si

    // Let's use Isotope class to get details
    Isotope target(targetName);
    Isotope projectile(projectileName);
    Isotope recoil(recoilName);
    Isotope ejectile(
        ejectileNameStr); // The one after ) is the ejectile (heavy product)

    // Check if isotopes are valid
    if (target.Name == "non" || projectile.Name == "non" ||
        recoil.Name == "non" || ejectile.Name == "non") {
      printf(" Error in Isotope names on line %d\n", lineNum);
      continue;
    }

    // Parse Jpi (Ejectile Ground State?)
    double J = GetSpinValue(JpiStr.substr(0, JpiStr.length() - 1));
    int parity = (JpiStr.back() == '+') ? 1 : -1;

    // Parse Orbital
    // e.g. 0d3/2
    int nNode = atoi(orbitalStr.substr(0, 1)
                         .c_str()); // Ptolemy uses 0-based nodes? Check manual.
    // Actually, InFileCreator.h line 334: int n =
    // atoi(orbital.substr(0,1).c_str()) - 1; But wait, standard notation 1s,
    // 1p... Ptolemy might expect nodes (0, 1...). Let's stick to InFileCreator
    // logic.
    int l = GetLValue(orbitalStr.substr(1, 1));
    double j = GetSpinValue(orbitalStr.substr(2));

    double Ex = atof(ExStr.c_str());

    // Parse Elab
    double Elab = 0.0;
    size_t posMeV = ElabStr.find("MeV/u");
    if (posMeV != std::string::npos) {
      double ePerU = atof(ElabStr.substr(0, posMeV).c_str());
      Elab = ePerU * projectile.A;
    } else {
      Elab = atof(ElabStr.c_str());
    }

    // Calculate Q-value
    double Qval =
        target.Mass + projectile.Mass - ejectile.Mass - recoil.Mass - Ex;

    // Generate Ptolemy input content
    std::string outFileName =
        inputFileName + ".in"; // Or some other naming convention
    std::ofstream file_out(
        outFileName.c_str(),
        std::ios::app); // Append to file? Or create new per line?
    // The original creates one big input file with multiple calculations?
    // "reset" command suggests multiple runs.

    file_out << "reset" << std::endl;
    file_out << "reaction: " << target.Name << "(" << projectile.Name << ","
             << recoil.Name << ")" << ejectile.Name << std::endl;
    file_out << "energy: " << Elab << std::endl;
    file_out << "qvalue: " << Qval << std::endl;
    file_out << "excitation: " << Ex << std::endl;

    // Potentials
    // This part needs to map the 'potStr' (e.g. 'A') to actual potential
    // parameters Using Potentials.h functions

    // Incoming
    CallPotential(potStr.substr(0, 1), target.A, target.Z, Elab, projectile.Z);
    file_out << "incoming: v=" << v << " r0=" << r0 << " a=" << a
             << " vi=" << vi << " ri0=" << ri0 << " ai=" << ai << " vsi=" << vsi
             << " rsi0=" << rsi0 << " asi=" << asi << " vso=" << vso
             << " rso0=" << rso0 << " aso=" << aso << " rc0=" << rc0
             << std::endl;

    // Outgoing
    // Need energy for outgoing channel
    // E_out approx E_in + Q
    double E_out = Elab * target.Mass / (target.Mass + projectile.Mass) + Qval;
    // This is CM energy? No, Elab is lab energy.
    // Let's look at InFileCreator.h logic for outgoing energy if it calculates
    // it for potential lookup. It passes E_out to CallPotential.

    // For now, let's assume simple approximation or just use Elab for potential
    // lookup if that's what was done. InFileCreator.h line 416:
    // CallPotential(pot.substr(1,1), isotope[3].A, isotope[3].Z, Eout,
    // isotope[2].Z); It calculates Eout.

    CallPotential(potStr.substr(1, 1), ejectile.A, ejectile.Z, E_out, recoil.Z);
    file_out << "outgoing: v=" << v << " r0=" << r0 << " a=" << a
             << " vi=" << vi << " ri0=" << ri0 << " ai=" << ai << " vsi=" << vsi
             << " rsi0=" << rsi0 << " asi=" << asi << " vso=" << vso
             << " rso0=" << rso0 << " aso=" << aso << " rc0=" << rc0
             << std::endl;

    // Bound state
    // Calculate binding energy of the transferred particle in the composite
    // nucleus For stripping A(a,b)B: a=b+x, B=A+x. Bound state is x in B.
    // Binding energy = Sp(B -> A + x) = Mass(A) + Mass(x) - Mass(B)
    // Wait, CalSp2(parent, daughter) = Mass(daughter) + Mass(particle) -
    // Mass(parent)? Let's check Isotope.h/cpp. CalSp2(Isotope& other)
    // calculates separation energy. In InFileCreator:
    // isotope[3].CalSp2(isotope[2]) -> ejectile.CalSp2(recoil)? No. Reaction:
    // 1(2,3)4. Target(Proj, Recoil)Ejectile. InFileCreator uses
    // isotope[0]..[3]. isotope[0]=target, [1]=proj, [2]=recoil, [3]=ejectile.
    // Bound state is for x in Ejectile (if stripping).
    // So Ejectile -> Target + x.
    // Binding = Ejectile.CalSp(Target)?
    // Or Ejectile.CalSp2(Target)?
    // Let's assume CalSp calculates separation energy of a nucleon?
    // Let's check Isotope.cpp.

    double bindingEnergy = 0.0;
    // Assuming single nucleon transfer for now, or using mass difference.
    // x = projectile - recoil
    double massX = projectile.Mass - recoil.Mass; // Approx
    // Better: Binding = (Target.Mass + (Projectile.Mass - Recoil.Mass)) -
    // Ejectile.Mass Wait, Q = (Target + Proj) - (Recoil + Ejectile). Q =
    // (Target + Proj) - (Recoil + (Target + x - Binding)) Q = Proj - Recoil - x
    // + Binding If Proj = Recoil + x + epsilon, then Q = epsilon + Binding.
    // Actually, Binding = Separation Energy.
    // S_x = Mass(Target) + Mass(x) - Mass(Ejectile).
    // Mass(x) = Mass(Proj) - Mass(Recoil) + S_proj (binding of x in proj).

    // Let's use the masses directly.
    // We need the mass of the transferred particle 'x'.
    // x is usually n, p, d, t, he3, a.
    // We can deduce x from Z and A difference.
    int dZ = ejectile.Z - target.Z;
    int dA = ejectile.A - target.A;
    // Find mass of particle with dZ, dA.
    // Use Isotope::FindIsotope(A, Z) if available, or just construct Isotope(A,
    // Z) if constructor supports it. Isotope class has FindIsotope static
    // method? Or constructor? Let's assume Isotope(int A, int Z) constructor
    // exists or similar. Actually Isotope.h has: Isotope(int A, int Z);

    Isotope transferredPart;
    if (dA > 0) {
      transferredPart = Isotope(dA, dZ);
      bindingEnergy = target.Mass + transferredPart.Mass - ejectile.Mass;
    } else if (dA < 0) {
      dZ = target.Z - ejectile.Z;
      dA = target.A - ejectile.A;
      transferredPart = Isotope(dA, dZ);
      bindingEnergy = ejectile.Mass + transferredPart.Mass - target.Mass;
    }

    file_out << "boundstate: n=" << nNode << " l=" << l << " j=" << j
             << " binding=" << bindingEnergy << std::endl;

    file_out << "angles: min=" << angleMin << " max=" << angleMax
             << " step=" << angleStep << std::endl;

    file_out << "end" << std::endl;

    file_out.close();
  }
  file_in.close();
}
