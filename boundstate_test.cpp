// boundstate_test.cpp — 16O(d,p)17O bound state: neutron in 16O core
// Computes phi = u/r for 0d5/2 (gs) and 1s1/2 (Ex=0.871 MeV)
// Uses DWBA::CalculateBoundState from bound.cpp
// Potential params matching Raphael dwba_zr.py defaults:
//   r0=1.10, a=0.65, Vso=-6, rso=1.25, aso=0.65, rc=1.30

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "include/dwba.h"

int main() {
    const double AMU   = 931.494061;
    const double HBARC = 197.32697;

    // ── Neutron (transferred nucleon) bound in 16O core ──────────────
    // mu = m_n * m_16O / (m_n + m_16O)  [AMU]
    double m_n  = 1.00866;   // neutron mass AMU
    double m_16O = 15.9949;  // 16O mass AMU
    double mu_AMU = m_n * m_16O / (m_n + m_16O);

    // Potential (matching Raphael defaults from dwba_zr.py)
    // r0=1.10, a=0.65, Vso=-6, rso=1.25, aso=0.65, rc=1.30
    // Note: V (depth) found by bisection — set dummy here
    ChannelPotential pot;
    pot.V    = -60.0;   // initial guess (will be overwritten by bisection)
    pot.R0   = 1.10;
    pot.A    = 0.65;
    pot.VI   = 0.0;
    pot.RI0  = 0.0;
    pot.AI   = 0.0;
    pot.VSO  = -6.0;
    pot.RSO0 = 1.25;
    pot.ASO  = 0.65;
    pot.RC0  = 1.30;

    // Grid
    int    NSTEPS  = 3000;
    double STEPSIZE = 0.01;   // 0.01 fm → 30 fm total

    struct OrbitalSpec {
        std::string label;
        std::string tag;
        int    n;    // nodes (radial)
        int    l;
        double j;
        double BE;   // binding energy (positive MeV)
    };

    std::vector<OrbitalSpec> orbs = {
        {"0d5/2 (gs, Ex=0)",      "d52", 0, 2, 2.5, 4.1431},
        {"1s1/2 (Ex=0.871 MeV)", "s12", 1, 0, 0.5, 3.2721},
    };

    DWBA dwba;

    for (auto& orb : orbs) {
        std::cout << "\n=== C++: " << orb.label << " ===" << std::endl;

        // Build Channel
        Channel ch;
        ch.Pot        = pot;
        ch.mu         = mu_AMU;
        ch.Projectile.Z = 0;   // neutron
        ch.Projectile.A = 1;
        ch.Target.Z   = 8;     // 16O
        ch.Target.A   = 16;
        ch.StepSize   = STEPSIZE;
        ch.MaxR       = NSTEPS * STEPSIZE;
        ch.NSteps     = NSTEPS;
        ch.RGrid.resize(NSTEPS);
        for (int i = 0; i < NSTEPS; i++) ch.RGrid[i] = (i+1)*STEPSIZE;
        ch.WaveFunction.resize(NSTEPS, {0,0});
        ch.V_real.resize(NSTEPS, 0);
        ch.V_imag.resize(NSTEPS, 0);
        ch.V_so_real.resize(NSTEPS, 0);
        ch.V_so_imag.resize(NSTEPS, 0);
        ch.V_coulomb.resize(NSTEPS, 0);

        dwba.CalcBoundState(ch, orb.n, orb.l, orb.j, orb.BE);

        // Write phi = u/r to file
        std::string fname = std::string("/tmp/cpp_wf_") + orb.tag + ".txt";
        std::ofstream ofs(fname);
        ofs << "# r(fm)  phi=u/r\n";
        for (int i = 0; i < NSTEPS; i++) {
            double r = (i+1) * STEPSIZE;
            double phi = ch.WaveFunction[i].real();
            ofs << std::setw(10) << std::fixed << std::setprecision(4) << r
                << "  " << std::setw(14) << std::scientific << std::setprecision(6) << phi << "\n";
        }
        ofs.close();
        std::cout << "  Saved: " << fname << std::endl;

        // Normalization check: integral phi^2 r^2 dr = 1
        double norm = 0.0;
        for (int i = 0; i < NSTEPS; i++) {
            double r   = (i+1) * STEPSIZE;
            double phi = ch.WaveFunction[i].real();
            norm += phi*phi * r*r * STEPSIZE;
        }
        std::cout << "  Norm integral phi^2*r^2*dr = " << norm << " (should be 1)" << std::endl;
        std::cout << "  V_sol = " << ch.Pot.V << " MeV" << std::endl;
    }

    return 0;
}
