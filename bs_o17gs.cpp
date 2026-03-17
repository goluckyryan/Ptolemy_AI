// bs_o17gs.cpp — Standalone 17O g.s. (0d5/2) bound state calculator
// Matches Ptolemy TARGET section: n=0, l=2, j=5/2, BE=4.143 MeV
// WS: V=60.0 (initial), R0=1.10, A=0.65, VSO=0, RC0=1.30
// Uses DWBA::CalcBoundState()

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "include/dwba.h"

int main() {
    const double AMU = 931.494061;

    // Neutron (Z=0, A=1) + 16O core (Z=8, A=16)
    double m_n  = 1.00866;
    double m_16O = 15.9949;
    double mu_AMU = m_n * m_16O / (m_n + m_16O);

    // Potential matching Ptolemy TARGET section (NO spin-orbit)
    ChannelPotential pot;
    pot.V    = -60.0;   // initial V (Ptolemy will bisect to find correct V)
    pot.R0   = 1.10;
    pot.A    = 0.65;
    pot.VI   = 0.0;
    pot.RI0  = 0.0;
    pot.AI   = 0.0;
    pot.VSO  = 0.0;    // NO spin-orbit in Ptolemy TARGET section
    pot.RSO0 = 1.25;
    pot.ASO  = 0.65;
    pot.VSOI = 0.0;
    pot.RSOI0 = 0.0;
    pot.ASOI = 0.0;
    pot.VSI  = 0.0;
    pot.RSI0 = 0.0;
    pot.ASI  = 0.0;
    pot.RC0  = 1.30;

    // Grid: 0.08125 fm step → 30 fm (matches Ptolemy step)
    double h = 0.08125;
    int    NSTEPS = static_cast<int>(30.0 / h) + 5;  // ~375 steps

    Channel ch;
    ch.Pot = pot;
    ch.mu  = mu_AMU;
    ch.Projectile.Z = 0;   // neutron
    ch.Projectile.A = 1;
    ch.Target.Z   = 8;     // 16O core
    ch.Target.A   = 16;
    ch.StepSize   = h;
    ch.MaxR       = NSTEPS * h;
    ch.NSteps     = NSTEPS;
    ch.RGrid.resize(NSTEPS);
    for (int i = 0; i < NSTEPS; i++) ch.RGrid[i] = (i+1)*h;
    ch.WaveFunction.resize(NSTEPS, {0,0});
    ch.V_real.resize(NSTEPS, 0);
    ch.V_imag.resize(NSTEPS, 0);
    ch.V_so_real.resize(NSTEPS, 0);
    ch.V_so_imag.resize(NSTEPS, 0);
    ch.V_coulomb.resize(NSTEPS, 0);

    DWBA dwba;
    dwba.CalcBoundState(ch, 0, 2, 2.5, 4.143);

    // Print r and phi=u/r (matching Ptolemy output format)
    // Ptolemy prints: step r phi
    // We print: r phi (for easy comparison)
    // phi_vec[i] is at r=i*h (0-indexed: i=0→r=0, i=1→r=h, ...)
    // ch.WaveFunction[i] = phi_vec[i], so print r=i*h for WaveFunction[i]
    // (Note: ch.RGrid[i] = (i+1)*h, so WaveFunction[i] is at r=i*h, NOT (i+1)*h)
    for (int i = 0; i < NSTEPS; i++) {
        double r = i * h;   // correct r: phi_vec[i] is at r=i*h
        double phi = ch.WaveFunction[i].real();
        if (r > 30.0) break;
        if (r < 1e-10) continue;  // skip r=0
        std::cout << std::fixed << std::setprecision(5) << r
                  << "  " << std::scientific << std::setprecision(6) << phi << "\n";
    }

    // Print summary to stderr
    double maxAmp = 0.0; double rPeak = 0.0;
    double norm = 0.0;
    for (int i = 0; i < NSTEPS; i++) {
        double r = (i+1) * h;
        double phi = ch.WaveFunction[i].real();
        norm += phi*phi * r*r * h;
        if (std::abs(phi) > maxAmp) { maxAmp = std::abs(phi); rPeak = r; }
    }
    std::cerr << "C++ V_sol=" << ch.Pot.V << " MeV\n";
    std::cerr << "C++ maxAmp=" << maxAmp << " at r=" << rPeak << " fm\n";
    std::cerr << "C++ norm=" << norm << "\n";

    return 0;
}
