// test_prjbs.cpp — Compare C++ projectile bound state (n in deuteron) vs Raphael
// n=0, l=0, j=0.5, BE=2.224 MeV (1s1/2), WS R0=1.00 A=0.50 Vso=0 RC=0
// Core = proton (Z=1, A=1), particle = neutron (Z=0, A=1)

#include "dwba.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

int main() {
    // Build a minimal Channel for n+p bound state
    Channel ch;
    ch.Projectile.Z = 0;  // neutron
    ch.Projectile.A = 1;
    ch.Projectile.Mass = 939.565;  // MeV
    ch.Target.Z = 1;               // proton (core)
    ch.Target.A = 1;
    ch.Target.Mass = 938.272;
    
    // Reduced mass (AMU): mu = mn*mp / (mn+mp) / 931.494
    double mn = 939.565, mp = 938.272, amu = 931.494;
    ch.mu = mn * mp / (mn + mp) / amu;  // ≈ 0.4997 AMU

    // Grid: same as Raphael (0 to 29.9 fm, h=0.1, 300 pts)
    ch.StepSize = 0.1;
    ch.NSteps   = 300;
    ch.MaxR     = 30.0;
    ch.RGrid.resize(300);
    for (int i = 0; i < 300; i++) ch.RGrid[i] = i * 0.1;
    ch.WaveFunction.resize(300);
    ch.V_real.resize(300); ch.V_imag.resize(300);
    ch.V_so_real.resize(300); ch.V_so_imag.resize(300);
    ch.V_coulomb.resize(300);

    // Potential: WS V0 to be fitted, R0=1.00 A=0.50 Vso=0 RC=0
    ch.Pot.R0   = 1.00;
    ch.Pot.A    = 0.50;
    ch.Pot.VSO  = 0.0;
    ch.Pot.RSO0 = 1.00;
    ch.Pot.ASO  = 0.50;
    ch.Pot.RC0  = 0.0;

    // Solve bound state: n=0 (ground state, 0 nodes), l=0, j=0.5, BE=2.224 MeV
    DWBA dwba;
    dwba.CalculateBoundState(ch, 0, 0, 0.5, 2.224);

    // Print wavefunction phi = u/r
    std::cout << "\n=== Projectile BS (n in deuteron, 1s1/2, BE=2.224 MeV) ===\n";
    std::cout << std::setw(8) << "r(fm)"
              << std::setw(16) << "phi=u/r (C++)"
              << std::endl;
    double norm = 0.0;
    for (int i = 0; i < 300; i++) {
        double r = i * 0.1;
        double phi = ch.WaveFunction[i].real();
        norm += phi * phi * 0.1;  // ∫phi²dr (trapez approx)
        if (i % 10 == 0 || i < 5)
            std::cout << std::setw(8) << std::fixed << std::setprecision(2) << r
                      << std::setw(16) << std::scientific << std::setprecision(6) << phi
                      << "\n";
    }
    std::cout << "\nNorm ∫phi²dr = " << norm << " (should be ~1.0)\n";
    std::cout << "V_sol = " << ch.Pot.V << " MeV\n";

    // Compute kappa for comparison
    double BE = 2.224;
    double hbarc = 197.32697;
    double mu_MeV = ch.mu * amu;
    double kappa = std::sqrt(2.0 * mu_MeV * BE) / hbarc;
    std::cout << "kappa = " << kappa << " fm^-1\n";
    std::cout << "Expected tail: phi ~ ANC * exp(-kappa*r) / r  at large r\n";
    
    // Print tail (r=8..15 fm)
    std::cout << "\nTail comparison:\n";
    std::cout << std::setw(8) << "r" << std::setw(16) << "phi(C++)"
              << std::setw(16) << "exp(-kappa*r)/r" << std::setw(12) << "ratio\n";
    for (int i = 80; i <= 150; i += 5) {
        double r = i * 0.1;
        double phi = ch.WaveFunction[i].real();
        double tail = std::exp(-kappa*r) / r;
        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << r
                  << std::setw(16) << std::scientific << std::setprecision(6) << phi
                  << std::setw(16) << tail
                  << std::setw(12) << std::fixed << std::setprecision(6) << phi/tail
                  << "\n";
    }

    // Save to file for comparison with Raphael
    std::ofstream fout("prjbs_cpp.txt");
    fout << "# r(fm)  phi=u/r\n";
    for (int i = 0; i < 300; i++)
        fout << std::fixed << std::setprecision(4) << (i*0.1) << "  "
             << std::scientific << std::setprecision(8) << ch.WaveFunction[i].real() << "\n";
    std::cout << "\nSaved to prjbs_cpp.txt\n";
    return 0;
}
