// reid_ivphi_compare.cpp
// Print current C++ IVPHI_P at sample rP values vs Fortran reference
// Compile and link against existing DWBA infrastructure
// Build: see below

#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include "include/dwba.h"

int main() {
    DWBA dwba;

    // Setup: match fr_o16dp parameters exactly
    // Incoming: d + 16O
    Channel Incoming;
    Incoming.Projectile = {1, 2, 2.0136};  // deuteron
    Incoming.Target     = {8, 16, 15.9949}; // 16O
    Incoming.Elab = 20.0;
    Incoming.Pot.V   = 88.9546; Incoming.Pot.R0  = 1.1489; Incoming.Pot.A   = 0.7508;
    Incoming.Pot.W   = 2.3480;  Incoming.Pot.RI0 = 1.3446; Incoming.Pot.AI  = 0.6030;
    Incoming.Pot.VSO = 3.5570;  Incoming.Pot.RSO0= 0.9720; Incoming.Pot.ASO = 1.0110;
    Incoming.Pot.RC0 = 1.3;
    Incoming.JSPS = 2;  // deuteron spin=1 → 2S=2
    dwba.SetupChannel(Incoming);
    dwba.WavSet(Incoming);

    // Outgoing: p + 17O
    Channel Outgoing;
    Outgoing.Projectile = {1, 1, 1.00728}; // proton
    Outgoing.Target     = {8, 17, 16.9991}; // 17O
    Outgoing.Pot.V   = 49.5434; Outgoing.Pot.R0  = 1.1462; Outgoing.Pot.A   = 0.6753;
    Outgoing.Pot.W   = 2.0611;  Outgoing.Pot.RI0 = 1.1462; Outgoing.Pot.AI  = 0.6753;
    Outgoing.Pot.VSO = 5.2956;  Outgoing.Pot.RSO0= 0.9338; Outgoing.Pot.ASO = 0.5900;
    Outgoing.Pot.RC0 = 1.3;
    Outgoing.JSPS = 1;  // proton spin=1/2 → 2S=1
    Outgoing.Elab = 0.0;
    dwba.SetupChannel(Outgoing);
    dwba.WavSet(Outgoing);

    // Projectile bound state: Reid deuteron (L=0, S-state)
    BoundState ProjectileBS;
    ProjectileBS.l = 0; ProjectileBS.j = 1.0; ProjectileBS.n = 0;
    ProjectileBS.BindingEnergy = 2.2247;
    ProjectileBS.Pot.V  = 1.0; ProjectileBS.Pot.R0 = 1.0; ProjectileBS.Pot.A = 0.5;
    ProjectileBS.Pot.RC0 = 1.3;

    // Build PrjBS_ch exactly as InelDc does
    Channel PrjBS_ch;
    PrjBS_ch.Pot    = ProjectileBS.Pot;
    PrjBS_ch.Target = Outgoing.Projectile;
    PrjBS_ch.Projectile.Z    = 1;  // neutron (d - p)
    PrjBS_ch.Projectile.A    = 1;
    PrjBS_ch.Projectile.Mass = 1.00866;
    double mb = Outgoing.Projectile.Mass;
    double mx = PrjBS_ch.Projectile.Mass;
    const double AMU_MEV = 931.494061;
    PrjBS_ch.mu   = mb * mx / (mb + mx) / AMU_MEV;
    PrjBS_ch.StepSize = Incoming.StepSize;
    PrjBS_ch.MaxR     = Incoming.MaxR;
    PrjBS_ch.NSteps   = Incoming.NSteps;
    PrjBS_ch.RGrid    = Incoming.RGrid;
    PrjBS_ch.WaveFunction.resize(Incoming.NSteps, {0,0});
    PrjBS_ch.V_real.resize(Incoming.NSteps, 0);
    PrjBS_ch.V_imag.resize(Incoming.NSteps, 0);
    PrjBS_ch.V_so_real.resize(Incoming.NSteps, 0);
    PrjBS_ch.V_so_imag.resize(Incoming.NSteps, 0);
    PrjBS_ch.V_coulomb.resize(Incoming.NSteps, 0);
    dwba.WavSet(PrjBS_ch);

    // Load current projectile WF (AV18 or WS fallback)
    bool loaded = dwba.LoadDeuteronWavefunction(PrjBS_ch,
        "/home/node/working/ptolemy_2019/Cpp_AI/data", "av18-phi-v");
    if (!loaded) {
        printf("# AV18 not loaded — using WS CalculateBoundState\n");
        dwba.CalculateBoundState(PrjBS_ch, ProjectileBS.n, ProjectileBS.l,
                                 ProjectileBS.j, ProjectileBS.BindingEnergy);
    } else {
        printf("# AV18 loaded\n");
    }

    // Rebuild V_real using WS potential (as InelDc does)
    for (int i = 0; i < PrjBS_ch.NSteps; ++i) {
        double r = i * PrjBS_ch.StepSize;
        if (r < 0.001) { PrjBS_ch.V_real[i] = 0.0; continue; }
        double vi, vso, vsoi, vc;
        dwba.EvaluatePotential(r, PrjBS_ch.Pot, PrjBS_ch.V_real[i], vi,
                               vso, vsoi, vc,
                               PrjBS_ch.Projectile.Z, PrjBS_ch.Target.Z,
                               PrjBS_ch.Target.A, PrjBS_ch.Projectile.A);
    }

    // Build IVPHI_P = phi_P * V_real
    int N = PrjBS_ch.NSteps;
    double h = PrjBS_ch.StepSize;
    std::vector<double> IVPHI_P(N, 0.0);
    for (int i = 1; i < N; ++i)
        IVPHI_P[i] = std::abs(PrjBS_ch.WaveFunction[i].real()) * PrjBS_ch.V_real[i];

    // Print comparison at sample points
    printf("\n# C++ vs Fortran IVPHI_P comparison\n");
    printf("# Fortran (Reid cubic spline + L2 norm, V_Reid):\n");
    printf("#   rP=0.5: -174.1  rP=1.0: 45.87  rP=1.5: 17.49  rP=2.0: 3.697\n");
    printf("#   rP=2.5: 0.844   rP=3.0: 0.245  rP=3.5: 0.090  rP=4.0: 0.039\n");
    printf("\n");
    printf("# %6s  %14s  %14s  %8s\n", "rP(fm)", "C++_IVPHI_P", "Fortran_ref", "Ratio");

    double fortran_ref[] = {-174.148, 45.868, 17.487, 3.697, 0.844, 0.245, 0.090, 0.039};
    double rP_vals[]     = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};

    for (int k = 0; k < 8; ++k) {
        double rP = rP_vals[k];
        int idx = (int)std::round(rP / h);
        if (idx >= N) { printf("  %6.2f  OUT_OF_RANGE\n", rP); continue; }
        double cpp_val = IVPHI_P[idx];
        double fort = fortran_ref[k];
        double ratio = (std::abs(fort) > 1e-8) ? cpp_val / fort : 0.0;
        printf("  %6.2f  %14.5e  %14.5e  %8.4f\n", rP, cpp_val, fort, ratio);
    }

    // Also print phi_P at sample points for extra diagnostic
    printf("\n# phi_P comparison:\n");
    printf("# %6s  %14s\n", "rP(fm)", "C++_phi_P");
    for (int k = 0; k < 8; ++k) {
        double rP = rP_vals[k];
        int idx = (int)std::round(rP / h);
        if (idx >= N) continue;
        double phi = std::abs(PrjBS_ch.WaveFunction[idx].real());
        printf("  %6.2f  %14.5e\n", rP, phi);
    }
    return 0;
}
