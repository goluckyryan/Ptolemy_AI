#include <iostream>
#include "potential_eval.h"
#include <cmath>
#include <algorithm>

// Compute Ptolemy R0MASS factor
// R0MASS = A_heavy^(1/3) + A_light^(1/3) when A_light > 2.5 (heavy-ion convention)
// R0MASS = A_heavy^(1/3)                  when A_light <= 2.5
// For alpha (Ap=4>2.5): R0MASS = At^(1/3)+Ap^(1/3) — matches Fortran elastic S-matrix
static double computeR0MASS(double A_target, double A_projectile) {
    double A_heavy = std::max(A_target, A_projectile);
    double A_light = std::min(A_target, A_projectile);
    if (A_light > 2.5) {
        return std::cbrt(A_heavy) + std::cbrt(A_light);
    } else {
        return std::cbrt(A_heavy);
    }
}

// Woods-Saxon form factor
double WoodsSaxon(double r, double V, double r0, double a, double R0MASS) {
    if (V == 0.0) return 0.0;
    double R = r0 * R0MASS;
    return V / (1.0 + std::exp((r - R) / a));
}

// Derivative of Woods-Saxon form factor (for Spin-Orbit)
// d/dr (1 / (1 + exp((r-R)/a))) = - (1/a) * exp(...) / (1 + exp(...))^2
double DerivWoodsSaxon(double r, double V, double r0, double a, double R0MASS) {
    if (V == 0.0) return 0.0;
    double R = r0 * R0MASS;
    double ex = std::exp((r - R) / a);
    double denom = 1.0 + ex;
    return - (V / a) * ex / (denom * denom);
}

void EvaluatePotential(double r, const ChannelPotential& pot, double& V_real, double& V_imag, 
                       double& V_so_real, double& V_so_imag, double& V_coulomb, 
                       double Z1, double Z2, double A_target, double A_projectile) {
    
    // Compute R0MASS once — used for ALL potential terms
    double R0MASS = computeR0MASS(A_target, A_projectile);
    
    // Real Central
    V_real = WoodsSaxon(r, pot.V, pot.R0, pot.A, R0MASS);
    
    // Imaginary Central (Volume)
    V_imag = WoodsSaxon(r, pot.VI, pot.RI0, pot.AI, R0MASS);
    
    // Imaginary Surface (Derivative Woods-Saxon)
    // Standard form: 4 * W_D * exp / (1+exp)^2
    if (pot.VSI != 0.0) {
        double R_si = pot.RSI0 * R0MASS;
        double ex_si = std::exp((r - R_si) / pot.ASI);
        double denom_si = 1.0 + ex_si;
        double surface_term = 4.0 * pot.VSI * ex_si / (denom_si * denom_si);
        V_imag += surface_term;
    }
    
    // Spin-Orbit (Real)
    // Ptolemy WOODSX type-2 uses: 2*VSO * (1/r) * d/dr[WS]
    // In WAVELJ: SDOTL * LSOR = SDOTL * (-H²/12E) * WOODSX(type=2)
    // SDOTL = 2*<L·S> for proton (JSPS=1). Physical: 4*<L·S>*VSO/r*d(WS)/dr.
    // So SDOTL * WOODSX = 2*<L·S> * 2*VSO/r*d(WS)/dr = 4*<L·S>*VSO/r*d(WS)/dr. ✓
    // Ptolemy convention: V_SO(r) = VSO * (1/r) * d(WS)/dr  (Thomas form, negative at surface)
    // In f(r): f_so = L·S * f_conv * V_SO(r)  where L·S = SDOTL computed in wavelj.cpp
    if (pot.VSO != 0.0 && r > 0.001) {
        double R_so = pot.RSO0 * R0MASS;
        double ex_so = std::exp((r - R_so) / pot.ASO);
        double denom_so = 1.0 + ex_so;
        double deriv = - (1.0 / pot.ASO) * ex_so / (denom_so * denom_so);
        V_so_real = pot.VSO * 2.0 * (1.0/r) * deriv;  // 2*VSO*Thomas(r), matches WOODSX ITYPE=2
    } else {
        V_so_real = 0.0;
    }
    
    // Spin-Orbit (Imaginary) — same factor 2 convention
    if (pot.VSOI != 0.0 && r > 0.001) {
        double R_soi = pot.RSOI0 * R0MASS;
        double ex_soi = std::exp((r - R_soi) / pot.ASOI);
        double denom_soi = 1.0 + ex_soi;
        double deriv = - (1.0 / pot.ASOI) * ex_soi / (denom_soi * denom_soi);
        V_so_imag = pot.VSOI * 2.0 * (1.0/r) * deriv;  // 2*VSOI*Thomas(r), matches WOODSX ITYPE=2
    } else {
        V_so_imag = 0.0;
    }
    
    // Coulomb — Uniformly charged sphere of radius Rc
    double Rc = pot.RC0 * R0MASS;
    const double e2 = 1.439965; // MeV fm (e^2 / 4pi eps0)
    double coul_const = Z1 * Z2 * 1.439965;
    
    if (r >= Rc) {
        V_coulomb = coul_const / r;
    } else {
        V_coulomb = (coul_const / (2.0 * Rc)) * (3.0 - (r * r) / (Rc * Rc));
    }
}
