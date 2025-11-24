#include <iostream>
#include "potential_eval.h"
#include <cmath>

// Woods-Saxon form factor
double WoodsSaxon(double r, double V, double r0, double a, double A_target) {
    if (V == 0.0) return 0.0;
    double R = r0 * std::pow(A_target, 1.0/3.0);
    return V / (1.0 + std::exp((r - R) / a));
}

// Derivative of Woods-Saxon form factor (for Spin-Orbit)
// d/dr (1 / (1 + exp((r-R)/a))) = - (1/a) * exp(...) / (1 + exp(...))^2
double DerivWoodsSaxon(double r, double V, double r0, double a, double A_target) {
    if (V == 0.0) return 0.0;
    double R = r0 * std::pow(A_target, 1.0/3.0);
    double ex = std::exp((r - R) / a);
    double denom = 1.0 + ex;
    return - (V / a) * ex / (denom * denom);
}

void EvaluatePotential(double r, const ChannelPotential& pot, double& V_real, double& V_imag, 
                       double& V_so_real, double& V_so_imag, double& V_coulomb, 
                       double Z1, double Z2, double A_target) {
    
    // Real Central
    V_real = WoodsSaxon(r, pot.V, pot.R0, pot.A, A_target);
    
    // Imaginary Central (Volume)
    V_imag = WoodsSaxon(r, pot.VI, pot.RI0, pot.AI, A_target);
    
    // Imaginary Surface (Derivative Woods-Saxon)
    // 4 * a * d/dr (WS) ? Or just 4 * W_D * exp / (1+exp)^2
    // Ptolemy manual says: 4 * ASI * (D/DR) 1/(1+XSI) * VSI
    // My DerivWoodsSaxon returns V * d/dr (form).
    // So if I pass VSI, I get VSI * d/dr.
    // But Ptolemy formula has 4 * ASI.
    // Let's check DerivWoodsSaxon implementation.
    // It returns - (V/a) * ...
    // So 4 * a * Deriv = -4 * V * ...
    // Standard surface absorption is usually positive W_D * ...
    // Ptolemy: VSI * 4 * ASI * (D/DR) ...
    // Let's implement Surface term explicitly.
    
    if (pot.VSI != 0.0) {
        double R_si = pot.RSI0 * std::pow(A_target, 1.0/3.0);
        double ex_si = std::exp((r - R_si) / pot.ASI);
        double denom_si = 1.0 + ex_si;
        double surface_term = 4.0 * pot.VSI * ex_si / (denom_si * denom_si); // This is positive
        V_imag += surface_term;
    }
    
    // Spin-Orbit (Real)
    // (VSO) * (1/r) * d/dr (WS) ?
    // Ptolemy: (VSO) * 4 * L.S * (1/R) * (D/DR) 1/(1+XSO) ?
    // Wait, the L.S factor is applied later in the radial equation.
    // Here we just return the radial part V_so(r).
    // Standard Thomas form: V_so * (hbar/m_pi c)^2 * (1/r) * d/dr (WS)
    // Ptolemy manual might have specific factors.
    // " (VSO+TAU*V) * 4*L.S * (1/R) * (D/DR) 1/(1+XSO) "
    // So the factor 4 is there.
    // And 1/r is there.
    // I will return the term without L.S.
    // V_so(r) = VSO * 4 * (1/r) * (D/DR) ...
    // But wait, D/DR of 1/(1+X) is -1/a * ...
    // So it's -4 * VSO / (r * a) * ...
    
    if (pot.VSO != 0.0 && r > 0.001) {
        double R_so = pot.RSO0 * std::pow(A_target, 1.0/3.0);
        double ex_so = std::exp((r - R_so) / pot.ASO);
        double denom_so = 1.0 + ex_so;
        double deriv = - (1.0 / pot.ASO) * ex_so / (denom_so * denom_so);
        V_so_real = pot.VSO * 4.0 * (1.0/r) * deriv; // Note: L.S factor to be applied later
        // Also check if factor 2 or 4 is used. Ptolemy says 4.
    } else {
        V_so_real = 0.0;
    }
    
    // Spin-Orbit (Imaginary)
    if (pot.VSOI != 0.0 && r > 0.001) {
        double R_soi = pot.RSOI0 * std::pow(A_target, 1.0/3.0);
        double ex_soi = std::exp((r - R_soi) / pot.ASOI);
        double denom_soi = 1.0 + ex_soi;
        double deriv = - (1.0 / pot.ASOI) * ex_soi / (denom_soi * denom_soi);
        V_so_imag = pot.VSOI * 4.0 * (1.0/r) * deriv;
    } else {
        V_so_imag = 0.0;
    }
    
    // Coulomb
    // Uniformly charged sphere of radius Rc
    double Rc = pot.RC0 * std::pow(A_target, 1.0/3.0);
    const double e2 = 1.439965; // MeV fm (e^2 / 4pi eps0 ?)
    // Wait, HBARC / 137.036 = 197.327 / 137.036 = 1.43996
    double coul_const = Z1 * Z2 * 1.439965;
    
    if (r >= Rc) {
        V_coulomb = coul_const / r;
    } else {
        V_coulomb = (coul_const / (2.0 * Rc)) * (3.0 - (r * r) / (Rc * Rc));
    }
}
