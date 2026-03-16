// elastic_sm148aa.cpp — Standalone elastic scattering for 148Sm(α,α) at 50 MeV
// Uses the existing WavSet/WavElj machinery from the DWBA project.
//
// Computes:
//   - S-matrix Re(S_L), Im(S_L) for L = 0..90
//   - Elastic DCS dσ/dΩ (mb/sr) for θ = 1°..180°
//
// Output files:
//   sm148_aa_smat_cpp.txt   — columns: L  Re(S)  Im(S)
//   sm148_aa_xsec_cpp.txt   — columns: theta_deg  dcs_mb_sr

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <string>

// Include project headers
#include "include/dwba.h"
#include "include/Potentials.h"
#include "include/Isotope.h"
#include "include/rcwfn.h"
#include "include/potential_eval.h"

// ============================================================
// Physical constants (matching the rest of the project)
// ============================================================
static const double AMU_MEV  = 931.494061;   // MeV/c² per amu
static const double HBARC    = 197.32697;    // MeV fm
static const double FINE     = 1.0 / 137.036; // fine-structure constant
static const double PI       = 3.14159265358979323846;

// ============================================================
// Compute elastic DCS (mb/sr) for one angle using partial-wave sum
// including full Coulomb amplitude
//
// f(θ) = f_C(θ) + f_N(θ)
//   f_C = -η / (2k sin²(θ/2)) * exp(2iσ₀ - 2iη ln(sin(θ/2))) * Γ(1+iη)/|Γ(1+iη)|
//   f_N = (1/2ik) Σ_L (2L+1) exp(2iσ_L) (S_L - 1) P_L(cosθ)
//
// For speed we compute the Coulomb amplitude using the exact formula.
// ============================================================

// Legendre polynomial via recurrence
static double LegendreP(int L, double x) {
    if (L == 0) return 1.0;
    if (L == 1) return x;
    double pm2 = 1.0, pm1 = x, p = 0.0;
    for (int l = 2; l <= L; ++l) {
        p = ((2*l - 1) * x * pm1 - (l - 1) * pm2) / l;
        pm2 = pm1;
        pm1 = p;
    }
    return p;
}

// Coulomb phase shift σ_L = arg Γ(L+1+iη) using Stirling / recurrence
static double CoulombPhase(int L, double eta) {
    // σ_L = σ_0 + Σ_{l=1}^{L} arctan(η/l)
    // σ_0 = arg Γ(1+iη)
    // Use series: arg Γ(1+iη) = -γ*η + Σ_{n=1}^{∞} [η/n - arctan(η/n)]
    // For our purpose, compute via recurrence from the Coulomb functions.
    // Simplest: sigma_L = sum_{n=1}^{L} atan(eta/n)  +  sigma_0
    // sigma_0: use the known series (Abramowitz 6.1.27)
    // sigma_0 = -Euler*eta + sum_{n=1}^inf [eta/n - atan(eta/n)]
    // But this converges slowly. Instead use lgamma:
    // sigma_L = Im[ ln Γ(L+1+iη) ]
    // We compute via the recurrence:  σ_L = σ_{L-1} + atan(eta/L)
    // Need σ_0 = Im ln Γ(1+iη)
    // Use the Lanczos approximation or the series relation.

    // Fast computation using lgamma_r or the recurrence-only approach:
    // σ_L = atan(eta/1) + atan(eta/2) + ... + atan(eta/L) + sigma_0
    // sigma_0 computed via arg(Γ(1+iη)) - use spouge / Stirling for large eta
    // For small eta: sigma_0 ≈ -gamma_euler * eta + eta^3/3*(1 - 1/4 + 1/9 - ...)
    // For our case eta ~ 5, just use the partial sum approach with enough terms.

    // Use the relation: Im[ln Γ(z)] computed by the Stirling series for |z| large
    // z = L+1+iη; for L >> eta this is accurate.
    // For small L+1, accumulate arctan terms from a reference sigma_0.

    // Practical: compute sigma_0 using the accurate formula from DLMF 5.9.3
    // log Γ(1+iη) = -γη + Σ_{k=2}^∞ (-1)^k ζ(k)/k * (iη)^k
    // Im part:  sigma_0 = -γη - η³/3*ζ(3) + ... convergence ok for eta<5
    // OR: use the product formula: Γ(1+iη) = iη * Γ(iη)
    // Γ(iη) = π/(η*sinh(πη)*|Γ(iη)|²)^0.5 ... messy.

    // EASIEST: compute sigma_0 via rcwfn (Coulomb wave functions at some rho)
    // Actually we'll just use the recurrence from sigma_0 = 0 for η=0 case,
    // but we have η ≠ 0 here.

    // Use the accurate simple recurrence: sigma_L = Im[ ln Γ(L+1+iη) ]
    // Compute by accumulating: sigma_L = sigma_0 + sum_{n=1}^{L} atan(eta/n)
    // sigma_0 from: Im[ln Γ(1+iη)] using the asymptotic series starting at large L0
    // and subtracting the accumulated sum.

    // Implementation: compute numerically with sufficient L0
    static const double euler_gamma = 0.5772156649015328606;
    // Lanczos approx for Im[ln Gamma(1+ieta)]:
    // use the Stirling series for Gamma(N+1+ieta) for large N, then subtract atan sum
    // sigma_0 = Im[ln Γ(N+1+iη)] - Σ_{n=1}^{N} atan(η/n)
    // Im[ln Γ(z)] for large |z|: Im[ln Γ(z)] ≈ (z-0.5)*arg(z) - Im(z)*ln|z|
    //   + Im[ -z + 0.5*ln(2π) + 1/(12z) - 1/(360z³) + ... ]
    // z = N+1+iη: arg(z) = atan(η/(N+1)), |z| = sqrt((N+1)²+η²)
    int N0 = 200;  // large enough
    double z_re = N0 + 1.0;
    double z_im = eta;
    double z2 = z_re*z_re + z_im*z_im;
    double ln_mod_z = 0.5 * std::log(z2);
    double arg_z    = std::atan2(z_im, z_re);
    // Im[ln Γ(z)] ≈ (z_re - 0.5)*arg_z - z_im*ln|z| - z_im
    //   + Im[1/(12z)] - Im[1/(360z³)] + ...
    double im_lngamma_N0 = (z_re - 0.5) * arg_z - z_im * ln_mod_z - z_im
                           - z_im / (12.0 * z2);  // first Stirling correction
    // Subtract sum_{n=1}^{N0} atan(eta/n)
    double sum_atan = 0.0;
    for (int n = 1; n <= N0; ++n)
        sum_atan += std::atan(eta / n);
    double sigma_0 = im_lngamma_N0 - sum_atan;

    // Now compute sigma_L
    double sigma_L = sigma_0;
    for (int n = 1; n <= L; ++n)
        sigma_L += std::atan(eta / n);
    return sigma_L;
}

// Elastic DCS in fm² using partial-wave Coulomb + nuclear amplitudes
static double ElasticDCS(double theta_deg,
                          double k, double eta,
                          const std::vector<std::complex<double>>& SMatrix,
                          const std::vector<double>& CoulPhase) {
    if (theta_deg < 0.01) return std::numeric_limits<double>::quiet_NaN();

    double theta = theta_deg * PI / 180.0;
    double cos_th = std::cos(theta);
    double sin_th2 = std::sin(theta / 2.0);  // sin(θ/2)

    int Lmax = (int)SMatrix.size() - 1;

    // Coulomb amplitude: f_C(θ) = -η/(2k sin²(θ/2)) * exp(i*phi_C)
    // phi_C = 2σ_0 - 2η*ln(sin(θ/2))
    double sigma0 = CoulPhase.size() > 0 ? CoulPhase[0] : CoulombPhase(0, eta);
    double phi_C = 2.0 * sigma0 - 2.0 * eta * std::log(sin_th2);
    double fc_mod = eta / (2.0 * k * sin_th2 * sin_th2);
    std::complex<double> f_C(-fc_mod * std::cos(phi_C), -fc_mod * std::sin(phi_C));
    // Note: standard form is -η/(2k sin²(θ/2)) exp(iφ_C) (negative sign)

    // Nuclear amplitude: f_N = (1/2ik) Σ_L (2L+1) e^{2iσ_L} (S_L - 1) P_L(cosθ)
    std::complex<double> f_N(0.0, 0.0);
    std::complex<double> two_ik(0.0, 2.0 * k);

    for (int L = 0; L <= Lmax; ++L) {
        double sigma_L = CoulPhase[L];
        std::complex<double> e2sig(std::cos(2.0 * sigma_L), std::sin(2.0 * sigma_L));
        std::complex<double> SL = SMatrix[L];
        double PL = LegendreP(L, cos_th);
        f_N += (double)(2 * L + 1) * e2sig * (SL - 1.0) * PL;
    }
    f_N /= two_ik;

    std::complex<double> f_tot = f_C + f_N;
    return std::norm(f_tot);  // |f|² in fm²
}

// ============================================================
// main
// ============================================================
int main() {
    std::cout << "=== Elastic scattering 148Sm(α,α) at 50 MeV (C++ code) ===\n\n";

    // --- Set up the incoming channel manually ---
    // 148Sm: Z=62, A=148; alpha: Z=2, A=4
    // Potential: V=65.5, R0=1.427, A=0.671 (volume WS, real)
    //            W=29.8, RI0=1.427, AI=0.671 (volume WS, imaginary)
    //            RC0=1.4 (Coulomb)
    // Ptolemy uses R0MASS = A_heavy^(1/3) + A_light^(1/3) for A_light > 2.5
    // R0MASS = 148^(1/3) + 4^(1/3) = 5.2945 + 1.5874 = 6.8819
    // => R_opt = 1.427 * 6.8819 = 9.8125 fm (Ptolemy says 9.8134 — matches)
    // => R_C   = 1.400 * 6.8819 = 9.6346 fm (Ptolemy says 9.6278 — close)

    // Build channel — non-relativistic kinematics matching Ptolemy/Raphael convention.
    // Ptolemy uses AMU masses (not binding-energy masses) for kinematics at these energies.
    double Ap_kin = 4.0, At_kin = 148.0;   // mass numbers
    double Z1 = 2.0, Z2 = 62.0;

    double Elab = 50.0;  // MeV lab
    double mu    = Ap_kin * At_kin / (Ap_kin + At_kin);   // reduced mass in AMU
    double Ecm   = Elab * At_kin / (Ap_kin + At_kin);     // CM energy (non-rel)
    double k     = std::sqrt(2.0 * mu * AMU_MEV * Ecm) / HBARC;  // fm^-1
    double eta   = Z1 * Z2 * mu * AMU_MEV / (137.036 * HBARC * k);

    std::cout << "Kinematics:\n";
    std::cout << "  Elab = " << Elab << " MeV\n";
    std::cout << "  Ecm  = " << Ecm  << " MeV\n";
    std::cout << "  k    = " << k    << " fm^-1\n";
    std::cout << "  eta  = " << eta  << "\n";
    std::cout << "  mu   = " << mu   << " AMU\n\n";

    // --- Build the Channel struct ---
    Channel ch;
    ch.Projectile.SetIso(4, 2);   // 4He
    ch.Target.SetIso(148, 62);    // 148Sm

    ch.Elab  = Elab;
    ch.Ecm   = Ecm;
    ch.k     = k;
    ch.eta   = eta;
    ch.mu    = mu;
    ch.JSPS  = 0;  // alpha spin = 0 (spin-0 projectile, no spin-orbit)

    // Potential: V=65.5, W=29.8, R0=RI0=1.427, A=AI=0.671, RC0=1.4
    ChannelPotential pot = {};
    pot.V   = 65.5;  pot.R0  = 1.427;  pot.A   = 0.671;
    pot.VI  = 29.8;  pot.RI0 = 1.427;  pot.AI  = 0.671;
    pot.VSO = 0.0;   pot.RSO0= 0.0;    pot.ASO = 0.0;
    pot.VSOI= 0.0;   pot.RSOI0=0.0;    pot.ASOI= 0.0;
    pot.VSI = 0.0;   pot.RSI0= 0.0;    pot.ASI = 0.0;
    pot.RC0 = 1.4;
    ch.Pot  = pot;

    // --- Set up radial grid (same as WavSet in wavelj.cpp) ---
    ch.StepSize = 0.1;
    ch.MaxR     = 30.0;
    ch.NSteps   = static_cast<int>(ch.MaxR / ch.StepSize) + 1;

    ch.RGrid.resize(ch.NSteps);
    ch.V_real.resize(ch.NSteps, 0.0);
    ch.V_imag.resize(ch.NSteps, 0.0);
    ch.V_so_real.resize(ch.NSteps, 0.0);
    ch.V_so_imag.resize(ch.NSteps, 0.0);
    ch.V_coulomb.resize(ch.NSteps, 0.0);

    for (int i = 0; i < ch.NSteps; ++i) {
        double r = i * ch.StepSize;
        ch.RGrid[i] = r;
        if (r < 0.001) continue;
        EvaluatePotential(r, ch.Pot,
                          ch.V_real[i], ch.V_imag[i],
                          ch.V_so_real[i], ch.V_so_imag[i], ch.V_coulomb[i],
                          Z1, Z2, 148.0, 4.0);
    }

    // --- Create a minimal DWBA object to use WavElj ---
    // We cannot call WavElj directly (it's private), so we replicate the
    // Numerov integration + S-matrix extraction inline here.
    // This mirrors wavelj.cpp exactly for spin-0 projectile (JSPS=0, JP=0 for all L).

    int Lmax = 90;
    std::vector<std::complex<double>> SMatrix(Lmax + 1);
    std::vector<double> CoulPhase(Lmax + 1);

    // Pre-compute Coulomb phases
    for (int L = 0; L <= Lmax; ++L)
        CoulPhase[L] = CoulombPhase(L, eta);

    std::cout << "Computing S-matrix for L = 0.." << Lmax << " ...\n";

    const double AMU_MEV2 = 931.494061;
    const double HBARC_L  = 197.32697;
    double f_conv = 2.0 * mu * AMU_MEV2 / (HBARC_L * HBARC_L);

    int N = ch.NSteps;
    double h  = ch.StepSize;
    double h2 = h * h;
    double h2_12 = h2 / 12.0;
    double k2 = k * k;

    const double BIGNUM = 1e30;
    const double STEPI  = 1e-10;

    // Pre-compute Coulomb functions F_L, G_L at both matching points for all L.
    // CRITICAL: call Rcwfn with MINL=0 so it uses the correct downward F recursion
    // and upward G recursion across all L values. Calling with MINL=MAXL=L gives
    // wrong G values because Steed's CF2 at a single L computes Q incorrectly for
    // large L (the normalization W is computed only once for LMIN and applied
    // uniformly to all L, but Q = -1/(F^2+G^2) varies with L).
    int n_match_pt = N - 3;
    int n1_pt = n_match_pt - 1;
    int n2_pt = n_match_pt;
    double rho1_pt = k * (n1_pt * h);
    double rho2_pt = k * (n2_pt * h);

    std::vector<double> FC1_all(Lmax + 2), FCP1_all(Lmax + 2), GC1_all(Lmax + 2), GCP1_all(Lmax + 2);
    std::vector<double> FC2_all(Lmax + 2), FCP2_all(Lmax + 2), GC2_all(Lmax + 2), GCP2_all(Lmax + 2);
    Rcwfn(rho1_pt, eta, 0, Lmax, FC1_all, FCP1_all, GC1_all, GCP1_all);
    Rcwfn(rho2_pt, eta, 0, Lmax, FC2_all, FCP2_all, GC2_all, GCP2_all);

    for (int L = 0; L <= Lmax; ++L) {
        // Build f(r) array
        double LL1_d = (double)L * (L + 1);
        std::vector<std::complex<double>> f(N + 2, 0.0);
        for (int i = 0; i <= N; ++i) {
            double r = (i == 0) ? 1e-10
                     : ((i < ch.NSteps) ? ch.RGrid[i] : i * h);
            if (r < 1e-10) r = 1e-10;
            int idx = std::min(i, N - 1);
            double Vr  = ch.V_real[idx];
            double Wi  = ch.V_imag[idx];
            double Vc  = ch.V_coulomb[idx];
            double LL1 = (double)L * (L + 1);

            double f_re = k2 - LL1 / (r * r) - f_conv * Vc + f_conv * Vr;
            double f_im = f_conv * Wi;
            f[i] = std::complex<double>(f_re, f_im);
        }

        // Numerov integration
        std::vector<std::complex<double>> u(N + 2, 0.0);
        u[0] = 0.0;
        u[1] = std::pow(h, L + 1);
        if (std::abs(u[1]) < 1e-300) u[1] = 1e-300;

        int ISTRT = 0;

        for (int i = 1; i < N; ++i) {
            auto t1 = 2.0 * (1.0 - 5.0 * h2_12 * f[i]) * u[i];
            auto t2 = (1.0 + h2_12 * f[i-1]) * u[i-1];
            auto dn = 1.0 + h2_12 * f[i+1];
            if (std::abs(dn) < 1e-300) dn = 1e-300;
            u[i+1] = (t1 - t2) / dn;

            double mag = std::abs(u[i+1]);
            if (mag > BIGNUM) {
                double inv_mag = 1.0 / mag;
                double threshold = STEPI * mag;
                int new_istrt = i + 1;
                for (int j = ISTRT; j <= i + 1; ++j) {
                    if (std::abs(u[j]) >= threshold) { new_istrt = j; break; }
                }
                for (int j = ISTRT; j < new_istrt; ++j) u[j] = 0.0;
                ISTRT = new_istrt;
                for (int j = ISTRT; j <= i + 1; ++j) u[j] *= inv_mag;
            }
        }

        // S-matrix extraction: two-point matching using pre-computed Coulomb functions.
        // F,G from Rcwfn called with MINL=0,MAXL=Lmax (correct upward G recursion).
        // u(r) = A*G_L(r) + B*F_L(r)  =>  S = (B + iA)/(B - iA)
        int n1 = n1_pt, n2 = n2_pt;

        double f1 = FC1_all[L], g1 = GC1_all[L];
        double f2 = FC2_all[L], g2 = GC2_all[L];

        std::complex<double> u1 = u[n1];
        std::complex<double> u2 = u[n2];

        // det = f2*g1 - f1*g2
        double det = f2 * g1 - f1 * g2;

        // A = G-coefficient, B = F-coefficient
        std::complex<double> A = (f2 * u1 - u2 * f1) / det;
        std::complex<double> B = (u2 * g1 - g2 * u1) / det;

        // S = (B + iA) / (B - iA)
        using namespace std::complex_literals;
        std::complex<double> Snum = B + 1i * A;
        std::complex<double> Sden = B - 1i * A;
        SMatrix[L] = Snum / Sden;
    }

    std::cout << "S-matrix computed.\n\n";

    // --- Print S-matrix table ---
    std::cout << std::setw(5) << "L"
              << std::setw(14) << "Re(S)"
              << std::setw(14) << "Im(S)"
              << std::setw(12) << "|S|" << "\n";
    std::cout << std::string(45, '-') << "\n";

    // Save S-matrix to file
    std::ofstream smat_file("sm148_aa_smat_cpp.txt");
    smat_file << "# 148Sm(a,a) 50 MeV S-matrix from C++ elastic code\n";
    smat_file << "# L  Re(S)  Im(S)\n";

    for (int L = 0; L <= Lmax; ++L) {
        double re = SMatrix[L].real();
        double im = SMatrix[L].imag();
        double mod = std::abs(SMatrix[L]);
        std::cout << std::setw(5) << L
                  << std::scientific << std::setprecision(6)
                  << std::setw(16) << re
                  << std::setw(16) << im
                  << std::fixed << std::setprecision(6)
                  << std::setw(12) << mod << "\n";
        smat_file << L << " " << std::scientific << std::setprecision(8)
                  << re << " " << im << "\n";
    }
    smat_file.close();
    std::cout << "\nS-matrix saved to: sm148_aa_smat_cpp.txt\n\n";

    // --- Compute DCS ---
    std::cout << "Computing DCS for θ = 1°..180°...\n";
    std::ofstream xsec_file("sm148_aa_xsec_cpp.txt");
    xsec_file << "# 148Sm(a,a) 50 MeV DCS from C++ elastic code\n";
    xsec_file << "# theta_deg  dcs_mb_sr\n";

    std::cout << std::setw(10) << "theta(deg)"
              << std::setw(18) << "DCS (mb/sr)" << "\n";
    std::cout << std::string(28, '-') << "\n";

    for (int ith = 1; ith <= 180; ++ith) {
        double theta = (double)ith;
        double dcs_fm2 = ElasticDCS(theta, k, eta, SMatrix, CoulPhase);
        double dcs_mb  = dcs_fm2 * 10.0;   // fm² → mb/sr
        std::cout << std::setw(10) << theta
                  << std::setw(18) << std::scientific << std::setprecision(4) << dcs_mb << "\n";
        xsec_file << theta << " " << std::scientific << std::setprecision(8) << dcs_mb << "\n";
    }
    xsec_file.close();
    std::cout << "\nDCS saved to: sm148_aa_xsec_cpp.txt\n";
    std::cout << "\n=== Done ===\n";
    return 0;
}
