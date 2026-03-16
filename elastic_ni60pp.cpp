// elastic_ni60pp.cpp
// Elastic scattering: 60Ni(p,p) at 30 MeV
// Spin-1/2 projectile: volume WS + surface WS + spin-orbit + Coulomb
// Two S-matrix values per L: S_{J=L+1/2} and S_{J=L-1/2}
// DCS (unpolarized): spin-averaged Coulomb + nuclear partial-wave sum

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include "include/rcwfn.h"

static constexpr double PI    = 3.14159265358979323846;
static constexpr double AMU   = 931.494061;   // MeV
static constexpr double HBARC = 197.32697;    // MeV·fm

// === Coulomb phase sigma_L = Im[lnGamma(L+1+i*eta)] ===
static double CoulombPhase(int L, double eta) {
    int N0 = 200;
    double sum_atan = 0;
    for (int n = 1; n <= N0; ++n) sum_atan += std::atan(eta / n);
    double x = N0 + 1, y = eta;
    double theta_N0 = (x - 0.5)*std::atan2(y, x) - y/2*std::log(x*x + y*y)
                    + y*(1.0 - std::log(2*PI)/2);
    double sig0 = theta_N0 - sum_atan;
    for (int n = 1; n <= L; ++n) sig0 += std::atan(eta / n);
    return sig0;
}

// === Legendre polynomial P_L(x) ===
static double LegendreP(int L, double x) {
    if (L == 0) return 1.0;
    if (L == 1) return x;
    double pm2 = 1.0, pm1 = x;
    for (int l = 2; l <= L; ++l) {
        double p = ((2*l - 1)*x*pm1 - (l - 1)*pm2) / l;
        pm2 = pm1; pm1 = p;
    }
    return pm1;
}

// === Associated Legendre P_L^1(x) / sin(theta) = dP_L/d(cos theta) ===
// We need (1/sin θ) * dP_L/dθ for the spin-orbit DCS terms
// P_L^1(x) = -sin(theta) * dP_L/d(cos theta)
// So: dP_L/dθ = -P_L^1(x)   where x = cos θ
// For the DCS formula with spin-1/2 (McIntyre form):
// dσ/dΩ = |f_0|² + |g|²
// where f_0 = Coulomb + (1/2ik)*Σ[(L+1)(e^{2iσ}S_{L+}-1) + L(e^{2iσ}S_{L-}-1)] P_L
//         g  = (1/2ik)*Σ[e^{2iσ}(S_{L-} - S_{L+})] * dP_L/dθ  (spin-flip)
// But for unpolarized DCS: dσ/dΩ = |f_0|² + |g|² where we sum spin states
// Equivalent to: dσ/dΩ = (1/2) Σ_m |<m'|f|m>|² summed over initial/final m states

// === Derivative of Legendre polynomial: P_L'(x) = d[P_L(x)]/dx ===
static double LegendrePprime(int L, double x) {
    if (L == 0) return 0.0;
    if (L == 1) return 1.0;
    // Recurrence: (1-x^2)*P_L'(x) = L*(P_{L-1}(x) - x*P_L(x))
    double PL   = LegendreP(L, x);
    double PLm1 = LegendreP(L-1, x);
    if (std::abs(1.0 - x*x) < 1e-12) {
        // Near poles: use limiting formula L*(L+1)/2 for x→±1
        return (x > 0 ? 1.0 : (L%2==0 ? 1.0 : -1.0)) * L*(L+1)/2.0;
    }
    return L*(PLm1 - x*PL)/(1.0 - x*x);
}

// === Elastic DCS for spin-1/2 projectile (unpolarized) ===
// Ref: Satchler "Direct Nuclear Reactions" Ch. 3; Thompson & Nunes "Nuclear Reactions for Astrophysics"
// f(θ) = f_C(θ) δ_{m'm} + f_N(θ)
// Spin-non-flip: f = f_C + A(θ)
// Spin-flip:     g = B(θ)
// dσ/dΩ = |f|² + |g|²
// A(θ) = (1/2ik) Σ_L [(L+1)e^{2iσ_L}(S_{L+}-1) + L·e^{2iσ_L}(S_{L-}-1)] P_L(cosθ)
//   (+ L term adds J=L-1/2, (L+1) adds J=L+1/2 weight for spin average)
// B(θ) = (1/2ik) Σ_L e^{2iσ_L}(S_{L+} - S_{L-}) dP_L/dθ  ← spin-flip
// For L=0: only J=1/2, so S_{0-}=1 by convention; (L+1)=1, L=0.
static double ElasticDCS_pp(double theta_deg, double k, double eta,
    const std::vector<std::complex<double>>& Sp,   // S_{J=L+1/2}
    const std::vector<std::complex<double>>& Sm,   // S_{J=L-1/2}
    const std::vector<double>& CoulPhase) {

    if (theta_deg < 0.01) return std::numeric_limits<double>::quiet_NaN();
    double theta  = theta_deg * PI / 180.0;
    double cos_th = std::cos(theta);
    double sin_th = std::sin(theta);
    double sin_th2 = std::sin(theta / 2.0);

    int Lmax = (int)Sp.size() - 1;

    // Coulomb amplitude (spin non-flip only)
    double sigma0 = CoulPhase[0];
    double phi_C  = 2.0*sigma0 - 2.0*eta*std::log(sin_th2);
    double fc_mod = eta / (2.0*k*sin_th2*sin_th2);
    std::complex<double> f_C(-fc_mod*std::cos(phi_C), -fc_mod*std::sin(phi_C));

    // Nuclear amplitudes A (non-flip) and B (flip)
    std::complex<double> A(0, 0), B(0, 0);
    std::complex<double> two_ik(0, 2.0*k);

    for (int L = 0; L <= Lmax; ++L) {
        double sL  = CoulPhase[L];
        std::complex<double> e2sig(std::cos(2*sL), std::sin(2*sL));
        double PL  = LegendreP(L, cos_th);
        // dP_L/dθ = -sin(θ) * dP_L/d(cosθ)
        double dPL = -sin_th * LegendrePprime(L, cos_th);

        // S_{L-} = 1 for L=0 (no J=L-1/2 state)
        std::complex<double> sp = Sp[L];
        std::complex<double> sm = (L == 0) ? std::complex<double>(1, 0) : Sm[L];

        // Spin-non-flip (McIntyre 1960):
        // A += e^{2iσ} * [(L+1)*(S+ - 1) + L*(S- - 1)] * P_L / (2ik)
        A += e2sig * ((double)(L+1)*(sp - 1.0) + (double)L*(sm - 1.0)) * PL;
        // Spin-flip:
        // B += e^{2iσ} * (S- - S+) * dP_L/dθ / (2ik)
        B += e2sig * (sm - sp) * dPL;
    }
    A /= two_ik;
    B /= two_ik;

    // Total unpolarized DCS: |f_C + A|² + |B|²
    std::complex<double> f_tot = f_C + A;
    return std::norm(f_tot) + std::norm(B);  // fm²
}

int main() {
    std::cout << "=== Elastic scattering 60Ni(p,p) at 30 MeV (C++ code) ===\n\n";

    // Relativistic kinematics matching Raphael (nuclear masses from AME)
    double Ap=1.0, At=60.0, Zp=1.0, Zt=28.0, Elab=30.0;
    // Proton mass; 60Ni: Z=28, N=32, BE/A=8.7808 MeV
    double mp_MeV = 938.272013, mn_MeV = 939.565346;
    double mass_a = mp_MeV;                               // proton (MeV/c^2)
    double mass_A = 28*mp_MeV + 32*mn_MeV - 8.7808*60;  // 60Ni  (MeV/c^2)
    double E_tot  = std::sqrt((mass_a+mass_A)*(mass_a+mass_A) + 2*mass_A*Elab);
    double Ecm    = E_tot - (mass_a + mass_A);
    double mu     = mass_a*mass_A/(mass_a+mass_A);        // reduced mass (MeV/c^2)
    double k      = std::sqrt(2.0*mu*Ecm)/HBARC;          // fm^-1
    double eta    = Zp*Zt*(HBARC/137.036)*k/(2.0*Ecm);
    double f_conv = 2.0*mu/(HBARC*HBARC);                 // MeV^-1 fm^-2
    double e2   = HBARC/137.036;

    std::cout << "Kinematics:\n";
    std::cout << "  Elab = " << Elab << " MeV\n";
    std::cout << "  Ecm  = " << Ecm  << " MeV\n";
    std::cout << "  k    = " << k    << " fm^-1\n";
    std::cout << "  eta  = " << eta  << "\n";
    std::cout << "  mu   = " << mu   << " AMU\n\n";

    // Grid — use smaller step to match Raphael's dr=0.05 fm accuracy
    double h = 0.05;
    int N = 601;   // r_max = 30 fm

    // Potential: V + iW (volume) + iWs (surface) + Vso (spin-orbit) + Coulomb
    // Volume: V=-47.937, W=-2.853; r0=1.20, a=0.669
    // Surface: Wi=-6.878; rs=1.28, as=0.550  (derivative WS)
    // Spin-orbit: Vso=-5.250, Wso=0.162; rso=1.02, aso=0.590
    // Coulomb: rc=1.258
    double V0=47.937, W0=2.853, r0=1.20, av=0.669;
    double Ws=6.878,  rs=1.28,  as=0.550;
    double Vso=5.250, Wso=-0.162, rso=1.02, aso=0.590;
    double RC0=1.258;

    // Radii (non-R0MASS for A_proj=1): R = r0 * At^{1/3}
    double Rv  = r0  * std::cbrt(At);
    double Rws = rs  * std::cbrt(At);
    double Rso = rso * std::cbrt(At);
    double Rc  = RC0 * std::cbrt(At);

    std::cout << "Potential radii:\n";
    std::cout << "  Rv=" << Rv << " fm, Rws=" << Rws << " fm, Rso=" << Rso << " fm, Rc=" << Rc << " fm\n\n";

    // Pre-compute potential on grid
    // Vr[i] = volume real, Wi_vol[i] = volume imag, Wsurf[i] = surface imag
    // Vc[i] = Coulomb; VsoRe/VsoIm = spin-orbit (computed per-step in Numerov)
    std::vector<double> Vr_grid(N), Wi_vol_grid(N), Wsurf_grid(N), Vc_grid(N);
    std::vector<double> Vso_re_grid(N), Vso_im_grid(N);   // SO derivative factor

    for (int i = 0; i < N; ++i) {
        double r = (i == 0) ? 1e-10 : i*h;
        // Volume WS: f(r) = 1/(1+exp((r-R)/a))
        double fvol = 1.0/(1+std::exp((r-Rv)/av));
        Vr_grid[i]    = V0*fvol;
        Wi_vol_grid[i] = W0*fvol;
        // Surface imaginary: matches Raphael WS_SurfacePot.output(x) = 4*Ws*exp/(1+exp)^2
        // Note: Raphael does NOT divide by a_s or include 1/r — it's 4*W_s * [-a*df/dr]
        // which equals 4*W_s * exp/(1+exp)^2 directly (no extra a or r factors)
        double exp_s = std::exp((r-Rws)/as);
        Wsurf_grid[i] = 4.0*Ws*exp_s/((1+exp_s)*(1+exp_s));  // MeV (negative = absorptive)
        // Coulomb
        Vc_grid[i] = (r < Rc) ?
            Zp*Zt*e2*(3.0 - r*r/(Rc*Rc))/(2.0*Rc) :
            Zp*Zt*e2/r;
        // Spin-orbit potential matching Raphael convention:
        // V_SO(r) = LS_val * 4*Vso * (-df_so/dr) / r
        // where f_so(r) = 1/(1+exp((r-Rso)/aso))
        // Raphael's SpinOrbit_Pot.output(x) = 4*Vso * exp/((1+exp)^2 * aso * x)
        // Note: no Thomas factor — Vso is used directly in MeV
        double exp_so = std::exp((r-Rso)/aso);
        double fso_deriv_over_r = (r > 1e-5) ?
            -exp_so/(aso*(1+exp_so)*(1+exp_so)*r) : 0.0;
        // Factor 4 matches Raphael's SpinOrbit_Pot.output
        Vso_re_grid[i] = 4.0*Vso*fso_deriv_over_r;   // per-unit-LS, sign: -df/dr
        Vso_im_grid[i] = 4.0*Wso*fso_deriv_over_r;
    }

    // Precompute Coulomb functions at matching points
    int Lmax = 25;
    int n2_pt = N - 3, n1_pt = n2_pt - 1;
    double rho1_pt = k*n1_pt*h, rho2_pt = k*n2_pt*h;
    std::vector<double> FC1(Lmax+2),FCP1(Lmax+2),GC1(Lmax+2),GCP1(Lmax+2);
    std::vector<double> FC2(Lmax+2),FCP2(Lmax+2),GC2(Lmax+2),GCP2(Lmax+2);
    Rcwfn(rho1_pt, eta, 0, Lmax, FC1, FCP1, GC1, GCP1);
    Rcwfn(rho2_pt, eta, 0, Lmax, FC2, FCP2, GC2, GCP2);

    // Precompute Coulomb phases
    std::vector<double> CoulPhase(Lmax+1);
    for (int L = 0; L <= Lmax; ++L) CoulPhase[L] = CoulombPhase(L, eta);

    // S-matrix for J=L+1/2 and J=L-1/2
    std::vector<std::complex<double>> Sp(Lmax+1), Sm(Lmax+1);

    std::cout << "Computing S-matrix for L=0.." << Lmax << "...\n";

    double k2 = k*k;
    double h2_12 = h*h/12.0;
    const double BIGNUM = 1e30;

    for (int L = 0; L <= Lmax; ++L) {
        double LL1 = (double)L*(L+1);

        // Compute S for both J = L+1/2 and J = L-1/2
        // L·S eigenvalue: (J(J+1) - L(L+1) - S(S+1))/2, S=1/2, S(S+1)=3/4
        // For J=L+1/2: LS = L/2
        // For J=L-1/2: LS = -(L+1)/2
        double LS_plus  = (L > 0) ?  (double)L/2.0    : 0.5;   // L=0: J=1/2=S only
        double LS_minus = -(double)(L+1)/2.0;

        for (int spin = 0; spin <= 1; ++spin) {
            // spin=0: J=L+1/2, spin=1: J=L-1/2
            if (L == 0 && spin == 1) { Sm[0] = {1, 0}; continue; }

            double LS_val = (spin == 0) ? LS_plus : LS_minus;

            // Build f(r) array for this L, J
            std::vector<std::complex<double>> f(N+2, 0.0);
            for (int i = 0; i <= N; ++i) {
                double r = (i == 0) ? 1e-10 : i*h;
                int idx = std::min(i, N-1);
                // Numerov f = k^2 - L(L+1)/r^2 - f_conv*V_total
                // Sign conventions (verified vs Raphael at r=5, L=3, J=3.5):
                // fr = k^2 - LL1/r^2 - f_conv*Vc + f_conv*Vr - f_conv*LS*Vso_re
                //   (Vso_re_grid is negative for attractive SO → minus flips to +attractive)
                // fi = f_conv*(Wi_vol - Wsurf - LS*Vso_im)
                //   (Wsurf_grid<0 → -Wsurf>0=absorptive; Vso_im_grid>0 → -LS*Vso_im<0)
                double fr = k2 - LL1/(r*r) - f_conv*Vc_grid[idx]
                             + f_conv*Vr_grid[idx]
                             - f_conv*LS_val*Vso_re_grid[idx];
                double fi = f_conv*Wi_vol_grid[idx]
                           + f_conv*Wsurf_grid[idx]
                           - f_conv*LS_val*Vso_im_grid[idx];
                f[i] = {fr, fi};
            }

            // Numerov integration
            std::vector<std::complex<double>> u(N+2, 0.0);
            u[0] = 0.0; u[1] = std::pow(h, L+1);
            if (std::abs(u[1]) < 1e-300) u[1] = 1e-300;
            for (int i = 1; i < N; ++i) {
                auto t1 = 2.0*(1.0 - 5.0*h2_12*f[i])*u[i];
                auto t2 = (1.0 + h2_12*f[i-1])*u[i-1];
                auto dn = 1.0 + h2_12*f[i+1];
                if (std::abs(dn) < 1e-300) dn = 1e-300;
                u[i+1] = (t1 - t2)/dn;
                if (std::abs(u[i+1]) > BIGNUM)
                    for (int j = 0; j <= i+1; ++j) u[j] /= std::abs(u[i+1]);
            }

            // Two-point S-matrix extraction
            double f1=FC1[L], g1=GC1[L], f2=FC2[L], g2=GC2[L];
            std::complex<double> u1=u[n1_pt], u2=u[n2_pt];
            double det = f2*g1 - f1*g2;
            std::complex<double> A = (f2*u1 - u2*f1)/det;
            std::complex<double> B = (u2*g1 - g2*u1)/det;
            using namespace std::complex_literals;
            auto S = (B + 1i*A)/(B - 1i*A);

            if (spin == 0) Sp[L] = S;
            else           Sm[L] = S;
        }
    }

    std::cout << "S-matrix computed.\n\n";

    // Print and save S-matrix
    std::cout << std::setw(4) << "L"
              << std::setw(14) << "Re(S+)"  << std::setw(14) << "Im(S+)"
              << std::setw(14) << "Re(S-)"  << std::setw(14) << "Im(S-)" << "\n";
    std::cout << std::string(56, '-') << "\n";

    std::ofstream smat_file("ni60pp_smat_cpp.txt");
    smat_file << "# 60Ni(p,p) 30 MeV S-matrix from C++\n";
    smat_file << "# L  Re(S+)  Im(S+)  Re(S-)  Im(S-)\n";
    for (int L = 0; L <= Lmax; ++L) {
        std::cout << std::setw(4) << L
                  << std::scientific << std::setprecision(5)
                  << std::setw(14) << Sp[L].real() << std::setw(14) << Sp[L].imag()
                  << std::setw(14) << Sm[L].real() << std::setw(14) << Sm[L].imag() << "\n";
        smat_file << L << " " << std::scientific << std::setprecision(8)
                  << Sp[L].real() << " " << Sp[L].imag() << " "
                  << Sm[L].real() << " " << Sm[L].imag() << "\n";
    }
    smat_file.close();
    std::cout << "\nS-matrix saved to: ni60pp_smat_cpp.txt\n\n";

    // Compute DCS
    std::cout << "Computing DCS...\n";
    std::ofstream xsec_file("ni60pp_xsec_cpp.txt");
    xsec_file << "# 60Ni(p,p) 30 MeV DCS from C++ (mb/sr)\n# theta_deg  dcs_mb_sr\n";
    for (int ith = 1; ith <= 180; ++ith) {
        double dcs_fm2 = ElasticDCS_pp(ith, k, eta, Sp, Sm, CoulPhase);
        double dcs_mb  = dcs_fm2 * 10.0;
        xsec_file << ith << " " << std::scientific << std::setprecision(6) << dcs_mb << "\n";
    }
    xsec_file.close();
    std::cout << "DCS saved to: ni60pp_xsec_cpp.txt\n";

    return 0;
}
