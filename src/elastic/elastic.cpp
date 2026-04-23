// elastic.cpp — Unified elastic scattering solver implementation
// Sign conventions verified against Raphael (Python) for:
//   148Sm(α,α) 50 MeV,  60Ni(p,p) 30 MeV,  60Ni(d,d) 60 MeV
#include "elastic.h"
#include "rcwfn.h"
// Forward-declare only — do NOT include math_utils.h here (JSymbols.h ODR issue)
double ClebschGordan(double J1, double m1, double J2, double m2, double J, double m);
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>

static constexpr double PI     = 3.14159265358979323846;
static constexpr double HBARC  = 197.326979;   // MeV·fm
static constexpr double MP_MEV = 938.272013;    // proton mass MeV
static constexpr double MN_MEV = 939.565346;    // neutron mass MeV
static constexpr double E2     = HBARC / 137.036;  // e² in MeV·fm = 1.43996

// ============================================================
// Construction / setup
// ============================================================
ElasticSolver::ElasticSolver()
    : Ap_(0), Zp_(0), At_(0), Zt_(0), Elab_(0),
      S_(0), h_(0.05), N_(601), Lmax_(-1),
      k_(0), eta_(0), mu_(0), Ecm_(0), f_conv_(0),
      rc0_(0), hasCoulomb_(false) {}

void ElasticSolver::SetSystem(int Ap, int Zp, int At, int Zt, double Elab) {
    Ap_ = Ap; Zp_ = Zp; At_ = At; Zt_ = Zt; Elab_ = Elab;
    // Determine projectile spin from mass number
    // A=1 (p,n) → S=1/2;  A=2 (d) → S=1;  A=3 (t,He3) → S=1/2;  A=4 (α) → S=0
    if      (Ap == 1) S_ = 0.5;
    else if (Ap == 2) S_ = 1.0;
    else if (Ap == 3) S_ = 0.5;
    else              S_ = 0.0;   // alpha and heavier ions treated as spin-0
}

void ElasticSolver::SetGrid(double h, int N) { h_ = h; N_ = N; }

void ElasticSolver::SetLmax(int Lmax) { Lmax_ = Lmax; }

void ElasticSolver::AddVolumeWS(std::complex<double> V, double r0, double a0) {
    pots_.push_back({kVolumeWS, V, r0, a0});
}
void ElasticSolver::AddSurfaceWS(std::complex<double> V, double r0, double a0) {
    pots_.push_back({kSurfaceWS, V, r0, a0});
}
void ElasticSolver::AddSpinOrbit(std::complex<double> Vso, double r0, double a0) {
    pots_.push_back({kSpinOrbit, Vso, r0, a0});
}
void ElasticSolver::AddCoulomb(double rc0) {
    rc0_ = rc0; hasCoulomb_ = true;
}
void ElasticSolver::ClearPotentials() {
    pots_.clear(); hasCoulomb_ = false; rc0_ = 0;
}

// ============================================================
// Nuclear mass (MeV/c²) — relativistic, matching Raphael
// For common projectiles use exact masses; targets use binding energy.
// ============================================================
double ElasticSolver::NuclearMass(int A, int Z) const {
    // Common light projectiles (exact AME values)
    if (A == 1 && Z == 1) return 938.272013;  // proton
    if (A == 1 && Z == 0) return 939.565346;  // neutron
    if (A == 2 && Z == 1) return 1875.612928; // deuteron (AME2020)
    if (A == 3 && Z == 1) return 2808.920930; // triton
    if (A == 3 && Z == 2) return 2808.391500; // He-3
    if (A == 4 && Z == 2) return 3727.379378; // He-4 (alpha)
    // For heavier nuclei: use Z*mp + N*mn - BE/A * A
    // We use a simple estimate: (A-Z)*MN + Z*MP - 8.5*A (rough BE ~8.5 MeV/A)
    // This is approximate but gives correct kinematics to ~0.1%
    // Better: look up in AME table via Isotope class, but keep it self-contained here
    double N = A - Z;
    // Use liquid-drop BE/A ≈ 15.84 - 18.33/A^{1/3} - 0.71*Z(Z-1)/A^{4/3} - 23.2*(A-2Z)²/A² (rough)
    // For simplicity, use the Isotope class data file if available, else estimate
    // Actually let's just use a table for common targets:
    struct { int A, Z; double mass_MeV; } table[] = {
        {12, 6,  11174.864},  // 12C
        {16, 8,  14895.080},  // 16O
        {28,14,  26060.054},  // 28Si
        {40,20,  37224.058},  // 40Ca  (ME=-34.846 MeV)
        {48,20,  44657.277},  // 48Ca  (ME=-44.215 MeV)
        {58,28,  53955.344},  // 58Ni
        {60,28,  55814.996},  // 60Ni (Z=28, N=32, BE/A=8.7808)
        {90,40,  83729.000},  // 90Zr
        {120,50,111876.000},  // 120Sn
        {148,62,137801.069},  // 148Sm (Z=62, N=86, BE/A=8.1521)
        {208,82,193729.000},  // 208Pb
    };
    for (auto& e : table)
        if (e.A == A && e.Z == Z) return e.mass_MeV;
    // Fallback: rough estimate
    return Z * MP_MEV + N * MN_MEV - 8.5 * A;
}

// ============================================================
// Kinematics (relativistic, matching Raphael)
// ============================================================
void ElasticSolver::CalcKinematics() {
    static constexpr double AMUMEV_PTL = 931.5016;  // Ptolemy AMUMEV
    static constexpr double HBARC_PTL  = 197.32858; // Ptolemy HBARC
    static constexpr double AFINE_PTL  = 137.03604; // Ptolemy AFINE
    static constexpr double E2_PTL     = HBARC_PTL / AFINE_PTL;

    if (ptolemyMass_) {
        // Ptolemy convention: A*AMUMEV + mass-excess correction (matching Fortran source line ~32371)
        // Mass excess (MeV): ME = nuclear_mass - A*u  (from AME; Fortran reads from input AMXGPT)
        // Electron mass correction: -Z * m_e  (Fortran: -IZ*(EMASS/AMUMEV)*AMUMEV = -IZ*EMASS)
        // For common cases, use hardcoded AME values (same as NuclearMass table but in AMU)
        static constexpr double EMASS_PTL = 0.511;  // MeV
        auto ptolMass = [&](int A, int Z) -> double {
            // Return mass in AMU = (AME_mass_MeV - Z*EMASS) / AMUMEV_PTL
            // Use exact nuclear masses from NuclearMass() which already exclude electrons
            return NuclearMass(A, Z) / AMUMEV_PTL;
        };
        double ma_amu = ptolMass(Ap_, Zp_);
        double mA_amu = ptolMass(At_, Zt_);
        mu_    = AMUMEV_PTL * ma_amu * mA_amu / (ma_amu + mA_amu);
        Ecm_   = Elab_ * mA_amu / (ma_amu + mA_amu);   // non-relativistic Ecm
        k_     = std::sqrt(2.0*mu_*Ecm_) / HBARC_PTL;
        eta_   = Zp_ * Zt_ * E2_PTL * k_ / (2.0 * Ecm_);
        f_conv_ = 2.0 * mu_ / (HBARC_PTL * HBARC_PTL);
    } else {
    double ma = NuclearMass(Ap_, Zp_);
    double mA = NuclearMass(At_, Zt_);
    double E_tot = std::sqrt((ma + mA)*(ma + mA) + 2.0*mA*Elab_);
    Ecm_   = E_tot - ma - mA;
    mu_    = ma * mA / (ma + mA);        // MeV/c^2
    k_     = std::sqrt(2.0*mu_*Ecm_) / HBARC;
    eta_   = Zp_ * Zt_ * E2 * k_ / (2.0 * Ecm_);
    f_conv_ = 2.0 * mu_ / (HBARC * HBARC);
    }

    // Auto Lmax: largest L where S-matrix differs appreciably from 1
    // Classical turning: rho_turn ≈ eta + sqrt(eta² + L(L+1))
    // Use Lmax where rho_turn < k * (r_nuclear + 2*a_max) to be safe
    if (Lmax_ < 0) {
        double r_nuclear = 1.5 * std::cbrt(At_);
        Lmax_ = (int)(k_ * r_nuclear + 5) + 10;   // generous estimate
    }
}

void ElasticSolver::PrintKinematics() const {
    std::cout << "Kinematics:\n"
              << "  Elab = " << Elab_ << " MeV\n"
              << "  Ecm  = " << Ecm_  << " MeV\n"
              << "  k    = " << k_    << " fm^-1\n"
              << "  eta  = " << eta_  << "\n"
              << "  mu   = " << mu_   << " MeV/c^2\n"
              << "  Lmax = " << Lmax_ << "\n"
              << "  Spin = " << S_    << "\n\n";
}

// ============================================================
// Coulomb phase σ_L = Im[lnΓ(L+1+iη)]
// ============================================================
double ElasticSolver::CoulombPhase(int L) const {
    // Stirling at N0=200 then recurse down
    int N0 = 200;
    double sum_atan = 0;
    for (int n = 1; n <= N0; ++n) sum_atan += std::atan(eta_ / n);
    double x = N0 + 1, y = eta_;
    // Stirling approximation for Im[ln Gamma(x + iy)]:
    // Im[(z-0.5)*ln(z) - z + 0.5*ln(2π)] = (x-0.5)*atan2(y,x) + y/2*ln(x²+y²) - y
    // (the 0.5*ln(2π) is real, contributes 0 to imaginary part)
    double theta_N0 = (x-0.5)*std::atan2(y, x) + y/2.0*std::log(x*x + y*y) - y;
    double sig = theta_N0 - sum_atan;
    for (int n = 1; n <= L; ++n) sig += std::atan(eta_ / n);
    return sig;
}

// ============================================================
// Legendre polynomials
// ============================================================
double ElasticSolver::LegendreP(int L, double x) {
    if (L == 0) return 1.0;
    if (L == 1) return x;
    double pm2 = 1.0, pm1 = x;
    for (int l = 2; l <= L; ++l) {
        double p = ((2*l-1)*x*pm1 - (l-1)*pm2) / l;
        pm2 = pm1; pm1 = p;
    }
    return pm1;
}
double ElasticSolver::LegendrePprime(int L, double x) {
    if (L == 0) return 0.0;
    if (L == 1) return 1.0;
    double PL = LegendreP(L, x), PLm1 = LegendreP(L-1, x);
    double sin2 = 1.0 - x*x;
    if (sin2 < 1e-12) return (x > 0 ? 1.0 : (L%2==0?1.0:-1.0)) * 0.5*L*(L+1);
    return L*(PLm1 - x*PL) / sin2;
}

// ============================================================
// Build potential arrays on radial grid
// ============================================================
static void BuildPotentialArrays(
    const std::vector<PotEntry>& pots, bool hasCoulomb, double rc0,
    int Ap, int Zp, int Zt, int At,
    double h, int N,
    std::vector<double>& Vr, std::vector<double>& Wi,
    std::vector<double>& Vc,
    std::vector<double>& VsoRe, std::vector<double>& VsoIm)
{
    Vr.assign(N, 0); Wi.assign(N, 0); Vc.assign(N, 0);
    VsoRe.assign(N, 0); VsoIm.assign(N, 0);

    // R0MASS convention — from Ptolemy source.mor lines 28572-28581 (R0DEFAULT):
    //   AMMORE = heavier mass, AMLESS = lighter mass
    //   R0MASS = AMMORE^(1/3)
    //   IF AMLESS > 2.5: R0MASS += AMLESS^(1/3)   (add projectile only if A > 2.5)
    // This means:
    //   p, n, d (Ap <= 2): R0MASS = At^(1/3)  only
    //   t, He3, alpha, heavier (Ap > 2.5, i.e. Ap >= 3): R0MASS = At^(1/3) + Ap^(1/3)
    double At13 = std::cbrt((double)At);
    double Ap13 = (Ap > 2.5) ? std::cbrt((double)Ap) : 0.0;
    double R0MASS = At13 + Ap13;
    double Rc   = hasCoulomb ? rc0 * R0MASS : 0.0;  // Coulomb follows same R0MASS convention

    for (int i = 0; i < N; ++i) {
        double r = (i == 0) ? 1e-10 : i * h;

        // Coulomb
        // If hasCoulomb: use uniform-sphere for r<Rc, point for r>=Rc.
        // If Rc==0 (rc0=0, e.g. Koning convention): pure point Coulomb (Raphael-compatible).
        if (hasCoulomb) {
            if (Rc > 0.0 && r < Rc)
                Vc[i] = Zp * Zt * E2 * (3.0 - r*r/(Rc*Rc)) / (2.0*Rc);
            else
                Vc[i] = Zp * Zt * E2 / r;
        }

        // Nuclear potentials
        for (const auto& p : pots) {
            double R0 = p.r0 * R0MASS;
            double ex = std::exp((r - R0) / p.a0);

            switch (p.type) {
            case kVolumeWS: {
                // V(r) = V0 / (1+exp)
                // Raphael: V0 is signed (negative attractive)
                // f_re contribution: -f_conv * Re(V0) * fWS → +f_conv*|ReV| for attract.
                // We store Vr as -Re(V)*f (so Vr>0 for attractive Re(V)<0)
                // Wi as -Im(V)*f (so Wi>0 for absorptive Im(V)<0)
                double fws = 1.0 / (1.0 + ex);
                Vr[i] += -p.V.real() * fws;  // > 0 if V.real() < 0 (attractive)
                Wi[i] += -p.V.imag() * fws;  // > 0 if V.imag() < 0 (absorptive)
                break;
            }
            case kSurfaceWS: {
                // Raphael: output = 4*V0*exp/(1+exp)^2  (no 1/a, no 1/r)
                // Surface is added directly to potential (not derivative in strict sense)
                // Typically V.real()=0, V.imag()<0 (absorptive surface)
                double fsurf = 4.0 * ex / ((1.0 + ex)*(1.0 + ex));
                Vr[i] += -p.V.real() * fsurf;
                Wi[i] += -p.V.imag() * fsurf;
                break;
            }
            case kSpinOrbit: {
                // Raphael: output = 4*Vso * exp / (a*(1+exp)^2 * r)  [per unit LS]
                // VsoRe/VsoIm store the radial factor (signed, per unit LS)
                // f_re contribution: -f_conv*LS*Re(V_so) → need: -f_conv*LS*(-|Vso|*factor)
                // store as: VsoRe = -Re(Vso) * factor  so: -f_conv*LS*VsoRe is correct
                if (r > 1e-5) {
                    double fso = 4.0 * ex / (p.a0 * (1.0 + ex)*(1.0 + ex) * r);
                    VsoRe[i] += -p.V.real() * fso;
                    VsoIm[i] += -p.V.imag() * fso;
                }
                break;
            }
            case kCoulomb:
                break;  // handled separately
            }
        }
    }
}

// ============================================================
// Numerov integration + two-point S-matrix extraction
// ============================================================
std::complex<double> ElasticSolver::RunNumerov(
    int L, double LS_val,
    const std::vector<double>& Vr, const std::vector<double>& Wi,
    const std::vector<double>& Vc,
    const std::vector<double>& VsoRe, const std::vector<double>& VsoIm,
    const std::vector<double>& FC1, const std::vector<double>& GC1,
    const std::vector<double>& FC2, const std::vector<double>& GC2,
    std::vector<std::complex<double>>* wf_out) const
{
    double k2 = k_*k_, h2_12 = h_*h_/12.0;
    double LL1 = (double)L*(L+1);
    // Matching points: back = N_-4, end = N_ (Fortran: NSTEP+1-NBAKCM, NSTEP+1)
    int n_back = N_ - 4, n_end = N_;

    // Build f(r) array
    // f_re = k² - L(L+1)/r² - f_conv*Vc + f_conv*(Vr + LS*VsoRe)
    // f_im = f_conv*(Wi + LS*VsoIm)
    // where Vr, Wi, VsoRe, VsoIm are all stored as POSITIVE for attractive/absorptive
    std::vector<std::complex<double>> f(N_ + 2, 0.0);
    for (int i = 0; i <= N_; ++i) {
        double r = (i == 0) ? 1e-10 : i * h_;
        int idx = std::min(i, N_-1);
        double fr = k2 - LL1/(r*r) - f_conv_*Vc[idx]
                  + f_conv_*(Vr[idx] + LS_val * VsoRe[idx]);
        double fi = f_conv_*(Wi[idx] + LS_val * VsoIm[idx]);
        f[i] = {fr, fi};
    }

    // Numerov integration: u[0]=0, u[1]=r^{L+1}
    std::vector<std::complex<double>> u(N_ + 2, 0.0);
    u[0] = 0.0;
    u[1] = std::pow(h_, L + 1);
    if (std::abs(u[1]) < 1e-300) u[1] = 1e-300;

    const double BIGNUM = 1e30;
    for (int i = 1; i < N_; ++i) {
        auto t1 = 2.0*(1.0 - 5.0*h2_12*f[i])*u[i];
        auto t2 = (1.0 + h2_12*f[i-1])*u[i-1];
        auto dn = 1.0 + h2_12*f[i+1];
        if (std::abs(dn) < 1e-300) dn = 1e-300;
        u[i+1] = (t1 - t2) / dn;
        if (std::abs(u[i+1]) > BIGNUM)
            for (int j = 0; j <= i+1; ++j) u[j] /= std::abs(u[i+1]);
    }

    // Two-point Wronskian matching: u = A*G_L + B*F_L
    // Use Fortran-matching points: back=N_-4, end=N_
    // FC1/GC1 = F,G at back; FC2/GC2 = F1,G1 at end
    double fb = FC1[L], gb = GC1[L];  // Coulomb at back (SUMMAX-4h)
    double fe = FC2[L], ge = GC2[L];  // Coulomb at end  (SUMMAX)
    std::complex<double> ub = u[n_back], ue = u[n_end];
    // Solve: ub = A*gb + B*fb,  ue = A*ge + B*fe
    double det = fe*gb - fb*ge;
    if (std::abs(det) < 1e-30) return {1.0, 0.0};
    std::complex<double> A = (fe*ub - ue*fb) / det;
    std::complex<double> B = (ue*gb - ge*ub) / det;
    using namespace std::complex_literals;
    std::complex<double> norm = B - 1i*A;
    if (wf_out) {
        *wf_out = u;
        if (std::abs(norm) > 1e-30)
            for (auto& v : *wf_out) v /= norm;
    }
    return (B + 1i*A) / norm;
}

// ============================================================
// Main S-matrix computation
// ============================================================
void ElasticSolver::CalcScatteringMatrix() {
    if (k_ <= 0) CalcKinematics();

    // Auto-extend grid to cover classical turning point for all L (matching Fortran WAVPOT).
    // Fortran: RMAX = MAX(RMAX, (eta + sqrt(eta^2 + L*(L+1))) / k)
    {
        double Rmax = N_ * h_;
        for (int L = 0; L <= Lmax_; ++L) {
            double LL1 = (double)L * (L + 1);
            double rTurn = (eta_ + std::sqrt(eta_*eta_ + LL1)) / k_;
            if (rTurn > Rmax) Rmax = rTurn;
        }
        int Nnew = (int)(Rmax / h_ + 0.5);
        if (Nnew > N_) N_ = Nnew;
    }

    // Build potential arrays
    std::vector<double> Vr, Wi, Vc, VsoRe, VsoIm;
    BuildPotentialArrays(pots_, hasCoulomb_, rc0_, Ap_, Zp_, Zt_, At_,
                         h_, N_, Vr, Wi, Vc, VsoRe, VsoIm);

    // Matching points: n_back = N_-4 (NSTEP+1-NBAKCM in Fortran, NBAKCM=4)
    //                  n_end  = N_   (NSTEP+1 in Fortran, 1-indexed)
    // Fortran SETFG: F,G at R1=SUMMAX-4h (back), F1,G1 at R2=SUMMAX (end)
    int n_back = N_ - 4, n_end = N_;
    double rho_back = k_ * n_back * h_, rho_end = k_ * n_end * h_;
    std::vector<double> FC1(Lmax_+2), FCP1(Lmax_+2), GC1(Lmax_+2), GCP1(Lmax_+2);
    std::vector<double> FC2(Lmax_+2), FCP2(Lmax_+2), GC2(Lmax_+2), GCP2(Lmax_+2);
    Rcwfn(rho_back, eta_, 0, Lmax_, FC1, FCP1, GC1, GCP1);  // F,G at back point
    Rcwfn(rho_end,  eta_, 0, Lmax_, FC2, FCP2, GC2, GCP2);  // F1,G1 at end point

    // Coulomb phases
    CoulPhase_.resize(Lmax_ + 1);
    for (int L = 0; L <= Lmax_; ++L)
        CoulPhase_[L] = CoulombPhase(L);

    // Allocate S-matrix and wavefunctions: [L][idx]
    int nJ = (int)(2*S_) + 1;  // number of J channels per L
    Smat_.assign(Lmax_ + 1, std::vector<std::complex<double>>(nJ, {1,0}));
    wfunc_.assign(Lmax_ + 1, std::vector<std::vector<std::complex<double>>>(nJ));

    bool hasSO = false;
    for (const auto& p : pots_) if (p.type == kSpinOrbit) { hasSO = true; break; }

    for (int L = 0; L <= Lmax_; ++L) {
        if (S_ == 0.0) {
            // Spin-0: single channel, no SO splitting
            Smat_[L][0] = RunNumerov(L, 0.0, Vr, Wi, Vc, VsoRe, VsoIm,
                                     FC1, GC1, FC2, GC2, &wfunc_[L][0]);
        } else {
            // Spin-S (S>0): compute all J channels even if no SO potential
            for (int idx = 0; idx < nJ; ++idx) {
                double J = L - S_ + idx;
                double Jmin = std::abs(L - S_);
                if (J < Jmin || J < 0) {
                    Smat_[L][idx] = {0,0};
                    wfunc_[L][idx].assign(N_+2, {0,0});
                    continue;
                }
                double LS_val = (J*(J+1) - L*(L+1) - S_*(S_+1)) / 2.0;
                Smat_[L][idx] = RunNumerov(L, LS_val, Vr, Wi, Vc, VsoRe, VsoIm,
                                            FC1, GC1, FC2, GC2, &wfunc_[L][idx]);
            }
        }
    }
}

// ============================================================
// DCS computation — general spin formula (matches Raphael exactly)
// ============================================================
std::complex<double> ElasticSolver::CoulombAmp(double theta_deg) const {
    double theta = theta_deg * PI / 180.0;
    double sin2  = std::sin(theta / 2.0);
    double phi_C = 2.0*CoulPhase_[0] - 2.0*eta_*std::log(sin2);
    double fmod  = eta_ / (2.0 * k_ * sin2 * sin2);
    return std::complex<double>(-fmod*std::cos(phi_C), -fmod*std::sin(phi_C));
}

// Associated Legendre P_l^m(x), m >= 0, Condon-Shortley convention
// P_l^m — standard associated Legendre WITHOUT Condon-Shortley phase.
// Matches Raphael (assLegendreP.py) convention: no (-1)^m factor.
// P_m^m = (1-2m) * P_{m-1}^{m-1} * sqrt(1-x^2)  (Raphael recurrence)
double ElasticSolver::AssocLegendreP(int l, int m, double x) {
    if (m < 0 || m > l) return 0.0;
    // Start from P_0^0 = 1 and build up using Raphael's recurrence
    double pmm = 1.0;  // P_0^0
    double omx2_sqrt = std::sqrt((1.0 - x)*(1.0 + x));  // sqrt(1-x^2)
    for (int i = 1; i <= m; ++i) {
        pmm = (1 - 2*i) * pmm * omx2_sqrt;  // P_i^i = (1-2i)*P_{i-1}^{i-1}*sqrt(1-x^2)
    }
    if (l == m) return pmm;
    double pmmp1 = x * (2*m + 1) * pmm;
    if (l == m+1) return pmmp1;
    double pll = 0;
    for (int ll = m+2; ll <= l; ++ll) {
        pll = ((2*ll - 1)*x*pmmp1 - (ll + m - 1)*pmm) / (ll - m);
        pmm = pmmp1; pmmp1 = pll;
    }
    return pll;
}

// Thin CG wrapper (math_utils.h → JSymbols.h)
double ElasticSolver::CG(double j1, double m1, double j2, double m2, double J, double M) {
    return ClebschGordan(j1, m1, j2, m2, J, M);
}

// ============================================================
// Wynn's epsilon algorithm — accelerates convergence of a slowly
// converging series.  Given the partial sum sequence S_0, S_1, ..., S_N,
// applies the Shanks transformation repeatedly via the epsilon table.
//
// The ε-table is built column by column:
//   col[-1] = 0  (auxiliary)
//   col[ 0] = S_n  (partial sums)
//   col[ k] = col[k-2][n+1] + 1 / (col[k-1][n+1] - col[k-1][n])
//
// Even-numbered columns (ε_0, ε_2, ε_4, ...) are the Shanks estimates.
// The last even-column element is returned as the accelerated estimate.
// Reference: Wynn (1956), Math. Comp. 10, 91.
// ============================================================
std::complex<double> ElasticSolver::WynnEpsilon(
    const std::vector<std::complex<double>>& S)
{
    // Faithful line-by-line port of Ptolemy's EPSLON subroutine
    // (source_annotated.f lines 13330-13524).  Uses the same anti-diagonal
    // walking, singular-rule, and loss-of-significance checks.
    //
    // Input:  S[0..N-1] = partial sums of the series.
    // Output: accelerated estimate of the series limit.
    //
    // The Fortran uses 1-based indexing; comments note [F:...] for
    // the original 1-based index.  C++ uses 0-based throughout.

    int N = (int)S.size();
    if (N < 5) return N > 0 ? S[N-1] : std::complex<double>(0,0);

    // Approximate absolute value (matches Fortran APXABS)
    auto apx = [](std::complex<double> z) -> double {
        return std::abs(z.real()) + std::abs(z.imag());
    };

    const double BIG = 1e38;
    const double SMLNUM = 1e-300;
    const double DEPS = 1e-5;
    const std::complex<double> CBIG(BIG, 0.0);

    // Work on a mutable copy (Fortran overwrites X in-place)
    std::vector<std::complex<double>> X(S);
    int n = N;
    double acc = std::max(1e-8, DEPS * DEPS);  // tighter first pass

    // Outer convergence loop (Fortran label 41)
    for (;;) {
        if (n <= 0) return S[N-1];
        std::complex<double> fin = X[n-1];
        if (n == 1) return fin;

        // Check convergence
        double derr = 0;
        for (int i = 0; i < n; i++)
            derr = std::max(derr, apx(fin - X[i]));
        derr /= (apx(fin) + SMLNUM);
        if (derr <= acc) return fin;
        acc = DEPS;  // relax for subsequent passes
        if (n < 6) return fin;

        // --- Initialize anti-diagonal (Fortran label 2) ---
        bool ISW1 = false, ISW2 = false;

        // W1 = 1/(X[3]-X[2])  [F: 1/(X(4)-X(3))]
        std::complex<double> W1 = CBIG;
        { auto d = X[3]-X[2]; if (apx(d)!=0) W1 = 1.0/d; }

        // W5 = 1/(X[1]-X[0])  [F: 1/(X(2)-X(1))]
        std::complex<double> W5 = CBIG;
        { auto d = X[1]-X[0]; if (apx(d)!=0) W5 = 1.0/d; }

        // W4 = 1/(X[2]-X[1])  [F: 1/(X(3)-X(2))]
        std::complex<double> W4, T, W2, W3;
        { auto d = X[2]-X[1];
          if (apx(d) == 0) {
              W4 = CBIG; T = X[1]; W2 = X[2]; W3 = CBIG;
          } else {
              W4 = 1.0/d;
              T = CBIG;
              { auto dd = W4-W5; if (apx(dd)!=0) T = X[1] + 1.0/dd; }
              auto diff = W1 - W4;
              if (apx(diff) == 0) {
                  W2 = CBIG;
                  ISW2 = (T.real() != BIG);
                  W3 = W4;
              } else {
                  W2 = X[2] + 1.0/diff;
                  auto dd = W2 - T;
                  if (apx(dd) == 0) { W3 = CBIG; }
                  else              { W3 = W4 + 1.0/dd; }
              }
          }
        }

        ISW1 = ISW2;
        ISW2 = false;
        int IMIN = 3;  // 0-based (Fortran IMIN=4)
        bool truncated = false;

        // --- Main loop (Fortran label 40: DO I=5,N) ---
        for (int i = 4; i < n; i++) {
            int iaus = i - IMIN - 1;  // output index [F: I-IMIN]
            std::complex<double> W6;
            W4 = CBIG;
            W5 = X[i-1];
            auto dX = X[i] - X[i-1];

            enum Action { STORE_W2, COMPUTE_28, TRUNCATE, RESET_IMIN };
            Action action = STORE_W2;  // default path

            if (apx(dX) != 0) {
                W4 = 1.0 / dX;
                if (W1.real() == BIG) {
                    // goto 25: W6=BIG, goto 26
                    W6 = CBIG;
                    action = STORE_W2;
                } else {
                    W6 = W4 - W1;
                    // Singular rule test
                    if (apx(W6) <= 1e-12 * apx(W4)) {
                        ISW2 = true;
                        if (apx(W6) == 0) {
                            W5 = CBIG;
                            W6 = W1;
                            if (W2.real() != BIG) { action = COMPUTE_28; }
                            else                  { action = STORE_W2; ISW2 = false; }
                        } else {
                            // Fall through to W5 = X[i-1]+1/W6
                            goto compute_w5;
                        }
                    } else {
                        compute_w5:
                        W5 = X[i-1] + 1.0 / W6;
                        // Loss of significance
                        if (apx(W5) < 1e-10 * apx(X[i-1])) {
                            if (apx(W5) != 0) { n = iaus; truncated = true; break; }
                            // W5 == 0: fall through to STORE_W2 path
                        }
                        // goto 24
                        auto dd = W5 - W2;
                        if (apx(dd) != 0) {
                            W6 = W1 + 1.0 / dd;
                            action = COMPUTE_28;
                        } else {
                            W6 = CBIG;
                            action = STORE_W2;
                            ISW2 = false;
                        }
                    }
                }
            } else {
                // dX == 0: goto 24 with W4=BIG, W5=X[i-1]
                auto dd = W5 - W2;
                if (apx(dd) != 0) {
                    W6 = W1 + 1.0 / dd;
                    action = COMPUTE_28;
                } else {
                    W6 = CBIG;
                    action = STORE_W2;
                    ISW2 = false;
                }
            }

            // Execute chosen action
            if (action == STORE_W2) {
                X[iaus] = W2;
            } else if (action == COMPUTE_28) {
                if (ISW1) {
                    // Singular rule (Fortran labels 28→32 or 31)
                    if (W2.real() == BIG) {
                        X[iaus] = W5 + T - X[i-2];
                    } else {
                        auto denom1 = W2 - W5;
                        auto denom2 = W2 - T;
                        auto denom3 = X[i-2] - W2;
                        // Guard against division by zero in singular rule
                        if (apx(denom1) == 0 || apx(denom2) == 0 || apx(denom3) == 0) {
                            X[iaus] = CBIG; IMIN = i;
                            goto next_iter;
                        }
                        auto w7 = W5/denom1 + T/denom2 + X[i-2]/denom3;
                        if (apx(w7 + 1.0) == 0) {
                            X[iaus] = CBIG; IMIN = i;
                            goto next_iter;
                        }
                        X[iaus] = w7 * W2 / (1.0 + w7);
                    }
                } else {
                    // Normal rule (Fortran label 33)
                    auto dd = W6 - W3;
                    if (apx(dd) == 0) {
                        X[iaus] = CBIG; IMIN = i;
                        goto next_iter;
                    }
                    X[iaus] = W2 + 1.0 / dd;
                    // Second loss-of-significance
                    if (apx(X[iaus]) < 1e-10 * apx(W2)) {
                        if (apx(X[iaus]) != 0) { n = iaus; truncated = true; break; }
                    }
                }
            }

            // label 37: check for BIG W2 → reset IMIN
            if (W2.real() == BIG) { X[iaus] = CBIG; IMIN = i; }

            next_iter:
            W1 = W4;
            T  = W2;
            W2 = W5;
            W3 = W6;
            ISW1 = ISW2;
            ISW2 = false;
        }

        // Fortran: N = N - IMIN (1-based), our IMIN is 0-based
        if (!truncated)
            n = n - (IMIN + 1);
        // Loop back to convergence check (label 41)
    }
}


// GMatrix element G(v, v0, L) = Σ_J CG(L,v0-v; S,v | J,v0) * CG(L,0; S,v0 | J,v0) * S_J(L)
// Then nuclear amp f_N(v,v0,θ) = (1/2ik) Σ_L (2L+1) * (-1)^(v0-v)
//   * sqrt(fact(L-|v0-v|)/fact(L+|v0-v|)) * P_L^|v0-v|(cosθ) * e^{2iσ_L} * G(v,v0,L)
// (Raphael NuclearScatteringAmp)
// When useWynn_=true: build partial sum sequence and apply Wynn epsilon extrapolation.
std::complex<double> ElasticSolver::NuclearAmp(double v, double v0, double theta_deg) const {
    double theta  = theta_deg * PI / 180.0;
    double cos_th = std::cos(theta);
    int    m      = (int)std::abs(v0 - v);  // |Δm|

    std::complex<double> sum(0,0);
    // Wynn epsilon: collect partial sums starting from L=m (first non-zero term).
    // Only record after each contributing term — skip zero-contribution L<m terms.
    std::vector<std::complex<double>> partial_sums;
    if (useWynn_) partial_sums.reserve(Lmax_ - m + 1);

    for (int L = 0; L <= Lmax_; ++L) {
        if (m > L) continue;  // P_L^m = 0 for m>L, skip entirely (don't pollute Wynn sequence)

        // Compute GMatrix element
        std::complex<double> Gmat(0,0);
        if (S_ == 0.0) {
            Gmat = Smat_[L][0] - (v == v0 ? 1.0 : 0.0);
        } else {
            int nJ = (int)(2*S_) + 1;
            for (int idx = 0; idx < nJ; ++idx) {
                double J = L - S_ + idx;
                if (J < 0) continue;  // skip unphysical channels (J must be non-negative)
                double cg1 = CG(L, v0 - v, S_, v,  J, v0);
                double cg2 = CG(L, 0.0,    S_, v0,  J, v0);
                if (std::isnan(cg1) || std::isnan(cg2)) continue;
                Gmat += cg1 * cg2 * Smat_[L][idx];
            }
            if (v == v0) Gmat -= 1.0;
        }

        double sL = CoulPhase_[L];
        std::complex<double> e2s(std::cos(2*sL), std::sin(2*sL));

        // sign = (-1)^(v0-v): use signed dv (matches Raphael convention)
        int dv_int = (int)std::round(v0 - v);
        double sign = (dv_int % 2 == 0) ? 1.0 : -1.0;  // (-1)^(v0-v), signed
        double fact_ratio = 1.0;
        for (int k = L - m + 1; k <= L + m; ++k) fact_ratio *= k;
        double norm = sign / std::sqrt(fact_ratio);

        double PLm = AssocLegendreP(L, m, cos_th);
        sum += (double)(2*L + 1) * norm * PLm * e2s * Gmat;

        if (useWynn_) partial_sums.push_back(sum / std::complex<double>(0, 2.0*k_));
    }

    if (useWynn_ && (int)partial_sums.size() >= 6) {
        return WynnEpsilon(partial_sums);
    }
    return sum / std::complex<double>(0, 2.0 * k_);
}

// ── Associated Legendre P_L^1(cosθ) = -sinθ * dP_L/d(cosθ) ──────────────────
// Recursion: P^1_0=0, P^1_1=-sinθ, (L-1)*P^1_L=(2L-1)*cosθ*P^1_{L-1}-L*P^1_{L-2}
static std::vector<double> PlmAssoc1(double cosT, int Lmax) {
    double sinT = std::sqrt(std::max(0.0, (1.0-cosT)*(1.0+cosT)));
    std::vector<double> P1(Lmax+1, 0.0);
    if (Lmax >= 1) P1[1] = -sinT;
    for (int L = 2; L <= Lmax; ++L)
        P1[L] = ((2.0*L-1.0)*cosT*P1[L-1] - (double)L*P1[L-2]) / (double)(L-1);
    return P1;
}

double ElasticSolver::DCSUnpolarized(double theta_deg) const {
    if (theta_deg < 0.01) return std::numeric_limits<double>::quiet_NaN();
    std::complex<double> fC = CoulombAmp(theta_deg);

    // Handbook §2.9 general formula for all spins (Eq. 2.9, §2.9.7 spin-1 case):
    // dσ/dΩ = (1/(2S+1)) Σ_{ν,ν0} |fc·δ(ν,ν0) + f(ν,ν0)|² × 10
    // For spin-1: all 9 terms f(1,1), f(1,0), f(1,-1), f(0,1), f(0,0), ...
    // f(ν,ν0) = (1/2ik) Σ_L √(4π(2L+1)) e^{2iσ_L} g_{νν0}(L)
    // g_{νν0}(L) = Y_{L,ν0-ν} Σ_J CG(L,0;S,ν0|J,ν0)·CG(L,ν0-ν;S,ν|J,ν0)·S^J_L - δ_{νν0}
    // This is implemented by NuclearAmp(v,v0,θ) using CG coefficients.
    double total = 0.0;
    int nSpin = (int)(2*S_) + 1;
    for (int iv = 0; iv < nSpin; ++iv) {
        double v = -S_ + iv;
        for (int iv0 = 0; iv0 < nSpin; ++iv0) {
            double v0 = -S_ + iv0;
            std::complex<double> fN = NuclearAmp(v, v0, theta_deg);
            std::complex<double> amp = fN + (v == v0 ? fC : std::complex<double>(0,0));
            total += std::norm(amp);
        }
    }
    total /= (2*S_ + 1);   // average over initial spin projections
    return total * 10.0;   // fm² → mb/sr
}


// ============================================================
// Output
// ============================================================
void ElasticSolver::PrintSMatrix(int Lmax_print) const {
    if (Lmax_print < 0) Lmax_print = Lmax_;
    int nJ = (int)(2*S_) + 1;
    std::cout << std::setw(4) << "L";
    for (int idx = 0; idx < nJ; ++idx) {
        double Jrel = -S_ + idx;  // J - L
        std::string lbl = (Jrel >= 0 ? "J=L+" : "J=L") + std::to_string((int)std::abs(Jrel));
        std::cout << std::setw(13) << ("Re("+lbl+")") << std::setw(13) << ("Im("+lbl+")");
    }
    std::cout << "\n" << std::string(4 + 26*nJ, '-') << "\n";
    for (int L = 0; L <= Lmax_print; ++L) {
        std::cout << std::setw(4) << L << std::scientific << std::setprecision(4);
        for (int idx = 0; idx < nJ; ++idx) {
            double J = L - S_ + idx;
            if (J < 0) { std::cout << std::setw(13) << "---" << std::setw(13) << "---"; continue; }
            std::cout << std::setw(13) << Smat_[L][idx].real()
                      << std::setw(13) << Smat_[L][idx].imag();
        }
        std::cout << "\n";
    }
}

void ElasticSolver::SaveSMatrix(const std::string& fname) const {
    std::ofstream f(fname);
    int nJ = (int)(2*S_) + 1;
    f << "# Elastic S-matrix  (S=" << S_ << ", " << nJ << " J-channels per L)\n";
    f << "# L";
    for (int idx = 0; idx < nJ; ++idx) {
        double Jrel = -S_ + idx;
        std::string lbl = (Jrel >= 0 ? "J=L+" : "J=L") + std::to_string((int)std::abs(Jrel));
        f << "  Re(" << lbl << ")  Im(" << lbl << ")";
    }
    f << "\n";
    for (int L = 0; L <= Lmax_; ++L) {
        f << L;
        for (int idx = 0; idx < nJ; ++idx) {
            double J = L - S_ + idx;
            if (J < 0) { f << "  0.0  0.0"; continue; }
            f << " " << std::scientific << std::setprecision(8)
              << Smat_[L][idx].real() << " " << Smat_[L][idx].imag();
        }
        f << "\n";
    }
}

void ElasticSolver::SaveDCS(const std::string& fname,
                             double th_min, double th_max, double dth) const {
    std::ofstream f(fname);
    f << "# Elastic DCS (mb/sr)\n# theta_deg  dcs_mb_sr\n";
    for (double th = th_min; th <= th_max + 1e-9; th += dth) {
        double dcs = DCSUnpolarized(th);
        if (dcs > 0)
            f << th << " " << std::scientific << std::setprecision(6) << dcs << "\n";
    }
}
