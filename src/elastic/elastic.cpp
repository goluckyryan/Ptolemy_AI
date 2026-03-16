// elastic.cpp — Unified elastic scattering solver implementation
// Sign conventions verified against Raphael (Python) for:
//   148Sm(α,α) 50 MeV  and  60Ni(p,p) 30 MeV
#include "elastic.h"
#include "rcwfn.h"
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
        {40,20,  37224.000},  // 40Ca
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
    double ma = NuclearMass(Ap_, Zp_);
    double mA = NuclearMass(At_, Zt_);
    double E_tot = std::sqrt((ma + mA)*(ma + mA) + 2.0*mA*Elab_);
    Ecm_   = E_tot - ma - mA;
    mu_    = ma * mA / (ma + mA);        // MeV/c^2
    k_     = std::sqrt(2.0*mu_*Ecm_) / HBARC;
    eta_   = Zp_ * Zt_ * E2 * k_ / (2.0 * Ecm_);
    f_conv_ = 2.0 * mu_ / (HBARC * HBARC);

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
    double theta_N0 = (x-0.5)*std::atan2(y, x) - y/2*std::log(x*x + y*y)
                    + y*(1.0 - std::log(2*PI)/2);
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
    int Zp, int Zt, int At,
    double h, int N,
    std::vector<double>& Vr, std::vector<double>& Wi,
    std::vector<double>& Vc,
    std::vector<double>& VsoRe, std::vector<double>& VsoIm)
{
    Vr.assign(N, 0); Wi.assign(N, 0); Vc.assign(N, 0);
    VsoRe.assign(N, 0); VsoIm.assign(N, 0);

    double At13 = std::cbrt((double)At);
    double Rc   = hasCoulomb ? rc0 * At13 : 0.0;

    for (int i = 0; i < N; ++i) {
        double r = (i == 0) ? 1e-10 : i * h;

        // Coulomb
        if (hasCoulomb && Rc > 0) {
            Vc[i] = (r < Rc)
                ? Zp * Zt * E2 * (3.0 - r*r/(Rc*Rc)) / (2.0*Rc)
                : Zp * Zt * E2 / r;
        }

        // Nuclear potentials
        for (const auto& p : pots) {
            double R0 = p.r0 * At13;
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
    const std::vector<double>& FC2, const std::vector<double>& GC2) const
{
    double k2 = k_*k_, h2_12 = h_*h_/12.0;
    double LL1 = (double)L*(L+1);
    int n1 = N_ - 4, n2 = N_ - 3;

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

    // Two-point matching: u = A*G_L + B*F_L → S = (B+iA)/(B-iA)
    double f1 = FC1[L], g1 = GC1[L], f2 = FC2[L], g2 = GC2[L];
    std::complex<double> u1 = u[n1], u2 = u[n2];
    double det = f2*g1 - f1*g2;
    if (std::abs(det) < 1e-30) return {1.0, 0.0};  // fallback: no interaction
    std::complex<double> A = (f2*u1 - u2*f1) / det;
    std::complex<double> B = (u2*g1 - g2*u1) / det;
    using namespace std::complex_literals;
    return (B + 1i*A) / (B - 1i*A);
}

// ============================================================
// Main S-matrix computation
// ============================================================
void ElasticSolver::CalcScatteringMatrix() {
    if (k_ <= 0) CalcKinematics();

    // Build potential arrays
    std::vector<double> Vr, Wi, Vc, VsoRe, VsoIm;
    BuildPotentialArrays(pots_, hasCoulomb_, rc0_, Zp_, Zt_, At_,
                         h_, N_, Vr, Wi, Vc, VsoRe, VsoIm);

    // Precompute Coulomb functions at both matching points (MINL=0 critical!)
    int n1 = N_ - 4, n2 = N_ - 3;
    double rho1 = k_ * n1 * h_, rho2 = k_ * n2 * h_;
    std::vector<double> FC1(Lmax_+2), FCP1(Lmax_+2), GC1(Lmax_+2), GCP1(Lmax_+2);
    std::vector<double> FC2(Lmax_+2), FCP2(Lmax_+2), GC2(Lmax_+2), GCP2(Lmax_+2);
    Rcwfn(rho1, eta_, 0, Lmax_, FC1, FCP1, GC1, GCP1);
    Rcwfn(rho2, eta_, 0, Lmax_, FC2, FCP2, GC2, GCP2);

    // Coulomb phases
    CoulPhase_.resize(Lmax_ + 1);
    for (int L = 0; L <= Lmax_; ++L)
        CoulPhase_[L] = CoulombPhase(L);

    // Allocate S-matrix
    Smat_.resize(Lmax_ + 1);
    for (auto& s : Smat_) s = {std::complex<double>(1,0), std::complex<double>(1,0)};

    bool hasSO = false;
    for (const auto& p : pots_) if (p.type == kSpinOrbit) { hasSO = true; break; }

    for (int L = 0; L <= Lmax_; ++L) {
        if (S_ == 0.0 || !hasSO) {
            // Spin-0 or no SO: single channel
            Smat_[L][0] = RunNumerov(L, 0.0, Vr, Wi, Vc, VsoRe, VsoIm,
                                     FC1, GC1, FC2, GC2);
            Smat_[L][1] = Smat_[L][0];
        } else {
            // Spin-1/2: two J channels per L
            for (int spin = 0; spin <= 1; ++spin) {
                if (L == 0 && spin == 1) { Smat_[0][1] = {1,0}; continue; }
                double J = (spin == 0) ? L + S_ : L - S_;
                double LS_val = (J*(J+1) - L*(L+1) - S_*(S_+1)) / 2.0;
                Smat_[L][spin] = RunNumerov(L, LS_val, Vr, Wi, Vc, VsoRe, VsoIm,
                                            FC1, GC1, FC2, GC2);
            }
        }
    }
}

// ============================================================
// DCS computation
// ============================================================
std::complex<double> ElasticSolver::CoulombAmp(double theta_deg) const {
    double theta = theta_deg * PI / 180.0;
    double sin2  = std::sin(theta / 2.0);
    double phi_C = 2.0*CoulPhase_[0] - 2.0*eta_*std::log(sin2);
    double fmod  = eta_ / (2.0 * k_ * sin2 * sin2);
    return std::complex<double>(-fmod*std::cos(phi_C), -fmod*std::sin(phi_C));
}

double ElasticSolver::DCSUnpolarized(double theta_deg) const {
    if (theta_deg < 0.01) return std::numeric_limits<double>::quiet_NaN();
    double theta  = theta_deg * PI / 180.0;
    double cos_th = std::cos(theta);
    double sin_th = std::sin(theta);

    std::complex<double> fC = CoulombAmp(theta_deg);

    // Nuclear amplitudes
    std::complex<double> A(0,0), B(0,0);
    std::complex<double> two_ik(0, 2.0*k_);

    for (int L = 0; L <= Lmax_; ++L) {
        double sL = CoulPhase_[L];
        std::complex<double> e2s(std::cos(2*sL), std::sin(2*sL));
        double PL  = LegendreP(L, cos_th);

        if (S_ == 0.0) {
            // Spin-0: single S-matrix
            A += (double)(2*L+1) * e2s * (Smat_[L][0] - 1.0) * PL;
        } else {
            // Spin-1/2: spin-averaged (McIntyre formula)
            // S_[L][0] = J=L+1/2, S_[L][1] = J=L-1/2 (for L>0)
            std::complex<double> sp = Smat_[L][0];
            std::complex<double> sm = (L == 0) ? std::complex<double>(1,0) : Smat_[L][1];
            double dPL = -sin_th * LegendrePprime(L, cos_th);
            A += e2s * ((double)(L+1)*(sp-1.0) + (double)L*(sm-1.0)) * PL;
            B += e2s * (sm - sp) * dPL;
        }
    }
    A /= two_ik;
    B /= two_ik;

    // dσ/dΩ = |f_C + A|² + |B|²  (fm²)
    double dcs_fm2 = std::norm(fC + A) + std::norm(B);
    return dcs_fm2 * 10.0;  // mb/sr
}

// ============================================================
// Output
// ============================================================
void ElasticSolver::PrintSMatrix(int Lmax_print) const {
    if (Lmax_print < 0) Lmax_print = Lmax_;
    if (S_ == 0.0) {
        std::cout << std::setw(4) << "L"
                  << std::setw(14) << "Re(S)" << std::setw(14) << "Im(S)"
                  << std::setw(12) << "|S|" << "\n";
        std::cout << std::string(44, '-') << "\n";
        for (int L = 0; L <= Lmax_print; ++L) {
            auto& s = Smat_[L][0];
            std::cout << std::setw(4) << L
                      << std::scientific << std::setprecision(5)
                      << std::setw(14) << s.real() << std::setw(14) << s.imag()
                      << std::fixed << std::setprecision(6)
                      << std::setw(12) << std::abs(s) << "\n";
        }
    } else {
        std::cout << std::setw(4) << "L"
                  << std::setw(13) << "Re(S+)" << std::setw(13) << "Im(S+)"
                  << std::setw(13) << "Re(S-)" << std::setw(13) << "Im(S-)" << "\n";
        std::cout << std::string(56, '-') << "\n";
        for (int L = 0; L <= Lmax_print; ++L) {
            std::cout << std::setw(4) << L
                      << std::scientific << std::setprecision(4)
                      << std::setw(13) << Smat_[L][0].real()
                      << std::setw(13) << Smat_[L][0].imag()
                      << std::setw(13) << Smat_[L][1].real()
                      << std::setw(13) << Smat_[L][1].imag() << "\n";
        }
    }
}

void ElasticSolver::SaveSMatrix(const std::string& fname) const {
    std::ofstream f(fname);
    f << "# Elastic S-matrix\n";
    if (S_ == 0.0) {
        f << "# L  Re(S)  Im(S)\n";
        for (int L = 0; L <= Lmax_; ++L)
            f << L << " " << std::scientific << std::setprecision(8)
              << Smat_[L][0].real() << " " << Smat_[L][0].imag() << "\n";
    } else {
        f << "# L  Re(S+)  Im(S+)  Re(S-)  Im(S-)\n";
        for (int L = 0; L <= Lmax_; ++L)
            f << L << " " << std::scientific << std::setprecision(8)
              << Smat_[L][0].real() << " " << Smat_[L][0].imag() << " "
              << Smat_[L][1].real() << " " << Smat_[L][1].imag() << "\n";
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
