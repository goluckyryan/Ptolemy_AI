// radint_test.cpp
// Radial integral for 16O(d,p)17O at 10 MeV/u
// I(L1,J1,L2,J2) = ∫ chi_out(r/mF) * phi_bs(r) * chi_in(r) dr * e^(iσ_L1) * e^(iσ_L2) * mF
//
// Loads bound state from /tmp/cpp_wf_d52.txt (precomputed by boundstate_test)
// Grid h=0.1 fm matches Raphael's dr=0.1 nStep=300 rStart=0
//
// ROOT CAUSES FIXED:
//  1. Coulomb phase: elastic.cpp CoulombPhase() Stirling formula had sign error on log term.
//     Fixed: theta_N0 = (x-0.5)*atan2(y,x) + y/2*ln(x²+y²) - y  [was: -y/2*ln - y*log(2π)/2]
//  2. Grid alignment: bs[i] = phi(r=i*h) for i=0..N-1 (r=0..29.9), matching wfI/wfO indexing.
//  3. Coordinate transform: Raphael evaluates chi_O at r/mF where mF=A_B/A_A=17/16.
//     Implemented via interp(rpos*mF, wfO)(rpos) which gives chi_O at rpos/mF.
//  4. No conjugate on chi_O: Raphael uses wf1*wf2 (NOT conj(wf2)), consistent with
//     prior-form DWBA where chi_f^(-*)=chi_f^(+)=wfO.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "include/elastic.h"

static const double HBARC = 197.32697, AMU = 931.494061;

// ── Build incoming d+16O solver ──────────────────────────────
static ElasticSolver BuildIncoming(int Lmax) {
    double m_d=2.01410, m_A=15.99491;
    double mu_amu = m_d*m_A/(m_d+m_A), Ecm=20.0*m_A/(m_d+m_A);
    double mu_MeV=mu_amu*AMU, k=std::sqrt(2.0*mu_MeV*Ecm)/HBARC;
    double eta=1.0*8.0/137.035999*mu_MeV/(HBARC*k);
    ElasticSolver s; s.SetTarget(16,8); s.SetProjectile(2,1);
    s.SetKinematics(k,eta,mu_amu); s.SetLmax(Lmax); s.SetGrid(0.1,30.0);
    s.AddVolumeWS({-88.955, 0.0},   1.149, 0.751);
    s.AddVolumeWS({ 0.0,   -2.348}, 1.345, 0.603);
    s.AddSurfaceWS({0.0,  -10.218}, 1.397, 0.687);
    s.AddSpinOrbit({-3.557, 0.0},   0.972, 1.011);
    s.AddCoulomb(1.303);
    return s;
}

// ── Build outgoing p+17O solver ──────────────────────────────
static ElasticSolver BuildOutgoing(int Lmax) {
    double m_p=1.00728, m_B=16.99913;
    double mu_amu=m_p*m_B/(m_p+m_B), Elab=20.85054, Ecm=Elab*m_B/(m_p+m_B);
    double mu_MeV=mu_amu*AMU, k=std::sqrt(2.0*mu_MeV*Ecm)/HBARC;
    double eta=1.0*8.0/137.035999*mu_MeV/(HBARC*k);
    ElasticSolver s; s.SetTarget(17,8); s.SetProjectile(1,1);
    s.SetKinematics(k,eta,mu_amu); s.SetLmax(Lmax); s.SetGrid(0.1,30.0);
    s.AddVolumeWS({-49.544, 0.0},   1.250, 0.650);
    s.AddVolumeWS({ 0.0,   -2.061}, 1.250, 0.650);
    s.AddSurfaceWS({0.0,   -7.670}, 1.250, 0.650);
    s.AddSpinOrbit({-5.296, 0.0},   1.250, 0.650);
    s.AddCoulomb(1.419);
    return s;
}

// ── Load bound state from file, subsample to target grid step ──
// File format: r(fm)  phi=u/r   (h_file=0.01 fm, h_target=0.1 fm)
// Returns bs[i] = phi(r=i*h_target) for i=0..N_target-1
// (index 0 → r=0, aligned with Raphael's rpos = [0, 0.1, ..., 29.9])
static std::vector<double> LoadBoundState(const std::string& fname,
                                          double h_target, int N_target) {
    std::ifstream f(fname);
    if (!f) { std::cerr << "Cannot open " << fname << std::endl; return {}; }
    std::vector<double> r_raw, phi_raw;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream ss(line);
        double r, phi; ss >> r >> phi;
        r_raw.push_back(r); phi_raw.push_back(phi);
    }
    double h_raw = r_raw[0];  // e.g. 0.01 fm

    // bs[0] = phi(r=0) = 0  (u(0)=0, phi=u/r → 0 as r→0)
    // bs[i] = phi(r=i*h_target) for i=1..N_target-1
    std::vector<double> bs(N_target, 0.0);
    bs[0] = 0.0;
    for (int i = 1; i < N_target; i++) {
        double r_target = i * h_target;
        int idx = (int)std::round(r_target / h_raw) - 1;
        idx = std::max(0, std::min(idx, (int)r_raw.size()-1));
        bs[i] = phi_raw[idx];
    }
    return bs;
}

// ── Linear interpolation of complex wavefunction ──────────────
// Interpolate wf (defined on grid r[i]=i*h, i=0..n-1) at position r_query.
// Mirrors Raphael's interp1d(rpos*mF, wfO)(rpos) → chi_O at rpos/mF.
static std::complex<double> InterpWF(
    const std::vector<std::complex<double>>& wf, double h, double r_query) {
    if (r_query <= 0.0) return wf[0];
    double fidx = r_query / h;
    int i0 = (int)fidx;
    int i1 = i0 + 1;
    if (i0 < 0) return wf[0];
    if (i0 >= (int)wf.size()) return {0.0, 0.0};
    if (i1 >= (int)wf.size()) return wf[i0];
    double t = fidx - i0;
    return wf[i0] * (1.0 - t) + wf[i1] * t;
}

// ── Simpson's rule (complex) ──────────────────────────────────
static std::complex<double> Simpson(
    const std::vector<std::complex<double>>& f, double dx, int n) {
    if (n < 3) return 0;
    int m = (n%2==1) ? n : n-1;
    std::complex<double> s(0,0);
    for (int i=0; i<m-2; i+=2) s += f[i] + 4.0*f[i+1] + f[i+2];
    s *= dx/3.0;
    if (n > m) s += (f[m-1]+f[n-1])*(dx/2.0);
    return s;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    const double h=0.1, rmax=30.0;
    const int N=(int)(rmax/h);  // 300 pts: r=0, 0.1, ..., 29.9 (matches Raphael rpos)
    const int maxL1=18, maxL2=20;

    // ZR coordinate transformation massFactor = A_B/A_A = 17/16
    // (since A_A=16 < A_B=17; and A_a_incoming=2 > A_a_outgoing=1 → Raphael scales outgoing grid)
    // Raphael: rpos_O_temp = rpos_I * mF; interp(rpos_O_temp, wfu_O)(r) = chi_O(r/mF)
    // Integral: ∫ chi_I(r) * chi_O(r/mF) * bs(r) dr * pf1 * pf2 * mF
    const double massFactor = 17.0 / 16.0;

    std::cout << "Solving incoming d+16O..." << std::flush;
    ElasticSolver dwI = BuildIncoming(maxL1); dwI.Solve();
    std::cout << " k=" << dwI.GetK() << " eta=" << dwI.GetEta() << std::endl;

    std::cout << "Solving outgoing p+17O..." << std::flush;
    ElasticSolver dwO = BuildOutgoing(maxL2); dwO.Solve();
    std::cout << " k=" << dwO.GetK() << " eta=" << dwO.GetEta() << std::endl;

    // Check Coulomb phases
    std::cout << "Coulomb phases (corrected):" << std::endl;
    std::cout << "  sigma_L=0 incoming = " << dwI.GetCoulombPhase(0)
              << " (expect -0.207079)" << std::endl;
    std::cout << "  sigma_L=0 outgoing = " << dwO.GetCoulombPhase(0)
              << " (expect -0.151647)" << std::endl;

    std::cout << "Loading bound state from /tmp/cpp_wf_d52.txt..." << std::flush;
    auto bs = LoadBoundState("/tmp/cpp_wf_d52.txt", h, N);
    if (bs.empty()) return 1;
    std::cout << " " << bs.size() << " pts (r=0.." << (N-1)*h << " fm)" << std::endl;

    // Norm check: ∫phi²r²dr ≈ 1
    double norm2=0;
    for(int i=1; i<N; i++) { double r=i*h; norm2+=bs[i]*bs[i]*r*r*h; }
    std::cout << "  Norm: ∫phi²r²dr = " << norm2 << " (expect ~1)" << std::endl;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\n=== Radial integrals (C++) vs Raphael ===" << std::endl;
    std::cout << "Format: {L1, J1, L2, J2, Re+ImI}" << std::endl;

    for (int dJ1 : {-1,0,1}) {
        for (int dJ2half : {-1,1}) {
            double dJ2 = dJ2half * 0.5;
            std::string tag1 = (dJ1==-1)?"L-1":(dJ1==0)?"L+0":"L+1";
            std::string tag2 = (dJ2half==-1)?"L-1/2":"L+1/2";
            std::cout << "\n======================= J1 = " << tag1
                      << ", J2 = " << tag2 << std::endl;

            for (int L1=0; L1<=maxL1; L1++) {
                double J1 = L1 + dJ1;
                if (J1 < 0 || J1 < std::abs(L1-1.0)) continue;
                int twoJ1 = (int)std::round(2*J1);
                const auto& wfI = dwI.GetWavefunction(L1, twoJ1);
                if (wfI.empty()) continue;
                double sL1 = dwI.GetCoulombPhase(L1);
                std::complex<double> pf1(std::cos(sL1), std::sin(sL1));

                std::cout << "{";
                for (int L2=std::max(0,L1-2); L2<=L1+2 && L2<=maxL2; L2++) {
                    double J2 = L2 + dJ2;
                    if (J2 < 0 || J2 < std::abs(L2-0.5)) {
                        std::cout << "{" << L1 << "," << J1 << "," << L2 << "," << J2 << ", 0.0},  ";
                        continue;
                    }
                    int twoJ2 = (int)std::round(2*J2);
                    const auto& wfO = dwO.GetWavefunction(L2, twoJ2);
                    if (wfO.empty()) {
                        std::cout << "{" << L1 << "," << J1 << "," << L2 << "," << J2 << ", empty},  ";
                        continue;
                    }
                    double sL2 = dwO.GetCoulombPhase(L2);
                    std::complex<double> pf2(std::cos(sL2), std::sin(sL2));

                    // Integration grid: r[i] = i*h, i=0..N-1 (matches Raphael rpos_I)
                    // Coordinate transform: chi_O evaluated at r/massFactor
                    //   = Raphael interp(rpos_I*mF, wfO)(rpos_I) = chi_O_orig(rpos_I/mF)
                    // No conjugate: Raphael uses wf1*wf2 (prior-form DWBA convention)
                    int npts = N;
                    std::vector<std::complex<double>> integrand(npts);
                    for (int i=0; i<npts; i++) {
                        double r = i * h;
                        double r_O = r / massFactor;  // chi_O at r/mF = r*16/17
                        std::complex<double> chi_I = (i < (int)wfI.size()) ? wfI[i] : std::complex<double>(0,0);
                        std::complex<double> chi_O = InterpWF(wfO, h, r_O);  // no conj
                        integrand[i] = chi_I * chi_O * bs[i];  // wf1 * wf2 (Raphael: no conj)
                    }

                    auto I = Simpson(integrand, h, npts) * pf1 * pf2 * massFactor;
                    std::cout << "{" << std::setw(2) << L1 << ", " << J1 << ", "
                              << std::setw(2) << L2 << ", " << J2 << ", "
                              << std::setw(10) << I.real() << "+"
                              << std::setw(10) << I.imag() << "I},  ";
                }
                std::cout << "}," << std::endl;
            }
        }
    }

    return 0;
}
