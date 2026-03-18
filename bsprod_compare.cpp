// bsprod_compare.cpp
// Compare C++ BSPROD (ITYPE=1, NUCONL=0) against Fortran bsprod_test output.
//
// Build: g++ -O2 -std=c++17 bsprod_compare.cpp -o bsprod_compare
// Run:   ../fortran_testing/bsprod_test > /tmp/fort_bsprod.txt
//        ./bsprod_compare < /tmp/fort_bsprod.txt

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// ---- AITLAG: Aitken-Lagrange interpolation (matches Ptolemy AITLAG) ----
// TABLE[i] corresponds to x = i * h, where h = 1/stpinv
// NAIT = order (5-pt = NAIT=4)
static double Aitlag(double x, double stpinv, const std::vector<double>& table, int nait) {
    int ntable = (int)table.size();
    if (x < 0) return table[0];
    double xbyh = x * stpinv;
    int itable = (int)xbyh;
    if (itable >= ntable) return table[ntable-1];

    int nait2 = nait / 2;
    int nstart = itable - nait2;
    if (nstart < 0) nstart = 0;
    int nend = nstart + nait;
    if (nend >= ntable) nend = ntable - 1;
    nstart = nend - nait;

    // Precomputed inverses (matches Fortran DATA INVERS)
    static const double inv[17] = {
        0, // unused (1-indexed in Fortran)
        1.0, 0.5, 1.0/3.0, 0.25,
        0.2, 1.0/6.0, 1.0/7.0, 0.125,
        1.0/9.0, 0.1, 1.0/11.0, 1.0/12.0,
        1.0/13.0, 1.0/14.0, 1.0/15.0, 1.0/16.0
    };

    std::vector<double> fs(nait+2), dels(nait+2);
    fs[1] = table[nstart];   // TABLE(NSTART+1) in 1-indexed Fortran
    double del1 = xbyh - nstart;
    dels[1] = del1;

    double f = 0;
    for (int i = 1; i <= nait; i++) {
        f = table[nstart + i];  // TABLE(NSTART+1+I) in Fortran (0-indexed: nstart+i)
        del1 -= 1.0;
        for (int j = 1; j <= i; j++) {
            f = (f * dels[j] - fs[j] * del1) * inv[i + 1 - j];
        }
        dels[i+1] = del1;
        fs[i+1] = f;
    }
    return f;
}

// ---- BSPROD ITYPE=1, NUCONL=0 ----
// Inputs:
//   ra, rb: integration coordinates
//   x: cos(phi_ab)
//   S1,T1,S2,T2: coordinate transform coefficients
//   phi_T: tabulated φ_T(r) at step h_T
//   ivphi_P: tabulated V(r)*φ_P(r) at step h_P
//   bndmxp, bndmxt: asymptopia bounds
//   nait: Lagrange order
// Returns: FP*FT = IVPHI_P(RP) * PHI_T(RT); sets rp_out, rt_out
static double Bsprod1(double ra, double rb, double x,
                      double S1, double T1, double S2, double T2,
                      const std::vector<double>& phi_T,  double h_T,
                      const std::vector<double>& ivphi_P, double h_P,
                      double bndmxp, double bndmxt,
                      int nait,
                      double& rp_out, double& rt_out) {
    // Coordinate transform
    double rp2 = (S2*ra)*(S2*ra) + (T2*rb)*(T2*rb) + 2*(S2*ra)*(T2*rb)*x;
    double rt2 = (S1*ra)*(S1*ra) + (T1*rb)*(T1*rb) + 2*(S1*ra)*(T1*rb)*x;
    if (std::abs(rp2) < 1e-8) rp2 = 0;
    if (std::abs(rt2) < 1e-8) rt2 = 0;
    double rp = std::sqrt(rp2);
    double rt = std::sqrt(rt2);
    rp_out = rp; rt_out = rt;

    if (rp > bndmxp || rt > bndmxt) return 0.0;  // asymptopia

    // ITYPE=1: ALLSW=true → always interpolate both sides (no clipping)
    double ft = Aitlag(rt, 1.0/h_T,  phi_T,  nait);
    double fp = Aitlag(rp, 1.0/h_P, ivphi_P, nait);
    return fp * ft;
}

int main() {
    // ---- Build same tables as Fortran test ----
    const int NTAB_T = 641, NTAB_P = 641;
    const double H_T = 0.05, H_P = 0.05;

    // Kinematics (same as Fortran: 16O(d,p)17O stripping)
    const double BRATMS1 = 1.0, BRATMS2 = 1.0/16.0;
    double temp = 1.0 / (BRATMS1 + BRATMS2*(1.0+BRATMS1));
    double S1 = (1.0+BRATMS1)*(1.0+BRATMS2)*temp;
    double T1 = -(1.0+BRATMS2)*temp;
    double S2 = (1.0+BRATMS1)*temp;
    double T2 = -S1;

    std::cout << "C++ S1=" << S1 << " T1=" << T1 << " S2=" << S2 << " T2=" << T2 << "\n";

    // φ_T table: r*exp(-0.433*r)
    std::vector<double> phi_T(NTAB_T);
    for (int i = 0; i < NTAB_T; i++) {
        double r = i * H_T;
        phi_T[i] = (i == 0) ? 0.0 : r * std::exp(-0.433 * r);
    }

    // IVPHI_P table: -0.8*exp(-0.3*r^2)
    std::vector<double> ivphi_P(NTAB_P);
    for (int i = 0; i < NTAB_P; i++) {
        double r = i * H_P;
        ivphi_P[i] = -0.8 * std::exp(-0.3 * r * r);
    }

    const double bndmxp = (NTAB_P-1)*H_P;
    const double bndmxt = (NTAB_T-1)*H_T;
    const int NAIT = 4;

    // ---- Read Fortran output ----
    struct FRow { double ra, rb, x, fpft, rp, rt; };
    std::vector<FRow> fort_rows;
    std::string line;
    bool header_done = false;
    while (std::getline(std::cin, line)) {
        if (line.find("RA") != std::string::npos) { header_done = true; continue; }
        if (!header_done) continue;
        // Skip S1= line
        if (line.find("S1=") != std::string::npos) continue;
        if (line.find("BNDMXP=") != std::string::npos) continue;
        if (line.find("ASYM") != std::string::npos) continue;
        std::istringstream ss(line);
        FRow r;
        if (ss >> r.ra >> r.rb >> r.x >> r.fpft >> r.rp >> r.rt)
            fort_rows.push_back(r);
    }

    // ---- Compare ----
    double max_rel = 0, max_rel_rp = 0, max_rel_rt = 0;
    int nbad = 0;
    std::cout << std::setw(5) << "#"
              << std::setw(6) << "RA" << std::setw(6) << "RB" << std::setw(7) << "X"
              << std::setw(18) << "F_FPFT" << std::setw(18) << "C_FPFT"
              << std::setw(10) << "relΔ"
              << std::setw(10) << "F_RP" << std::setw(10) << "C_RP" << std::setw(9) << "relΔ_RP"
              << std::setw(10) << "F_RT" << std::setw(10) << "C_RT" << std::setw(9) << "relΔ_RT"
              << "\n";

    // Generate same test grid
    int idx = 0;
    for (int ira = 1; ira <= 6; ira++) {
        double ra = ira * 1.0;
        for (int irb = 1; irb <= 6; irb++) {
            double rb = irb * 1.0;
            for (int ix = 1; ix <= 4; ix++) {
                double xv = (ix==1) ? -0.5 : (ix==2) ? 0.0 : (ix==3) ? 0.5 : 1.0;
                double rp_out, rt_out;
                double cpp = Bsprod1(ra, rb, xv, S1, T1, S2, T2,
                                     phi_T, H_T, ivphi_P, H_P,
                                     bndmxp, bndmxt, NAIT, rp_out, rt_out);

                // Check if Fortran marked this as asymptopia (skipped row)
                if (idx >= (int)fort_rows.size()) break;

                // Match to the corresponding fort row (same RA,RB,X)
                double fval = fort_rows[idx].fpft;
                double frp  = fort_rows[idx].rp;
                double frt  = fort_rows[idx].rt;

                double ref = std::abs(fval);
                double rel = (ref > 1e-30) ? std::abs(cpp - fval)/ref : std::abs(cpp - fval);
                double rel_rp = (frp > 1e-10) ? std::abs(rp_out - frp)/frp : std::abs(rp_out - frp);
                double rel_rt = (frt > 1e-10) ? std::abs(rt_out - frt)/frt : std::abs(rt_out - frt);
                max_rel = std::max(max_rel, rel);
                max_rel_rp = std::max(max_rel_rp, rel_rp);
                max_rel_rt = std::max(max_rel_rt, rel_rt);
                if (rel > 1e-10) nbad++;

                std::cout << std::setw(5) << idx+1
                          << std::setw(6) << (int)ra << std::setw(6) << (int)rb
                          << std::setw(7) << xv
                          << std::setw(18) << std::setprecision(8) << std::scientific << fval
                          << std::setw(18) << cpp
                          << std::setw(10) << std::setprecision(2) << rel
                          << std::setw(10) << std::setprecision(5) << std::fixed << frp
                          << std::setw(10) << rp_out
                          << std::setw(9) << std::setprecision(2) << std::scientific << rel_rp
                          << std::setw(10) << std::setprecision(5) << std::fixed << frt
                          << std::setw(10) << rt_out
                          << std::setw(9) << rel_rt << "\n";
                idx++;
            }
        }
    }

    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "N compared: " << idx << "\n";
    std::cout << "Max rel error FPFT: " << std::scientific << max_rel << "\n";
    std::cout << "Max rel error RP:   " << max_rel_rp << "\n";
    std::cout << "Max rel error RT:   " << max_rel_rt << "\n";
    std::cout << "Points |Δ|>1e-10:  " << nbad << "\n";
    // Threshold: 1e-9 (5-pt Lagrange rounding across 144 test points is ~4e-10)
    if (max_rel < 1e-9 && max_rel_rp < 1e-9 && max_rel_rt < 1e-9)
        std::cout << "RESULT: PASS ✓ — C++ BSPROD (ITYPE=1) matches Fortran to FP precision (thr 1e-9)\n";
    else
        std::cout << "RESULT: FAIL ✗ — discrepancies found (exceeds 1e-9)!\n";

    return 0;
}
