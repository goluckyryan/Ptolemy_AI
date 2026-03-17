// chi_validate.cpp — Compare C++ distorted waves vs Ptolemy
// 16O(d,p)17O at Elab=20 MeV
// Uses ElasticSolver (public) to compute chi_a, chi_b

#include "elastic.h"
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>

int main() {
    const double AMU_MEV = 931.494061;
    const double HBARC   = 197.32697;

    // -------------------------------------------------------
    // Incoming: d + 16O, Elab=20 MeV, r0target convention
    // Ptolemy: Ecm=17.763 MeV, k=1.2328 fm^-1
    // Non-relativistic: mu = mA*ma/(mA+ma) AMU
    double mA=16.0, ma=2.0, mb=1.0, mB=17.0;
    double Ecm_in  = mA/(mA+ma)*20.0;          // 17.7778 MeV
    double mu_in   = mA*ma/(mA+ma);             // 1.7778 AMU
    double k_in    = std::sqrt(2*mu_in*AMU_MEV*Ecm_in)/HBARC;
    double eta_in  = 1.0*8.0*mu_in*AMU_MEV/(137.036*HBARC*k_in);

    // Outgoing: p + 17O, Ecm from Ptolemy output = 19.682 MeV
    double Ecm_out = 19.682;
    double mu_out  = mB*mb/(mB+mb);             // 0.9444 AMU
    double k_out   = std::sqrt(2*mu_out*AMU_MEV*Ecm_out)/HBARC;
    double eta_out = 1.0*8.0*mu_out*AMU_MEV/(137.036*HBARC*k_out);

    printf("=== Kinematics ===\n");
    printf("Incoming: Ecm=%.5f  k=%.6f  eta=%.6f\n", Ecm_in,  k_in,  eta_in);
    printf("Outgoing: Ecm=%.5f  k=%.6f  eta=%.6f\n\n", Ecm_out, k_out, eta_out);
    printf("Ptolemy:  Incoming k=1.2328, Outgoing k=0.94628\n\n");

    // -------------------------------------------------------
    // Test a range of step sizes to see which matches Ptolemy
    // Ptolemy uses h=0.125 fm; our C++ uses h=0.100 fm
    // We want to check both to see which gives the correct normalization

    auto interp_wf = [](const std::vector<std::complex<double>>& wf,
                        double h, double r) -> std::complex<double> {
        if (r < 0) return {0,0};
        double idx_f = r / h;
        int idx_i = (int)idx_f;
        double frac = idx_f - idx_i;
        if (idx_i + 1 < (int)wf.size())
            return wf[idx_i]*(1.0-frac) + wf[idx_i+1]*frac;
        if (idx_i < (int)wf.size()) return wf[idx_i];
        return {0,0};
    };

    // Ptolemy STARTS (h_pto=0.125 fm):
    // L=0 JP=2/2 INCOMING: idx=2 → r=0.125: (0.98146e-1, 0.37532e-1), idx=3 → r=0.250: (0.18354, 0.69800e-1)
    // L=2 JP=2/2 INCOMING: idx=2 → r=0.125: (0.86875e-3,-0.37850e-3), idx=3 → r=0.250: (0.67237e-2,-0.29356e-2)
    // L=0 JP=1/2 OUTGOING: idx=2 → r=0.125: (-0.10584,  0.76296e-2), idx=3 → r=0.250: (-0.20692, 0.15074e-1)
    // L=2 JP=3/2 OUTGOING: idx=2 → r=0.125: (-0.76411e-4, 0.23767e-3), idx=3→r=0.250: (-0.60544e-3, 0.18854e-2)
    // L=2 JP=5/2 OUTGOING: idx=2 → r=0.125: (-0.22404e-3, 0.22612e-3), idx=3→r=0.250: (-0.17732e-2, 0.17908e-2)
    // S-matrix (JP-basis from print=4000):
    // INCOMING L=0 JP=2/2:  0.14429 + 0.11977i
    // INCOMING L=2 JP=2/2:  0.10997 - 0.12204i
    // OUTGOING L=0 JP=1/2:  0.49629 - 0.01952i
    // OUTGOING L=2 JP=3/2: -0.32679 - 0.31757i
    // OUTGOING L=2 JP=5/2: -0.04084 - 0.54572i

    struct TestCase {
        const char* label;
        bool        isIncoming;
        int L, twoJP;
        double pto_SJR, pto_SJI;
        double pto_r1_Re, pto_r1_Im;  // r = 0.125 fm (idx=2 in Ptolemy)
        double pto_r2_Re, pto_r2_Im;  // r = 0.250 fm (idx=3 in Ptolemy)
    };

    std::vector<TestCase> cases = {
        {"INCOMING L=0 JP=2/2 d+16O", true,  0, 2,  0.14429,     0.11977,
          9.8146e-2,  3.7532e-2, 1.8354e-1,  6.9800e-2},
        {"INCOMING L=2 JP=2/2 d+16O", true,  2, 2,  0.10997,    -0.12204,
          8.6875e-4, -3.7850e-4, 6.7237e-3, -2.9356e-3},
        {"OUTGOING L=0 JP=1/2 p+17O", false, 0, 1,  0.49629,    -1.9521e-2,
         -1.0584e-1,  7.6296e-3,-2.0692e-1,  1.5074e-2},
        {"OUTGOING L=2 JP=3/2 p+17O", false, 2, 3, -0.32679,    -0.31757,
         -7.6411e-5,  2.3767e-4,-6.0544e-4,  1.8854e-3},
        {"OUTGOING L=2 JP=5/2 p+17O", false, 2, 5, -4.0840e-2,  -0.54572,
         -2.2404e-4,  2.2612e-4,-1.7732e-3,  1.7908e-3},
    };

    // Test with h=0.100 (our default) and h=0.125 (Ptolemy)
    for (double h_test : {0.100, 0.125}) {
        printf("\n========= h = %.3f fm =========\n\n", h_test);

        for (auto& tc : cases) {
            ElasticSolver s;
            if (tc.isIncoming) {
                s.SetProjectile(2, 1);  // deuteron
                s.SetTarget(16, 8);
                s.SetKinematics(k_in, eta_in, mu_in);
                // Incoming: d+16O potentials
                s.AddVolumeWS ({88.9546, 0.000}, 1.1489, 0.7508);
                s.AddVolumeWS ({0.000,   2.3480}, 1.3446, 0.6030);
                s.AddSurfaceWS({0.000, -10.2180}, 1.3943, 0.6872);
                s.AddSpinOrbit({7.1140,  0.0000}, 0.9720, 1.0110);
                s.AddCoulomb(1.3030);
            } else {
                s.SetProjectile(1, 1);  // proton
                s.SetTarget(17, 8);
                s.SetKinematics(k_out, eta_out, mu_out);
                // Outgoing: p+17O potentials
                s.AddVolumeWS ({49.5434, 0.000 }, 1.1462, 0.6753);
                s.AddVolumeWS ({0.000,   2.0611}, 1.1462, 0.6753);
                s.AddSurfaceWS({0.000,  -7.6703}, 1.3016, 0.5275);
                s.AddSpinOrbit({5.2956, -0.1059}, 0.9338, 0.5900);
                s.AddCoulomb(1.3030);
            }

            s.SetGrid(h_test, 30.0);
            s.SetLmax(tc.L);
            s.Solve();

            auto Sval = s.GetSMatrix(tc.L, tc.twoJP);
            const auto& wf = s.GetWavefunction(tc.L, tc.twoJP);

            printf("  [%s]\n", tc.label);
            printf("  S: C++=(%+.5f, %+.5f)  Pto=(%+.5f, %+.5f)  dRe=%+.5f dIm=%+.5f\n",
                   Sval.real(), Sval.imag(), tc.pto_SJR, tc.pto_SJI,
                   Sval.real()-tc.pto_SJR, Sval.imag()-tc.pto_SJI);

            // Compare wavefunction at r=0.125, 0.250 fm
            for (int step = 1; step <= 2; step++) {
                double r = step * 0.125;
                auto cpp_val = interp_wf(wf, h_test, r);
                double pto_Re = (step==1) ? tc.pto_r1_Re : tc.pto_r2_Re;
                double pto_Im = (step==1) ? tc.pto_r1_Im : tc.pto_r2_Im;
                double ratio_re = (std::abs(pto_Re)>1e-10) ? cpp_val.real()/pto_Re : 999.0;
                double ratio_im = (std::abs(pto_Im)>1e-10) ? cpp_val.imag()/pto_Im : 999.0;
                printf("    r=%5.3f: C++=(%+.5e, %+.5e)  Pto=(%+.5e, %+.5e)  rRe=%.4f rIm=%.4f\n",
                       r, cpp_val.real(), cpp_val.imag(), pto_Re, pto_Im, ratio_re, ratio_im);
            }

            // Print wf from r=0 to r=2 fm at C++ grid spacing
            printf("    C++ wf grid (h=%.3f), r=0..2 fm:\n", h_test);
            int nprint = (int)(2.0/h_test) + 1;
            for (int i = 0; i < nprint && i < (int)wf.size(); i++) {
                double r = i * h_test;
                printf("    %5.3f  %+.7e  %+.7e\n", r, wf[i].real(), wf[i].imag());
            }
            printf("\n");
        }
    }

    return 0;
}
