// chi_validate2.cpp — Use SetSystem instead of SetKinematics
// Try to reproduce Ptolemy S-matrix for d+16O and p+17O

#include "elastic.h"
#include <cmath>
#include <cstdio>
#include <complex>

int main() {
    // === INCOMING: d + 16O at Elab=20 MeV ===
    // Ptolemy INCOMING S-matrix (from print=4000):
    //   L=0 JP=2/2: 0.14429 + 0.11977i
    //   L=2 JP=2/2: 0.10997 - 0.12204i
    printf("=== INCOMING: d + 16O at Elab=20 MeV ===\n");
    printf("Ptolemy: k=1.2328, step=0.125 fm\n\n");

    for (double h : {0.100, 0.125}) {
        printf("-- h=%.3f --\n", h);
        ElasticSolver s;
        s.SetSystem(2, 1, 16, 8, 20.0);
        // Potentials from o16dp_gs_new.in INCOMING block
        s.AddVolumeWS ({88.9546, 0.0000}, 1.1489, 0.7508);
        s.AddVolumeWS ({0.0,     2.3480}, 1.3446, 0.6030);
        s.AddSurfaceWS({0.0,   -10.2180}, 1.3943, 0.6872);
        s.AddSpinOrbit({7.1140,  0.0000}, 0.9720, 1.0110);
        s.AddCoulomb(1.3030);
        s.SetGrid(h, 30.0);
        s.SetLmax(4);
        s.Solve();

        printf("  k=%.6f  eta=%.6f\n", s.GetK(), s.GetEta());

        for (int L=0; L<=4; L+=2) {
            for (int twoJ : {2*L-2, 2*L, 2*L+2}) {
                if (twoJ < 0) continue;
                auto S = s.GetSMatrix(L, twoJ);
                if (std::isnan(S.real())) continue;
                printf("  L=%d JP=%d/2: S=(%+.5f, %+.5f)\n", L, twoJ, S.real(), S.imag());
            }
        }
        printf("\n");
    }

    // === OUTGOING: p + 17O ===
    // Ptolemy: Ecm=19.682 MeV, k=0.94628
    // In SetSystem, Elab_p = what value gives Ecm_p = 19.682?
    // Ecm = At/(At+Ap) * Elab → Elab = Ecm * (At+Ap)/At = 19.682 * 18/17 = 20.848 MeV ✓
    // (Ptolemy printed: E LAB = 20.848 MeV)
    printf("\n=== OUTGOING: p + 17O at Elab=20.848 MeV ===\n");
    printf("Ptolemy: k=0.94628, step=0.125 fm\n");
    printf("Ptolemy S-matrix:\n");
    printf("  L=0 JP=1/2:  0.49629 - 0.01952i\n");
    printf("  L=2 JP=3/2: -0.32679 - 0.31757i\n");
    printf("  L=2 JP=5/2: -0.04084 - 0.54572i\n\n");

    for (double h : {0.100, 0.125}) {
        printf("-- h=%.3f --\n", h);
        ElasticSolver s;
        s.SetSystem(1, 1, 17, 8, 20.848);
        // Potentials from o16dp_gs_new.in OUTGOING block
        s.AddVolumeWS ({49.5434, 0.0000}, 1.1462, 0.6753);
        s.AddVolumeWS ({0.0,     2.0611}, 1.1462, 0.6753);
        s.AddSurfaceWS({0.0,    -7.6703}, 1.3016, 0.5275);
        s.AddSpinOrbit({5.2956, -0.1059}, 0.9338, 0.5900);
        s.AddCoulomb(1.3030);
        s.SetGrid(h, 30.0);
        s.SetLmax(4);
        s.Solve();

        printf("  k=%.6f  eta=%.6f\n", s.GetK(), s.GetEta());

        for (int L=0; L<=4; L+=2) {
            for (int twoJ : {2*L-1, 2*L+1}) {
                if (twoJ < 0) continue;
                auto S = s.GetSMatrix(L, twoJ);
                if (std::isnan(S.real())) continue;
                printf("  L=%d JP=%d/2: S=(%+.5f, %+.5f)\n", L, twoJ, S.real(), S.imag());
            }
        }
        printf("\n");
    }

    // Now print the wavefunction for L=0 to check normalization
    printf("\n=== Wavefunction check: L=0 incoming (d+16O) ===\n");
    printf("Ptolemy STARTS (h=0.125): idx=2 → r=0.125: (0.09815, 0.03753)\n");
    printf("                          idx=3 → r=0.250: (0.18354, 0.06980)\n\n");

    for (double h : {0.100, 0.125}) {
        printf("h=%.3f:\n", h);
        ElasticSolver s;
        s.SetSystem(2, 1, 16, 8, 20.0);
        s.AddVolumeWS ({88.9546, 0.0000}, 1.1489, 0.7508);
        s.AddVolumeWS ({0.0,     2.3480}, 1.3446, 0.6030);
        s.AddSurfaceWS({0.0,   -10.2180}, 1.3943, 0.6872);
        s.AddSpinOrbit({7.1140,  0.0000}, 0.9720, 1.0110);
        s.AddCoulomb(1.3030);
        s.SetGrid(h, 30.0);
        s.SetLmax(2);
        s.Solve();

        const auto& wf = s.GetWavefunction(0, 2);  // L=0, JP=2/2
        printf("  r(fm)   Re(chi)         Im(chi)\n");
        for (int i = 0; i <= (int)(2.0/h); i++) {
            if (i >= (int)wf.size()) break;
            printf("  %.3f  %+.7e  %+.7e\n", i*h, wf[i].real(), wf[i].imag());
        }
        printf("\n");
    }

    return 0;
}
