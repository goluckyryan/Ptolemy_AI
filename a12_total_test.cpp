// a12_total_test.cpp — compute total A12(phi_T, phi_ab) for Li=1, Lo=1 or 3, Lx=2
// and compare directly with what we expect from Ptolemy's terms
#include "dwba.h"
#include "math_utils.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <tuple>

int main() {
    // 16O(d,p)17O gs: lT=2, lP=0, Lx=2
    int lT = 2, lP = 0, Lx = 2;

    // Test A12 for Li=1 with different Lo values
    for (int Li = 1; Li <= 1; Li++) {
        for (int Lo : {1, 3}) {
            DWBA dwba;
            auto terms = dwba.ComputeA12Terms(Li, Lo, Lx, lT, lP);

            printf("\n=== A12 terms: Li=%d Lo=%d Lx=%d (lT=%d lP=%d) ===\n", Li, Lo, Lx, lT, lP);
            printf("  Number of terms: %zu\n", terms.size());
            for (auto& [MT, MU, coeff] : terms) {
                printf("  MT=%2d MU=%2d coeff=%+.6f\n", MT, MU, coeff);
            }

            // Evaluate A12 at several (phi_T, phi_ab) angles
            printf("  A12 values at key angles:\n");
            for (double phi_T : {0.0, M_PI/4, M_PI/2}) {
                for (double phi_ab : {0.0, M_PI/4, M_PI/2}) {
                    double val = dwba.EvalA12(terms, phi_T, phi_ab);
                    printf("    phi_T=%.4f phi_ab=%.4f -> A12=%.8f\n", phi_T, phi_ab, val);
                }
            }
        }
    }

    // Now: what does Ptolemy compute for the ANSWER terms?
    // Ptolemy printed for Li=1, MT=-2, MU=1:
    //   LO=1: TEMP=0.0000  stored=0.0000
    //   LO=3: TEMP=-0.20963 stored=-0.41926
    // From our a12_debug, for MT=-2, MU=1 (only):
    //   LO=3: stored=-0.41926 ✓
    //   LO=1: stored=0.0000 ✓
    // For MT=0, MU=1:
    //   LO=1: stored=+0.13693
    //   LO=3: stored=-0.16771
    // For MT=2, MU=1:
    //   LO=1: stored=+0.41079
    //   LO=3: stored=-0.08385
    //
    // A12VL[LO=1] total = 0 + 0.13693 + 0.41079 = 0.54772
    // A12VL[LO=3] total = -0.41926 + (-0.16771) + (-0.08385) = -0.67082
    //
    // Our ComputeA12Terms for Li=1, Lo=1: sum of all MT,MU contributions
    // The stored values go into A12VL array indexed by LO.
    // EvalA12 for Lo=1: cos(MT*phi_T - MU*phi_ab) * coeff
    // At phi_T=0, phi_ab=0: sum of all coeffs

    printf("\n--- Expected from summing Ptolemy ANSWER terms ---\n");
    printf("For Li=1, Lo=1, Lx=2: sum at phi_T=0,phi_ab=0 = 0.0 + 0.13693 + 0.41079 = %.5f\n",
           0.0 + 0.13693 + 0.41079);
    printf("For Li=1, Lo=3, Lx=2: sum at phi_T=0,phi_ab=0 = -0.41926 + (-0.16771) + (-0.08385) = %.5f\n",
           -0.41926 + (-0.16771) + (-0.08385));

    return 0;
}
