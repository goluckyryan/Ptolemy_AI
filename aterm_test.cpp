// aterm_test.cpp — Standalone ATERM validation for 16O(d,p)17O g.s.
//
// Tests the ATERM formula in isolation against Ptolemy print=2 output.
//
// Build:
//   g++ -O2 -std=c++17 -Iinclude aterm_test.cpp \
//       src/dwba/math_utils.cpp -o aterm_test && ./aterm_test
//
// Inputs (hardcoded for 16O(d,p)17O g.s., lT=2, jT=2.5, lP=0, jP=0.5):
//   Target residual:    17O g.s.  J=5/2+   SpinResidual = 2.5
//   Target initial:     16O       J=0+     SpinTarget    = 0.0
//   Neutron in target:  lBT=2, jBT=5/2
//   Neutron in proj:    lBP=0, jBP=1/2
//   jx = 1/2 (neutron)
//   SPAMP = 0.97069 (AV18 S-state amplitude)
//   SPAMT = 1.0

#include "math_utils.h"
#include <cmath>
#include <cstdio>

int main() {
    // ---- Quantum numbers for 16O(d,p)17O g.s. ----
    int    lBT = 2;        // neutron orbital L in 17O g.s.
    double jBT = 2.5;      // neutron j in 17O g.s. (d5/2)
    int    lBP = 0;        // neutron orbital L in deuteron (S-state)
    double jBP = 0.5;      // neutron j in deuteron (s1/2)
    double jx  = 0.5;      // neutron spin

    // Nuclear spins
    double SpinTarget   = 0.0;   // 16O: J=0
    double SpinResidual = 2.5;   // 17O g.s.: J=5/2

    int JBIGA = (int)std::round(2.0 * SpinTarget);    // 0
    int JBIGB = (int)std::round(2.0 * SpinResidual);  // 5

    // Spectroscopic amplitudes
    double SPAMP = 0.97069;   // AV18 S-state amplitude
    double SPAMT = 1.0;

    // Allowed Lx: |lBT - lBP| <= Lx <= lBT + lBP → only Lx=2
    int LxMin = std::abs(lBT - lBP);
    int LxMax = lBT + lBP;

    printf("=== ATERM standalone test: 16O(d,p)17O g.s. ===\n");
    printf("  lBT=%d jBT=%.1f  lBP=%d jBP=%.1f  jx=%.1f\n", lBT, jBT, lBP, jBP, jx);
    printf("  JBIGA=%d (2*J_16O)  JBIGB=%d (2*J_17O)\n", JBIGA, JBIGB);
    printf("  SPAMP=%.5f  SPAMT=%.1f\n", SPAMP, SPAMT);
    printf("  Allowed Lx: %d..%d\n\n", LxMin, LxMax);

    // ---- Ptolemy ATERM formula (source.mor ~25631) ----
    // TEMP = sqrt((JBIGB+1)/(JBIGA+1))
    // ATERM = TEMP * sqrt(2*Lx+1) * SPAMP * SPAMT * RACAH(2*lBT, 2*jBT, 2*lBP, 2*jBP, 2*jx, 2*Lx)
    // Sign: ITEST = (JX_doubled - JBP_doubled + 2*(lBP+lBT)) / 2 + 1
    //        if ITEST is odd → flip sign
    //
    // Note: RACAH here is Ptolemy's convention = (-1)^(a/2+b/2+c/2+d/2) * 6j{a/2,b/2,e/2;d/2,c/2,f/2}
    // with args (a=2lBT, b=2jBT, c=2lBP, d=2jBP, e=2jx, f=2Lx)
    // which means 6j{lBT, jBT, jx; jBP, lBP, Lx}

    printf("%-4s  %-10s  %-10s  %-10s  %-12s  %-12s\n",
           "Lx", "TEMP", "sqrt(2Lx+1)", "6j×RACAH_phase", "ATERM(raw)", "ATERM(signed)");
    printf("%-4s  %-10s  %-10s  %-10s  %-12s  %-12s\n",
           "----", "----------", "-----------", "--------------", "------------", "-------------");

    for (int Lx = LxMin; Lx <= LxMax; Lx += 2) {
        double TEMP = std::sqrt((JBIGB + 1.0) / (JBIGA + 1.0));

        // 6-j symbol: {lBT, jBT, jx; jBP, lBP, Lx}
        double sj = SixJ((double)lBT, jBT, jx, jBP, (double)lBP, (double)Lx);

        // RACAH phase: (-1)^((2*lBT + 2*jBT + 2*lBP + 2*jBP)/2)
        int twoj_sum = 2*lBT + (int)(2*jBT + 0.5) + 2*lBP + (int)(2*jBP + 0.5);
        double racah_phase = ((twoj_sum / 2) % 2 == 0) ? 1.0 : -1.0;
        double RACAH_val = racah_phase * sj;

        double ATERM_raw = TEMP * std::sqrt(2.0*Lx + 1.0) * SPAMP * SPAMT * RACAH_val;

        // Sign convention (Ptolemy source.mor 25636-25639)
        int JX_doubled  = (int)std::round(2 * jx);       // = 1
        int JBP_doubled = (int)std::round(2 * jBP);      // = 1
        int ITEST_aterm = JX_doubled - JBP_doubled + 2*(lBP + lBT);
        int ITEST_parity = (ITEST_aterm / 2 + 1);

        double ATERM_signed = ATERM_raw;
        if (ITEST_parity % 2 != 0) ATERM_signed = -ATERM_signed;

        printf("Lx=%-2d  TEMP=%8.5f  sqrt=%8.5f  6j=%10.6f  RACAH=%10.6f  "
               "ATERM_raw=%10.6f  ATERM_signed=%10.6f\n",
               Lx, TEMP, std::sqrt(2.0*Lx+1.0), sj, RACAH_val, ATERM_raw, ATERM_signed);

        // ---- Intermediate breakdown ----
        printf("       ITEST_aterm=%d  ITEST_parity=%d → sign %s\n",
               ITEST_aterm, ITEST_parity,
               (ITEST_parity % 2 != 0) ? "FLIPPED" : "unchanged");
        printf("       6-j args: {%d, %.1f, %.1f; %.1f, %d, %d}\n\n",
               lBT, jBT, jx, jBP, lBP, Lx);
    }

    // ---- Also print SFROMI normalisation factor (for reference) ----
    // SFROMI uses: TEMP_sfromi = FACTOR * ATERM(Lx) / sqrt(2*Li+1)
    // where FACTOR = 2*sqrt(ka*kb/(Ea*Eb))
    // For 16O(d,p)17O at Elab=20 MeV (from Ptolemy print=2 header):
    //   ka = 1.2297 fm^-1 (Incoming.k in C++)
    //   kb = 0.9907 fm^-1 (Outgoing.k)
    //   Ea = 17.234 MeV (Incoming.Ecm)
    //   Eb = 15.479 MeV (Outgoing.Ecm)
    double ka = 1.2297, kb = 0.9907;
    double Ea = 17.234, Eb = 15.479;
    double FACTOR_sfromi = 2.0 * std::sqrt(ka * kb / (Ea * Eb));
    printf("=== SFROMI kinematic FACTOR = %.6f ===\n", FACTOR_sfromi);
    printf("    (= 2*sqrt(ka*kb/(Ea*Eb)) = 2*sqrt(%.4f*%.4f/(%.4f*%.4f)))\n", ka, kb, Ea, Eb);

    // Print full SFROMI norm for Lx=2, each Li
    printf("\n%-4s  %-10s  %-12s  %-12s\n", "Li", "ATERM_sgn", "SFROMI_norm(Lx=2)", "=FACTOR*ATERM/sqrt(2Li+1)");
    // Re-compute ATERM_signed for Lx=2 (should be the only one)
    {
        int Lx = 2;
        double TEMP = std::sqrt((JBIGB + 1.0) / (JBIGA + 1.0));
        double sj = SixJ((double)lBT, jBT, jx, jBP, (double)lBP, (double)Lx);
        int twoj_sum = 2*lBT + (int)(2*jBT+0.5) + 2*lBP + (int)(2*jBP+0.5);
        double racah_phase = ((twoj_sum/2) % 2 == 0) ? 1.0 : -1.0;
        double ATERM_raw = TEMP * std::sqrt(2.0*Lx+1.0) * SPAMP * SPAMT * racah_phase * sj;
        int JX_d = 1, JBP_d = (int)std::round(2*jBP);
        int ITEST_a = JX_d - JBP_d + 2*(lBP + lBT);
        int ITEST_p = (ITEST_a / 2 + 1);
        double ATERM_signed = (ITEST_p % 2 != 0) ? -ATERM_raw : ATERM_raw;

        for (int Li = 0; Li <= 5; Li++) {
            double sfromi_norm = FACTOR_sfromi * std::abs(ATERM_signed) / std::sqrt(2.0*Li + 1.0);
            printf("Li=%-2d  ATERM=%9.6f  SFROMI_norm=%12.6f\n",
                   Li, ATERM_signed, sfromi_norm);
        }
    }

    printf("\n=== Expected from Ptolemy (read from print=2 output) ===\n");
    printf("  Li=0 JPI=2/2 Lo=2 JPO=3/2 Lx=2: Ptolemy S = (0.1212, 0.04796)\n");
    printf("  Li=0 JPI=2/2 Lo=2 JPO=5/2 Lx=2: Ptolemy S = (0.1157, 0.1201)\n");
    printf("  → Compare ATERM_signed value above to Ptolemy's ATERM printout (if any)\n");
    printf("  → Ptolemy doesn't print ATERM directly, but you can back-calculate:\n");
    printf("     ATERM = S_ptolemy / (radial_integral * SFROMI_norm / ATERM)\n");

    return 0;
}
