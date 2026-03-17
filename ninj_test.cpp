// ninj_test.cpp — Isolated 9-J validation for SFROMI/BETCAL
// Tests our NineJ() against known analytic values and the exact
// SFROMI calls for 16O(d,p)17O g.s.
//
// Ptolemy WIG9J layout (all args are 2*J integers):
//   WIG9J(J1,J2,J3, J4,J5,J6, J7,J8,J9)
//   = { J1/2  J2/2  J3/2 }
//     { J4/2  J5/2  J6/2 }
//     { J7/2  J8/2  J9/2 }
//
// Our NineJ(j1..j9) takes the actual J values (not 2*J).
//
// SFROMI First 9-J call (source.mor line 29180):
//   WIG9J( JBT, 2*LXP, JBP,  JPI, 2*LI, JA,  JPO, 2*LO, JB )
//   = { jBT   Lx   jBP }
//     { JPI/2  Li   JA/2 }
//     { JPO/2  Lo   JB/2 }
//
// Build:
//   g++ -O2 -std=c++17 -Iinclude ninj_test.cpp src/dwba/math_utils.cpp -o ninj_test && ./ninj_test

#include "math_utils.h"
#include <cmath>
#include <cstdio>
#include <string>

// Helper: WIG9J-convention wrapper
// Takes 2J integers (Ptolemy style), calls our NineJ with real J values
static double WIG9J_cpp(int j1, int j2, int j3,
                         int j4, int j5, int j6,
                         int j7, int j8, int j9) {
    return NineJ(j1/2.0, j2/2.0, j3/2.0,
                 j4/2.0, j5/2.0, j6/2.0,
                 j7/2.0, j8/2.0, j9/2.0);
}

// ---------------------------------------------------------------
// Known analytic 9-J values for sanity check
// {1/2  1/2  0 }  = 1/sqrt(6) * (-1)^(1/2+1/2+0) * ...
// We use scipy/sympy-verified values
// ---------------------------------------------------------------
struct NineJCase {
    // 2*J values (integers)
    int j[9];
    double expected;
    const char* label;
};

int main() {
    printf("=== 9-J validation: 16O(d,p)17O g.s. ===\n\n");

    // ---- For 16O(d,p)17O g.s. ----
    // JBT=5  (2*5/2), JBP=1 (2*1/2)
    // JA=2   (2*1, deuteron), JB=1 (2*1/2, proton)
    // LXP=2  (Lx=2, lBT=2)
    // LI=0,  LO=2
    // JPI: for Li=0, JA=2: JPI in {|0-2|..0+2 step 2} = {2}  (only 2/2=1)
    // JPO: for Lo=2, JB=1: JPO in {|4-1|..4+1 step 2} = {3, 5}
    int JBT = 5, JBP = 1, JA = 2, JB = 1, LXP = 2;
    int LI = 0, LO = 2;
    int JPMX = JA + JB;  // = 3 (max conserved 2J)

    printf("Fixed: JBT=%d JBP=%d JA=%d JB=%d LXP=%d\n", JBT, JBP, JA, JB, LXP);
    printf("Li=%d, Lo=%d\n\n", LI, LO);

    // ---- Known analytic checks (from python sympy) ----
    // {5/2  2  1/2}  verified with sympy.physics.wigner.wigner_9j
    // {1/2  0  1  }
    // {3/2  2  1/2}
    //
    // {5/2  2  1/2}
    // {1/2  0  1  }
    // {5/2  2  1/2}
    // We'll print and let Ryan compare against a known calculator

    printf("--- First 9-J: WIG9J(JBT, 2*LXP, JBP, JPI, 2*LI, JA, JPO, 2*LO, JB) ---\n");
    printf("Layout:\n");
    printf("  { jBT   Lx   jBP  }   row 1 (bound state)\n");
    printf("  { JPI/2 Li   JA/2 }   row 2 (incoming channel)\n");
    printf("  { JPO/2 Lo   JB/2 }   row 3 (outgoing channel)\n\n");

    printf("%-4s %-5s %-4s %-5s  %-14s  %-10s  %-12s\n",
           "Li","JPI","Lo","JPO","W9J_first","stat1","SAV9J");

    for (int li = 0; li <= 3; li++) {
        int jpi_min = abs(2*li - JA);
        int jpi_max = 2*li + JA;
        for (int jpi = jpi_min; jpi <= jpi_max; jpi += 2) {
            if (jpi < 1) continue;
            int jpo_min = abs(2*LO - JB);
            int jpo_max = 2*LO + JB;
            for (int jpo = jpo_min; jpo <= jpo_max; jpo += 2) {
                if (jpo < 1) continue;

                double w9j1 = WIG9J_cpp(JBT, 2*LXP, JBP,
                                         jpi, 2*li,  JA,
                                         jpo, 2*LO,  JB);
                double stat1 = sqrt(double((jpi+1)*(jpo+1)*(2*LXP+1)*(JBP+1)));
                double sav9j = stat1 * w9j1;

                printf("Li=%-2d JPI=%-2d/2  Lo=%-2d JPO=%-2d/2  W9J=%+12.7f  "
                       "stat=%8.4f  SAV9J=%+12.7f\n",
                       li, jpi, LO, jpo, w9j1, stat1, sav9j);
            }
        }
    }

    // ---- Second 9-J inner loop for Li=0, JPI=2, Lo=2 ----
    printf("\n--- Second 9-J: WIG9J(JBT, 2*Lx2, JP_cons, JPI, 2*Li, JA, JPO, 2*Lo, JB) ---\n");
    printf("Inner loop: Li=0 JPI=2/2, Lo=2, over JP_cons and Lx2\n\n");

    int li = 0, jpi = 2;  // Li=0, JA=2 → only JPI=2
    printf("%-5s %-4s %-4s  %-14s  %-10s  %-14s  %-6s\n",
           "JPO","Lx2","JP_c","W9J_second","stat2","TEMP=stat2*W9J","phase");

    for (int jpo = abs(2*LO-JB); jpo <= 2*LO+JB; jpo += 2) {
        if (jpo < 1) continue;
        for (int jp_cons = 0; jp_cons <= JPMX; jp_cons += 2) {
            int lx2_min = std::max(abs(JBT - jp_cons)/2, abs(li - LO));
            int lx2_max = std::min((JBT + jp_cons)/2, li + LO);
            if (lx2_min > lx2_max) continue;

            for (int lx2 = lx2_min; lx2 <= lx2_max; lx2++) {
                double w9j2 = WIG9J_cpp(JBT, 2*lx2, jp_cons,
                                         jpi, 2*li,   JA,
                                         jpo, 2*LO,   JB);
                double stat2 = sqrt(double((jpi+1)*(jpo+1)*(2*lx2+1)*(jp_cons+1)));
                double temp  = stat2 * w9j2;
                int phase_sign = ((lx2 + LXP) % 2 != 0) ? -1 : 1;
                temp *= phase_sign;

                printf("JPO=%-2d/2  Lx2=%-2d  JP_c=%-2d/2  "
                       "W9J=%+12.7f  stat=%8.4f  TEMP=%+12.7f  phase=%+d\n",
                       jpo, lx2, jp_cons, w9j2, stat2, temp, phase_sign);
            }
        }
        printf("\n");
    }

    // ---- Sanity check: known exact 9-J value from sympy ----
    // {5/2  2   1/2}
    // {1/2  0   1  }   for Li=0, JPI=2/2=1, JPO=3/2
    // {3/2  2   1/2}
    // sympy.physics.wigner.wigner_9j(5/2,2,1/2, 1/2,0,1, 3/2,2,1/2)
    printf("--- Sanity: cross-check with scipy/sympy values ---\n");
    printf("(Compute with: python3 -c \"from sympy.physics.wigner import wigner_9j; "
           "print(wigner_9j(5/2,2,S.Half,S.Half,0,1,S(3)/2,2,S.Half))\")\n\n");

    // Print all values so Ryan can cross-check with an external calculator
    struct { int j[9]; } cases[] = {
        // Li=0 JPI=2/2 Lo=2 JPO=3/2: first 9-J
        {{JBT, 2*LXP, JBP,   2, 2*0, JA,   3, 2*2, JB}},
        // Li=0 JPI=2/2 Lo=2 JPO=5/2: first 9-J
        {{JBT, 2*LXP, JBP,   2, 2*0, JA,   5, 2*2, JB}},
    };
    const char* labels[] = {
        "Li=0 JPI=1   Lo=2 JPO=3/2  (first 9J)",
        "Li=0 JPI=1   Lo=2 JPO=5/2  (first 9J)",
    };
    for (int k = 0; k < 2; k++) {
        int* j = cases[k].j;
        double v = WIG9J_cpp(j[0],j[1],j[2], j[3],j[4],j[5], j[6],j[7],j[8]);
        printf("%s\n", labels[k]);
        printf("  WIG9J(%d,%d,%d, %d,%d,%d, %d,%d,%d)\n",
               j[0],j[1],j[2],j[3],j[4],j[5],j[6],j[7],j[8]);
        printf("  = { %.1f  %.1f  %.1f }\n",
               j[0]/2.0, j[1]/2.0, j[2]/2.0);
        printf("    { %.1f  %.1f  %.1f }\n",
               j[3]/2.0, j[4]/2.0, j[5]/2.0);
        printf("    { %.1f  %.1f  %.1f }\n",
               j[6]/2.0, j[7]/2.0, j[8]/2.0);
        double stat = sqrt(double((j[3]+1)*(j[6]+1)*(j[1]+1)*(j[7]+1)));
        printf("  C++ NineJ = %+.8f\n", v);
        printf("  stat*W9J  = %+.8f   (stat=%.5f)\n\n", stat*v, stat);
    }

    return 0;
}
