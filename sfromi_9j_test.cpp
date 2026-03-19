// sfromi_9j_test.cpp
// C++ replica of sfromi_9j_test.f — SFROMI 9-J branch
// for 16O(d,p)17O g.s., JA=2, JB=1, JBP=1, JBT=5, LXP=2
// Uses same LCG seed and identical 9-J loop structure.
// Prints accumulated S-matrix slots to stdout.
// Compare against Fortran: /tmp/sfromi_9j_fort.txt

#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>
#include <cstdint>

// Include NineJ from the Cpp_AI library
#include "include/math_utils.h"

// Constants (16O(d,p)17O)
static const int JA   = 2;   // 2*S_deuteron
static const int JB   = 1;   // 2*S_proton
static const int JBP  = 1;   // 2*j_neutron_in_deuteron (0s1/2)
static const int JBT  = 5;   // 2*j_neutron_in_17O (0d5/2)
static const int LXP  = 2;   // transferred L
static const int JPBASE = 1; // |JA-JB| = |2-1| = 1
static const int LMAX = 8;
static const int LDMAX = 16;

static const double FACTOR = 0.11100;
static const double ATERM  = -1.681284;

// LCG — Fortran 32-bit signed convention
static int32_t lcg_state = 98765;
static double lcg_next() {
    lcg_state = (int32_t)(1664525LL * (int64_t)lcg_state + 1013904223LL);
    return (double)lcg_state * 4.656612873e-10;
}

// S-matrix indexed by (JP/2, LX, IDEL=LO-LI+LDMAX)
struct SKey {
    int jp2, lx, idel;
    bool operator<(const SKey &o) const {
        if (jp2 != o.jp2) return jp2 < o.jp2;
        if (lx  != o.lx ) return lx  < o.lx;
        return idel < o.idel;
    }
};
std::map<SKey, double> SMATR, SMATI;

int main() {
    printf("# SFROMI 9-J test (C++): 16O(d,p)17O g.s.\n");
    printf("# JA=%d, JB=%d, JBP=%d, JBT=%d, LXP=%d\n", JA,JB,JBP,JBT,LXP);
    printf("#\n");

    int JPMX = JBT;  // conservative upper bound for JP loop

    for (int LI = 0; LI <= LMAX; LI++) {
        for (int LO = 0; LO <= LMAX; LO++) {

            // Random XI
            double XI_R = lcg_next() * 0.2;
            double XI_I = lcg_next() * 0.2;

            // SFROMI core
            double TEMP = FACTOR * ATERM / std::sqrt(double(2*LI+1));
            int ITEST = LI + LO + 2*LXP + 1;
            if (ITEST % 4 >= 2) TEMP = -TEMP;
            double TEMPR = TEMP * XI_R;
            double TEMPI = TEMP * XI_I;
            double SR = TEMPR, SI = TEMPI;
            if (ITEST % 2 != 0) { SR = -TEMPI; SI = TEMPR; }

            // 9-J branch: loop over JPI, JPO (SOSWS active)
            int JPIMN = std::abs(2*LI - JA);
            int JPIMX = 2*LI + JA;
            int JPOMN = std::abs(2*LO - JB);
            int JPOMX = 2*LO + JB;

            for (int JPI = JPIMN; JPI <= JPIMX; JPI += 2) {
                for (int JPO = JPOMN; JPO <= JPOMX; JPO += 2) {

                    // First 9-J
                    double w9j1 = NineJ(JBT/2.0, (double)LXP, JBP/2.0,
                                        JPI/2.0, (double)LI,  JA/2.0,
                                        JPO/2.0, (double)LO,  JB/2.0);
                    if (std::isnan(w9j1) || w9j1 == 0.0) continue;

                    double SAV9J = std::sqrt(double((JPI+1)*(JPO+1)*(2*LXP+1)*(JBP+1))) * w9j1;
                    if (SAV9J == 0.0) continue;

                    double TEMPR2 = SAV9J * SR;
                    double TEMPI2 = SAV9J * SI;

                    // Inner loop: JP and LX
                    for (int JP = JPBASE; JP <= JPMX; JP += 2) {
                        int LXMN2 = std::max(std::abs(JBT-JP)/2, std::abs(LO-LI));
                        int LXMX2 = std::min((JBT+JP)/2, LO+LI);
                        if (LXMN2 > LXMX2) continue;

                        for (int LX = LXMN2; LX <= LXMX2; LX++) {
                            double TEMP2;
                            if (LX == LXP && JP == JBP) {
                                TEMP2 = SAV9J;
                            } else {
                                double w9j2 = NineJ(JBT/2.0, (double)LX,  JP/2.0,
                                                    JPI/2.0, (double)LI,  JA/2.0,
                                                    JPO/2.0, (double)LO,  JB/2.0);
                                if (std::isnan(w9j2) || w9j2 == 0.0) continue;
                                TEMP2 = std::sqrt(double((JPI+1)*(JPO+1)*(2*LX+1)*(JP+1))) * w9j2;
                            }
                            if (TEMP2 == 0.0) continue;

                            if ((LX + LXP) % 2 != 0) TEMP2 = -TEMP2;

                            int IDEL = LO - LI + LDMAX;
                            int JP2  = JP / 2;
                            SKey k{JP2, LX, IDEL};
                            SMATR[k] += TEMP2 * TEMPR2;
                            SMATI[k] += TEMP2 * TEMPI2;
                        }
                    }
                }
            }
        }
    }

    // Print all non-zero slots
    printf("# JP  LX  dLOLI    SMATR           SMATI\n");
    for (auto &[k, v] : SMATR) {
        double sr = v, si = SMATI[k];
        if (std::abs(sr) + std::abs(si) > 1e-20) {
            printf("%4d%4d%4d  %18.10E  %18.10E\n",
                   2*k.jp2, k.lx, k.idel - LDMAX, sr, si);
        }
    }

    return 0;
}
