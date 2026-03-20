// a12_debug.cpp — compare our A12 terms vs Ptolemy annotated output
// For Li=1, Lo=1, Lx=2 (16O d,p ground state case)
#include "dwba.h"
#include "math_utils.h"
#include <cstdio>
#include <cmath>

// Re-implement A12 with full debug output matching Ptolemy format
// OUTTMP: MP, MT, LX, xlam_T, xlam_P, outtmp
// OUTTER: MP, MT, MU, LI, LO, xlam_Li, sqrt(2LO+1), outter
// ANSWER: MP, MT, MU, LX, LO, LOMNMN, LXMIN, XLOTMP, OUTTMP_LX, TTT, TEMP

static double xlam_correct(int L, int am) {
    am = std::abs(am);
    if (am > L) return 0.0;
    if ((L + am) % 2 != 0) return 0.0;
    if (L == 0) return 1.0;

    static double cache[20][10];
    double outter = 1.0;
    int m_cur = 0;
    cache[0][0] = 1.0;
    for (int ll = 1; ll <= L; ll++) {
        m_cur = 1 - m_cur;
        outter *= std::sqrt((double)(ll + m_cur - 1) / (double)(ll + m_cur));
        if (m_cur == 1) outter = -outter;
        cache[ll][m_cur / 2] = outter;
        for (int mm = m_cur + 2; mm <= ll; mm += 2) {
            int idx = mm / 2;
            cache[ll][idx] = -cache[ll][idx-1] *
                std::sqrt((double)(ll-mm+2)*(double)(ll+mm-1) /
                         ((double)(ll+mm)*(double)(ll-mm+1)));
        }
    }
    return cache[L][am / 2];
}

int main() {
    // 16O(d,p)17O: lT=2 (neutron orbital), lP=0 (proton s-wave)
    // We run Li=1 case
    int lT = 2, lP = 0;
    int LXMIN = 2, LXMAX = 2;  // for 16O(d,p) gs, lT=2 so Lx=2

    for (int Li = 1; Li <= 1; Li++) {
        printf("\n=== Li=%d ===\n", Li);
        int LLI = 2*Li;
        double XN = 0.5 * std::sqrt((double)(2*Li+1) * (double)(2*lT+1) * (double)(2*lP+1));

        // Compute LOMNMN
        int LOMNMN = std::abs(Li - LXMIN);
        LOMNMN += (lT + lP + Li + LOMNMN) % 2;
        int LOMXMX = Li + LXMAX;
        LOMXMX -= (lT + lP + Li + LOMXMX) % 2;
        printf("LOMNMN=%d LOMXMX=%d XN=%.6f\n", LOMNMN, LOMXMX, XN);

        // HALFSW: false for our case (lT=2 even, lP=0 even)
        bool HALFSW = false;
        int MTMIN = -lT, MPMIN = -lP;
        int MUSTRT = Li % 2;  // =1 for Li=1
        if (HALFSW) MUSTRT = -Li;

        // MT loop
        for (int MP = MPMIN; MP <= lP; MP += 2) {
            for (int MT = MTMIN; MT <= lT; MT += 2) {
                int MX = MT + MP;
                if (std::abs(MX) > LXMAX) continue;

                double xlam_T = xlam_correct(lT, std::abs(MT));
                double xlam_P = xlam_correct(lP, std::abs(MP));

                int IP3 = 0, IP4 = 0;
                if (MT < 0) IP3 = lT;
                if (MP < 0) IP4 = lP;

                // OUTTMP for each LX
                for (int LX = LXMIN; LX <= LXMAX; LX += 1) {
                    double outtmp = xlam_T * xlam_P * XN *
                        ThreeJ((double)lT, (double)MT, (double)lP, (double)MP, (double)LX, (double)(-MX));
                    if (std::abs(outtmp) < 1e-15) continue;

                    printf("OUTTMP  MP=%2d MT=%2d LX=%d  xlam_T=%.5f xlam_P=%.5f outtmp=%.5f\n",
                           MP, MT, LX, xlam_T, xlam_P, outtmp);

                    // MU loop (HALFSW=false: MU from MUSTRT to LI step 2)
                    for (int MU = MUSTRT; MU <= Li; MU += 2) {
                        int MX_minus_MU = MX - MU;
                        if (std::abs(MX_minus_MU) > LOMXMX) continue;

                        double xlam_Li = xlam_correct(Li, std::abs(MU));
                        if (std::abs(xlam_Li) < 1e-15) continue;

                        // OUTTER loop over LO
                        for (int LO = LOMNMN; LO <= LOMXMX; LO += 2) {
                            if (LO < std::abs(MX_minus_MU)) continue;
                            double xlam_Lo = xlam_correct(LO, std::abs(MX_minus_MU));
                            double outter = xlam_Li * xlam_Lo * std::sqrt(2.0*LO+1.0);

                            // ITES sign
                            int IP1 = 0, IP2 = 0;
                            if (MU < 0) IP1 = Li;
                            if (MX_minus_MU < 0) IP2 = LOMNMN;
                            int ITES = IP4 + IP3 + IP2 + IP1;
                            if (ITES % 2 != 0) outter = -outter;

                            printf("  OUTTER MP=%2d MT=%2d MU=%2d Li=%d LO=%d  xlam_Li=%.5f sqrt(2LO+1)=%.5f outter=%.5f\n",
                                   MP, MT, MU, Li, LO, xlam_Li, std::sqrt(2.0*LO+1.0), outter);

                            // TTT inner 3-J
                            double ttt = ThreeJ((double)Li, (double)MU,
                                               (double)LO, (double)MX_minus_MU,
                                               (double)LX, (double)(-MX));

                            double TEMP = outter * outtmp * ttt;
                            double doubling = (HALFSW || MU != 0) ? 2.0 : 1.0;

                            printf("    ANSWER MP=%2d MT=%2d MU=%d LX=%d LO=%d LOMNMN=%d: "
                                   "outter=%.5f outtmp=%.5f TTT=%.5f TEMP=%.5f -> stored=%.5f\n",
                                   MP, MT, MU, LX, LO, LOMNMN,
                                   outter, outtmp, ttt, TEMP, TEMP*doubling);
                        }
                    }
                }
            }
        }
    }
    return 0;
}
