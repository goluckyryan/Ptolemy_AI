// betcal_faithful_test.cpp
// C++ faithful replica of betcal_faithful_test.f
// Replicates BETCAL + AMPCAL pipeline for 16O(d,p)17O
// using identical random SMAG/SPHASE (same LCG seed=12345)
// and identical JTOCS structure (JP=1,3; LX=2; LDEL=-2..2).
//
// Compare output DCS against Fortran reference:
//   0°:  F^2=659.38  DCS=1446.2 mb/sr
//  10°:  F^2= 71.39  DCS= 156.6 mb/sr
//  30°:  F^2= 52.87  DCS= 115.9 mb/sr
//  45°:  F^2= 43.79  DCS=  96.0 mb/sr
//  90°:  F^2= 26.09  DCS=  57.2 mb/sr
//
// FIXED BUGS (compared to previous version):
//   1. TEMPS indexing: was NMX*(lx-LXMN) + NMLX*(li-LIMN)/2 + mx
//                  → correct: NMX*(lx-LXMN + NMLX*(li-LIMN)/2) + mx
//      (The factor of NMX must multiply the whole inner index, not just lx term.)
//   2. ClebschGordan argument order: was ClebschGordan(li, lx, lo, 0, mx, mx)
//                                  → correct: ClebschGordan(li, 0, lx, mx, lo, mx)
//      (Fortran CLEBSH(2*LI, 2*LX, 0, 2*MX, 2*LO, 2*MX) = <LI,0; LX,MX | LO,MX>)
//      (C++ signature: ClebschGordan(J1, m1, J2, m2, J, m))

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <complex>
#include <vector>

#include "include/math_utils.h"

// ---- Parameters ----
static const int JA   = 2;    // 2*S_deuteron
static const int JB   = 1;    // 2*S_proton (unused in this test)
static const int LMN  = 0;
static const int LMX  = 20;
static const int LSKP = 1;
static const int LXMN = 2;
static const int LXMX = 2;
static const int LDELMX = 2;
static const int LXONLY = 2;
static const int JPBASE = 1;   // |JA-JB| = 1
static const int JPMX   = 3;   // JA+JB = 3
static const int NLO = LMX - LMN + 1;  // 21
static const int NSPL = 10;    // 2 JP values * 5 LDEL values

static const double UK   = 1.23280;
static const double ETA_A = 7.416;
static const double ETA_B = 0.612;
static const double PI = 3.14159265358979320;
static const double DEGREE = PI / 180.0;

// ---- LCG (32-bit signed, Fortran convention) ----
static int32_t lcg_state = 12345;
static float lcg_next() {
    lcg_state = (int32_t)(1664525LL * (int64_t)lcg_state + 1013904223LL);
    return (float)((double)lcg_state * 4.656612873e-10);
}

// ---- JTOCS struct ----
struct JTOCSSlot {
    int ldel, lx, jp;
    bool active;
};

// ---- Coulomb phase: sigma_L = sum_{k=1}^{L} atan(eta/k) ----
// sigma[0]=0, sigma[L] = sum_{k=1}^{L} atan(eta/k)
// Matches Fortran: SIGIN(L+1) = SIGIN(L) + ATAN(ETA/L), SIGIN(1)=0
// So SIGIN(LI+1) in Fortran = sigin[LI] in C++.
std::vector<double> coulomb_phases(int lmax, double eta) {
    std::vector<double> sig(lmax+2, 0.0);
    for (int l = 1; l <= lmax+1; l++)
        sig[l] = sig[l-1] + std::atan2(eta, (double)l);
    return sig;
}

int main() {
    printf("# betcal_faithful_test.cpp — C++ replica (FIXED)\n");
    printf("# JP range: 1,3  LX=2  LDEL=-2..2  LMX=%d\n", LMX);

    // ---- Step 1: Build JTOCS ----
    // Order: JP=1,3 (outer); LDEL=-2,-1,0,1,2 (inner) — same as Fortran test
    std::vector<JTOCSSlot> jtocs;
    for (int jp = JPBASE; jp <= JPMX; jp += 2) {
        for (int ldel = -LDELMX; ldel <= LDELMX; ldel++) {
            jtocs.push_back({ldel, LXONLY, jp, true});
        }
    }
    int nact = (int)jtocs.size();
    printf("# NSPL_ACT=%d\n", nact);

    // Print JTOCS for verification
    for (int k = 0; k < nact; k++) {
        int mx = (jtocs[k].ldel + jtocs[k].lx + 1) / 2;  // C++ integer div (matches Fortran for non-negative)
        // Note: Fortran integer division truncates toward zero (same as C++ for positive)
        // For negative: Fortran (-1)/2 = 0, C++ (-1)/2 = 0. For (-3)/2: F=−1, C++=−1. OK.
        printf("#  k=%2d  ldel=%2d  lx=%d  jp=%d  MX(2ndpass)=%d\n",
               k+1, jtocs[k].ldel, jtocs[k].lx, jtocs[k].jp, mx);
    }

    // ---- Step 2: Random SMAG/SPHASE (same LCG as Fortran) ----
    // Fortran loops: outer LI=LMN..LMX, inner K=1..NACT
    // SMAG(K, ILI) and SPHASE(K, ILI) — stored C-style as [k][ili]
    std::vector<std::vector<float>> SMAG(nact, std::vector<float>(NLO));
    std::vector<std::vector<float>> SPHASE(nact, std::vector<float>(NLO));

    for (int li = LMN; li <= LMX; li++) {
        int ili = li - LMN;
        for (int k = 0; k < nact; k++) {
            SMAG[k][ili]   = std::abs(lcg_next() * 0.5f);
            SPHASE[k][ili] = lcg_next() * (float)PI;
        }
    }

    printf("# SMAG[0][0]=%.6f  SPHASE[0][0]=%.6f\n", SMAG[0][0], SPHASE[0][0]);
    printf("# SMAG[2][4]=%.6f  SPHASE[2][4]=%.6f\n", SMAG[2][4], SPHASE[2][4]);

    // ---- Step 3: Coulomb phases ----
    // sigin[L] = sum_{k=1}^{L} atan(ETA_A/k)  (sigin[0]=0)
    // Fortran: SIGIN(LI+1) = sigin[LI]  and  SIGOT(LO+1) = sigot[LO]
    auto sigin = coulomb_phases(LMX+1, ETA_A);
    auto sigot = coulomb_phases(LMX+1, ETA_B);

    // Verify: Fortran prints SIGIN(1)=sigin[0]=0, SIGIN(3)=sigin[2]
    printf("# SIGIN[0]=%.6f  SIGIN[2]=%.6f\n", sigin[0], sigin[2]);
    printf("# SIGOT[0]=%.6f  SIGOT[2]=%.6f\n", sigot[0], sigot[2]);

    // ---- Step 4: BETCAL — build BETAS[k][lo_idx] ----
    //
    // Ptolemy BETCAL algorithm (reaction mode):
    //   IODD = MOD(LDELMX, 2) = 0  → FACTOR = +0.5/UK
    //   NMX  = LXMX+1 = 3
    //   NMLX = LXMX-LXMN+1 = 1
    //   LBASE = LMN = 0
    //   LOMN = LMN = 0  (LSKP=1, IODD=0)
    //
    //   Outer loop LO = LOMN..LMX:
    //     ILO = LO - LBASE + 1  (1-based in Fortran → 0-based in C++ as lo-LMN)
    //     LIMN = LO - LDELMX
    //     LIMX = LO + LDELMX
    //     MXMX = MIN(LXMX, LO)
    //
    //     Build TEMPS(I) for each (LI, LX, MX):
    //       LI in [LIMN..LIMX step 2], LI >= LBASE
    //       LX in [max(|LI-LO|, LXMN) .. LXMX]
    //       MXZ = MOD(LX+LI-LO, 2)  [parity of LI+LX-LO]
    //       I(0-based) = NMX*(LX-LXMN + NMLX*(LI-LIMN)/2)
    //       TEMPS[I+MX] = FACTOR*(2*LI+1)*CG(LI,0; LX,MX | LO,MX)   for MX=MXZ..LX
    //
    //     Main KOFFS loop (LIPREV = -100000 initially):
    //       For k=0..nact-1 (0-based):
    //         LX = jtocs[k].lx, LDEL = jtocs[k].ldel, LI = LO - LDEL
    //         If LI >= LIPREV: UPDATE MXZ = LX+LDEL (clamp to 0), KOFFZ = k - max(0,MXZ)
    //         LIPREV = LI  (always)
    //         If LI out of [LMN..LMX]: skip accumulation
    //         Load S: AMAG=SMAG[k][LI-LMN], PHASE=SPHASE[k][LI-LMN]+sigin[LI]+sigot[LO]
    //           SMATR = AMAG*sin(PHASE), SMATI = -AMAG*cos(PHASE)
    //           (This is the Ptolemy "divide by i" convention for transfer amplitudes)
    //         TEMPS_idx(0-based) = NMX*(LX-LXMN + NMLX*(LI-LIMN)/2)
    //         For MX=MXZ..LX:
    //           BETAS[KOFFZ+MX][ilo] += TEMPS[TEMPS_idx+MX] * (SMATR + i*SMATI)
    //
    //     Second pass (sqrt-factorial correction):
    //       TEMPS[0] = 1.0
    //       For MX=1..MXMX: TEMPS[MX] = TEMPS[MX-1] / sqrt((LO+MX)*(LO-MX+1))
    //       For k=0..nact-1:
    //         MX = (LDEL + LX + 1) / 2  [integer division, Fortran trunc-toward-zero]
    //         BETAS[k][ilo] *= TEMPS[MX]

    double FACTOR = 0.5 / UK;
    // IODD = MOD(LDELMX, 2) = MOD(2,2) = 0 → FACTOR = +0.5/UK (positive)

    int NMX  = LXMX + 1;          // = 3
    int NMLX = LXMX - LXMN + 1;  // = 1

    // BETAS[k][lo_idx] = complex amplitude for slot k at LO index
    // This matches Fortran BETAS(2, NSPL, NLO) with BETAS(1,k,ilo) = real part
    std::vector<std::vector<std::complex<double>>> BETAS(
        nact, std::vector<std::complex<double>>(NLO, 0.0));

    // TEMPS size: Fortran NTEMPS = (LXMX+1)*(LXMX-LXMN+1)*(LDELMX+1) = 3*1*3 = 9
    // But the actual max index accessed is NMX*(NMX-1 + NMLX*(2*LDELMX/2)) + LXMX
    //   = 3*(0 + 1*2) + 2 = 8  (0-based), so size 9 is sufficient.
    //   But for LO=LMX=20: LIMN=18, LIMX=22 (capped to LMX=20). (LI-LIMN)/2 can be 0 or 1.
    //   Max idx = 3*(0+1*1)+2 = 5. Actually we need NTEMPS=500 as in Fortran test (safe).
    int NTEMPS = 500;

    for (int lo = LMN; lo <= LMX; lo++) {
        int ilo = lo - LMN;  // 0-based index
        int LIMN = lo - LDELMX;
        // int LIMX = lo + LDELMX;  // not used directly
        int MXMX = std::min(LXMX, lo);

        // Build TEMPS: one entry per (LX, LI, MX) triple
        std::vector<double> TEMPS(NTEMPS, 0.0);

        // LI loop: step 2 to match parity (same as Fortran "DO LI = LIMN, LIMX, 2")
        for (int li = LIMN; li <= lo + LDELMX; li += 2) {
            if (li < LMN) continue;  // Fortran: IF (LI .LT. LBASE) GO TO 189
            if (li > LMX) continue;
            int LX1 = std::max(std::abs(li - lo), LXMN);
            for (int lx = LX1; lx <= LXMX; lx++) {
                int mxz = (lx + li - lo) % 2;
                if (mxz < 0) mxz += 2;
                // 0-based TEMPS index: NMX*(lx-LXMN + NMLX*(li-LIMN)/2)
                // Fortran 1-based: I = 1 + NMX*(lx-LXMN + NMLX*(li-LIMN)/2)
                // then accesses TEMPS(I+MX) = TEMPS_0based[NMX*(...) + MX]
                int base = NMX * (lx - LXMN + NMLX * (li - LIMN) / 2);
                for (int mx = mxz; mx <= lx; mx++) {
                    // Fortran: CLEBSH(2*LI, 2*LX, 0, 2*MX, 2*LO, 2*MX)
                    // = CG(j1=LI, m1=0, j2=LX, m2=MX | j3=LO, m3=MX)
                    // C++ signature: ClebschGordan(J1, m1, J2, m2, J, m)
                    double cg = ClebschGordan((double)li, 0.0,
                                              (double)lx, (double)mx,
                                              (double)lo, (double)mx);
                    TEMPS[base + mx] = FACTOR * (2*li + 1) * cg;
                }
            }
        }

        // Main KOFFS loop
        int LIPREV = -100000;
        int MXZ_cur = 0, KOFFZ_cur = 0;

        for (int koffs = 0; koffs < nact; koffs++) {
            if (!jtocs[koffs].active) continue;
            int lx   = jtocs[koffs].lx;
            int ldel = jtocs[koffs].ldel;
            int li   = lo - ldel;

            // Fortran: IF (LI .LT. LIPREV) GO TO 210
            //   (at 210: LIPREV = LI)
            // → update MXZ/KOFFZ when LI >= LIPREV (not decreasing)
            if (li >= LIPREV) {
                MXZ_cur   = lx + ldel;  // MXZ = LX + JTOCS(1,KOFFS) = LX + LDEL
                KOFFZ_cur = koffs - std::max(0, MXZ_cur);
            }
            LIPREV = li;  // always update (label 210)

            if (li < LMN || li > LMX) continue;

            // Load S from SMAG/SPHASE with Coulomb phases
            // Fortran: PHASE = SPHASE(KOFFS,I) + SIGIN(LI+1) + SIGOT(LO+1)
            //   where I = LI - LBASE + 1  (1-based index into SMAG/SPHASE)
            // C++: ili_li = li - LMN (0-based), sigin[li] = SIGIN(LI+1)
            int ili_li = li - LMN;
            float amag = SMAG[koffs][ili_li];
            if (amag == 0.0f) continue;
            float phase_f = SPHASE[koffs][ili_li];
            double phase = (double)phase_f + sigin[li] + sigot[lo];

            // Fortran BETCAL "divide by i" for reaction S-matrix:
            //   SMATR = AMAG * DSIN(PHASE)
            //   SMATI = -AMAG * DCOS(PHASE)
            // These form the complex: SMATR + i*SMATI = AMAG*(sin(φ) - i*cos(φ))
            //   = AMAG * (-i) * (cos(φ) + i*sin(φ)) = AMAG * e^{iφ} / i
            double smatr = (double)amag * std::sin(phase);
            double smati = -(double)amag * std::cos(phase);
            std::complex<double> S_val(smatr, smati);

            // TEMPS index (0-based): NMX*(lx-LXMN + NMLX*(li-LIMN)/2)
            int temps_base = NMX * (lx - LXMN + NMLX * (li - LIMN) / 2);

            // MX loop from MXZ_cur to lx; accumulate into BETAS[KOFFZ+MX]
            for (int mx = std::max(0, MXZ_cur); mx <= lx; mx++) {
                double t = TEMPS[temps_base + mx];
                if (t == 0.0) continue;
                int kidx = KOFFZ_cur + mx;
                if (kidx >= 0 && kidx < nact)
                    BETAS[kidx][ilo] += t * S_val;
            }
        }

        // Second pass: sqrt-factorial correction on BETAS
        // Fortran:
        //   TEMPS(1) = 1.0
        //   DO MX=1,MXMX: TEMPS(MX+1) = TEMPS(MX)/sqrt((DLO+MX)*(DLO-MX+1))
        //   Then for each KOFFS:
        //     MX = (JTOCS(1,KOFFS) + LX + 1) / 2  [Fortran integer division]
        //     BETAS(...,KOFFS,ILO) *= TEMPS(MX+1)
        std::vector<double> sqfac(LXMX + 2, 1.0);
        for (int mx = 1; mx <= MXMX; mx++) {
            double denom = ((double)lo + mx) * ((double)lo - mx + 1.0);
            sqfac[mx] = sqfac[mx-1] / std::sqrt(denom);
        }
        // For mx > MXMX (= MIN(LXMX, LO)): sqfac[mx] is unset.
        // But Fortran only computes up to MXMX too. If MX > MXMX, TEMPS(MX+1) is 0 (cleared at top of LO loop).
        // However TEMPS was zeroed at start, but sqfac here is a separate vector.
        // Handle: if mx > lo (MX > LO), the sqrt has (LO-MX+1) <= 0 → BETAS = 0.
        // Just ensure sqfac[mx]=0 for mx > MXMX:
        for (int mx = MXMX + 1; mx <= LXMX + 1; mx++)
            sqfac[mx] = 0.0;

        for (int koffs = 0; koffs < nact; koffs++) {
            if (!jtocs[koffs].active) continue;
            int lx   = jtocs[koffs].lx;
            int ldel = jtocs[koffs].ldel;
            // Fortran integer division: (LDEL + LX + 1) / 2
            // C++ truncates toward zero, which matches Fortran for non-negative (LDEL+LX+1 >= 0 here).
            // LDEL in [-2..2], LX=2: LDEL+LX+1 in [1..5] → always positive. Safe.
            int mx = (ldel + lx + 1) / 2;
            if (mx < 0 || mx > LXMX + 1) continue;
            BETAS[koffs][ilo] *= sqfac[mx];
        }
    }

    // ---- Print some BETAS for cross-check vs Fortran ----
    printf("#\n");
    printf("# BETAS(Re+Im, k=1..10, ILO=1) [Fortran 1-based ILO=1 = LO=0]:\n");
    for (int k = 0; k < nact; k++) {
        printf("#  k=%2d  BETAS(1,k,1)=%10.5f  BETAS(2,k,1)=%10.5f\n",
               k+1, BETAS[k][0].real(), BETAS[k][0].imag());
    }
    printf("# Expected from Fortran ILO=1 (LO=0):\n");
    printf("#  k=1: -0.00421 + -0.04060i  (JP=1/2, MX=0)\n");
    printf("#  k=6:  0.24275 + -0.17778i  (JP=3/2, MX=0)\n");
    printf("#  k=2..5, k=7..10: all zero\n");

    printf("#\n");
    printf("# BETAS(Re, k=1..10, ILO=1..5):\n");
    for (int ilo = 0; ilo < 5; ilo++) {
        printf("# ILO=%d:", ilo+1);
        for (int k = 0; k < nact; k++)
            printf(" %9.5f", BETAS[k][ilo].real());
        printf("\n");
    }

    // ---- Step 5: AMPCAL + DCS ----
    //
    // Ptolemy AMPCAL:
    //   For each KOFFS slot k:
    //     MX = (LDEL + LX + 1) / 2  [same as second-pass formula]
    //     LPLM = MX*(2*LMX+1-MX)/2  (Fortran 1-based offset → C++ 0-based same formula)
    //     Innermost:
    //       For LO = max(LMN, MX) .. LMX:
    //         F[k] += PLM(LO + LPLM) * BETAS[k][ilo]
    //   where PLM is the Ptolemy-convention Legendre array:
    //     PLM(1..LMX+1) = P_L(x) for L=0..LMX  (MX=0 sector, 1-indexed)
    //     PLM(LPLM+L+1) = P_L^MX(x) for L=MX..LMX
    //   In 0-based C++:
    //     PLM[L]       = P_L(x)   for MX=0, L=0..LMX
    //     PLM[LPLM+L]  = P_L^MX(x) for MX>0, L=MX..LMX
    //
    // Ptolemy PLMSUB convention: uses Condon-Shortley phase convention
    //   (includes (-1)^MX factor, same as std::assoc_legendre).
    //
    // In AMPCAL Fortran:
    //   F(1,K) += PLM(LO+LPLM) * BETAS(1,K,ILO)   [using 1-based PLM array]
    //   (No separate real/imag split here — F is real, BETAS has real+imag parts)
    // Wait: AMPCAL builds F(2,NSPL) where F(1,K)=real part, F(2,K)=imag part.
    // It loops over LO summing PLM*BETAS into both parts.
    //
    // Then DCS = sum_k (F(1,k)^2 + F(2,k)^2) / (uk^2 * (JA+1)) * 10

    printf("#\n");
    printf("# Angle    F^2_sum_raw        DCS_mb_sr\n");
    printf("# Fortran reference:\n");
    printf("#   0.0:  F^2=659.38  DCS=1446.2\n");
    printf("#  10.0:  F^2= 71.39  DCS= 156.6\n");
    printf("#  30.0:  F^2= 52.87  DCS= 115.9\n");
    printf("#  45.0:  F^2= 43.79  DCS=  96.0\n");
    printf("#  90.0:  F^2= 26.09  DCS=  57.2\n");
    printf("#\n");
    printf("# C++ output:\n");

    double angles[] = {0.0, 10.0, 30.0, 45.0, 90.0};

    // PLM array layout (Ptolemy PLMSUB, 0-based C++ indexing):
    //   LPLM(MX) = MX*(2*LMX+1-MX)/2  (same formula as Fortran)
    //   PLM[LPLM(MX) + L] = P_L^MX(cos_theta) for L=MX..LMX
    //   PLM[L] = P_L(cos_theta) for L=0..LMX  [MX=0 sector: LPLM(0)=0]
    // Total size: LPLM(LXMX+1) + LMX + 1 = LXMX*(2*LMX+1-LXMX)/2 + LMX + 1
    int PLM_SIZE = LMX + 1 + LXMX * (2*LMX + 1 - LXMX) / 2 + 10;
    std::vector<double> PLM(PLM_SIZE, 0.0);

    for (int iang = 0; iang < 5; iang++) {
        double angle = angles[iang];
        double an    = DEGREE * angle;
        double costh = std::cos(an);
        double sinth = std::sin(an);

        // Build PLM array (Ptolemy PLMSUB convention)
        // MX=0: standard Legendre P_L(x) for L=0..LMX
        std::fill(PLM.begin(), PLM.end(), 0.0);
        PLM[0] = 1.0;
        if (LMX >= 1) PLM[1] = costh;
        for (int l = 2; l <= LMX; l++)
            PLM[l] = ((2*l - 1) * costh * PLM[l-1] - (l-1) * PLM[l-2]) / l;

        // MX>0: associated Legendre P_L^MX(x) with Condon-Shortley phase
        // Ptolemy PLMSUB: starts from P_MX^MX and recurs upward.
        // P_MX^MX = (-1)^MX * (2*MX-1)!! * sin^MX(theta)
        // The Condon-Shortley convention includes (-1)^MX.
        double fac = 1.0;
        for (int imx = 1; imx <= std::min(LXMX, LMX); imx++) {
            // fac = P_imx^imx = (-1)^imx * (2*imx-1)!! * sin^imx(theta)
            fac *= -(double)(2*imx - 1) * sinth;
            int LPLM = imx * (2*LMX + 1 - imx) / 2;
            // P_imx^imx
            PLM[LPLM + imx] = fac;
            // P_{imx+1}^imx
            if (imx + 1 <= LMX)
                PLM[LPLM + imx + 1] = costh * (2*imx + 1) * fac;
            // P_L^imx for L > imx+1
            for (int l = imx + 2; l <= LMX; l++) {
                PLM[LPLM + l] = ((2*l - 1) * costh * PLM[LPLM + l - 1]
                                 - (l + imx - 1) * PLM[LPLM + l - 2]) / (double)(l - imx);
            }
        }

        // Compute F[k] = sum_LO PLM[LPLM(MX)+LO] * BETAS[k][ilo]
        std::vector<std::complex<double>> F(nact, 0.0);

        for (int koffs = 0; koffs < nact; koffs++) {
            if (!jtocs[koffs].active) continue;
            int lx   = jtocs[koffs].lx;
            int ldel = jtocs[koffs].ldel;
            // MX from the second-pass formula (same as AMPCAL MX formula)
            int mx   = (ldel + lx + 1) / 2;
            if (mx < 0 || mx > lx) continue;

            int LPLM = mx * (2*LMX + 1 - mx) / 2;
            int lomin = std::max(LMN, mx);

            for (int lo = lomin; lo <= LMX; lo++) {
                double plmval = PLM[LPLM + lo];
                int ilo = lo - LMN;
                F[koffs] += plmval * BETAS[koffs][ilo];
            }
        }

        // FSQ = sum_k |F[k]|^2
        double fsq = 0.0;
        for (int k = 0; k < nact; k++)
            fsq += std::norm(F[k]);

        // DCS = FSQ / (uk^2 * (JA+1)) * 10  [mb/sr]
        double dcs = fsq / (UK*UK * (JA+1)) * 10.0;
        printf("%8.3f  %18.8E  %18.8E\n", angle, fsq, dcs);
    }

    return 0;
}
