// wig9j.h — Faithful C++ port of Ptolemy's WIG9J (S. Pieper, 1974)
// Computes the Wigner 9-J coefficient using the fused Racah sum method.
// Input: 9 integers that are TWICE the corresponding j values.
// Reference: J.G. Wills, Computer Physics Communications 2, 38 (1971)

#ifndef WIG9J_H
#define WIG9J_H

#include <cmath>
#include <algorithm>

// Log-factorial table (initialized on first use)
static double _wig9j_logfac[2001];
static bool   _wig9j_logfac_init = false;
static const int _wig9j_maxlf = 2000;

static void _wig9j_init_logfac() {
    if (_wig9j_logfac_init) return;
    _wig9j_logfac[0] = 0.0;
    for (int i = 1; i <= _wig9j_maxlf; i++)
        _wig9j_logfac[i] = _wig9j_logfac[i-1] + std::log((double)i);
    _wig9j_logfac_init = true;
}

// Log(n!) — uses table for small n, lgamma for large n
static double _wig9j_lf(int n) {
    if (n < 0) return -1e300;  // invalid
    if (n <= _wig9j_maxlf) return _wig9j_logfac[n];
    return std::lgamma((double)(n + 1));
}

// Wigner 9-J symbol.  Arguments are 2*j (integers).
//
//   { j1  j2  j3 }
//   { j4  j5  j6 }     =  WIG9J(2*j1, 2*j2, ..., 2*j9)
//   { j7  j8  j9 }
//
static double WIG9J_int(int J1, int J2, int J3, int J4, int J5, int J6,
                        int J7, int J8, int J9) {
    _wig9j_init_logfac();

    // Rearrange J's for the 3 internal 6-J symbols:
    //   { J1S(I)    J2S(I)    J3S(I)  }
    //   { J1S(I+1)  J2S(I+2)    X     }    for I = 1, 2, 3 (cyclic)
    int J1S[6], J2S[6], J3S[6];
    J1S[1] = J1;  J1S[2] = J6;  J1S[3] = J8;
    J2S[1] = J2;  J2S[2] = J4;  J2S[3] = J9;
    J3S[1] = J3;  J3S[2] = J5;  J3S[3] = J7;
    // Cyclic extension
    J1S[4] = J1S[1]; J1S[5] = J1S[2];
    J2S[4] = J2S[1]; J2S[5] = J2S[2];
    J3S[4] = J3S[1]; J3S[5] = J3S[2];

    // Determine the range of 2*X (the summation variable)
    int xstart2 = std::abs(J2S[1] - J1S[2]);  // 2*xstart
    int xend2   = J2S[1] + J1S[2];
    int xsum2   = xend2;
    for (int i = 2; i <= 3; i++) {
        xstart2 = std::max(xstart2, std::abs(J2S[i] - J1S[i+1]));
        int xpiece = J2S[i] + J1S[i+1];
        xend2 = std::min(xend2, xpiece);
        xsum2 += xpiece;
    }
    if (xend2 < xstart2) return 0.0;  // triangle violation

    xsum2 -= xend2;

    // Extract half-integer part and convert 2*X to integer part
    int xhalf = (xstart2 & 1);  // 1 if x is half-integer, 0 if integer
    int xstart = xstart2 >> 1;  // integer part of x
    int xend   = xend2 >> 1;
    int xsum   = (xsum2 >> 1) + 1;
    double dxhalf = (xhalf == 1) ? 1.0 : 0.0;

    // Setup AS, BS arrays for each of the 3 internal 6-J's
    // All indices are in units of integer parts (original 2j values divided by 2 where appropriate)
    int AS[5][4], BS[3][4];  // [j][i], 1-indexed
    double DAS[5][4], DBS[3][4];
    int MN1S[4], MN2S[4], MX2S[4], RS_arr[4], IDS_arr[3][4];
    double DMN1S[4], DMN2S[4], DMX2S[4], DRS[4];

    for (int i = 1; i <= 3; i++) {
        // A1 = (J1S(I) + J2S(I) - J3S(I)) / 2
        AS[1][i] = (J1S[i] + J2S[i] - J3S[i]) >> 1;
        // A2 = (J1S(I+1) + J2S(I+2) - J3S(I)) / 2
        AS[2][i] = (J1S[i+1] + J2S[i+2] - J3S[i]) >> 1;
        // A3 = (J1S(I) + J2S(I+2)) / 2   (X component left out)
        AS[3][i] = (J1S[i] + J2S[i+2]) >> 1;
        // A4 = (J2S(I) + J1S(I+1)) / 2
        AS[4][i] = (J2S[i] + J1S[i+1]) >> 1;

        // Triangle check: AS[1] and AS[2] must be non-negative
        if (AS[1][i] < 0 || AS[2][i] < 0) return 0.0;

        // B1 = J2S(I+2) - A2 - A3  (X component left out)
        BS[1][i] = J2S[i+2] - AS[2][i] - AS[3][i];
        // B2 = J1S(I) - A1 - A3
        BS[2][i] = J1S[i] - AS[1][i] - AS[3][i];
    }

    // Setup derived quantities for each 6-J
    for (int i = 1; i <= 3; i++) {
        MN1S[i] = std::min(AS[1][i], AS[2][i]);
        MN2S[i] = std::min(AS[3][i], AS[4][i]);
        MX2S[i] = -std::min(BS[1][i], BS[2][i]);
        RS_arr[i] = 1 + xhalf + AS[3][i] + AS[4][i];
        for (int j = 1; j <= 4; j++) DAS[j][i] = (double)AS[j][i];
        DBS[1][i] = (double)BS[1][i];
        DBS[2][i] = (double)BS[2][i];
        DRS[i] = 1.0 + dxhalf + DAS[3][i] + DAS[4][i];
        DMN1S[i] = (double)MN1S[i];
        DMN2S[i] = (double)MN2S[i];
        DMX2S[i] = (double)MX2S[i];
        IDS_arr[1][i] = AS[2][i] + BS[2][i];
        IDS_arr[2][i] = AS[1][i] + BS[1][i];
    }

    // Compute the log of outside factors (factorials that completely factor out)
    double sumout = 0.0;
    if (xsum <= _wig9j_maxlf) {
        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 2; j++) {
                int idx_j1 = (j == 1) ? i : i + 1;  // J1S index: I for j=1, I+1 for j=2
                int idx_j2 = (j == 1) ? i : i + 2;  // J2S index: I for j=1, I+2 for j=2
                // Fortran: J1S(I-1+J) and J2S(I-2+2*J)
                // j=1: J1S(I), J2S(I)     → idx_j1=I, idx_j2=I
                // j=2: J1S(I+1), J2S(I+2) → idx_j1=I+1, idx_j2=I+2  ← matches!
                sumout += _wig9j_lf(AS[j][i])
                        - _wig9j_lf(1 + J3S[i] + AS[j][i])
                        + _wig9j_lf(-AS[j][i] + J1S[idx_j1])
                        + _wig9j_lf(-AS[j][i] + J2S[idx_j2]);
            }
        }
    } else {
        // Use lgamma for large arguments
        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 2; j++) {
                int idx_j1 = (j == 1) ? i : i + 1;
                int idx_j2 = (j == 1) ? i : i + 2;
                sumout += std::lgamma(1.0 + DAS[j][i])
                        - std::lgamma(2.0 + J3S[i] + AS[j][i])
                        + std::lgamma(1.0 - AS[j][i] + J1S[idx_j1])
                        + std::lgamma(1.0 - AS[j][i] + J2S[idx_j2]);
            }
        }
    }
    sumout *= 0.5;

    // Main 9-J loop over X (integer part)
    double result = 0.0;
    double DX = (double)(xstart - 1);

    for (int X = xstart; X <= xend; X++) {
        DX += 1.0;

        double sumlog = sumout;
        int isign = 0;
        double SUMS[4];  // 1-indexed: SUMS[1..3]

        // Process each of the 3 internal 6-J's
        bool triangle_fail = false;
        for (int i = 1; i <= 3; i++) {
            int MX2 = MX2S[i] - X;
            int MN2 = MN2S[i] - X;
            int NMAX = MN1S[i];
            double DNMAX = DMN1S[i];
            double DNMIN = 0.0;
            int NMIN = 0;

            // Determine range of Racah sum
            if (MN2 < NMAX) {
                NMAX = MN2;
                DNMAX = DMN2S[i] - DX;
            }
            if (MX2 > 0) {
                NMIN = MX2;
                DNMIN = DMX2S[i] - DX;
            }

            double SUM = 1.0;
            isign += NMIN;

            // Triangle check: NMAX must be >= NMIN
            if (NMAX < NMIN) {
                triangle_fail = true;
                break;
            }

            if (NMAX > NMIN) {
                // Compute the Racah sum using Wills' factorized recurrence
                double D1 = DAS[1][i] - (DNMAX - 1);
                double D2 = DAS[2][i] - (DNMAX - 1);
                double D8 = DRS[i] - (DNMAX - 1);
                double D3 = DAS[3][i] - (DNMAX - 1 + DX);
                double D4 = DAS[4][i] - (DNMAX - 1 + DX);
                double D5 = -DNMAX;
                double D6 = -DBS[1][i] - (DX + DNMAX);
                double D7 = -DBS[2][i] - (DX + DNMAX);

                double E1 = D1 * D2 * D3 * D4;
                double E5 = D5 * D6 * D7 * D8;

                int NMAXM2 = NMAX - 2;
                SUM = 1.0 + E1 / E5;

                if (NMAXM2 >= NMIN) {
                    double F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1;
                    double F5 = (D5+1)*(D6+1)*(D7+1)*(D8+1) - E5;
                    double G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1;
                    double G5 = (D5+2)*(D6+2)*(D7+2)*(D8+2) - F5 - F5 - E5;
                    double H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1;
                    double H5 = (D5+3)*(D6+3)*(D7+3)*(D8+3) - 3*(G5+F5) - E5;

                    for (int N = NMIN; N <= NMAXM2; N++) {
                        E1 += F1; E5 += F5;
                        F1 += G1; F5 += G5;
                        G1 += H1; G5 += H5;
                        H1 += 24.0; H5 += 24.0;
                        SUM = 1.0 + SUM * (E1 / E5);
                    }
                }
            }
            // else NMAX == NMIN: SUM stays 1.0

            SUMS[i] = SUM;

            // Accumulate X-dependent log-factorial terms
            if (RS_arr[i] <= _wig9j_maxlf) {
                sumlog += _wig9j_lf(-X + AS[4][i])
                        + _wig9j_lf(X + IDS_arr[1][i])
                        + _wig9j_lf(X + IDS_arr[2][i])
                        - _wig9j_lf(1 + xhalf + X + AS[4][i])
                        + _wig9j_lf(-NMIN + RS_arr[i])
                        - _wig9j_lf(NMIN)
                        - _wig9j_lf(NMIN + X + BS[1][i])
                        - _wig9j_lf(NMIN + X + BS[2][i])
                        - _wig9j_lf(-NMIN + AS[1][i])
                        - _wig9j_lf(-NMIN + AS[2][i])
                        - _wig9j_lf(-NMIN - X + AS[3][i])
                        - _wig9j_lf(-NMIN - X + AS[4][i]);
            } else {
                sumlog += std::lgamma(1.0 - DX + DAS[4][i])
                        + std::lgamma(1.0 + DX + DAS[2][i] + DBS[2][i])
                        + std::lgamma(1.0 + DX + DAS[1][i] + DBS[1][i])
                        - std::lgamma(2.0 + dxhalf + DX + DAS[4][i])
                        + std::lgamma(1.0 - DNMIN + DRS[i])
                        - std::lgamma(1.0 + DNMIN)
                        - std::lgamma(1.0 + DNMIN + DX + DBS[1][i])
                        - std::lgamma(1.0 + DNMIN + DX + DBS[2][i])
                        - std::lgamma(1.0 - DNMIN + DAS[1][i])
                        - std::lgamma(1.0 - DNMIN + DAS[2][i])
                        - std::lgamma(1.0 - DNMIN - DX + DAS[3][i])
                        - std::lgamma(1.0 - DNMIN - DX + DAS[4][i]);
            }
        }

        if (triangle_fail) return 0.0;

        // Multiply the 3 6-J parts together with phase and (2x+1) weight
        double sump = SUMS[1] * SUMS[2] * std::exp(sumlog) * SUMS[3]
                     * (1.0 + dxhalf + 2.0 * DX);

        if (isign & 1) sump = -sump;
        result += sump;
    }

    return result;
}

// Wrapper: takes j values (double, possibly half-integer) instead of 2j
static double WIG9J(double j1, double j2, double j3,
                    double j4, double j5, double j6,
                    double j7, double j8, double j9) {
    return WIG9J_int((int)std::round(2*j1), (int)std::round(2*j2), (int)std::round(2*j3),
                     (int)std::round(2*j4), (int)std::round(2*j5), (int)std::round(2*j6),
                     (int)std::round(2*j7), (int)std::round(2*j8), (int)std::round(2*j9));
}

#endif // WIG9J_H
