// coulin.cpp — Coulomb integral routines ported from Ptolemy Fortran
// RCASYM + CLINTS + COULIN
// Computes ∫(R,∞) F_out(r)*F_in(r)/r^N dr by Belling's method + L-recursion

#include "coulin.h"
#include "rcwfn.h"
#include "math_utils.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>

static const double PI = 3.14159265358979323846;
static const double ACCUR_RCWFN = 1.0e-14;

namespace {

struct CoulombPoint {
    double F = 0.0, Fp = 0.0, G = 0.0, Gp = 0.0;
    double Z = 0.0, dZ = 0.0, phase = 0.0;
};

struct AsymState {
    bool active = false;
    double rho0I = 0.0, rho0F = 0.0;
    int nasyI = 0, nasyF = 0;
    std::vector<double> zI, zF;
};

int EvalCoulombDirect(double rho, double eta, int L, CoulombPoint &out) {
    std::vector<double> FC, FCP, GC, GCP;
    int iret = Rcwfn(rho, eta, L, L, FC, FCP, GC, GCP, ACCUR_RCWFN);
    if (iret != 0) return iret;
    out.F = FC[L];
    out.Fp = FCP[L];
    out.G = GC[L];
    out.Gp = GCP[L];
    out.Z = 1.0 / (out.F * out.F + out.G * out.G);
    out.dZ = 0.0;
    out.phase = 0.0;
    return 0;
}

int EvalCoulombAsym(double rho, double eta, int L, double sigL,
                    double accura, int nmax,
                    CoulombPoint &out, std::vector<double> &zCoeffs, int &ntz) {
    double zeta, phi, dzeta, F, Fp, G, Gp;
    int iret = Rcasym(L, eta, rho, sigL, zeta, phi, dzeta, F, Fp, G, Gp,
                      zCoeffs, accura, nmax, ntz);
    if (iret != 0) return iret;
    out.F = F;
    out.Fp = Fp;
    out.G = G;
    out.Gp = Gp;
    out.Z = zeta;
    out.dZ = dzeta;
    out.phase = phi;
    return 0;
}

void EvalFromAsymSeries(double rho, double eta, int L, double sigL,
                        double rho0, const std::vector<double> &coeffs, int nasy,
                        double accura, CoulombPoint &out) {
    double RI = rho0 / rho;
    double Z = 1.0;
    double dZ = 0.0;
    double PHS = 0.0;
    double WI = RI;
    double DI = 1.0;
    double term2 = 1.0;

    for (int i = 2; i <= nasy; i++) {
        double term = WI * coeffs[i - 1];
        WI *= RI;
        Z += term;
        dZ -= DI * term;
        if (i > 2) PHS += term / (DI - 1.0);
        DI += 1.0;
        if (std::fabs(term) + term2 < accura * Z) {
            break;
        }
        term2 = std::fabs(term);
    }

    PHS = rho - eta * std::log(2.0 * rho) + sigL - 0.5 * PI * L - rho0 * PHS;
    dZ /= rho;

    out.Z = Z;
    out.dZ = dZ;
    out.phase = PHS;
    out.F = std::sin(PHS) / std::sqrt(Z);
    out.G = std::cos(PHS) / std::sqrt(Z);
    out.Fp = Z * out.G - (0.5 * dZ / Z) * out.F;
    out.Gp = -Z * out.F - (0.5 * dZ / Z) * out.G;
}

} // namespace

// ============================================================
// RCASYM: Asymptotic series in 1/rho for Coulomb functions
// Fortran: source_annotated.f lines 30474-30644
// Z_coeffs is output array of series coefficients (resized internally)
// Returns 0 on success, -5 on convergence failure
// ============================================================
int Rcasym(int L, double eta, double rho, double sigL,
           double &zeta, double &phi, double &dzeta,
           double &F, double &Fp, double &G, double &Gp,
           std::vector<double> &Z_coeffs,
           double eps, int nmax, int &ntz) {

    double FLL = (double)L * (L + 1);
    double etasq = eta * eta;
    double rhosq = rho * rho;
    double BF1 = 2.0 * (eta / rho);
    double BF2 = FLL / rhosq;

    // Turning point check
    double rhot;
    if (eta > 0)       rhot = eta * (1.0 + std::sqrt(1.0 + FLL / etasq));
    else if (eta == 0) rhot = (L == 0) ? 0.5 : std::sqrt(FLL);
    else               rhot = std::fabs(eta) * (std::sqrt(1.0 + FLL / etasq) - 1.0);

    if (rho <= rhot) {
        return -5;  // inside turning point
    }

    // Allocate work arrays (1-based in Fortran, 0-based here)
    Z_coeffs.resize(nmax + 1, 0.0);
    std::vector<double> S(nmax + 1, 0.0);
    std::vector<double> DZSQ(nmax + 1, 0.0);
    std::vector<double> ZINV(nmax + 1, 0.0);

    // Initialize (Fortran index 1 → C++ index 0)
    Z_coeffs[0] = 1.0;
    S[0] = 1.0;
    DZSQ[0] = 0.0;
    ZINV[0] = 1.0;
    zeta = 1.0;
    dzeta = 0.0;
    phi = rho - eta * std::log(2.0 * rho) + sigL - 0.5 * L * PI;
    int KSER = 2;
    ntz = 1;

    // Iterative solution for zeta coefficients
    // Fortran I goes 2..NMAX (1-based), C++ i goes 1..nmax-1 (0-based)
    for (int i = 1; i < nmax; i++) {
        // i corresponds to Fortran I = i+1
        int I_f = i + 1;  // Fortran 1-based index
        Z_coeffs[i] = 0.0;
        int jtop = (I_f + 1) / 2;  // Fortran (I+1)/2
        S[i] = 0.0;
        DZSQ[i] = 0.0;
        double CF = 0.0;

        if (I_f > 2) {
            // S[i] = sum over j=2..jtop of factor*Z[j-1]*Z[I_f-j]
            for (int j = 2; j <= jtop; j++) {
                double factor = 2.0;
                if (j == (I_f - j + 1)) factor = 1.0;
                S[i] += factor * Z_coeffs[j - 1] * Z_coeffs[I_f - j];
            }
        }

        // CF = sum over j=1..I_f-1 of Z[j-1]*S[I_f-j]
        for (int j = 1; j <= I_f - 1; j++) {
            CF += Z_coeffs[j - 1] * S[I_f - j];
        }

        double BF = -BF1 * Z_coeffs[i - 1];  // Z(I-1) → Z_coeffs[i-1]
        if (I_f > 2) BF -= BF2 * Z_coeffs[i - 2];  // Z(I-2) → Z_coeffs[i-2]

        // ZINV computation for I > 5
        if (I_f > 5) {
            int iloc = I_f - 4;  // Fortran ILOC = I-4
            int iloc_c = iloc - 1;  // 0-based
            ZINV[iloc_c] = 0.0;
            for (int j = 1; j <= iloc - 1; j++) {
                ZINV[iloc_c] -= ZINV[j - 1] * Z_coeffs[iloc - j];
            }
        }

        double TF = 0.0;
        if (I_f > 4) {
            // DZSQ sum
            for (int j = 3; j <= jtop; j++) {
                double factor = 2.0 / rhosq;
                if (j == (I_f - j + 1)) factor = 1.0 / rhosq;
                DZSQ[i] += factor * (j - 2) * (I_f - j - 1) * Z_coeffs[j - 2] * Z_coeffs[I_f - j - 1];
            }
            // TF sum
            int jjtop = I_f - 4;
            for (int j = 1; j <= jjtop; j++) {
                TF += ZINV[j - 1] * DZSQ[I_f - j];
            }
        }

        double D2F = 0.0;
        if (I_f > 3) D2F = (I_f - 2) * (I_f - 3) * (Z_coeffs[i - 2] / rhosq);

        Z_coeffs[i] = (-CF + BF + 0.75 * TF - 0.5 * D2F) * 0.5;

        ntz++;
        S[i] += 2.0 * Z_coeffs[i];
        zeta += Z_coeffs[i];

        if (std::fabs(Z_coeffs[i]) < eps) KSER--;

        // dzeta update: dzeta -= (I-2)*Z(I-1)/rho
        dzeta -= (I_f - 2) * Z_coeffs[i - 1] / rho;
        if (I_f == 2 || std::fabs((I_f - 2) * Z_coeffs[i - 1] / rho) >= eps) {
            // don't decrement KSER
        } else {
            KSER--;
        }

        // phi update
        if (I_f > 2) {
            phi -= (Z_coeffs[i] * rho) / (double)(I_f - 2);
        }

        if (KSER <= 0) break;  // converged

        // Check series hasn't started diverging (for i >= 7 in 0-based, I_f >= 8)
        if (I_f >= 8) {
            if (std::fabs(Z_coeffs[i]) > std::fabs(Z_coeffs[i - 6])) {
                return -5;  // series turning
            }
        }
    }

    if (KSER > 0) {
        return -5;  // didn't converge within nmax terms
    }

    // Compute F, G, F', G' from zeta, phi, dzeta
    double rtz = std::sqrt(zeta);
    F = std::sin(phi) / rtz;
    G = std::cos(phi) / rtz;
    double factor = (0.5 * dzeta) / (zeta * zeta);
    Fp = zeta * (G - factor * F);
    Gp = -zeta * (F + factor * G);

    return 0;
}

// ============================================================
// CLINTS: Coulomb integral for one (LI, LF) pair
// Fortran: source_annotated.f lines 7369-7870
// Current port chunk: setup + direct/asymptotic Coulomb evaluation path
// ============================================================
int Clints(double rlower, double etaI, double etaF, double kI, double kF,
           double sigI, double sigF, double accura, double &asmult,
           double &ffint, double &fgint, double &gfint, double &ggint,
           int N, int LI, int LF, int nterms, int npts) {

    std::vector<double> pts, wts;
    GaussLegendre(npts, -1.0, 1.0, pts, wts);

    ffint = fgint = gfint = ggint = 0.0;

    int ISYN = 1;
    int IPIECE = 0;
    int ICHI = 0;

    double EPSILO = asmult * accura;
    double EPSIL = EPSILO;
    double DELTA = PI / std::max(kI, kF);
    double RVAL = rlower;
    double RTURN = std::max((etaI + std::sqrt(etaI * etaI + LI * (LI + 1.0))) / kI,
                            (etaF + std::sqrt(etaF * etaF + LF * (LF + 1.0))) / kF);

    AsymState asym;
    CoulombPoint inPt, outPt;

    double ROI = kI * RVAL;
    double ROF = kF * RVAL;

    if (RVAL >= 1.2 * RTURN) {
        int nmax = std::max(LI, 4 * nterms);
        int iret = EvalCoulombAsym(ROI, etaI, LI, sigI, accura, nmax,
                                   inPt, asym.zI, asym.nasyI);
        if (iret == 0) {
            iret = EvalCoulombAsym(ROF, etaF, LF, sigF, accura, nmax,
                                   outPt, asym.zF, asym.nasyF);
            if (iret == 0) {
                asym.active = true;
                asym.rho0I = ROI;
                asym.rho0F = ROF;
            }
        }
    }

    if (!asym.active) {
        int iret = EvalCoulombDirect(ROI, etaI, LI, inPt);
        if (iret != 0) return iret;
        iret = EvalCoulombDirect(ROF, etaF, LF, outPt);
        if (iret != 0) return iret;
    }

    if (asym.active) {
        EvalFromAsymSeries(ROI, etaI, LI, sigI, asym.rho0I, asym.zI, asym.nasyI, accura, inPt);
        EvalFromAsymSeries(ROF, etaF, LF, sigF, asym.rho0F, asym.zF, asym.nasyF, accura, outPt);
    }

    double TI = inPt.F * inPt.Fp + inPt.G * inPt.Gp;
    double TF = outPt.F * outPt.Fp + outPt.G * outPt.Gp;
    double CHI = kI * inPt.Z + ISYN * kF * outPt.Z;
    double tiny = 1.0e-8 * kI * inPt.Z;

    double CUMC = 0.0, CUMS = 0.0, CIP = 0.0, SIP = 0.0;

    while (true) {
        while (true) {
            bool needNumeric = false;
            if (std::fabs(CHI) < tiny) {
                needNumeric = true;
            } else if (std::fabs(inPt.G) > 1000.0) {
                needNumeric = true;
            } else {
                double F = 1.0 / (std::sqrt(inPt.Z * outPt.Z) * std::pow(RVAL, N));
                double A0 = F / CHI;
                double Rtest = (((double)N / RVAL) + 2.0 * (kI * inPt.Z * std::fabs(TI) + kF * outPt.Z * std::fabs(TF))) / std::fabs(CHI);
                double EPS = EPSIL;
                if (ISYN == -1) EPS *= (1.0 + std::fabs(CUMC / A0));
                double DB = Rtest * Rtest * Rtest;
                if (DB * DB <= EPS) break;
                needNumeric = true;
            }

            if (!needNumeric) break;
            if (ISYN == -1) DELTA = PI / (std::fabs(CHI) + (double)N / RVAL);

            IPIECE++;
            double CI = 0.0, SI = 0.0;

            for (int i = 0; i < npts; i++) {
                double R = RVAL + 0.5 * DELTA + 0.5 * DELTA * pts[i];
                double roi = kI * R;
                double rof = kF * R;
                double zfac = 0.5 * DELTA * wts[i] / std::pow(R, N);

                CoulombPoint pin, pout;
                if (asym.active) {
                    EvalFromAsymSeries(roi, etaI, LI, sigI, asym.rho0I, asym.zI, asym.nasyI, accura, pin);
                    EvalFromAsymSeries(rof, etaF, LF, sigF, asym.rho0F, asym.zF, asym.nasyF, accura, pout);
                } else {
                    int iret = EvalCoulombDirect(roi, etaI, LI, pin);
                    if (iret != 0) return iret;
                    iret = EvalCoulombDirect(rof, etaF, LF, pout);
                    if (iret != 0) return iret;
                }

                if (ISYN == 1) {
                    ffint += zfac * pout.F * pin.F;
                    fgint += zfac * pout.F * pin.G;
                    gfint += zfac * pout.G * pin.F;
                    ggint += zfac * pout.G * pin.G;
                } else {
                    if (asym.active) {
                        CI += (zfac / std::sqrt(pin.Z * pout.Z)) * std::cos(pin.phase - pout.phase);
                        SI += (zfac / std::sqrt(pin.Z * pout.Z)) * std::sin(pin.phase - pout.phase);
                    } else {
                        CI += zfac * (pin.G * pout.G + pin.F * pout.F);
                        SI += zfac * (pin.F * pout.G - pin.G * pout.F);
                    }
                }
            }

            RVAL += DELTA;

            double ROI2 = kI * RVAL;
            double ROF2 = kF * RVAL;
            if (asym.active) {
                EvalFromAsymSeries(ROI2, etaI, LI, sigI, asym.rho0I, asym.zI, asym.nasyI, accura, inPt);
                EvalFromAsymSeries(ROF2, etaF, LF, sigF, asym.rho0F, asym.zF, asym.nasyF, accura, outPt);
            } else {
                int iret = EvalCoulombDirect(ROI2, etaI, LI, inPt);
                if (iret != 0) return iret;
                iret = EvalCoulombDirect(ROF2, etaF, LF, outPt);
                if (iret != 0) return iret;
            }

            TI = inPt.F * inPt.Fp + inPt.G * inPt.Gp;
            TF = outPt.F * outPt.Fp + outPt.G * outPt.Gp;
            CHI = kI * inPt.Z + ISYN * kF * outPt.Z;
            tiny = 1.0e-8 * kI * inPt.Z;

            if (ISYN == -1) {
                CUMC += CI;
                CUMS += SI;
                if (std::fabs(CI) < accura * std::fabs(ffint + 0.5 * (CUMC - CIP)) &&
                    std::fabs(SI) < accura * std::fabs(gfint + 0.5 * (CUMS + SIP))) {
                    ffint += 0.5 * (CUMC - CIP);
                    ggint += 0.5 * (CUMC + CIP);
                    fgint += 0.5 * (SIP - CUMS);
                    gfint += 0.5 * (SIP + CUMS);
                    return 0;
                }
            }

            if (IPIECE > 2000) return -1;
        }

        double PHSI = std::atan2(inPt.F, inPt.G);
        double PHSF = std::atan2(outPt.F, outPt.G);
        double DZI = -2.0 * kI * inPt.Z * inPt.Z * TI;
        double DZF = -2.0 * kF * outPt.Z * outPt.Z * TF;
        double DCHI = kI * DZI + ISYN * kF * DZF;
        double D1 = DZI / inPt.Z + DZF / outPt.Z;
        double Famp = 1.0 / (std::sqrt(inPt.Z * outPt.Z) * std::pow(RVAL, N));
        double C1 = (double)N / RVAL + 0.5 * D1;
        double DF = -Famp * C1;
        double WI = 1.0 - (2.0 * etaI) / (kI * RVAL) - (LI * (LI + 1.0)) / ((kI * RVAL) * (kI * RVAL));
        double WF = 1.0 - (2.0 * etaF) / (kF * RVAL) - (LF * (LF + 1.0)) / ((kF * RVAL) * (kF * RVAL));
        double DTI = kI * (inPt.Fp * inPt.Fp + inPt.Gp * inPt.Gp - WI / inPt.Z);
        double DTF = kF * (outPt.Fp * outPt.Fp + outPt.Gp * outPt.Gp - WF / outPt.Z);
        double D2ZI = -2.0 * kI * inPt.Z * (2.0 * DZI * TI + inPt.Z * DTI);
        double D2ZF = -2.0 * kF * outPt.Z * (2.0 * DZF * TF + outPt.Z * DTF);
        double D2CHI = kI * D2ZI + ISYN * kF * D2ZF;
        double D12 = (DZI / inPt.Z) * (DZI / inPt.Z) + (DZF / outPt.Z) * (DZF / outPt.Z);
        double D2 = D2ZI / inPt.Z + D2ZF / outPt.Z;
        double D2F = (Famp / RVAL - DF) * C1 - (Famp * D1) / (2.0 * RVAL) - 0.5 * Famp * (D2 - D12);

        double TIMS = 1.0 / (RVAL * CHI);
        double PHANG = PHSI + ISYN * PHSF;
        double SN = std::sin(PHANG);
        double CSN = std::cos(PHANG);
        double A0 = Famp / CHI;
        double CILAST = 0.0, SILAST = 0.0;
        double relEven = 0.0, relOdd = 0.0;
        bool converged = false;

        for (int NIT = 1; NIT <= 2; NIT++) {
            double FACTOR = Famp / CHI;
            double EVEN = A0;
            double AM1 = A0;
            double B = 1.0;
            double ODD = 0.0;
            double SIGN = 1.0;
            bool oddsw = true;
            double DB = 0.0;
            double AM2 = 0.0;
            double CI = 0.0, SI = 0.0;
            bool localConv = false;

            for (int I = 1; I <= nterms; I++) {
                FACTOR *= TIMS;
                double C2 = DF / Famp - (I * DCHI) / CHI;
                double BM1 = B;
                B = RVAL * DB + (C2 * RVAL - I + 1.0) * B;
                if (NIT > 1) {
                    DB = (C2 * RVAL - I + 2.0) * DB +
                         (C2 - RVAL * ((DF / Famp) * (DF / Famp) - I * (DCHI / CHI) * (DCHI / CHI)) +
                          RVAL * (D2F / Famp - I * D2CHI / CHI)) * BM1;
                }
                double A = B * FACTOR;

                if (oddsw) {
                    ODD += SIGN * A;
                    SIGN = -SIGN;
                } else {
                    EVEN += SIGN * A;
                }
                oddsw = !oddsw;

                if (I > 1) {
                    if (std::fabs(A) > std::fabs(AM2)) {
                        localConv = false;
                        break;
                    }
                    if (std::fabs(A) < accura * (std::fabs(A0) + std::min(std::fabs(CUMC), std::fabs(CUMS)))) {
                        CI = -EVEN * SN - ODD * CSN;
                        SI = -ODD * SN + EVEN * CSN;
                        if (NIT == 2) {
                            relEven = std::fabs((CILAST - CI) / (CUMC + CI));
                            relOdd = std::fabs((SILAST - SI) / (CUMS + SI));
                        }
                        CILAST = CI;
                        SILAST = SI;
                        localConv = true;
                        break;
                    }
                }

                AM2 = AM1;
                AM1 = A;
                if (I == nterms) localConv = false;
            }

            if (!localConv) {
                converged = false;
                break;
            }
            converged = true;
        }

        if (!converged || std::max(relEven, relOdd) * std::max(relEven, relOdd) > accura) {
            // Fortran label 700: reduce ASMULT
            double Rtest = (((double)N / RVAL) + 2.0 * (kI * inPt.Z * std::fabs(TI) + kF * outPt.Z * std::fabs(TF))) / std::fabs(CHI);
            double EPS = EPSIL;
            double Rnew = std::min(0.5 * EPSIL, 0.75 * std::pow(Rtest, 6.0) * (EPSIL / EPS));
            asmult = asmult * Rnew / EPSIL;
            EPSIL = Rnew;
            ICHI++;
            if (ICHI > 20) return -1;
            continue;
        }

        if (ISYN == -1) {
            ffint += 0.5 * (CUMC - CILAST);
            ggint += 0.5 * (CUMC + CILAST);
            fgint += 0.5 * (SILAST - CUMS);
            gfint += 0.5 * (SILAST + CUMS);
            return 0;
        }

        CIP = CILAST;
        SIP = SILAST;
        ISYN = -1;
        IPIECE = 0;
        CUMS = 0.0;
        CUMC = 0.0;
        EPSIL = EPSILO;
        CHI = kI * inPt.Z + ISYN * kF * outPt.Z;
        tiny = 1.0e-8 * kI * inPt.Z;
    }
}

// ============================================================
// COULIN: Coulomb integrals by L-recursion for all (LIN,LOUT) pairs
// Fortran: source_annotated.f lines 9319-9797
// Current port chunk: storage layout + CLINTS seeding + R=0 inward FF accumulation
// ============================================================
int Coulin(int N, int maxdel, int lmin, int lmax,
           double etaOut, double kOut, const std::vector<double> &sigOut,
           double etaIn, double kIn, const std::vector<double> &sigIn,
           double R, bool allSw,
           CoulinResult &result,
           double accura, double asmult, int nterms, int npts) {

    int MXDEL = std::max(2, maxdel);
    int LMAXMX = lmax + MXDEL;
    int LMINMN = std::max(0, lmin - MXDEL);
    if (result.ldldim <= MXDEL) return -10;

    int NTRM = std::max(nterms, lmax / 4);
    double ACC = accura * std::exp(10.0 - 5.0 * N - 2.0 * std::fabs(etaIn - etaOut)
                                   - (double)(lmax - lmin) / 100.0);
    ACC = std::max(ACC, std::max(1.0e-12, accura * 0.01));

    double RMINS[2] = {R, R};
    double ACCS[2] = {ACC, 0.01 * accura};
    double RTURNS[2] = {0.0, 0.0};

    if (R <= 0.0) {
        ACCS[0] = 0.1 * accura;
        ACCS[1] = 0.1 * ACC;
        if (allSw) return -12;

        RTURNS[0] = std::max((etaIn + std::sqrt(etaIn * etaIn + LMINMN * (LMINMN + 1.0))) / kIn,
                             (etaOut + std::sqrt(etaOut * etaOut + LMINMN * (LMINMN + 1.0))) / kOut);
        RTURNS[1] = std::max((etaIn + std::sqrt(etaIn * etaIn + LMAXMX * (LMAXMX + 1.0))) / kIn,
                             (etaOut + std::sqrt(etaOut * etaOut + LMAXMX * (LMAXMX + 1.0))) / kOut);
        RMINS[0] = 1.4 * RTURNS[0];
        RMINS[1] = 1.4 * RTURNS[1];
    }

    result.nils = 2 * (lmax - lmin) + 1;
    result.FF.assign(result.ldldim * result.nils, 0.0);
    result.FG.assign(result.ldldim * result.nils, 0.0);
    result.GF.assign(result.ldldim * result.nils, 0.0);
    result.GG.assign(result.ldldim * result.nils, 0.0);

    std::vector<double> starts(result.ldldim * 2 * 4 * 2, 0.0);
    auto SIDX = [&](int id, int is, int comp, int ii) -> double& {
        return starts[(((id * 2 + is) * 4 + comp) * 2 + ii)];
    };
    auto FFAT = [&](int id, int il) -> double& { return result.FF[id + il * result.ldldim]; };
    auto FGAT = [&](int id, int il) -> double& { return result.FG[id + il * result.ldldim]; };
    auto GFAT = [&](int id, int il) -> double& { return result.GF[id + il * result.ldldim]; };
    auto GGAT = [&](int id, int il) -> double& { return result.GG[id + il * result.ldldim]; };
    // Safe read-only accessor: returns 0.0 for out-of-bounds IL
    // Fortran equivalent: FF(id, 0) accesses memory before array start, which is 0
    // in the seeded region (zero-initialized). C++ needs explicit bounds guard.
    auto FF_safe = [&](int id, int il) -> double {
        if (il < 0 || il >= result.nils || id < 0 || id >= result.ldldim) return 0.0;
        return result.FF[id + il * result.ldldim];
    };
    auto FG_safe = [&](int id, int il) -> double {
        if (il < 0 || il >= result.nils || id < 0 || id >= result.ldldim) return 0.0;
        return result.FG[id + il * result.ldldim];
    };
    auto GF_safe = [&](int id, int il) -> double {
        if (il < 0 || il >= result.nils || id < 0 || id >= result.ldldim) return 0.0;
        return result.GF[id + il * result.ldldim];
    };
    auto GG_safe = [&](int id, int il) -> double {
        if (il < 0 || il >= result.nils || id < 0 || id >= result.ldldim) return 0.0;
        return result.GG[id + il * result.ldldim];
    };

    int LS = 2 * lmin;
    int ISTOPR = std::max(N - 2 * lmin, 2);
    if (allSw) ISTOPR = std::max(MXDEL + 1 - 2 * lmin, ISTOPR);
    int ISTOP = ISTOPR;

    for (int II = 0; II < 2; II++) {
        for (int I = 1; I <= ISTOP; I++) {
            int IL = LS - 2 * lmin + I;
            int MD = MXDEL + ((MXDEL + IL) % 2);
            int IS = (I == ISTOP) ? 1 : 0;

            for (int ID = 1; ID <= MD; ID++) {
                int LIN = (LS + I + MXDEL) / 2 + ID - 1 - MXDEL;
                int LOUT = LS + I - 1 - LIN;
                if (LIN < 0 || LOUT < 0) continue;

                double A = asmult;
                double ff = 0.0, fg = 0.0, gf = 0.0, gg = 0.0;
                int iret = Clints(RMINS[II], etaIn, etaOut, kIn, kOut,
                                  sigIn.at(LIN), sigOut.at(LOUT),
                                  ACCS[II], A, ff, fg, gf, gg,
                                  N, LIN, LOUT, NTRM, npts);
                if (iret != 0) return iret;

                SIDX(ID - 1, IS, 0, II) = ff;
                SIDX(ID - 1, IS, 1, II) = fg;
                SIDX(ID - 1, IS, 2, II) = gf;
                SIDX(ID - 1, IS, 3, II) = gg;
            }

            for (int ID = 1; ID <= MD; ID++) {
                FFAT(ID - 1, IL - 1) = SIDX(ID - 1, IS, 0, II);
                if (!allSw) continue;
                FGAT(ID - 1, IL - 1) = SIDX(ID - 1, IS, 1, II);
                GFAT(ID - 1, IL - 1) = SIDX(ID - 1, IS, 2, II);
                GGAT(ID - 1, IL - 1) = SIDX(ID - 1, IS, 3, II);
            }
        }
        LS = 2 * lmax - 1;
        ISTOP = 2;
    }

    if (R == 0.0) {
        double DELTA = 3.1 / kIn;
        ISTOP = ISTOPR;
        int LMN = LMINMN;
        int LMX = lmin + (ISTOP + MXDEL + 1) / 2;
        LS = 2 * lmin;

        std::vector<double> pts, wts;
        GaussLegendre(npts, -1.0, 1.0, pts, wts);

        for (int II = 0; II < 2; II++) {
            double RTOP = RMINS[II];
            while (true) {
                double RBOT = std::max(0.0, RTOP - DELTA);
                if (RBOT >= RTOP) break;
                double B = RTOP - RBOT;
                double lastR = RBOT;
                bool done = false;

                for (int IPT = 0; IPT < npts; IPT++) {
                    double RVAL = RTOP - 0.5 * B * (1.0 + pts[IPT]);
                    lastR = RVAL;

                    std::vector<double> FI, FIP, GI, GIP;
                    int iret = Rcwfn(kIn * RVAL, etaIn, LMN, LMX, FI, FIP, GI, GIP, 1.0e-8);
                    if (iret != 0 && iret != 2) return iret;
                    if (std::fabs(FI.at(LMN)) < accura && RVAL < RTURNS[II]) {
                        done = true;
                        break;
                    }

                    std::vector<double> FO, FOP, GO, GOP;
                    iret = Rcwfn(kOut * RVAL, etaOut, LMN, LMX, FO, FOP, GO, GOP, 1.0e-8);
                    if (iret != 0 && iret != 2) return iret;
                    if (std::fabs(FO.at(LMN)) < accura && RVAL < RTURNS[II]) {
                        done = true;
                        break;
                    }

                    double X = 0.5 * B * wts[IPT] / std::pow(RVAL, N);
                    for (int I = 1; I <= ISTOP; I++) {
                        int IL = LS - 2 * lmin + I;
                        int MD = MXDEL + ((MXDEL + IL) % 2);
                        int IS = (I == ISTOP) ? 1 : 0;
                        for (int ID = 1; ID <= MD; ID++) {
                            int LIN = (LS + I + MXDEL) / 2 + ID - 1 - MXDEL;
                            int LOUT = LS + I - 1 - LIN;
                            if (LIN < 0 || LOUT < 0) continue;
                            double A = X * FI.at(LIN) * FO.at(LOUT);
                            FFAT(ID - 1, IL - 1) += A;
                            SIDX(ID - 1, IS, 0, II) = FFAT(ID - 1, IL - 1);
                        }
                    }
                }

                if (done) {
                    RMINS[II] = lastR;
                    LS = 2 * lmax - 1;
                    ISTOP = 2;
                    LMN = (LS - MXDEL - 1) / 2;
                    LMX = lmax + MXDEL;
                    break;
                }
                RTOP = RBOT;
            }
        }
    }

    int LSUMMN = 2 * lmin + 1 + ISTOPR - 2;
    int LSUMMX = 2 * lmax - 1;
    int MMXDEL = -MXDEL;
    double DN = (double)N;

    // __float128 arrays for quad-precision R=0 downward recursion
    int nFF = result.ldldim * result.nils;
    std::vector<__float128> FF128(nFF, (__float128)0.0);
    auto FFAT128 = [&](int id, int il) -> __float128& { return FF128[id + il * result.ldldim]; };
    auto FF_safe128 = [&](int id, int il) -> __float128 {
        if (il < 0 || il >= result.nils || id < 0 || id >= result.ldldim) return (__float128)0.0;
        return FF128[id + il * result.ldldim];
    };
    // Sync FF128 from result.FF (inner loop has modified result.FF)
    for (int i = 0; i < nFF; i++) FF128[i] = (__float128)result.FF[i];

    // STABILITY FIX: Save pre-recursion seed+inner values for il=0,1 (the ISTOPR=2 seed rows).
    // These come from Clints + inner loop and are accurate at double precision.
    // The downward recursion overwrites them with values propagated from top seeds,
    // which can be wrong due to floating-point amplification.
    // After recursion, restore il=0,1 to preserve the accurate seed values.
    std::vector<double> saved_seeds;
    if (R == 0.0) {
        int ISTOPR_save = std::max(N - 2 * lmin, 2);
        if (allSw) ISTOPR_save = std::max(MXDEL + 1 - 2 * lmin, ISTOPR_save);
        for (int il = 0; il < ISTOPR_save && il < result.nils; il++)
            for (int id = 0; id < result.ldldim; id++)
                saved_seeds.push_back(result.FF[id + il * result.ldldim]);
    }

    if (R == 0.0) {
        // Recurse downward for FF
        double DLI = 0, DLO = 0;
        int ID = 0, ID2 = 0, IL = 0;

        for (int LS2 = LSUMMN; LS2 <= LSUMMX; LS2++) {
            int LSUM = LSUMMX + LSUMMN - LS2;
            for (int LDEL = MMXDEL; LDEL <= MXDEL; LDEL++) {
                int LIN = LSUM - LDEL;
                if (LIN < 0) continue;
                LIN = LIN / 2;
                int LOUT = (LSUM + LDEL) / 2;
                if (LSUM != LIN + LOUT) continue;
                int LOUTP = LOUT - 1;
                if (std::abs(LIN - LOUTP) > MXDEL) continue;
                if (LOUTP < 0) continue;
                if (std::abs(LOUT + 1 - LIN) > MXDEL) continue;

                DLI = (double)LIN;
                DLO = (double)LOUT;
                ID = (LIN - LOUT + MXDEL) + 2;
                IL = LIN + LOUT - 2 * lmin + 1;
                ID2 = (ID + 1) / 2;
                ID = ID / 2;

                double E = (2.0 * DLO + 1.0) / ((DLI + DLO - DN + 2.0) * std::sqrt(DLO * DLO + etaOut * etaOut));
                double A = etaOut * (DLI - DN + 2.0) / (DLO + 1.0)
                           - (etaIn * (kIn / kOut)) * DLO / (DLI + 1.0);
                double B = (kIn / kOut) * (DLO / (DLI + 1.0)) * std::sqrt((DLI + 1.0) * (DLI + 1.0) + etaIn * etaIn);
                double C = DLO * (DLO - DLI + DN - 1.0) * std::sqrt((DLO + 1.0) * (DLO + 1.0) + etaOut * etaOut)
                           / ((DLO + 1.0) * (2.0 * DLO + 1.0));

                // Recurse in __float128 (quad precision) for numerical stability
                __float128 E128 = (__float128)E, A128 = (__float128)A, B128 = (__float128)B, C128 = (__float128)C;
                FFAT128(ID2-1, IL-2) = E128*(A128*FFAT128(ID-1,IL-1) + C128*FF_safe128(ID2-2,IL)
                                           + B128*FF_safe128(ID2-1,IL));
                // Also update double array for reading by other parts
                FFAT(ID2-1, IL-2) = (double)FFAT128(ID2-1, IL-2);
            }

            // Reduce LIN for the last case
            if (DLI <= 0.0) continue;

            double E = (2.0*DLI+1.0)/((DLI+DLO-DN+2.0)*std::sqrt(DLI*DLI+etaIn*etaIn));
            double A = etaIn*(DLO-DN+2.0)/(DLI+1.0) - (etaOut*(kOut/kIn))*DLI/(DLO+1.0);
            double B = (kOut/kIn)*(DLI/(DLO+1.0))*std::sqrt((DLO+1.0)*(DLO+1.0)+etaOut*etaOut);
            double C = DLI*(DLI-DLO+DN-1.0)*std::sqrt((DLI+1.0)*(DLI+1.0)+etaIn*etaIn)
                       /((DLI+1.0)*(2.0*DLI+1.0));
            __float128 E2=(__float128)E,A2=(__float128)A,B2=(__float128)B,C2=(__float128)C;
            FFAT128(ID2-2, IL-2) = E2*(A2*FFAT128(ID-1,IL-1) + B2*FF_safe128(ID2-2,IL)
                                      + C2*FF_safe128(ID2-1,IL));
            FFAT(ID2-2, IL-2) = (double)FFAT128(ID2-2, IL-2);
        }
        // Copy __float128 result back to double array
        for (int i = 0; i < nFF; i++) result.FF[i] = (double)FF128[i];
        // STABILITY FIX: Restore il=0,1 seed rows with accurate pre-recursion values.
        // If the recursion reproduced them correctly, this is a no-op.
        // If the recursion gave wrong values (due to top-seed noise amplification),
        // this restores the accurate Clints+inner-loop values.
        if (!saved_seeds.empty()) {
            int ISTOPR_save = std::max(N - 2 * lmin, 2);
            if (allSw) ISTOPR_save = std::max(MXDEL + 1 - 2 * lmin, ISTOPR_save);
            int k = 0;
            for (int il = 0; il < ISTOPR_save && il < result.nils; il++)
                for (int id = 0; id < result.ldldim; id++)
                    result.FF[id + il * result.ldldim] = saved_seeds[k++];
        }
    } else {
        // R > 0: get inhomogeneous terms (Coulomb wavefunctions at R)
        std::vector<double> FI, FIP, GI_wf, GIP;
        int iret = Rcwfn(kIn * R, etaIn, LMINMN, LMAXMX, FI, FIP, GI_wf, GIP, 1.0e-14);
        if (iret != 0) return iret;

        std::vector<double> FO, FOP, GO_wf, GOP;
        iret = Rcwfn(kOut * R, etaOut, LMINMN, LMAXMX, FO, FOP, GO_wf, GOP, 1.0e-14);
        if (iret != 0) return iret;

        // Recurse upwards on FG, GF, GG, and FF when R > 0
        double DLI = 0, DLO = 0;
        int LINLS = 0, LOUTLS = 0;
        int ID = 0, ID2 = 0, IL = 0;

        for (int LSUM = LSUMMN; LSUM <= LSUMMX; LSUM++) {
            for (int LDEL = MMXDEL; LDEL <= MXDEL; LDEL++) {
                int LIN = (LSUM + LDEL) / 2;
                int LOUT = (LSUM - LDEL) / 2;
                if (LSUM != LIN + LOUT) continue;
                int LOUTP = LOUT + 1;
                if (std::abs(LIN - LOUTP) > MXDEL) continue;
                if (std::abs(LOUT - 1 - LIN) > MXDEL) continue;

                DLI = (double)LIN;
                DLO = (double)LOUT;
                LINLS = LIN;
                LOUTLS = LOUT;
                ID = (LIN - LOUT + MXDEL) + 2;
                IL = LIN + LOUT - 2 * lmin + 1;
                ID2 = (ID + 1) / 2;
                ID = ID / 2;

                double E = (2.0 * DLO + 1.0) / ((DLI + DLO + DN) * std::sqrt((DLO + 1.0) * (DLO + 1.0) + etaOut * etaOut));
                double A = etaOut * (DLI + DN - 1.0) / DLO
                           - (etaIn * (kIn / kOut)) * (DLO + 1.0) / DLI;
                double B = (kIn / kOut) * ((DLO + 1.0) / DLI) * std::sqrt(DLI * DLI + etaIn * etaIn);
                double C2 = (DLO + 1.0) * (DLO - DLI - DN + 1.0) * std::sqrt(DLO * DLO + etaOut * etaOut)
                            / (DLO * (2.0 * DLO + 1.0));
                double D = (DLO + 1.0) / (kOut * std::pow(R, N));

                // Fortran: FF(ID2-1, IL+1) — 1-based → 0-based: (ID2-2, IL)
                // IL-2 may be negative at boundary; use safe accessor → returns 0
                FFAT(ID2 - 2, IL) = E * (A * FFAT(ID - 1, IL - 1) + C2 * FF_safe(ID2 - 1, IL - 2)
                                         + B * FF_safe(ID2 - 2, IL - 2) + D * FO[LOUT] * FI[LIN]);
                if (!allSw) continue;
                FGAT(ID2 - 2, IL) = E * (A * FGAT(ID - 1, IL - 1) + C2 * FG_safe(ID2 - 1, IL - 2)
                                         + B * FG_safe(ID2 - 2, IL - 2) + D * FO[LOUT] * GI_wf[LIN]);
                GFAT(ID2 - 2, IL) = E * (A * GFAT(ID - 1, IL - 1) + C2 * GF_safe(ID2 - 1, IL - 2)
                                         + B * GF_safe(ID2 - 2, IL - 2) + D * GO_wf[LOUT] * FI[LIN]);
                GGAT(ID2 - 2, IL) = E * (A * GGAT(ID - 1, IL - 1) + C2 * GG_safe(ID2 - 1, IL - 2)
                                         + B * GG_safe(ID2 - 2, IL - 2) + D * GO_wf[LOUT] * GI_wf[LIN]);
            }

            // Increase LIN for last case to complete generation of LSUM+1 from LSUM
            double E = (2.0 * DLI + 1.0) / ((DLI + DLO + DN) * std::sqrt((DLI + 1.0) * (DLI + 1.0) + etaIn * etaIn));
            double A = etaIn * (DLO + DN - 1.0) / DLI
                       - (etaOut * (kOut / kIn)) * (DLI + 1.0) / DLO;
            double B = (kOut / kIn) * ((DLI + 1.0) / DLO) * std::sqrt(DLO * DLO + etaOut * etaOut);
            double C2 = (DLI + 1.0) * (DLI - DLO - DN + 1.0) * std::sqrt(DLI * DLI + etaIn * etaIn)
                        / (DLI * (2.0 * DLI + 1.0));
            double D = (DLI + 1.0) / (kIn * std::pow(R, N));

            FFAT(ID2 - 1, IL) = E * (A * FFAT(ID - 1, IL - 1) + B * FF_safe(ID2 - 1, IL - 2)
                                     + C2 * FF_safe(ID2 - 2, IL - 2) + D * FO[LOUTLS] * FI[LINLS]);
            if (allSw) {
                FGAT(ID2 - 1, IL) = E * (A * FGAT(ID - 1, IL - 1) + B * FG_safe(ID2 - 1, IL - 2)
                                         + C2 * FG_safe(ID2 - 2, IL - 2) + D * FO[LOUTLS] * GI_wf[LINLS]);
                GFAT(ID2 - 1, IL) = E * (A * GFAT(ID - 1, IL - 1) + B * GF_safe(ID2 - 1, IL - 2)
                                         + C2 * GF_safe(ID2 - 2, IL - 2) + D * GO_wf[LOUTLS] * FI[LINLS]);
                GGAT(ID2 - 1, IL) = E * (A * GGAT(ID - 1, IL - 1) + B * GG_safe(ID2 - 1, IL - 2)
                                         + C2 * GG_safe(ID2 - 2, IL - 2) + D * GO_wf[LOUTLS] * GI_wf[LINLS]);
            }
        }
    }

    return 0;
}
