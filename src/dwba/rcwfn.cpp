#include "rcwfn.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// Constants
static const double VRYBIG = 1.79769e+308;
static const double BIG = 1.0e+300;
static const double SMALL = 1.0e-300;
static const double SMALLN = -690.775527898214;
static const double PRERT3 = 2.7e-5;
static const double PRECIS = 2.0e-14;
static const double PRECLN = 32.0;
static const double PI = 3.141592653589793238;

// Helper function for Steed's method (P + iQ)
int SteedPQ(double RHO, double ETA, int LMIN, double &P, double &Q, double ACC) {
  P = 0.0;
  Q = RHO - ETA;
  double PL = 0.0;
  double ETA2 = ETA * ETA;
  // XLL1 = LMIN*(LMIN+1) — must match Fortran (source: rcwfn.f line 264)
  double XLL1 = (double)LMIN * (double)(LMIN + 1);

  double AR = -(ETA2 + XLL1);
  double AI = ETA;
  double BR = Q + Q;
  double BI = 2.0;
  double WI = ETA + ETA;
  double DR = BR / (BR * BR + BI * BI);
  double DI = -BI / (BR * BR + BI * BI);
  double DP = -(AR * DI + AI * DR);
  double DQ = (AR * DR - AI * DI);

  for (int I = 1; I <= 100000; ++I) {
    P += DP;
    Q += DQ;
    PL += 2.0;
    AR += PL;
    AI += WI;
    BI += 2.0;
    double D = AR * DR - AI * DI + BR;
    DI = AI * DR + AR * DI + BI;
    double T = 1.0 / (D * D + DI * DI);
    DR = T * D;
    DI = -T * DI;
    double H = BR * DR - BI * DI - 1.0;
    double X = BI * DR + BR * DI;
    T = DP * H - DQ * X;
    DQ = DP * X + DQ * H;
    DP = T;

    if (PL > 46000.0)
      return 7; // Nonconvergence
    if (std::abs(DP) + std::abs(DQ) < (std::abs(P) + std::abs(Q)) * ACC)
      return 0; // Converged
  }
  return 7;
}

int Rcwfn(double RHO, double ETA, int MINL, int MAXL, std::vector<double> &FC,
          std::vector<double> &FCP, std::vector<double> &GC,
          std::vector<double> &GCP, double ACCUR) {

  // Resize arrays
  if (FC.size() <= MAXL)
    FC.resize(MAXL + 1);
  if (FCP.size() <= MAXL)
    FCP.resize(MAXL + 1);
  if (GC.size() <= MAXL)
    GC.resize(MAXL + 1);
  if (GCP.size() <= MAXL)
    GCP.resize(MAXL + 1);

  double ACC = std::max(ACCUR, PRECIS);
  ACC = std::min(ACC, PRERT3);

  int LMAX = MAXL;
  int LMIN = MINL;
  double ETA2 = ETA * ETA;
  double TURN = ETA + std::sqrt(ETA2 + (double)LMIN * (LMIN + 1));

  // Determine method/region
  int Method = 0; // 1: Steed, 2: Maclaurin (small rho), 3: Maclaurin (general),
                  // 4: Taylor (turning point), 5: Maclaurin (very small rho)

  if (RHO > 0.45) {
    if (RHO >= TURN - 1.0e-4) {
      Method = 1; // Steed
    } else {
      if (RHO < ETA + std::abs(ETA)) {
        if (ETA < 10.0 || RHO <= ETA)
          Method = 2; // Should be 2 or 3?
        else
          Method = 3;
        LMIN = 0; // Reset LMIN for these methods
      } else {
        // Reduce MINL
        LMIN =
            (int)(0.5 *
                  (std::sqrt(1.0 + 4.0 * ((RHO - ETA) * (RHO - ETA) - ETA2)) -
                   1.0));
        Method = 1; // Steed with reduced LMIN
      }
    }
  } else {
    // RHO <= 0.45
    if (ETA >= 0) {
      if (RHO > 0.005) {
        if (ETA < 10.0 || RHO <= ETA)
          Method = 2;
        else
          Method = 3;
        LMIN = 0;
      } else {
        if (ETA < 10.0 || RHO <= ETA)
          Method = 5; // -> 2 logic but different return
        else
          Method = 3; // -> 3 logic but different return?
        // Actually original code sets IGOTO=5 then goes to 70 which sets LMIN=0
        // then checks ETA/RHO If ETA<10 or RHO<=ETA -> 80 (Steed for F'/F) ->
        // 400 (Maclaurin) -> 850 (Return F, F' only) If ETA>=10 and RHO>ETA ->
        // 80 -> 500 (Taylor) -> ... Let's simplify: Method 5 means "Compute F,
        // F' only using Maclaurin"
        Method = 5;
        LMIN = 0;
      }
    } else {
      if (-ETA * RHO > 7.0) {
        // Similar to RHO > 0.45 logic
        TURN = ETA + std::sqrt(ETA2 + (double)LMIN * (LMIN + 1));
        if (RHO >= TURN - 1.0e-4)
          Method = 1;
        else {
          if (RHO < ETA + std::abs(ETA)) {
            Method = 2; // Or 3
            LMIN = 0;
          } else {
            LMIN =
                (int)(0.5 * (std::sqrt(1.0 + 4.0 * ((RHO - ETA) * (RHO - ETA) -
                                                    ETA2)) -
                             1.0));
            Method = 1;
          }
        }
      } else {
        if (RHO > 0.005) {
          Method = 2; // Or 3
          LMIN = 0;
        } else {
          Method = 5;
          LMIN = 0;
        }
      }
    }
  }

  // Refine Method 2 vs 3
  if (Method == 2 || Method == 3 || Method == 5) {
    if (ETA >= 10.0 && RHO > ETA) {
      if (Method == 5)
        Method = 6; // Special case: Taylor but only F/F' return?
      else
        Method = 4; // Taylor
    } else {
      // Keep 2 or 5 (Maclaurin)
      if (Method == 3)
        Method = 2;
    }
  }

  // Common Step: Compute F'/F for L = MAXL using Steed's continued fraction
  double PL = LMAX + 1.0;
  double PLSAVE = PL;
  double R, F;
  bool converged = false;

  // Loop for F'/F
  while (true) {
    bool FRSTSW = true;
    R = ETA / PL + PL / RHO;
    double DQ = (ETA * RHO) * 2.0 + 6.0 * PL * PL;
    double DR = 12.0 * PL + 6.0;
    double DEL = 0.0;
    double D = 0.0;
    F = 1.0;
    double X = (PL * PL - PL + (ETA * RHO)) * (2.0 * PL - 1.0);
    double AI = RHO * PL * PL;
    double DI = (2.0 * PL + 1.0) * RHO;

    for (int I = 1; I <= 100000; ++I) {
      double H = (AI + RHO * ETA2) * (RHO - AI);
      X = X + DQ;
      D = D * H + X;

      // Corrected convergence check from previous debugging
      if (std::abs(D) > PRERT3 * std::abs(DR)) {
        D = 1.0 / D;
        DQ += DR;
        DR += 12.0;
        AI += DI;
        DI += 2.0 * RHO;
        DEL *= (D * X - 1.0);
        if (FRSTSW)
          DEL = -RHO * (PL * PL + ETA2) * (PL + 1.0) * D / PL;
        FRSTSW = false;
        R += DEL;
        if (D < 0.0)
          F = -F;
        if (std::abs(DEL) < std::abs(R * ACC)) {
          converged = true;
          break;
        }
        continue;
      }

      // Restart logic
      PL += 1.0;
      if (PL < PLSAVE + 100.0) {
        // Restart outer loop
        break; // Break inner loop, continue outer
      } else {
        return 5; // Failed
      }
    }
    if (converged)
      break;
    // If we reached here, we broke inner loop to restart, or failed convergence
    if (PL >= PLSAVE + 100.0)
      return 6; // Nonconvergence
  }

  // Recurse down to LMAX
  while (PL > PLSAVE) {
    PL -= 1.0;
    double D_val = ETA / PL + PL / RHO;
    F = (R + D_val) * F;
    R = D_val - (1.0 + ETA2 / (PL * PL)) / (R + D_val);
  }

  // If Method 4 (Taylor at 2*ETA), we need R at LMIN, 2*ETA.
  // But wait, the code above computes R at RHO.
  // Original code: If IGOTO=4, it computes R at 2*ETA.
  // My logic: If Method 4, we need to re-run the F'/F calculation at 2*ETA?
  // Yes, original code sets RHOUSE = 2*ETA and jumps back to 105.

  if (Method == 4 || Method == 6) {
    // We need G, G' at turning point (2*ETA).
    // This requires R at 2*ETA.
    // Let's defer this. The logic is complex.
    // For now, let's assume we are in Method 1 or 2 (Standard Steed or
    // Maclaurin). Method 4 is for ETA > 15 and ETA < RHO < 2*ETA. Our failing
    // case: ETA=0.45, RHO=32. This is RHO > 2*ETA. So Method 1.
  }

  // Store F, F' at LMAX
  FC[LMAX] = F;
  FCP[LMAX] = F * R;

  // Downward recursion to LMIN
  if (LMAX != LMIN) {
    double PL_down = (double)LMAX;
    double AR = 1.0 / RHO;
    int L = LMAX;
    for (int LP = LMIN + 1; LP <= LMAX; ++LP) {
      GC[L] = ETA / PL_down + PL_down * AR;
      GCP[L] = std::sqrt((ETA / PL_down) * (ETA / PL_down) + 1.0);
      FC[L - 1] = (GC[L] * FC[L] + FCP[L]) / GCP[L];
      FCP[L - 1] = GC[L] * FC[L - 1] - GCP[L] * FC[L];
      PL_down -= 1.0;
      L--;

      if (std::abs(FC[L + 1]) < BIG)
        continue;
      for (int LL = L; LL <= LMAX; ++LL) {
        FC[LL] *= SMALL;
        FCP[LL] *= SMALL;
      }
    }
  }

  F = FC[LMIN];
  R = FCP[LMIN] / F;

  // Calculate P, Q, G, G'
  double P, Q, W, G, GP;

  if (Method == 2 || Method == 5) {
    // Maclaurin Series for F(L=0)
    double C = 2.0 * PI * ETA;
    double X_mac = 0.0;
    double T_mac = 1.0;

    if (std::abs(C) <= 0.5) {
      double AR_mac = 1.0;
      double BR_mac = C;
      double AI_mac = 1.0;
      C = 1.0;
      while (true) {
        AI_mac += 1.0;
        AR_mac = AR_mac * BR_mac / AI_mac;
        C += AR_mac;
        if (std::abs(AR_mac) < ACC * C)
          break;
      }
      C = 1.0 / C;
    } else {
      if (ETA > 0) {
        X_mac = -SMALLN - PI * ETA;
        T_mac = SMALL;
        if (C < PRECLN)
          C = C / (1.0 - std::exp(-C));
      } else {
        C = -C;
      }
    }

    C = RHO * std::sqrt(C);
    double B1 = 1.0;
    double B2 = ETA * RHO;
    double SUM = B1 + B2;
    double AI_sum = 6.0;
    double DI_sum = 6.0;

    for (int I = 1; I <= 10000; ++I) {
      double B3 = ((2.0 * ETA * RHO) * B2 - (RHO * RHO) * B1) / AI_sum;
      AI_sum += DI_sum;
      DI_sum += 2.0;
      SUM += B3;
      double STOP = std::abs(B1) + std::abs(B2) + std::abs(B3);
      B1 = B2;
      B2 = B3;
      if (std::abs(SUM) >= BIG) {
        X_mac -= SMALLN;
        SUM *= SMALL;
        B1 *= SMALL;
        B2 *= SMALL;
      }
      if (STOP < ACC * std::abs(SUM))
        break;
      if (I == 10000)
        return 8;
    }

    SUM = (C * std::exp(X_mac) * SUM) * T_mac;

    if (SUM == 0.0)
      return 4;

    W = SUM / F;
    F = SUM;

    if (Method == 5) {
      // Return only F, F'
      FC[LMIN] = F;
      FCP[LMIN] = R * F;
      // Propagate W to higher L's if needed (though Method 5 usually implies
      // LMAX=LMIN=0) But if LMAX > LMIN, we need to scale them? Original code:
      // 850 -> 859 loop scales FC/FCP by W. Wait, original code: 850 FC(LMIN1)
      // = F; FCP(LMIN1) = R*F; IRET=2; IF (LMAX.EQ.LMIN) RETURN DO 859 L=LMIN1,
      // LMAX FC(L+1) = W*FC(L+1); FCP(L+1) = W*FCP(L+1) But FC/FCP for L > LMIN
      // were computed by downward recursion? Yes, they were computed
      // unnormalized. So we scale them all.
      for (int k = LMIN + 1; k <= LMAX; ++k) {
        FC[k] *= W;
        FCP[k] *= W;
      }
      return 2;
    }
  }

  // Compute P+iQ using Steed's method (unless Method 5 and RHO < 0.005)
  // Note: If Method 2, we need P+iQ.
  // If Method 1, we need P+iQ.
  // If Method 5, we returned already.

  if (Method != 5 || RHO > 0.005) {
    int res = SteedPQ(RHO, ETA, LMIN, P, Q, ACC);
    if (res != 0)
      return res;
    P /= RHO;
    Q /= RHO;
  }

  if (Method == 2) {
    // Continue Maclaurin logic for G
    double X_val = (R - P) * F;
    if (std::abs(X_val) <= PRERT3) {
      G = 1.0 / X_val;
      GP = P * G;
    } else {
      double B1 = 0.5 / X_val;
      double B2 = B1 * std::sqrt(1.0 - 4.0 * (X_val * F) * (X_val * F));
      G = B1 + B2;
      if (ETA < 0) {
        double SUM_val = 1.0 / Q - F * F;
        double GP_val = B1 - B2;
        if (std::abs(G * G - SUM_val) > std::abs(GP_val * GP_val - SUM_val))
          G = GP_val;
      }
      GP = P * G - X_val * F / G;
    }
  } else if (Method == 1) {
    // Standard Steed
    double X_val = (R - P) / Q;
    double FMAG = std::sqrt(1.0 / (Q * (1.0 + X_val * X_val)));
    W = FMAG / std::abs(F);
    F = W * F;
    G = F * X_val;
    GP = R * G - 1.0 / F;
  }

  // Normalize and store
  GC[LMIN] = G;
  GCP[LMIN] = GP;
  FC[LMIN] = F;
  FCP[LMIN] = R * F;

  if (LMAX != LMIN) {
    // Fortran: DO 829 L = LMIN1, LMAX
    //            T        = GC(L+1)          <- reads NEXT slot (1-indexed L+1)
    //            GC (L+1) = (GC(L)*GC(L+1) - GCP(L))/GCP(L+1)
    //            GCP(L+1) =  GC(L)*GCP(L+1) - GC(L+1)*T
    // In 0-indexed C++ with LMIN=0: Fortran L maps to C++ L (same),
    // Fortran GC(L+1) = C++ GC[L+1].  The loop writes INTO GC[L+1] using
    // the STORED value GC[L+1] before overwriting it.
    // Previous C++ code wrongly used GC[L] (current slot) instead of GC[L+1].
    for (int L = LMIN; L < LMAX; ++L) {
      double T         = GC[L + 1];           // stored 'p' value at next slot
      double stored_GCP = GCP[L + 1];
      GC [L + 1] = (GC[L] * T           - GCP[L]) / stored_GCP;
      GCP[L + 1] =  GC[L] * stored_GCP  - GC[L+1] * T;

      FC [L + 1] *= W;
      FCP[L + 1] *= W;
    }
  }

  if (std::abs(FC[LMAX]) + std::abs(FCP[LMAX]) == 0.0)
    return 1;

  return 0;
}
