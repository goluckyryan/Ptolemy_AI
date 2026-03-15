// xsectn.cpp — DWBA::XSectn() extracted from dwba.cpp
// Implements BETCAL + AMPCAL + cross-section output (Ptolemy XSECTN port)
//
// DO NOT MODIFY dwba.cpp — this file is the canonical home for XSectn().

#include "dwba.h"
#include "math_utils.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <tuple>

void DWBA::XSectn() {
  // ============================================================
  // DWBA Cross Section — Ptolemy SFROMI (double 9-J) + BETCAL + AMPCAL + XSECTN
  //
  // Quantum numbers (doubled, integers) for 33Si(d,p)34Si:
  //   JA  = 2*j_deuteron = 2   (spin of incoming projectile a)
  //   JB  = 2*j_proton   = 1   (spin of outgoing ejectile b)
  //   JBT = 2*j_target_neutron = 3   (j of neutron orbit in 33Si, 0d3/2)
  //   JBP = 2*j_proj_neutron   = 1   (j of neutron orbit in deuteron, 0s1/2)
  //   JT  = 2*J_target   = 3   (33Si ground state J=3/2)
  //   JB_res = 2*J_residual = 0 (34Si ground state J=0)
  //
  // SFROMI (Ptolemy source.mor ~29150):
  //   For each radial integral I(Li, Lo, Lx):
  //     Loop JPI = |2*Li - JA| .. 2*Li+JA  step 2  (incoming J = L +/- 1/2 for spin-1 deuteron)
  //     Loop JPO = |2*Lo - JB| .. 2*Lo+JB  step 2  (outgoing J = L +/- 1/2 for proton)
  //       SAV9J = sqrt((JPI+1)(JPO+1)(2Lx+1)(JBP+1))
  //               * W9J( JBT,  2*Lx, JBP,
  //                      JPI,  2*Li, JA,
  //                      JPO,  2*Lo, JB )
  //       if SAV9J == 0: skip
  //       Loop JP_cons = JPBASE .. JPMX  step 2   [conserved total J]
  //         Loop Lx2 = triangle(Li,Lo) ∩ [|JBT-JP_cons|/2 .. (JBT+JP_cons)/2]
  //           TEMP = sqrt((JPI+1)(JPO+1)(2*Lx2+1)(JP_cons+1))
  //                  * W9J( JBT,  2*Lx2, JP_cons,
  //                         JPI,  2*Li,  JA,
  //                         JPO,  2*Lo,  JB )
  //           if LX2==LXP && JP==JBP && LASI==LI && LASO==LO: TEMP=SAV9J (skip 2nd 9J)
  //           if (Lx2+Lx) odd: TEMP = -TEMP
  //           S(JP_cons, JT, Lx2, LDEL=Lo-Li) += SAV9J * TEMP * I_phased
  //
  // BETCAL (Ptolemy ~3480):
  //   For each (JP_cons, JT, Lx2, LDEL) KOFFS entry:
  //     Li = Lo - LDEL  (LDEL = Lo - Li)
  //     MXZ_second = Lx2 + LDEL   (NOT parity-based!)
  //     For Mx = max(0,MXZ_second) .. Lx2:
  //       cg = CG(Li, 0; Lx2, Mx | Lo, Mx)
  //       BETAS(JP_cons, Lx2, Mx, Lo) += FACTOR_BET*(2Li+1)*cg * S_phased * e^i(sigma_Li+sigma_Lo)
  //
  // AMPCAL + XSECTN:
  //   F(JP_cons, Lx2, Mx, theta) = sum_Lo BETAS(...,Lo) * sf(Lo,Mx) * P_Lo^Mx(cos_theta)
  //   dSigma = 10 * sum_{JP_cons,Lx2,Mx} FMNEG(Mx) * |F|^2

  double ka    = Incoming.k;
  double kb    = Outgoing.k;
  (void)kb;  // used in cross-section formula via kinematic prefactor (absorbed in FACTOR_BET)

  // -----------------------------------------------
  // Reaction-generic doubled quantum numbers (read from DWBA object)
  //   JA  = 2*spin(projectile)    e.g. deuteron=2, proton=1
  //   JB  = 2*spin(ejectile)      e.g. proton=1, deuteron=2
  //   JBT = 2*j(neutron in target bound state)
  //   JBP = 2*j(neutron in projectile bound state)
  //   JT  = 2*J(target nucleus)
  // -----------------------------------------------
  const int JA  = Incoming.JSPS;                             // 2*j_projectile
  const int JB  = Outgoing.JSPS;                             // 2*j_ejectile
  const int JBT = (int)std::round(2.0 * TargetBS.j);        // 2*j_neutron_in_target
  const int JBP = (int)std::round(2.0 * ProjectileBS.j);    // 2*j_neutron_in_proj
  const int JT  = (int)std::round(2.0 * SpinTarget);        // 2*J_target_nucleus
  (void)JT;

  // JP_conserved (2*J_total) runs from |JA-JB| to JA+JB (Ptolemy ANGSET line ~28894)
  // JPBASE = |JA - JB| = |2 - 1| = 1
  // JPMX   = JA + JB   = 2 + 1   = 3
  const int JPBASE = std::abs(JA - JB);  // = 1
  const int JPMX   = JA + JB;            // = 3

  // -----------------------------------------------
  // SFROMI: build S_koffs(JP_cons, Lx2, Lo, Li) via double 9-J
  // -----------------------------------------------
  using KOFFS_key = std::tuple<int,int,int,int>; // (JP_cons, Lx2, Lo, Li)
  std::map<KOFFS_key, std::complex<double>> S_koffs;

  int lT = TargetBS.l;
  int lP = ProjectileBS.l;
  int LxMax_bs = lT + lP;

  // Determine whether spin-orbit coupling is active
  bool SOSWT = (Incoming.Pot.VSO != 0.0 || Outgoing.Pot.VSO != 0.0);

  if (!SOSWT) {
    // Simple case: no spin-orbit; S stored directly without 9-J coupling
    for (const auto &elem : TransferSMatrix) {
      int Lx = elem.Lx;
      int Li = elem.Li;
      int Lo = elem.Lo;
      if ((Li + Lo + Lx) % 2 != 0) continue;
      if (Lx != LxMax_bs) continue;
      S_koffs[{JBP, LxMax_bs, Lo, Li}] += elem.S;
    }
  } else {
    for (const auto &elem : TransferSMatrix) {
      int Lx  = elem.Lx;
      int Li  = elem.Li;
      int Lo  = elem.Lo;
      int JPI = elem.JPI;   // 2*J_incoming (already J-split in InelDc)
      int JPO = elem.JPO;   // 2*J_outgoing

      // Parity filter: (Li + Lo + Lx) must be even
      if ((Li + Lo + Lx) % 2 != 0) continue;

      // First 9-J: W9J(JBT, 2*Lx, JBP, JPI, 2*Li, JA, JPO, 2*Lo, JB)
      double w9j1 = NineJ(JBT/2.0, (double)Lx,  JBP/2.0,
                          JPI/2.0, (double)Li,   JA/2.0,
                          JPO/2.0, (double)Lo,   JB/2.0);

      if (std::abs(w9j1) < 1e-14) continue;

      double stat1 = std::sqrt(double((JPI+1)*(JPO+1)*(2*Lx+1)*(JBP+1)));
      double SAV9J = stat1 * w9j1;
      if (std::abs(SAV9J) < 1e-14) continue;

      // Second 9-J loop: over conserved total J (JP_cons) and transfer multipolarity (Lx2)
      for (int JP_cons = JPBASE; JP_cons <= JPMX; JP_cons += 2) {
        // Triangle constraints for Lx2:
        //   from row (JBT, 2*Lx2, JP_cons): |JBT-JP_cons|/2 <= Lx2 <= (JBT+JP_cons)/2
        //   from triangle(Li, Lo, Lx2):      |Li-Lo| <= Lx2 <= Li+Lo
        int Lx2_mn_jbt = std::abs(JBT - JP_cons) / 2;
        int Lx2_mx_jbt = (JBT + JP_cons) / 2;
        int Lx2_mn_lo  = std::abs(Li - Lo);
        int Lx2_min = std::max(Lx2_mn_jbt, Lx2_mn_lo);
        int Lx2_max = std::min(Lx2_mx_jbt, Li + Lo);
        if (Lx2_min > Lx2_max) continue;

        for (int Lx2 = Lx2_min; Lx2 <= Lx2_max; Lx2++) {
          // ANGSET parity restriction (Ptolemy SETSPT formula):
          //   NDEL = Lx2 + 1 - MOD(Lx2, 2)
          int NDEL_ang = Lx2 + 1 - (Lx2 % 2);
          if (NDEL_ang <= 0) continue;

          double TEMP;
          // Diagonal shortcut: if second 9-J matches first exactly
          if (Lx2 == Lx && JP_cons == JBP) {
            TEMP = SAV9J;
          } else {
            // Second 9-J: W9J(JBT, 2*Lx2, JP_cons, JPI, 2*Li, JA, JPO, 2*Lo, JB)
            double w9j2 = NineJ(JBT/2.0, (double)Lx2, JP_cons/2.0,
                                JPI/2.0, (double)Li,   JA/2.0,
                                JPO/2.0, (double)Lo,   JB/2.0);
            double stat2 = std::sqrt(double((JPI+1)*(JPO+1)*(2*Lx2+1)*(JP_cons+1)));
            TEMP = stat2 * w9j2;
          }
          if (std::abs(TEMP) < 1e-14) continue;

          // Phase factor: (-1)^(Lx + Lx2)
          if ((Lx + Lx2) % 2 != 0) TEMP = -TEMP;

          // Accumulate S_koffs
          auto contrib = SAV9J * TEMP * elem.S;
          S_koffs[{JP_cons, Lx2, Lo, Li}] += contrib;
        }
      }
    }
  } // end SOSWT block

  // -----------------------------------------------
  // Coulomb phases sigma_L  (cumulative Coulomb phase shifts)
  // -----------------------------------------------
  int Lmax_col = 25;
  std::vector<double> sig_a(Lmax_col + 1, 0.0);
  std::vector<double> sig_b(Lmax_col + 1, 0.0);
  for (int L = 1; L <= Lmax_col; L++) {
    sig_a[L] = sig_a[L-1] + std::atan2(Incoming.eta, (double)L);
    sig_b[L] = sig_b[L-1] + std::atan2(Outgoing.eta, (double)L);
  }
  auto get_sig_a = [&](int L) { return sig_a[std::min(L, Lmax_col)]; };
  auto get_sig_b = [&](int L) { return sig_b[std::min(L, Lmax_col)]; };

  // -----------------------------------------------
  // BETCAL — build BETAS array
  //
  // For each (JP_cons, Lx2, Lo, Li) in S_koffs:
  //   LDEL = Lo - Li
  //   NDEL = Lx2 + 1 - MOD(Lx2 + LBP + LBT, 2)   [Ptolemy SETSPT]
  //   LDELMN = 1 - NDEL
  //   MXZ_loop = Lx2 + LDEL   (starting Mx_loop)
  //   For Mx_loop = max(0, MXZ_loop) .. Lx2:
  //     official_Mx = Mx_loop   (direct from Ptolemy AMPCAL MX formula)
  //     cg = CG(Li, 0; Lx2, Mx_loop | Lo, Mx_loop)
  //     BETAS[JP_cons, Lx2, official_Mx, Lo] += FACTOR_BET*(2Li+1)*cg * e^i(sig_a(Li)+sig_b(Lo)) * S_val
  //
  // FACTOR_BET = 0.5/ka  (positive, since LDELMN=-2 → IODD=0 for 33Si case)
  // -----------------------------------------------

  // Compute LDELMN_primary dynamically from quantum numbers
  int NDEL_primary = LxMax_bs + 1 - ((LxMax_bs + lP + lT) % 2);
  int LDELMN_primary = 1 - NDEL_primary;
  int IODD_global = std::abs(LDELMN_primary) % 2;
  double FACTOR_BET = (IODD_global != 0) ? (-0.5 / ka) : (0.5 / ka);  // = +0.5/ka

  using BETA_key = std::tuple<int,int,int,int>; // (JP_cons, Lx2, official_Mx, Lo)
  std::map<BETA_key, std::complex<double>> BETAS;

  for (auto &[k, S_val] : S_koffs) {
    auto [JP_cons, Lx2, Lo, Li] = k;
    int LDEL = Lo - Li;

    // Per-Lx2 NDEL/LDELMN from Ptolemy SETSPT formula:
    //   NDEL = Lx2 + 1 - MOD(Lx2 + LBP + LBT, 2)
    //   LDELMN = 1 - NDEL
    int LBP_bc = lP;
    int LBT_bc = lT;
    int NDEL_lx   = Lx2 + 1 - ((Lx2 + LBP_bc + LBT_bc) % 2);
    int LDELMN_lx = 1 - NDEL_lx;

    if (NDEL_lx <= 0) continue;

    // Range check
    if (LDEL < LDELMN_lx || LDEL > LDELMN_lx + 2*(NDEL_lx-1)) continue;
    if ((LDEL - LDELMN_lx) % 2 != 0) continue;

    // Coulomb phase factor: e^i(sigma_Li + sigma_Lo)
    double phase_angle = get_sig_a(Li) + get_sig_b(Lo);
    std::complex<double> phase_fac(std::cos(phase_angle), std::sin(phase_angle));
    std::complex<double> S_phased = phase_fac * S_val;

    // BETCAL Mx_loop range: [max(0, Lx2+LDEL), Lx2]
    int MXZ_loop = Lx2 + LDEL;

    for (int Mx_loop = std::max(0, MXZ_loop); Mx_loop <= Lx2; Mx_loop++) {
      if (Lo < Mx_loop) continue;  // CG vanishes: m cannot exceed L

      // official_Mx = Mx_loop (see Ptolemy AMPCAL line ~389 derivation in header)
      int official_Mx = Mx_loop;
      if (official_Mx < 0 || official_Mx > Lx2) continue;

      // Clebsch-Gordan coefficient: CG(Li, 0; Lx2, Mx_loop | Lo, Mx_loop)
      double cg = ClebschGordan((double)Li,  0.0,
                                (double)Lx2, (double)Mx_loop,
                                (double)Lo,  (double)Mx_loop);
      if (!std::isfinite(cg)) continue;
      if (std::abs(cg) < 1e-14) continue;

      double TEMPS = FACTOR_BET * (2.0*Li + 1.0) * cg;
      BETAS[{JP_cons, Lx2, official_Mx, Lo}] += TEMPS * S_phased;
    }
  }

  // -----------------------------------------------
  // AMPCAL + XSECTN — angular distribution
  //
  // F(JP_cons, Lx2, Mx, theta) = sum_Lo BETAS[JP_cons,Lx2,Mx,Lo]
  //                               * sf(Lo,Mx) * P_Lo^Mx(cos_theta)
  //
  // sf(Lo, Mx) = prod_{n=1}^{Mx} 1/sqrt((Lo+n)(Lo-n+1))
  //   [sqrt_factorial correction for Mx>0 amplitudes; see AGENT_FINDINGS.md]
  //
  // PLM: Ptolemy PLMSUB uses Condon-Shortley phase convention.
  //   std::assoc_legendre(l, m, x) also includes (-1)^m — same convention.
  //   std::legendre(l, x) = P_l(x) for m=0 (no phase).
  //
  // dSigma/dOmega = 10 * sum_{JP_cons,Lx2,Mx} FMNEG(Mx) * |F(theta)|^2
  //   Factor 10 converts fm^2 -> mb (1 fm^2 = 10 mb).
  //   FMNEG = 1 for Mx=0; FMNEG = 2 for Mx>0 (time-reversal degeneracy).
  // -----------------------------------------------
  std::cout << "\nAngle (deg)    dSigma/dOmega (mb/sr)\n";

  for (double theta_deg = AngleMin; theta_deg <= AngleMax + 1e-9; theta_deg += AngleStep) {
    double theta_rad  = theta_deg * M_PI / 180.0;
    double cos_theta  = std::cos(theta_rad);

    double dSigma = 0.0;

    // Collect unique (JP_cons, Lx2, Mx) KOFFS slots from BETAS
    std::set<std::tuple<int,int,int>> koffs_set;
    for (auto &[k, b] : BETAS) {
      auto [jp, lx2, mx, lo] = k;
      koffs_set.insert({jp, lx2, mx});
    }

    for (auto &[jp, lx2, mx] : koffs_set) {
      double FMNEG = (mx == 0) ? 1.0 : 2.0;

      std::complex<double> F_amp(0.0, 0.0);

      for (auto &[k2, b] : BETAS) {
        auto [jp2, lx3, mx2, lo] = k2;
        if (jp2 != jp || lx3 != lx2 || mx2 != mx) continue;
        if (lo < mx) continue;  // guard against invalid PLM

        // sf: sqrt_factorial correction (Ptolemy BETCAL second-pass / AMPCAL correction)
        // sf = prod_{n=1}^{mx} 1 / sqrt( (lo+n)*(lo-n+1) )
        double sf = 1.0;
        for (int n = 1; n <= mx; n++) {
          double denom = double(lo + n) * double(lo - n + 1);
          if (denom <= 0.0) { sf = 0.0; break; }
          sf /= std::sqrt(denom);
        }
        if (sf == 0.0) continue;

        // Associated Legendre polynomial P_lo^mx(cos_theta)
        // Both std::legendre (m=0) and std::assoc_legendre (m>0) match Ptolemy PLMSUB
        double PLo_Mx;
        if (mx == 0) {
          PLo_Mx = std::legendre(lo, cos_theta);
        } else if (mx > lo) {
          continue;
        } else {
          // std::assoc_legendre includes (-1)^m Condon-Shortley phase — same as Ptolemy
          PLo_Mx = std::assoc_legendre(lo, mx, cos_theta);
        }

        auto contrib = b * sf * PLo_Mx;
        if (std::isfinite(contrib.real()) && std::isfinite(contrib.imag()))
          F_amp += contrib;
      }

      // FMNEG * |F|^2 contribution; factor 10 applied outside inner loop
      dSigma += 10.0 * FMNEG * std::norm(F_amp);
    }

    std::cout << std::fixed << std::setprecision(1) << theta_deg
              << "           " << std::scientific << std::setprecision(4)
              << dSigma << "\n";
  }
}
