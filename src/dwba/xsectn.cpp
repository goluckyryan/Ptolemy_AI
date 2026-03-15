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

  std::printf("[XSectn] Generic QNs: JA=%d JB=%d JBT=%d JBP=%d JT=%d\n",
              JA, JB, JBT, JBP, JT);

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
    std::cout << "\n[SFROMI] No-SO path: using S directly (JP_cons=" << JBP << "/2, Lx2=" << LxMax_bs << ")\n";
    for (const auto &elem : TransferSMatrix) {
      int Lx = elem.Lx;
      int Li = elem.Li;
      int Lo = elem.Lo;
      if ((Li + Lo + Lx) % 2 != 0) continue;
      if (Lx != LxMax_bs) continue;
      S_koffs[{JBP, LxMax_bs, Lo, Li}] += elem.S;
    }
    std::cout << "  KOFFS entries: " << S_koffs.size() << "\n";
  } else {
    std::cout << "\n[SFROMI] Applying double 9-J (SOSWT=TRUE) coupling...\n";
    int n9j_nonzero = 0;

    // Debug: print first few TSM entries
    {
      int cnt = 0;
      for (const auto &elem : TransferSMatrix) {
        if ((elem.Li + elem.Lo + elem.Lx) % 2 != 0) continue;
        std::cout << "TSM: Li=" << elem.Li << " Lo=" << elem.Lo
                  << " Lx=" << elem.Lx
                  << " JPI=" << elem.JPI << " (J=" << elem.JPI/2.0 << ")"
                  << " JPO=" << elem.JPO << " (J=" << elem.JPO/2.0 << ")"
                  << " |S|=" << std::abs(elem.S) << std::endl;
        if (++cnt >= 20) { std::cout << "  ...(showing first 20 TSM)\n"; break; }
      }
    }

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

      // Debug print for representative entries
      if (Li==2 && Lo==2 && Lx==2 && std::abs(elem.S) > 1e-10) {
        std::cout << "NineJ args: " << JBT/2.0 << " " << Lx << " " << JBP/2.0
                  << " " << JPI/2.0 << " " << Li << " " << JA/2.0
                  << " " << JPO/2.0 << " " << Lo << " " << JB/2.0
                  << " = " << w9j1 << std::endl;
      }
      if (std::abs(w9j1) < 1e-14) continue;

      double stat1 = std::sqrt(double((JPI+1)*(JPO+1)*(2*Lx+1)*(JBP+1)));
      double SAV9J = stat1 * w9j1;
      if (std::abs(SAV9J) < 1e-14) continue;
      n9j_nonzero++;

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

          // Debug accumulation for specific key
          if (Li==2 && Lo==2 && Lx==2 && JP_cons==1 && Lx2==2) {
            std::printf("  [ACCUM] Li=%d Lo=%d Lx=%d JPI=%d JPO=%d JP=%d Lx2=%d SAV9J=%.4f TEMP=%.4f |elemS|=%.4e\n",
              Li, Lo, Lx, JPI, JPO, JP_cons, Lx2, SAV9J, TEMP, std::abs(elem.S));
          }
        }
      }
    }
    std::cout << "  9-J non-zero terms: " << n9j_nonzero
              << ", KOFFS entries: " << S_koffs.size() << "\n";

    // Debug: print S_koffs for key entry
    {
      auto key = std::make_tuple(1, 2, 2, 2); // JP=1, Lx2=2, Lo=2, Li=2
      auto it = S_koffs.find(key);
      if (it != S_koffs.end()) {
        std::printf("  [DEBUG] S_koffs[JP=1,Lx2=2,Lo=2,Li=2] = %.5e + %.5e*i, |S|=%.5e, arg=%.3f\n",
          it->second.real(), it->second.imag(), std::abs(it->second), std::arg(it->second));
      } else {
        std::printf("  [DEBUG] S_koffs[JP=1,Lx2=2,Lo=2,Li=2] = NOT FOUND\n");
      }
    }
  } // end SOSWT block

  // -----------------------------------------------
  // S_koffs vs Ptolemy comparison table
  // -----------------------------------------------
  static const double ptol_mag[31][3] = {
    // LDEL:     -2          0          +2
    /* LI=0 */  {0.0,       0.0,       0.031541},
    /* LI=1 */  {0.0,       0.028498,  0.021680},
    /* LI=2 */  {0.020288,  0.019032,  0.017769},
    /* LI=3 */  {0.028015,  0.017730,  0.0062396},
    /* LI=4 */  {0.026046,  0.017251,  0.015631},
    /* LI=5 */  {0.027465,  0.010069,  0.011666},
    /* LI=6 */  {0.030268,  0.025254,  0.0071561},
    /* LI=7 */  {0.027662,  0.021121,  0.0052177},
    /* LI=8 */  {0.047739,  0.015263,  0.0034708},
    /* LI=9 */  {0.029899,  0.0090015, 0.0021234},
    /* LI=10 */ {0.016492,  0.0048676, 0.0011820},
    /* LI=11 */ {0.0090187, 0.0026025, 0.00064164},
    /* LI=12 */ {0.0049538, 0.0013955, 0.00034699},
    /* LI=13 */ {0.0027331, 0.00075197,0.00018802},
    /* LI=14 */ {0.0015128, 0.00040714,0.00010228},
    /* LI=15 */ {0.00083940,0.00022140,0.000055873},
    {0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},
  };
  static const double ptol_phase[31][3] = {
    /* LI=0 */  {0.0,    0.0,    -2.236},
    /* LI=1 */  {0.0,    0.325,  -3.835},
    /* LI=2 */  {3.538, -0.455,  -4.403},
    /* LI=3 */  {2.172, -2.011,  -2.169},
    /* LI=4 */  {1.265, -2.645,  -2.182},
    /* LI=5 */  {-0.680,-0.754,  -2.110},
    /* LI=6 */  {-1.574,-1.007,  -1.778},
    /* LI=7 */  {-0.262,-1.198,  -1.544},
    /* LI=8 */  {-0.977,-1.316,  -1.468},
    /* LI=9 */  {-1.365,-1.463,  -1.505},
    /* LI=10 */ {-1.499,-1.529,  -1.542},
    /* LI=11 */ {-1.545,-1.555,  -1.559},
    /* LI=12 */ {-1.562,-1.565,  -1.566},
    /* LI=13 */ {-1.568,-1.568,  -1.569},
    /* LI=14 */ {-1.570,-1.570,  -1.570},
    /* LI=15 */ {-1.570,-1.570,  -1.570},
    {0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},{0},
  };
  std::printf("\n=== S_koffs vs Ptolemy (JP=1/2, Lx2=2) comparison ===\n");
  std::printf("  LI  LDEL  |S|_cpp    |S|_ptol   ratio   arg_cpp  arg_ptol  darg\n");
  for (auto &[k, v] : S_koffs) {
    auto [jp, lx2, lo, li] = k;
    if (jp != 1 || lx2 != 2) continue;
    int LDEL = lo - li;
    if (LDEL < -2 || LDEL > 2 || (LDEL % 2 != 0)) continue;
    int ldel_idx = (LDEL + 2) / 2; // 0,1,2 for LDEL=-2,0,+2
    if (li < 0 || li > 15) continue;
    double mag_ptol  = ptol_mag[li][ldel_idx];
    double phase_ptol= ptol_phase[li][ldel_idx];
    double mag_cpp   = std::abs(v);
    double arg_cpp   = std::arg(v);
    double ratio     = (mag_ptol > 1e-10) ? mag_cpp/mag_ptol : 0.0;
    double darg      = (mag_ptol > 1e-10) ? arg_cpp - phase_ptol : 0.0;
    if (mag_ptol > 1e-10) {
      std::printf("  Li=%2d LDEL=%+d |S|_cpp=%.5e |S|_ptol=%.5e ratio=%.3f arg_cpp=%+.3f arg_ptol=%+.3f darg=%+.3f\n",
                  li, LDEL, mag_cpp, mag_ptol, ratio, arg_cpp, phase_ptol, darg);
    }
  }
  {
    int cnt_diag = 0, cnt_offdiag = 0;
    for (auto &[k, v] : S_koffs) {
      auto [jp, lx2, lo, li] = k;
      if (li == lo) cnt_diag++; else cnt_offdiag++;
    }
    std::printf("  Total S_koffs=%zu (%d diag, %d offdiag)\n",
                S_koffs.size(), cnt_diag, cnt_offdiag);
  }

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

  // Determine IODD from primary (full) LDELMN
  //   LBP_bf = lP (actual orbital L of projectile constituent)
  //   LBT_bf = lT (actual orbital L of target constituent)
  {
    int LBP_bf = lP;
    int LBT_bf = lT;
    int NDEL_bf = LxMax_bs + 1 - ((LxMax_bs + LBP_bf + LBT_bf) % 2);
    int LDELMN_bf = 1 - NDEL_bf;
    std::printf("  [BETCAL] LBP=%d LBT=%d NDEL_primary=%d LDELMN_primary=%d\n",
                LBP_bf, LBT_bf, NDEL_bf, LDELMN_bf);
    [[maybe_unused]] int IODD_check = std::abs(LDELMN_bf) % 2;
  }

  // Compute LDELMN_primary dynamically from quantum numbers
  int NDEL_primary = LxMax_bs + 1 - ((LxMax_bs + lP + lT) % 2);
  int LDELMN_primary = 1 - NDEL_primary;
  int IODD_global = std::abs(LDELMN_primary) % 2;
  std::printf("  [BETCAL] NDEL_primary=%d LDELMN_primary=%d IODD=%d\n",
              NDEL_primary, LDELMN_primary, IODD_global);
  double FACTOR_BET = (IODD_global != 0) ? (-0.5 / ka) : (0.5 / ka);  // = +0.5/ka
  std::printf("  [BETCAL] FACTOR_BET=%.4e (IODD=%d, ka=%.4f)\n",
              FACTOR_BET, IODD_global, ka);

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

    // Debug print for off-diagonal Lx2=2 entries
    if (Lx2 == 2 && Lo != Li) {
      std::printf("  [BETCAL DBG] JP=%d/2 Lx2=%d Lo=%d Li=%d LDEL=%d NDEL=%d LDELMN=%d range=[%d,%d]\n",
                  JP_cons, Lx2, Lo, Li, LDEL, NDEL_lx, LDELMN_lx,
                  LDELMN_lx, LDELMN_lx + 2*(NDEL_lx-1));
    }

    // Range check
    if (LDEL < LDELMN_lx || LDEL > LDELMN_lx + 2*(NDEL_lx-1)) {
      if (Lx2 == 2 && Lo != Li)
        std::printf("    -> FILTERED OUT by LDEL range\n");
      continue;
    }
    if ((LDEL - LDELMN_lx) % 2 != 0) {
      if (Lx2 == 2 && Lo != Li)
        std::printf("    -> FILTERED OUT by step parity\n");
      continue;
    }

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

      if (Lx2 == 2) {
        std::printf("    BETAS[JP=%d/2,Lx2=%d,Mx=%d,Lo=%d] += TEMPS=%.4e * |S|=%.4e CG=%.4e"
                    "  (Li=%d,LDEL=%d,Mx_loop=%d)\n",
                    JP_cons, Lx2, official_Mx, Lo, TEMPS,
                    std::abs(S_phased), cg, Li, LDEL, Mx_loop);
      }
    }
  }

  std::cout << "\nBETAS entries: " << BETAS.size() << "\n";
  {
    int cnt = 0;
    for (auto &[k, b] : BETAS) {
      auto [jp, lx2, mx, lo] = k;
      std::printf("  BETAS(JP=%d/2,Lx2=%d,Mx=%d,Lo=%d)=(%.4e,%.4e) |b|=%.4e\n",
                  jp, lx2, mx, lo, b.real(), b.imag(), std::abs(b));
      if (++cnt >= 12) { std::printf("  ... (%zu total)\n", BETAS.size()); break; }
    }
  }

  // -----------------------------------------------
  // BETCAL (JP,Lx2,Mx) breakdown and spot-check dumps
  // -----------------------------------------------
  {
    std::map<std::tuple<int,int,int>, int> mx_count;
    for (auto &[k, b] : BETAS) {
      auto [jp, lx2, mx, lo] = k;
      mx_count[{jp, lx2, mx}]++;
    }
    std::printf("\nBETAS (JP,Lx2,Mx) breakdown:\n");
    for (auto &[k, n] : mx_count) {
      auto [jp, lx2, mx] = k;
      std::printf("  JP=%d/2 Lx2=%d Mx=%d : %d entries\n", jp, lx2, mx, n);
    }

    std::printf("\nBETAS JP=1/2 Lx2=2 Mx=0 (for comparison with Ptolemy):\n");
    for (auto &[k, b] : BETAS) {
      auto [jp, lx2, mx, lo] = k;
      if (jp==1 && lx2==2 && mx==0)
        std::printf("  Lo=%-2d : %+.5e %+.5ei\n", lo, b.real(), b.imag());
    }
    std::printf("\nBETAS JP=1/2 Lx2=2 Mx=1:\n");
    for (auto &[k, b] : BETAS) {
      auto [jp, lx2, mx, lo] = k;
      if (jp==1 && lx2==2 && mx==1)
        std::printf("  Lo=%-2d : %+.5e %+.5ei\n", lo, b.real(), b.imag());
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
  std::cout << "\nAngle (deg)    dSigma/dOmega (mb/sr)   [DWBA C++ 9J]\n";

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

    // Debug spot-checks vs Ptolemy reference values
    if (theta_deg == 0.0) {
      std::cout << "[DEBUG 0deg] dSigma_9J=" << dSigma << " mb/sr\n";
      std::cout << "[DEBUG 0deg] Ptolemy ref = 1.863 mb/sr\n";
      std::cout << "[DEBUG 0deg] Ratio C++/Ptol = " << dSigma/1.863 << "\n";
    }
    if (theta_deg == 15.0)
      std::cout << "[DEBUG 15deg] dSigma_9J=" << dSigma << " mb/sr  (Ptol: 2.535)\n";
    if (theta_deg == 30.0)
      std::cout << "[DEBUG 30deg] dSigma_9J=" << dSigma << " mb/sr  (Ptol: 0.905)\n";

    std::cout << std::fixed << std::setprecision(1) << theta_deg
              << "           " << std::scientific << std::setprecision(4)
              << dSigma << "\n";
  }
}
