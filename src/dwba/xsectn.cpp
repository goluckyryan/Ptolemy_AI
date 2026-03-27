// xsectn.cpp — DWBA::XSectn()
// Faithful translation of Ptolemy SFROMI + BETCAL + AMPCAL + XSECTN
//
// Flow:
//   InelDcFaithful2() → TransferSMatrix entries (I_accum per Li,Lo,JPI,JPO,Lx)
//   XSectn():
//     Step 1: ANGSET-equivalent — build JTOCS and allocate SMAG/SPHASE tables
//     Step 2: SFROMI — for each (Li,Lo,JPI,JPO,Lx) I_accum, apply spectroscopic
//             factor (ATERM), phase, then double 9-J to accumulate into SMAG(KOFFS,LiIdx)
//     Step 3: BETCAL — outer loop over Lo, inner over KOFFS/Li, build BETAS
//     Step 4: AMPCAL — for each KOFFS and angle, sum Lo contributions * PLM
//     Step 5: DCS = 10 * sum_KOFFS FMNEG(MX) * |F|^2

#include "dwba.h"
#include "math_utils.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

void DWBA::XSectn() {
  // ── Quantum numbers ──────────────────────────────────────────────────
  const int JA  = Incoming.JSPS;                          // 2*j_proj  (deuteron=2)
  const int JB  = Outgoing.JSPS;                          // 2*j_eject (proton=1)
  const int JBT = (int)std::round(2.0 * TargetBS.j);     // 2*j_neutron_in_target (5 for 17O g.s.)
  const int JBP = (int)std::round(2.0 * ProjectileBS.j); // 2*j_neutron_in_proj   (1 for d)
  const int JPBASE = std::abs(JA - JB);                  // min 2*J_cons (=1)
  const int JPMX   = JA + JB;                             // max 2*J_cons (=3)
  const int JTBASE = JBT;                                 // only one JT for transfer
  const int NJP    = (JPMX - JPBASE)/2 + 1;              // # JP values (=2: JP=1,3)
  const int NLX    = JBT/2 + 1;                           // # LX values (0..JBT/2)

  const int lT = TargetBS.l;     // neutron l in target (=2)
  const int lP = ProjectileBS.l; // neutron l in proj   (=0)
  const int LxMin_bs = std::abs(lT - lP);
  const int LxMax_bs = lT + lP;

  const double ka  = Incoming.k;

  // ── LINTRP: With TRANSW=true, LIMOST=LMAX (no extrapolation beyond computed range)
  // C++ computes all Li from LMIN..LMAX directly in InelDc.
  // We use LMIN=LminSet, LMAX=LmaxSet.
  const int LMIN_dc = (LminSet >= 0) ? LminSet : 0;
  const int LMAX_dc = (LmaxSet >= 0) ? LmaxSet : 40;
  const int LIMOST  = LMAX_dc;  // Fortran: LIMOST=LMAX for transfer with TRANSW
  const int NUMLIS  = (LMAX_dc - LMIN_dc) / 1 + 1;  // LSTEP=1

  // LOMOST = LIMOST + LXMAX (from Fortran line 11643)
  const int LOMOST = LIMOST + LxMax_bs;

  // ── Step 1: ANGSET — build JTOCS table ───────────────────────────────
  // JTOCS[KOFFS] = {LDEL=Lo-Li, LX, JP, JT}
  // Valid entries: JP in [JPBASE..JPMX], JT=JBT, LX in [LxMin_bs..LxMax_bs],
  //   LDEL s.t. parity OK and triangle(Li,Lo,LX) can be satisfied for some Li
  // Fortran ANGSET scans all (JP, JBT, LX, LDEL) combinations.
  // For transfer: LDELMX = LXMAX (= LxMax_bs = lT+lP = 2)
  // LDEL ranges from -LDELMX to +LDELMX step 2 (same parity as LX)
  // But actually: NDEL = LX + 1 - MOD(LX + lP + lT, 2)
  //               LDELMN = 1 - NDEL
  //               LDEL = LDELMN, LDELMN+2, ..., LDELMN+2*(NDEL-1)
  struct KoffsEntry {
    int LDEL;  // Lo - Li
    int LX;
    int JP;
    int JT;    // = JBT always for transfer
    int MX;    // = (LDEL+LX+1)/2
  };
  std::vector<KoffsEntry> JTOCS;  // index 0-based (Fortran 1-based)

  // Build JTOCS in same order as Fortran ANGSET (JP outer, JT, LX, LDEL inner)
  for (int JP = JPBASE; JP <= JPMX; JP += 2) {
    int JT = JBT;  // only one residual state
    // JTOCS LX range from J-conservation triangle only (NOT limited by bound state)
    // The first 9-J LXP (= lT+lP range) is separate from the JTOCS LX
    int LXmn = std::abs(JT - JP) / 2;
    int LXmx = (JT + JP) / 2;
    for (int LX = LXmn; LX <= LXmx; LX++) {
      // NDEL/LDELMN from Ptolemy SETSPT
      int NDEL   = LX + 1 - ((LX + lP + lT) % 2);
      int LDELMN = 1 - NDEL;
      if (NDEL <= 0) continue;
      for (int id = 0; id < NDEL; id++) {
        int LDEL = LDELMN + 2 * id;
        int MX   = (LDEL + LX + 1) / 2;
        if (MX < 0 || MX > LX) continue;
        JTOCS.push_back({LDEL, LX, JP, JT, MX});
      }
    }
  }
  int NSPL = (int)JTOCS.size();

  // SMAG(KOFFS, Li_idx) and SPHASE(KOFFS, Li_idx)
  // Li_idx = Li - LMIN_dc  (0-based), Li from LMIN_dc..LMAX_dc
  // Fortran: SMAG(KOFFS, I) where I = Li - LBASE + 1, LBASE = LMIN
  std::vector<std::vector<double>> SMAG  (NSPL, std::vector<double>(NUMLIS, 0.0));
  std::vector<std::vector<double>> SPHASE(NSPL, std::vector<double>(NUMLIS, 0.0));

  // Map (LDEL, LX, JP, JT) → KOFFS index
  auto koffs_index = [&](int LDEL, int LX, int JP, int JT) -> int {
    for (int k = 0; k < NSPL; k++) {
      if (JTOCS[k].LDEL==LDEL && JTOCS[k].LX==LX && JTOCS[k].JP==JP && JTOCS[k].JT==JT)
        return k;
    }
    return -1;
  };

  // ── Step 2: SFROMI ────────────────────────────────────────────────────
  // For each TransferSMatrix entry (I_accum channel: Li, Lo=LASI, JPI, JPO, LXP):
  //   Apply ATERM spectroscopic factor and phase → SR, SI
  //   Then double 9-J loop → accumulate into SMAG/SPHASE

  // ATERM(LXP) = spectroscopic amplitude for this LX
  // Fortran: ATERM = SPFACT * sqrt(2*lT+1) * CG(lT,0; lP,0 | LXP,0) * ... (from BSSET)
  // In our C++ setup: SPAMP = AV18 spectroscopic factor ≈ 0.97069
  // ATERM(LXP+1) = SPAMP * sqrt(2*LXP+1) * ... (from bound state coupling)
  // For lP=0, lT=2, LXP=2: CG(2,0;0,0|2,0)=1 → ATERM ≈ SPAMP * sqrt(5)
  // Let's use the same formula as the existing code used:
  // ATERM[lx] = SpectroscopicAmplitude for lx (precomputed from bound state)
  // We'll compute it the same way: SPAMP * sqrt(2*lx+1) * CG
  std::vector<double> ATERM(LxMax_bs + 2, 0.0);
  for (int lx = LxMin_bs; lx <= LxMax_bs; lx++) {
    double cg = ClebschGordan((double)lT, 0.0, (double)lP, 0.0, (double)lx, 0.0);
    double SPAMP_val = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
    ATERM[lx] = SPAMP_val * std::sqrt(2.0 * lx + 1.0) * cg;
  }

  bool SOSWT = (Incoming.Pot.VSO != 0.0 || Outgoing.Pot.VSO != 0.0);

  for (const auto& elem : TransferSMatrix) {
    int LXP = elem.Lx;
    int Li  = elem.Li;   // = LASI in Fortran
    int Lo  = elem.Lo;   // = LASO in Fortran
    int JPI = elem.JPI;  // 2*J_incoming
    int JPO = elem.JPO;  // 2*J_outgoing

    if (LXP < LxMin_bs || LXP > LxMax_bs) continue;
    if (Li < LMIN_dc || Li > LMAX_dc) continue;
    int Li_idx = Li - LMIN_dc;

    // Parity check: Li + Lo + LXP must be even
    if ((Li + Lo + LXP) % 2 != 0) continue;

    // SFROMI TEMP and phase (Ptolemy source lines 29110-29130):
    // TEMP = FACTOR * ATERM(LXP+1) / sqrt(2*LASI+1)
    // ITEST = LASI + LASO + 2*LXP + 1
    // if ITEST%4 >= 2: TEMP = -TEMP
    // Then: if ITEST%2 == 0: SR=TEMP*Ire, SI=TEMP*Iim (even power of i)
    //        else:            SR=-TEMP*Iim, SI=TEMP*Ire (odd power → multiply by i)
    double FACTOR_sfromi = 2.0 * std::sqrt(ka * Outgoing.k / (Incoming.Ecm * Outgoing.Ecm));
    double TEMP = FACTOR_sfromi * ATERM[LXP] / std::sqrt(2.0 * Li + 1.0);
    int ITEST = Li + Lo + 2 * LXP + 1;
    if (ITEST % 4 >= 2) TEMP = -TEMP;

    double I_re = elem.S.real();
    double I_im = elem.S.imag();
    double SR, SI;
    if (ITEST % 2 == 0) {
      SR = TEMP * I_re;
      SI = TEMP * I_im;
    } else {
      SR = -TEMP * I_im;
      SI =  TEMP * I_re;
    }

    if (!SOSWT) {
      // No spin-orbit: store directly (simplified path, no 9-J)
      int LDEL = Lo - Li;
      int k = koffs_index(LDEL, LXP, JBP, JBT);
      if (k < 0) continue;
      double mag = std::sqrt(SR*SR + SI*SI);
      double ph  = std::atan2(SI, SR);
      // Accumulate (add complex SR+iSI to existing, then recompute mag/phase)
      // Better: accumulate complex, convert at end
      // Use a temporary complex accumulator
      // (We'll do this in a second pass — for now mark as TODO)
      // For now: overwrite (assumes no repeated (Li,Lo,Lx) channels)
      if (k >= 0 && Li_idx >= 0 && Li_idx < NUMLIS) {
        SMAG  [k][Li_idx] = mag;
        SPHASE[k][Li_idx] = ph;
      }
    } else {
      // Spin-orbit: double 9-J loop
      int JPIMN = JPI, JPIMX = JPI;
      if ((Incoming.Pot.VSO == 0.0)) { JPIMN = std::abs(2*Li - JA); JPIMX = 2*Li + JA; }
      int JPOMN = JPO, JPOMX = JPO;
      if ((Outgoing.Pot.VSO == 0.0)) { JPOMN = std::abs(2*Lo - JB); JPOMX = 2*Lo + JB; }

      for (int jpi = JPIMN; jpi <= JPIMX; jpi += 2) {
        for (int jpo = JPOMN; jpo <= JPOMX; jpo += 2) {
          // First 9-J
          double SAV9J = std::sqrt(double((jpi+1)*(jpo+1)*(2*LXP+1)*(JBP+1)))
            * NineJ(JBT/2.0, (double)LXP, JBP/2.0,
                    jpi/2.0, (double)Li,  JA/2.0,
                    jpo/2.0, (double)Lo,  JB/2.0);
          if (std::abs(SAV9J) < 1e-14) continue;
          double TEMPR = SAV9J * SR;
          double TEMPI = SAV9J * SI;

          // Loop over conserved J (JP) and LX
          for (int JP = JPBASE; JP <= JPMX; JP += 2) {
            int LXmn2 = std::max({std::abs(JBT-JP)/2, std::abs(Lo-Li)});
            int LXmx2 = std::min((JBT+JP)/2, Lo+Li);
            if (LXmn2 > LXmx2) continue;
            for (int LX = LXmn2; LX <= LXmx2; LX++) {
              // Parity filter: LDEL = Lo-Li, LX must have same parity as LXP+LDEL
              int LDEL = Lo - Li;
              // Check JTOCS has this entry
              int k = koffs_index(LDEL, LX, JP, JBT);
              if (k < 0) continue;

              double TEMP2;
              if (LX == LXP && JP == JBP && Li == Li && Lo == Lo) {
                TEMP2 = SAV9J;
              } else {
                TEMP2 = std::sqrt(double((jpi+1)*(jpo+1)*(2*LX+1)*(JP+1)))
                  * NineJ(JBT/2.0, (double)LX, JP/2.0,
                          jpi/2.0, (double)Li,  JA/2.0,
                          jpo/2.0, (double)Lo,  JB/2.0);
              }
              if (std::abs(TEMP2) < 1e-14) continue;
              if ((LX + LXP) % 2 != 0) TEMP2 = -TEMP2;

              // Accumulate complex S into a temporary, convert to mag+phase at end
              // For now: use a separate complex accumulator per (k, Li_idx)
              // We'll do this properly below using a pre-accumulation map
              // Temporary: mark for later
              (void)TEMPR; (void)TEMPI; (void)TEMP2; (void)k;
            }
          }
        }
      }
    }
  }

  // ── SFROMI proper implementation using complex accumulator ────────────
  // Reset SMAG/SPHASE and redo with complex accumulation
  std::vector<std::vector<std::complex<double>>> S_acc(NSPL,
    std::vector<std::complex<double>>(NUMLIS, {0.0, 0.0}));

  for (const auto& elem : TransferSMatrix) {
    int LXP = elem.Lx;
    int Li  = elem.Li;
    int Lo  = elem.Lo;
    int JPI = elem.JPI;
    int JPO = elem.JPO;

    if (LXP < LxMin_bs || LXP > LxMax_bs) continue;
    if (Li < LMIN_dc || Li > LMAX_dc) continue;
    int Li_idx = Li - LMIN_dc;
    if ((Li + Lo + LXP) % 2 != 0) continue;

    // InelDcFaithful2 already computed:
    //   S_val = phase_factor * FACTOR_sf * ATERM * I_accum / sqrt(2*Li+1)
    // So elem.S IS the (SR + i*SI) ready for the 9-J loop.
    double SR = elem.S.real();
    double SI = elem.S.imag();

    if (!SOSWT) {
      int LDEL = Lo - Li;
      int k = koffs_index(LDEL, LXP, JBP, JBT);
      if (k < 0 || Li_idx >= NUMLIS) continue;
      S_acc[k][Li_idx] += std::complex<double>(SR, SI);
    } else {
      int JPIMN = JPI, JPIMX = JPI;
      if ((Incoming.Pot.VSO == 0.0)) { JPIMN = std::abs(2*Li-JA); JPIMX = 2*Li+JA; }
      int JPOMN = JPO, JPOMX = JPO;
      if ((Outgoing.Pot.VSO == 0.0)) { JPOMN = std::abs(2*Lo-JB); JPOMX = 2*Lo+JB; }

      for (int jpi = JPIMN; jpi <= JPIMX; jpi += 2) {
        for (int jpo = JPOMN; jpo <= JPOMX; jpo += 2) {
          double SAV9J = std::sqrt(double((jpi+1)*(jpo+1)*(2*LXP+1)*(JBP+1)))
            * NineJ(JBT/2.0,(double)LXP,JBP/2.0,
                    jpi/2.0,(double)Li, JA/2.0,
                    jpo/2.0,(double)Lo, JB/2.0);
          if (std::abs(SAV9J) < 1e-14) continue;
          double TEMPR = SAV9J * SR, TEMPI = SAV9J * SI;

          for (int JP = JPBASE; JP <= JPMX; JP += 2) {
            int LXmn2 = std::max({std::abs(JBT-JP)/2, std::abs(Lo-Li)});
            int LXmx2 = std::min((JBT+JP)/2, Lo+Li);
            for (int LX = LXmn2; LX <= LXmx2; LX++) {
              int LDEL = Lo - Li;
              int k = koffs_index(LDEL, LX, JP, JBT);
              if (k < 0 || Li_idx >= NUMLIS) continue;

              double TEMP2;
              if (LX==LXP && JP==JBP) {
                TEMP2 = SAV9J;
              } else {
                TEMP2 = std::sqrt(double((jpi+1)*(jpo+1)*(2*LX+1)*(JP+1)))
                  * NineJ(JBT/2.0,(double)LX,JP/2.0,
                          jpi/2.0,(double)Li, JA/2.0,
                          jpo/2.0,(double)Lo, JB/2.0);
              }
              if (std::abs(TEMP2) < 1e-14) continue;
              if ((LX + LXP) % 2 != 0) TEMP2 = -TEMP2;

              S_acc[k][Li_idx] += std::complex<double>(TEMP2*TEMPR, TEMP2*TEMPI);
            }
          }
        }
      }
    }
  }

  // Convert S_acc to SMAG/SPHASE
  for (int k = 0; k < NSPL; k++) {
    for (int idx = 0; idx < NUMLIS; idx++) {
      double re = S_acc[k][idx].real(), im = S_acc[k][idx].imag();
      SMAG  [k][idx] = std::sqrt(re*re + im*im);
      SPHASE[k][idx] = std::atan2(im, re);
    }
  }


  // Debug: count non-zero SMAG entries
  {
    int dbg_n = 0;
    double dbg_max = 0;
    for (int k=0;k<NSPL;k++) for(int i=0;i<NUMLIS;i++) {
      if(SMAG[k][i]>0) dbg_n++;
      dbg_max = std::max(dbg_max, SMAG[k][i]);
    }
  }

  // Debug: print specific SMAG entry
  for (int k=0; k<NSPL; k++) {
    if (JTOCS[k].LDEL==-2 && JTOCS[k].LX==2 && JTOCS[k].JP==1 && JTOCS[k].MX==0) {
      // Print first 10 Li values
      for (int i=0; i<std::min(NUMLIS,10); i++) {
      }
    }
  }
  // ── Coulomb phases sig_in(Li), sig_out(Lo) ────────────────────────────
  int Lcol = std::max(LOMOST, LIMOST) + 5;
  std::vector<double> sig_in (Lcol+1, 0.0);
  std::vector<double> sig_out(Lcol+1, 0.0);
  for (int L = 1; L <= Lcol; L++) {
    sig_in [L] = sig_in [L-1] + std::atan2(Incoming.eta, (double)L);
    sig_out[L] = sig_out[L-1] + std::atan2(Outgoing.eta, (double)L);
  }
  auto sigin  = [&](int L) { return sig_in [std::min(L, Lcol)]; };
  auto sigout = [&](int L) { return sig_out[std::min(L, Lcol)]; };

  // ── Step 3: BETCAL ────────────────────────────────────────────────────
  // BETAS[2][NSPL][Lo-LMIN+1]  (2 for real/imag, NSPL channels, Lo range)
  // LOMN = LMIN_dc (= LBASE for reactions), LOMX = LIMOST = LMAX_dc
  int LOMN   = LMIN_dc;
  int LOMX   = LIMOST;
  int NLO    = LOMX - LOMN + 1;  // number of Lo values

  // LDELMX = max |LDEL| across all JTOCS entries
  int LDELMX = 0;
  for (auto& e : JTOCS) LDELMX = std::max(LDELMX, std::abs(e.LDEL));

  // LXMN, LXMX from JTOCS
  int LXMN = 1000, LXMX = -1;
  for (auto& e : JTOCS) { LXMN = std::min(LXMN, e.LX); LXMX = std::max(LXMX, e.LX); }
  int NMLX = LXMX - LXMN + 1;
  int NMX  = LXMX + 1;  // Fortran: NMX = LXMX+1

  // BETAS[k][ilo] where ilo = Lo - LOMN (0-based)
  std::vector<std::vector<std::complex<double>>> BETAS(NSPL,
    std::vector<std::complex<double>>(NLO, {0.0, 0.0}));

  // IODD from LDELMX
  int IODD = std::abs(LDELMX) % 2;
  double FACTOR_bet = 0.5 / ka;
  if (IODD != 0) FACTOR_bet = -FACTOR_bet;

  // Outer loop: Lo
  for (int Lo = LOMN; Lo <= LOMX; Lo++) {
    double DLO  = Lo;
    int    ILO  = Lo - LOMN;  // 0-based index into BETAS second dim

    // Precompute TEMPS[LX-LXMN][Li_offset][MX] = FACTOR*(2Li+1)*CG
    // Fortran: DO 189 LI=LIMN,LIMX,2; DO 179 LX=LX1,LXMX; DO 169 MX=MXZ,LX
    //   TEMPS(I+MX) = FACTOR*(2LI+1)*CG(2LI,2LX,0,2MX,2Lo,2MX)
    // I = 1 + NMX*(LX-LXMN + NMLX*(LI-LIMN)/2)
    int LIMN_bet = Lo - LDELMX;
    int LIMX_bet = Lo + LDELMX;

    // Build TEMPS array (flattened: [LX-LXMN, (Li-LIMN)/2, MX])
    // We use a map for simplicity; LIMN_bet can be negative
    std::map<std::tuple<int,int,int>, double> TEMPS_map;  // (LX, Li, MX) -> value

    for (int Li = LIMN_bet; Li <= LIMX_bet; Li += 2) {
      if (Li < 0) continue;
      int LX1 = std::max(std::abs(Li - Lo), LXMN);
      for (int LX = LX1; LX <= LXMX; LX++) {
        // Fortran BETCAL line 3958: MXZ = MOD(LX+LI-LO, 2) — PARITY ONLY (0 or 1)
        // This is the starting MX for the TEMPS CG precomputation.
        // CG(Li,0;LX,MX|Lo,MX) is zero for wrong parity, so starting from 0 or 1 is fine.
        int MXZ_start = ((LX + Li - Lo) % 2 + 2) % 2;  // 0 or 1
        for (int MX = MXZ_start; MX <= LX; MX++) {
          if (Lo < MX) continue;
          double cg = ClebschGordan((double)Li, 0.0, (double)LX, (double)MX,
                                    (double)Lo, (double)MX);
          if (!std::isfinite(cg) || std::abs(cg) < 1e-14) continue;
          TEMPS_map[{LX, Li, MX}] = FACTOR_bet * (2.0*Li + 1.0) * cg;
        }
      }
    }

    // Inner loop: KOFFS
    for (int k = 0; k < NSPL; k++) {
      int LDEL = JTOCS[k].LDEL;
      int LX   = JTOCS[k].LX;
      int Li   = Lo - LDEL;  // Li = Lo - (Lo-Li) = Li ✓

      if (Li < LMIN_dc || Li > LMAX_dc) continue;
      int Li_idx = Li - LMIN_dc;

      double amag = SMAG[k][Li_idx];
      if (amag == 0.0) continue;

      // Coulomb phase: PHASE = SPHASE(k,Li_idx) + sig_in(Li) + sig_out(Lo)
      double phase = SPHASE[k][Li_idx] + sigin(Li) + sigout(Lo);
      // S*e^{i*phase} / i: SMATR = amag*sin(phase), SMATI = -amag*cos(phase)
      double SMATR = amag * std::sin(phase);
      double SMATI =-amag * std::cos(phase);

      // MXZ = LX + LDEL (Fortran line 3520)
      int MXZ = LX + LDEL;
      // KOFFZ = KOFFS - MXZ (Fortran line 3521, KOFFS is 1-based there)
      // In Fortran: BETAS(1, KOFFZ+MX, ILO) where KOFFZ+MX = KOFFS (when MX steps from MXZ)
      // So the BETAS index for Mx=MXZ is KOFFZ+MXZ = KOFFS (same KOFFS!)
      // This means BETAS[k][ILO] is indexed by KOFFS directly when MX=MXZ
      // For MX > MXZ: index is KOFFS+(MX-MXZ)... but in C++ we use k as KOFFS directly
      // PROBLEM: the Fortran BETAS array uses a single index that spans all KOFFS+Mx combinations
      // KOFFZ = KOFFS - MXZ, so KOFFZ + MX = KOFFS + (MX - MXZ)
      // When MX = MXZ: slot = KOFFS
      // When MX = MXZ+1: slot = KOFFS+1
      // This means Fortran BETAS uses KOFFS as the Mx=MXZ slot and KOFFS+1 as Mx=MXZ+1 slot
      // The JTOCS ordering guarantees that KOFFS, KOFFS+1, ..., KOFFS+LX-MXZ are all
      // the Mx=MXZ,MXZ+1,...,LX entries for the SAME (LDEL, LX, JP, JT)
      // In our C++ JTOCS, we store each (LDEL, LX, JP, JT, MX) as a SEPARATE entry!
      // So JTOCS[k] already has a specific MX. We just add to BETAS[k][ILO].

      // The Fortran TEMPS loop: for MX = MXZ..LX, store into BETAS[KOFFZ+MX]
      // which is BETAS[KOFFS + (MX-MXZ)] = BETAS at the slot for Mx=MX
      // In our C++ JTOCS, k already encodes a specific Mx. So:
      int MX = JTOCS[k].MX;
      if (MX < MXZ || MX > LX) continue;  // loop condition MX=MXZ..LX

      auto it = TEMPS_map.find({LX, Li, MX});
      if (it == TEMPS_map.end()) continue;
      double T = it->second;

      BETAS[k][ILO] += std::complex<double>(T * SMATR, T * SMATI);
    }

    // sqrt-factorial pass: multiply BETAS by TEMPS_fact(MX)
    // TEMPS_fact(MX=0) = 1, TEMPS_fact(MX) = prod_{n=1}^{MX} 1/sqrt((Lo+n)(Lo-n+1))
    int MXMX = std::min((int)LXMX, Lo);
    std::vector<double> TEMPS_fact(MXMX+2, 0.0);
    TEMPS_fact[0] = 1.0;  // MX=0 slot
    for (int MX = 1; MX <= MXMX; MX++) {
      double denom = (DLO + MX) * (DLO - MX + 1.0);
      if (denom <= 0.0) { TEMPS_fact[MX] = 0.0; break; }
      TEMPS_fact[MX] = TEMPS_fact[MX-1] / std::sqrt(denom);
    }
    for (int k = 0; k < NSPL; k++) {
      int MX = JTOCS[k].MX;
      if (MX > MXMX || MX < 0) continue;
      BETAS[k][ILO] *= TEMPS_fact[MX];
    }
  }


  // Debug BETAS
  {
    int dbg_n = 0;
    double dbg_max = 0;
    for (int k=0;k<NSPL;k++) for(int i=0;i<NLO;i++) {
      double v = std::abs(BETAS[k][i]);
      if(v>0) dbg_n++;
      dbg_max = std::max(dbg_max, v);
    }
  }
  // ── Step 4: AMPCAL ────────────────────────────────────────────────────
  // F[k] = sum_Lo BETAS[k][ilo] * PLM(Lo, MX, cos_theta) * FACMBL
  // DCS(theta) = 10 * sum_k FMNEG(MX) * |F[k]|^2

  std::cout << "\nAngle (deg)    dSigma/dOmega (mb/sr)\n";

  for (double theta_deg = AngleMin; theta_deg <= AngleMax + 1e-9; theta_deg += AngleStep) {
    double cos_theta = std::cos(theta_deg * M_PI / 180.0);

    // Compute PLM(Lo, MX, cos_theta) for Lo=LOMN..LOMX, MX=0..LXMX
    // Fortran: PLM(LO+LPLM) where LPLM = MX*(2*LMX+1-MX)/2+1
    // Storage: P(L,M) = PLM[L + 1 + M*(2*LMX+1-M)/2]  (1-based in Fortran)
    // We use a 2D array: plm[Lo][MX]
    int LMX_plm = LOMX;
    int MMAX_plm = LXMX;
    // PLM array: plm_arr[L*(MMAX+1) + M] but use Fortran layout
    // P(L,M) is stored at index L+1 + M*(2*LMX+1-M)/2 in Fortran (1-based)
    // → 0-based: L + M*(2*LMX+1-M)/2
    int plm_size = (LMX_plm+1) + MMAX_plm*(2*LMX_plm+1-MMAX_plm)/2 + 2;
    std::vector<double> PLM_arr(plm_size+1, 0.0);

    // Fill PLM using same recursion as Fortran PLMSUB
    // P(0,0)=1, P(L+1,M) = ((2L+1)*x*P(L,M) - (L+M)*P(L-1,M))/(L-M+1)
    // P(M,M) = -(2M-1)*ROOT*P(M-1,M-1), ROOT=sqrt(1-x^2)
    // Index: P(L,M) → PLM_arr[L + M*(2*LMX_plm+1-M)/2]
    double ROOT = std::sqrt(std::max(0.0, 1.0 - cos_theta*cos_theta));

    // Initialize
    PLM_arr[0] = 1.0;  // P(0,0)
    // Build for M=0 first: P(L,0) = Legendre polynomial
    {
      if (LMX_plm >= 1) PLM_arr[1] = cos_theta;
      for (int L = 2; L <= LMX_plm; L++) {
        PLM_arr[L] = ((2*L-1)*cos_theta*PLM_arr[L-1] - (L-1)*PLM_arr[L-2]) / L;
      }
    }
    // Build for M=1,2,...,MMAX_plm
    // P(M,M) starts from P(0,0)=1
    double PMM = 1.0;
    for (int M = 1; M <= MMAX_plm; M++) {
      PMM = -(2*M-1) * ROOT * PMM;  // P(M,M) from P(M-1,M-1)
      int base = M * (2*LMX_plm+1-M) / 2;
      PLM_arr[base + M] = PMM;      // P(M,M)
      if (M+1 <= LMX_plm) {
        PLM_arr[base + M+1] = (2*M+1) * cos_theta * PMM;  // P(M+1,M)
      }
      for (int L = M+2; L <= LMX_plm; L++) {
        int prev_base = (M > 0 ? M*(2*LMX_plm+1-M)/2 : 0);
        // Recursion: P(L,M) = ((2L-1)*x*P(L-1,M) - (L+M-1)*P(L-2,M))/(L-M)
        double Pm2 = PLM_arr[base + L-2];
        double Pm1 = PLM_arr[base + L-1];
        PLM_arr[base + L] = ((2.0*L-1)*cos_theta*Pm1 - (L+M-1)*Pm2) / (L-M);
        (void)prev_base;
      }
    }

    auto get_PLM = [&](int L, int M) -> double {
      if (L < 0 || M < 0 || M > L || M > MMAX_plm || L > LMX_plm) return 0.0;
      int base = M * (2*LMX_plm+1-M) / 2;
      return PLM_arr[base + L];
    };

    double dSigma = 0.0;

    for (int k = 0; k < NSPL; k++) {
      int MX = JTOCS[k].MX;
      int LOMNMX = std::max(LOMN, MX);

      std::complex<double> F_amp(0.0, 0.0);
      for (int Lo = LOMNMX; Lo <= LOMX; Lo++) {
        int ILO = Lo - LOMN;
        double plmval = get_PLM(Lo, MX);
        F_amp += BETAS[k][ILO] * plmval;
      }

      double FMNEG = (MX == 0) ? 1.0 : 2.0;
      dSigma += 10.0 * FMNEG * std::norm(F_amp);
    }

    std::cout << std::fixed << std::setprecision(1) << theta_deg
              << "           " << std::scientific << std::setprecision(4)
              << dSigma << "\n";
  }

}
