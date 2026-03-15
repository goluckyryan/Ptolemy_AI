// a12.cpp — DWBA::ComputeA12Terms() and DWBA::EvalA12()
// Angular coupling kernel for DWBA transfer integral.
// Extracted from ineldc.cpp (A12 precomputation and evaluation).

#include "dwba.h"
#include "math_utils.h"
#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

// xlam(L, |M|) = Wigner d^L_{|M|,0}(pi/2) — verified by standalone Fortran test.
// Ptolemy source.mor lines 1619-1644 (HALFSW=FALSE, even-M path).
// Recurrence: start OUTTER=1, m_cur=0.
//   For each LL=1..L: m_cur = 1-m_cur (alternates 0→1→0→1...)
//     OUTTER *= sqrt((LL+m_cur-1)/(LL+m_cur))
//     if m_cur==1: OUTTER = -OUTTER  [sign convention]
//     xlam(LL, m_cur) = OUTTER
//     For MM = m_cur+2, m_cur+4, ..., LL:
//       xlam(LL, MM) = -xlam(LL, MM-2) * sqrt((LL-MM+2)*(LL+MM-1)/((LL+MM)*(LL-MM+1)))
// This exactly reproduces d^L_{M,0}(pi/2).
static double xlam_correct(int L, int am) {
  am = std::abs(am);
  if (am > L) return 0.0;
  if ((L + am) % 2 != 0) return 0.0;  // must have same parity

  // Build the table up to level L using Ptolemy recurrence
  double outter = 1.0;
  int m_cur = 0;
  static double xlam_cache[20][20];  // [ll][mm/2]
  xlam_cache[0][0] = 1.0;  // xlam(0,0)=1

  for (int ll = 1; ll <= L; ll++) {
    m_cur = 1 - m_cur;   // alternate 0→1→0→...
    outter *= std::sqrt((double)(ll + m_cur - 1) / (double)(ll + m_cur));
    if (m_cur == 1) outter = -outter;
    xlam_cache[ll][m_cur / 2] = outter;  // store xlam(ll, m_cur)
    // Higher M by inner recursion (step 2)
    for (int mm = m_cur + 2; mm <= ll; mm += 2) {
      int idx = mm / 2;
      xlam_cache[ll][idx] = -xlam_cache[ll][idx - 1] *
        std::sqrt((double)(ll - mm + 2) * (double)(ll + mm - 1) /
                 ((double)(ll + mm) * (double)(ll - mm + 1)));
    }
  }
  return xlam_cache[L][am / 2];
}

std::vector<std::tuple<int,int,double>> DWBA::ComputeA12Terms(
    int Li, int Lo, int Lx, int lT, int lP) {
  // Ptolemy A12 structure (HALFSW=FALSE, i.e. LBT and LBP both even):
  //   XN  = 0.5 * sqrt((2Li+1)*(2LBT+1)*(2LBP+1))
  //   MT  loops from -LBT to +LBT step 2 (EVEN values only!)
  //   MX  = MT + MP  (MP=0 for LBP=lP=0)
  //   OUTTMP(LX,MT) = XLAM(LBT,MT) * XLAM(LBP,MP=0) * XN * ThreeJ(LBT,MT,LBP,0,LX,-MX)
  //   MU  loops from MOD(LI,2) to LI step 2 (same parity as Li)
  //   MUPOS = |MX-MU|
  //   OUTTER = XLAM(LI,MU) * XLAM(LO,MUPOS) * sqrt(2*LO+1)
  //   TTT = ThreeJ(LI,MU,LO,MX-MU,LX,-MX)
  //   A12 += OUTTMP * OUTTER * TTT * doubling_MU   [doubled for MU!=0]

  std::vector<std::tuple<int,int,double>> A12_terms;

  double XN_a12 = 0.5 * std::sqrt((2.0*Li+1.0)*(2.0*lT+1.0)*(2.0*lP+1.0));

  // Compute LOMNMN: minimum Lo with right parity
  int min_Lo_tri = std::abs(Li - Lx);
  int lo_parity = Li % 2;
  int LOMNMN = min_Lo_tri;
  if (LOMNMN % 2 != lo_parity) LOMNMN++;

  for (int MT = -lT; MT <= lT; MT += 2) {
    int Mx = MT;  // MX = MT + MP, MP=0
    double xlam_T = xlam_correct(lT, std::abs(MT));
    double xlam_P = xlam_correct(lP, 0);
    double outtmp = xlam_T * xlam_P * XN_a12
                  * ThreeJ((double)lT, (double)MT, (double)lP, 0.0, (double)Lx, (double)(-Mx));
    if (std::abs(outtmp) < 1e-15) continue;

    // MU: same parity as Li, from (Li%2) to Li step 2 (always non-negative)
    for (int MU = (Li % 2); MU <= Li; MU += 2) {
      int MX_minus_MU = Mx - MU;
      if (std::abs(MX_minus_MU) > Lo) continue;
      double xlam_Li = xlam_correct(Li, std::abs(MU));
      double xlam_Lo = xlam_correct(Lo, std::abs(MX_minus_MU));
      double outter = xlam_Li * xlam_Lo * std::sqrt(2.0*Lo + 1.0);
      // ITES sign correction (Ptolemy lines 1838-1853):
      int ITES = 0;
      if (MT < 0) ITES += lT;
      if (MX_minus_MU < 0) ITES += LOMNMN;
      if (ITES % 2 != 0) outter = -outter;
      // TTT inner ThreeJ:
      double ttt = ThreeJ((double)Li, (double)MU, (double)Lo, (double)MX_minus_MU,
                          (double)Lx, (double)(-Mx));
      double A12_val = outtmp * outter * ttt;
      if (std::abs(A12_val) < 1e-15) continue;
      // Doubling for MU!=0
      double doubling = (MU != 0) ? 2.0 : 1.0;
      A12_terms.push_back({MT, MU, A12_val * doubling});
    }
  }

  return A12_terms;
}

double DWBA::EvalA12(
    const std::vector<std::tuple<int,int,double>>& A12_terms,
    double phi_T_angle, double phi_ab) {
  double A12_val = 0.0;
  for (auto &[MT_k, MU_k, coeff] : A12_terms) {
    double arg = MT_k * phi_T_angle - MU_k * phi_ab;
    A12_val += coeff * std::cos(arg);
  }
  return A12_val;
}
