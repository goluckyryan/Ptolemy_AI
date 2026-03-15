// wavelj.cpp — Distorted wave setup and Numerov integration
//
// Extracted from dwba.cpp (original lines 342–567):
//   DWBA::WavSet  (lines 342–376)
//   DWBA::WavElj  (lines 377–567)
//
// These functions port Ptolemy's WAVSET and WAVELJ subroutines to C++.
// WAVSET allocates and fills the radial potential grid for a channel.
// WavElj performs the outward Numerov integration for partial wave (L, Jp),
// extracts the S-matrix element via Wronskian matching, and normalizes the
// stored wavefunction to the Coulomb asymptotic form.

#include "dwba.h"
#include "potential_eval.h"
#include "rcwfn.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// DWBA::WavSet
//
// Allocates the radial grid and evaluates the optical potential at each point.
// Matches Ptolemy's WAVSET: sets step size H = 0.1 fm, MaxR = 30 fm.
// EvaluatePotential fills V_real, V_imag, V_so_real, V_so_imag, V_coulomb.
// ---------------------------------------------------------------------------
void DWBA::WavSet(Channel &ch) {
  ch.StepSize = 0.1;
  ch.MaxR = 30.0;
  ch.NSteps = static_cast<int>(ch.MaxR / ch.StepSize) + 1;

  ch.RGrid.resize(ch.NSteps);
  ch.V_real.resize(ch.NSteps);
  ch.V_imag.resize(ch.NSteps);
  ch.V_so_real.resize(ch.NSteps);
  ch.V_so_imag.resize(ch.NSteps);
  ch.V_coulomb.resize(ch.NSteps);

  double A_target = ch.Target.A;
  double A_projectile = ch.Projectile.A;

  for (int i = 0; i < ch.NSteps; ++i) {
    double r = i * ch.StepSize;
    ch.RGrid[i] = r;

    if (r < 0.001) {
      ch.V_real[i] = 0.0;
      ch.V_imag[i] = 0.0;
      ch.V_so_real[i] = 0.0;
      ch.V_so_imag[i] = 0.0;
      ch.V_coulomb[i] = 0.0;
      continue;
    }

    EvaluatePotential(r, ch.Pot, ch.V_real[i], ch.V_imag[i], ch.V_so_real[i],
                      ch.V_so_imag[i], ch.V_coulomb[i], ch.Projectile.Z,
                      ch.Target.Z, A_target, A_projectile);
  }
}

// ---------------------------------------------------------------------------
// DWBA::WavElj
//
// Outward Numerov integration matching Ptolemy's WAVELJ subroutine.
//
// Key insight from Ptolemy source (source.mor ~30793):
//   When |u| exceeds BIGNUM during integration, Ptolemy scales all stored
//   values by 1/|u| AND ZEROS OUT the inner region (where |u| < STEPI
//   after rescaling). This means the stored wavefunction is only nonzero
//   from some ISTRT outward.  The normalization alpha = chi_last / u_last
//   then correctly maps these values to the Coulomb-normalized form.
//
// The critical fix vs previous code:
//   We must also zero the inner region after rescaling, so that the
//   inner values (which grew exponentially due to centrifugal barrier)
//   do NOT contaminate the final alpha-normalized wavefunction.
//
// Spin-orbit coupling factor (Ptolemy WAVELJ convention):
//   SDOTL = 0.25*(JP*(JP+2) - JSPS*(JSPS+2)) - L*(L+1)
//   spin_dot_L = SDOTL / JSPS
// where JP = Jp (doubled integer), JSPS = ch.JSPS (doubled integer, 2*spin).
// For proton (JSPS=1):  JP=2L+1 → spin_dot_L = L     (J=L+1/2);
//                        JP=2L-1 → spin_dot_L = -(L+1) (J=L-1/2)
// For deuteron (JSPS=2): JP=2L+2 → spin_dot_L = L     (J=L+1);
//                         JP=2L   → spin_dot_L = -1    (J=L);
//                         JP=2L-2 → spin_dot_L = -(L+1) (J=L-1)
// ---------------------------------------------------------------------------
void DWBA::WavElj(Channel &ch, int L, int Jp) {
  int N = ch.NSteps;
  double h = ch.StepSize;
  double h2 = h * h;
  double h2_12 = h2 / 12.0;

  ch.WaveFunction.assign(N + 1, std::complex<double>(0.0, 0.0));

  // --- Spin-orbit coupling factor ---
  int JSP_ch = (ch.JSPS > 0) ? ch.JSPS : 1;  // default to 1 if not set

  // Ptolemy WAVELJ: if SOSWS && JP outside valid range, return (invalid JP)
  // JP must be in [|2L-JSPS|, 2L+JSPS] with same parity as JSPS
  int JP_min_valid = std::abs(2 * L - JSP_ch);
  int JP_max_valid = 2 * L + JSP_ch;
  if (Jp < JP_min_valid || Jp > JP_max_valid || ((Jp + JSP_ch) % 2 != 0)) {
    // Invalid JP for this channel — zero out wavefunction and return
    ch.WaveFunction.assign(N + 1, std::complex<double>(0.0, 0.0));
    return;
  }

  double DL2 = (double)L * (L + 1);
  double SDOTL_raw =
      0.25 * double(Jp * (Jp + 2) - JSP_ch * (JSP_ch + 2)) - DL2;
  double spin_dot_L = (JSP_ch > 0) ? SDOTL_raw / JSP_ch : 0.0;

  // Physical constants (same as dwba.cpp globals)
  const double AMU_MEV  = 931.494061;
  const double HBARC_L  = 197.32697;
  double f_conv = 2.0 * ch.mu * AMU_MEV / (HBARC_L * HBARC_L);
  double k2 = ch.k * ch.k;

  // --- Build f(r) array (complex: absorptive imaginary potential) ---
  // Schrodinger equation: u'' + f(r)*u = 0
  // f(r) = k² - l(l+1)/r² - f_conv*V_coulomb + f_conv*V_real
  //        + spin_dot_L * f_conv * V_so_real
  //        + i*(f_conv*V_imag + spin_dot_L * f_conv * V_so_imag)
  // EvaluatePotential returns positive depths: V_real>0, W>0, Vc>0.
  // Nuclear potential is attractive → adds to k² with + sign.
  // Coulomb is repulsive → subtracts from k² with - sign.
  std::vector<std::complex<double>> f(N + 2, 0.0);
  for (int i = 0; i <= N; ++i) {
    // For i >= NSteps, extrapolate; i=0 handled with r=epsilon
    double r = (i == 0) ? 1e-10
               : ((i < ch.NSteps) ? ch.RGrid[i] : i * ch.StepSize);
    if (r < 1e-10) r = 1e-10;

    int idx = std::min(i, N - 1);
    double Vr   = ch.V_real[idx];
    double Wi   = ch.V_imag[idx];
    double Vc   = ch.V_coulomb[idx];
    double Vso  = ch.V_so_real[idx];
    double Vsoi = ch.V_so_imag[idx];
    double LL1  = (double)L * (L + 1);

    double f_re = k2 - LL1 / (r * r) - f_conv * Vc + f_conv * Vr
                  + spin_dot_L * f_conv * Vso;
    double f_im = f_conv * Wi + spin_dot_L * f_conv * Vsoi;
    f[i] = std::complex<double>(f_re, f_im);

    // Debug output for L=0 at key radii
    if (L == 0 && (i == 10 || i == 30 || i == 50 || i == 80 || i == 100)) {
      std::printf("  WavElj L=0 i=%3d r=%5.2f: Vr=%8.3f Wi=%7.3f Vc=%6.3f "
                  "f_re=%9.5f f_im=%+9.5f\n",
                  i, r, Vr, Wi, Vc, f_re, f_im);
    }
  }

  // --- Standard Numerov integration (Ptolemy WAVELJ form) ---
  // Recurrence: (1 + h²/12·f_{i+1}) u_{i+1}
  //           = 2*(1 - 5h²/12·f_i)*u_i - (1 + h²/12·f_{i-1})*u_{i-1}
  //
  // Initial condition: u(0)=0, u(1)=h^{L+1}  (regular solution at origin)
  const double BIGNUM = 1e30;   // overflow threshold
  const double STEPI  = 1e-10;  // inner cutoff after rescaling

  std::vector<std::complex<double>> u(N + 2, 0.0);
  u[0] = 0.0;
  u[1] = std::pow(h, L + 1);
  if (std::abs(u[1]) < 1e-300) u[1] = 1e-300;

  int ISTRT = 0;  // first nonzero index (tracks inner zeroing after rescale)

  for (int i = 1; i < N; ++i) {
    std::complex<double> t1 = 2.0 * (1.0 - 5.0 * h2_12 * f[i]) * u[i];
    std::complex<double> t2 = (1.0 + h2_12 * f[i - 1]) * u[i - 1];
    std::complex<double> dn = 1.0 + h2_12 * f[i + 1];
    if (std::abs(dn) < 1e-300) dn = 1e-300;
    u[i + 1] = (t1 - t2) / dn;

    double mag = std::abs(u[i + 1]);
    if (mag > BIGNUM) {
      // Rescale: multiply all stored values by 1/mag
      // AND zero out the inner region where |u * (1/mag)| < STEPI
      double inv_mag = 1.0 / mag;
      double threshold = STEPI * mag;  // |u| must exceed this to survive

      // Find new ISTRT: scan forward from old ISTRT to find first nonzero
      int new_istrt = i + 1;  // default: everything up to now is zeroed
      for (int j = ISTRT; j <= i + 1; ++j) {
        if (std::abs(u[j]) >= threshold) {
          new_istrt = j;
          break;
        }
      }

      // Zero the inner region
      for (int j = ISTRT; j < new_istrt; ++j) u[j] = 0.0;
      ISTRT = new_istrt;

      // Scale the surviving region
      for (int j = ISTRT; j <= i + 1; ++j) u[j] *= inv_mag;
    }
  }

  // --- S-matrix extraction: handbook single-point method ---
  // Match at index n_match = N-3 so 5-point stencil fits within [0..N]
  int n_match = N - 3;
  double R_match = n_match * h;
  double rho_match = ch.k * R_match;

  // Get Coulomb functions and derivatives at R_match
  std::vector<double> FC(L + 2), FCP(L + 2), GC(L + 2), GCP(L + 2);
  Rcwfn(rho_match, ch.eta, L, L, FC, FCP, GC, GCP);
  double FL  = FC[L];
  double GL  = GC[L];
  double FLP = FCP[L] * ch.k;   // dF/dr = (dF/drho)*(drho/dr) = FP * k
  double GLP = GCP[L] * ch.k;   // same for G

  // 5-point stencil for u'(R_match)
  std::complex<double> u_nm2 = u[n_match - 2];
  std::complex<double> u_nm1 = u[n_match - 1];
  std::complex<double> u_np1 = u[n_match + 1];
  std::complex<double> u_np2 = u[n_match + 2];
  std::complex<double> u_prime = (u_nm2 - 8.0*u_nm1 + 8.0*u_np1 - u_np2) / (12.0 * h);

  std::complex<double> u_m = u[n_match];
  double ur = u_m.real(), ui = u_m.imag();
  double upr = u_prime.real(), upi = u_prime.imag();

  // Wronskians: A = W(F,u), B = W(u,G)
  double Ar = FLP * ur - upr * FL;
  double Ai = FLP * ui - upi * FL;
  double Br = upr * GL - GLP * ur;
  double Bi = upi * GL - GLP * ui;

  // S = (B + iA) / (B - iA)
  double num_r = Br - Ai;
  double num_i = Bi + Ar;
  double den_r = Br + Ai;
  double den_i = Bi - Ar;
  double den   = den_r*den_r + den_i*den_i;

  double SJR = (den > 1e-60) ? (num_r*den_r + num_i*den_i) / den : 0.0;
  double SJI = (den > 1e-60) ? (num_i*den_r - num_r*den_i) / den : 0.0;

  if (ch.SMatrix.size() <= (size_t)L) ch.SMatrix.resize(L + 1);
  ch.SMatrix[L] = std::complex<double>(SJR, SJI);

  // Debug output for selected partial waves
  if (L == 12 || L == 7 || L == 5 || L == 2 || L == 0) {
    std::printf("  [WavElj L=%2d JP=%d] FL=%.5f GL=%.5f FLP=%.5f GLP=%.5f\n",
                L, Jp, FL, GL, FLP, GLP);
    std::printf("               ur=%.4e ui=%.4e upr=%.4e upi=%.4e\n",
                ur, ui, upr, upi);
    std::printf("               Ar=%.4e Ai=%.4e Br=%.4e Bi=%.4e den=%.4e\n",
                Ar, Ai, Br, Bi, den);
    std::printf("               |S|=%.5f SJR=%.5f SJI=%.5f\n",
                std::sqrt(SJR * SJR + SJI * SJI), SJR, SJI);
    std::printf("               R_match=%.3f k=%.4f eta=%.4f\n",
                R_match, ch.k, ch.eta);
  }

  // --- Normalization (Ptolemy source.mor lines 30913–30916) ---
  // A1n = 0.5*(F*(1+SJR) + SJI*G)
  // A2n = 0.5*(G*(1-SJR) + SJI*F)
  // alpha = (u[n_match] · A1n + conj-part A2n) / |u[n_match]|²
  // Uses Coulomb functions at matching radius and u at matching point.
  double A1n  = 0.5 * (FL * (1.0 + SJR) + SJI * GL);
  double A2n  = 0.5 * (GL * (1.0 - SJR) + SJI * FL);
  double den2 = ur * ur + ui * ui;  // |u[n_match]|²

  if (den2 > 1e-60) {
    double alpha_r = (ur * A1n + ui * A2n) / den2;
    double alpha_i = (ur * A2n - ui * A1n) / den2;
    std::complex<double> alpha(alpha_r, alpha_i);

    for (int i = 0; i <= N; ++i)
      ch.WaveFunction[i] = u[i] * alpha;
  }

  if (std::isnan(SJR) || std::isnan(SJI))
    std::cerr << "WavElj NaN: L=" << L << " k=" << ch.k
              << " eta=" << ch.eta << std::endl;
}
