// wavelj.cpp — Distorted wave setup and Numerov integration
//
// Single-point Wronskian implementation (shared between FR and ZR paths).
// The ZR reference copy is preserved in wavelj_ZR.cpp.
//
// NOTE: The 2-point Wronskian (Fortran NBAKCM=4) was investigated but found
// ill-conditioned at R~30 fm with h=0.125 fm (separation = 0.5 fm << wavelength).
// Fortran uses h~0.15 fm with matching at ~24 fm where the formula is well-conditioned.
// The 6.4° S-matrix phase error is caused by a potential V(R) discrepancy (~1.5% at R~3 fm)
// that must be fixed in the potential computation, not the Wronskian method.
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
// Step size follows Ptolemy's formula:
//   - Scattering channel (k > 0 set externally): h = min(2π/k, 1.0) / STEPSPER
//   - Bound state: h = min(1/kappa, A_diffuseness) / STEPSPER  (set in InelDc)
// Default: h = 0.1 fm (Ptolemy scattering h ≈ 0.125 fm at typical energies;
//   caller sets ch.StepSize before calling WavSet to override).
// EvaluatePotential fills V_real, V_imag, V_so_real, V_so_imag, V_coulomb.
// ---------------------------------------------------------------------------

// Ptolemy step size for a scattering channel:
//   h = min(2π/k, 1.0) / STEPSPER  (STEPSPER=8 from dpsb parameterset)
static double PtolemyScatStepSize(double k_fm, int stepsper = 8) {
  if (k_fm < 1e-10) return 0.1;
  double lambda = 2.0 * M_PI / k_fm;
  return std::min(lambda, 1.0) / stepsper;
}

// Ptolemy step size for a bound state channel:
//   h = min(1/kappa, A_diff) / STEPSPER
static double PtolemyBoundStepSize(double kappa_fm, double A_diff, int stepsper = 8) {
  if (kappa_fm < 1e-10) return 0.1;
  double asym_range = std::min(1.0 / kappa_fm, A_diff);
  return asym_range / stepsper;
}

void DWBA::WavSet(Channel &ch) {
  // Use caller-set StepSize if already specified (> 0); otherwise default 0.1 fm.
  // Scattering channels call WavSet with ch.k already set → caller can pre-set StepSize.
  // Bound state channels override StepSize in InelDc before calling WavSet.
  if (ch.StepSize <= 0.0) ch.StepSize = 0.1;
  // Use caller-set MaxR if > 0 (setup.cpp computes it from turning point);
  // otherwise default to 30.0 fm (bound state channels).
  if (ch.MaxR <= 0.0) ch.MaxR = 30.0;
  // Ptolemy: NSTEP = RMAX/STEPSZ + 0.5  (Fortran integer truncation = rounding)
  ch.NSteps = static_cast<int>(ch.MaxR / ch.StepSize + 0.5);

  ch.RGrid.resize(ch.NSteps);
  ch.V_real.resize(ch.NSteps);
  ch.V_imag.resize(ch.NSteps);
  ch.V_so_real.resize(ch.NSteps);
  ch.V_so_imag.resize(ch.NSteps);
  ch.V_coulomb.resize(ch.NSteps);
  ch.WaveFunction.assign(ch.NSteps + 2, {0.0, 0.0});  // pre-allocate for CalculateBoundState

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
  bool is_in = (&ch == &Incoming);
  if (L <= 3) {
    fprintf(stderr, "WAVELJ %s L=%d JP=%d N=%d h=%.6f Rmax=%.3f\n",
            is_in?"IN":"OUT", L, Jp, N, h, N*h);
  }
  double h2 = h * h;
  double h2_12 = h2 / 12.0;

  // +4 guard zeros at end for 5-point Lagrange interpolation in InelDc
  ch.WaveFunction.assign(N + 5, std::complex<double>(0.0, 0.0));

  // --- Spin-orbit coupling factor ---
  int JSP_ch = (ch.JSPS > 0) ? ch.JSPS : 1;  // default to 1 if not set

  // Ptolemy WAVELJ: if SOSWS && JP outside valid range, return (invalid JP)
  // JP must be in [|2L-JSPS|, 2L+JSPS] with same parity as JSPS
  int JP_min_valid = std::abs(2 * L - JSP_ch);
  int JP_max_valid = 2 * L + JSP_ch;
  if (Jp < JP_min_valid || Jp > JP_max_valid || ((Jp + JSP_ch) % 2 != 0)) {
    // Invalid JP for this channel — zero out wavefunction and return
    ch.WaveFunction.assign(N + 5, std::complex<double>(0.0, 0.0));
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
                  - spin_dot_L * f_conv * Vso;
    double f_im = f_conv * Wi - spin_dot_L * f_conv * Vsoi;
    f[i] = std::complex<double>(f_re, f_im);

    // Fortran stores WAVR(I+4) = h²/12 × f(r) = total Numerov potential
    // C++ f[i] = k² - L(L+1)/r² - f_conv×Vc + f_conv×Vr - sdotL×f_conv×Vso + i×(...)
    // Fortran WAVR(I+4) = VREAL(I+1) + DLSQ×VCENT(I+1) + SDOTL×ALLOC(LSOR+I)
    // where VREAL = -(h/ℏc)²×μ/6 × V_nuclear, VCENT = -(h/ℏc)²×μ/6 × 1/r²
    // So WAVR(I+4) = h²/(12E) × f(r) ... but with different sign convention
    // Actually: h²_12 = h² * μ / (6 ℏc²) = (StepSize/ℏc)² × μ_MeV / 6
    // And WAVR(I+4) = 1 + h²_12 × f(r) in Numerov form: W_i = 1 + h²/12 × f_i
    // No wait, Fortran WAVR = -h²/(12E) × V, so it IS h²_12 but with V/E factor
    // Let me just print the Numerov W(r) = 1 + h²/12 × f(r) for direct comparison
    if (L == 4 && Jp == 10 && i >= 2) {
      // Fortran: WAVR(I+4) is stored BEFORE Numerov integration
      // It's the quantity h²/12 × [V_nuc + L(L+1)/r² + sdotL×VSO - k²×12/h² ...]
      // Actually Fortran stores: -h²μ/(6ℏc²) × [V(r) + L(L+1)/r² × ℏc²/μ + sdotL×Vso]
      // = -(StepSize)²×μ/(6×ℏc²) × potentials
      // The Numerov form is WAVR(I+4) = total potential in Numerov units
      // In C++: h²_12 × f(r) gives the same quantity
      double W_re = h2_12 * f_re;  // This should match Fortran WAVR(I+4)
      double W_im = h2_12 * f_im;
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
  // Initial conditions: u[0]=0, u[1]=h^(L+1) (real regular solution near origin)
  // Note: Fortran uses u[1] = THISR*(1+i) — complex IC. For a complex optical potential,
  // this gives a slightly different WF phase at matching radius but the same S-matrix
  // (the Wronskian normalization removes the overall complex scale).
  u[0] = std::complex<double>(0.0, 0.0);
  u[1] = std::complex<double>(std::pow(h, L + 1), 0.0);
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
      // Rescale: multiply ALL stored values by 1/mag to prevent overflow.
      // DO NOT zero inner points here — that happens after normalization.
      // (Fortran WAVELJ: zeroing only after final ALPHAR/ALPHAI computation)
      double inv_mag = 1.0 / mag;
      for (int j = ISTRT; j <= i + 1; ++j) u[j] *= inv_mag;
      // Track ISTRT: find first nonzero (inner region is naturally small)
      double threshold = STEPI * inv_mag;  // already-rescaled threshold
      for (int j = ISTRT; j <= i + 1; ++j) {
        if (std::abs(u[j]) >= threshold) { ISTRT = j; break; }
      }
    }
  }


  // --- S-matrix extraction: Wronskian matching (5-point stencil) ---
  // Match at N-3 (5-point stencil fits within [0..N])
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

  // Print S-matrix
  {
    double mag_S = std::sqrt(SJR*SJR + SJI*SJI);
    double phase_S = std::atan2(SJI, SJR);
    bool is_in = (&ch == &Incoming);
    fprintf(stderr, "SMAT_%s L=%2d JP=%2d |S|=%.6f phase=%.4f SJR=%.6f SJI=%.6f\n",
      is_in ? "IN" : "OUT", L, Jp, mag_S, phase_S, SJR, SJI);
  }

  // --- Normalization (Ptolemy source.mor lines 30913–30916) ---
  // Ptolemy normalizes at NSTEP using F, G at R_max.
  // ALPHAR = (u_N·A1n + v_N·A2n) / (u_N² + v_N²)
  // ALPHAI = (u_N·A2n - v_N·A1n) / (u_N² + v_N²)
  // ur, ui already defined from u_m = u[n_match]
  double A1n  = 0.5 * (FL * (1.0 + SJR) + SJI * GL);
  double A2n  = 0.5 * (GL * (1.0 - SJR) + SJI * FL);
  double den2 = ur * ur + ui * ui;  // |u[NSTEP]|²

  if (den2 > 1e-60) {
    double alpha_r = (ur * A1n + ui * A2n) / den2;
    double alpha_i = (ur * A2n - ui * A1n) / den2;
    std::complex<double> alpha(alpha_r, alpha_i);

    double AA1 = std::abs(alpha_r) + std::abs(alpha_i);
    // Fortran DO 629: after normalization, zero points where |u_norm| < STEPI
    // BB1 = STEPI / AA1 — threshold for the UNnormalized u array
    double BB1 = (AA1 > 1e-300) ? STEPI / AA1 : 1e300;

    for (int i = 0; i <= N; ++i) {
      double ui_mag = std::abs(u[i].real()) + std::abs(u[i].imag());
      if (ui_mag < BB1) {
        ch.WaveFunction[i] = {0.0, 0.0};
      } else {
        ch.WaveFunction[i] = u[i] * alpha;
      }
    }
  }

  // --- Coulomb asymptotic extension beyond NSTEP (Ptolemy WAVELJ NSTP2S logic) ---
  // Ptolemy extends chi to NSTP2S = RIMAX/h+3.5 (RIMAX can exceed MaxR=30fm).
  // Integration domain reaches RI,RO up to SUMMAX=30.4fm > MaxR.
  // Without extension, 5-pt Lagrange stencil reads guard zeros → large error.
  // Here we extend with the exact Coulomb asymptotic:
  //   chi_b(r) = 0.5*[(1+S)*F_L + i*(1-S)*G_L] / (normalized to our convention)
  // Since we already know alpha = normalization constant, and the Coulomb functions,
  // just call Rcwfn at each extension point and apply alpha.
  {
    // Fortran NSTP2S = RIMAX/h + 3.5 → INT = RIMAX/h + 3 (for h=0.125, RIMAX=30: NSTP2S=243)
    // Extra points = NSTP2S - NSTEPS = (RIMAX/h+3) - RIMAX/h = 3.
    // This gives max chi R = NSTP2S * h = (RIMAX/h+3)*h = RIMAX+3h = 30.375 fm (→ Fortran "30.4").
    const int NEXTRA = static_cast<int>(ch.MaxR / h + 3.5) - static_cast<int>(ch.MaxR / h);
    int old_size = static_cast<int>(ch.WaveFunction.size());
    int need_size = N + NEXTRA + 4;
    if (need_size > old_size) ch.WaveFunction.resize(need_size, {0.0, 0.0});

    if (den2 > 1e-60) {
      double alpha_r = (ur * A1n + ui * A2n) / den2;
      double alpha_i = (ur * A2n - ui * A1n) / den2;
      std::complex<double> alpha(alpha_r, alpha_i);

      std::vector<double> FC_ext(L + 2, 0.0), FCP_ext(L + 2, 0.0);
      std::vector<double> GC_ext(L + 2, 0.0), GCP_ext(L + 2, 0.0);

      for (int I = N + 1; I <= N + NEXTRA; ++I) {
        double r_ext   = I * h;
        double rho_ext = ch.k * r_ext;
        Rcwfn(rho_ext, ch.eta, L, L, FC_ext, FCP_ext, GC_ext, GCP_ext);
        double FL_ext = FC_ext[L];
        double GL_ext = GC_ext[L];
        // Unnormalized chi at r_ext: same asymptotic as at matching point
        // u_unnorm = 0.5*((1+S)*F + i*(1-S)*G) at large r
        double A1_ext = 0.5 * (FL_ext * (1.0 + SJR) + SJI * GL_ext);
        double A2_ext = 0.5 * (GL_ext * (1.0 - SJR) + SJI * FL_ext);
        // Normalized: chi = (A1+iA2) * alpha
        std::complex<double> chi_ext(A1_ext, A2_ext);
        ch.WaveFunction[I] = chi_ext * alpha;
      }
    }
  }

}
