// ineldc.cpp — DWBA::GrdSet() and DWBA::InelDc() (transfer integral)
// Extracted from dwba.cpp lines 568-1697
// File-local CubMap helper reproduced here (static, also in dwba.cpp).

#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <tuple>

// Constants
const double HBARC = 197.32697; // MeV fm
const double AMU = 931.494;     // MeV/c^2
const double FINE_STRUCTURE = 1.0 / 137.035999;

// ---------------------------------------------------------------

// CubMap: Ptolemy CUBMAP subroutine ported to C++
// Maps GL points on [-1,1] into [xlo,xhi] with concentration near xmid.
// MAPTYP: 0=linear, 1=cubic-sinh, 2=rational-sinh, 3=linear-sinh
// gamma=0 or GAMPHI~1e-6 → nearly linear (rational-sinh degenerates)
// Called with GL points already in args[], wts[].
// ---------------------------------------------------------------
static void CubMap(int maptyp, double xlo, double xmid, double xhi, double gamma,
                   std::vector<double>& args, std::vector<double>& wts) {
  int npts = (int)args.size();
  // arcsinh(gamma)
  double tau = (gamma > 1e-6) ? std::log(gamma + std::sqrt(gamma*gamma + 1.0))
                                : gamma * (1.0 - gamma*gamma/6.0);
  double xlen = xhi - xlo;
  double xadd = xlo + xhi;

  if (maptyp == 0) {
    // Linear map to [xlo,xhi]
    for (int i = 0; i < npts; i++) {
      args[i] = xlo + 0.5*xlen*(args[i] + 1.0);
      wts[i] *= 0.5*xlen;
    }
  } else if (maptyp == 1) {
    // Cubic-sinh map, piles points near xmid
    xmid = std::max(xmid, xlo + xlen/7.0);
    xmid = std::min(xmid, 0.5*xadd);
    double A = 0.5*xadd - xmid;
    double B = 0.5*xlen;
    double C = 0.5*xadd;
    for (int i = 0; i < npts; i++) {
      double tu = tau * args[i];
      double xi = std::sinh(tu) / gamma;
      args[i] = A*(xi*xi - 1.0)*(xi + 1.0) + B*xi + C;
      wts[i] *= (tau/gamma) * std::cosh(tu) * ((3.0*xi - 1.0)*(xi + 1.0)*A + B);
    }
  } else if (maptyp == 2) {
    // Rational-sinh map (Ptolemy default for sum/phi)
    // For XMID=0.5*(XLO+XHI), this reduces to linear map.
    // For GAMPHI→0 (tau→gamma), denom→XLEN, → linear map.
    double A = -xmid * xlen;
    double B = xlen;
    double C = xmid*xadd - 2.0*xlo*xhi;
    double D = xadd - 2.0*xmid;
    for (int i = 0; i < npts; i++) {
      double tu = tau * args[i];
      double sh = std::sinh(tu);
      double denom = B - (D/gamma)*sh;
      args[i] = (-A + (C/gamma)*sh) / denom;
      wts[i] *= (tau/gamma) * std::cosh(tu) * ((B*C - A*D) / (denom*denom));
    }
  } else if (maptyp == 3) {
    // Linear-sinh map
    for (int i = 0; i < npts; i++) {
      double tu = tau * args[i];
      args[i] = xlo + 0.5*xlen * (std::sinh(tu)/gamma + 1.0);
      wts[i] *= 0.5*xlen * (tau/gamma) * std::cosh(tu);
    }
  }
}


void DWBA::InelDc() {
  // ===========================================================
  // Finite-Range DWBA integration
  // Produces TransferSMatrix: S[Lx][{Li,Lo}] = radial integral
  // projected onto angular momentum transfer Lx
  //
  // The radial integral is:
  //   I(Lx, Li, Lo) = int dra int drb
  //                   chi_b*(rb) * K_Lx(ra,rb) * chi_a(ra)
  //
  // where the kernel is:
  //   K_Lx(ra,rb) = int_{-1}^{1} d(cos_theta)
  //                 phi_T(rx) * V(rp) * phi_P(rp) * P_Lx(cos_theta)
  // ===========================================================

  // Masses
  double mA = Incoming.Target.Mass;
  double ma = Incoming.Projectile.Mass;
  double mb = Outgoing.Projectile.Mass;
  double mB = Outgoing.Target.Mass;

  const double AMU_MEV = 931.494061;

  // True mass of transferred particle (recover neutron mass from binding)
  double mx_kinematic = ma - mb;
  double mx = mx_kinematic + ProjectileBS.BindingEnergy / AMU_MEV;

  // Coordinate transformation coefficients (Ptolemy GRDSET convention)
  // r_x = S1*ra + T1*rb  (coordinate of x in the target frame)
  // r_p = S2*ra + T2*rb  (coordinate of x in the projectile frame)
  double ratio_xA = mx / mA;  // BRATMS(2) = x/BIGA
  double ratio_xb = mx / mb;  // BRATMS(1) = x/b
  // Ptolemy GRDSET (source.mor lines 15875-15901) for STRIPPING (ISTRIP=+1):
  // BRATMS(1) = x/b = ratio_xb, BRATMS(2) = x/BIGA = ratio_xA
  // TEMP = 1/(BRATMS(1) + BRATMS(2)*(1+BRATMS(1)))
  // S1 = (1+BRATMS(1))*(1+BRATMS(2))*TEMP
  // T1 = -(1+BRATMS(2))*TEMP
  // S2 = (1+BRATMS(1))*TEMP
  // T2 = -S1
  double denom = ratio_xb + ratio_xA * (1.0 + ratio_xb);  // CORRECT: BRATMS(1) + BRATMS(2)*(1+BRATMS(1))
  double S1 = (1.0 + ratio_xb) * (1.0 + ratio_xA) / denom;  // ≈ 1.940
  double T1 = -(1.0 + ratio_xA) / denom;                     // ≈ -0.970 (CORRECT)
  double S2 = (1.0 + ratio_xb) / denom;                      // ≈ 1.883 (CORRECT)
  double T2 = -S1;                                            // ≈ -1.940 (T2 = -S1)
  // JACOB = S1^3 confirmed by standalone Fortran test (test_rirowts.f).
  // Ptolemy GRDSET line 15882: JACOB = S1**3 for stripping.
  // No compensating /JACOB found anywhere in source.mor.
  // See AGENT_FINDINGS.md section "GRDSET RIROWTS Standalone Fortran Validation".
  double JACOB_grdset = S1*S1*S1;  // Ptolemy GRDSET: JACOB = S1^3 for stripping

  // -------------------------------------------------------
  // Build bound state channels and solve for wave functions
  // -------------------------------------------------------

  // Target Bound State: neutron x orbiting target core A
  Channel TgtBS_ch;
  TgtBS_ch.Pot    = TargetBS.Pot;
  TgtBS_ch.Target = Incoming.Target;
  TgtBS_ch.Projectile.Z    = Incoming.Projectile.Z - Outgoing.Projectile.Z;
  TgtBS_ch.Projectile.A    = Incoming.Projectile.A - Outgoing.Projectile.A;
  TgtBS_ch.Projectile.Mass = mx;
  TgtBS_ch.mu   = mA * mx / (mA + mx) / AMU_MEV;   // AMU
  TgtBS_ch.StepSize = Incoming.StepSize;
  TgtBS_ch.MaxR     = Incoming.MaxR;
  TgtBS_ch.NSteps   = Incoming.NSteps;
  TgtBS_ch.RGrid    = Incoming.RGrid;
  TgtBS_ch.WaveFunction.resize(Incoming.NSteps);
  TgtBS_ch.V_real.resize(Incoming.NSteps);
  TgtBS_ch.V_imag.resize(Incoming.NSteps);
  TgtBS_ch.V_so_real.resize(Incoming.NSteps);
  TgtBS_ch.V_so_imag.resize(Incoming.NSteps);
  TgtBS_ch.V_coulomb.resize(Incoming.NSteps);

  WavSet(TgtBS_ch);
  CalculateBoundState(TgtBS_ch, TargetBS.n, TargetBS.l, TargetBS.j, TargetBS.BindingEnergy);

  // Projectile Bound State: neutron x orbiting ejectile b (proton)
  Channel PrjBS_ch;
  PrjBS_ch.Pot    = ProjectileBS.Pot;
  PrjBS_ch.Target = Outgoing.Projectile;
  PrjBS_ch.Projectile.Z    = Incoming.Projectile.Z - Outgoing.Projectile.Z;
  PrjBS_ch.Projectile.A    = Incoming.Projectile.A - Outgoing.Projectile.A;
  PrjBS_ch.Projectile.Mass = mx;
  PrjBS_ch.mu   = mb * mx / (mb + mx) / AMU_MEV;   // AMU
  PrjBS_ch.StepSize = Incoming.StepSize;
  PrjBS_ch.MaxR     = Incoming.MaxR;
  PrjBS_ch.NSteps   = Incoming.NSteps;
  PrjBS_ch.RGrid    = Incoming.RGrid;
  PrjBS_ch.WaveFunction.resize(Incoming.NSteps);
  PrjBS_ch.V_real.resize(Incoming.NSteps);
  PrjBS_ch.V_imag.resize(Incoming.NSteps);
  PrjBS_ch.V_so_real.resize(Incoming.NSteps);
  PrjBS_ch.V_so_imag.resize(Incoming.NSteps);
  PrjBS_ch.V_coulomb.resize(Incoming.NSteps);

  WavSet(PrjBS_ch);

  // Projectile bound state: for the deuteron (n+p), use the tabulated Reid soft-core
  // S-state wavefunction from reid-phi-v (same data directory as mass table).
  // Ptolemy uses WAVEFUNCTION=REID linkule — it never runs BOUND for the projectile.
  // CalculateBoundState fails for n+p (matching radius too small for A_core=1).
  // For other projectile BSs (not deuteron), fall back to CalculateBoundState.
  bool reidLoaded = false;
  if (Incoming.Projectile.A == 2 && Incoming.Projectile.Z == 1 &&
      ProjectileBS.l == 0) {
    // Deuteron S-state (L=0): load AV18 wavefunction (same as Ptolemy "wavefunction av18")
    // AV18 is the reference for o16dp_gs_new.in; Reid is available as "reid-phi-v"
    reidLoaded = LoadDeuteronWavefunction(PrjBS_ch,
        "/home/node/working/ptolemy_2019/Cpp_AI/data", "av18-phi-v");
  }
  if (!reidLoaded) {
    // Fallback: solve WS bound state via Numerov (works for heavy cores)
    CalculateBoundState(PrjBS_ch, ProjectileBS.n, ProjectileBS.l,
                        ProjectileBS.j, ProjectileBS.BindingEnergy);
  }

  // Rebuild PrjBS_ch.V_real using the solved potential depth V_sol (set by CalculateBoundState).
  // WavSet filled V_real with the initial pot.V; after bound state solve, pot.V = V_sol.
  // We need the correct V_np(r) for the POST form vertex.
  {
    for (int i = 0; i < PrjBS_ch.NSteps; ++i) {
      double r = i * PrjBS_ch.StepSize;
      if (r < 0.001) { PrjBS_ch.V_real[i] = 0.0; continue; }
      EvaluatePotential(r, PrjBS_ch.Pot, PrjBS_ch.V_real[i], PrjBS_ch.V_imag[i],
                        PrjBS_ch.V_so_real[i], PrjBS_ch.V_so_imag[i],
                        PrjBS_ch.V_coulomb[i], PrjBS_ch.Projectile.Z,
                        PrjBS_ch.Target.Z, PrjBS_ch.Target.A, PrjBS_ch.Projectile.A);
    }
  }

  // Override projectile WF with tabulated values (e.g. AV18) if loaded
  if (ProjectileWFLoaded) {
    double h_prj = PrjBS_ch.StepSize;
    double h_tab = ProjectileWFGridH;
    for (int I = 0; I < PrjBS_ch.NSteps; I++) {
      double r = I * h_prj;
      // Interpolate tabulated WF onto the calculation grid
      double idx_f = r / h_tab;
      int    idx_i = static_cast<int>(idx_f);
      double frac  = idx_f - idx_i;
      double phi   = 0.0;
      if (idx_i + 1 < (int)ProjectileWFTable.size()) {
        phi = ProjectileWFTable[idx_i] * (1.0 - frac)
            + ProjectileWFTable[idx_i+1] * frac;
      } else if (idx_i < (int)ProjectileWFTable.size()) {
        phi = ProjectileWFTable[idx_i];
      }
      // Do NOT multiply by SPAM here — SPAM is already in ATERM
      PrjBS_ch.WaveFunction[I] = std::complex<double>(phi, 0.0);
    }
  }

  // Sanity check and find cutoff radius for bound states
  // Beyond this radius phi < 1e-6 * peak, contribution is negligible
  double maxT = 0, maxP = 0;
  int cutT = TgtBS_ch.NSteps - 1, cutP = PrjBS_ch.NSteps - 1;
  for (int i = 1; i < TgtBS_ch.NSteps; ++i)
    maxT = std::max(maxT, std::abs(TgtBS_ch.WaveFunction[i].real()));
  for (int i = 1; i < PrjBS_ch.NSteps; ++i)
    maxP = std::max(maxP, std::abs(PrjBS_ch.WaveFunction[i].real()));
  for (int i = TgtBS_ch.NSteps - 1; i > 0; --i) {
    if (std::abs(TgtBS_ch.WaveFunction[i].real()) > 1e-5 * maxT) { cutT = i + 5; break; }
  }
  for (int i = PrjBS_ch.NSteps - 1; i > 0; --i) {
    if (std::abs(PrjBS_ch.WaveFunction[i].real()) > 1e-5 * maxP) { cutP = i + 5; break; }
  }
  cutT = std::min(cutT, TgtBS_ch.NSteps - 1);
  cutP = std::min(cutP, PrjBS_ch.NSteps - 1);

  // -------------------------------------------------------
  // Determine allowed Lx (angular momentum transfer)
  // Triangle rule: |lT - lP| <= Lx <= lT + lP
  // -------------------------------------------------------
  int lT = TargetBS.l;
  int lP = ProjectileBS.l;
  int LxMin = std::abs(lT - lP);
  int LxMax_bs = lT + lP;

  // -------------------------------------------------------
  // Partial wave loops
  // -------------------------------------------------------
  int Lmax = 15;  // Increased; Lmax=20 needed for full accuracy

  TransferSMatrix.clear();



  // Precompute the bound state potential on the grid.
  // Vertex potential grid. Currently uses target BS potential V_nA at rx.
  // Ptolemy BSPROD with IVRTEX=1 uses (V_np*phi_P)(rp) but with Ryan's T1/S2 values,
  // rp is large (~5 fm at integrand peak) so V_np(rp)≈0 — gives wrong result (3e-3 mb/sr).
  // The correct Ptolemy behavior for this reaction needs further investigation.
  // See AGENT_FINDINGS.md for full analysis of the 2.33x discrepancy.
  std::vector<double> V_bx_grid(TgtBS_ch.NSteps, 0.0);
  for (int i = 1; i < TgtBS_ch.NSteps; ++i) {
    V_bx_grid[i] = TgtBS_ch.V_real[i]; // target BS potential V_nA
  }

  // Ptolemy BSPROD ITYPE=2 clipping of (V*phi) products.
  // For PRIOR form: clips (V_nA * phi_T)(rx) — target vertex product.
  // For POST form:  clips (V_np * phi_P)(rp) — projectile vertex product.
  // Clipping rule: f'(r) = f(r) for r >= r_peak; f'(r) = f_max for r < r_peak.
  // Where r_peak = location of max(|f(r)|) and f_max = max(|f(r)|).
  //
  // For PRIOR form: the vertex array is V_nA stored in TgtBS_ch.V_real.
  //   IVPHI_T[i] = TgtBS_ch.WaveFunction[i] * TgtBS_ch.V_real[i] = phi_T * V_nA
  // For POST form: the vertex array is V_np stored in PrjBS_ch.V_real.
  //   IVPHI_P[i] = PrjBS_ch.WaveFunction[i] * PrjBS_ch.V_real[i] = phi_P * V_np
  // The "other" WF (not at vertex) is NOT clipped in the PRIOR/POST vertex product,
  // but IS clipped separately: phi'_T or phi'_P.
  //
  // Build IVPHI arrays:
  int NSteps_common = TgtBS_ch.NSteps;  // = PrjBS_ch.NSteps = 301
  double h_common = TgtBS_ch.StepSize;  // = 0.1 fm
  // USECORE correction parameters for POST form (NUCONL=3, stripping, IVRTEX=1):
  // For stripping (d,p): ISC=2 (outgoing channel), IOUTSW=TRUE
  //   RNSCAT = R0_out * mB^(1/3)   [outgoing p-B scattering radius]
  //   RNCORE = RNSCAT * (ma^(1/3) + mb^(1/3)) / (mB^(1/3) + mb^(1/3))
  //   VOPT = -V0_out (negative of outgoing real depth)
  //   AOPT = A_out  (outgoing diffuseness)
  // Core correction: V_core(r_core) = VOPT / (1 + exp((r_core - RNCORE) / AOPT))
  // where r_core is the core-core separation, evaluated at r_core = (RNCORE/RNSCAT) * rb
  // This is added per (ra,rb) in the inner loop for POST form.
  // Use integer mass numbers (A) for radius scaling, not masses in MeV
  int A_A = Incoming.Target.A;       // 33Si: A=33  (AMBGA = target A)
  int A_b = Outgoing.Projectile.A;   // proton: A=1  (AMB = ejectile b)
  int A_B = Outgoing.Target.A;       // 34Si: A=34  (AMBGB = residual B)
  // BSPROD stripping/projectile-vertex: ISC=2 (outgoing channel)
  // Ptolemy BSSET (source.mor ~4829):
  //   RNSCAT = RSCTS(2) = R0_out * A_B^(1/3)
  //   RNCORE = RNSCAT * (AMBGA3 + AMB3) / (AMBGB3 + AMB3)
  //         = RNSCAT * (A_A^(1/3) + A_b^(1/3)) / (A_B^(1/3) + A_b^(1/3))
  // BUG FIX: previously used A_a (projectile=deuteron=2) instead of A_A (target=33Si=33)
  double RNSCAT_post = Outgoing.Pot.R0 * std::pow((double)A_B, 1.0/3.0);
  double RNCORE_post = RNSCAT_post * (std::pow((double)A_A, 1.0/3.0) + std::pow((double)A_b, 1.0/3.0))
                                   / (std::pow((double)A_B, 1.0/3.0) + std::pow((double)A_b, 1.0/3.0));
  double VOPT_post   = -Outgoing.Pot.V;  // negative of outgoing real depth (MeV)
  double AOPT_post   = Outgoing.Pot.A;
  double CORE_SCALE  = RNCORE_post / RNSCAT_post;  // r_core = CORE_SCALE * rb
  auto V_core_post = [&](double r_core) -> double {
    if (r_core < 1e-6) return VOPT_post;
    return VOPT_post / (1.0 + std::exp((r_core - RNCORE_post) / AOPT_post));
  };

  std::vector<double> IVPHI_T(NSteps_common, 0.0);  // phi_T * V_nA (PRIOR vertex)
  std::vector<double> IVPHI_P(NSteps_common, 0.0);  // phi_P * V_np (POST vertex)
  double IVPHI_T_max = 0, IVPHI_P_max = 0;
  int idx_T_vert_peak = 1, idx_P_vert_peak = 1;

  // Vertex potential:
  //   Target vertex: phi_T(r) * V_WS_nA(r) — always WS bound-state potential ✅
  //   Projectile vertex: phi_P(r) * V_WS_np(r)
  //     WAVEFUNCTION=REID replaces only the bound-state solver for phi_P(r).
  //     The vertex potential V_np is still the WS from the PROJECTILE input section
  //     (V=60, R=1.0, A=0.5). Block 1 of reid-phi-v is NOT used by Ptolemy.
  //     V_real[] is already set to this WS by WavSet/EvaluatePotential. ✅
  std::cout << "  Vertex: V_WS(V=60,R=1,A=0.5)*phi_Reid (projectile), V_WS*phi_T (target)" << std::endl;

  for (int i = 1; i < NSteps_common; ++i) {
    IVPHI_T[i] = std::abs(TgtBS_ch.WaveFunction[i].real()) * TgtBS_ch.V_real[i];
    // Use VPhiProduct (Reid V_eff * phi_S) if available; else fall back to WS V * phi
    IVPHI_P[i] = !PrjBS_ch.VPhiProduct.empty()
                 ? PrjBS_ch.VPhiProduct[i]
                 : std::abs(PrjBS_ch.WaveFunction[i].real()) * PrjBS_ch.V_real[i];
    if (IVPHI_T[i] > IVPHI_T_max) { IVPHI_T_max = IVPHI_T[i]; idx_T_vert_peak = i; }
    if (IVPHI_P[i] > IVPHI_P_max) { IVPHI_P_max = IVPHI_P[i]; idx_P_vert_peak = i; }
  }
  double r_T_vert_peak = idx_T_vert_peak * h_common;
  double r_P_vert_peak = idx_P_vert_peak * h_common;

  // Report vertex peaks for diagnostic comparison vs Ptolemy BSSET output
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  IVPHI_T peak = " << IVPHI_T_max << " at r=" << r_T_vert_peak << " fm" << std::endl;
  std::cout << "  IVPHI_P peak = " << IVPHI_P_max << " at r=" << r_P_vert_peak << " fm" << std::endl;

  // DIAGNOSTIC: print IVPHI_P and phi_P at sample rP values vs Fortran reference
  {
    double rP_vals[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    double fort_ref[] = {-174.148, 45.868, 17.487, 3.697, 0.844, 0.245, 0.090, 0.039};
    fprintf(stderr, "\n# IVPHI_P comparison: C++ vs Fortran (Reid cubic+norm, V_Reid)\n");
    fprintf(stderr, "# %6s  %14s  %14s  %8s  %14s\n",
            "rP(fm)", "C++_IVPHI_P", "Fort_ref", "Ratio", "C++_phi_P");
    for (int k = 0; k < 8; ++k) {
      double rP = rP_vals[k];
      int idx = (int)std::round(rP / h_common);
      if (idx >= NSteps_common) continue;
      double cpp_ivphi = IVPHI_P[idx];
      double cpp_phi   = std::abs(PrjBS_ch.WaveFunction[idx].real());
      double fort      = fort_ref[k];
      double ratio     = (std::abs(fort) > 1e-8) ? cpp_ivphi / fort : 0.0;
      fprintf(stderr, "  %6.2f  %14.5e  %14.5e  %8.4f  %14.5e\n",
              rP, cpp_ivphi, fort, ratio, cpp_phi);
    }
    fprintf(stderr, "\n");
  }

  // The "other" WF (non-vertex side): clipped separately
  auto findWFPeak = [&](const std::vector<std::complex<double>> &wf, int nsteps)
      -> std::pair<double, int> {
    double phi_max = 0; int idx_peak = 1;
    for (int i = 1; i < nsteps; ++i) {
      double v = std::abs(wf[i].real());
      if (v > phi_max) { phi_max = v; idx_peak = i; }
    }
    return {phi_max, idx_peak};
  };
  auto [phi_T_max, idx_T_peak] = findWFPeak(TgtBS_ch.WaveFunction, NSteps_common);
  auto [phi_P_max, idx_P_peak] = findWFPeak(PrjBS_ch.WaveFunction, NSteps_common);
  double r_T_peak = idx_T_peak * h_common;
  double r_P_peak = idx_P_peak * h_common;

  // Interpolation helper (linear, returns real WF value)
  auto Interpolate = [&](const std::vector<std::complex<double>> &wf,
                         const std::vector<double> &grid,
                         int nsteps, double stepsize, double maxr,
                         double r) -> double {
    if (r >= maxr || r < 0) return 0.0;
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return wf[ii].real() * (1.0 - frac) + wf[ii + 1].real() * frac;
  };

  // Clipped interpolation for plain phi (non-vertex side):
  // phi'(r) = phi_max for r < r_peak; phi(r) for r >= r_peak
  // Ptolemy Pass 2 uses ITYPE=1 (PHI V PHI), NO clipping.
  // The bound state WFs are phi=u/r, always positive for 0-node states.
  auto InterpolateClipped = [&](const std::vector<std::complex<double>> &wf,
                                const std::vector<double> &/*grid*/,
                                int nsteps, double stepsize, double maxr,
                                double r, double /*phi_max*/, double /*r_peak*/) -> double {
    if (r >= maxr || r < 0) return 0.0;
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return wf[ii].real() * (1.0 - frac) + wf[ii + 1].real() * frac;
  };

  // Clipped IVPHI interpolation for vertex product (V*phi at vertex side):
  // IVPHI'(r) = IVPHI_max for r < r_vert_peak; IVPHI(r) for r >= r_vert_peak
  // InterpolateIVPHI: interpolate V(r)*phi(r) product WITHOUT clipping
  // Ptolemy Pass 2 uses ITYPE=1 (PHI V PHI), no clipping.
  auto InterpolateIVPHI = [&](const std::vector<double> &ivphi,
                               double stepsize, int nsteps, double maxr,
                               double r, double /*ivphi_max*/, double /*r_vert_peak*/) -> double {
    if (r >= maxr || r < 0) return 0.0;
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return ivphi[ii] * (1.0 - frac) + ivphi[ii + 1] * frac;
  };

  auto InterpolateV = [&](const std::vector<double> &v,
                          double stepsize, int nsteps, double maxr,
                          double r) -> double {
    if (r >= maxr || r < 0) return 0.0;
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return v[ii] * (1.0 - frac) + v[ii + 1] * frac;
  };

  // Determine effective grid cutoff for distorted waves contribution.
  // The integrand is negligible where BOTH bound states are ~0.
  // The max extent is driven by the larger of the two bound state cutoffs.
  // For the outer integral (ra), the chi_a(ra) that matters is up to ~cutT/S1 or ~cutP/S2
  // Use full grid (30 fm) to capture the long tail of the bound states
  // (Target BS cutoff can be ~22 fm — must not truncate!)
  int N_int = Incoming.NSteps; // full 30 fm range

  bool debugPrinted = false;

  // Spin of projectile/ejectile (doubled integers) for J-split DW
  const int JA_dw = 2;  // 2 * spin_deuteron = 2*1 = 2
  const int JB_dw = 1;  // 2 * spin_proton   = 2*(1/2) = 1

  for (int Li = 0; Li <= Lmax; ++Li) {
    // --- J-split incoming distorted waves for all JPI ---
    // JPI ranges from |2*Li - JA_dw| to 2*Li + JA_dw in steps of 2
    int JPI_min = std::abs(2*Li - JA_dw);
    int JPI_max = 2*Li + JA_dw;
    std::map<int, std::vector<std::complex<double>>> chi_a_byJPI;
    for (int JPI = JPI_min; JPI <= JPI_max; JPI += 2) {
      WavElj(Incoming, Li, JPI);
      chi_a_byJPI[JPI] = Incoming.WaveFunction;
    }

    // For each outgoing partial wave Lo, compute ALL J-split JPO values and integrals.
    // Lo ranges over all values allowed by parity and triangle conditions.
    // JPO ranges from |2*Lo - JB_dw| to 2*Lo + JB_dw in steps of 2 (J-split).

    for (int Lo = 0; Lo <= Lmax; Lo++) {
      // Pre-compute all J-split outgoing waves for this Lo
      std::map<int, std::vector<std::complex<double>>> chi_b_byJPO;
      int JPO_min = std::abs(2*Lo - JB_dw);
      int JPO_max = 2*Lo + JB_dw;
      for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
        if (JPO < 1) continue;  // must be positive
        WavElj(Outgoing, Lo, JPO);
        chi_b_byJPO[JPO] = Outgoing.WaveFunction;
      }
      
      // For each allowed Lx, compute the radial integral — separate per (JPI, JPO)
      for (int Lx = LxMin; Lx <= LxMax_bs; Lx += 2) {
        // Triangle rule check: |Li - Lo| <= Lx <= Li + Lo
        if (Lx < std::abs(Li - Lo)) continue;
        if (Lx > Li + Lo) continue;
        // Parity check: (Li + Lo + Lx) must be even
        if ((Li + Lo + Lx) % 2 != 0) continue;

      // ---- Precompute A12 angular coupling coefficients for (Li, Lo, Lx) ----
      // Extracted to DWBA::ComputeA12Terms() in a12.cpp
      std::vector<std::tuple<int,int,double>> A12_terms = ComputeA12Terms(Li, Lo, Lx, lT, lP);

      // J-split loop: compute I(Li,Lo,Lx,JPI,JPO) for each (JPI,JPO) pair
      for (auto &[JPI, chi_a] : chi_a_byJPI) {
        for (auto &[JPO, chi_b] : chi_b_byJPO) {

        // ---------------------------------------------------------------
        // GL quadrature in (U,V) = (sum, dif) coordinates — Ptolemy GRDSET
        // U = (ri+ro)/2, V = ri-ro
        // ri = U + V/2, ro = U - V/2
        // Domain: triangular — U in [0, SUMMAX], V in [-2U, +2U]
        // Jacobian d(ri)d(ro) = d(U)d(V): |∂(ri,ro)/∂(U,V)| = 1 (exact)
        // ---------------------------------------------------------------
        std::complex<double> Integral(0.0, 0.0);

        // Chi interpolation at arbitrary r: linear interpolation on uniform grid
        // Defined here so it's accessible in both USE_ZR and normal (#else) paths.
        double h_interp = Incoming.StepSize;
        auto interp_chi_b = [&](double r) -> std::complex<double> {
          if (r <= 0 || r >= Outgoing.MaxR) return {0.0, 0.0};
          double idx_f = r / h_interp;
          int ii = static_cast<int>(idx_f);
          if (ii >= (int)chi_b.size() - 1) return {0.0, 0.0};
          double frac = idx_f - ii;
          return chi_b[ii] * (1.0 - frac) + chi_b[ii+1] * frac;
        };

#ifdef USE_ZR
        // ---------------------------------------------------------------
        // ZERO-RANGE (ZR) approximation — PHYSICALLY CORRECT VERSION
        //
        // In ZR, V_np(rp)*phi_P(rp) = D0 * delta^3(rp).
        // The delta function forces rp=0, which (at phi_ab=0) means:
        //   ra/rb = S1_c/S2_c  (where _c = CORRECT Ptolemy values)
        //   rx = ra  (neutron coordinate = deuteron coordinate exactly)
        //
        // This collapses the 3D integral to 1D over ra:
        //   I_ZR(Li,Lo,Lx) = D0 * A12(phi=0) *
        //       ∫_0^∞ chi_a(ra) * phi_T(ra) * chi_b*(rb=zr_scale*ra)
        //              * ra * rb * J_zr * dra
        //
        // Coordinate geometry (Ptolemy GRDSET BRATMS convention):
        //   BRATMS(1) = mx/mb (= x/b, neutron/proton ≈ 1.001)
        //   BRATMS(2) = mx/mA (= x/A, neutron/33Si ≈ 0.0306)
        //   S1_c = (1+BRATMS1)*(1+BRATMS2)/denom_c   ≈ 1.941
        //   T1_c = -(1+BRATMS2)/denom_c               ≈ -0.970
        //   S2_c = (1+BRATMS1)/denom_c                ≈ 1.883
        //   T2_c = -S1_c                              ≈ -1.941
        //   denom_c = BRATMS1 + BRATMS2*(1+BRATMS1)
        //
        // ZR geometry: phi_ab=0, rp=0 → rb = (S2_c/S1_c)*ra ≈ (A/(A+1))*ra ≈ 0.970*ra
        //                                  rx = ra (exact: S1_c + T1_c*(S2_c/S1_c) = 1)
        //
        // No Li==Lo restriction in ZR! The phi_ab=0 limit allows all Li,Lo pairs.
        // At phi_ab=0: phi_T_angle = acos((T1_c*rb + S1_c*ra)/rx)
        //   = acos((T1_c*zr_scale + S1_c)*ra / ra) = acos(S1_c + T1_c*zr_scale)
        //   = acos(1.0) = 0  (since rx/ra = 1 exactly)
        //
        // D0 = -120.1 MeV·fm^{3/2} (ZR constant, from AGENT_FINDINGS.md)
        // ---------------------------------------------------------------
        {
          const double D0_ZR = -120.1;  // MeV·fm^{3/2}

          // Compute CORRECT S1_c/S2_c/T1_c for ZR geometry
          // Note: current code has T1/S2 SWAPPED. Use correct Ptolemy convention here.
          double BRATMS1 = mx / mb;           // x/b = neutron/proton
          double BRATMS2 = mx / mA;           // x/A = neutron/33Si
          double denom_c  = BRATMS1 + BRATMS2*(1.0 + BRATMS1);
          double S1_c = (1.0 + BRATMS1)*(1.0 + BRATMS2)/denom_c;
          double T1_c = -(1.0 + BRATMS2)/denom_c;
          double S2_c = (1.0 + BRATMS1)/denom_c;
          // ZR scale: at phi=0, rp=0 → rb = (S2_c/S1_c)*ra
          double zr_scale = S2_c / S1_c;  // ≈ A/(A+1) ≈ 0.970

          double A12_at_zero = EvalA12(A12_terms, 0.0, 0.0);

          double h_zr    = Incoming.StepSize;
          int    N_zr    = TgtBS_ch.NSteps;
          std::complex<double> I_1D(0.0, 0.0);

          for (int i = 1; i < N_zr - 1; ++i) {
            double ra = i * h_zr;
            double rb = zr_scale * ra;

            std::complex<double> ca_r = chi_a[i];
            double phi_T_r = TgtBS_ch.WaveFunction[i].real();
            // Ptolemy uses chi_b * chi_a (NOT conj(chi_b) * chi_a): DWR=Re(chi_b*chi_a)
            std::complex<double> cb_r = interp_chi_b(rb);  // no conjugate

            I_1D += ca_r * phi_T_r * cb_r * h_zr;
          }

          Integral = D0_ZR * A12_at_zero * I_1D;
        }

#else
        double h = Incoming.StepSize;  // kept for chi interpolation

        // ---------------------------------------------------------------
        // Quadrature parameters:
        //   Sum (U): Large GL grid for chi oscillations (direct quadrature)
        //   Dif (V): Large GL grid (symmetric)
        //   Phi:     CUBMAP adaptive PHI0 per-(ra,rb) (Ptolemy MAPPHI)
        // ---------------------------------------------------------------
        const int NPSUM = 40;    // GL points for U  (dpsb: NPSUM=40)
        const int NPDIF = 40;    // GL points for V  (dpsb: NPDIF=40)
        const int NPPHI = 20;    // GL points for phi (dpsb: NPPHI=20)
        const double SUMMAX = 30.0;
        const double GAMPHI = 1.0e-6;  // → linear phi map
        const double PHIMID = 0.5;    // midpoint in [0,1]
        const double DWCUT = 1.0e-3;  // Ptolemy default DWCUTOFF
        const int LOOKST = 250;       // test points for phi scan
        const int NPHIAD = 4;         // extra phi points margin
        const double DXV = 2.0 / ((double)LOOKST * LOOKST);
        const int IXTOPZ = LOOKST + 1;

        // Build CUBMAP GL points for phi on [0,1] (scaled by PHI0 later)
        // CUBMAP(MAPPHI=2, XLO=0, XMID=0.5, XHI=1, GAMPHI=1e-6) → nearly linear
        std::vector<double> phi_pts(NPPHI), phi_wts(NPPHI);
        GaussLegendre(NPPHI, -1.0, 1.0, phi_pts, phi_wts);
        CubMap(2, 0.0, PHIMID, 1.0, GAMPHI, phi_pts, phi_wts);
        // phi_pts[i] ∈ [0,1], phi_wts[i] weights for integration over [0,1]
        // In Pass 2: PHI = PHI0 * phi_pts[i], DPHI = PHI0 * phi_wts[i] * sin(PHI)

        // Sum/dif GL on standard grids (direct quadrature)
        std::vector<double> xi_s(NPSUM), wi_s(NPSUM);
        GaussLegendre(NPSUM, -1.0, 1.0, xi_s, wi_s);

        // Chi interpolation at arbitrary r: linear interpolation on uniform grid
        auto interp_chi_a = [&](double r) -> std::complex<double> {
          if (r <= 0 || r >= Incoming.MaxR) return {0.0, 0.0};
          double idx_f = r / h;
          int ii = static_cast<int>(idx_f);
          if (ii >= (int)chi_a.size() - 1) return {0.0, 0.0};
          double frac = idx_f - ii;
          return chi_a[ii] * (1.0 - frac) + chi_a[ii+1] * frac;
        };
        // Note: interp_chi_b already defined above (before #ifdef USE_ZR) — reuse it here.

        // Helper: evaluate BSPROD-like integrand at (ra, rb, x=cos_phi)
        // Returns PRIOR: |IVPHI_T(rx)| * |phi_P(rp)|  (same as used in AngKernel)
        auto bsprod_val = [&](double ra, double rb, double x) -> double {
          double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*x;
          double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*x;
          if (rx2 < 0) rx2 = 0;
          if (rp2 < 0) rp2 = 0;
          double rx = std::sqrt(rx2);
          double rp = std::sqrt(rp2);
#ifndef USE_PRIOR_FORM
          {
            double phi_T  = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                               TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                               TgtBS_ch.MaxR, rx, phi_T_max, r_T_peak);
            double ivphi_P_nc = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                              h_common * NSteps_common, rp,
                                              IVPHI_P_max, r_P_vert_peak);
            // Ptolemy RCORE formula
            double dS2 = S1 - S2, dT2 = T1 - T2;
            double rc2 = dS2*dS2*ra*ra + dT2*dT2*rb*rb + 2.0*dS2*dT2*ra*rb*x;
            double rc = std::sqrt(std::max(rc2, 0.0));
            double fc = 1.0/(1.0+std::exp((rc - RNCORE_post)/AOPT_post));
            double fs = 1.0/(1.0+std::exp((rb - RNSCAT_post)/AOPT_post));
            double dnu = VOPT_post * (fc - fs);
            double veff = InterpolateV(PrjBS_ch.V_real, PrjBS_ch.StepSize,
                                       PrjBS_ch.NSteps, PrjBS_ch.MaxR, rp);
            double fac = (std::abs(veff) > 1e-6) ? (1.0 + dnu/veff) : 1.0;
            double ivphi_P = fac * ivphi_P_nc;
            return std::abs(phi_T * ivphi_P);
          }
#else
          double ivphi_T = InterpolateIVPHI(IVPHI_T, h_common, NSteps_common,
                                             h_common * NSteps_common, rx,
                                             IVPHI_T_max, r_T_vert_peak);
          double phi_P  = InterpolateClipped(PrjBS_ch.WaveFunction, PrjBS_ch.RGrid,
                                             PrjBS_ch.NSteps, PrjBS_ch.StepSize,
                                             PrjBS_ch.MaxR, rp, phi_P_max, r_P_peak);
          return std::abs(ivphi_T * phi_P);
#endif
        };

        // ── WVWMAX scan & SUMMIN search (Ptolemy GRDSET steps 1&2) ─────────
        // Ptolemy uses BSPROD ITYPE=3 (R*PSI*phi'*V*phi'*PSI*R) which includes
        // the distorted waves chi_a, chi_b. The chi suppression at small r
        // naturally limits WVWMAX to the physical peak and sets SUMMIN to the
        // inner integration boundary.
        //
        // Our bsprod_val doesn't include chi — we compute WVWMAX and SUMMIN
        // using a separate helper that evaluates |chi_a(U)| * bsprod_val(U,U,x) * |chi_b(U)|.
        auto bsprod_with_chi = [&](double U, double x) -> double {
          if (U < 1e-6) return 0.0;
          // chi_a at ra=U, chi_b at rb=U (same U for the diagonal scan)
          double ca_mag = std::abs(interp_chi_a(U));
          double cb_mag = std::abs(interp_chi_b(U));
          // Include ra*rb factor (Ptolemy ITYPE=3: R*PSI*...*PSI*R)
          return U * ca_mag * bsprod_val(U, U, x) * U * cb_mag;
        };

        // Step 1: find WVWMAX by scanning U from SUMMAX/2 down to 0
        double WVWMAX_pre = 0.0;
        for (double U_scan = 0.5 * SUMMAX; U_scan > 0.05; U_scan -= 0.2) {
          double v1 = bsprod_with_chi(U_scan,  1.0);
          double v2 = bsprod_with_chi(U_scan, -1.0);
          WVWMAX_pre = std::max(WVWMAX_pre, std::max(v1, v2));
        }
        // RVRLIM threshold
        double RVRLIM = DWCUT * std::max(WVWMAX_pre, 1.0e-30);

        // Step 2: find SUMMIN — first U where integrand exceeds RVRLIM
        double SUMMIN = 0.0;
        {
          double U_s = 0.0;
          while (U_s <= SUMMAX) {
            double v1 = bsprod_with_chi(U_s,  1.0);
            double v2 = bsprod_with_chi(U_s, -1.0);
            if (std::max(v1, v2) >= RVRLIM) {
              SUMMIN = std::max(0.0, U_s - 0.2);
              break;
            }
            U_s += 0.2;
          }
        }
        std::cout << std::scientific << std::setprecision(3)
                  << "  WVWMAX=" << WVWMAX_pre << "  RVRLIM=" << RVRLIM
                  << "  SUMMIN=" << SUMMIN << " fm" << std::endl;

        // Pre-allocate dif GL arrays (reused per U, resized below)
        std::vector<double> xi_d_base(NPDIF), wi_d_base(NPDIF);
        GaussLegendre(NPDIF, -1.0, 1.0, xi_d_base, wi_d_base);

        // Outer loop: U = (ri+ro)/2 in [SUMMIN, SUMMAX] from GL on [-1,1]
        for (int is = 0; is < NPSUM; ++is) {
          // Map GL point from [-1,1] onto [SUMMIN, SUMMAX]
          double U_range = SUMMAX - SUMMIN;
          double U = SUMMIN + U_range * (xi_s[is] + 1.0) / 2.0;
          double WOW = U_range * wi_s[is] / 2.0;

          if (U < 1e-6) continue;

          // Inner loop: V = ri-ro in [-2U, +2U] via GL on [-1,1]
          // Use plain GL (large NPDIF for convergence)
          for (int id = 0; id < NPDIF; ++id) {
            double V = 2.0 * U * xi_d_base[id];       // V in [-2U, +2U]
            double DIFWT = 2.0 * U * wi_d_base[id];   // dif weight

            double ra = U + V / 2.0;  // ri = U + V/2
            double rb = U - V / 2.0;  // ro = U - V/2

            // Domain constraint: ri>0, ro>0 (should always hold for |V|<=2U)
            if (ra < 1e-6 || rb < 1e-6) continue;

            // Interpolate distorted waves at non-grid (ra, rb) values
            // Ptolemy INELDC: DWR=Re(chi_b*chi_a), DWI=Im(chi_b*chi_a) — NO conjugate on chi_b
            std::complex<double> ca      = interp_chi_a(ra);
            std::complex<double> cb_conj = interp_chi_b(rb);  // no conjugate (matches Ptolemy)

            // ---------------------------------------------------------------
            // ADAPTIVE PHI0 per (ra, rb): Ptolemy two-pass approach
            //
            // Pass 1: Scan phi from 0 upward via X = 1 - DXV*(II-1)^2
            //         (X=1 → phi=0, X decreases as phi increases)
            //         Find where |BSPROD| drops below ULIM for 2 consecutive pts
            // Pass 2: Set PHI0 = acos(X0) where X0=1-DXV*(IEND-1)^2
            //         Integrate from 0 to PHI0 using NPPHI GL points
            // ---------------------------------------------------------------
            double ULIM = RVRLIM / (std::max(1.0, ra) * std::max(1.0, rb));

            // Pass 1: scan phi test points to find IEND
            int IEND = IXTOPZ;  // default: full range if never drops below cutoff
            double WOW_prev = 0.0;
            for (int II = 1; II <= IXTOPZ; ++II) {
              double X = 1.0 - DXV * (double)(II-1) * (double)(II-1);
              if (X < -1.0) X = -1.0;
              double fifo = bsprod_val(ra, rb, X);
              if (II == 1) {
                WOW_prev = fifo;
              } else {
                if (fifo < ULIM && WOW_prev < ULIM) {
                  // Both this and previous point below threshold → stop here
                  IEND = II - 1 + NPHIAD;
                  IEND = std::min(IEND, IXTOPZ);
                  IEND = std::max(IEND, 2);
                  break;
                }
                WOW_prev = fifo;
              }
            }
            // If we never broke (full range needed): IEND stays at IXTOPZ
            // IEND=IXTOPZ means integrate full [0,pi]
            IEND = std::max(IEND, 2);

            // Compute PHI0 from IEND
            double X0 = 1.0 - DXV * (double)(IEND-1) * (double)(IEND-1);
            X0 = std::max(-1.0, std::min(1.0, X0));
            double PHI0 = std::acos(X0);  // integration limit [0, PHI0]

            // Angular kernel: integrate over phi_ab from 0 to PHI0
            // using NPPHI GL points mapped by CUBMAP(MAPPHI=2, 0, 0.5, 1, 1e-6)
            // phi = PHI0 * phi_pts[k], dphi_weighted = PHI0 * phi_wts[k] * sin(phi)
            double AngKernel = 0.0;

            for (int k = 0; k < NPPHI; ++k) {
              double phi_frac = phi_pts[k];  // ∈ [0,1]
              double phi_ab = PHI0 * phi_frac;
              double phi_weight = PHI0 * phi_wts[k] * std::sin(phi_ab);  // = DPHI in Ptolemy

              double x = std::cos(phi_ab);   // cos(phi_ab)

              // Compute rx and rp magnitudes via law of cosines
              double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*x;
              double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*x;
              if (rx2 < 0) rx2 = 0;
              if (rp2 < 0) rp2 = 0;
              double rx = std::sqrt(rx2);
              double rp = std::sqrt(rp2);

              // Ptolemy BSPROD ITYPE=2 integrand:
              // PRIOR: FPFT = IVPHI_T'(rx) * phi_P'(rp)
              //   where IVPHI_T = phi_T * V_nA (clipped at its max)
              //   and phi_P' = phi_P clipped at its max
              // POST:  FPFT = phi_T'(rx) * IVPHI_P'(rp)
              //   where IVPHI_P = phi_P * V_np (clipped at its max)
              //   and phi_T' = phi_T clipped at its max
              // Then Vbx = 1.0 (already included in IVPHI)
#ifndef USE_PRIOR_FORM
              // POST form: vertex at rp (USEPROJECTILE + USECORE, Ptolemy default NUCONL=3)
              // Ptolemy BSPROD NUCONL=3 formula (source.mor ~4640-4660):
              //   FPFT_nuconly = phi_T(rx) * [V_np(rp)*phi_P(rp)]
              //   RCORE = sqrt(((S1-S2)*ra)^2 + ((T1-T2)*rb)^2 + 2*(S1-S2)*(T1-T2)*ra*rb*x)
              //   RSCAT = rb (IOUTSW=TRUE for stripping/projective vertex)
              //   DELVNU = VOPT * (WS(RCORE, RNCORE, AOPT) - WS(RSCAT, RNSCAT, AOPT))
              //   VEFF = V_np(rp)  [the vertex potential at rp]
              //   FACTOR = 1 + DELVNU / VEFF
              //   FPFT = FACTOR * FPFT_nuconly
              // POST form: phi_T at rx, vertex V_np at rp
              double phi_T  = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                                 TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                                 TgtBS_ch.MaxR, rx, phi_T_max, r_T_peak);
              // NUCONLY vertex: V_np(rp) * phi_P(rp)
              double ivphi_P_nuconly = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                                h_common * NSteps_common, rp,
                                                IVPHI_P_max, r_P_vert_peak);
              
              // USECORE: Ptolemy BSPROD NUCONL=3 formula
              // FACTOR = 1 + DELV/VEFF, where DELV = DELVNU (nuclear)
              // FPFT = FACTOR * V_np(rp) * phi_P(rp) * phi_T(rx)
              double ivphi_P;
              double VEFF = InterpolateV(PrjBS_ch.V_real, PrjBS_ch.StepSize,
                                         PrjBS_ch.NSteps, PrjBS_ch.MaxR, rp);
              if (std::abs(VEFF) > 1e-10) {
                // Ptolemy RCORE formula (BSPROD source.mor ~line 4860)
                double dS = S1 - S2;
                double dT = T1 - T2;
                double rcore2 = dS*dS*ra*ra + dT*dT*rb*rb + 2.0*dS*dT*ra*rb*x;
                if (rcore2 < 0) rcore2 = 0;
                double r_core = std::sqrt(rcore2);
                double r_scat = rb;
                double fcore = 1.0 / (1.0 + std::exp((r_core - RNCORE_post) / AOPT_post));
                double fscat = 1.0 / (1.0 + std::exp((r_scat - RNSCAT_post) / AOPT_post));
                double DELVNU = VOPT_post * (fcore - fscat);
                double FACTOR = 1.0 + DELVNU / VEFF;
                ivphi_P = FACTOR * ivphi_P_nuconly;
              } else {
                ivphi_P = ivphi_P_nuconly;
              }
              double phi_P  = 1.0;   // Already included in ivphi_P
#ifdef NO_USECORE_DBG
              double Vbx    = ivphi_P_nuconly;
#else
              double Vbx    = ivphi_P;
#endif
#else
              // PRIOR form: vertex at rx
              double ivphi_T = InterpolateIVPHI(IVPHI_T, h_common, NSteps_common,
                                                 h_common * NSteps_common, rx,
                                                 IVPHI_T_max, r_T_vert_peak);
              double phi_T  = 1.0;   // Already included in ivphi_T
              double Vbx    = ivphi_T;  // = (V_nA * phi_T)'(rx)
              double phi_P  = InterpolateClipped(PrjBS_ch.WaveFunction, PrjBS_ch.RGrid,
                                                 PrjBS_ch.NSteps, PrjBS_ch.StepSize,
                                                 PrjBS_ch.MaxR, rp, phi_P_max, r_P_peak);
#endif

              // Angle of rx in the integration coordinate frame
              // Ptolemy line 16670: TEMP = (T1*RO + S1*RI*X)/RT
              double cos_phiT = (rx2 > 1e-30) ? (T1*rb + S1*ra*x) / rx : 1.0;
              cos_phiT = std::max(-1.0, std::min(1.0, cos_phiT));
              double phi_T_angle = std::acos(cos_phiT);

              // A12 angular coupling kernel — same formula as rectangular
              double A12_val = EvalA12(A12_terms, phi_T_angle, phi_ab);

              // phi_weight = PHI0 * phi_wts[k] * sin(phi) = DPHI in Ptolemy GRDSET
              AngKernel += phi_weight * phi_T * Vbx * phi_P * A12_val;
            }

            // Accumulate: RIROWTS = JACOB * ri * ro * WOW * DIFWT
            // Jacobian d(ri)d(ro)/d(U)d(V) = 1 (exact, verified analytically)
            // WOW = wi_s[is] (CUBMAP sum weight), DIFWT = wi_d[id] (CUBMAP dif weight)
            auto contrib = ca * cb_conj * AngKernel * JACOB_grdset * ra * rb * WOW * DIFWT;
            Integral += contrib;

            // DEBUG: print integrand at diagonal (V=0, ra=rb=U) for Li=0,Lo=2,Lx=2
#ifdef DEBUG_INELDC
            if (Li == 0 && Lo == 2 && Lx == 2 && JPO == 3) {
              double V = ra - rb;  // = U + V/2 - (U - V/2) = V
              if (std::abs(V) < 0.5 * ra * 0.1) {  // near-diagonal (V < 5% of 2U)
                fprintf(stderr, "DBG_HKernel Li=%d Lo=%d Lx=%d JPO=%d: "
                        "ra=%.4f rb=%.4f AngKernel=%.6e "
                        "ca=(%.4e,%.4e) cb=(%.4e,%.4e) "
                        "contrib=(%.4e,%.4e)\n",
                        Li, Lo, Lx, JPO, ra, rb, AngKernel,
                        ca.real(), ca.imag(), cb_conj.real(), cb_conj.imag(),
                        contrib.real(), contrib.imag());
              }
            }
#endif
          }
        }

#endif  // USE_ZR

        // Phase factor from Ptolemy SFROMI (source.mor ~line 29130):
        //   ITEST = LASI + LASO + 2*LXP + 1
        //   IF MOD(ITEST,4) >= 2: TEMP = -TEMP  (real sign flip)
        //   IF MOD(ITEST,2) == 1: multiply by i  (imaginary rotation)
        // Combined: phase = i^ITEST = i^(Li+Lo+2*Lx+1)
        int ITEST = ((Li + Lo + 2*Lx + 1) % 4 + 4) % 4;
        std::complex<double> phase_factor;
        switch (ITEST) {
          case 0: phase_factor = { 1.0,  0.0}; break;
          case 1: phase_factor = { 0.0,  1.0}; break;
          case 2: phase_factor = {-1.0,  0.0}; break;
          case 3: phase_factor = { 0.0, -1.0}; break;
        }

        Integral *= phase_factor;

        // Ptolemy SFROMI normalization (source.mor ~29130):
        //   TEMP = FACTOR * ATERM(LXP+1) / DSQRT(2*LASI+1)
        //   S_sfromi = TEMP * I * i^ITEST
        // where:
        //   FACTOR = 2*sqrt(AKI*AKO/(ES1*ES2))  [kinematic factor]
        //   ATERM(Lx) = sqrt((2jB+1)/(2jA+1)) * sqrt(2Lx+1) * SPAMP*SPAMT * RACAH(...)
        //     RACAH(2*LBT, JBT, 2*LBP, JBP, JX, 2*LX) with all args doubled
        //   LASI = Li (incident L), LASO = Lo
        //
        // Compute FACTOR_sfromi from kinematics
        double ka_sf = Incoming.k;
        double kb_sf = Outgoing.k;
        double Ea_sf = Incoming.Ecm;
        double Eb_sf = Outgoing.Ecm;
        double FACTOR_sfromi = 2.0 * std::sqrt(ka_sf * kb_sf / (Ea_sf * Eb_sf));
        // Factor of 2: Ptolemy SFROMI uses FACTOR = 2*sqrt(ki*kf/(Ei*Ef))
        // (see source.mor line 29025 comment: FACTOR = KINEMATIC FACTOR = 2*SQRT(KI*KF/EI*EF))
        
        // Compute ATERM generically (same formula as ineldc_zr.cpp).
        // Ptolemy source.mor ~25631:
        //   TEMP = sqrt((JBIGB+1)/(JBIGA+1))
        //   ATERM = TEMP * sqrt(2*Lx+1) * SPAMP * SPAMT * RACAH(2*LBT, JBT, 2*LBP, JBP, JX, 2*LX)
        //   where RACAH(a,b,c,d,e,f) = (-1)^((a+b+c+d)/2) * {a/2,b/2,e/2; d/2,c/2,f/2}  (6-j)
        //   Sign: ITEST = (JX - JBP + 2*(LBP+LBT))/2 + 1; if odd, flip ATERM
        double ATERM_val = 0.0;
        if (Lx >= std::abs(TargetBS.l - ProjectileBS.l) &&
            Lx <= TargetBS.l + ProjectileBS.l) {
          double jT_bs = TargetBS.j;          // j of neutron in target bound state
          double jP_bs = ProjectileBS.j;       // j of neutron in projectile bound state
          double jx = 0.5;                     // neutron spin

          // SixJ{lBT, jBT, jX; jBP, lBP, Lx}
          double sj = SixJ((double)TargetBS.l, jT_bs, jx,
                           jP_bs, (double)ProjectileBS.l, (double)Lx);

          // Phase: (-1)^((2*lBT + 2*jBT + 2*lBP + 2*jBP)/2) = RACAH phase
          int twoj_sum = 2*TargetBS.l + (int)(2*jT_bs + 0.5)
                       + 2*ProjectileBS.l + (int)(2*jP_bs + 0.5);
          double sign_val = ((twoj_sum / 2) % 2 == 0) ? 1.0 : -1.0;
          double RACAH_val = sign_val * sj;

          // Nuclear spins from DWBA object (generic)
          int JBIGA = (int)std::round(2.0 * SpinTarget);    // 2*J(target nucleus)
          int JBIGB = (int)std::round(2.0 * SpinResidual);  // 2*J(residual nucleus)
          double TEMP_aterm = std::sqrt((JBIGB + 1.0) / (JBIGA + 1.0));

          // Spectroscopic amplitudes
          double SPAMP = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
          double SPAMT = 1.0;

          ATERM_val = TEMP_aterm * std::sqrt(2.0*Lx + 1.0) * SPAMP * SPAMT * RACAH_val;

          // Sign convention (Ptolemy source.mor line 25636-25639):
          // ITEST = JX - JBP + 2*(LBP+LBT)   [all doubled integers]
          int JX_doubled  = 1;  // neutron spin × 2
          int JBP_doubled = (int)(2 * ProjectileBS.j);
          int ITEST_aterm = JX_doubled - JBP_doubled + 2*(ProjectileBS.l + TargetBS.l);
          if ((ITEST_aterm / 2 + 1) % 2 != 0) ATERM_val = -ATERM_val;
        }
        
        // Total SFROMI factor: FACTOR_sfromi * ATERM * (1/sqrt(2*Li+1))
        double sfromi_norm = FACTOR_sfromi * std::abs(ATERM_val) / std::sqrt(2.0 * Li + 1.0);
        // Apply sign of ATERM to the integral (ATERM sign affects overall sign)
        if (ATERM_val < 0) {
          Integral = -Integral;
        }

        // Store: (Lx, Li, Lo, JPI, JPO) — SFROMI will apply 9-J coupling
        auto S_pre9j = Integral * sfromi_norm;
        TransferSMatrix.push_back({Lx, Li, Lo, JPI, JPO, S_pre9j});

        // Print pre-9J S-matrix (compare vs Ptolemy PRINT=2 "S-MATRIX BEFORE 9-J" table)
        fprintf(stderr, "PRE9J Li=%2d JPI=%2d/2  Lo=%2d JPO=%2d/2  Lx=%d  "
                "S=(%9.4e, %9.4e)  |S|=%9.4e\n",
                Li, JPI, Lo, JPO, Lx,
                S_pre9j.real(), S_pre9j.imag(), std::abs(S_pre9j));

        } // end JPO loop
        } // end JPI loop
      } // end Lx loop
    } // end Lo loop
  } // end Li loop


}
