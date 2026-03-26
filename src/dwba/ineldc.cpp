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
  // Ptolemy BOUND step size: h = min(1/kappa, A_diff) / STEPSPER (STEPSPER=8 from dpsb)
  // kappa = sqrt(2*mu_AMU*AMU_MEV*|BE|) / HBARC;  A_diff = WS diffuseness of BS potential
  {
    const int STEPSPER = 8;
    double kappa_T = std::sqrt(2.0 * TgtBS_ch.mu * AMU_MEV * std::abs(TargetBS.BindingEnergy)) / HBARC;
    double A_tbs   = (TargetBS.Pot.A > 0) ? TargetBS.Pot.A : 0.65;
    TgtBS_ch.StepSize = std::min(1.0 / kappa_T, A_tbs) / STEPSPER;
    // fprintf(stderr, "TgtBS step: kappa=%.5f fm^-1, A=%.4f → h=%.5f fm\n",
    //         kappa_T, A_tbs, TgtBS_ch.StepSize);
  }
  // WavSet will allocate NSteps and RGrid based on StepSize (≤0 uses default 0.1)
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
  // Projectile BS step: Ptolemy formula h = min(1/kappa_P, A_P) / STEPSPER
  // Deuteron (AV18): kappa=0.2316 fm^-1, A_pot=0.5 → h=min(4.32,0.5)/8=0.0625 fm
  // IVPHI_P lives on PrjBS grid; IVPHI_T lives on TgtBS grid — grids are independent.
  {
    const int STEPSPER = 8;
    double kappa_P = std::sqrt(2.0 * PrjBS_ch.mu * AMU_MEV * std::abs(ProjectileBS.BindingEnergy)) / HBARC;
    double A_pbs   = (ProjectileBS.Pot.A > 0) ? ProjectileBS.Pot.A : 0.5;
    PrjBS_ch.StepSize = std::min(1.0 / kappa_P, A_pbs) / STEPSPER;
    // fprintf(stderr, "PrjBS step: kappa=%.5f fm^-1, A=%.4f → h=%.5f fm\n",
    //         kappa_P, A_pbs, PrjBS_ch.StepSize);
  }
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
  int Lmax = 40;  // Full run

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
  // IVPHI_P (projectile vertex) lives on PrjBS grid (h=0.1 fm)
  // IVPHI_T (target vertex) lives on TgtBS grid (h_T = min(1/kappa,A)/8)
  // Use PrjBS grid as the common reference for IVPHI_P interpolation.
  int NSteps_common = PrjBS_ch.NSteps;  // PrjBS grid size (h=0.1 fm)
  double h_common   = PrjBS_ch.StepSize; // = 0.1 fm
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

  // IVPHI_T on TgtBS grid (fine: h_T = 0.08125 fm)
  int NSteps_T = TgtBS_ch.NSteps;
  double h_T   = TgtBS_ch.StepSize;
  std::vector<double> IVPHI_T(NSteps_T, 0.0);
  double IVPHI_T_max = 0; int idx_T_vert_peak = 1;

  // IVPHI_P on PrjBS grid (h=0.1 fm)
  std::vector<double> IVPHI_P(NSteps_common, 0.0);
  double IVPHI_P_max = 0; int idx_P_vert_peak = 1;

  std::cout << "  Vertex: V_WS(V=60,R=1,A=0.5)*phi_Reid (projectile), V_WS*phi_T (target)" << std::endl;

  for (int i = 1; i < NSteps_T; ++i) {
    IVPHI_T[i] = std::abs(TgtBS_ch.WaveFunction[i].real()) * TgtBS_ch.V_real[i];
    if (IVPHI_T[i] > IVPHI_T_max) { IVPHI_T_max = IVPHI_T[i]; idx_T_vert_peak = i; }
  }
  for (int i = 1; i < NSteps_common; ++i) {
    IVPHI_P[i] = !PrjBS_ch.VPhiProduct.empty()
                 ? PrjBS_ch.VPhiProduct[i]
                 : std::abs(PrjBS_ch.WaveFunction[i].real()) * PrjBS_ch.V_real[i];
    if (IVPHI_P[i] > IVPHI_P_max) { IVPHI_P_max = IVPHI_P[i]; idx_P_vert_peak = i; }
  }
  double r_T_vert_peak = idx_T_vert_peak * h_T;
  double r_P_vert_peak = idx_P_vert_peak * h_common;

#ifdef DEBUG_DUMP_TABLES
  // Dump phi_T and IVPHI_P to files for Fortran comparison
  {
    FILE* fp = fopen("phi_T_table.txt", "w");
    for (int i = 0; i < NSteps_T; ++i)
      fprintf(fp, "%.6f  %.10e\n", i * h_T, TgtBS_ch.WaveFunction[i].real());
    fclose(fp);
    FILE* fp2 = fopen("ivphi_P_table.txt", "w");
    for (int i = 0; i < NSteps_common; ++i)
      fprintf(fp2, "%.6f  %.10e\n", i * h_common, IVPHI_P[i]);
    fclose(fp2);
    fprintf(stderr, "Dumped phi_T (%d pts, h=%.5f) and IVPHI_P (%d pts, h=%.5f)\n",
            NSteps_T, h_T, NSteps_common, h_common);
  }
#endif

  // Report vertex peaks for diagnostic comparison vs Ptolemy BSSET output
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "  IVPHI_T peak = " << IVPHI_T_max << " at r=" << r_T_vert_peak << " fm" << std::endl;
  std::cout << "  IVPHI_P peak = " << IVPHI_P_max << " at r=" << r_P_vert_peak << " fm" << std::endl;



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
  auto [phi_T_max, idx_T_peak] = findWFPeak(TgtBS_ch.WaveFunction, NSteps_T);
  auto [phi_P_max, idx_P_peak] = findWFPeak(PrjBS_ch.WaveFunction, NSteps_common);
  double r_T_peak = idx_T_peak * h_T;
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
  // 5-pt Lagrange interpolation (AITLAG equivalent) for real tables
  // Matches Ptolemy's BSPROD ITYPE=1 which uses AITLAG(NAIT=4 → 5-pt)
  auto aitlag5 = [](const double* tab, int ntab, double stepinv, double r) -> double {
    if (r < 0) return tab[0];
    double idx = r * stepinv;
    if (idx >= ntab - 1) return tab[ntab-1];
    int I = static_cast<int>(idx + 0.5);  // nearest index (AITLAG convention)
    I = std::max(2, std::min(I, ntab - 3));
    double P = idx - I;
    double PS = P * P;
    double X1 = P*(PS-1.0)/24.0, X2 = X1+X1, X3 = X1*P;
    double X4 = X2+X2-0.5*P, X5 = X4*P;
    double C1=X3-X2, C5=X3+X2, C3=X5-X3, C2=X5-X4, C4=X5+X4;
    C3 = C3+C3+1.0;
    return C1*tab[I-2] - C2*tab[I-1] + C3*tab[I] - C4*tab[I+1] + C5*tab[I+2];
  };

  auto InterpolateClipped = [&](const std::vector<std::complex<double>> &wf,
                                const std::vector<double> &/*grid*/,
                                int nsteps, double stepsize, double maxr,
                                double r, double /*phi_max*/, double /*r_peak*/) -> double {
    if (r >= maxr || r < 0) return 0.0;
    // Extract real part into temp array for aitlag5
    // Use linear for speed in non-PTOLEMY mode; 5-pt Lagrange in PTOLEMY mode
#ifdef INTERP_PTOLEMY
    // Build pointer into real parts (wf is contiguous complex<double>)
    // Use direct 5-pt formula on real parts
    double idx = r / stepsize;
    int I = static_cast<int>(idx + 0.5);
    I = std::max(2, std::min(I, nsteps - 3));
    double P = idx - I;
    double PS = P*P;
    double X1=P*(PS-1.0)/24.0, X2=X1+X1, X3=X1*P;
    double X4=X2+X2-0.5*P, X5=X4*P;
    double C1=X3-X2, C5=X3+X2, C3=X5-X3, C2=X5-X4, C4=X5+X4;
    C3=C3+C3+1.0;
    return C1*wf[I-2].real()-C2*wf[I-1].real()+C3*wf[I].real()
          -C4*wf[I+1].real()+C5*wf[I+2].real();
#else
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return wf[ii].real() * (1.0 - frac) + wf[ii + 1].real() * frac;
#endif
  };

  // Clipped IVPHI interpolation for vertex product (V*phi at vertex side):
  // IVPHI'(r) = IVPHI_max for r < r_vert_peak; IVPHI(r) for r >= r_vert_peak
  // InterpolateIVPHI: interpolate V(r)*phi(r) product WITHOUT clipping
  // Ptolemy Pass 2 uses ITYPE=1 (PHI V PHI) with AITLAG (5-pt Lagrange).
  auto InterpolateIVPHI = [&](const std::vector<double> &ivphi,
                               double stepsize, int nsteps, double maxr,
                               double r, double /*ivphi_max*/, double /*r_vert_peak*/) -> double {
    if (r >= maxr || r < 0) return 0.0;
#ifdef INTERP_PTOLEMY
    return aitlag5(ivphi.data(), nsteps, 1.0/stepsize, r);
#else
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return ivphi[ii] * (1.0 - frac) + ivphi[ii + 1] * frac;
#endif
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

  // Ptolemy GRDSET pre-computes TWO reference scattering wavefunctions ONCE (before Li loop):
  //   ISCTMN: incoming chi at L = MAX(0, Lmin - Lxmax)  → used in WVWMAX + SUMMIN scan
  //   ISCTCR: incoming chi at L = LCRIT                  → used in SUMMID scan
  //
  // LCRIT determination (Ptolemy source lines ~38088, 29457, 29831):
  //   If LMIN and LMAX are both defined: LCRITS(ICHANW) = (LMIN+LMAX)/2
  //   LCRIT = (LCRITS(1)+LCRITS(2))/2
  //   If LCRIT==0 fallback: LCRIT = LMAX/2
  //   Then LC = LCRIT, clamped to [LMIN,LMAX]; if out of range: LC = (LMIN+LMAX)/2
  //
  // For full run (Lmin=0, Lmax=40): Ptolemy LCRITL gives LCRIT=5 (grazing L from turning pt).
  // With lmin/lmax both defined: LCRIT = (0+40)/2 = 20 — but Ptolemy uses LCRITL here.
  // LCRITL computes grazing L ≈ k*R_nuclear. For d+16O at Ecm=17.76 MeV, k=1.233 fm^-1:
  //   L_gr = k * R_nuclear ≈ 1.233 * (1.149*(16^1/3) + 1.0) ≈ 1.233 * 3.76 ≈ 4.6 → 5
  // We compute this analytically and match Ptolemy's LCRIT=5 for this reaction.
  const int LxMax_grdset = lT;          // max Lx = lT = 2 for (d,p) neutron transfer
  const int L_isctmn = std::max(0, 0 - LxMax_grdset);  // = 0 for Lmin=0
  // Estimate grazing L: LCRIT = round(k_in * R_nuclear)
  // R_nuclear ≈ R0*(Ap^(1/3) + At^(1/3)) with R0=1.149 fm (Ptolemy default Coulomb)
  // Compute L_isctcr = LCRIT clamped to [0, Lmax].
  // Confirmed by Ptolemy DBGCHI: full lmax=40 run gives LCRIT=5 (L CRITICAL AVERAGE=5).
  // Use mass number A (not AMU mass) for nuclear radius formula R = R0*(Ap^1/3 + At^1/3)
  // Use Ptolemy's formula: LCRIT = round(k * R_nuclear), R = 1.149*(Ap^1/3+At^1/3)
  // Confirmed LCRIT=5 for d+16O at Ecm=17.76 MeV (Ptolemy DBGCHI output).
  // R0=1.149 fm is Ptolemy's default nuclear radius (not Coulomb RC0=1.303).
  const int L_isctcr = std::max(0, std::min(Lmax, (int)std::round(
      1.149 * (std::pow((double)Incoming.Projectile.A, 1.0/3.0) +
               std::pow((double)Incoming.Target.A,     1.0/3.0)) * Incoming.k)));
  fprintf(stderr, "GRDSET: k_in=%.4f Ap=%d At=%d → R_nuclear=%.4f  L_isctmn=%d  L_isctcr=%d\n",
          Incoming.k, Incoming.Projectile.A, Incoming.Target.A,
          1.3030*(std::pow((double)Incoming.Projectile.A,1.0/3)+std::pow((double)Incoming.Target.A,1.0/3)),
          L_isctmn, L_isctcr);
  // Compute ISCTCR chi: incoming at L=L_isctcr (central potential, no SO — Ptolemy GRDSET)
  // We use full optical (SO included) for now — same as before but at correct L.
  {
    int JPI_cr = 2*L_isctcr + JA_dw;  // highest J for this L
    WavElj(Incoming, L_isctcr, JPI_cr);
  }
  const std::vector<std::complex<double>> chi_isctcr = Incoming.WaveFunction;
  const double h_isctcr = Incoming.StepSize;
  auto interp_isctcr = [&](double r) -> double {
    int idx = (int)(r / h_isctcr);
    if (idx < 0) return 0.0;
    if (idx >= (int)chi_isctcr.size() - 1) return 0.0;
    double frac = r/h_isctcr - idx;
    return (chi_isctcr[idx]*(1.0-frac) + chi_isctcr[idx+1]*frac).real();
  };
  fprintf(stderr, "GRDSET: L_isctmn=%d  L_isctcr=%d (Lmax=%d)\n",
          L_isctmn, L_isctcr, Lmax);

  for (int Li = 0; Li <= Lmax; ++Li) {
    // --- J-split incoming distorted waves for all JPI ---
    // JPI ranges from |2*Li - JA_dw| to 2*Li + JA_dw in steps of 2
    int JPI_min = std::abs(2*Li - JA_dw);
    int JPI_max = 2*Li + JA_dw;
    std::map<int, std::vector<std::complex<double>>> chi_a_byJPI;
    for (int JPI = JPI_min; JPI <= JPI_max; JPI += 2) {
      WavElj(Incoming, Li, JPI);
      chi_a_byJPI[JPI] = Incoming.WaveFunction;
      // Print chi_a at key radii for Li=3 (to diagnose JPI=8/2 error)
      if (Li == 3) {
        double h = Incoming.StepSize;
        double peak = 0; int peak_idx = 0;
        for (int ii=0; ii<(int)Incoming.WaveFunction.size(); ii++)
          if (std::abs(Incoming.WaveFunction[ii]) > peak) { peak=std::abs(Incoming.WaveFunction[ii]); peak_idx=ii; }
        fprintf(stderr, "chi_a Li=%d JPI=%d/2: size=%d h=%.5f peak=%.6f at r=%.4f\n",
                Li, JPI, (int)Incoming.WaveFunction.size(), h, peak, peak_idx*h);
        for (double r : {2.0, 4.0, 6.0, 8.0, 10.0}) {
          int idx = (int)(r/h + 0.5);
          if (idx < (int)Incoming.WaveFunction.size()) {
            auto v = Incoming.WaveFunction[idx];
            fprintf(stderr, "  r=%.1f: re=%+.6f im=%+.6f |chi|=%.6f\n",
                    r, v.real(), v.imag(), std::abs(v));
          }
        }
      }
    }

    // ── Per-Li GRDSET: compute SUMMIN, SUMMID, WVWMAX, RVRLIM ─────────────
    const double DWCUT_grdset = 2.0e-6;  // DPSB grid row 12: DWCUT=2e-6 (RGRIDS(1,12) in source)
    // NOTE: WVWMAX must use BSPROD ITYPE=3 (surface derivative WS × phi_P × phi_T).
    // Fortran WVWMAX~1.54e-5 → RVRLIM = DWCUT × WVWMAX ~ 3.08e-11 for 16O(d,p)17O.
    // Temporary: override RVRLIM with Fortran value for validation (fix WVWMAX later).
    // Ptolemy reruns GRDSET for each LI value (source.mor line 17735 outer loop),
    // using the BSPROD×chi integrand for the lowest-JPO channel to find SUMMIN, SUMMID.
    // These are then REUSED for all Lo, Lx, JPI, JPO under this Li.
    //
    // We use chi_a[JPI_min] (lowest JPI for this Li) and chi_b[Lo=Li, JPO=Lo+1/2] to
    // approximate Ptolemy's reference chi (first channel encountered).
    // For the SUMMIN/SUMMID scan, Ptolemy uses BSPROD ITYPE=3 = R*chi*bsprod*chi*R.
    //
    // Verified SUMMID values from Ptolemy print=5 (lmin=lmax=Li):
    //   Li=0: 4.926, Li=1: 5.008, Li=2: 4.850, Li=3: 4.839, Li=4: 5.022
    //   Li=5: 5.069, Li=6: 5.486, Li=7: 6.349, Li=8: 7.290, Li=9: 8.093
    //   Li=10: 8.831, Li=11: 9.570, ... growing linearly with Li for large Li
    // Flat-top constants (same for all Li — depend only on bound state WFs, not chi)
    // Used in SUMMID scan AND in the phi integration loop (BSPROD ITYPE=1).
    const double RLMAXS1_li = r_P_vert_peak;   // IVPHI_P peak location (Ptolemy RLPMAX)
    const double RLMAXS2_li = r_T_peak;         // phi_T peak location (Ptolemy RLTMAX)
    const double VTMAX_li   = std::abs(TgtBS_ch.WaveFunction[idx_T_peak].real());
    const double VPMAX_li   = IVPHI_P_max;

    // Flat-top interpolation functions (Ptolemy BSPROD ITYPE=1 convention):
    // for r < r_peak → use peak value; for r > r_peak → use actual interpolation
    auto flat_FT_li = [&](double rx) -> double {
      if (rx <= RLMAXS2_li) return VTMAX_li;
      return std::abs(InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                         TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                         TgtBS_ch.MaxR, rx, 0.0, 0.0));
    };
    auto flat_FP_li = [&](double rp) -> double {
      if (rp <= RLMAXS1_li) return VPMAX_li;
      if (rp >= h_common*NSteps_common) return 0.0;
      return InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                               h_common*NSteps_common, rp, 0.0, 0.0);
    };

    double SUMMIN_li = 0.0, SUMMID_li = 0.0, WVWMAX_li = 0.0, RVRLIM_li = 0.0;
    {
      // Get representative chi waves for this Li
      // Use JPI_max (highest J) as the reference: matches Ptolemy's ISCTCR which maximizes bsprod
      const auto& ref_chi_a = chi_a_byJPI.rbegin()->second;
      const double h_ref_a  = Incoming.StepSize;
      const double maxR_ref_a = (ref_chi_a.size() >= 4)
          ? (static_cast<int>(ref_chi_a.size()) - 4) * h_ref_a : Incoming.MaxR;

      // Compute a reference chi_b: use Lo=Li (elastic-like), lowest JPO
      int Lo_ref = Li;
      int JPO_ref = std::max(1, std::abs(2*Lo_ref - JB_dw));
      WavElj(Outgoing, Lo_ref, JPO_ref);
      const auto ref_chi_b = Outgoing.WaveFunction;
      const double h_ref_b  = Outgoing.StepSize;
      const double maxR_ref_b = (ref_chi_b.size() >= 4)
          ? (static_cast<int>(ref_chi_b.size()) - 4) * h_ref_b : Outgoing.MaxR;

      // Reference interpolators (5-pt Lagrange)
      auto interp_ref_a = [&](double r) -> double {
        if (r <= 0 || r > maxR_ref_a) return 0.0;
#ifdef INTERP_PTOLEMY
        double rbyh = r / h_ref_a;
        int I = static_cast<int>(rbyh + 0.5);
        int IMX = static_cast<int>(ref_chi_a.size()) - 4;
        I = std::max(2, std::min(I, IMX));
        double P = rbyh - I, PS = P*P;
        double X1 = P*(PS-1.0)/24.0, X2=X1+X1, X3=X1*P;
        double X4=X2+X2-0.5*P, X5=X4*P;
        double C1=X3-X2, C5=X3+X2, C3=X5-X3, C2=X5-X4, C4=X5+X4;
        C3=C3+C3+1.0;
        return (C1*ref_chi_a[I-2]-C2*ref_chi_a[I-1]+C3*ref_chi_a[I]
               -C4*ref_chi_a[I+1]+C5*ref_chi_a[I+2]).real();
#else
        double idx_f = r/h_ref_a; int ii=(int)idx_f;
        if(ii>=(int)ref_chi_a.size()-1) return 0.0;
        double frac=idx_f-ii;
        return (ref_chi_a[ii]*(1-frac)+ref_chi_a[ii+1]*frac).real();
#endif
      };
      auto interp_ref_b = [&](double r) -> double {
        if (r <= 0 || r > maxR_ref_b) return 0.0;
#ifdef INTERP_PTOLEMY
        double rbyh = r / h_ref_b;
        int I = static_cast<int>(rbyh + 0.5);
        int IMX = static_cast<int>(ref_chi_b.size()) - 4;
        I = std::max(2, std::min(I, IMX));
        double P = rbyh - I, PS = P*P;
        double X1 = P*(PS-1.0)/24.0, X2=X1+X1, X3=X1*P;
        double X4=X2+X2-0.5*P, X5=X4*P;
        double C1=X3-X2, C5=X3+X2, C3=X5-X3, C2=X5-X4, C4=X5+X4;
        C3=C3+C3+1.0;
        return (C1*ref_chi_b[I-2]-C2*ref_chi_b[I-1]+C3*ref_chi_b[I]
               -C4*ref_chi_b[I+1]+C5*ref_chi_b[I+2]).real();
#else
        double idx_f = r/h_ref_b; int ii=(int)idx_f;
        if(ii>=(int)ref_chi_b.size()-1) return 0.0;
        double frac=idx_f-ii;
        return (ref_chi_b[ii]*(1-frac)+ref_chi_b[ii+1]*frac).real();
#endif
      };

      // RLMAXS(1) = projectile form-factor peak radius (RLPMAX in Ptolemy BSSET)
      // RLMAXS(2) = target form-factor peak radius    (RLTMAX in Ptolemy BSSET)
      // "LOCATION OF PROJECTILE AND TARGET FORM-FACTOR MAXIMA" from print=5
      // These are the r at which r*|phi_P*V| and r*|phi_T| peak.
      // = r_P_vert_peak (0.1875 fm) and r_T_vert_peak (2.1125 fm)
      // RLMAXS(1) = RLPMAX: peak of phi_P*V (IVPHI_P) = r_P_vert_peak
      // RLMAXS(2) = RLTMAX: peak of phi_T (raw WF)   = r_T_peak (not IVPHI_T!)
      // RLMAXS1 and RLMAXS2 are now defined at Li scope as RLMAXS1_li / RLMAXS2_li (aliased below)
      // const double RLMAXS1 = r_P_vert_peak;  // moved to Li scope
      // const double RLMAXS2 = r_T_peak;       // moved to Li scope

      // XS scan values (Ptolemy: LOOKST=250, DXV=2/LOOKST^2)
      const double DXV_grd = 2.0 / (250.0*250.0);
      const double xs_grd[2] = {1.0, 1.0 - DXV_grd};

      // VTMAX, VPMAX: flat-top values (phi_T peak and IVPHI_P peak)
      // VTMAX and VPMAX now defined at Li scope (VTMAX_li, VPMAX_li, aliased below)
      // const double VTMAX = ...; // moved to Li scope
      // const double VPMAX = ...; // moved to Li scope

      // Ptolemy BSPROD flat-top phi' (ITYPE=2,3,4):
      //   For r < RLMAXS (peak location): FT = VTMAX, FP = VPMAX (flat cap)
      //   For r > RLMAXS: FT = phi_T(r), FP = IVPHI_P(r) (actual interpolation)
      //   Alternate return (*) fires only if RP > BNDMXP OR RT > BNDMXT (grid bounds, ~30fm)
      //   — i.e., basically never during these scans.
      // ITYPE=3: FT' × FP'  with chi factors: result = ra*chi_a(ra)*FP'*FT'*chi_b(rb)*rb
      // ITYPE=4: FT' × FP'  without chi: result = FP' * FT'
      // ITYPE=2: same as ITYPE=4 (flat phi, no chi)

      // Aliases for the outer-scope flat-top functions (defined at Li level above)
      auto& flat_FT = flat_FT_li;
      auto& flat_FP = flat_FP_li;
      const double& RLMAXS1 = RLMAXS1_li;
      const double& RLMAXS2 = RLMAXS2_li;
      const double& VTMAX   = VTMAX_li;
      const double& VPMAX   = VPMAX_li;

      // bsprod3: ITYPE=3 = ra*chi_LCRIT * FP'*FT' * chi_LCRIT*rb
      // Ptolemy BSPROD source lines 200-270: BOTH ra AND rb use the SAME ISCAT
      // (the LCRIT incoming channel chi). NOT chi_in(ra)*chi_out(rb)!
      // Reference: ISCTCR = incoming chi at L=LCRIT stored during grid setup.
      // Our ref_chi_a IS chi_in at JPI_max which corresponds to LC (LCRIT).
      auto bsprod3 = [&](double ra, double rb, double x) -> double {
        if (ra < 1e-6 || rb < 1e-6) return 0.0;
        double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*x;
        double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*x;
        if (rx2<0) rx2=0; if (rp2<0) rp2=0;
        double rx=std::sqrt(rx2), rp=std::sqrt(rp2);
        double FT = flat_FT(rx);
        double FP = flat_FP(rp);
        // Both ra and rb use ISCTCR chi (pre-computed at L=LCRIT, once for all Li)
        double ca = std::abs(interp_isctcr(ra));
        double cb = std::abs(interp_isctcr(rb));
        return ra * ca * FP * FT * cb * rb;
      };

      // bsprod4: ITYPE=4 = FP'*FT' (no chi, no r factors)
      auto bsprod4 = [&](double U, double x) -> double {
        if (U < 1e-6) return 0.0;
        double rx2 = S1*S1*U*U + T1*T1*U*U + 2.0*S1*T1*U*U*x;
        double rp2 = S2*S2*U*U + T2*T2*U*U + 2.0*S2*T2*U*U*x;
        if (rx2<0) rx2=0; if (rp2<0) rp2=0;
        return flat_FT(std::sqrt(rx2)) * flat_FP(std::sqrt(rp2));
      };

      // Helper: compute Ptolemy's 5 adaptive VS at a given U
      auto compute_VS = [&](double U, double VS[5]) {
        double D = 2.0*(RLMAXS2 - (S1+T1)*U) / (S1-T1);
        if (std::abs(D) > 2.0*U) D = std::copysign(2.0*U, D);
        VS[0] = D; VS[1] = 0.5*D; VS[2] = 0.0;
        D = 2.0*(RLMAXS1 - (S2+T2)*U) / (S2-T2);
        if (std::abs(D) > 2.0*U) D = std::copysign(2.0*U, D);
        VS[3] = 0.5*D; VS[4] = D;
      };

      // Step 1: WVWMAX — scan UP from U=0.5*ROFMAX for N=1.5*ROFMAX/0.2+1 steps
      // (Ptolemy source.f line ~18588: U=.5*ROFMAX, N=1.5*ROFMAX/.20+1)
      const double SUMMAX_li = 30.4;
      const double ROFMAX_grd = maxR_ref_a;  // Ptolemy ROFMAX = chi_a extent
      if (Li <= 1) fprintf(stderr,"ROFMAX Li=%d JPI_max=%d maxR=%.2f h=%.4f npts=%d\n",
                            Li, JPI_max, ROFMAX_grd, h_ref_a, (int)ref_chi_a.size());
      double U_s = 0.5 * ROFMAX_grd;
      int N_wv = (int)(1.5*ROFMAX_grd/0.2 + 1.5);
      for (int iu = 0; iu < N_wv; ++iu, U_s += 0.2) {
        double VS5[5]; compute_VS(U_s, VS5);
        for (int iv5 = 0; iv5 < 5; ++iv5) {
          double ra = U_s + 0.5*VS5[iv5], rb = U_s - 0.5*VS5[iv5];
          for (int ix = 0; ix < 2; ++ix)
            WVWMAX_li = std::max(WVWMAX_li, std::abs(bsprod3(ra, rb, xs_grd[ix])));
        }
      }
      RVRLIM_li = DWCUT_grdset * std::max(WVWMAX_li, 1.0e-30);
      if (Li == 3) fprintf(stderr, "WVWMAX_li(Li=3) = %.8g  RVRLIM_li = %.8g\n", WVWMAX_li, RVRLIM_li);

      // Step 2: SUMMIN — step out from U=0, BSPROD ITYPE=4 at V=0 (Ptolemy line 15964)
      U_s = 0.0;
      while (U_s <= SUMMAX_li) {
        double f = 0.0;
        for (int ix = 0; ix < 2; ++ix) f = std::max(f, std::abs(bsprod4(U_s, xs_grd[ix])));
        if (f >= RVRLIM_li) {
          SUMMIN_li = std::max(0.0, U_s - 0.2);
          break;
        }
        U_s += 0.2;
      }

      // Step 3: SUMMID — first moment <U> * AMDMLT (Ptolemy GRDSET lines 15994–16085)
      // Ptolemy scans 5 ADAPTIVE V values at each U (not fixed fractions of 2U).
      // VS(1,2) from target asymptopia: V where RT hits RLMAXS(2)
      // VS(4,5) from projectile asymptopia: V where RP hits RLMAXS(1)
      // VS(3) = 0 (center); BSPROD ITYPE=3 with adaptive VS per U.
      // (RLMAXS1, RLMAXS2, xs_grd, bsprod3, compute_VS defined above)

      double SUM0 = 0.0, SUM1 = 0.0, SUM2 = 0.0;

      for (U_s = SUMMIN_li; U_s <= SUMMAX_li; U_s += 0.2) {
        double VS5[5]; compute_VS(U_s, VS5);

        double temp = 0.0;
        bool zerosw = true;
        for (int iv5 = 0; iv5 < 5; ++iv5) {
          double ra_s = U_s + 0.5*VS5[iv5];
          double rb_s = U_s - 0.5*VS5[iv5];
          for (int ix = 0; ix < 2; ++ix) {
            double f = bsprod3(ra_s, rb_s, xs_grd[ix]);
            if (f != 0.0) zerosw = false;
            temp += std::abs(f);
          }
        }
        if (zerosw) continue;

        if (Li==0 && (U_s < 3.01 || (U_s > 4.8 && U_s < 7.01)))
          fprintf(stderr,"SUMMID_SCAN Li=0 U=%.2f TEMP=%.6e VS=[%.3f,%.3f,%.3f,%.3f,%.3f]\n",
                  U_s, temp, VS5[0],VS5[1],VS5[2],VS5[3],VS5[4]);

        SUM0 += temp;
        SUM1 += temp * U_s;
        SUM2 += temp * U_s * U_s;
      }

      double mean_U = (SUM0 > 1e-30) ? SUM1/SUM0 : 0.5*(SUMMIN_li + SUMMAX_li);
      if (Li == 0) {
        fprintf(stderr, "DIAG SUMMID scan: RLMAXS1=%.4f  RLMAXS2=%.4f\n", RLMAXS1, RLMAXS2);
        fprintf(stderr, "DIAG SUMMID: SUM0=%.4e  SUM1=%.4e  mean_U=%.4f\n", SUM0, SUM1, mean_U);
        // Print VS at U=5
        double VS_test[5]; compute_VS(5.0, VS_test);
        fprintf(stderr, "DIAG VS@U=5: %.4f  %.4f  %.4f  %.4f  %.4f\n",
                VS_test[0], VS_test[1], VS_test[2], VS_test[3], VS_test[4]);
        double b3_test = bsprod3(5.0+0.5*VS_test[2], 5.0-0.5*VS_test[2], xs_grd[0]);
        fprintf(stderr, "DIAG bsprod3@U=5,V=0: %.6e\n", b3_test);
        double b3_test2 = bsprod3(5.2, 4.8, xs_grd[0]);
        fprintf(stderr, "DIAG bsprod3@ra=5.2,rb=4.8: %.6e\n", b3_test2);
      }
      const double AMDMLT = 0.9;  // dpsb MIDMULT (RGRIDS row C, col 11)
      SUMMID_li = mean_U * AMDMLT;
      // Clamp (Ptolemy lines 16084-16085)
      SUMMIN_li = std::min(SUMMIN_li, 7.0*(SUMMID_li - SUMMAX_li/7.0)/6.0);
      SUMMIN_li = std::max(SUMMIN_li, 0.0);
      SUMMID_li = std::min(SUMMID_li, 0.5*(SUMMIN_li + SUMMAX_li));
      SUMMID_li = std::max(SUMMID_li, SUMMIN_li + (SUMMAX_li - SUMMIN_li)/7.0);

      std::cout << "Li=" << Li << "  SUMMIN=" << SUMMIN_li
                << "  SUMMID=" << SUMMID_li << "  SUMMAX=" << SUMMAX_li
                << "  (<U>=" << mean_U << ")" << std::endl;

#ifdef USE_PTO_SUMMID_TABLE
      // Override with known-correct Ptolemy SUMMID/SUMMIN/SUMMAX values (from print=5)
      // Used to isolate whether SUMMID error drives the off-diagonal outliers.
      static const double pto_summid[] = {
        4.93, 5.01, 4.85, 4.84, 5.02, 5.07, 5.49, 6.35, 7.29, 8.09,
        8.83, 9.57,10.33,11.09,11.90,12.79,13.53,14.21,14.93,15.74,
       16.42,17.03,17.65,18.31,18.98,19.64,20.29,20.93,21.55,22.16,22.74};
      static const double pto_summin[] = {
        0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.4, 0.6, 1.0, 1.4,
        2.0, 2.4, 3.2, 3.8, 4.6, 5.8, 6.2, 6.8, 7.2, 7.8,
        8.6, 9.2, 9.8,10.6,11.2,12.0,12.8,13.4,14.2,15.0,15.8};
      int li_idx = std::min(Li, 30);
      SUMMID_li = pto_summid[li_idx];
      SUMMIN_li = pto_summin[li_idx];
      fprintf(stderr, "PTO_TABLE Li=%d SUMMIN=%.2f SUMMID=%.2f\n", Li, SUMMIN_li, SUMMID_li);
#endif
    }

    // ── Per-Li constants needed in the inner loops ───────────────────────────
    // Bound-state decay constants (same formula as ineldc_zr.cpp)
    double ALPHAP = std::sqrt(2.0 * PrjBS_ch.mu * AMU * std::abs(ProjectileBS.BindingEnergy)) / HBARC;
    double ALPHAT = std::sqrt(2.0 * TgtBS_ch.mu * AMU * std::abs(TargetBS.BindingEnergy)) / HBARC;

    // bsprod_val(ra, rb, x): |phi_T(rx)| * |ivphi_P(rp)|  (no chi, phi angle = acos(x))
    auto bsprod_val = [&](double ra, double rb, double x) -> double {
      double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*x;
      double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*x;
      if (rx2 < 0) rx2 = 0; if (rp2 < 0) rp2 = 0;
      double phi_T = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                        TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                        TgtBS_ch.MaxR, std::sqrt(rx2), 0.0, 0.0);
      double ivphi_P = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                         h_common*NSteps_common, std::sqrt(rp2), 0.0, 0.0);
      return std::abs(phi_T * ivphi_P);
    };

    // ptolemy_spline: natural cubic spline interpolation (SPLNCB+INTRPC)
    // Captures the ptolemy_splncb and ptolemy_intrpc lambdas defined inside the old inner block.
    // We implement directly here using the same algorithm.
    auto ptolemy_spline = [&](const std::vector<double>& X_in,
                               const std::vector<double>& Y_in,
                               const std::vector<double>& X_out,
                               std::vector<double>& Y_out) {
      int N = (int)X_in.size();
      int M = (int)X_out.size();
      Y_out.assign(M, 0.0);
      if (N < 2) return;
      // Natural cubic spline (not-a-knot from Ptolemy SPLNCB with leading/trailing trim)
      // Simplified version: standard natural cubic spline
      std::vector<double> h_sp(N-1), alpha(N), l(N), mu(N), z(N);
      for (int i=0;i<N-1;i++) h_sp[i]=X_in[i+1]-X_in[i];
      for (int i=1;i<N-1;i++)
        alpha[i]=3.0/h_sp[i]*(Y_in[i+1]-Y_in[i])-3.0/h_sp[i-1]*(Y_in[i]-Y_in[i-1]);
      l[0]=1; mu[0]=0; z[0]=0;
      for (int i=1;i<N-1;i++){
        l[i]=2*(X_in[i+1]-X_in[i-1])-h_sp[i-1]*mu[i-1];
        if(std::abs(l[i])<1e-30) l[i]=1e-30;
        mu[i]=h_sp[i]/l[i]; z[i]=(alpha[i]-h_sp[i-1]*z[i-1])/l[i];
      }
      l[N-1]=1; z[N-1]=0;
      std::vector<double> c(N,0),b(N,0),d(N,0);
      for (int j=N-2;j>=0;j--){
        c[j]=z[j]-mu[j]*c[j+1];
        b[j]=(Y_in[j+1]-Y_in[j])/h_sp[j]-h_sp[j]*(c[j+1]+2*c[j])/3.0;
        d[j]=(c[j+1]-c[j])/(3.0*h_sp[j]);
      }
      // Evaluate at X_out points
      for (int k=0;k<M;k++){
        double x=X_out[k];
        // find segment
        int seg=0;
        for (int i=0;i<N-1;i++) if(x>=X_in[i]) seg=i;
        double dx=x-X_in[seg];
        Y_out[k]=Y_in[seg]+b[seg]*dx+c[seg]*dx*dx+d[seg]*dx*dx*dx;
      }
    };

    // ── FAITHFUL PTOLEMY ARCHITECTURE ────────────────────────────────────────────
    // Ptolemy INELDC loop structure (source.mor lines 17870-18200):
    //   1. Build full (Lo,Lx) list for this Li → IHMAX pairs
    //   2. Precompute chi_b for all (Lo,JPO) under this Li
    //   3. DO IV=1,NPDIF (V-slice):
    //        DO IU=1,NPSUM (H-computation U-grid):
    //          DO II=1,NPPHI (phi):
    //            Accumulate H[IH] for ALL IH (Lo,Lx) simultaneously
    //          SMHVL[IH][IU] = H[IH] * RIOEX
    //        Spline SMHVL[IH][*] → SMIVL[IH][*] for each IH
    //        DO IU=1,NPSUMI (chi-integration):
    //          FOR each (JPI,JPO,IH): I(JPI,JPO,Lo,Lx) += SMIVL[IH][IU]*chi_a*chi_b*TERM
    //
    // KEY: H is computed once per (Li,IV,IU) for ALL Lo simultaneously.
    //      GRDSET (SUMMIN/SUMMID/PHI0) is per-Li and already computed above.
    // ─────────────────────────────────────────────────────────────────────────────

    // ── Step 1: Build (Lo,Lx) pair list for this Li ──────────────────────────────
    struct LoLxPair { int Lo, Lx; };
    std::vector<LoLxPair> lolx_list;
    {
      for (int Lx = LxMin; Lx <= LxMax_bs; Lx += 2) {
        for (int Lo = 0; Lo <= Lmax; Lo++) {
          if (Lx < std::abs(Li - Lo)) continue;
          if (Lx > Li + Lo) continue;
          if ((Li + Lo + Lx) % 2 != 0) continue;
          lolx_list.push_back({Lo, Lx});
        }
      }
    }
    int IHMAX = (int)lolx_list.size();
    if (IHMAX == 0) continue;  // no valid (Lo,Lx) for this Li

    // ── Step 2: Precompute chi_b for all unique (Lo,JPO) ─────────────────────────
    // Map (Lo,JPO) → chi_b wavefunction
    struct LoJpo { int Lo, JPO;
      bool operator<(const LoJpo& o) const {
        return Lo!=o.Lo ? Lo<o.Lo : JPO<o.JPO; } };
    std::map<LoJpo, std::vector<std::complex<double>>> chi_b_cache;
    {
      std::set<int> unique_Lo;
      for (auto& [Lo, Lx] : lolx_list) unique_Lo.insert(Lo);
      for (int Lo : unique_Lo) {
        int JPO_min = std::max(1, std::abs(2*Lo - JB_dw));
        int JPO_max = 2*Lo + JB_dw;
        for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
          WavElj(Outgoing, Lo, JPO);
          chi_b_cache[{Lo, JPO}] = Outgoing.WaveFunction;
        }
      }
    }
    const double h_chi_b = Outgoing.StepSize;
    const double chi_b_maxR_all = (Outgoing.NSteps > 4)
        ? (Outgoing.NSteps - 4) * h_chi_b : Outgoing.MaxR;

    // ── Step 3: Precompute A12 terms for all IHMAX (Lo,Lx) pairs ─────────────────
    std::vector<std::vector<std::tuple<int,int,double>>> A12_all(IHMAX);
    for (int IH = 0; IH < IHMAX; ++IH) {
      A12_all[IH] = ComputeA12Terms(Li, lolx_list[IH].Lo, lolx_list[IH].Lx, lT, lP);
      if (Li == 3 && IH == 0) {
        fprintf(stderr, "CPP_A12_TERMS Li=%d IH=%d Lo=%d Lx=%d lT=%d lP=%d:\n",
                Li, IH, lolx_list[IH].Lo, lolx_list[IH].Lx, lT, lP);
        for (auto &[MT_k, MU_k, coeff] : A12_all[IH])
          fprintf(stderr, "  MT=%2d MU=%2d coeff=%+.8f\n", MT_k, MU_k, coeff);
      }
    }

    // ── Quadrature parameters ─────────────────────────────────────────────────────
    const int NPSUM = 40;
    const int NPDIF = 40;
    const int NPPHI = 20;
    const double SUMMAX = 30.4;  // Ptolemy SUMMAX (overridden per Li below for high Li)
    const double GAMSUM = 2.0;
    const double GAMDIF = 12.0;
    const double GAMPHI = 1.0e-6;
    const double PHIMID_frac = 0.20;
    const double DWCUT = 2.0e-6;
    const int LOOKST = 250;
    const int NPHIAD = 4;
    const double DXV = 2.0 / ((double)LOOKST * LOOKST);
    double AKI = Incoming.k;
    double AKO = Outgoing.k;

    // ── SUMMAX varies at high Li (from Ptolemy per-Li runs) ─────────────────────
    static const double SUMMAX_table2[] = {
      30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4,
      30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.4, 30.8, 31.4,
      31.8, 31.8, 31.8, 31.8, 31.8, 31.8, 31.8, 31.8, 31.8, 31.8, 31.8
    };
    const double SUMMAX_eff = SUMMAX_table2[std::min(Li, 30)];

    // Phi GL base points [0,1]
    std::vector<double> phi_pts(NPPHI), phi_wts(NPPHI);
    GaussLegendre(NPPHI, -1.0, 1.0, phi_pts, phi_wts);
    CubMap(2, 0.0, PHIMID_frac, 1.0, GAMPHI, phi_pts, phi_wts);

    // H-computation U-grid (NPSUM)
    std::vector<double> xi_s(NPSUM), wi_s(NPSUM);
    GaussLegendre(NPSUM, -1.0, 1.0, xi_s, wi_s);
    CubMap(2, SUMMIN_li, SUMMID_li, SUMMAX_eff, GAMSUM, xi_s, wi_s);

    // Chi-integration U-grid (NPSUMI)
    int NPSUMI = (int)((SUMMAX_eff - SUMMIN_li) * 8.0 * (AKI + AKO) / (4.0*M_PI));
    if (NPSUMI < NPSUM) NPSUMI = NPSUM;
    std::vector<double> xi_si(NPSUMI), wi_si(NPSUMI);
    GaussLegendre(NPSUMI, -1.0, 1.0, xi_si, wi_si);
    CubMap(2, SUMMIN_li, SUMMID_li, SUMMAX_eff, GAMSUM, xi_si, wi_si);

    // Chi_a interpolation helper
    const double h_chi_a = Incoming.StepSize;

    // ── GRDSET DO 489: per-IU V-range scan (faithful Ptolemy replica) ────────────
    // For each U on the H-grid, scan BSPROD ITYPE=2 (phi'*V*phi' only, no chi)
    // to find VMAX and VMIN fractions in [0,1] (relative to VLEN=2U).
    // For U < 1 fm: skip scan; use VMIN=VMAX=1 → VABSMIN=VABSMAX=2U (clip to ±2U).
    // DV = 1/LOOKST step size for the scan.
    // RVRLIM_li = DWCUT * WVWMAX (precomputed above).
    //
    // bsprod2(U, vfrac, xs): phi'_T(rx) * phi'_P_vphi(rp) at (ra=U+vfrac*U, rb=U-vfrac*U, x=xs)
    // Uses XS(1)=1.0, XS(2)=1-DXV ≈ 0.999968 (Ptolemy GRDSET line 16209)
    const double XS1 = 1.0;
    const double XS2 = 1.0 - DXV;

    auto bsprod2 = [&](double ra, double rb, double xs) -> double {
      // BSPROD ITYPE=2: PHI'_T(rx) * PHI'_P_vphi(rp)  (flat-top bound states)
      // PHI'(r) = phi(r) for r >= r_peak; phi_max for r < r_peak
      double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*xs;
      double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*xs;
      if (rx2 < 0) rx2 = 0; if (rp2 < 0) rp2 = 0;
      double rx = std::sqrt(rx2), rp = std::sqrt(rp2);
      // phi_T flat-top (ITYPE=2: clipped)
      double phi_T_val;
      {
        double rv = rx / h_T;
        int ii = (int)rv; if (ii >= NSteps_T-1) { phi_T_val = 0.0; }
        else {
          double frac = rv - ii;
          double raw = TgtBS_ch.WaveFunction[ii].real()*(1-frac)
                     + TgtBS_ch.WaveFunction[ii+1].real()*frac;
          phi_T_val = (rx < r_T_peak) ? std::abs(TgtBS_ch.WaveFunction[idx_T_peak].real()) : raw;
        }
      }
      // IVPHI_P flat-top (ITYPE=2: clipped at peak)
      double ivphi_val;
      {
        double rv = rp / h_common;
        int ii = (int)rv; if (ii >= NSteps_common-1) { ivphi_val = 0.0; }
        else {
          double frac = rv - ii;
          double raw = IVPHI_P[ii]*(1-frac) + IVPHI_P[ii+1]*frac;
          ivphi_val = (rp < r_P_peak) ? IVPHI_P_max : raw;
        }
      }
      return std::abs(phi_T_val * ivphi_val);
    };

    // Faithful Ptolemy GRDSET DO 489 V-range scan as a reusable lambda.
    // Scans downward from VVAL=min(1, prev+3*DV) to find BSPROD cutoff.
    // For U < 1: Ptolemy skips scan (GO TO 465) and uses VMIN/VMAX from previous IU.
    // vmax_prev/vmin_prev track across IU iterations.
    const double DV_scan = 1.0 / (double)LOOKST;
    auto scan_vrange = [&](double U, double& vmax_prev, double& vmin_prev,
                            double& vmax_out, double& vmin_out) {
      if (U < 1.0) {
        // Ptolemy GO TO 465: inherit VMAX/VMIN from previous IU (already in vmax_prev/vmin_prev)
        vmax_out = vmax_prev;
        vmin_out = vmin_prev;
        return;
      }
      // Scan VMAX (positive V side): RI=U+VVAL*U, RO=U-VVAL*U
      double VVAL = std::min(1.0, vmax_prev + 3.0*DV_scan);
      vmax_out = DV_scan;
      while (true) {
        if (VVAL <= 0.5*DV_scan) { vmax_out = DV_scan; break; }
        double RI = U + VVAL * U;
        double RO = U - VVAL * U;
        if (RO < 0) RO = 0;
        double ULIM_loc = RVRLIM_li / std::max(1.0e-2, RI * RO);
        double f1 = std::abs(bsprod2(RI, RO, XS1));
        double f2 = std::abs(bsprod2(RI, RO, XS2));
        if (f1 > ULIM_loc || f2 > ULIM_loc) { vmax_out = std::min(1.0, VVAL); break; }
        VVAL -= DV_scan;
      }
      // Scan VMIN (negative V side): RI=U-VVAL*U, RO=U+VVAL*U
      VVAL = std::min(1.0, vmin_prev + 3.0*DV_scan);
      vmin_out = DV_scan;
      while (true) {
        if (VVAL <= 0.5*DV_scan) { vmin_out = DV_scan; break; }
        double RI = U - VVAL * U;
        double RO = U + VVAL * U;
        if (RI < 0) RI = 0;
        double ULIM_loc = RVRLIM_li / std::max(1.0e-2, RI * RO);
        double f1 = std::abs(bsprod2(RI, RO, XS1));
        double f2 = std::abs(bsprod2(RI, RO, XS2));
        if (f1 > ULIM_loc || f2 > ULIM_loc) { vmin_out = std::min(1.0, VVAL); break; }
        VVAL -= DV_scan;
      }
      vmax_prev = vmax_out;
      vmin_prev = vmin_out;
      if (Li==3 && U >= 1.0 && U <= 1.10) {
        fprintf(stderr, "SCAN_IU Li3 U=%.4f  vmax_out=%.4f vmin_out=%.4f\n", U, vmax_out, vmin_out);
        // Show BSPROD at VVAL=1.0 (full range)
        double RI_t = U + 1.0*U, RO_t = U - 1.0*U;
        double ULIM_t = RVRLIM_li / std::max(1.0e-2, RI_t * RO_t);
        double f1 = std::abs(bsprod2(RI_t, RO_t, XS1));
        fprintf(stderr, " at VVAL=1.0: RI=%.4f RO=%.4f ULIM=%.4e BSPROD2=%.4e (above=%s)\n",
                RI_t, RO_t, ULIM_t, f1, (f1>ULIM_t)?"YES":"NO");
      }
    };

    // Per-IU VMIN/VMAX fractions [0,1] for H-grid (NPSUM)
    std::vector<double> vmin_frac(NPSUM, 1.0), vmax_frac(NPSUM, 1.0);
    {
      double vmax_prev = 1.0, vmin_prev = 1.0;
      for (int IU = 0; IU < NPSUM; ++IU) {
        scan_vrange(xi_s[IU], vmax_prev, vmin_prev, vmax_frac[IU], vmin_frac[IU]);
      }
    }

    // Per-IU VMIN/VMAX absolute for chi-grid (NPSUMI)
    // Ptolemy uses polynomial interpolation (LSQPOL) from H-grid to chi-grid.
    // We instead scan BSPROD directly at each chi-grid U (same RVRLIM/threshold).
    // This avoids polynomial fitting complexity and gives exact Ptolemy-equivalent ranges.
    std::vector<double> vmin_abs_chi(NPSUMI), vmax_abs_chi(NPSUMI);
    {
      double vmax_prev = 1.0, vmin_prev = 1.0;
      for (int IU = 0; IU < NPSUMI; ++IU) {
        double U = xi_si[IU];
        double vmax_f, vmin_f;
        scan_vrange(U, vmax_prev, vmin_prev, vmax_f, vmin_f);
        vmax_abs_chi[IU] = std::min(vmax_f * 2.0 * U,  2.0*U);   // +VABSMAX
        vmin_abs_chi[IU] = std::max(-vmin_f * 2.0 * U, -2.0*U);  // -VABSMIN (stored as negative)
      }
    }

    // ── Polynomial fit to V-ranges (Ptolemy GRDSET Stage 3, LSQPOL, NVPOLY=3) ────
    // Ptolemy fits cubic polynomial to LVMIN and LVMAX (absolute V) vs U (H-grid).
    // This allows evaluation at chi-grid (NPSUMI) points (NPLYSW=false path).
    // We fit vmin_abs and vmax_abs (absolute V boundaries) using least-squares cubic.
    // LVMIN[IU] = -vmin_frac[IU]*2U (negative abs boundary)
    // LVMAX[IU] = +vmax_frac[IU]*2U (positive abs boundary)
    // Polynomial: p(U) = c[0] + c[1]*U + c[2]*U^2 + c[3]*U^3  (Horner: c[3]+U*(c[2]+U*(c[1]+U*c[0])))
    // Ptolemy stores coefficients in Horner order starting from highest degree.
    // Weighted least squares cubic polynomial fit (Ptolemy LSQPOL replica)
    // Weight = 1/(VMAX-VMIN)^2 per Ptolemy ALLOC(LVWTS-1+IU)=1/(VMAX-VMIN)^2
    // vmin_abs_arr and vmax_abs_arr are the absolute V-boundaries at H-grid points
    auto fit_cubic = [&](const std::vector<double>& x, const std::vector<double>& y,
                          const std::vector<double>& vmin_arr, const std::vector<double>& vmax_arr,
                          std::vector<double>& c) {
      int N = x.size();
      std::vector<double> ATA(16, 0.0), ATy(4, 0.0);
      for (int i = 0; i < N; ++i) {
        double U = x[i];
        // Ptolemy weight = 1/(VMAX-VMIN)^2; VMAX-VMIN = vmax_abs - vmin_abs
        double width = vmax_arr[i] - vmin_arr[i];
        double w = (width > 1e-12) ? 1.0 / (width * width) : 0.0;
        double basis[4] = {1.0, U, U*U, U*U*U};
        for (int j = 0; j < 4; ++j) {
          ATy[j] += w * basis[j] * y[i];
          for (int k = 0; k < 4; ++k)
            ATA[j*4+k] += w * basis[j] * basis[k];
        }
      }
      // Solve via Gaussian elimination with partial pivoting
      std::vector<double> A(ATA), b(ATy);
      for (int col = 0; col < 4; ++col) {
        int pivot = col;
        for (int row = col+1; row < 4; ++row)
          if (std::abs(A[row*4+col]) > std::abs(A[pivot*4+col])) pivot = row;
        std::swap(b[col], b[pivot]);
        for (int k = 0; k < 4; ++k) std::swap(A[col*4+k], A[pivot*4+k]);
        double diag = A[col*4+col];
        if (std::abs(diag) < 1e-30) continue;
        for (int k = col; k < 4; ++k) A[col*4+k] /= diag;
        b[col] /= diag;
        for (int row = 0; row < 4; ++row) {
          if (row == col) continue;
          double fac = A[row*4+col];
          for (int k = col; k < 4; ++k) A[row*4+k] -= fac * A[col*4+k];
          b[row] -= fac * b[col];
        }
      }
      c = b;  // c[0]=const, c[1]=linear, c[2]=quad, c[3]=cubic
    };

    // Build absolute V-boundary arrays for H-grid
    std::vector<double> vmin_abs_h(NPSUM), vmax_abs_h(NPSUM);
    for (int IU = 0; IU < NPSUM; ++IU) {
      double U = xi_s[IU];
      vmin_abs_h[IU] = -vmin_frac[IU] * 2.0 * U;  // negative (LVMIN)
      vmax_abs_h[IU] = +vmax_frac[IU] * 2.0 * U;  // positive (LVMAX)
    }
    std::vector<double> poly_vmin(4), poly_vmax(4);
    fit_cubic(xi_s, vmin_abs_h, vmin_abs_h, vmax_abs_h, poly_vmin);  // weighted by 1/(VMAX-VMIN)^2
    fit_cubic(xi_s, vmax_abs_h, vmin_abs_h, vmax_abs_h, poly_vmax);

    // Helper: evaluate cubic polynomial at U
    auto eval_poly = [](const std::vector<double>& c, double U) {
      return c[0] + U*(c[1] + U*(c[2] + U*c[3]));
    };

    // Diagnostic: print H-grid raw values AND polynomial V-ranges for chi-grid
    if (Li == 3) {
      fprintf(stderr, "=== H-grid vmin/vmax_abs (before poly fit, Li=3) ===\n");
      for (int IU = 0; IU < NPSUM; ++IU) {
        double U = xi_s[IU];
        fprintf(stderr, "H_VR IU=%2d U=%7.4f  vmin_abs=%8.5f vmax_abs=%8.5f\n",
                IU+1, U, vmin_abs_h[IU], vmax_abs_h[IU]);
      }
    }
    // Diagnostic: print polynomial-interpolated V-ranges for chi-grid IU=1..10
    if (Li == 3) {
      fprintf(stderr, "=== Poly V-ranges for chi-grid (Li=3) ===\n");
      for (int IU = 0; IU < std::min(NPSUMI, 10); ++IU) {
        double U = xi_si[IU];
        double vmin_p = eval_poly(poly_vmin, U);
        double vmax_p = eval_poly(poly_vmax, U);
        // VMID: use midpoint (simplified; Ptolemy uses moment scan DO 559)
        double VMID_p = 0.5*(vmin_p + vmax_p);
        fprintf(stderr, "CPP_POLY_VR IU=%2d U=%7.4f  vmin=%.5f vmax=%.5f vmid=%.5f\n",
                IU+1, U, vmin_p, vmax_p, VMID_p);
      }
    }

    // ── Step 4: Integral accumulators: I[IH][JPI][JPO] ───────────────────────────
    struct IntKey { int IH, JPI, JPO;
      bool operator<(const IntKey& o) const {
        if (IH!=o.IH) return IH<o.IH;
        if (JPI!=o.JPI) return JPI<o.JPI;
        return JPO<o.JPO; } };
    std::map<IntKey, std::complex<double>> I_accum;

    // ── Step 5: MAIN IV LOOP (Ptolemy DO 859 IV=1,NPDIF) ─────────────────────────
    // Ptolemy: outer IV, inner IU. Each IU has its own VMIN/VMAX → per-IU CUBMAP.
    // The GL base points gl_v_frac are mapped per-IU via CUBMAP(VMIN_IU,VMID_IU,VMAX_IU).
    // We replicate this: for each IU, run CUBMAP to get NPDIF points, then pick IV-th one.

    // Precompute per-IU CUBMAP output (NPDIF points + weights) for H-grid
    // Stored as dif_pts_h[IU][IV], dif_wts_h[IU][IV]
    std::vector<std::vector<double>> dif_pts_h(NPSUM, std::vector<double>(NPDIF)),
                                     dif_wts_h(NPSUM, std::vector<double>(NPDIF));
    std::vector<double> RIOEX_h(NPSUM * NPDIF, 0.0);  // indexed [IV*NPSUM + IU]

    for (int IU = 0; IU < NPSUM; ++IU) {
      double U = xi_s[IU];
      // VABSMAX = vmax_frac * 2U;  VABSMIN = -vmin_frac * 2U
      double VABSMAX =  vmax_frac[IU] * 2.0 * U;
      double VABSMIN = -vmin_frac[IU] * 2.0 * U;
      // Clip to ±2U (Ptolemy DO 689 line 16410)
      VABSMAX = std::min(VABSMAX,  2.0*U);
      VABSMIN = std::max(VABSMIN, -2.0*U);
      if (VABSMAX <= VABSMIN + 1e-12) {
        // Degenerate range — fill with zero
        std::fill(dif_pts_h[IU].begin(), dif_pts_h[IU].end(), 0.0);
        std::fill(dif_wts_h[IU].begin(), dif_wts_h[IU].end(), 0.0);
        continue;
      }
      // VMID: first moment of BSPROD2 weight (simplified: use midpoint for now)
      // TODO: compute proper VMID from DO 559 moment scan if needed
      double VMID = 0.5 * (VABSMIN + VABSMAX);
      // Clamp VMID to [VABSMIN+0.3*width, VABSMAX-0.3*width]
      double width = VABSMAX - VABSMIN;
      VMID = std::max(VABSMIN + 0.3*width, std::min(VABSMAX - 0.3*width, VMID));

      // SYNE: if VMID <= midpoint → positive orientation (SYNE=+1)
      //       else swap sign and reflect (SYNE=-1)
      double mid = 0.5*(VABSMIN + VABSMAX);
      double vmin_c = VABSMIN, vmax_c = VABSMAX, vmid_c = VMID;
      double SYNE = 1.0;
      if (vmid_c > mid) {
        // Reflect
        SYNE = -1.0;
        vmid_c = -vmid_c;
        double tmp = vmax_c; vmax_c = -vmin_c; vmin_c = -tmp;
      }

      // CUBMAP: map NPDIF GL points to [vmin_c, vmid_c, vmax_c]
      std::vector<double> pts(NPDIF), wts(NPDIF);
      GaussLegendre(NPDIF, -1.0, 1.0, pts, wts);
      CubMap(1, vmin_c, vmid_c, vmax_c, GAMDIF, pts, wts);

      // Store in correct SYNE order (Ptolemy DO 689: IPLUNK = NPDIF-IV if SYNE<0)
      for (int IV = 0; IV < NPDIF; ++IV) {
        int IV_store = (SYNE < 0) ? (NPDIF - 1 - IV) : IV;
        double VVAL = pts[IV] * SYNE;     // actual V value (un-reflect)
        double WT   = wts[IV];
        double ra = U + 0.5*VVAL;
        double rb = U - 0.5*VVAL;
        dif_pts_h[IU][IV_store] = VVAL;
        dif_wts_h[IU][IV_store] = WT;
        // RIOEX = exp(+ALPHAP*RP + ALPHAT*RT) at (ra,rb)
        // Ptolemy GRDSET: RP = SQRT(1+(S1*RI+T1*RO)^2), RT = SQRT(1+(S2*RI+T2*RO)^2)
        // Note: S1/T1 → target coord (rx), S2/T2 → proj coord (rp)
        // Ptolemy's variable naming is inverted vs physics but we must match it exactly.
        double RP0 = std::sqrt(1.0 + (S1*ra+T1*rb)*(S1*ra+T1*rb));  // S1,T1 = target
        double RT0 = std::sqrt(1.0 + (S2*ra+T2*rb)*(S2*ra+T2*rb));  // S2,T2 = proj
        RIOEX_h[IV_store*NPSUM + IU] = std::exp(ALPHAP*RP0 + ALPHAT*RT0);
      }
    }

    // Precompute per-IU CUBMAP for chi-integration grid (NPSUMI × NPDIF)
    // Uses polynomial-interpolated VMIN/VMAX at NPSUMI U-points (simplified: same scan)
    std::vector<std::vector<double>> dif_pts_chi(NPSUMI, std::vector<double>(NPDIF)),
                                     dif_wts_chi(NPSUMI, std::vector<double>(NPDIF));
    std::vector<double> LWIO_chi(NPSUMI * NPDIF, 0.0);  // JACOB*RI*RO*WOW*DIFWT*exp_neg

    for (int IU = 0; IU < NPSUMI; ++IU) {
      double U = xi_si[IU];
      double WOW = wi_si[IU];  // Ptolemy SMIVL (U-grid weight after spline)
      // V-range at this U: use weighted cubic polynomial fit (Ptolemy GRDSET Stage 3 LSQPOL)
      // Weight = 1/(VMAX-VMIN)^2 emphasizes large-U points with smaller integration ranges.
      double VABSMIN = eval_poly(poly_vmin, U);  // negative value
      double VABSMAX = eval_poly(poly_vmax, U);  // positive value
      // Clip to physical limits ±2U
      VABSMAX = std::min(VABSMAX,  2.0*U);
      VABSMIN = std::max(VABSMIN, -2.0*U);
      if (VABSMAX <= VABSMIN + 1e-12) {
        std::fill(dif_pts_chi[IU].begin(), dif_pts_chi[IU].end(), 0.0);
        std::fill(dif_wts_chi[IU].begin(), dif_wts_chi[IU].end(), 0.0);
        continue;
      }
      double width = VABSMAX - VABSMIN;
      double VMID = std::max(VABSMIN + 0.3*width, std::min(VABSMAX - 0.3*width,
                             0.5*(VABSMIN+VABSMAX)));
      double mid = 0.5*(VABSMIN+VABSMAX);
      double vmin_c=VABSMIN, vmax_c=VABSMAX, vmid_c=VMID, SYNE=1.0;
      if (vmid_c > mid) {
        SYNE=-1.0; vmid_c=-vmid_c;
        double tmp=vmax_c; vmax_c=-vmin_c; vmin_c=-tmp;
      }
      if (Li == 3 && IU == 4)  // IU=5 (0-indexed=4)
        fprintf(stderr, "CPP_CHI_VRANGE IU=5 U=%.6g VABSMIN=%.6g VABSMAX=%.6g VMID=%.6g SYNE=%.1f vmin_c=%.6g vmid_c=%.6g vmax_c=%.6g\n",
                U, VABSMIN, VABSMAX, VMID, SYNE, vmin_c, vmid_c, vmax_c);
      std::vector<double> pts(NPDIF), wts(NPDIF);
      GaussLegendre(NPDIF, -1.0, 1.0, pts, wts);
      CubMap(1, vmin_c, vmid_c, vmax_c, GAMDIF, pts, wts);

      for (int IV = 0; IV < NPDIF; ++IV) {
        int IV_store = (SYNE < 0) ? (NPDIF-1-IV) : IV;
        double VVAL = pts[IV] * SYNE;
        double WT   = wts[IV];
        double ra = U + 0.5*VVAL;
        double rb = U - 0.5*VVAL;
        dif_pts_chi[IU][IV_store] = VVAL;
        dif_wts_chi[IU][IV_store] = WT;
        // LWIO = JACOB * RI * RO * WOW * DIFWT * exp(-ALPHAP*RP - ALPHAT*RT)
        // Ptolemy: RP uses S1/T1, RT uses S2/T2 (same convention as RIOEX)
        double RP_c = std::sqrt(1.0 + (S1*ra+T1*rb)*(S1*ra+T1*rb));
        double RT_c = std::sqrt(1.0 + (S2*ra+T2*rb)*(S2*ra+T2*rb));
        double JACOB_chi = S1*S1*S1;
        double LWIO_val = JACOB_chi * ra * rb * WOW * WT * std::exp(-ALPHAP*RP_c - ALPHAT*RT_c);
        LWIO_chi[IV_store*NPSUMI + IU] = LWIO_val;
        if (Li == 3 && IU == 4 && IV == 0) {  // IU=5,IV=1 in Fortran (0-indexed)
          fprintf(stderr, "CPP_LWIO IU=%d IV=%d\n", IU+1, IV+1);
          fprintf(stderr, " JACOB= %.8g\n", JACOB_chi);
          fprintf(stderr, "     U= %.8g\n", U);
          fprintf(stderr, "  VVAL= %.8g\n", VVAL);
          fprintf(stderr, "    RI= %.8g\n", ra);
          fprintf(stderr, "    RO= %.8g\n", rb);
          fprintf(stderr, "   WOW= %.8g\n", WOW);
          fprintf(stderr, " DIFWT= %.8g\n", WT);
          fprintf(stderr, "  TEMP= %.8g\n", std::exp(-ALPHAP*RP_c - ALPHAT*RT_c));
          fprintf(stderr, "  LWIO= %.8g\n", LWIO_val);
          fprintf(stderr, " ALPHAP=%.8g ALPHAT=%.8g RP=%.8g RT=%.8g\n",
                  ALPHAP, ALPHAT, RP_c, RT_c);
        }
      }
    }

    fprintf(stderr, "Li=%d  SUMMIN=%.2f  SUMMID=%.2f  SUMMAX=%.2f  NPSUMI=%d\n",
            Li, SUMMIN_li, SUMMID_li, SUMMAX_eff, NPSUMI);

    // ── Step 5: MAIN DOUBLE LOOP: outer IV (DO 859), inner IU (DO 549) ──────────
    // Ptolemy structure:
    //   DO 859 IV = 1, NPDIF:           ← outer V loop
    //     DO 549 IU = 1, NPSUM:         ← H-integral accumulation on H-grid
    //       SMHVL[IH][IU] += phi_loop contributions
    //   DO 609 IH: spline SMHVL → SMIVL  ← ONCE after IV loop
    //   DO 789 IU = 1, NPSUMI:           ← chi integral (ONCE, not per IV)
    //     I[IH,JPI,JPO] += TERM * SMIVL[IH][IU] * chi_b * chi_a
    //
    // Faithful Ptolemy structure:
    //   DO 859 IV:
    //     DO 549 IU: SMHVL[IH][IU] = H(IU,IV)*RIOEX(IV,IU)  ← SET (not accumulate)
    //     DO 609 IH: spline SMHVL → SMIVL
    //     DO 789 IU: I[IH,JPI,JPO] += LWIO(IV,IU)*SMIVL[IH][IU]*DW(IV,IU)
    // The outer IV sum in I_accum is done by looping and adding each IV's contribution.

    std::vector<std::vector<double>> SMHVL(IHMAX, std::vector<double>(NPSUM, 0.0));

    for (int IV = 0; IV < NPDIF; ++IV) {
      // ── H-computation IU loop (DO 549): SET SMHVL for this IV ─────────────────
      // Ptolemy: SMHVL[IH][IU] = H(IU,IV)*RIOEX(IV,IU) — NOT accumulated across IV
      for (int IH = 0; IH < IHMAX; ++IH)
        std::fill(SMHVL[IH].begin(), SMHVL[IH].end(), 0.0);

      for (int IU = 0; IU < NPSUM; ++IU) {
        double U = xi_s[IU];
        if (U < 1e-6) continue;

        double V    = dif_pts_h[IU][IV];
        double DIFWT = dif_wts_h[IU][IV];
        if (std::abs(DIFWT) < 1e-30) continue;

        double ra = U + 0.5 * V;
        double rb = U - 0.5 * V;
        if (ra < 1e-6 || rb < 1e-6) continue;

        // RIOEX: use precomputed value from per-IU CUBMAP
        double RIOEX = RIOEX_h[IV * NPSUM + IU];

        // ── PHI0 scan (Ptolemy GRDSET pass-1 phi cutoff, source.mor 16840-16865) ────
        // Scan X = cos(phi) from 1 (phi=0) toward -1 (phi=π).
        // Stop when BOTH current AND previous |bsprod| < ULIM (two consecutive points).
        // IEND = II-1 + NPHIAD (capped at IXTOPZ = LOOKST+1).
        // Ptolemy Fortran: ULIM = RVRLIM / (max(1,RI) * max(1,RO))
        double ULIM = RVRLIM_li / (std::max(1.0, ra) * std::max(1.0, rb));
        int IEND = LOOKST + 1;  // default: use full range (PHI0=π)
        double WOW = std::abs(bsprod_val(ra, rb, 1.0));  // II=1: always keep
        for (int II = 2; II <= LOOKST + 1; ++II) {
          double X = 1.0 - DXV * (double)(II-1)*(double)(II-1);
          if (X < -1.0) X = -1.0;
          double fifo = std::abs(bsprod_val(ra, rb, X));
          if (fifo < ULIM && WOW < ULIM) {
            // Both consecutive points below ULIM → cutoff here
            IEND = std::min(II - 1 + NPHIAD, LOOKST + 1);
            IEND = std::max(IEND, 2);
            break;
          }
          WOW = fifo;
        }
        double X0  = 1.0 - DXV * (double)(IEND-1)*(double)(IEND-1);
        X0 = std::max(-1.0, std::min(1.0, X0));
        double PHI0 = std::acos(X0);

        // ── Phi integral: accumulate H[IH] for ALL IH simultaneously ───────────
        // Ptolemy DO 489 II=1,NPPHI  — one loop, all IH at once (innermost = IH loop)
        std::vector<double> HINT(IHMAX, 0.0);   // Ptolemy SALLOC(LHSM1+IH)

        for (int k = 0; k < NPPHI; ++k) {
          double PHI  = PHI0 * phi_pts[k];
          double DPHI = PHI0 * phi_wts[k] * std::sin(PHI);
          if (DPHI == 0.0) continue;
          double X    = std::cos(PHI);

          // BSPROD at (ra, rb, PHI)
          double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*X;
          double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*X;
          if (rx2 < 0) rx2 = 0; if (rp2 < 0) rp2 = 0;
          double rx = std::sqrt(rx2), rp = std::sqrt(rp2);

          // Ptolemy BSPROD ITYPE=1 (Pass 2): ALLSW=TRUE → raw interpolation, NO flat-top.
          // Both FT and FP use raw aitlag for all r (no peak-clip).
          // phi_T_val = target BS WF at rx; ivphi_val = V_np*phi_P at rp.
          double phi_T_val = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                                TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                                TgtBS_ch.MaxR, rx, 0.0, 0.0);
          double ivphi_val = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                              h_common*NSteps_common, rp, 0.0, 0.0);
          double PVPDX = DPHI * phi_T_val * ivphi_val;
          if (std::abs(PVPDX) < 1e-30) continue;

          // PHIT = acos((T1*rb + S1*ra*cos(PHI)) / rx)
          double cos_phiT = (rx2 > 1e-30) ? (T1*rb + S1*ra*X) / rx : 1.0;
          cos_phiT = std::max(-1.0, std::min(1.0, cos_phiT));
          double PHIT = std::acos(cos_phiT);

          // DIAGNOSTIC: at Li=3, IV=0, IU=10 — print all phi points for IH=0
          if (Li==3 && IV==0 && IU==10) {
            double A12_val0 = EvalA12(A12_all[0], PHIT, PHI);
            fprintf(stderr,"CPP_PHILOOP k=%2d PHI=%.6f PHIT=%.6f DPHI=%.6e phi_T=%.6e ivphi=%.6e PVPDX=%.6e A12[0]=%.6e term=%.6e\n",
                    k+1, PHI, PHIT, DPHI, phi_T_val, ivphi_val, PVPDX, A12_val0, PVPDX*A12_val0);
          }

          // ── Innermost loop: IH (all (Lo,Lx) pairs) — Ptolemy DO 459/469 ──────
          for (int IH = 0; IH < IHMAX; ++IH) {
            double A12_val = EvalA12(A12_all[IH], PHIT, PHI);
            if (Li==1 && IV==0 && IU==1 && k<3)
              fprintf(stderr,"A12_DETAIL IH=%d IU=%d k=%d PHI=%.5f PHIT=%.5f A12=%.6e PVPDX=%.6e term=%.6e\n",
                      IH, IU+1, k+1, PHI, PHIT, A12_val, PVPDX, PVPDX*A12_val);
            HINT[IH] += PVPDX * A12_val;
          }
        }
        // End phi loop

        // SET SMHVL[IH][IU] = H[IH] * RIOEX (Ptolemy DO 509 — assignment, not accumulation)
        // Reset to zero above at top of IV loop, then set here for each IU
        for (int IH = 0; IH < IHMAX; ++IH)
          SMHVL[IH][IU] = HINT[IH] * RIOEX;
        // STEP_A equivalent: Li=3, IV=0 — print HINT[0..1], RIOEX, SMHVL[0..1]
        if (Li==3 && IV==0 && IHMAX>=2)
          fprintf(stderr,"CPP_STEP_A LI3 IU=%2d U=%12.5e HINT1=%12.5e HINT2=%12.5e RIOEX=%12.5e SMHVL1=%12.5e SMHVL2=%12.5e\n",
                  IU+1, U, HINT[0], HINT[1], RIOEX, SMHVL[0][IU], SMHVL[1][IU]);

        // DIAGNOSTIC: print H values for Li=1, IV=0 (most-negative V side), all IU
        // Ptolemy HGRD_LI1 prints at IV=1 (most-negative V) on the NPSUM grid
        // Our IV=0 is the most-negative-V end (SYNE positive)
        if (Li==1 && IV==0 && IU < NPSUM && IHMAX>=2) {
          fprintf(stderr,"CPP_H IU=%2d U=%.7f V=%.7f ra=%.6f rb=%.6f  H[0]=%.6e H[1]=%.6e RIOEX=%.4e PHI0=%.5f IEND=%d\n",
                  IU+1, U, V, ra, rb, HINT[0], HINT[1], RIOEX, PHI0, IEND);
        }

      }  // End H-computation IU loop (Ptolemy DO 549)

    // ── DO 609 + DO 789: spline cumulative SMHVL → SMIVL, then chi-integrate at this IV ──
    // (Ptolemy: this runs INSIDE DO 859, once per IV, using cumulative SMHVL up to this IV)

      // ── SPLNCB + INTRPC: spline SMHVL[IH][NPSUM] → SMIVL[IH][NPSUMI] ──────────
      // (Ptolemy DO 609 IH=1,IHMAX)
      std::vector<std::vector<double>> SMIVL(IHMAX, std::vector<double>(NPSUMI, 0.0));

      // Ptolemy spline: same validated ptolemy_spline lambda from earlier
      // Build a standalone spline lambda (reuse the same code as before)
      auto do_spline = [&](const std::vector<double>& xi_in,
                           const std::vector<double>& yi_in,
                           const std::vector<double>& xi_out,
                           std::vector<double>& yi_out) {
        // Thin wrapper: uses the ptolemy_splncb + ptolemy_intrpc helpers
        // which are defined earlier in this function scope.
        // We replicate the call pattern used in the old H_smhvl spline block.
        std::vector<double> ywork = yi_in; // copy
        int N = (int)xi_in.size();
        int M = (int)xi_out.size();
        yi_out.assign(M, 0.0);
        if (N < 2) return;

        // Natural cubic spline (Ptolemy SPLNCB+INTRPC via ptolemy_spline)
        // ptolemy_spline is already defined as a lambda earlier in InelDc() scope
        // but we're in a new block — call it via the capture.
        ptolemy_spline(xi_in, yi_in, xi_out, yi_out);
      };

      for (int IH = 0; IH < IHMAX; ++IH)
        do_spline(xi_s, SMHVL[IH], xi_si, SMIVL[IH]);

      // STEP_B equivalent: SMIVL[IH=0..1][IU=0..9] at Li=3, IV=0
      if (Li==3 && IV==0 && IHMAX>=2) {
        for (int iu=0; iu<std::min(10,NPSUMI); ++iu)
          fprintf(stderr,"CPP_STEP_B LI3 IV=0 IU=%3d Usi=%.4f SMIVL[1]=%.6e SMIVL[2]=%.6e\n",
                  iu+1, xi_si[iu], SMIVL[0][iu], SMIVL[1][iu]);
      }

      // ── Chi integration: Ptolemy DO 789 (IU), at current IV of DO 859 ──────────────
      // Uses current IV's LWIO_chi and (ra_chi, rb_chi) grid.
      // SMIVL is the cumulative spline up to this IV step.
      for (int IU = 0; IU < NPSUMI; ++IU) {
        double U_chi  = xi_si[IU];
        double WGT_U  = wi_si[IU];
        if (U_chi < 1e-6) continue;

        // LWIO(IV,IU): precomputed = JACOB*ra*rb*WOW*DIFWT*exp_neg
        double LWIO = LWIO_chi[IV * NPSUMI + IU];
        double TERM = LWIO;
        // STEP_C equivalent: LWIO/TERM at Li=3, IV=0, IU=0..4
        if (Li==3 && IV==0 && IU<5)
          fprintf(stderr,"CPP_STEP_C LI3 IV=0 IU=%2d IPLUNK=%5d TERM=%12.5e\n",
                  IU+1, IV*NPSUMI+IU+1, TERM);

        double ra_chi = U_chi + 0.5 * dif_pts_chi[IU][IV];
        double rb_chi = U_chi - 0.5 * dif_pts_chi[IU][IV];
        if (ra_chi < 1e-6 || rb_chi < 1e-6) continue;

        // Check DWMAX (Ptolemy line 18130: skip if all DW products small)
        // Compute max |chi_a * chi_b| across all (JPI,Lo,JPO)
        // For now: just check one representative chi_a
        // (Ptolemy uses the NDW stored DW products — we do per-JPI below)

        // Loop over all (JPI, JPO, IH) combinations
        for (auto& [JPI_key, chi_a] : chi_a_byJPI) {
          // chi_a at ra_chi
          double chi_a_maxR = (chi_a.size() >= 4)
              ? (static_cast<int>(chi_a.size()) - 4) * h_chi_a : Incoming.MaxR;
          std::complex<double> ca;
          {
            double r = ra_chi;
            if (r <= 0 || r > chi_a_maxR) { ca = {0,0}; }
            else {
#ifdef INTERP_PTOLEMY
              double rbyh = r / h_chi_a;
              int I = std::max(2, std::min((int)(rbyh+0.5), (int)chi_a.size()-4));
              double P = rbyh-I, PS=P*P;
              double X1=P*(PS-1.)/24.,X2=X1+X1,X3=X1*P,X4=X2+X2-.5*P,X5=X4*P;
              double C1=X3-X2,C5=X3+X2,C3=(X5-X3)*2.+1.,C2=X5-X4,C4=X5+X4;
              ca = C1*chi_a[I-2]-C2*chi_a[I-1]+C3*chi_a[I]-C4*chi_a[I+1]+C5*chi_a[I+2];
#else
              int ii=(int)(r/h_chi_a); if(ii>=(int)chi_a.size()-1){ca={0,0};}
              else { double f=r/h_chi_a-ii; ca=chi_a[ii]*(1-f)+chi_a[ii+1]*f; }
#endif
            }
          }

          for (int IH = 0; IH < IHMAX; ++IH) {
            int Lo = lolx_list[IH].Lo;

            int JPO_min = std::max(1, std::abs(2*Lo - JB_dw));
            int JPO_max = 2*Lo + JB_dw;
            for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
              auto cb_it = chi_b_cache.find({Lo, JPO});
              if (cb_it == chi_b_cache.end()) continue;
              const auto& chi_b = cb_it->second;

              // chi_b at rb_chi (5-pt Lagrange or linear)
              std::complex<double> cb;
              {
                double r = rb_chi;
                if (r <= 0 || r > chi_b_maxR_all) { cb = {0,0}; }
                else {
#ifdef INTERP_PTOLEMY
                  double rbyh = r / h_chi_b;
                  int I = std::max(2, std::min((int)(rbyh+0.5), (int)chi_b.size()-4));
                  double P=rbyh-I,PS=P*P;
                  double X1=P*(PS-1.)/24.,X2=X1+X1,X3=X1*P,X4=X2+X2-.5*P,X5=X4*P;
                  double C1=X3-X2,C5=X3+X2,C3=(X5-X3)*2.+1.,C2=X5-X4,C4=X5+X4;
                  cb = C1*chi_b[I-2]-C2*chi_b[I-1]+C3*chi_b[I]-C4*chi_b[I+1]+C5*chi_b[I+2];
#else
                  int ii=(int)(r/h_chi_b); if(ii>=(int)chi_b.size()-1){cb={0,0};}
                  else { double f=r/h_chi_b-ii; cb=chi_b[ii]*(1-f)+chi_b[ii+1]*f; }
#endif
                }
              }

              // Ptolemy INELDC line 18147: DWR = Re(chi_b)*Re(chi_a) - Im(chi_b)*Im(chi_a)
              //                             DWI = Re(chi_b)*Im(chi_a) + Im(chi_b)*Re(chi_a)
              // i.e. DW = chi_b * chi_a  (not conj)
              std::complex<double> DW = cb * ca;
              if (std::abs(DW) < 1e-30) continue;

              // Contribution: TERM * SMIVL[IH][IU] * DW
              I_accum[{IH, JPI_key, JPO}] += TERM * SMIVL[IH][IU] * DW;
            }
          }
        }
      }  // End chi IU loop (Ptolemy DO 789) — chi integration for this IV

      // STEP_E equivalent: print running I_accum[IH=0] after each IV (matches Ptolemy II=1)
      if (Li == 3 && !I_accum.empty()) {
        auto it0 = I_accum.begin();
        fprintf(stderr, "CPP_STEP_E LI3 after IV=%2d I_accum[0] re=%.6e im=%.6e\n",
                IV, it0->second.real(), it0->second.imag());
      }

  }  // End IV loop (Ptolemy DO 859) — SMHVL accumulated, chi integrated at each IV

    // ── Step 6: Apply SFROMI factors and store ────────────────────────────────────
    for (int IH = 0; IH < IHMAX; ++IH) {
      int Lo = lolx_list[IH].Lo;
      int Lx = lolx_list[IH].Lx;

      // ITEST phase (Ptolemy SFROMI): i^(Li+Lo+2*Lx+1)
      int ITEST_val = ((Li + Lo + 2*Lx + 1) % 4 + 4) % 4;
      std::complex<double> phase_factor;
      switch (ITEST_val) {
        case 0: phase_factor = { 1.0,  0.0}; break;
        case 1: phase_factor = { 0.0,  1.0}; break;
        case 2: phase_factor = {-1.0,  0.0}; break;
        default:phase_factor = { 0.0, -1.0}; break;
      }

      // ATERM (same formula as before — independent of JPI/JPO)
      double ATERM_val = 0.0;
      if (Lx >= std::abs(TargetBS.l - ProjectileBS.l) &&
          Lx <= TargetBS.l + ProjectileBS.l) {
        double jT_bs = TargetBS.j, jP_bs = ProjectileBS.j, jx = 0.5;
        double sj = SixJ((double)TargetBS.l, jT_bs, jx,
                         jP_bs, (double)ProjectileBS.l, (double)Lx);
        int twoj_sum = 2*TargetBS.l + (int)(2*jT_bs+0.5)
                     + 2*ProjectileBS.l + (int)(2*jP_bs+0.5);
        double sign_val = ((twoj_sum/2)%2==0) ? 1.0 : -1.0;
        double RACAH_val = sign_val * sj;
        int JBIGA = (int)std::round(2.0*SpinTarget);
        int JBIGB = (int)std::round(2.0*SpinResidual);
        double TEMP_aterm = std::sqrt((JBIGB+1.0)/(JBIGA+1.0));
        // SPAMP: Ptolemy's SPAMP = bare_amplitude (e.g. AV18=0.97069)
        // No extra sqrt(2S+1) factor needed — Ptolemy's SPAMP is already the bare value.
        // The sqrt(3) factor seen earlier was a red herring from RACAH sign confusion.
        // Verified by mini Fortran test: RACAH(4,5,0,1,1,4)=+0.31623 (not -0.18257)
        // ATERM = sqrt(6)*sqrt(5)*0.97069*0.31623 = +1.6813 before ITEST flip.
        double SPAMP = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
        double SPAMT = 1.0;  // target spectroscopic amplitude (default)
        ATERM_val = TEMP_aterm * std::sqrt(2.0*Lx+1.0) * SPAMP * SPAMT * RACAH_val;
        // Sign: Ptolemy BETCAL lines 29890-29892:
        //   ITEST = JX - JBP + 2*(LBP+LBT)
        //   ITEST = ITEST/2 + 1          ← integer divide then +1
        //   IF(MOD(ITEST,2) != 0) flip   ← check odd after transform
        int JX_d  = 1;  // 2*jx = 1
        int JBP_d = (int)std::round(2.0*ProjectileBS.j);
        int ITEST_a = JX_d - JBP_d + 2*(ProjectileBS.l + TargetBS.l);
        ITEST_a = ITEST_a/2 + 1;   // Ptolemy two-step transform
        if (ITEST_a % 2 != 0) ATERM_val = -ATERM_val;
        if (Li == 3) {
          fprintf(stderr, "ATERM_DBG Li=%d Lx=%d TEMP_aterm=%.5e SPAMP=%.5f RACAH=%.5e ATERM=%.5e\n",
                  Li, Lx, std::sqrt((JBIGB+1.0)/(JBIGA+1.0)) * std::sqrt(2.0*Lx+1.0),
                  ProjectileWFLoaded ? ProjectileWFSpam : 0.97069,
                  SixJ((double)TargetBS.l, TargetBS.j, 0.5, ProjectileBS.j,
                       (double)ProjectileBS.l, (double)Lx) *
                  (((2*TargetBS.l+(int)(2*TargetBS.j+0.5)+2*ProjectileBS.l+(int)(2*ProjectileBS.j+0.5))/2%2==0)?1.0:-1.0),
                  ATERM_val);
        }
      }

      double FACTOR_sfromi = 2.0 * std::sqrt(Incoming.k * Outgoing.k /
                                              (Incoming.Ecm * Outgoing.Ecm));

      for (auto& [JPI_key, chi_a] : chi_a_byJPI) {
        int JPO_min = std::max(1, std::abs(2*Lo - JB_dw));
        int JPO_max = 2*Lo + JB_dw;
        for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
          auto it = I_accum.find({IH, JPI_key, JPO});
          if (it == I_accum.end()) continue;

          std::complex<double> Integral = it->second * phase_factor;

          double sfromi_norm = FACTOR_sfromi * std::abs(ATERM_val) / std::sqrt(2.0*Li+1.0);
          if (ATERM_val < 0) Integral = -Integral;

          auto S_pre9j = Integral * sfromi_norm;
          TransferSMatrix.push_back({Lx, Li, Lo, JPI_key, JPO, S_pre9j});

          // DIAGNOSTIC: print raw I_accum for Li=3, Lo=3 to compare JPI dependence
          if (Li==3 && Lo==3 && (JPO==5 || JPO==7)) {
            fprintf(stderr, "IACCUM_DIAG Li=%d JPI=%d/2 Lo=%d JPO=%d/2 Lx=%d  "
                    "I_raw=(%10.4e,%10.4e) ATERM=%.5f phase=(%g,%g) S=(%10.4e,%10.4e)\n",
                    Li, JPI_key, Lo, JPO, Lx,
                    it->second.real(), it->second.imag(),
                    ATERM_val, phase_factor.real(), phase_factor.imag(),
                    S_pre9j.real(), S_pre9j.imag());
          }

#ifdef DEBUG_PRE9J
          fprintf(stderr, "PRE9J Li=%2d JPI=%2d/2  Lo=%2d JPO=%2d/2  Lx=%d  "
                  "S=(%9.4e, %9.4e)  |S|=%9.4e\n",
                  Li, JPI_key, Lo, JPO, Lx,
                  S_pre9j.real(), S_pre9j.imag(), std::abs(S_pre9j));
#endif
        }
      }
    }
    // ── End faithful INELDC body for this Li ──────────────────────────────────────
  } // end Li loop


}
