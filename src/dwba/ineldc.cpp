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
  int Lmax = 40;  // Match Ptolemy lmax=40

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

    // ── Per-Li GRDSET: compute SUMMIN, SUMMID, WVWMAX, RVRLIM ─────────────
    const double DWCUT_grdset = 2.0e-6;  // DWCUTOFF (same as DWCUT inside inner loop)
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
    double SUMMIN_li = 0.0, SUMMID_li = 0.0, WVWMAX_li = 0.0, RVRLIM_li = 0.0;
    {
      // Get representative chi waves for this Li
      const auto& ref_chi_a = chi_a_byJPI.begin()->second;
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

      // bsprod_with_ref_chi: |chi_a(U)| * bsprod(U,U,x) * |chi_b(U)| * U * U
      // (Ptolemy BSPROD ITYPE=3 at V=0: ra=rb=U)
      auto bsprod_ref = [&](double U, double x) -> double {
        if (U < 1e-6) return 0.0;
        double rx2 = S1*S1*U*U + T1*T1*U*U + 2.0*S1*T1*U*U*x;
        double rp2 = S2*S2*U*U + T2*T2*U*U + 2.0*S2*T2*U*U*x;
        if (rx2 < 0) rx2 = 0; if (rp2 < 0) rp2 = 0;
        double rx = std::sqrt(rx2), rp = std::sqrt(rp2);
        double phi_T = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                          TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                          TgtBS_ch.MaxR, rx, 0.0, 0.0);
        double ivphi_P = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                           h_common*NSteps_common, rp, 0.0, 0.0);
        double ca = std::abs(interp_ref_a(U));
        double cb = std::abs(interp_ref_b(U));
        return U * ca * std::abs(phi_T * ivphi_P) * cb * U;
      };

      // Step 1: WVWMAX — scan U from SUMMAX/2 down (Ptolemy line 15924)
      const double SUMMAX_li = 30.4;
      for (double U_s = 0.5*SUMMAX_li; U_s > 0.05; U_s -= 0.2)
        WVWMAX_li = std::max(WVWMAX_li, std::max(bsprod_ref(U_s,1.0), bsprod_ref(U_s,-1.0)));
      RVRLIM_li = DWCUT_grdset * std::max(WVWMAX_li, 1.0e-30);

      // Step 2: SUMMIN — step out from 0 (Ptolemy line 15964)
      double U_s = 0.0;
      while (U_s <= SUMMAX_li) {
        if (std::max(bsprod_ref(U_s,1.0), bsprod_ref(U_s,-1.0)) >= RVRLIM_li) {
          SUMMIN_li = std::max(0.0, U_s - 0.2);
          break;
        }
        U_s += 0.2;
      }

      // Step 3: SUMMID — first moment <U> * AMDMLT (Ptolemy line 16070)
      // Scan V=0 (diagonal) + 4 V fractions, 2 x values (Ptolemy scans 5 V's × 2 x's)
      double SUM0 = 0.0, SUM1 = 0.0;
      const double vfracs[5] = {-0.8, -0.4, 0.0, 0.4, 0.8};
      const double xscan[2] = {1.0 - 2.0/(100.0*100.0), -1.0 + 2.0/(100.0*100.0)};
      for (U_s = SUMMIN_li; U_s <= SUMMAX_li; U_s += 0.2) {
        double temp = 0.0;
        for (double vf : vfracs) {
          double ra_s = U_s*(1+vf), rb_s = U_s*(1-vf);
          if (ra_s < 1e-6 || rb_s < 1e-6) continue;
          for (double xs : xscan) {
            // Full (ra,rb) bsprod with chi (Ptolemy ITYPE=3)
            double rx2 = S1*S1*ra_s*ra_s + T1*T1*rb_s*rb_s + 2.0*S1*T1*ra_s*rb_s*xs;
            double rp2 = S2*S2*ra_s*ra_s + T2*T2*rb_s*rb_s + 2.0*S2*T2*ra_s*rb_s*xs;
            if (rx2<0) rx2=0; if (rp2<0) rp2=0;
            double rx=std::sqrt(rx2), rp=std::sqrt(rp2);
            double phi_T = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                              TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                              TgtBS_ch.MaxR, rx, 0.0, 0.0);
            double ivphi_P = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                               h_common*NSteps_common, rp, 0.0, 0.0);
            double ca = std::abs(interp_ref_a(ra_s));
            double cb = std::abs(interp_ref_b(rb_s));
            temp += ra_s * ca * std::abs(phi_T*ivphi_P) * cb * rb_s;
          }
        }
        SUM0 += temp; SUM1 += temp * U_s;
      }
      double mean_U = (SUM0 > 1e-30) ? SUM1/SUM0 : 0.5*(SUMMIN_li + SUMMAX_li);
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
    for (int IH = 0; IH < IHMAX; ++IH)
      A12_all[IH] = ComputeA12Terms(Li, lolx_list[IH].Lo, lolx_list[IH].Lx, lT, lP);

    // ── Quadrature parameters ─────────────────────────────────────────────────────
    const int NPSUM = 40;
    const int NPDIF = 40;
    const int NPPHI = 20;
    const double SUMMAX = 30.4;
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

    // NPSUMI (chi-integration grid)
    int NPSUMI = (int)((SUMMAX - SUMMIN_li) * 8.0 * (AKI + AKO) / (4.0*M_PI));
    if (NPSUMI < NPSUM) NPSUMI = NPSUM;

    // Phi GL base points [0,1] (CubMap concentrates near 0 per Ptolemy)
    std::vector<double> phi_pts(NPPHI), phi_wts(NPPHI);
    GaussLegendre(NPPHI, -1.0, 1.0, phi_pts, phi_wts);
    CubMap(2, 0.0, PHIMID_frac, 1.0, GAMPHI, phi_pts, phi_wts);

    // H-computation U-grid (NPSUM) and chi-integration U-grid (NPSUMI)
    std::vector<double> xi_s(NPSUM), wi_s(NPSUM);
    GaussLegendre(NPSUM, -1.0, 1.0, xi_s, wi_s);
    CubMap(2, SUMMIN_li, SUMMID_li, SUMMAX, GAMSUM, xi_s, wi_s);

    std::vector<double> xi_si(NPSUMI), wi_si(NPSUMI);
    GaussLegendre(NPSUMI, -1.0, 1.0, xi_si, wi_si);
    CubMap(2, SUMMIN_li, SUMMID_li, SUMMAX, GAMSUM, xi_si, wi_si);

    // V-grid GL base (NPDIF)
    std::vector<double> gl_dif_base(NPDIF), gl_dif_wts_base(NPDIF);
    GaussLegendre(NPDIF, -1.0, 1.0, gl_dif_base, gl_dif_wts_base);

    // Chi_a interpolation helper
    const double h_chi_a = Incoming.StepSize;

    // ── Step 4: Integral accumulators: I[IH][JPI][JPO] ───────────────────────────
    // Key: (IH index into lolx_list, JPI, JPO)
    struct IntKey { int IH, JPI, JPO;
      bool operator<(const IntKey& o) const {
        if (IH!=o.IH) return IH<o.IH;
        if (JPI!=o.JPI) return JPI<o.JPI;
        return JPO<o.JPO; } };
    std::map<IntKey, std::complex<double>> I_accum;

    // ── Step 5: MAIN IV LOOP (Ptolemy DO 859 IV=1,NPDIF) ─────────────────────────
    for (int IV = 0; IV < NPDIF; ++IV) {
      double gl_v_frac = gl_dif_base[IV];
      double gl_v_wt   = gl_dif_wts_base[IV];

      // ── H-computation IU loop (Ptolemy DO 549 IU=1,NPSUM) ──────────────────────
      // SMHVL[IH][IU] = H[IH] * RIOEX   (IHMAX × NPSUM array)
      std::vector<std::vector<double>> SMHVL(IHMAX, std::vector<double>(NPSUM, 0.0));

      for (int IU = 0; IU < NPSUM; ++IU) {
        double U = xi_s[IU];
        if (U < 1e-6) { continue; }

        // V range: flat ±2U (symmetric, Ptolemy IRECT=1 before LSQPOL trim)
        double VLEN = 2.0 * U;
        double V    = gl_v_frac * VLEN;
        double DIFWT = gl_v_wt * VLEN;

        double ra = U + V * 0.5;
        double rb = U - V * 0.5;
        if (ra < 1e-6 || rb < 1e-6) continue;

        // RIOEX = exp(+ALPHAP*RP + ALPHAT*RT) at phi=0 (x=1)
        // Ptolemy source.mor line 16443-16444:
        //   RP = SQRT( 1 + (S1*RI+T1*RO)**2 )   ← +1 inside sqrt (regularizes singularity)
        //   RT = SQRT( 1 + (S2*RI+T2*RO)**2 )
        // Ptolemy's S1,T1 → projectile coord; S2,T2 → target coord
        // Our code: S1,T1 → target (rx); S2,T2 → projectile (rp)
        // So: Ptolemy RP (projectile) = sqrt(1 + (S1_pto*RI+T1_pto*RO)^2) = sqrt(1 + (S2*ra+T2*rb)^2) [our S2,T2]
        //     Ptolemy RT (target)     = sqrt(1 + (S2_pto*RI+T2_pto*RO)^2) = sqrt(1 + (S1*ra+T1*rb)^2) [our S1,T1]
        // ALPHAP goes with RP (projectile BS), ALPHAT goes with RT (target BS) — unchanged
        double RP0 = std::sqrt(1.0 + (S2*ra + T2*rb)*(S2*ra + T2*rb));
        double RT0 = std::sqrt(1.0 + (S1*ra + T1*rb)*(S1*ra + T1*rb));
        double RIOEX = std::exp(ALPHAP * RP0 + ALPHAT * RT0);

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

          // ── Innermost loop: IH (all (Lo,Lx) pairs) — Ptolemy DO 459/469 ──────
          for (int IH = 0; IH < IHMAX; ++IH) {
            double A12_val = EvalA12(A12_all[IH], PHIT, PHI);
            HINT[IH] += PVPDX * A12_val;
          }
        }
        // End phi loop

        // Store SMHVL[IH][IU] = H[IH] * RIOEX (Ptolemy DO 509)
        for (int IH = 0; IH < IHMAX; ++IH)
          SMHVL[IH][IU] = HINT[IH] * RIOEX;

        // DIAGNOSTIC: print first 10 H values for Li=0, IV=0
        if (Li==0 && IV==0 && IU < 10 && IHMAX>0)
          fprintf(stderr,"DIAG_H IU=%2d U=%.4f ra=%.4f rb=%.4f  HINT[0]=%.4e RIOEX=%.4e SMHVL[0]=%.4e\n",
                  IU, U, ra, rb, HINT[0], RIOEX, SMHVL[0][IU]);

      }  // End H-computation IU loop (Ptolemy DO 549)

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

      // ── Chi integration: IU=1..NPSUMI (Ptolemy DO 789) ─────────────────────────
      // TERM = exp(-ALPHAP*RP - ALPHAT*RT) * DIFWT * wi_si[IU]  (Ptolemy LWIO * GL wt)
      // I(IH,JPI,JPO) += TERM * chi_a(ra) * SMIVL[IH][IU] * chi_b*(rb)
      for (int IU = 0; IU < NPSUMI; ++IU) {
        double U_chi  = xi_si[IU];
        double WGT_U  = wi_si[IU];
        if (U_chi < 1e-6) continue;

        // V at this chi-grid U: same IV fraction, recomputed at U_chi
        double VLEN_chi = 2.0 * U_chi;
        double V_chi    = gl_v_frac * VLEN_chi;
        double DIFWT    = gl_v_wt * VLEN_chi;

        double ra_chi = U_chi + V_chi * 0.5;
        double rb_chi = U_chi - V_chi * 0.5;
        if (ra_chi < 1e-6 || rb_chi < 1e-6) continue;

        // LWIO = JACOB * ra * rb * exp(-ALPHAP*RP - ALPHAT*RT) * DIFWT
        // (matches old code: TERM = JACOB_grdset * ra * rb * WOW * DIFWT * exp_neg)
        // JACOB_grdset = S1^3 (Ptolemy GRDSET line 15882: JACOB = S1**3 for stripping)
        // ra*rb comes from RIROWTS (Ptolemy chi-grid Jacobian factor)
        // Ptolemy source.mor line 16531-16532: same sqrt(1+x^2) formula
        double RP_chi = std::sqrt(1.0 + (S2*ra_chi + T2*rb_chi)*(S2*ra_chi + T2*rb_chi));
        double RT_chi = std::sqrt(1.0 + (S1*ra_chi + T1*rb_chi)*(S1*ra_chi + T1*rb_chi));
        double JACOB_chi = S1*S1*S1;  // S1^3 (stripping Jacobian)
        double LWIO = JACOB_chi * ra_chi * rb_chi * std::exp(-ALPHAP * RP_chi - ALPHAT * RT_chi) * DIFWT;
        double TERM = LWIO * WGT_U;

        // DIAGNOSTIC: print chi-grid TERM + SMIVL for Li=0, IV=19, IU=9
        if (Li==0 && IV==19 && IU==9)
          fprintf(stderr,"DIAG_CHI IV=19 IU=9 U_chi=%.4f ra_chi=%.4f rb_chi=%.4f  DIFWT=%.4e LWIO=%.4e TERM=%.4e SMIVL[0]=%.4e\n",
                  U_chi, ra_chi, rb_chi, DIFWT, LWIO, TERM, SMIVL[0][IU]);

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
      }  // End chi IU loop (Ptolemy DO 789)

    }  // End IV loop (Ptolemy DO 859)

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
        double SPAMP = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
        ATERM_val = TEMP_aterm * std::sqrt(2.0*Lx+1.0) * SPAMP * 1.0 * RACAH_val;
        int JX_d=(int)1, JBP_d=(int)(2*ProjectileBS.j);
        int ITEST_a = JX_d - JBP_d + 2*(ProjectileBS.l + TargetBS.l);
        if ((ITEST_a/2+1)%2!=0) ATERM_val = -ATERM_val;
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
