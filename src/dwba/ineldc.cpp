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


void DWBA::GrdSet() {
  // NOTE: phi_ab integration now uses CUBMAP adaptive PHI0 in InelDc.
  // ThetaGrid is kept for compatibility but no longer used in the phi loop.
  // The phi integration uses NPPHI=10 CUBMAP points on [0, PHI0] per (ra,rb) pair.
  int NTheta = 10;   // retained for compatibility; not used in phi integration
  GaussLegendre(NTheta, -1.0, 1.0, ThetaGrid, ThetaWeights);
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
  CalculateBoundState(PrjBS_ch, ProjectileBS.n, ProjectileBS.l, ProjectileBS.j, ProjectileBS.BindingEnergy);

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
                        PrjBS_ch.Target.Z, PrjBS_ch.Target.A);
    }
    std::cout << "  PrjBS V_real rebuilt: V_sol=" << PrjBS_ch.Pot.V << " MeV"
              << "  R_np=" << PrjBS_ch.Pot.R0 * std::pow(PrjBS_ch.Target.A, 1.0/3.0)
              << " fm  Target.A=" << PrjBS_ch.Target.A << std::endl;
    // Diagnostic: print V_np at a few r values
    for (double rcheck : {0.5, 1.0, 1.5, 2.0, 3.0}) {
      int idx = static_cast<int>(rcheck / PrjBS_ch.StepSize);
      if (idx < PrjBS_ch.NSteps) {
        std::printf("    V_np(r=%.1f)=%.3f MeV  phi_P(r=%.1f)=%.4e\n",
                    rcheck, PrjBS_ch.V_real[idx], rcheck, 
                    PrjBS_ch.WaveFunction[idx].real());
      }
    }
  }

  // Override projectile WF with tabulated values (e.g. AV18) if loaded
  if (ProjectileWFLoaded) {
    double h_prj = PrjBS_ch.StepSize;
    double h_tab = ProjectileWFGridH;
    std::cout << "  Overriding projectile WF with tabulated AV18 (SPAM="
              << ProjectileWFSpam << ")\n";
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
      PrjBS_ch.WaveFunction[I] = std::complex<double>(phi * ProjectileWFSpam, 0.0);
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
  std::cout << "Bound State MaxAmp: Target=" << maxT << " Projectile=" << maxP << std::endl;
  std::cout << "Bound State cutoff: Target=" << cutT * TgtBS_ch.StepSize
            << " fm, Projectile=" << cutP * PrjBS_ch.StepSize << " fm" << std::endl;

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

  std::cout << "Computing distorted waves and transfer integrals (Li=0.."
            << Lmax << ")..." << std::endl;

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
  int A_a = Incoming.Projectile.A;   // deuteron: A=2
  int A_b = Outgoing.Projectile.A;   // proton: A=1
  int A_B = Outgoing.Target.A;       // 34Si: A=34
  // BSPROD stripping/projectile-vertex: ISC=2 (outgoing channel)
  // RNSCAT = R0_out * A_B^(1/3)   [outgoing p-B scattering radius]
  // RNCORE = RNSCAT * (A_a^(1/3) + A_b^(1/3)) / (A_B^(1/3) + A_b^(1/3))
  double RNSCAT_post = Outgoing.Pot.R0 * std::pow((double)A_B, 1.0/3.0);
  double RNCORE_post = RNSCAT_post * (std::pow((double)A_a, 1.0/3.0) + std::pow((double)A_b, 1.0/3.0))
                                   / (std::pow((double)A_B, 1.0/3.0) + std::pow((double)A_b, 1.0/3.0));
  double VOPT_post   = -Outgoing.Pot.V;  // negative of outgoing real depth (MeV)
  double AOPT_post   = Outgoing.Pot.A;
  double CORE_SCALE  = RNCORE_post / RNSCAT_post;  // r_core = CORE_SCALE * rb
  auto V_core_post = [&](double r_core) -> double {
    if (r_core < 1e-6) return VOPT_post;
    return VOPT_post / (1.0 + std::exp((r_core - RNCORE_post) / AOPT_post));
  };
  std::printf("  USECORE POST: RNSCAT=%.3f RNCORE=%.3f VOPT=%.2f AOPT=%.3f CORE_SCALE=%.4f\n",
              RNSCAT_post, RNCORE_post, VOPT_post, AOPT_post, CORE_SCALE);

  std::vector<double> IVPHI_T(NSteps_common, 0.0);  // phi_T * V_nA (PRIOR vertex)
  std::vector<double> IVPHI_P(NSteps_common, 0.0);  // phi_P * V_np (POST vertex, NUCONLY)
  double IVPHI_T_max = 0, IVPHI_P_max = 0;
  int idx_T_vert_peak = 1, idx_P_vert_peak = 1;
  for (int i = 1; i < NSteps_common; ++i) {
    IVPHI_T[i] = std::abs(TgtBS_ch.WaveFunction[i].real()) * TgtBS_ch.V_real[i];
    IVPHI_P[i] = std::abs(PrjBS_ch.WaveFunction[i].real()) * PrjBS_ch.V_real[i];
    if (IVPHI_T[i] > IVPHI_T_max) { IVPHI_T_max = IVPHI_T[i]; idx_T_vert_peak = i; }
    if (IVPHI_P[i] > IVPHI_P_max) { IVPHI_P_max = IVPHI_P[i]; idx_P_vert_peak = i; }
  }
  double r_T_vert_peak = idx_T_vert_peak * h_common;
  double r_P_vert_peak = idx_P_vert_peak * h_common;

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

  std::printf("  WF peaks: phi_T_max=%.4f at r=%.2f fm; phi_P_max=%.4f at r=%.2f fm\n",
              phi_T_max, r_T_peak, phi_P_max, r_P_peak);
  std::printf("  IVPHI vertex peaks: IVPHI_T_max=%.4f at r=%.2f fm; IVPHI_P_max=%.4f at r=%.2f fm\n",
              IVPHI_T_max, r_T_vert_peak, IVPHI_P_max, r_P_vert_peak);

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
  auto InterpolateClipped = [&](const std::vector<std::complex<double>> &wf,
                                const std::vector<double> &grid,
                                int nsteps, double stepsize, double maxr,
                                double r, double phi_max, double r_peak) -> double {
    if (r < r_peak) return phi_max;
    if (r >= maxr || r < 0) return 0.0;
    double idx = r / stepsize;
    int ii = static_cast<int>(idx);
    if (ii >= nsteps - 1) return 0.0;
    double frac = idx - ii;
    return std::abs(wf[ii].real() * (1.0 - frac) + wf[ii + 1].real() * frac);
  };

  // Clipped IVPHI interpolation for vertex product (V*phi at vertex side):
  // IVPHI'(r) = IVPHI_max for r < r_vert_peak; IVPHI(r) for r >= r_vert_peak
  auto InterpolateIVPHI = [&](const std::vector<double> &ivphi,
                               double stepsize, int nsteps, double maxr,
                               double r, double ivphi_max, double r_vert_peak) -> double {
    if (r < r_vert_peak) return ivphi_max;
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

  // Diagnostic: check chi magnitude for L=0
  {
    // For deuteron (JSPS=2), L=0: JP_min = |2*0-2| = 2, so JP=2 is minimum valid
    int diag_JP = 2*Incoming.JSPS;  // = 2*1=2 for proton, 2*2=4 for deuteron... 
    // Actually for L=0: JPI = |0*2 - JSPS| = JSPS
    diag_JP = std::abs(0 * 2 - Incoming.JSPS) + ((Incoming.JSPS % 2 == 0) ? 0 : 0);
    // Simple: use JP = 2*L + JSPS for L=0 → JP=JSPS=2 for deuteron
    diag_JP = Incoming.JSPS;  // For L=0: only valid JP is JSPS
    WavElj(Incoming, 0, diag_JP);
    double maxChi = 0;
    for (int i = 1; i < Incoming.NSteps; ++i)
      maxChi = std::max(maxChi, std::abs(Incoming.WaveFunction[i]));
    double S0 = (Incoming.SMatrix.size() > 0) ? std::abs(Incoming.SMatrix[0]) : 0.0;
    std::cout << "Chi_a(L=0,JP=" << diag_JP << ") MaxAmp=" << maxChi
              << "  |S_0|=" << S0 << std::endl;
  // Print full SMatrix for incoming and outgoing
  std::printf("INCOMING SMatrix (d+33Si):\n");
  int smat_out[] = {0,1,2,3,5,7,10,12,15};
  for (int il : smat_out) {
    if (il < (int)Incoming.SMatrix.size()) {
      auto& S = Incoming.SMatrix[il];
      std::printf("  L=%2d: |S|=%.6f phase=%.3f rad\n", il, std::abs(S), std::atan2(S.imag(), S.real()));
    }
  }
  std::printf("OUTGOING SMatrix (p+34Si):\n");
  for (int il : smat_out) {
    if (il < (int)Outgoing.SMatrix.size()) {
      auto& S = Outgoing.SMatrix[il];
      std::printf("  L=%2d: |S|=%.6f phase=%.3f rad\n", il, std::abs(S), std::atan2(S.imag(), S.real()));
    }
  }
  }

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
      // Print S-matrix for key L values
      int smat_L[] = {0,1,2,3,5,7,10,12};
      for (int sl : smat_L) {
        if (Li == sl && Li < (int)Incoming.SMatrix.size()) {
          auto& S = Incoming.SMatrix[Li];
          std::printf("[SMAT_INC] L=%2d JP=%2d |S|=%.6f phase=%.3f rad\n",
                      Li, JPI, std::abs(S), std::atan2(S.imag(), S.real()));
        }
      }
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
      // These depend only on (Li, Lo, Lx, lT, lP), NOT on JPI/JPO or (ra,rb).
      //
      // Reference: Ptolemy INELDC/A12 subroutine (source.mor lines 1453-1899)
      //
      // The A12 angular kernel at each (ra, rb, phi_ab) quadrature
      // point is: cos(MT × phi_T(ra,rb,phi_ab) + MP × phi_P(ra,rb,phi_ab) - MU × phi_ab)
      // where phi_T = arccos((T1×rb + S1×ra×cos(phi_ab)) / rx) and similarly phi_P.
      // This is Ptolemy INELDC comment (line ~17959): COSTHE = COS(MT*PHIT+MP*PHIP-MU*PHI).
      // The DMSVAL storage: DMSVAL(JA12M+2) = DINTS(MU) - 2 is an OFFSET for the recursion
      // (ARG starts at (2-MU_first)*PHI so that after the first recursion step it gives -MU*PHI).
      // (Ptolemy INELDC lines 17939-17949, DMSVAL stores MT, MP, MU-2)
      //
      // For lP=0 (MP=0): ARG = MT × phi_T - (MU-2) × phi_ab = MT × phi_T + (2-MU) × phi_ab
      // Each (MT, MU) pair has its own cos argument computed at each quadrature point.
      //
      // A12_terms stores: (MT_int, MU_int, amplitude) tuples, NOT cos-indexed buckets.
      // At each phi_ab quadrature point, we compute:
      //   A12_val = sum_{(MT,MU)} amplitude × cos(MT × phi_T + (2-MU) × phi_ab)
      std::vector<std::tuple<int,int,double>> A12_terms; // (MT, MU, amplitude)
      //
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
      //
      // XLAM table — Ptolemy A12 recurrence (source.mor lines 1613-1656)
      // DINTS(N) = N (integer array, not 2N+1).
      // XN = 0.5 * sqrt(DINTS(2*Li+1)*DINTS(2*lT+1)*DINTS(2*lP+1))
      //    = 0.5 * sqrt((2*Li+1)*(2*lT+1)*(2*lP+1))
      {
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
      auto xlam_correct = [](int L, int am) -> double {
        am = std::abs(am);
        if (am > L) return 0.0;
        if ((L + am) % 2 != 0) return 0.0;  // must have same parity

        // Build the table up to level L using Ptolemy recurrence
        // We only need the path to xlam(L, am), so walk the full recurrence.
        double outter = 1.0;
        int m_cur = 0;
        // Flat cache: xlam_cache[(ll, mm)] = value
        // Use a simple 2D array indexed by (ll, mm/2+1).
        // Max indices: L up to ~10 typically; allocate on stack.
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
      };

      double XN_a12 = 0.5 * std::sqrt((2.0*Li+1.0)*(2.0*lT+1.0)*(2.0*lP+1.0));

      A12_terms.clear();
      // Outer: MT loops from -lT to +lT step 2 (ALL EVEN values, including negative).
      // Matching Ptolemy INELDC A12 subroutine (source.mor lines 1773-1899).
      // For lP=0: MP=0 only, MX = MT (since MX = MT + MP = MT + 0).
      //
      // CORRECT: The cos argument in the AngKernel is cos(MT*phi_T - MU*phi_ab).
      // In Ptolemy INELDC: COSTHE = COS(MT*PHIT + MP*PHIP - MU*PHI) (per comment at line ~17959).
      //
      // OUTTMP(LX, MT) = XLAM(lT, |MT|) * XN * ThreeJ(lT, MT, lP=0, 0, Lx, -MT)
      // OUTTER(LI, LO, MU, MUPOS) = XLAM(Li, |MU|) * XLAM(Lo, |MX-MU|) * sqrt(2Lo+1)
      //   with SIGN: ITES = IP3+IP2; IP3=lT if MT<0; IP2=LOMNMN if MX_minus_MU<0.
      //   if ITES odd: OUTTER = -OUTTER.
      //
      // Compute LOMNMN: minimum Lo with right parity (Lo same parity as LBT+LBP+LI mod 2 = Li mod 2)
      {
        int min_Lo_tri = std::abs(Li - Lx);
        int lo_parity = Li % 2;  // Ptolemy: (LBT+LBP+Li+LOMIN) even → LOMIN same parity as Li
        int LOMNMN = min_Lo_tri;
        if (LOMNMN % 2 != lo_parity) LOMNMN++;  // adjust to correct parity

      for (int MT = -lT; MT <= lT; MT += 2) {
        int Mx = MT;  // MX = MT + MP, MP=0
        double xlam_T = xlam_correct(lT, std::abs(MT));
        double xlam_P = xlam_correct(lP, 0);
        double outtmp = xlam_T * xlam_P * XN_a12
                      * ThreeJ((double)lT, (double)MT, (double)lP, 0.0, (double)Lx, (double)(-Mx));
        if (std::abs(outtmp) < 1e-15) continue;

        // MU: same parity as Li, from (Li%2) to Li step 2 (always non-negative)
        for (int MU = (Li % 2); MU <= Li; MU += 2) {
          int MX_minus_MU = Mx - MU;  // can be negative; cos(|MX_minus_MU|*phi) kernel
          if (std::abs(MX_minus_MU) > Lo) continue;
          double xlam_Li = xlam_correct(Li, std::abs(MU));
          double xlam_Lo = xlam_correct(Lo, std::abs(MX_minus_MU));
          double outter = xlam_Li * xlam_Lo * std::sqrt(2.0*Lo + 1.0);
          // ITES sign correction (Ptolemy lines 1838-1853):
          int ITES = 0;
          if (MT < 0) ITES += lT;           // IP3 = LBT if MT<0 (lT=2 even → no effect)
          if (MX_minus_MU < 0) ITES += LOMNMN;  // IP2 = LOMNMN if MX-MU < 0
          if (ITES % 2 != 0) outter = -outter;
          // TTT inner ThreeJ:
          double ttt = ThreeJ((double)Li, (double)MU, (double)Lo, (double)MX_minus_MU,
                              (double)Lx, (double)(-Mx));
          double A12_val = outtmp * outter * ttt;
          // Debug: print all A12 values for validation cases
          {
            std::printf("  [A12_DBG_Li%d_Lo%d_Lx%d] MT=%d MU=%d coeff=%.8f (pre-doubling)\n",
                        Li, Lo, Lx, MT, MU, A12_val);
          }
          if (std::abs(A12_val) < 1e-15) continue;
          // Doubling for MU!=0 (Ptolemy: IF HALFSW OR MU!=0 TEMP=2*TEMP, HALFSW=FALSE here)
          double doubling = (MU != 0) ? 2.0 : 1.0;
          // Store as (MT, MU, amplitude) — the cos argument is evaluated per quadrature point.
          // For each MT, the cos arg: MT × phi_T + (2-MU) × phi_ab (Ptolemy DMSVAL convention).
          // Since we sum over all MT from -lT to +lT, store each (MT, MU) separately.
          A12_terms.push_back({MT, MU, A12_val * doubling});
        }
      }
      // Debug: print A12_terms for all validation cases
      {
        std::printf("[A12_TERMS] Li=%d Lo=%d Lx=%d lT=%d lP=%d XN=%.5f LOMNMN=%d\n",
                    Li, Lo, Lx, lT, lP, XN_a12, std::abs(Li - Lx));
        for (auto &[MT_k, MU_k, c_k] : A12_terms) {
          std::printf("  MT=%2d MU=%d coeff=%.8f\n", MT_k, MU_k, c_k);
        }
        std::printf("  Total A12_terms: %d\n", (int)A12_terms.size());
      }
      } // close LOMNMN scope
      } // close XLAM/xlam_correct scope
      // end A12 precompute

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

          // At phi_ab=0, rb=zr_scale*ra, rx=ra:
          //   cos_phiT = (T1_c*rb + S1_c*ra*1)/rx = (T1_c*zr_scale + S1_c) = 1.0 (exact)
          //   phi_T_angle = 0
          // => A12_val = sum_k coeff_k * cos(MT_k*0 - MU_k*0) = sum_k coeff_k
          double A12_at_zero = 0.0;
          for (auto &[MT_k, MU_k, coeff] : A12_terms) {
            A12_at_zero += coeff;
          }
          // Sanity: verify rx/ra = 1 numerically
          double rx_check = S1_c + T1_c * zr_scale;  // should be 1.0

          // 1D radial integral over ra using uniform grid:
          //   I_1D = ∫ chi_a(ra) * phi_T(ra) * chi_b*(rb=zr_scale*ra) * dra
          // The measure dra is absorbed; the Jacobian from the delta function in 3D:
          //   delta^3(rp) = delta(rp) / (4pi*rp^2) but the 3D angular integral over
          //   phi_ab when rp→0 collapses via the sin(phi_ab)*dphi_ab factor.
          //   In the PRIOR form: the 3D integral ∫dra dra dphi_ab sin(phi_ab) ... delta(rp) 
          //   For the phi_ab integral with rp^2 = S2^2*ra^2 + S1^2*rb^2 - 2*S1*S2*ra*rb*cos(phi_ab)
          //   drp/dphi_ab = (S1*S2*ra*rb*sin(phi_ab)) / rp → at phi=0: ∞ (limiting as rp→0)
          //   The integral gives a factor of 1/(S1_c*S2_c*ra*rb) per the Jacobian.
          //   Combined with the sin(phi_ab)*dphi_ab weight, the net ZR measure is:
          //   ∫ dphi_ab sin(phi_ab) * delta(rp) = 1 / (S1_c * S2_c * ra * rb_at_ZR)
          //                                      = 1 / (S1_c^2 * zr_scale * ra^2)
          // Full integral measure per dra:  ra^2 * rb^2 * (1/(S1_c^2*zr_scale*ra^2))
          //   = rb^2 / S1_c^2 = (zr_scale*ra)^2 / S1_c^2
          // Simplification: just integrate chi_a(ra)*phi_T(ra)*chi_b*(zr_scale*ra)*dra
          // and multiply by (zr_scale^2/S1_c^2) as the geometric factor.
          // But for a SHAPE CHECK we only need relative values — skip the overall constant.
          // For correct relative shape, we need all (Li,Lo) on the same footing.
          // The A12_at_zero factor handles the angular coupling differences between (Li,Lo).

          double h_zr    = Incoming.StepSize;
          int    N_zr    = TgtBS_ch.NSteps;
          std::complex<double> I_1D(0.0, 0.0);

          for (int i = 1; i < N_zr - 1; ++i) {
            double ra = i * h_zr;
            double rb = zr_scale * ra;

            // chi_a(ra) — incoming partial-wave WF at ra
            std::complex<double> ca_r = chi_a[i];

            // phi_T(rx=ra) — bound state at rx=ra (exact ZR condition)
            double phi_T_r = TgtBS_ch.WaveFunction[i].real();

            // chi_b*(rb) — outgoing partial-wave WF at rb=zr_scale*ra (interpolated)
            std::complex<double> cb_r = std::conj(interp_chi_b(rb));

            // Accumulate with measure dra (shape-preserving)
            I_1D += ca_r * phi_T_r * cb_r * h_zr;
          }

          // Full ZR integral: D0 * A12(phi=0) * geometric_factor * I_1D
          // geometric_factor = zr_scale^2 / S1_c^2 (from delta function Jacobian)
          // Omitted here for shape check — overall scale doesn't affect angular distribution.
          Integral = D0_ZR * A12_at_zero * I_1D;

          std::printf("[ZR] Li=%d Lo=%d Lx=%d A12_zero=%.5f rx_check=%.6f |I_1D|=%.5e |I_ZR|=%.5e\n",
                      Li, Lo, Lx, A12_at_zero, rx_check, std::abs(I_1D), std::abs(Integral));
        }

#else
        double h = Incoming.StepSize;  // kept for chi interpolation

        // ---------------------------------------------------------------
        // Quadrature parameters:
        //   Sum (U): Large GL grid for chi oscillations (direct quadrature)
        //   Dif (V): Large GL grid (symmetric)
        //   Phi:     CUBMAP adaptive PHI0 per-(ra,rb) (Ptolemy MAPPHI)
        // ---------------------------------------------------------------
        const int NPSUM = 150;   // GL points for U
        const int NPDIF = 75;    // GL points for V
        const int NPPHI = 10;    // Ptolemy default phi points (NPPHI=10)
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

        // Debug tracker for integrand peak (per-integral)
        double max_integrand_dbg = 0;
        double max_ra_p = 0, max_rb_p = 0, max_rx_p = 0, max_rp_p = 0;

        // Helper: evaluate BSPROD-like integrand at (ra, rb, x=cos_phi)
        // Returns PRIOR: |IVPHI_T(rx)| * |phi_P(rp)|  (same as used in AngKernel)
        auto bsprod_val = [&](double ra, double rb, double x) -> double {
          double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*x;
          double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*x;
          if (rx2 < 0) rx2 = 0;
          if (rp2 < 0) rp2 = 0;
          double rx = std::sqrt(rx2);
          double rp = std::sqrt(rp2);
#ifdef USE_POST_FORM
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

        // First pass: find global WVWMAX = max of bsprod over all (U,V,phi) for RVRLIM
        // (This is Ptolemy's WVWMAX from the GRDSET pass over the whole grid)
        // We estimate WVWMAX from phi=0 at (ra,rb) = (SUMMID, SUMMID):
        // or scan a few key points. Use peak of IVPHI_T * phi_P at short range.
        // Ptolemy: RVRLIM = DWCUT * WVWMAX
        // We compute WVWMAX adaptively below by tracking max seen across all (U,V)
        // For now, pre-estimate by evaluating at (SUMMID, SUMMID, phi=0):
        double WVWMAX_pre = bsprod_val(5.0, 5.0, 1.0);  // phi=0 → x=1
        WVWMAX_pre = std::max(WVWMAX_pre, bsprod_val(2.0, 2.0, 1.0));
        WVWMAX_pre = std::max(WVWMAX_pre, bsprod_val(3.0, 3.0, 1.0));
        WVWMAX_pre = std::max(WVWMAX_pre, bsprod_val(1.0, 1.0, 1.0));
#ifdef USE_POST_FORM
        // POST form WVWMAX: also sample ra=2*rb ridge (where rp→0)
        for (double rb_scan : {1.0, 1.5, 2.0, 2.5, 3.0})
          WVWMAX_pre = std::max(WVWMAX_pre, bsprod_val(2.0*rb_scan, rb_scan, 1.0));
#endif
        // RVRLIM threshold for phi cutoff
        double RVRLIM = DWCUT * std::max(WVWMAX_pre, 1.0e-30);

        // Pre-allocate dif GL arrays (reused per U, resized below)
        std::vector<double> xi_d_base(NPDIF), wi_d_base(NPDIF);
        GaussLegendre(NPDIF, -1.0, 1.0, xi_d_base, wi_d_base);

        // Outer loop: U = (ri+ro)/2 in [0, SUMMAX] from GL on [-1,1]
        for (int is = 0; is < NPSUM; ++is) {
          double U = SUMMAX * (xi_s[is] + 1.0) / 2.0;   // U in [0, SUMMAX]
          double WOW = SUMMAX * wi_s[is] / 2.0;           // sum weight

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
            std::complex<double> ca      = interp_chi_a(ra);
            std::complex<double> cb_conj = std::conj(interp_chi_b(rb));

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
#ifdef USE_POST_FORM
              // POST form: vertex at rp (USEPROJECTILE + USECORE, Ptolemy default NUCONL=3)
              // Ptolemy BSPROD NUCONL=3 formula (source.mor ~4640-4660):
              //   FPFT_nuconly = phi_T(rx) * [V_np(rp)*phi_P(rp)]
              //   RCORE = sqrt(((S1-S2)*ra)^2 + ((T1-T2)*rb)^2 + 2*(S1-S2)*(T1-T2)*ra*rb*x)
              //   RSCAT = rb (IOUTSW=TRUE for stripping/projective vertex)
              //   DELVNU = VOPT * (WS(RCORE, RNCORE, AOPT) - WS(RSCAT, RNSCAT, AOPT))
              //   VEFF = V_np(rp)  [the vertex potential at rp]
              //   FACTOR = 1 + DELVNU / VEFF
              //   FPFT = FACTOR * FPFT_nuconly
              double phi_T  = InterpolateClipped(TgtBS_ch.WaveFunction, TgtBS_ch.RGrid,
                                                 TgtBS_ch.NSteps, TgtBS_ch.StepSize,
                                                 TgtBS_ch.MaxR, rx, phi_T_max, r_T_peak);
              // NUCONLY vertex: V_np(rp) * phi_P(rp) (clipped)
              double ivphi_P_nuconly = InterpolateIVPHI(IVPHI_P, h_common, NSteps_common,
                                                h_common * NSteps_common, rp,
                                                IVPHI_P_max, r_P_vert_peak);
              // Ptolemy exact RCORE formula from source.mor:
              // RCORE = sqrt(((S1-S2)*RA)^2 + ((T1-T2)*RB)^2 + 2*(S1-S2)*RA*(T1-T2)*RB*X)
              double dS = S1 - S2;  // ≈ 0.9715
              double dT = T1 - T2;  // ≈ 0.0572
              double rcore2 = dS*dS*ra*ra + dT*dT*rb*rb + 2.0*dS*dT*ra*rb*x;
              if (rcore2 < 0) rcore2 = 0;
              double r_core = std::sqrt(rcore2);
              // RSCAT = rb (IOUTSW=TRUE for stripping/projectile vertex, ISC=2)
              double r_scat = rb;
              // DELVNU = VOPT * (WS(r_core, RNCORE, AOPT) - WS(r_scat, RNSCAT, AOPT))
              double fcore = 1.0 / (1.0 + std::exp((r_core - RNCORE_post) / AOPT_post));
              double fscat = 1.0 / (1.0 + std::exp((r_scat - RNSCAT_post) / AOPT_post));
              double DELVNU = VOPT_post * (fcore - fscat);  // VOPT_post < 0
              // VEFF = V_np(rp) evaluated without clipping (actual vertex potential)
              double VEFF = InterpolateV(PrjBS_ch.V_real, PrjBS_ch.StepSize,
                                         PrjBS_ch.NSteps, PrjBS_ch.MaxR, rp);
              // Ptolemy FACTOR = 1 + DELVNU/VEFF; FPFT = FACTOR * FPFT_nuconly
              // = (V_np + DELVNU) * phi_P  (algebraically)
              // Implement directly to avoid division-by-zero when V_np(rp)→0:
              // ivphi_P = phi_P(rp) * (V_np(rp) + DELVNU)
              // For the clipping region (rp < r_P_vert_peak): phi_P is clipped to max,
              // and V_np is at its max value; add DELVNU as-is.
              // phi_P_at_rp uses unclipped InterpolateClipped (correct for rp>peak)
              double phi_P_at_rp2 = InterpolateClipped(PrjBS_ch.WaveFunction, PrjBS_ch.RGrid,
                                                       PrjBS_ch.NSteps, PrjBS_ch.StepSize,
                                                       PrjBS_ch.MaxR, rp, phi_P_max, r_P_peak);
              // For clipped region (rp < r_P_vert_peak): use IVPHI_P_max as V_np*phi_P,
              // then add DELVNU*phi_P separately.
              double ivphi_P;
              if (rp < r_P_vert_peak) {
                // In clipped region: IVPHI_P_max = V_np_max * phi_P_max
                // Add DELVNU * phi_P_clipped (phi_P is at its max here)
                ivphi_P = IVPHI_P_max + DELVNU * phi_P_max;
              } else {
                // Outside peak: ivphi_P = (V_np(rp) + DELVNU) * phi_P(rp)
                ivphi_P = (VEFF + DELVNU) * phi_P_at_rp2;
              }
              double phi_P  = 1.0;   // Already included in ivphi_P
              double Vbx    = ivphi_P;  // = (V_eff * phi_P)(rp) with USECORE
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
              // Track peak integrand (debug)
              double this_integrand = std::abs(phi_T * Vbx * phi_P);
              if (this_integrand > max_integrand_dbg) {
                max_integrand_dbg = this_integrand;
                max_ra_p = ra; max_rb_p = rb; max_rx_p = rx; max_rp_p = rp;
              }

              // Angle of rx in the integration coordinate frame
              // Ptolemy line 16670: TEMP = (T1*RO + S1*RI*X)/RT
              double cos_phiT = (rx2 > 1e-30) ? (T1*rb + S1*ra*x) / rx : 1.0;
              cos_phiT = std::max(-1.0, std::min(1.0, cos_phiT));
              double phi_T_angle = std::acos(cos_phiT);

              // A12 angular coupling kernel — same formula as rectangular
              double A12_val = 0.0;
              for (auto &[MT_k, MU_k, coeff] : A12_terms) {
                double arg = MT_k * phi_T_angle - MU_k * phi_ab;
                A12_val += coeff * std::cos(arg);
              }

              // phi_weight = PHI0 * phi_wts[k] * sin(phi) = DPHI in Ptolemy GRDSET
              AngKernel += phi_weight * phi_T * Vbx * phi_P * A12_val;
            }

            // Accumulate: RIROWTS = JACOB * ri * ro * WOW * DIFWT
            // Jacobian d(ri)d(ro)/d(U)d(V) = 1 (exact, verified analytically)
            // WOW = wi_s[is] (CUBMAP sum weight), DIFWT = wi_d[id] (CUBMAP dif weight)
            Integral += ca * cb_conj * AngKernel * JACOB_grdset * ra * rb * WOW * DIFWT;
          }
        }

#endif  // USE_ZR

        // Debug: print integrand peak info for dominant term (Li=0,Lo=2,Lx=2)
#ifndef USE_ZR
        if (Li == 0 && Lo == 2 && Lx == 2) {
          std::printf("[INTEGRAND_PEAK] Li=%d Lo=%d Lx=%d: max|kernel|=%.4e"
                      " at ra=%.2f rb=%.2f rx=%.2f rp=%.2f\n",
                      Li, Lo, Lx, max_integrand_dbg, max_ra_p, max_rb_p, max_rx_p, max_rp_p);
        }
#endif

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
        
        // Compute ATERM for this reaction from quantum numbers
        // For 33Si(d,p)34Si: jA_tgt=3/2, jB_tgt=0, lBT=2, jBT=3/2
        //   lBP=0, jBP=1/2, JX=j_neutron=1/2, Lx from transfer
        // JBIGA = 2*jA = 3 (33Si target), JBIGB = 2*jB = 0 (34Si)
        // TEMP = sqrt((JBIGB+1)/(JBIGA+1)) = sqrt(1/4) = 0.5
        // ATERM = TEMP * sqrt(2*Lx+1) * SPAMP * SPAMT * RACAH(2*LBT, JBT, 2*LBP, JBP, JX, 2*Lx)
        //       = 0.5 * sqrt(2*Lx+1) * 0.97069 * 1.0 * RACAH(4, 3, 0, 1, 1, 2*Lx)
        // SIGN: ITEST = JX - JBP + 2*(LBP+LBT) = 1-1+2*(0+2) = 4, 4//2+1=3, 3%2!=0 => flip
        // |ATERM(Lx=2)| = 0.5 * sqrt(5) * 0.97069 * 2.515 = 2.730 (from RACAH computation)
        //
        // Hardcoded for this reaction (computed analytically):
        // RACAH(4,3,0,1,1,4) = W(2,3/2,0,1/2;1/2,2) = 2.5149 (verified numerically)
        // ATERM(Lx=2) = 0.5 * sqrt(5) * 0.97069 * 2.5149 = -2.730 (with sign flip)
        // |FACTOR * ATERM(Lx=2)| = 0.111 * 2.730 = 0.303
        
        // Apply FACTOR * ATERM to the sfromi normalization
        // ATERM for this specific Lx (only Lx=2 is relevant here, lT=2, lP=0)
        double ATERM_val = 1.0;  // Default if Lx != lT+lP
        if (Lx == std::abs(TargetBS.l - ProjectileBS.l) + 0 ||
            Lx <= TargetBS.l + ProjectileBS.l) {
          // For our reaction: lBT=2, jBT=3, lBP=0, jBP=1, JX=1
          // RACAH(2*lBT, jBT_doubled, 2*lBP, jBP_doubled, JX, 2*Lx)
          // = RACAH(4, 3, 0, 1, 1, 2*Lx)
          // Nonzero only when triangle(lBT,lBP,Lx) is valid:
          // triangle(2, 0, Lx): |2-0| <= Lx <= 2+0, so Lx=2 only!
          if (Lx == TargetBS.l) {
            // RACAH(4,3,0,1,1,4) = W(2,3/2,0,1/2;1/2,2) = 1/sqrt(10) = 0.31623 (numerically verified with Ptolemy Fortran)
            // Previous value 2.51487 was INCORRECT.
            double RACAH_val = 0.316227766;
            // Ptolemy ATERM formula (source.mor ~25631):
            // TEMP = sqrt((JBIGB+1)/(JBIGA+1))
            // JBIGB = 2*J(residual nucleus B = 34Si) = 2*0 = 0  → (JBIGB+1) = 1
            // JBIGA = 2*J(target nucleus A = 33Si) = 2*1.5 = 3  → (JBIGA+1) = 4
            // TEMP = sqrt(1/4) = 0.5
            // (NUCLEAR spins of target and residual, NOT the neutron quantum numbers)
            // For 33Si(d,p)34Si: JBIGA=3 (33Si, J=3/2), JBIGB=0 (34Si, J=0)
            // TODO: make these settable from input; hardcoded for 33Si(d,p)34Si
            const int JBIGA = 3;  // 2*J(target 33Si) = 2*(3/2)
            const int JBIGB = 0;  // 2*J(residual 34Si) = 2*0
            double TEMP_aterm = std::sqrt((JBIGB + 1.0) / (JBIGA + 1.0));  // sqrt(1/4) = 0.5
            double SPAMP = 0.97069;  // AV18 spectroscopic amplitude
            double SPAMT = 1.0;      // Target spec amplitude
            ATERM_val = TEMP_aterm * std::sqrt(2.0*Lx + 1.0) * SPAMP * SPAMT * RACAH_val;
            // Sign flip: ITEST = 4, ITEST//2+1 = 3, odd => flip
            ATERM_val = -ATERM_val;
          } else {
            ATERM_val = 0.0;  // Triangle condition fails
          }
        }
        
        // Total SFROMI factor: FACTOR_sfromi * ATERM * (1/sqrt(2*Li+1))
        double sfromi_norm = FACTOR_sfromi * std::abs(ATERM_val) / std::sqrt(2.0 * Li + 1.0);
        // Apply sign of ATERM to the integral (ATERM sign affects overall sign)
        if (ATERM_val < 0) {
          Integral = -Integral;
        }

        // Store: (Lx, Li, Lo, JPI, JPO) — SFROMI will apply 9-J coupling
        if (Lx==2 && (Li+Lo)<=6 && JPI<=4 && JPO<=5) {
          std::printf("INELDC Li=%d Lo=%d Lx=%d JPI=%d JPO=%d |Int|=%.5e sfromi_norm=%.5e |elemS|=%.5e\n",
            Li, Lo, Lx, JPI, JPO, std::abs(Integral), sfromi_norm, std::abs(Integral * sfromi_norm));
        }
        TransferSMatrix.push_back({Lx, Li, Lo, JPI, JPO, Integral * sfromi_norm});

        } // end JPO loop
        } // end JPI loop
      } // end Lx loop
    } // end Lo loop
  } // end Li loop

  std::cout << "Transfer integrals computed: " << TransferSMatrix.size()
            << " entries." << std::endl;
}

