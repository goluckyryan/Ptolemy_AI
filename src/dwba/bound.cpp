// bound.cpp — DWBA::CalculateKinematics and DWBA::CalculateBoundState
//
// Extracted from dwba.cpp (lines 181-248 and 2268-2487).
// Do NOT edit dwba.cpp directly — this file holds the canonical implementations.

#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// Module-level constants (duplicated from dwba.cpp for standalone compilation)
static const double HBARC_B  = 197.32697;  // MeV fm
static const double AMU_B    = 931.494;    // MeV/c^2
static const double FINE_B   = 1.0 / 137.035999;

// ---------------------------------------------------------------
// DWBA::CalculateKinematics
//
// Fills Incoming and Outgoing Channel kinematics:
//   Ecm, k (fm^-1), eta (Sommerfeld), mu (AMU), JSPS
//
// Uses relativistic kinematics for Ecm / momenta, then
// non-relativistic reduced mass for the potential.
// ---------------------------------------------------------------
void DWBA::CalculateKinematics() {
  // Incoming Channel (a + A)
  double ma = Incoming.Projectile.Mass;
  double mA = Incoming.Target.Mass;

  double E_tot_a = ma + Incoming.Elab;
  double p_a = std::sqrt(E_tot_a * E_tot_a - ma * ma);

  double s_in = (E_tot_a + mA) * (E_tot_a + mA) - p_a * p_a;
  double E_cm_tot = std::sqrt(s_in);
  Incoming.Ecm = E_cm_tot - ma - mA;

  double p_cm_in = p_a * mA / E_cm_tot;
  Incoming.k = p_cm_in / HBARC_B;

  // Coulomb parameter
  double Z1 = Incoming.Projectile.Z;
  double Z2 = Incoming.Target.Z;
  double E1_cm = std::sqrt(p_cm_in * p_cm_in + ma * ma);
  double E2_cm = std::sqrt(p_cm_in * p_cm_in + mA * mA);
  Incoming.eta = Z1 * Z2 * FINE_B * E1_cm * E2_cm / (p_cm_in * E_cm_tot);

  // Reduced mass for incoming channel (non-relativistic approx for potential)
  Incoming.mu = ma * mA / (ma + mA) / AMU_B;  // in AMU

  // Outgoing Channel (b + B)
  double mb = Outgoing.Projectile.Mass;
  double mB = Outgoing.Target.Mass;

  double mB_star = mB + Ex;
  double E_cm_tot_out = E_cm_tot;

  double s_out = E_cm_tot_out * E_cm_tot_out;
  double p_cm_out2 = (s_out - (mb + mB_star) * (mb + mB_star)) *
                     (s_out - (mb - mB_star) * (mb - mB_star)) / (4 * s_out);

  if (p_cm_out2 < 0) {
    std::cerr << "Error: Channel closed (below threshold)." << std::endl;
    p_cm_out2 = 0;
  }

  double p_cm_out = std::sqrt(p_cm_out2);
  Outgoing.k = p_cm_out / HBARC_B;

  double Z3 = Outgoing.Projectile.Z;
  double Z4 = Outgoing.Target.Z;
  double E3_cm = std::sqrt(p_cm_out2 + mb * mb);
  double E4_cm = std::sqrt(p_cm_out2 + mB_star * mB_star);

  Outgoing.eta = Z3 * Z4 * FINE_B * E3_cm * E4_cm / (p_cm_out * E_cm_tot_out);
  Outgoing.Ecm = E_cm_tot_out - mb - mB_star;

  // Reduced mass for outgoing channel
  Outgoing.mu = mb * mB / (mb + mB) / AMU_B;  // in AMU (use ground state mB)

  // Set projectile spin (JSPS = 2*spin) for each channel
  // For 33Si(d,p)34Si: incoming = deuteron (s=1, JSPS=2), outgoing = proton (s=1/2, JSPS=1)
  // Detect from Projectile mass: deuteron ~2 AMU, proton ~1 AMU
  Incoming.JSPS = (Incoming.Projectile.A >= 2) ? 2 : 1;  // deuteron=2, proton=1
  Outgoing.JSPS = (Outgoing.Projectile.A >= 2) ? 2 : 1;  // proton=1

  std::cout << "Kinematics:" << std::endl;
  std::cout << "  Incoming: k=" << Incoming.k << " eta=" << Incoming.eta
            << " Ecm=" << Incoming.Ecm << std::endl;
  std::cout << "  Outgoing: k=" << Outgoing.k << " eta=" << Outgoing.eta
            << " Ecm=" << Outgoing.Ecm << std::endl;
}

// ---------------------------------------------------------------
// DWBA::CalculateBoundState
//
// Ptolemy BOUND subroutine — exact translation (validated 2026-03-13)
//
// KEY conventions matching Ptolemy BOUND:
//   Array index I → r = I * h   (0-indexed: I=0 → r=0, I=1 → r=h)
//   Numerov: u[I+1] = (2 + h²*(κ²-AFAC*X+DL2/(h·I)²))·u[I] - u[I-1]
//   where X = V·V2[I] + V1[I]  (= negative total potential)
//   V2[I] = +WS(r) form factor (WOODSX type=1, V=-1)
//   V1[I] = -(Coulomb + ALS·f_SO)  (WOODSX type=2, V=-1 for SO)
//   ALS = VSO·[(0.25·(JP(JP+2)-JSP(JSP+2))-DL2)/JSP]  (σ·L coupling)
//   f_SO[I] = (2/(ASO·r))·exp/((1+exp)²)   (WOODSX type=2, V=-1)
//
// Normalization: phi = u/r stored in ch.WaveFunction (real part only)
//   Inner (0..NMATCH): scaled by SUMIN_sc = 1/sqrt(h·(SUMIN+XX²·SUMOUT))
//   Outer (NMATCH..NSTEPS): scaled by SUMOUT_sc = XX·SUMIN_sc
//   where XX = u_inner[NMATCH] / u_outer[NMATCH]
// ---------------------------------------------------------------
void DWBA::CalculateBoundState(Channel &ch, int n, int l, double j,
                               double bindingEnergy) {
  const double AMU_MEV = 931.494061;
  const double HBARC_L = 197.32697;
  double mu_MeV = ch.mu * AMU_MEV;
  double AFAC   = 2.0 * mu_MeV / (HBARC_L * HBARC_L);
  double kappa  = std::sqrt(2.0 * mu_MeV * bindingEnergy) / HBARC_L;
  double DL2    = l * (l + 1.0);

  int    NSTEPS = ch.NSteps;
  double h      = ch.StepSize;
  double h2     = h * h;

  std::cout << "Bound State (Ptolemy BOUND): E=" << -bindingEnergy
            << " MeV, kappa=" << kappa << " fm^-1"
            << ", L=" << l << ", J=" << j
            << ", mu=" << ch.mu << " AMU" << std::endl;

  // ---- Build V1 and V2 potential arrays (0-indexed, size NSTEPS+5) ----
  // Ptolemy BOUND: V2[I] = WS(r)/1  (WOODSX type=1, depth=-1 → returns +f_WS)
  //                V1[I] = -(Coulomb[I]) + ALS·f_SO[I]
  // For neutron: Coulomb=0 unless Z_proj>0

  // ALS: Ptolemy sigma·L factor (JP=2j, JSP=2s for projectile spin=1/2)
  int JP  = static_cast<int>(std::round(2.0 * j));
  int JSP = 1;  // 2*s, s=1/2 for neutron
  double ALS_raw = (0.25*(JP*(JP+2) - JSP*(JSP+2)) - DL2) / JSP;
  double ALS = ch.Pot.VSO * ALS_raw;

  // Coulomb constants (neutron: Z_proj=0 → all zero)
  double ZP = ch.Projectile.Z;
  double ZT = ch.Target.Z;
  double RC = ch.Pot.RC0 * std::pow(ch.Target.A, 1.0/3.0);
  double CI1 = ZP * ZT * HBARC_L * 1.5 / (137.036 * RC);
  double CI2 = ZP * ZT * HBARC_L * 0.5 / (137.036 * RC * RC * RC);
  double C0  = ZP * ZT * HBARC_L / 137.036;

  // Radii use ch.Target.A (the core nucleus mass)
  double R_nuc = ch.Pot.R0   * std::pow(ch.Target.A, 1.0/3.0);
  double R_SO  = ch.Pot.RSO0 * std::pow(ch.Target.A, 1.0/3.0);

  std::vector<double> V1(NSTEPS+5, 0.0);
  std::vector<double> V2(NSTEPS+5, 0.0);

  for (int I = 0; I <= NSTEPS+1; I++) {
    double r = I * h;

    // Coulomb (Ptolemy: V1 = -Coulomb)
    double VC = 0.0;
    if (ZP * ZT != 0.0) {
      VC = (r < RC) ? (CI1 - CI2*r*r) : (C0/r);
    }
    V1[I] = -VC;

    // Central WS: WOODSX type=1, V=-1 → V2[I] = 1/(1+exp((r-R)/A))
    double ex_c = (r > 1e-15) ? std::exp((r - R_nuc) / ch.Pot.A) : 1e100;
    V2[I] = 1.0 / (1.0 + ex_c);

    // SO: WOODSX type=2, V=-1 → f_SO = (2/(ASO·r))·exp/(1+exp)²
    if (ch.Pot.VSO != 0.0 && r > 1e-15) {
      double ex_s = std::exp((r - R_SO) / ch.Pot.ASO);
      double den  = 1.0 + ex_s;
      double f_so = (2.0 / (ch.Pot.ASO * r)) * ex_s / (den * den);
      V1[I] += ALS * f_so;   // Ptolemy: V1 += ALS * LWF1
    }
  }

  // ---- NMTOP and NMBOT (Ptolemy: NMTOP = max(R,A)/h + 0.5, clamp >=20) ----
  int NMTOP = static_cast<int>(std::max(R_nuc, ch.Pot.A) / h + 0.5);
  if (NMTOP < 20) NMTOP = 20;
  if (NMTOP > NSTEPS - 5) NMTOP = NSTEPS - 5;
  int NMBOT = std::max(NMTOP - 20, 3);

  // ---- V search via bisection (Ptolemy BOUND searches V to match derivatives) ----
  // We find V such that PHI = DROUT - DRIN = 0

  // Lambda: run one evaluation, return PHI and optionally build phi=u/r
  auto eval_bound = [&](double V0, bool build_phi, std::vector<double>* phi_out) -> double {

    // --- Outer (inward) integration ---
    std::vector<double> wf(NSTEPS+5, 0.0);
    wf[NSTEPS]   = std::exp(-kappa * NSTEPS * h);
    wf[NSTEPS-1] = std::exp(-kappa * (NSTEPS-1) * h);

    for (int II = NSTEPS-1; II >= NMBOT; II--) {
      double ALL = (II > 0) ? DL2/(h*II)/(h*II) : 0.0;
      double X   = V0*V2[II] + V1[II];
      double NF  = 2.0 + h2*(kappa*kappa - AFAC*X + ALL);
      wf[II-1]   = NF*wf[II] - wf[II+1];
    }

    // Find NMATCH
    int    NMATCH = NMTOP;
    double XXMAX  = 0.0;
    for (int II = NMTOP; II >= NMBOT; II--) {
      if (std::abs(wf[II]) > XXMAX) { XXMAX = std::abs(wf[II]); NMATCH = II; }
    }

    double SUMOUT = 0.5*(wf[NMATCH]*wf[NMATCH] + wf[NSTEPS-1]*wf[NSTEPS-1]);
    for (int II = NMATCH+1; II <= NSTEPS-2; II++) SUMOUT += wf[II]*wf[II];
    double VALOUT = wf[NMATCH];
    double DROUT  = 0.5/h * (wf[NMATCH+1] - wf[NMATCH-1]);

    std::vector<double> wf_outer;
    if (build_phi) wf_outer = wf;

    // --- Inner (forward) integration ---
    std::fill(wf.begin(), wf.end(), 0.0);
    wf[0] = 0.0;
    double TEMP = -AFAC*(V1[1] + V0*V2[1]) + kappa*kappa;
    wf[1] = h * (1.0 + h2*TEMP/6.0);

    double SUMIN = 0.0;
    int    nodes = 0;
    for (int II = 1; II <= NMATCH; II++) {
      double ALL = (II > 0) ? DL2/(h*II)/(h*II) : 0.0;
      double X   = V0*V2[II] + V1[II];
      double NF  = 2.0 + h2*(kappa*kappa - AFAC*X + ALL);
      wf[II+1]   = NF*wf[II] - wf[II-1];
      SUMIN += wf[II]*wf[II];
      if (wf[II]*wf[II-1] < 0.0) nodes++;
    }

    if (std::abs(VALOUT) < 1e-30) return 1e30;
    double XX2    = wf[NMATCH] / VALOUT;
    SUMIN        -= 0.5*wf[NMATCH]*wf[NMATCH];
    double norm_d = h*(SUMIN + XX2*XX2*SUMOUT);
    if (norm_d <= 0.0) return 1e30;
    double nf     = 1.0 / std::sqrt(norm_d);

    // PHI = DROUT_scaled - DRIN
    double DRIN   = nf   * 0.5/h * (wf[NMATCH+1]   - wf[NMATCH-1]);
    double DROUTs = XX2*nf * DROUT;
    double PHI    = DROUTs - DRIN;

    if (build_phi && phi_out) {
      phi_out->assign(NSTEPS+2, 0.0);
      double SUMIN_sc  = nf;
      double SUMOUT_sc = XX2 * nf;
      // phi = u/r; at r=0 for L=0 use Ptolemy's special formula
      if (l == 0) (*phi_out)[0] = SUMIN_sc;
      for (int II = 1; II <= NMATCH+1 && II < NSTEPS+1; II++) {
        (*phi_out)[II] = SUMIN_sc * wf[II] / (h*II);
      }
      for (int II = NMATCH+2; II <= NSTEPS; II++) {
        (*phi_out)[II] = SUMOUT_sc * wf_outer[II] / (h*II);
      }
    }
    return PHI;
  };

  // Scan for sign change then bisect
  double V_lo = -1.0, V_hi = -1.0, phi_lo = 0.0, phi_hi = 0.0, prev_phi = 1e30;
  for (double Vt = 10.0; Vt <= 200.0; Vt += 2.0) {
    double phi = eval_bound(Vt, false, nullptr);
    if (prev_phi != 1e30 && phi*prev_phi < 0.0 && V_lo < 0.0) {
      V_lo = Vt - 2.0; phi_lo = prev_phi;
      V_hi = Vt;       phi_hi = phi;
    }
    prev_phi = phi;
  }
  if (V_lo < 0.0) {
    std::cerr << "  ERROR: CalculateBoundState: no V solution found in [10,200] MeV\n";
    return;
  }

  // Bisect to convergence
  for (int it = 0; it < 100; it++) {
    double Vm = 0.5*(V_lo + V_hi);
    double pm = eval_bound(Vm, false, nullptr);
    if (phi_lo * pm < 0.0) { V_hi = Vm; phi_hi = pm; }
    else                   { V_lo = Vm; phi_lo = pm; }
    if (std::abs(pm) < 1e-10) break;
  }
  // suppress unused-variable warning for phi_hi
  (void)phi_hi;

  double V_sol = 0.5*(V_lo + V_hi);
  std::cout << "  Final V=" << V_sol << " MeV (n=" << n << " node(s))" << std::endl;

  // Build the normalized wavefunction phi = u/r
  std::vector<double> phi_vec;
  eval_bound(V_sol, true, &phi_vec);

  // Store in ch.WaveFunction (real part = phi = u/r)
  // ch.WaveFunction is indexed 0..NSTEPS-1; phi_vec is 0-indexed (I→r=I*h)
  ch.Pot.V = V_sol;
  for (int I = 0; I < NSTEPS && I < (int)phi_vec.size(); I++) {
    ch.WaveFunction[I] = std::complex<double>(phi_vec[I], 0.0);
  }
  double maxAmp = 0.0;
  for (int I = 0; I < NSTEPS; I++)
    maxAmp = std::max(maxAmp, std::abs(ch.WaveFunction[I].real()));
  std::cout << "  BS: maxAmp=" << maxAmp << "  kappa=" << kappa << " fm^-1" << std::endl;
}
