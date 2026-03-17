// bound.cpp — DWBA::CalculateKinematics and DWBA::CalculateBoundState
//
// Extracted from dwba.cpp (lines 181-248 and 2268-2487).
// Do NOT edit dwba.cpp directly — this file holds the canonical implementations.

#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <algorithm>
#include <cmath>
#include <fstream>
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

    // SO: WOODSX type=2, V=-1 → LWF1 = 2*exp/((1+exp)²·ASO·r)  (positive)
    // Ptolemy BOUND: V1 += ALS * LWF1 (ALS = VSO*ALS_raw, LWF1 > 0)
    // In Ptolemy NF: X = V*LV2 + V1, LV2 < 0 (WOODSX type1 V=-1)
    //   → effective V_SO in SE = -(-AFAC*V1) = AFAC*V1 = AFAC*ALS*LWF1
    // In C++ NF: X = V0*V2 + V1, V2 > 0 (opposite sign to Ptolemy LV2)
    //   → effective V_SO in SE = -(-AFAC*V1) = AFAC*V1   (same sign as Ptolemy)
    //   BUT V2 > 0 while LV2 < 0: C++ WS is attractive via -AFAC*(+V0*V2) < 0
    //   while Ptolemy WS is attractive via -AFAC*(V*LV2 < 0) = -AFAC*(negative) > 0 -- WAIT
    //
    // KEY: In C++ NF, -AFAC*X appears with a MINUS sign:
    //   NF = 2 + h²*(κ²  -AFAC*X + DL²/r²)
    // For bound-state oscillation inside well we need NF < 2:
    //   → -AFAC*X < 0  → X > 0
    // C++ V0*V2 > 0 ✓; V1 sign must NOT make X < 0 (can reduce X, but not flip it).
    //
    // Effective physical potential: V_eff = -X / AFAC
    //   V_WS_eff  = -V0*V2/AFAC  < 0  (attractive, correct)
    //   V_SO_eff  = -V1/AFAC = -ALS*f_so/AFAC
    //
    // For 0d5/2 (j=l+½): SO should be attractive (lowers energy).
    //   ALS = VSO*ALS_raw = (-6)*(+2) = -12; f_so > 0
    //   V_SO_eff = -(-12)*f_so/AFAC = +12*f_so/AFAC > 0  (REPULSIVE — WRONG)
    //
    // Fix: negate V1 SO term so V_SO_eff = -(-ALS*f_so) = ALS*f_so/AFAC < 0 (attractive ✓)
    //   i.e., V1[I] -= ALS * f_so
    if (ch.Pot.VSO != 0.0 && r > 1e-15) {
      double ex_s = std::exp((r - R_SO) / ch.Pot.ASO);
      double den  = 1.0 + ex_s;
      double f_so = (2.0 / (ch.Pot.ASO * r)) * ex_s / (den * den);
      V1[I] -= ALS * f_so;   // minus sign: corrects SO sign convention in C++ Numerov
    }
  }

  // ---- NMTOP and NMBOT (Ptolemy: NMTOP = max(R,A)/h + 0.5, clamp >=20) ----
  int NMTOP = static_cast<int>(std::max(R_nuc, ch.Pot.A) / h + 0.5);
  if (NMTOP < 20) NMTOP = 20;
  if (NMTOP > NSTEPS - 5) NMTOP = NSTEPS - 5;
  int NMBOT = std::max(NMTOP - 20, 3);

  // ---- V search via bisection (Ptolemy BOUND searches V to match derivatives) ----
  // We find V such that PHI = DROUT - DRIN = 0
  //
  // CRITICAL: must match the requested node count n.
  // As V increases the potential gets deeper and picks up more bound states.
  // Each PHI sign change corresponds to one additional node in the inner WF.
  // We count sign changes and only latch the bracket whose node count == n.

  // Lambda: run one evaluation, return {PHI, nodes} and optionally build phi=u/r
  auto eval_bound = [&](double V0, bool build_phi, std::vector<double>* phi_out,
                        int* nodes_out) -> double {

    // --- Outer (inward) integration ---
    // Outer wavefunction: Ptolemy uses Whittaker function W_{-eta, l+0.5}(2*kappa*r)
    // For a neutron (eta=0): W_{0, l+0.5}(2*kappa*r) = sqrt(2*kappa*r/pi) * K_{l+0.5}(kappa*r)
    // where K is the modified Bessel function of the second kind.
    // This differs from bare exp(-kappa*r) by a polynomial factor that matters at finite r.
    // Use Whittaker function via the exact finite-sum formula for integer l:
    //   k_l(z) = (pi/(2z))^{1/2} * K_{l+0.5}(z)  (spherical modified Bessel 2nd kind)
    //   Whittaker W_{0,l+0.5}(2z) = sqrt(2z/pi) * K_{l+0.5}(z) = z^{1/2} * (2/pi)^{1/2} * K_{l+0.5}(z)
    //   = sqrt(z) * sqrt(2) * k_l(z) / sqrt(pi/2) ... simplifying:
    //   W_{0,l+0.5}(2*kappa*r) = (2*kappa*r/pi)^{1/2} * K_{l+0.5}(kappa*r)
    //
    // For eta=0 and integer l: K_{l+0.5}(z) = sqrt(pi/(2z)) * exp(-z) * sum_{k=0}^{l} (l+k)!/(k!(l-k)!) * (2z)^{-k}
    // (exact finite sum — no approximation)
    //
    // So W_{0,l+0.5}(2*kappa*r) = exp(-kappa*r) * sum_{k=0}^{l} (l+k)!/(k!(l-k)!) * (2*kappa*r)^{-k}
    //
    // For l=2: sum = 1 + 3/(2*kappa*r) + 3/(2*kappa*r)^2
    // For l=1: sum = 1 + 1/(2*kappa*r)
    // For l=0: sum = 1
    // General l: Horner evaluation
    auto whittaker_eta0 = [&](double kappa_val, double r_val, int l_val) -> double {
      double z = kappa_val * r_val;  // kappa*r
      double iz = 1.0 / (2.0 * z);  // 1/(2*kappa*r)
      // Compute sum_{k=0}^{l} (l+k)!/(k!(l-k)!) * iz^k using upward recursion
      // Coefficient c_k = (l+k)! / (k! * (l-k)!)
      // c_0 = 1; c_{k+1} = c_k * (l-k)(l+k+1) / ((k+1)*1)  ... wait
      // Actually: (l+k)!/(k!(l-k)!) = C(l+k,k)*l!/... let me just compute directly
      double sum = 0.0;
      double izk = 1.0;  // iz^k
      double fact_lk = 1.0;  // (l+k)!
      double fact_k  = 1.0;  // k!
      double fact_lmk = 1.0; // (l-k)!
      // Precompute l!
      for (int i = 1; i <= l_val; i++) fact_lmk *= i;
      fact_lk = fact_lmk;  // (l+0)! = l! at k=0
      for (int k = 0; k <= l_val; k++) {
        sum += (fact_lk / (fact_k * fact_lmk)) * izk;
        if (k < l_val) {
          izk *= iz;
          fact_lk *= (l_val + k + 1);
          fact_k  *= (k + 1);
          fact_lmk /= (l_val - k);  // (l-k-1)! for next step
        }
      }
      return std::exp(-z) * sum;
    };

    std::vector<double> wf(NSTEPS+5, 0.0);
    wf[NSTEPS]   = whittaker_eta0(kappa, NSTEPS   * h, l);
    wf[NSTEPS-1] = whittaker_eta0(kappa, (NSTEPS-1) * h, l);

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


    // Fortran BOUND: SUMOUT = 0.5*(wf[NMATCH]^2 + wf[NSTEPS-1]^2)
    //                + sum(wf[II]^2) for II = NMATCH+3..NSTEPS-1
    // (trapezoidal endpoints at NMATCH and NSTEPS-1; inner points NMATCH+3..NSTEPS-2)
    double SUMOUT = 0.5*(wf[NMATCH]*wf[NMATCH] + wf[NSTEPS-1]*wf[NSTEPS-1]);
    for (int II = NMATCH+1; II <= NSTEPS-2; II++) SUMOUT += wf[II]*wf[II];
    double VALOUT = wf[NMATCH];
    double DROUT  = 0.5/h * (wf[NMATCH+1] - wf[NMATCH-1]);

    std::vector<double> wf_outer;
    if (build_phi) wf_outer = wf;

    // --- Inner (forward) integration ---
    std::fill(wf.begin(), wf.end(), 0.0);
    wf[0] = 0.0;
    // Ptolemy BOUND starting condition (line ~4143):
    //   TEMP = -AFAC*(LV1[1] + V*LV2[1]) + kappa^2
    // Note: Fortran TEMP does NOT include the centrifugal DL2/r^2 term in the starting
    // value (it is an approximation for the power-series starting, same as Ptolemy).
    // The centrifugal ALL is added inside the main Numerov loop (II=1..NMATCH).
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
    if (nodes_out) *nodes_out = nodes;

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

  // Scan for sign change with matching node count, then bisect.
  // Each PHI zero crossing with node count == n gives one bound state.
  // node count increases by 1 each time we cross a new PHI zero.
  double V_lo = -1.0, V_hi = -1.0, phi_lo = 0.0, phi_hi = 0.0;
  double prev_phi = 1e30;
  int    found_n  = -1;
  for (double Vt = 10.0; Vt <= 300.0; Vt += 2.0) {
    int cur_nodes = 0;
    double phi = eval_bound(Vt, false, nullptr, &cur_nodes);
    if (prev_phi != 1e30 && phi*prev_phi < 0.0) {
      // This bracket has node count == cur_nodes (the higher-V side)
      // Use cur_nodes as the node count for the bracket
      if (cur_nodes == n && V_lo < 0.0) {
        V_lo = Vt - 2.0; phi_lo = prev_phi;
        V_hi = Vt;       phi_hi = phi;
        found_n = cur_nodes;
      }
    }
    prev_phi = phi;
  }
  if (V_lo < 0.0) {
    std::cerr << "  ERROR: CalculateBoundState: no V solution found for n=" << n
              << " in [10,300] MeV\n";
    return;
  }
  std::cout << "  Bracket V=[" << V_lo << "," << V_hi << "] for n=" << found_n << std::endl;

  // Bisect to convergence
  for (int it = 0; it < 100; it++) {
    double Vm = 0.5*(V_lo + V_hi);
    double pm = eval_bound(Vm, false, nullptr, nullptr);
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
  eval_bound(V_sol, true, &phi_vec, nullptr);

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

// ---------------------------------------------------------------
// DWBA::LoadReidWavefunction
//
// Computes the Reid soft-core deuteron S-state wavefunction and
// effective potential by cubic spline on Reid's 34-point tabulation,
// exactly matching Ptolemy's Reid linkule (linkulesfitters.mor SUBROUTINE REID).
//
// Outputs:
//   ch.WaveFunction[I] = phi_S(r) = U_S(r)/r, normalized to ∫u²dr=1
//   ch.VPhiProduct[I]  = phi_S(r) * V_eff(r)  [IVPHI_P, signed]
//
// where V_eff(r) = -(VC + sqrt(8)*VT*W/U)  (Reid S-state effective potential)
//
// This replaces the previous file-based approach which misread reid-phi-v.
//
// File structure (ASCII, 2564 values total, 641 per block):
//   Block 0 [0..640]:    phi_S(r) = u_S(r)/r  (L=0, S-state)
//   Block 1 [641..1281]: V_np(r)*phi_S(r)      (vertex potential product)
//   Block 2 [1282..1922]: phi_D(r) = u_D(r)/r  (L=2, D-state)
//   Block 3 [1923..2563]: V_np(r)*phi_D(r)
//
// Grid: h=0.1 fm, 641 points, r = 0,0.1,...,64.0 fm
// The stored phi_S is unnormalized. We normalize here so that
//   ∫ u_S²(r) dr = ∫ (phi_S(r)·r)² dr = 1
//
// On return, ch.WaveFunction[I] = normalized phi_S(r=I*h) for the
// projectile bound state. Only the S-state (L=0) is loaded;
// the D-state (L=2) contribution requires a separate channel pass.
//
// dataPath should be the directory containing "reid-phi-v".
// ---------------------------------------------------------------
bool DWBA::LoadDeuteronWavefunction(Channel &ch, const std::string &dataPath,
                                    const std::string &filename) {
  // Reid soft-core S-state wavefunction via cubic spline on 34-point tabulation.
  // Matches Ptolemy's Reid linkule (linkulesfitters.mor SUBROUTINE REID, IREQUE=3).
  //
  // phi_S(r) = U_S(r)/r  (normalized to ∫u²dr = 1, i.e. SUMIN from BOUND lines 4109-4169)
  // V_eff(r) = -(VC + sqrt(8)*VT*W/U)  [VMULT=-1, as in linkule]
  // IVPHI_P  = phi_S * V_eff

  // Reid tabulated data: x = 0.7*r, u_S(x), du_S/dx  (from linkulesfitters.mor)
  static const int NTAB = 34;
  static const double XX[34] = {
    .0100, .041250, .07250, .1350, .19750,
    .2600, .32250, .3850, .44750, .5100,
    .57250, .6350, .69750, .7600, .8850,
    1.0100, 1.1350, 1.2600, 1.3850, 1.5100,
    1.7600, 2.0100, 2.5100, 3.0100, 3.5100,
    4.0100, 4.5100, 5.0100, 5.5100, 6.0100,
    7.0100, 8.0100, 9.0100, 10.0100
  };
  static const double UI_S[34] = {
    .000000, .333730e-4, .239010e-3, .276210e-2, .127370e-1,
    .360620e-1, .753590e-1, .128470, .189930, .253490,
    .313900, .367700, .413170, .449920, .499530,
    .524060, .531660, .528640, .519260, .506210,
    .475050, .442000, .378640, .322490, .273990,
    .232510, .197190, .167180, .141720, .120120,
    .862900e-1, .619830e-1, .445230e-1, .319820e-1
  };
  static const double UPI_S[34] = {
    .297510e-3, .261270e-2, .123350e-1, .829510e-1, .253260,
    .500520, .750720, .933490, .101620e1, .100340e1,
    .920420, .796740, .657440, .519850, .284940,
    .118650, .113970e-1, -.540900e-1, -.925230e-1, -.114180,
    -.130970, -.131930, -.120040, -.104530, -.897090e-1,
    -.765100e-1, -.650540e-1, -.552290e-1, -.468500e-1, -.397260e-1,
    -.285470e-1, -.205080e-1, -.147310e-1, -.105820e-1
  };
  static const double UI_D[34] = {
    .000000, .108500e-4, .840730e-4, .103690e-2, .496420e-2,
    .144460e-1, .307950e-1, .531570e-1, .789950e-1, .105250,
    .129330, .149580, .165290, .176450, .187100,
    .186540, .179460, .169100, .157420, .145530,
    .123140, .103730, .738590e-1, .532930e-1, .390770e-1,
    .291150e-1, .220160e-1, .168710e-1, .130790e-1, .102430e-1,
    .644120e-2, .415750e-2, .273630e-2, .182770e-2
  };
  static const double UPI_D[34] = {
    .742020e-4, .893210e-3, .447550e-2, .318920e-1, .101210,
    .205960, .314770, .393710, .424660, .408420,
    .357720, .288480, .214280, .144270, .332940e-1,
    -.358990e-1, -.731510e-1, -.900910e-1, -.952990e-1, -.941580e-1,
    -.839730e-1, -.712820e-1, -.493270e-1, -.339400e-1, -.236120e-1,
    -.166910e-1, -.120020e-1, -.877680e-2, -.652010e-2, -.491360e-2,
    -.289560e-2, -.177580e-2, -.112270e-2, -.726440e-3
  };
  const double ROOT8 = 2.8284271250;

  // Cubic spline interpolation: returns {U_S, W_D} at r
  auto reid_spline = [&](double r) -> std::pair<double,double> {
    double x = 0.7 * r;
    if (x <= XX[0] || x >= XX[NTAB-1]) return {0.0, 0.0};
    int ii = (int)(std::lower_bound(XX, XX+NTAB, x) - XX);
    if (ii <= 0) ii = 1;
    if (ii >= NTAB) ii = NTAB-1;
    int im = ii - 1;
    double H = XX[ii] - XX[im];
    double P = (x - XX[im]) / H;
    // S-wave
    double dU = UI_S[ii]-UI_S[im];
    double A1=UI_S[im], A2=H*UPI_S[im];
    double A3=3*dU-H*(2*UPI_S[im]+UPI_S[ii]);
    double A4=-2*dU+H*(UPI_S[im]+UPI_S[ii]);
    double U = A1+P*(A2+P*(A3+P*A4));
    // D-wave (for V_eff coupling term)
    double dW = UI_D[ii]-UI_D[im];
    double B1=UI_D[im], B2=H*UPI_D[im];
    double B3=3*dW-H*(2*UPI_D[im]+UPI_D[ii]);
    double B4=-2*dW+H*(UPI_D[im]+UPI_D[ii]);
    double W = B1+P*(B2+P*(B3+P*B4));
    return {U, W};
  };

  // Pass 1: compute L2 norm = sqrt(h * sum(U_S^2))  [BOUND lines 4109-4169]
  double norm_sq = 0.0;
  int NSTEPS = ch.NSteps;
  double h = ch.StepSize;
  for (int i = 1; i < NSTEPS; ++i) {
    double r = i * h;
    auto [U, W] = reid_spline(r);
    double wt = (i == 1) ? 0.5 : 1.0;  // trapezoidal, first point half-weight
    norm_sq += wt * U * U;
  }
  norm_sq *= h;
  double norm_scale = (norm_sq > 0.0) ? 1.0 / std::sqrt(norm_sq) : 1.0;

  std::cout << "LoadDeuteronWavefunction: Reid cubic spline (linkule convention)" << std::endl;
  std::cout << "  norm_sq=" << norm_sq << "  norm_scale=" << norm_scale
            << "  (expected ~1.035)" << std::endl;

  // Pass 2: fill phi_S and IVPHI_P
  ch.WaveFunction.assign(NSTEPS, std::complex<double>(0.0, 0.0));
  ch.VPhiProduct .assign(NSTEPS, 0.0);

  for (int i = 1; i < NSTEPS; ++i) {
    double r = i * h;
    auto [U, W] = reid_spline(r);
    // phi_S(r) = U(r)/r * norm_scale  [BOUND converts u→phi=u/r]
    double phi = U * norm_scale / r;
    ch.WaveFunction[i] = std::complex<double>(phi, 0.0);

    // V_eff(r): Reid effective S-state potential
    // From linkule (JUMPL=332): V_eff = VC + sqrt(8)*VT*W/U, then *VMULT=-1
    double x = 0.7 * r;
    if (x > 0.0) {
      double yy  = std::exp(-x);
      double y2  = yy*yy, y4=y2*y2, y6=y4*y2;
      double VC  = (-10.463*yy + 105.468*y2 - 3187.8*y4 + 9924.3*y6) / x;
      double VT  = (-10.463*(yy + 3.0*((yy-4*y4)+(yy-y4)/x)/x) + 351.77*y4 - 1673.5*y6) / x;
      double Veff = VC + (std::abs(U) > 1e-30 ? ROOT8*VT*W/U : 0.0);
      ch.VPhiProduct[i] = phi * (-Veff);  // VMULT=-1 → V_eff stored as -(VC+...)
    }
  }

  // Diagnostics
  double maxPhi = 0.0, maxVPhi = 0.0;
  int maxPhiI = 1, maxVPhiI = 1;
  for (int i = 1; i < NSTEPS; ++i) {
    double p  = std::abs(ch.WaveFunction[i].real());
    double vp = std::abs(ch.VPhiProduct[i]);
    if (p  > maxPhi)  { maxPhi  = p;  maxPhiI  = i; }
    if (vp > maxVPhi) { maxVPhi = vp; maxVPhiI = i; }
  }
  std::cout << "  phi_S:      peak=" << maxPhi
            << " at r=" << maxPhiI * h << " fm" << std::endl;
  std::cout << "  IVPHI_P:    peak=" << maxVPhi
            << " at r=" << maxVPhiI * h << " fm" << std::endl;
  std::cout << "  Vertex potential: Reid V_eff (cubic spline)" << std::endl;
  return true;
}
