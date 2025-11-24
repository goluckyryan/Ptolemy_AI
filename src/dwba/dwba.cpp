#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <cmath>
#include <iomanip>
#include <iostream>

// Constants
const double HBARC = 197.32697; // MeV fm
const double AMU = 931.494;     // MeV/c^2
const double FINE_STRUCTURE = 1.0 / 137.035999;

DWBA::DWBA() {
  AngleMin = 0.0;
  AngleMax = 180.0;
  AngleStep = 1.0;
  Ex = 0.0;

  // Initialize Bound States with safe defaults
  TargetBS.n = 0;
  TargetBS.l = 0;
  TargetBS.j = 0.0;
  TargetBS.BindingEnergy = 0.0;
  ProjectileBS.n = 0;
  ProjectileBS.l = 0;
  ProjectileBS.j = 0.0;
  ProjectileBS.BindingEnergy = 0.0;
}

DWBA::~DWBA() {}

void DWBA::SetReaction(const std::string &target, const std::string &projectile,
                       const std::string &ejectile, const std::string &recoil) {
  // Incoming.Target = Isotope(target);
  // Incoming.Projectile = Isotope(projectile);
  // Outgoing.Ejectile = Isotope(ejectile);
  // Outgoing.Recoil = Isotope(recoil);

  // Assuming A(a,b)B reaction:
  // Incoming: Target A, Projectile a
  // Outgoing: Ejectile B, Recoil b (light product)
  // Wait, usually recoil is the light product? No, recoil is usually the heavy
  // residual nucleus. Let's check standard notation A(a,b)B Target A,
  // Projectile a -> Ejectile b (light), Residual B (heavy) In InputGenerator:
  // target(projectile, recoil)ejectile
  // target=206Hg, proj=d, recoil=p, ejectile=207Hg
  // So recoil is the light product (p), ejectile is the heavy product (207Hg).
  // My naming in InputGenerator was:
  // targetName, projectileName, recoilName (light), ejectileName (heavy)

  // Let's align with that.
  // Incoming.Target = Isotope(target);
  // Incoming.Projectile = Isotope(projectile);
  // Outgoing.Target = Isotope(ejectile); // Residual nucleus acts as "target"
  // for outgoing? Outgoing.Projectile = Isotope(recoil); // Light product acts
  // as "projectile" for outgoing?

  // Actually, for DWBA, we have Incoming Channel (a+A) and Outgoing Channel
  // (b+B). Incoming: Projectile a, Target A Outgoing: Projectile b, Target B
  // (Residual)

  Incoming.Projectile = Isotope(projectile);
  Incoming.Target = Isotope(target);

  Outgoing.Projectile = Isotope(recoil); // Light outgoing particle
  Outgoing.Target = Isotope(ejectile);   // Heavy residual nucleus
}

void DWBA::SetEnergy(double Elab) { Incoming.Elab = Elab; }

void DWBA::SetExcitation(double Ex) { this->Ex = Ex; }

void DWBA::SetAngles(double min, double max, double step) {
  AngleMin = min;
  AngleMax = max;
  AngleStep = step;
}

void DWBA::SetIncomingPotential(const ChannelPotential &pot) {
  Incoming.Pot = pot;
}

void DWBA::SetOutgoingPotential(const ChannelPotential &pot) {
  Outgoing.Pot = pot;
}

void DWBA::SetTargetBoundState(int n, int l, double j, double bindingEnergy,
                               const ChannelPotential &pot) {
  TargetBS.n = n;
  TargetBS.l = l;
  TargetBS.j = j;
  TargetBS.BindingEnergy = bindingEnergy;
  TargetBS.Pot = pot;
}

void DWBA::SetProjectileBoundState(int n, int l, double j, double bindingEnergy,
                                   const ChannelPotential &pot) {
  ProjectileBS.n = n;
  ProjectileBS.l = l;
  ProjectileBS.j = j;
  ProjectileBS.BindingEnergy = bindingEnergy;
  ProjectileBS.Pot = pot;
}

void DWBA::CalculateKinematics() {
  // Incoming Channel (a + A)
  double ma = Incoming.Projectile.Mass;
  double mA = Incoming.Target.Mass;

  // Relativistic kinematics? Ptolemy uses relativistic kinematics.
  // E_cm^2 = (E_a + E_A)^2 - (p_a + p_A)^2
  // Lab frame: A at rest. E_A = mA, p_A = 0.
  // E_a = ma + Elab (Kinetic energy) -> Wait, Elab usually means kinetic energy
  // of projectile. Total energy of projectile E_tot_a = ma + Incoming.Elab
  // p_a^2 = E_tot_a^2 - ma^2

  double E_tot_a = ma + Incoming.Elab;
  double p_a = std::sqrt(E_tot_a * E_tot_a - ma * ma);

  double s_in = (E_tot_a + mA) * (E_tot_a + mA) - p_a * p_a;
  double E_cm_tot = std::sqrt(s_in);
  Incoming.Ecm =
      E_cm_tot - ma - mA; // Kinetic energy in CM? Or total CM energy?
  // Usually Ecm refers to kinetic energy available.

  // Reduced mass for wave number calculation
  // Relativistic reduced mass?
  // k^2 = 2*mu*E / hbar^2 (Non-relativistic)
  // Relativistic k: p_cm / hbar

  // Momentum in CM:
  // p_cm = p_lab * mA / E_cm_tot
  double p_cm_in = p_a * mA / E_cm_tot;
  Incoming.k = p_cm_in / HBARC;

  // Coulomb parameter eta = Z1*Z2*e^2 / (hbar*v)
  // eta = Z1*Z2 * (e^2/hbar*c) * (c/v)
  // v = p/E (relativistic) -> v_rel in CM?
  // Standard formula: eta = Z1*Z2 * alpha * (E_tot / p_cm) ?
  // Let's check Ptolemy manual or source.
  // In WAVSET/WAVPOT, eta is used.

  // From source.f:
  // ETA = ZP*ZT*AFINE * (E/P) ? No, let's check.
  // In KINEMA (if it exists) or similar.

  // Standard relativistic eta:
  // eta = Z1 Z2 alpha * (E1 E2 + p^2) / (p * E_tot) is for something else?

  // Let's use the definition: eta = alpha * Z1 * Z2 / beta_rel
  // v_rel = |v1 - v2|
  // In CM, v1 = p/E1, v2 = -p/E2.
  // v_rel = p/E1 + p/E2 = p * (E1+E2)/(E1*E2).
  // So eta = Z1 * Z2 * alpha * E1 * E2 / (p * (E1+E2)) = Z1 * Z2 * alpha * E1 *
  // E2 / (p * E_tot).

  double Z1 = Incoming.Projectile.Z;
  double Z2 = Incoming.Target.Z;
  double E1_cm = std::sqrt(p_cm_in * p_cm_in + ma * ma);
  double E2_cm = std::sqrt(p_cm_in * p_cm_in + mA * mA);
  Incoming.eta =
      Z1 * Z2 * FINE_STRUCTURE * E1_cm * E2_cm / (p_cm_in * E_cm_tot);

  // Outgoing Channel (b + B)
  double mb = Outgoing.Projectile.Mass;
  double mB = Outgoing.Target.Mass;

  // Q-value
  double Qval = (ma + mA) - (mb + mB); // Mass difference
  // Ex is excitation of residual nucleus B?
  // Effective Q = Qval - Ex.

  double E_cm_out = E_cm_tot - Ex; // Total energy in CM for outgoing?
  // Wait, E_cm_tot is conserved.
  // But masses changed.
  // E_cm_tot_out = E_cm_tot.
  // E_cm_tot_out = sqrt(p_out^2 + mb^2) + sqrt(p_out^2 + mB^2) + Ex?
  // No, Ex is internal energy of B. So mB_star = mB + Ex.

  double mB_star = mB + Ex;
  double E_cm_tot_out = E_cm_tot; // Energy conservation

  // Need to find p_out such that sqrt(p^2 + mb^2) + sqrt(p^2 + mB_star^2) =
  // E_cm_tot_out s_out = E_cm_tot_out^2 p_out^2 = (s_out - (mb + mB_star)^2) *
  // (s_out - (mb - mB_star)^2) / (4 * s_out)

  double s_out = E_cm_tot_out * E_cm_tot_out;
  double p_cm_out2 = (s_out - (mb + mB_star) * (mb + mB_star)) *
                     (s_out - (mb - mB_star) * (mb - mB_star)) / (4 * s_out);

  if (p_cm_out2 < 0) {
    std::cerr << "Error: Channel closed (below threshold)." << std::endl;
    p_cm_out2 = 0;
  }

  double p_cm_out = std::sqrt(p_cm_out2);
  Outgoing.k = p_cm_out / HBARC;

  double Z3 = Outgoing.Projectile.Z;
  double Z4 = Outgoing.Target.Z;
  double E3_cm = std::sqrt(p_cm_out2 + mb * mb);
  double E4_cm = std::sqrt(p_cm_out2 + mB_star * mB_star);

  Outgoing.eta =
      Z3 * Z4 * FINE_STRUCTURE * E3_cm * E4_cm / (p_cm_out * E_cm_tot_out);

  // Store outgoing Ecm (kinetic energy in CM)
  Outgoing.Ecm = E_cm_tot_out - mb - mB_star;

  std::cout << "Kinematics:" << std::endl;
  std::cout << "Incoming: k=" << Incoming.k << " eta=" << Incoming.eta
            << " Ecm=" << Incoming.Ecm << std::endl;
  std::cout << "Outgoing: k=" << Outgoing.k << " eta=" << Outgoing.eta
            << " Ecm=" << Outgoing.Ecm << std::endl;
}

void DWBA::Calculate() {
  CalculateKinematics();
  PrintParameters();

  // Setup grids and potentials
  WavSet(Incoming);
  WavSet(Outgoing);

  // Setup Integration Grid
  GrdSet();

  // Perform Finite Range Integration
  InelDc();

  // Calculate Cross Sections
  XSectn();
}

void DWBA::PrintParameters() {
  std::cout
      << "================================================================="
      << std::endl;
  // Constants
  const double AMU_MEV = 931.494061;
  const double HBARC = 197.32697;

  // Projectile BS
  // ProjectileBS
  // .mu; // Need to ensure this is set. It's set in InelDc, but not here?
  // Wait, ProjectileBS is a BoundState struct, it doesn't have mu.
  // ProjectileBSChannel has mu. But that's local to InelDc.
  // We need to calculate mu here.
  // Projectile BS: Particle x + Core b.
  // x = ma - mb.
  double ma = Incoming.Projectile.Mass; // In MeV
  double mb = Outgoing.Projectile.Mass; // In MeV
  double mx = ma - mb;
  double mu_pbs = mb * mx / (mb + mx); // In MeV
  // Kappa = sqrt(2 * mu * E) / hbarc. mu in MeV, E in MeV.
  double kappa_pbs =
      std::sqrt(2.0 * mu_pbs * (ProjectileBS.BindingEnergy)) / HBARC;

  std::cout << "        PROJECTILE BOUND STATE PARAMETERS" << std::endl;
  std::cout << "E = " << std::fixed << std::setprecision(4)
            << -ProjectileBS.BindingEnergy << " MEV"
            << "     KAPPA = " << std::setprecision(5) << kappa_pbs
            << std::endl;
  std::cout << " PROJECTILE MASS = " << std::fixed << std::setprecision(2)
            << mx / AMU_MEV << " AMU"
            << "     TARGET MASS = " << mb / AMU_MEV << " AMU"
            << "     REDUCED MASS = " << std::setprecision(2) << mu_pbs
            << " MEV/C**2" << std::endl;
  std::cout << " L = " << ProjectileBS.l << "  N = " << ProjectileBS.n
            << std::endl;
  std::cout << " POTENTIAL: V=" << std::setprecision(4) << ProjectileBS.Pot.V
            << " R0=" << ProjectileBS.Pot.R0 << " A=" << ProjectileBS.Pot.A
            << std::endl;
  std::cout << "            VSO=" << ProjectileBS.Pot.VSO
            << " RSO0=" << ProjectileBS.Pot.RSO0
            << " ASO=" << ProjectileBS.Pot.ASO << std::endl;
  std::cout << "            RC0=" << ProjectileBS.Pot.RC0 << std::endl;

  std::cout << std::endl;

  // Target BS
  // Target BS: Particle x + Core A.
  double mA = Incoming.Target.Mass;    // In MeV
  double mu_tbs = mA * mx / (mA + mx); // In MeV
  double kappa_tbs = std::sqrt(2.0 * mu_tbs * (TargetBS.BindingEnergy)) / HBARC;

  std::cout << "        TARGET BOUND STATE PARAMETERS" << std::endl;
  std::cout << "E = " << std::fixed << std::setprecision(4)
            << -TargetBS.BindingEnergy << " MEV"
            << "     KAPPA = " << std::setprecision(5) << kappa_tbs
            << std::endl;
  std::cout << " PROJECTILE MASS = " << std::fixed << std::setprecision(2)
            << mx / AMU_MEV << " AMU"
            << "     TARGET MASS = " << mA / AMU_MEV << " AMU"
            << "     REDUCED MASS = " << std::setprecision(2) << mu_tbs
            << " MEV/C**2" << std::endl;
  std::cout << " L = " << TargetBS.l << "  N = " << TargetBS.n << std::endl;
  std::cout << " POTENTIAL: V=" << std::setprecision(4) << TargetBS.Pot.V
            << " R0=" << TargetBS.Pot.R0 << " A=" << TargetBS.Pot.A
            << std::endl;
  std::cout << "            VSO=" << TargetBS.Pot.VSO
            << " RSO0=" << TargetBS.Pot.RSO0 << " ASO=" << TargetBS.Pot.ASO
            << std::endl;
  std::cout << "            RC0=" << TargetBS.Pot.RC0 << std::endl;

  std::cout << std::endl;
  std::cout << "        OPTICAL MODEL SCATTERING FOR THE INCOMING CHANNEL"
            << std::endl;
  std::cout << "E LAB = " << Incoming.Elab << " MEV" << std::endl;
  std::cout << " POTENTIAL: V=" << Incoming.Pot.V << " R0=" << Incoming.Pot.R0
            << " A=" << Incoming.Pot.A << std::endl;
  std::cout << "            VI=" << Incoming.Pot.VI
            << " RI0=" << Incoming.Pot.RI0 << " AI=" << Incoming.Pot.AI
            << std::endl;
  std::cout << "            VSI=" << Incoming.Pot.VSI
            << " RSI0=" << Incoming.Pot.RSI0 << " ASI=" << Incoming.Pot.ASI
            << std::endl;
  std::cout << "            VSO=" << Incoming.Pot.VSO
            << " RSO0=" << Incoming.Pot.RSO0 << " ASO=" << Incoming.Pot.ASO
            << std::endl;
  std::cout << "            RC0=" << Incoming.Pot.RC0 << std::endl;

  std::cout << std::endl;
  std::cout << "        OPTICAL MODEL SCATTERING FOR THE OUTGOING CHANNEL"
            << std::endl;
  std::cout << " POTENTIAL: V=" << Outgoing.Pot.V << " R0=" << Outgoing.Pot.R0
            << " A=" << Outgoing.Pot.A << std::endl;
  std::cout << "            VI=" << Outgoing.Pot.VI
            << " RI0=" << Outgoing.Pot.RI0 << " AI=" << Outgoing.Pot.AI
            << std::endl;
  std::cout << "            VSI=" << Outgoing.Pot.VSI
            << " RSI0=" << Outgoing.Pot.RSI0 << " ASI=" << Outgoing.Pot.ASI
            << std::endl;
  std::cout << "            VSO=" << Outgoing.Pot.VSO
            << " RSO0=" << Outgoing.Pot.RSO0 << " ASO=" << Outgoing.Pot.ASO
            << std::endl;
  std::cout << "            RC0=" << Outgoing.Pot.RC0 << std::endl;
  std::cout
      << "================================================================="
      << std::endl;
}

void DWBA::WavSet(Channel &ch) {
  // Determine grid parameters
  // Ptolemy uses heuristics or input. For now, use reasonable defaults or
  // simple heuristic. StepSize: 0.1 fm is standard. MaxR: R0 + 4*a + 10 fm? Or
  // based on Coulomb turning point? Let's use a fixed large range for now,
  // e.g., 20 fm, step 0.1. InFileCreator uses 0.05 or 0.1 usually.

  ch.StepSize = 0.1;
  ch.MaxR = 30.0; // Sufficient for most low energy reactions
  ch.NSteps = static_cast<int>(ch.MaxR / ch.StepSize) + 1;

  ch.RGrid.resize(ch.NSteps);
  ch.V_real.resize(ch.NSteps);
  ch.V_imag.resize(ch.NSteps);
  ch.V_so_real.resize(ch.NSteps);
  ch.V_so_imag.resize(ch.NSteps);
  ch.V_coulomb.resize(ch.NSteps);

  double A_target = ch.Target.A; // Mass number for radius scaling

  for (int i = 0; i < ch.NSteps; ++i) {
    double r = i * ch.StepSize;
    ch.RGrid[i] = r;

    if (r < 0.001) {
      // Avoid singularity at r=0
      ch.V_real[i] = 0.0; // Or potential at 0
      ch.V_imag[i] = 0.0;
      ch.V_so_real[i] = 0.0;
      ch.V_so_imag[i] = 0.0;
      ch.V_coulomb[i] = 0.0; // Diverges?
      // Coulomb potential for sphere is finite at 0: (3 - 0) / 2Rc ...
      // EvaluatePotential handles r < Rc correctly.
    }

    EvaluatePotential(r, ch.Pot, ch.V_real[i], ch.V_imag[i], ch.V_so_real[i],
                      ch.V_so_imag[i], ch.V_coulomb[i], ch.Projectile.Z,
                      ch.Target.Z, A_target);
  }
}

void DWBA::WavPot(Channel &ch) {
  // Allocate S-matrix array
  // Size depends on Lmax. For now, assume a max L or dynamic.
  // Let's say we support up to L=100.
  int MaxL = 100; // Should be dynamic based on Lmax
  ch.SMatrix.resize(MaxL + 1);
}

#include "rcwfn.h"

void DWBA::WavElj(Channel &ch, int L, int Jp) {
  // Solve radial equation using Numerov method
  // u''(r) + f(r) u(r) = 0
  // f(r) = k^2 - V(r) - L(L+1)/r^2
  // Numerov: (1 + h^2/12 f_{n+1}) u_{n+1} = 2(1 - 5h^2/12 f_n) u_n - (1 +
  // h^2/12 f_{n-1}) u_{n-1}

  int N = ch.NSteps;
  double h = ch.StepSize;
  double h2 = h * h;
  double h2_12 = h2 / 12.0;

  // Resize WaveFunction vector
  ch.WaveFunction.resize(N + 1);

  // Initial conditions
  // u(0) = 0
  // u(h) = C * h^(L+1)
  ch.WaveFunction[0] = 0.0;
  double start_val = std::pow(h, L + 1);
  // Avoid underflow/overflow if L is large
  if (start_val < 1e-30)
    start_val = 1e-30;
  ch.WaveFunction[1] = start_val;

  // Precompute f(r)
  std::vector<std::complex<double>> f(N + 1);
  double k2 = ch.k * ch.k;

  // Jp is 2*J.
  double spin_dot_L = 0.0;
  if (Jp == 2 * L + 1) {
    spin_dot_L = 0.5 * L;
  } else if (Jp == 2 * L - 1) {
    spin_dot_L = -0.5 * (L + 1);
  }

  // Constants for potential conversion
  // 2*mu/hbar^2
  // mu is in AMU. V is in MeV.
  // We need f(r) in fm^-2.
  // 1 AMU = 931.494 MeV/c^2
  // hbar*c = 197.327 MeV*fm
  // 2 * mu * c^2 / (hbar*c)^2 * V
  // = 2 * 931.494 / (197.327)^2 * V
  // = 0.0478 * V
  const double AMU_MEV = 931.494061;
  const double HBARC = 197.32697;
  double f_conv = 2.0 * ch.mu * AMU_MEV / (HBARC * HBARC);

  for (int i = 1; i <= N; ++i) {
    // double r = i * h; // Already have ch.RGrid[i]
    // Potential V(r)

    std::complex<double> V_central(ch.V_real[i], ch.V_imag[i]);
    V_central += ch.V_coulomb[i];

    std::complex<double> V_spin(ch.V_so_real[i], ch.V_so_imag[i]);

    std::complex<double> V_total = V_central + spin_dot_L * V_spin;

    // f(r) = k^2 - 2*mu/hbar^2 * V - L(L+1)/r^2
    // Note: ch.RGrid[i] should be i*h.
    double r = ch.RGrid[i];
    if (r < 1e-10)
      r = 1e-10; // Avoid division by zero

    f[i] = k2 - f_conv * V_total - (double)L * (L + 1) / (r * r);
  }

  // Integration
  for (int i = 1; i < N; ++i) {
    std::complex<double> term1 =
        2.0 * (1.0 - 5.0 * h2_12 * f[i]) * ch.WaveFunction[i];
    std::complex<double> term2 =
        (1.0 + h2_12 * f[i - 1]) * ch.WaveFunction[i - 1];
    std::complex<double> term3_inv = 1.0 / (1.0 + h2_12 * f[i + 1]);

    ch.WaveFunction[i + 1] = (term1 - term2) * term3_inv;

    // Check for overflow
    if (std::abs(ch.WaveFunction[i + 1]) > 1e100) {
      // Renormalize
      for (int j = 0; j <= i + 1; ++j) {
        ch.WaveFunction[j] *= 1e-100;
      }
    }
  }

  // S-matrix extraction
  double R_last = ch.MaxR;
  double R_prev = ch.MaxR - h;

  std::vector<double> FC(L + 2), FCP(L + 2), GC(L + 2), GCP(L + 2);
  std::vector<double> FC1(L + 2), FCP1(L + 2), GC1(L + 2), GCP1(L + 2);

  if (L == 0) {
    std::cout << "WavElj Debug: R_last=" << R_last << " R_prev=" << R_prev
              << " k=" << ch.k << " eta=" << ch.eta
              << " rho_last=" << ch.k * R_last << std::endl;
  }

  int ret1 = Rcwfn(ch.k * R_last, ch.eta, L, L, FC, FCP, GC, GCP);
  int ret2 = Rcwfn(ch.k * R_prev, ch.eta, L, L, FC1, FCP1, GC1, GCP1);

  if (ret1 != 0 || ret2 != 0) {
    std::cout << "Rcwfn Error: L=" << L << " k=" << ch.k << " eta=" << ch.eta
              << " ret1=" << ret1 << " ret2=" << ret2 << std::endl;
  }

  double F = FC[L];
  double G = GC[L];
  double F1 = FC1[L];
  double G1 = GC1[L];

  std::complex<double> u_last = ch.WaveFunction[N];
  std::complex<double> u_prev = ch.WaveFunction[N - 1];

  double wr_last = u_last.real();
  double wi_last = u_last.imag();
  double wr_prev = u_prev.real();
  double wi_prev = u_prev.imag();

  double A1 = wr_last * F1 + wi_last * G1 - wr_prev * F - wi_prev * G;
  double A2 = -wr_last * G1 + wi_last * F1 + wr_prev * G - wi_prev * F;
  double CR = -wr_last * F1 + wi_last * G1 + wr_prev * F - wi_prev * G;
  double CI = -wr_last * G1 - wi_last * F1 + wr_prev * G + wi_prev * F;

  double den = A1 * A1 + A2 * A2;
  double SJR = (CR * A1 + CI * A2) / den;
  double SJI = (CI * A1 - CR * A2) / den;

  // Ensure SMatrix is sized
  if (ch.SMatrix.size() <= L) {
    ch.SMatrix.resize(L + 1);
  }
  ch.SMatrix[L] = std::complex<double>(SJR, SJI);

  // Normalization
  double norm_A1 = 0.5 * (F * (1.0 + SJR) + SJI * G);
  double norm_A2 = 0.5 * (G * (1.0 - SJR) + SJI * F);

  double norm_den = wr_last * wr_last + wi_last * wi_last;
  double alpha_r = (wr_last * norm_A1 + wi_last * norm_A2) / norm_den;
  double alpha_i = (wr_last * norm_A2 - wi_last * norm_A1) / norm_den;

  std::complex<double> alpha(alpha_r, alpha_i);

  for (auto &val : ch.WaveFunction) {
    val *= alpha;
  }

  if (std::isnan(SJR) || std::isnan(SJI)) {
    std::cout << "WavElj NaN: L=" << L << " k=" << ch.k << " eta=" << ch.eta
              << std::endl;
    std::cout << "u_last=" << u_last << " u_prev=" << u_prev << std::endl;
    std::cout << "F=" << F << " G=" << G << " F1=" << F1 << " G1=" << G1
              << std::endl;
  }
}

void DWBA::GrdSet() {
  // Set up integration grid for finite range DWBA
  // We need a grid in r_in (Incoming.RGrid), r_out (Outgoing.RGrid), and theta.
  // Incoming and Outgoing grids are already set in WavSet.
  // We just need to define the angular grid.
  // Gauss-Legendre quadrature for cos(theta) from -1 to 1.

  int NTheta = 40; // Fixed number of points for now
  GaussLegendre(NTheta, -1.0, 1.0, ThetaGrid, ThetaWeights);

  // Convert cos(theta) to theta?
  // No, we usually integrate over d(cos theta).
  // But we might need theta for Legendre polynomials P_L(cos theta).
  // ThetaGrid contains x = cos(theta).
}

void DWBA::InelDc() {
  // Finite Range DWBA Integration
  // T = \int d^3 r_a \int d^3 r_b \chi_b^*(r_b) <\psi_B | V | \psi_A>
  // \chi_a(r_a)

  // 1. Calculate Mass Coefficients (S1, T1, S2, T2)
  // Based on GRDSET logic
  // A(a,b)B. a = b + x. B = A + x.
  // r_a = vector A->a
  // r_b = vector B->b (approx A->b)
  // r_x = vector A->x (Target Bound State radius)
  // r_p = vector b->x (Projectile Bound State radius)

  // Masses
  double mA = Incoming.Target.Mass;
  double ma = Incoming.Projectile.Mass;
  double mb = Outgoing.Projectile.Mass;
  double mB = Outgoing.Target.Mass;
  double mx = ma - mb; // Approx mass of transferred particle

  // Coordinate transformations
  // r_x = (ma/mx) r_a - (mb/mx) r_b ? No.
  // Let's use Ptolemy definitions from GRDSET
  // S1, T1 for Target BS (r_x)
  // S2, T2 for Projectile BS (r_p)
  // r_x = S1 * r_a + T1 * r_b
  // r_p = S2 * r_a + T2 * r_b

  // From GRDSET:
  // TEMP = 1 / ( BRATMS(1) + BRATMS(2) * (1+BRATMS(1)) )
  // BRATMS(1) = x/A, BRATMS(2) = x/b (for stripping)
  // S1 = (1 + x/A) * (1 + x/b) * TEMP
  // T1 = - (1 + x/b) * TEMP
  // S2 = (1 + x/A) * TEMP
  // T2 = -S1

  double ratio_xA = mx / mA;
  double ratio_xb = mx / mb;
  double temp = 1.0 / (ratio_xA + ratio_xb * (1.0 + ratio_xA));

  double S1 = (1.0 + ratio_xA) * (1.0 + ratio_xb) * temp;
  double T1 = -(1.0 + ratio_xb) * temp;
  double S2 = (1.0 + ratio_xA) * temp;
  double T2 = -S1; // For stripping

  // 2. Calculate Bound States

  // Refine mx (transferred particle mass) to include binding energy
  // ma = mb + mx - BindingEnergy -> mx = ma - mb + BindingEnergy
  // This recovers the actual neutron mass instead of (d-p) mass difference.
  const double AMU_MEV = 931.494061;
  mx = (ma - mb) + ProjectileBS.BindingEnergy / AMU_MEV;

  // Target Bound State (x in B)
  Channel TargetBSChannel = Incoming; // Copy properties
  TargetBSChannel.Pot = TargetBS.Pot;
  TargetBSChannel.mu = mA * mx / (mA + mx); // Reduced mass of x+A

  // Fix Target/Projectile for BS
  // Core is Target A (Incoming.Target).
  // Ptolemy uses Core Mass (33) for radius (R=4.009).
  // So we use Incoming.Target.
  TargetBSChannel.Target = Incoming.Target;
  // Particle is x (transferred particle).
  TargetBSChannel.Projectile.Z = Incoming.Projectile.Z - Outgoing.Projectile.Z;
  TargetBSChannel.Projectile.A = Incoming.Projectile.A - Outgoing.Projectile.A;
  TargetBSChannel.Projectile.Mass = mx;

  TargetBSChannel.NSteps = Incoming.NSteps;
  TargetBSChannel.StepSize = Incoming.StepSize;
  TargetBSChannel.RGrid = Incoming.RGrid;
  TargetBSChannel.WaveFunction.resize(Incoming.NSteps);
  TargetBSChannel.V_real.resize(Incoming.NSteps); // Resize potential arrays
  TargetBSChannel.V_imag.resize(Incoming.NSteps);
  TargetBSChannel.V_so_real.resize(Incoming.NSteps);
  TargetBSChannel.V_so_imag.resize(Incoming.NSteps);
  TargetBSChannel.V_coulomb.resize(Incoming.NSteps);

  // Setup potential arrays for bound state
  WavSet(TargetBSChannel);
  CalculateBoundState(TargetBSChannel, TargetBS.n, TargetBS.l, TargetBS.j,
                      TargetBS.BindingEnergy);

  // Projectile Bound State (x in a)
  Channel ProjectileBSChannel = Incoming;
  ProjectileBSChannel.Pot = ProjectileBS.Pot;
  ProjectileBSChannel.mu = mb * mx / (mb + mx); // Reduced mass of x+b

  // Fix Target/Projectile for BS
  // Core is b (Ejectile).
  ProjectileBSChannel.Target = Outgoing.Projectile;
  // Particle is x.
  ProjectileBSChannel.Projectile.Z =
      Incoming.Projectile.Z - Outgoing.Projectile.Z;
  ProjectileBSChannel.Projectile.A =
      Incoming.Projectile.A - Outgoing.Projectile.A;
  ProjectileBSChannel.Projectile.Mass = mx;

  ProjectileBSChannel.NSteps = Incoming.NSteps; // Use same grid for simplicity
  ProjectileBSChannel.StepSize = Incoming.StepSize;
  ProjectileBSChannel.RGrid = Incoming.RGrid;
  ProjectileBSChannel.WaveFunction.resize(Incoming.NSteps);
  ProjectileBSChannel.V_real.resize(Incoming.NSteps);
  ProjectileBSChannel.V_imag.resize(Incoming.NSteps);
  ProjectileBSChannel.V_so_real.resize(Incoming.NSteps);
  ProjectileBSChannel.V_so_imag.resize(Incoming.NSteps);
  ProjectileBSChannel.V_coulomb.resize(Incoming.NSteps);

  WavSet(ProjectileBSChannel);
  CalculateBoundState(ProjectileBSChannel, ProjectileBS.n, ProjectileBS.l,
                      ProjectileBS.j, ProjectileBS.BindingEnergy);

  // Debug: Check Bound State
  std::cout << "Target BS: Binding=" << TargetBS.BindingEnergy
            << " Norm=" << TargetBSChannel.WaveFunction[0] << std::endl;
  std::cout << "Projectile BS: Binding=" << ProjectileBS.BindingEnergy
            << " Norm=" << ProjectileBSChannel.WaveFunction[0] << std::endl;

  // 3. Loop over Partial Waves
  // For simplicity, let's just do L=0 to Lmax
  int Lmax = 5; // Reduced for testing speed (was 20)

  TransferSMatrix.clear();

  std::cout << "Calculating Transfer Integrals..." << std::endl;

  for (int La = 0; La <= Lmax; ++La) {
    // Calculate Incoming Distorted Wave
    WavElj(Incoming, La, 2 * La + 1); // Spin zero approx for now

    // For each La, we have a range of Lb allowed by selection rules
    // L_trans = |La - Lb| (approx)
    // L_trans is determined by bound state angular momenta l_T and l_P.
    // Vector sum: l_T = l_P + L_trans (or similar)
    // Let's iterate Lb
    for (int Lb = std::abs(La - TargetBS.l - ProjectileBS.l);
         Lb <= La + TargetBS.l + ProjectileBS.l; ++Lb) {
      if (Lb < 0)
        continue;

      // Calculate Outgoing Distorted Wave
      WavElj(Outgoing, Lb, 2 * Lb + 1);

      // 4. Perform Integration
      std::complex<double> Integral(0.0, 0.0);

      for (int i = 0; i < Incoming.NSteps; ++i) {
        double ra = Incoming.RGrid[i];
        if (ra < 1e-3)
          continue;
        std::complex<double> chi_a = Incoming.WaveFunction[i];

        for (int j = 0; j < Outgoing.NSteps; ++j) {
          double rb = Outgoing.RGrid[j];
          if (rb < 1e-3)
            continue;
          std::complex<double> chi_b = Outgoing.WaveFunction[j];

          // Angular Integration
          std::complex<double> AngInt(0.0, 0.0);

          for (size_t k = 0; k < ThetaGrid.size(); ++k) {
            double cos_theta = ThetaGrid[k];
            double weight = ThetaWeights[k];

            // Calculate r_x and r_p lengths
            double rx2 = S1 * S1 * ra * ra + T1 * T1 * rb * rb +
                         2 * S1 * T1 * ra * rb * cos_theta;
            double rp2 = S2 * S2 * ra * ra + T2 * T2 * rb * rb +
                         2 * S2 * T2 * ra * rb * cos_theta;

            if (rx2 < 0)
              rx2 = 0;
            if (rp2 < 0)
              rp2 = 0;

            double rx = std::sqrt(rx2);
            double rp = std::sqrt(rp2);

            // Interpolate Bound States
            // Simple linear interpolation
            double val_T = 0.0;
            double val_P = 0.0;

            // Helper lambda for interpolation
            auto Interpolate = [&](const Channel &ch, double r) -> double {
              if (r >= ch.MaxR)
                return 0.0;
              if (r < 0)
                return 0.0;
              double idx = r / ch.StepSize;
              int idx_i = static_cast<int>(idx);
              if (idx_i >= ch.NSteps - 1)
                return 0.0;
              double frac = idx - idx_i;
              return ch.WaveFunction[idx_i].real() * (1.0 - frac) +
                     ch.WaveFunction[idx_i + 1].real() * frac;
            };

            val_T = Interpolate(TargetBSChannel, rx);
            val_P = Interpolate(ProjectileBSChannel, rp);

            // Evaluate Interaction V(rp)
            // Prior form: V_bx(rp).
            // We can use EvaluatePotential or just use the bound state
            // potential depth? The bound state potential is V_bx. But
            // EvaluatePotential calculates V at a point. We need to use the
            // ProjectileBS.Pot parameters.
            double V_real, V_imag, V_so_r, V_so_i, V_c;
            EvaluatePotential(rp, ProjectileBS.Pot, V_real, V_imag, V_so_r,
                              V_so_i, V_c, 0, 0, 0); // Z=0 for nuclear only

            double Interaction = V_real; // Real central part

            // Form Factor
            double FormFactor = val_T * Interaction * val_P;

            // Angular part (Legendre Polynomials?)
            // We need to project onto Y_La and Y_Lb.
            // This is complex.
            // For now, assume monopole transfer (L=0) implies angle
            // independent? No, finite range makes it angle dependent. We
            // integrate F(ra, rb, theta) * P_L(cos theta)? The transfer
            // amplitude involves summing over M.

            // Let's assume a simplified kernel K(ra, rb) = \int d(cos theta)
            // F(...) But we need to couple La and Lb. For L=0 transfer, La =
            // Lb. Let's assume we are calculating the radial integral for a
            // specific multipole. For now, just integrate F.

            AngInt += FormFactor * weight;
          }

          Integral += chi_a * std::conj(chi_b) * AngInt * (ra * ra) *
                      (rb * rb) * (Incoming.StepSize * Outgoing.StepSize);
        }
      }

      // ... (calculation of Integral)

      // Debug: Check for NaN
      if (std::isnan(Integral.real()) || std::isnan(Integral.imag())) {
        std::cout << "NaN Integral at La=" << La << " Lb=" << Lb << std::endl;
        std::cout << "Incoming S: " << Incoming.SMatrix[La] << std::endl;
        std::cout << "Outgoing S: " << Outgoing.SMatrix[Lb] << std::endl;
      }

      TransferSMatrix.push_back({La, Lb, Integral});
    }
  }
}

void DWBA::XSectn() {
  // Calculate Cross Section
  // dsigma/dOmega = (mu_a mu_b / (2 pi hbar^2)^2) * (kb/ka) * |T|^2
  // T = \sum_{La, Lb} S_{La, Lb} * GeometricFactors * Y_{Lb} * Y_{La}^*

  // Constants
  const double HBARC = 197.32697;
  const double AMU = 931.494;

  // Kinematic factors
  double mu_a = Incoming.mu; // Reduced mass a+A
  double mu_b = Outgoing.mu; // Reduced mass b+B
  double ka = Incoming.k;
  double kb = Outgoing.k;

  // Factor = (mu_a * mu_b) / (2 * pi * hbar^2)^2 * (kb / ka)
  // mu is in AMU. We need energy units?
  // T-matrix usually has units of MeV^-1 fm^3?
  // Let's check standard DWBA formulas.
  // dsigma/dOmega (mb/sr) = 10 * (mu_a mu_b / (2 pi hbar^2)^2) * (kb/ka) *
  // |T|^2 if T is in MeV fm^3. My Integral in InelDc is \int chi_b^* V chi_a
  // d^3r. V is in MeV. chi are dimensionless. r^3 is fm^3. So Integral is MeV
  // fm^3. 2*pi*hbar^2 = 2*pi * (hbar*c)^2 / c^2? No. hbar^2/2mu = 20.736 MeV
  // fm^2 / mu(amu). So 2*pi*hbar^2 = 4*pi*mu * (hbar^2/2mu).

  // Let's use the prefactor:
  // C = (mu_a * mu_b) / (4 * pi^2 * hbar^4) * (kb/ka)
  // hbar^4 = (hbar c)^4 / c^4.
  // mu * c^2 is mass in MeV.
  // C = (mu_a c^2 * mu_b c^2) / (4 * pi^2 * (hbar c)^4) * (kb/ka)

  double muc2_a = mu_a * AMU;
  double muc2_b = mu_b * AMU;
  double hbarc4 = std::pow(HBARC, 4);
  double prefactor =
      (muc2_a * muc2_b) / (4.0 * M_PI * M_PI * hbarc4) * (kb / ka);

  // Convert to mb/sr. 1 fm^2 = 10 mb.
  prefactor *= 10.0;

  std::cout << "Calculating Cross Sections..." << std::endl;

  // Debug: Print Transfer S-Matrix
  std::cout << "Transfer S-Matrix:" << std::endl;
  std::cout << "La  Lb  Re(S)        Im(S)        Abs(S)" << std::endl;
  for (const auto &elem : TransferSMatrix) {
    std::cout << std::setw(3) << elem.La << " " << std::setw(3) << elem.Lb
              << " " << std::scientific << std::setprecision(4) << elem.S.real()
              << " " << elem.S.imag() << " " << std::abs(elem.S) << std::endl;
  }

  std::cout << "Angle (deg)   dSigma/dOmega (mb/sr)" << std::endl;

  for (double theta_deg = AngleMin; theta_deg <= AngleMax;
       theta_deg += AngleStep) {
    double theta_rad = theta_deg * M_PI / 180.0;
    double cos_theta = std::cos(theta_rad);

    std::complex<double> Amplitude(0.0, 0.0);

    for (const auto &elem : TransferSMatrix) {
      int La = elem.La;
      int Lb = elem.Lb;
      std::complex<double> S = elem.S;

      // We need Y_{La, 0}(theta, 0) and Y_{Lb, 0}(theta, 0)
      // assuming m=0 transfer dominates or we sum over m properly.
      // For L=0 transfer, La=Lb.
      // Y_{L,0} = sqrt((2L+1)/4pi) P_L(cos theta)

      // Legendre Polynomials
      // We can use std::legendre (C++17) or implement it.
      // math_utils doesn't have it yet.
      // But C++17 has std::legendre in <cmath>.

      double P_La = std::legendre(La, cos_theta);
      double P_Lb = std::legendre(Lb, cos_theta);

      double Y_La = std::sqrt((2.0 * La + 1.0) / (4.0 * M_PI)) * P_La;
      double Y_Lb = std::sqrt((2.0 * Lb + 1.0) / (4.0 * M_PI)) * P_Lb;

      // Phase factor i^(La - Lb)?
      // Usually included in S or expansion.
      // Let's assume S includes radial phases.
      // We sum coherent contributions.

      // Geometric factor for angular momentum coupling?
      // For L=0 transfer, it's simple.
      // For L>0, we need Clebsch-Gordan.
      // Let's assume simplified sum for now:
      // T ~ sum S * Y_Lb * Y_La^*

      Amplitude += S * Y_Lb * Y_La;
    }

    double dSigma = prefactor * std::norm(Amplitude);

    std::cout << std::fixed << std::setprecision(1) << theta_deg
              << "           " << std::scientific << std::setprecision(4)
              << dSigma << std::endl;
  }
}

void DWBA::CalculateBoundState(Channel &ch, int n, int l, double j,
                               double bindingEnergy) {
  double targetEnergy = -bindingEnergy;
  double mu = ch.mu;

  // Constants
  const double AMU_MEV = 931.494061;
  const double HBARC = 197.32697;
  // ch.mu is in MeV, so we don't need to multiply by AMU_MEV
  double factor = 2.0 * mu / (HBARC * HBARC);
  double kappa = std::sqrt(2.0 * mu * bindingEnergy) / HBARC;

  // Grid
  int N = ch.NSteps;
  double h = ch.StepSize;
  double h2 = h * h;
  double h12 = h2 / 12.0;

  // Matching Radius
  // Use R0 * A^(1/3) + 2.0 fm as a heuristic
  double R_nuc = ch.Pot.R0 * std::pow(ch.Target.A, 1.0 / 3.0);
  int MatchIdx = static_cast<int>((R_nuc + 2.0) / h);
  if (MatchIdx >= N - 2)
    MatchIdx = N - 5;
  if (MatchIdx < 5)
    MatchIdx = 5;

  std::cout << "Bound State Calculation: E=" << targetEnergy << " MeV, L=" << l
            << ", J=" << j << ", MatchR=" << MatchIdx * h << " fm"
            << ", mu=" << mu << ", factor=" << factor << std::endl;

  // Iteration to find V
  double V_depth = ch.Pot.V;
  double V_step = 5.0;

  double LS = 0.0;
  if (j > 0) {
    double s = 0.5;
    LS = (j * (j + 1) - l * (l + 1) - s * (s + 1)) / 2.0;
  }

  std::vector<double> u_out(N);
  std::vector<double> u_in(N);
  std::vector<double> f(N);
  int nodes = 0;

  for (int iter = 0; iter < 500; ++iter) {
    // 1. Calculate Potential and f(r)
    for (int i = 0; i < N; ++i) {
      double r = i * h; // Use explicit grid
      if (r < 1e-5)
        r = 1e-5;

      double V_r, V_i, V_so, V_soi, V_c;
      ChannelPotential pot = ch.Pot; // Local copy for evaluation
      pot.V = V_depth;
      EvaluatePotential(r, pot, V_r, V_i, V_so, V_soi, V_c, ch.Projectile.Z,
                        ch.Target.Z, ch.Target.A);

      double V_tot = -V_r + V_so * LS + V_c; // Attractive V_r
      f[i] = factor * (targetEnergy - V_tot) -
             static_cast<double>(l * (l + 1)) / (r * r);
    }

    // 2. Outward Integration (0 -> MatchIdx)
    std::fill(u_out.begin(), u_out.end(), 0.0);
    u_out[0] = 0.0;
    u_out[1] = std::pow(h, l + 1) * 1e-5;

    nodes = 0;
    for (int i = 1; i < MatchIdx; ++i) {
      double term1 = 2.0 * (1.0 - 5.0 * h12 * f[i]) * u_out[i];
      double term2 = (1.0 + h12 * f[i - 1]) * u_out[i - 1];
      double term3 = (1.0 + h12 * f[i + 1]);
      u_out[i + 1] = (term1 - term2) / term3;

      if (u_out[i + 1] * u_out[i] < 0)
        nodes++;

      if (std::abs(u_out[i + 1]) > 1e10) {
        for (int k = 0; k <= i + 1; ++k)
          u_out[k] *= 1e-10;
      }
    }

    // 3. Inward Integration (N-1 -> MatchIdx)
    std::fill(u_in.begin(), u_in.end(), 0.0);
    u_in[N - 1] = std::exp(-kappa * (N - 1) * h);
    u_in[N - 2] = std::exp(-kappa * (N - 2) * h);

    for (int i = N - 2; i > MatchIdx; --i) {
      double term1 = 2.0 * (1.0 - 5.0 * h12 * f[i]) * u_in[i];
      double term2 = (1.0 + h12 * f[i + 1]) * u_in[i + 1];
      double term3 = (1.0 + h12 * f[i - 1]);
      u_in[i - 1] = (term1 - term2) / term3;

      if (std::abs(u_in[i - 1]) > 1e10) {
        for (int k = N - 1; k >= i - 1; --k)
          u_in[k] *= 1e-10;
      }
    }

    // 4. Matching
    double scale = u_out[MatchIdx] / u_in[MatchIdx];
    for (int i = MatchIdx; i < N; ++i)
      u_in[i] *= scale;

    // Calculate Logarithmic Derivatives
    double term1_out = 2.0 * (1.0 - 5.0 * h12 * f[MatchIdx]) * u_out[MatchIdx];
    double term2_out = (1.0 + h12 * f[MatchIdx - 1]) * u_out[MatchIdx - 1];
    double term3_out = (1.0 + h12 * f[MatchIdx + 1]);
    double u_out_next = (term1_out - term2_out) / term3_out;

    double d_out_val = (u_out_next - u_out[MatchIdx - 1]) / (2.0 * h);

    double term1_in = 2.0 * (1.0 - 5.0 * h12 * f[MatchIdx]) * u_in[MatchIdx];
    double term2_in = (1.0 + h12 * f[MatchIdx + 1]) * u_in[MatchIdx + 1];
    double term3_in = (1.0 + h12 * f[MatchIdx - 1]);
    double u_in_prev = (term1_in - term2_in) / term3_in;

    double d_in_val = (u_in[MatchIdx + 1] - u_in_prev) / (2.0 * h);

    double log_deriv_out = d_out_val / u_out[MatchIdx];
    double log_deriv_in = d_in_val / u_in[MatchIdx];

    double Diff = log_deriv_out - log_deriv_in;

    if (iter % 50 == 0) {
      std::cout << "BS Iter " << iter << ": V=" << V_depth << " Nodes=" << nodes
                << " Diff=" << Diff << std::endl;
    }

    // 5. Update V
    if (nodes != n) {
      if (nodes > n)
        V_depth -= V_step;
      else
        V_depth += V_step;
      V_step *= 0.8;
      if (V_step < 0.1)
        V_step = 0.1;
    } else {
      static double update_factor = 2.0;
      static double prev_Diff = 0.0;

      // Reset factor if we just entered the correct node region
      if (std::abs(prev_Diff) < 1e-10)
        update_factor = 2.0;

      if (iter > 0 && Diff * prev_Diff < 0) {
        update_factor *= 0.5;
      }
      prev_Diff = Diff;

      V_depth += Diff * update_factor;

      if (std::abs(Diff) < 1e-5)
        break;
    }
  }

  std::cout << "Final Bound State V: " << V_depth << " MeV" << std::endl;

  // Stitch and Normalize
  std::vector<double> u(N);
  for (int i = 0; i < MatchIdx; ++i)
    u[i] = u_out[i];
  for (int i = MatchIdx; i < N; ++i)
    u[i] = u_in[i];

  double norm = 0.0;
  for (int i = 0; i < N; ++i)
    norm += u[i] * u[i];
  norm *= h;

  if (norm > 0) {
    double scale = 1.0 / std::sqrt(norm);
    double max_amp = 0.0;
    for (int i = 0; i < N; ++i) {
      ch.WaveFunction[i] = u[i] * scale;
      if (std::abs(ch.WaveFunction[i].real()) > max_amp)
        max_amp = std::abs(ch.WaveFunction[i].real());
    }
    ch.Pot.V = V_depth; // Update potential in channel
    std::cout << "BS Final: V=" << V_depth << " Norm=" << norm
              << " MaxAmp=" << max_amp << std::endl;
  } else {
    std::cerr << "Error: Bound state normalization failed (norm=0)."
              << std::endl;
  }
}
