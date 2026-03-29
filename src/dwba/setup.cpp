// setup.cpp — DWBA constructor, destructor, setters, Calculate, PrintParameters
// Extracted from dwba.cpp (lines 1–180, 249–341)

#include "dwba.h"
#include "math_utils.h"
#include "potential_eval.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

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

// ---------------------------------------------------------------
// Constructor / Destructor
// ---------------------------------------------------------------

DWBA::DWBA() {
  AngleMin = 0.0;
  AngleMax = 180.0;
  AngleStep = 0.5;
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

// ---------------------------------------------------------------
// Setters
// ---------------------------------------------------------------

void DWBA::SetReaction(const std::string &target, const std::string &projectile,
                       const std::string &ejectile, const std::string &recoil) {
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

void DWBA::SetTargetSpin(double J)    { SpinTarget   = J; }
void DWBA::SetResidualSpin(double J)  { SpinResidual = J; }

void DWBA::SetProjectileWFFile(const std::string &filename, double grid_h, double spam) {
  // Load tabulated phi=u/r from file with columns: step  r  phi
  // Ptolemy prints on a subsampled grid; we interpolate onto grid_h
  // File format: "   STEP    RADIUS    WAVEFUNCTION"
  std::ifstream f(filename);
  if (!f.is_open()) {
    return;
  }
  std::vector<double> r_tab, phi_tab;
  int step_dummy; double r_val, phi_val;
  while (f >> step_dummy >> r_val >> phi_val) {
    r_tab.push_back(r_val);
    phi_tab.push_back(phi_val);
  }
  f.close();

  // Interpolate onto uniform grid with step grid_h, from 0 to r_tab.back()
  int N = static_cast<int>(r_tab.back() / grid_h) + 2;
  ProjectileWFTable.assign(N, 0.0);
  for (int I = 0; I < N; I++) {
    double r = I * grid_h;
    // Linear interpolation
    auto it = std::lower_bound(r_tab.begin(), r_tab.end(), r);
    if (it == r_tab.end()) { ProjectileWFTable[I] = 0.0; continue; }
    if (it == r_tab.begin()) { ProjectileWFTable[I] = phi_tab[0]; continue; }
    int j2 = std::distance(r_tab.begin(), it);
    int j1 = j2 - 1;
    double frac = (r - r_tab[j1]) / (r_tab[j2] - r_tab[j1]);
    ProjectileWFTable[I] = phi_tab[j1] + frac*(phi_tab[j2]-phi_tab[j1]);
  }
  ProjectileWFGridH  = grid_h;
  ProjectileWFSpam   = spam;
  ProjectileWFLoaded = true;
  std::cout << "Loaded projectile WF: " << N << " points, h=" << grid_h
            << " fm, SPAM=" << spam << "\n";
}

void DWBA::SetProjectileWFFromFile(const std::vector<std::pair<double,double>>& wf_data,
                                   double h_ext, double spam) {
  if (wf_data.empty()) return;
  // Interpolate onto uniform grid
  int N = static_cast<int>(wf_data.back().first / h_ext) + 2;
  ProjectileWFTable.assign(N, 0.0);
  for (int I = 0; I < N; I++) {
    double r = I * h_ext;
    // Find bracket
    int j2 = -1;
    for (int k = 1; k < (int)wf_data.size(); k++) {
      if (wf_data[k].first >= r) { j2 = k; break; }
    }
    if (j2 <= 0) { ProjectileWFTable[I] = wf_data[0].second; continue; }
    int j1 = j2 - 1;
    double frac = (r - wf_data[j1].first) / (wf_data[j2].first - wf_data[j1].first);
    ProjectileWFTable[I] = wf_data[j1].second + frac*(wf_data[j2].second - wf_data[j1].second);
  }
  ProjectileWFGridH = h_ext;
  ProjectileWFSpam = spam;
  ProjectileWFLoaded = true;
  std::cout << "Loaded projectile WF from parsed data: " << N << " points, h=" << h_ext
            << " fm, SPAM=" << spam << "\n";
}

// ---------------------------------------------------------------
// Calculate — main orchestration function
// ---------------------------------------------------------------

void DWBA::Calculate() {
  CalculateKinematics();

  // Default nuclear spins if not explicitly set
  if (SpinTarget < -0.5) {
    // Heuristic: even-even nucleus → J=0, else use bound state j
    bool target_ee = (Incoming.Target.A % 2 == 0) && (Incoming.Target.Z % 2 == 0);
    SpinTarget = target_ee ? 0.0 : TargetBS.j;
    std::cout << "[DWBA] SpinTarget defaulted to " << SpinTarget
              << (target_ee ? " (even-even)" : " (from TargetBS.j)") << "\n";
  }
  if (SpinResidual < -0.5) {
    bool residual_ee = (Outgoing.Target.A % 2 == 0) && (Outgoing.Target.Z % 2 == 0);
    SpinResidual = residual_ee ? 0.0 : TargetBS.j;
    std::cout << "[DWBA] SpinResidual defaulted to " << SpinResidual
              << (residual_ee ? " (even-even)" : " (from TargetBS.j)") << "\n";
  }

  PrintParameters();

  // Setup grids and potentials
  // Ptolemy step size for scattering channels: h = min(2π/k, 1.0) / STEPSPER (STEPSPER=8)
  const int STEPSPER = 8;
  Incoming.StepSize = std::min(2.0 * M_PI / Incoming.k, 1.0) / STEPSPER;
  Outgoing.StepSize = std::min(2.0 * M_PI / Outgoing.k, 1.0) / STEPSPER;

  // Ptolemy WAVSET: extend asymptopia to classical turning point for large L
  // RMAX = max(ASYMPT, (eta + sqrt(eta^2 + LMAX*(LMAX+1))) / k)
  // Condition: NUMFIT==0 && LMAX!=NOTDEF && SCTASY<0 (always true for DWBA with set LMAX)
  //
  // CRITICAL: In Ptolemy, the 'asymptopia' PARAMETERSET keyword goes to SCTASY
  // (the GRDSET integration range), NOT to WAVSET's ASYMPT.
  // WAVSET's ASYMPT comes from BNDASY (default 20 fm), set by the bound-state
  // section BEFORE the scattering channels. Since the bound state call to SETBNDS
  // sets ASYMPT = BNDASY = 20, and subsequent calls skip (ASYMPT != UNDEF),
  // WAVSET always uses BNDASY = 20 as the base asymptopia.
  {
    // Ptolemy SETBNDS (source.f line 33041-33042):
    //   IF (ASYMPT .NE. UNDEF) RETURN  — user keyword overrides
    //   ASYMPT = BNDASY           (if IBND != 0, i.e. bound state)
    //   ASYMPT = ABS(SCTASY)      (if IBND == 0, i.e. scattering)
    //
    // For the ACTUAL scattering wavefunctions (used by WavElj for S-matrix extraction),
    // we use the user's asymptopia to ensure Numerov extends far enough.
    // The SCAN chi (GRDSET BSPROD) uses a SEPARATE, shorter range (ABS(SctAsySet))
    // — computed in grdset_ineldc_faithful.cpp.
    double ASYMPT = (AsymptopiaSet > 0) ? AsymptopiaSet : std::abs(SctAsySet);
    int LMAX_eff = (LmaxSet >= 0) ? LmaxSet : 40;

    // Incoming channel turning point
    double tp_in = (Incoming.eta + std::sqrt(Incoming.eta * Incoming.eta
                    + (double)LMAX_eff * (LMAX_eff + 1))) / Incoming.k;
    Incoming.MaxR = std::max(ASYMPT, tp_in);
    // Snap to grid: NSTEP = RMAX/h + 0.5, then RMAX = NSTEP*h
    int nstep_in = static_cast<int>(Incoming.MaxR / Incoming.StepSize + 0.5);
    Incoming.MaxR = nstep_in * Incoming.StepSize;

    // Outgoing channel turning point
    double tp_out = (Outgoing.eta + std::sqrt(Outgoing.eta * Outgoing.eta
                     + (double)LMAX_eff * (LMAX_eff + 1))) / Outgoing.k;
    Outgoing.MaxR = std::max(ASYMPT, tp_out);
    int nstep_out = static_cast<int>(Outgoing.MaxR / Outgoing.StepSize + 0.5);
    Outgoing.MaxR = nstep_out * Outgoing.StepSize;

    std::cout << "[WAVSET] Incoming: ASYMPTOPIA = " << Incoming.MaxR
              << " fm, NSTEP = " << nstep_in
              << ", StepSize = " << Incoming.StepSize << " fm\n";
    std::cout << "[WAVSET] Outgoing: ASYMPTOPIA = " << Outgoing.MaxR
              << " fm, NSTEP = " << nstep_out
              << ", StepSize = " << Outgoing.StepSize << " fm\n";
  }

  WavSet(Incoming);
  WavSet(Outgoing);

  // Setup Angular Integration Grid
  GrdSet();

  // Perform Finite Range Integration -> fills TransferSMatrix with S(Lx, Li, Lo)
  InelDc();

  // Calculate Cross Sections from Transfer S-Matrix
  XSectn();
}

void DWBA::CalculateZR() {
  CalculateKinematics();

  // Default nuclear spins if not explicitly set
  if (SpinTarget < -0.5) {
    bool target_ee = (Incoming.Target.A % 2 == 0) && (Incoming.Target.Z % 2 == 0);
    SpinTarget = target_ee ? 0.0 : TargetBS.j;
    std::cout << "[DWBA] SpinTarget defaulted to " << SpinTarget
              << (target_ee ? " (even-even)" : " (from TargetBS.j)") << "\n";
  }
  if (SpinResidual < -0.5) {
    bool residual_ee = (Outgoing.Target.A % 2 == 0) && (Outgoing.Target.Z % 2 == 0);
    SpinResidual = residual_ee ? 0.0 : TargetBS.j;
    std::cout << "[DWBA] SpinResidual defaulted to " << SpinResidual
              << (residual_ee ? " (even-even)" : " (from TargetBS.j)") << "\n";
  }

  PrintParameters();

  // Setup grids and potentials (Ptolemy step size formula)
  const int STEPSPER_ZR = 8;
  Incoming.StepSize = std::min(2.0 * M_PI / Incoming.k, 1.0) / STEPSPER_ZR;
  Outgoing.StepSize = std::min(2.0 * M_PI / Outgoing.k, 1.0) / STEPSPER_ZR;
  WavSet(Incoming);
  WavSet(Outgoing);

  // Setup Angular Integration Grid (still needed for ThetaGrid)
  GrdSet();

  // Zero-Range transfer integral (replaces InelDc)
  InelDcZR();

  // Cross sections (same as FR)
  XSectn();
}

// ---------------------------------------------------------------
// PrintParameters
// ---------------------------------------------------------------

void DWBA::PrintParameters() {
  std::cout
      << "================================================================="
      << std::endl;
  const double AMU_MEV = 931.494061;
  const double HBARC_L = 197.32697;

  double ma = Incoming.Projectile.Mass;
  double mb = Outgoing.Projectile.Mass;
  double mA = Incoming.Target.Mass;
  double mx = ma - mb;
  double mu_pbs = mb * mx / (mb + mx);
  double kappa_pbs = std::sqrt(2.0 * mu_pbs * ProjectileBS.BindingEnergy) / HBARC_L;

  std::cout << "        PROJECTILE BOUND STATE PARAMETERS" << std::endl;
  std::cout << "E = " << std::fixed << std::setprecision(4)
            << -ProjectileBS.BindingEnergy << " MEV"
            << "     KAPPA = " << std::setprecision(5) << kappa_pbs << std::endl;
  std::cout << " PROJECTILE MASS = " << std::fixed << std::setprecision(2)
            << mx / AMU_MEV << " AMU"
            << "     TARGET MASS = " << mb / AMU_MEV << " AMU"
            << "     REDUCED MASS = " << std::setprecision(2) << mu_pbs
            << " MEV/C**2" << std::endl;
  std::cout << " L = " << ProjectileBS.l << "  N = " << ProjectileBS.n << std::endl;
  std::cout << " POTENTIAL: V=" << std::setprecision(4) << ProjectileBS.Pot.V
            << " R0=" << ProjectileBS.Pot.R0 << " A=" << ProjectileBS.Pot.A << std::endl;
  std::cout << "            VSO=" << ProjectileBS.Pot.VSO
            << " RSO0=" << ProjectileBS.Pot.RSO0 << " ASO=" << ProjectileBS.Pot.ASO << std::endl;
  std::cout << "            RC0=" << ProjectileBS.Pot.RC0 << std::endl;
  std::cout << std::endl;

  double mu_tbs = mA * mx / (mA + mx);
  double kappa_tbs = std::sqrt(2.0 * mu_tbs * TargetBS.BindingEnergy) / HBARC_L;

  std::cout << "        TARGET BOUND STATE PARAMETERS" << std::endl;
  std::cout << "E = " << std::fixed << std::setprecision(4)
            << -TargetBS.BindingEnergy << " MEV"
            << "     KAPPA = " << std::setprecision(5) << kappa_tbs << std::endl;
  std::cout << " PROJECTILE MASS = " << std::fixed << std::setprecision(2)
            << mx / AMU_MEV << " AMU"
            << "     TARGET MASS = " << mA / AMU_MEV << " AMU"
            << "     REDUCED MASS = " << std::setprecision(2) << mu_tbs
            << " MEV/C**2" << std::endl;
  std::cout << " L = " << TargetBS.l << "  N = " << TargetBS.n << std::endl;
  std::cout << " POTENTIAL: V=" << std::setprecision(4) << TargetBS.Pot.V
            << " R0=" << TargetBS.Pot.R0 << " A=" << TargetBS.Pot.A << std::endl;
  std::cout << "            VSO=" << TargetBS.Pot.VSO
            << " RSO0=" << TargetBS.Pot.RSO0 << " ASO=" << TargetBS.Pot.ASO << std::endl;
  std::cout << "            RC0=" << TargetBS.Pot.RC0 << std::endl;
  std::cout << std::endl;

  std::cout << "        OPTICAL MODEL SCATTERING FOR THE INCOMING CHANNEL" << std::endl;
  std::cout << "E LAB = " << Incoming.Elab << " MEV" << std::endl;
  std::cout << " POTENTIAL: V=" << Incoming.Pot.V << " R0=" << Incoming.Pot.R0
            << " A=" << Incoming.Pot.A << std::endl;
  std::cout << "            VI=" << Incoming.Pot.VI << " RI0=" << Incoming.Pot.RI0
            << " AI=" << Incoming.Pot.AI << std::endl;
  std::cout << "            VSO=" << Incoming.Pot.VSO << " RSO0=" << Incoming.Pot.RSO0
            << " ASO=" << Incoming.Pot.ASO << std::endl;
  std::cout << "            RC0=" << Incoming.Pot.RC0 << std::endl;
  std::cout << std::endl;

  std::cout << "        OPTICAL MODEL SCATTERING FOR THE OUTGOING CHANNEL" << std::endl;
  std::cout << " POTENTIAL: V=" << Outgoing.Pot.V << " R0=" << Outgoing.Pot.R0
            << " A=" << Outgoing.Pot.A << std::endl;
  std::cout << "            VI=" << Outgoing.Pot.VI << " RI0=" << Outgoing.Pot.RI0
            << " AI=" << Outgoing.Pot.AI << std::endl;
  std::cout << "            VSO=" << Outgoing.Pot.VSO << " RSO0=" << Outgoing.Pot.RSO0
            << " ASO=" << Outgoing.Pot.ASO << std::endl;
  std::cout << "            RC0=" << Outgoing.Pot.RC0 << std::endl;
  std::cout
      << "================================================================="
      << std::endl;
}

// ComputeDCS: run BETCAL/AMPCAL/XSECTN on an externally-set TransferSMatrix
// Call after SetupChannels() and SetTransferSMatrix()
void DWBA::ComputeDCS() {
  XSectn();
}
