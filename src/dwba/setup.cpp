// setup.cpp — DWBA constructor, destructor, setters, Calculate, PrintParameters
// Extracted from dwba.cpp (lines 1–180, 249–341)

#include "dwba.h"
#include "ptolemy_mass_table.h"
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
  // INELOCA parameterset: STEPSPER=15; default (transfer/dpsb): STEPSPER=8
  // STEPSPER from Fortran PARAM subroutine:
  // INELOCA*: STPSPR=15 (source.f:28049 RGRIDS(3,N)=15)  
  // No PARAMETERSET, inelastic: STPSPR=10 (Fortran default for STPSPR, source.f:14132)
  // Transfer/elastic (dpsb): STPSPR=8 (from GRIDEL/RGRIDS)
  // STEPSPER from PARAMETERSET (Fortran RGRIDS(3,N) = STPSPR)
  // alpha3=25, DPSB=8, INELOCA=15, default=8
  int STEPSPER;
  if (ParameterSet.find("INELOCA") != std::string::npos) STEPSPER = 15;
  else if (BELx > 0.0) STEPSPER = 10;  // inelastic, no PARAMETERSET
  else STEPSPER = (int)GrdSTEPSPR;     // from PARAMETERSET (alpha3=25, DPSB=8)
  Incoming.StepSize = std::min(2.0 * M_PI / Incoming.k, 1.0) / STEPSPER;
  Outgoing.StepSize = std::min(2.0 * M_PI / Outgoing.k, 1.0) / STEPSPER;

  // ── WAVSET NSTEP: exact Fortran algorithm (source.f lines 32500-32502, 37793-37805) ──
  //
  // Ptolemy SETCHN (line 32500-32502) sets ASYMPT for scattering channels:
  //   ASYMPT = ABS(SCTASY)    (for IBND == 0, i.e., scattering)
  // where SCTASY comes from the DPSB preset (default -20 fm, negative = allow L-adjust).
  //
  // The user keyword 'asymptopia=50' sets FLOAT(9) which goes to GRDSET's integration
  // range (SUMMAX), NOT to WAVSET's ASYMPT. These are separate concepts:
  //   - WAVSET ASYMPT: controls Numerov integration range for chi (S-matrix extraction)
  //   - AsymptopiaSet: controls GRDSET radial integral range (SUMMAX)
  //
  // Ptolemy WAVSET (line 37793-37805) then computes RMAX:
  //   RMAX = ASYMPT                                             (base)
  //   IF (SCTASY < 0)  RMAX = max(RMAX, turning_point(LMAX))   (L-extend)
  //   NSTEP = RMAX/STEPSZ + 0.5                                (snap to grid)
  //   RMAX = NSTEP * STEPSZ
  //
  // The chi extension beyond NSTEP (NSTP2S) is handled separately in wavelj.cpp
  // for the GRDSET integration that needs chi up to SUMMAX > RMAX.
  {
    // Base asymptopia for scattering: ABS(SCTASY)
    // SCTASY from DPSB preset = -20.0 (negative = allow L-adjust)
    double ASYMPT_base = std::abs(SctAsySet);  // = 20.0 fm
    int LMAX_eff = (LmaxSet >= 0) ? LmaxSet : 40;
    // For inelastic without PARAMETERSET: no L-adjustment (SCTASY > 0 in Fortran default)
    // Fortran INELOCA sets SCTASY=-(GRIDIN(5,N)) which is -20 (allows adjustment)
    // Without PARAMETERSET, no SCTASY override → use ASYMPT_base = 20, no adjust
    bool allow_L_adjust = (SctAsySet < 0) && (ParameterSet.find("INELOCA") != std::string::npos || BELx == 0.0);

    // Incoming channel
    // For INELASTIC (INELOCA): NSTEP = ASYMPT_base/h (Fortran SCTASY=20fm)
    //   Coulomb extension to SUMMAX is done inside wavelj.cpp (NSTP2S)
    // For TRANSFER (BELx==0, dpsb): NSTEP = max(ASYMPT_base, turning_point)/h
    //   because the S-matrix matching must be beyond the Coulomb turning point
    double RMAX_in = ASYMPT_base;
    bool is_inelastic_ineloca = (ParameterSet.find("INELOCA") != std::string::npos && BELx > 0.0);
    if (allow_L_adjust && !is_inelastic_ineloca) {
      // Transfer: include Coulomb turning point for correct S-matrix extraction
      double tp_in = (Incoming.eta + std::sqrt(Incoming.eta * Incoming.eta
                      + (double)LMAX_eff * (LMAX_eff + 1))) / Incoming.k;
      RMAX_in = std::max(RMAX_in, tp_in);
    }
    int nstep_in = static_cast<int>(RMAX_in / Incoming.StepSize + 0.5);
    Incoming.MaxR = nstep_in * Incoming.StepSize;
    Incoming.NSteps = nstep_in;

    // Outgoing channel
    double RMAX_out = ASYMPT_base;
    if (allow_L_adjust && !is_inelastic_ineloca) {
      double tp_out = (Outgoing.eta + std::sqrt(Outgoing.eta * Outgoing.eta
                       + (double)LMAX_eff * (LMAX_eff + 1))) / Outgoing.k;
      RMAX_out = std::max(RMAX_out, tp_out);
    }
    int nstep_out = static_cast<int>(RMAX_out / Outgoing.StepSize + 0.5);
    Outgoing.MaxR = nstep_out * Outgoing.StepSize;
    Outgoing.NSteps = nstep_out;

    std::cout << "[WAVSET] Incoming: ASYMPT_base=" << ASYMPT_base
              << " RMAX=" << Incoming.MaxR
              << " fm, NSTEP=" << nstep_in
              << ", h=" << Incoming.StepSize << " fm\n";
    std::cout << "[WAVSET] Outgoing: ASYMPT_base=" << ASYMPT_base
              << " RMAX=" << Outgoing.MaxR
              << " fm, NSTEP=" << nstep_out
              << ", h=" << Outgoing.StepSize << " fm\n";
  }

  WavSet(Incoming);
  WavSet(Outgoing);

  // ── Inelastic collective branch ──
  if (BELx > 0.0) {
    InelDcCollective();
    return;
  }

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
  const double AMU_MEV = 931.5016;   // Ptolemy AMUMEV
  const double HBARC_L = 197.32858;  // Ptolemy HBARC

  // Masses from Ptolemy AME2003 mass excesses (matching SETCHN)
  const double EMASS_S = 0.511;
  auto ptolemy_mass_s = [&](int Z, int A) -> double {
      double MX = PtolemyMass::MassExcess_MeV(Z, A);
      return (A + MX/AMU_MEV - Z*(EMASS_S/AMU_MEV)) * AMU_MEV;
  };
  double ma = ptolemy_mass_s((int)Incoming.Projectile.Z, Incoming.Projectile.A);
  double mb = ptolemy_mass_s((int)Outgoing.Projectile.Z, Outgoing.Projectile.A);
  double mA = ptolemy_mass_s((int)Incoming.Target.Z, Incoming.Target.A);
  int Zx_s = std::abs((int)Incoming.Projectile.Z - (int)Outgoing.Projectile.Z);
  int Ax_s = std::abs(Incoming.Projectile.A - Outgoing.Projectile.A);
  double mx = ptolemy_mass_s(Zx_s, Ax_s);
  // Fortran BOUND: AMP=incoming(a), AMT=transferred(x) → mu = ma*mx/(ma+mx)
  // NOT mb*mx/(mb+mx). For (p,d): mu = m_proton*m_neutron/(m_proton+m_neutron) ≈ 469 MeV.
  double mu_pbs = ma * mx / (ma + mx);
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
