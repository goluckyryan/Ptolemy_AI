#ifndef DWBA_H
#define DWBA_H

#include "Isotope.h"
#include "Potentials.h"
#include <complex>
#include <string>
#include <tuple>
#include <vector>

// Structure to hold potential parameters for a channel
struct ChannelPotential {
  double V, R0, A;
  double VI, RI0, AI;
  double VSO, RSO0, ASO;
  double VSOI, RSOI0, ASOI;
  double VSI, RSI0, ASI;
  double RC0;
  // Tensor potentials...
};

// Structure to hold channel kinematics and properties
struct Channel {
  Isotope Projectile, Target, Ejectile, Recoil;
  double Elab, Ecm, Q;
  double k;   // Wave number
  double eta; // Coulomb parameter
  double mu;  // Reduced mass
  ChannelPotential Pot;
  int JSPS = 1;  // 2*spin of projectile (deuteron=2, proton=1)

  // Wave functions
  std::vector<std::complex<double>> WaveFunction; // Radial wave function phi = u/r
  std::vector<double> RGrid;

  // Potentials on grid
  std::vector<double> V_real;
  std::vector<double> V_imag;
  std::vector<double> V_so_real;
  std::vector<double> V_so_imag;
  std::vector<double> V_coulomb;

  // Precomputed vertex product V(r)*phi(r) — used in BSPROD-style integrals.
  // For a plain WS bound state: filled by InelDc as V_real[i]*WaveFunction[i].real().
  // For the Reid deuteron projectile: loaded from reid-phi-v Block 1 (precomputed V_Reid*phi_S).
  // Empty means "not set" — InelDc will compute V_real*WaveFunction on the fly.
  std::vector<double> VPhiProduct; // V(r)*phi(r) on the calculation grid

  // Grid parameters
  double StepSize;
  double MaxR;
  int NSteps;

  // S-matrix elements
  std::vector<std::complex<double>> SMatrix;
  std::vector<double> CoulombPhase;
};

class DWBA {
public:
  DWBA();
  ~DWBA();

  // Configuration
  void SetReaction(const std::string &target, const std::string &projectile,
                   const std::string &ejectile, const std::string &recoil);
  void SetEnergy(double Elab);
  void SetProjectileWFFile(const std::string &filename, double grid_h, double spam);
  void SetExcitation(double Ex);
  void SetAngles(double min, double max, double step);
  void SetIncomingPotential(const ChannelPotential &pot);
  void SetOutgoingPotential(const ChannelPotential &pot);

  struct BoundState {
    int n;
    int l;
    double j;
    double BindingEnergy;
    // Potential parameters for bound state (usually same as core-particle
    // potential)
    ChannelPotential Pot;
  };

  void SetTargetBoundState(int n, int l, double j, double bindingEnergy,
                           const ChannelPotential &pot);
  void SetProjectileBoundState(int n, int l, double j, double bindingEnergy,
                               const ChannelPotential &pot);
  void SetProjectileWFFromFile(const std::vector<std::pair<double,double>>& wf_data,
                               double h_ext, double spam);

  // Nuclear spins (needed for generic XSectn / ATERM)
  void SetTargetSpin(double J);       // J of target nucleus A (e.g. 33Si: 1.5, 16O: 0)
  void SetResidualSpin(double J);     // J of residual nucleus B (e.g. 34Si: 0, 17O: 2.5)

  // Main calculation
  void Calculate();
  void CalculateZR();  // Zero-Range DWBA
  void PrintParameters();

  friend class DWBATest;  // for wftest_main.cpp validation
private:
  // Internal state
  Channel Incoming, Outgoing;
  BoundState TargetBS, ProjectileBS;
  double AngleMin, AngleMax, AngleStep;
  double Ex;

  // Nuclear spins (set by parser or defaulted in Calculate/CalculateZR)
  double SpinTarget   = -1.0;   // J of target nucleus A (-1 = not set, use heuristic)
  double SpinResidual = -1.0;   // J of residual nucleus B (-1 = not set, use heuristic)

  // Integration Grid
  std::vector<double> ThetaGrid;
  std::vector<double> ThetaWeights;

  // Transfer S-matrix elements (indexed by La)
  struct TransferSMatrixElement {
    int Lx;
    int Li;
    int Lo;
    int JPI;  // 2*J_incoming (doubled)
    int JPO;  // 2*J_outgoing (doubled)
    std::complex<double> S;
  };
  std::vector<TransferSMatrixElement> TransferSMatrix;

  // Projectile WF table (optional: loaded from Ptolemy-output file)
  std::vector<double> ProjectileWFTable;
  double ProjectileWFGridH = 0.0;
  double ProjectileWFSpam  = 1.0;
  bool   ProjectileWFLoaded = false;

  // Methods corresponding to Fortran subroutines
  void WavSet(Channel &ch);
  void WavPot(Channel &ch);
  void WavElj(Channel &ch, int L, int Jp);
  void GrdSet();
  void InelDc();
  void InelDcZR();  // Zero-Range transfer integral
  std::vector<std::tuple<int,int,double>> ComputeA12Terms(int Li, int Lo, int Lx, int lT, int lP);
  double EvalA12(const std::vector<std::tuple<int,int,double>>& A12_terms,
                 double phi_T_angle, double phi_ab);
  void XSectn();

  // Helpers
  void CalculateKinematics();
  void CalculateBoundState(Channel &ch, int n, int l, double j,
                           double bindingEnergy);
  // Load a tabulated deuteron bound-state WF from a Ptolemy linkule file.
  // filename = "reid-phi-v" or "av18-phi-v" (in dataPath directory).
  // Populates ch.WaveFunction (phi_S = u/r, normalized) and
  // ch.VPhiProduct (V_np(r)*phi_S(r), same normalization — signed, not abs).
  bool LoadDeuteronWavefunction(Channel &ch, const std::string &dataPath,
                                const std::string &filename = "av18-phi-v");

public:
  // Public wrapper for standalone bound state tests
  void CalcBoundState(Channel &ch, int n, int l, double j, double bindingEnergy) {
    CalculateBoundState(ch, n, l, j, bindingEnergy);
  }
  // Public accessor for transfer S-matrix (before 9-J coupling)
  const std::vector<TransferSMatrixElement>& GetTransferSMatrix() const { return TransferSMatrix; }
  void Integrate(const Channel &ch, int L,
                 std::vector<std::complex<double>> &wf);
};

#endif
