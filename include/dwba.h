#ifndef DWBA_H
#define DWBA_H

#include "Isotope.h"
#include "Potentials.h"
#include <complex>
#include <string>
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

  // Wave functions
  std::vector<std::complex<double>> WaveFunction; // Radial wave function
  std::vector<double> RGrid;

  // Potentials on grid
  std::vector<double> V_real;
  std::vector<double> V_imag;
  std::vector<double> V_so_real;
  std::vector<double> V_so_imag;
  std::vector<double> V_coulomb;

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

  // Main calculation
  void Calculate();
  void PrintParameters();

private:
  // Internal state
  Channel Incoming, Outgoing;
  BoundState TargetBS, ProjectileBS;
  double AngleMin, AngleMax, AngleStep;
  double Ex;

  // Integration Grid
  std::vector<double> ThetaGrid;
  std::vector<double> ThetaWeights;

  // Transfer S-matrix elements (indexed by La)
  // We might need a more complex structure if we want to store by (La, Lb).
  // For now, let's store a struct or map?
  // Or just a vector of structs {La, Lb, S}.
  struct TransferSMatrixElement {
    int La;
    int Lb;
    std::complex<double> S;
  };
  std::vector<TransferSMatrixElement> TransferSMatrix;

  // Methods corresponding to Fortran subroutines
  void WavSet(Channel &ch);
  void WavPot(Channel &ch);
  void WavElj(Channel &ch, int L, int Jp);
  void GrdSet();
  void InelDc();
  void XSectn();

  // Helpers
  void CalculateKinematics();
  void CalculateBoundState(Channel &ch, int n, int l, double j,
                           double bindingEnergy);
  void Integrate(const Channel &ch, int L,
                 std::vector<std::complex<double>> &wf);
};

#endif
