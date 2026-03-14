// grdset.cpp — DWBA::GrdSet() (integration grid setup)
// Extracted from ineldc.cpp

#include "dwba.h"
#include "math_utils.h"

void DWBA::GrdSet() {
  // NOTE: phi_ab integration now uses CUBMAP adaptive PHI0 in InelDc.
  // ThetaGrid is kept for compatibility but no longer used in the phi loop.
  // The phi integration uses NPPHI=10 CUBMAP points on [0, PHI0] per (ra,rb) pair.
  int NTheta = 10;   // retained for compatibility; not used in phi integration
  GaussLegendre(NTheta, -1.0, 1.0, ThetaGrid, ThetaWeights);
}
