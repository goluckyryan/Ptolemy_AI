#ifndef POTENTIAL_EVAL_H
#define POTENTIAL_EVAL_H

#include "dwba.h"

// Function to evaluate potential at a given radius
void EvaluatePotential(double r, const ChannelPotential& pot, double& V_real, double& V_imag, 
                       double& V_so_real, double& V_so_imag, double& V_coulomb, 
                       double Z1, double Z2, double mu);

#endif
