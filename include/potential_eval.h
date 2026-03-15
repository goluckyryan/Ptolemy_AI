#ifndef POTENTIAL_EVAL_H
#define POTENTIAL_EVAL_H

#include "dwba.h"

// Function to evaluate potential at a given radius
// A_target: mass number of target, A_projectile: mass number of projectile
// Ptolemy R0DEFAULT convention:
//   R0MASS = A_heavy^(1/3) + A_light^(1/3)  when A_light > 2.5
//   R0MASS = A_heavy^(1/3)                   when A_light <= 2.5
void EvaluatePotential(double r, const ChannelPotential& pot, double& V_real, double& V_imag, 
                       double& V_so_real, double& V_so_imag, double& V_coulomb, 
                       double Z1, double Z2, double A_target, double A_projectile = 1.0);

#endif
