#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <string>

// Global variables for potential parameters
extern double v, r0, a;
extern double vi, ri0, ai;
extern double vsi, rsi0, asi;
extern double vso, rso0, aso;
extern double vsoi, rsoi0, asoi, rc0;

void PrintPotential();
std::string potentialRef(std::string name);
void CallPotential(std::string name, int A, int Z, double E, int Zproj);

#endif
