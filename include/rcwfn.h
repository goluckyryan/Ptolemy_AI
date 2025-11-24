#ifndef RCWFN_H
#define RCWFN_H

#include <complex>
#include <vector>

/**
 * @brief Computes regular and irregular Coulomb wavefunctions.
 *
 * Translated from Ptolemy's RCWFN subroutine.
 *
 * @param RHO   The value of the radial coordinate (kr).
 * @param ETA   The value of the charge parameter (Sommerfeld parameter).
 * @param MINL  The minimum value of L.
 * @param MAXL  The maximum value of L.
 * @param FC    Output array for Regular Coulomb functions F_L (size MAXL+1).
 * @param FCP   Output array for derivatives F'_L (size MAXL+1).
 * @param GC    Output array for Irregular Coulomb functions G_L (size MAXL+1).
 * @param GCP   Output array for derivatives G'_L (size MAXL+1).
 * @param ACCUR Desired accuracy (e.g., 1e-14).
 * @return int  Return code (0 = success).
 */
int Rcwfn(double RHO, double ETA, int MINL, int MAXL, std::vector<double> &FC,
          std::vector<double> &FCP, std::vector<double> &GC,
          std::vector<double> &GCP, double ACCUR = 1e-14);

#endif // RCWFN_H
