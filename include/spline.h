#pragma once
// spline.h — Ptolemy SPLNCB + INTRPC exact C++ port
// Source: fortlib.mor

#include <vector>

// Natural cubic spline coefficients (Ptolemy SPLNCB)
// Y(x) = Y[i] + del*B[i] + del^2*C[i] + del^3*D[i]  for X[i] <= x <= X[i+1]
// del = x - X[i]
// Boundary: natural (C=0 at both ends)
void Splncb(int N, const double* X, const double* Y, double* B, double* C, double* D);
void Splncb(int N, const std::vector<double>& X, const std::vector<double>& Y,
            std::vector<double>& B, std::vector<double>& C, std::vector<double>& D);

// Cubic spline evaluator (Ptolemy INTRPC)
// Evaluates spline from Splncb at NPTS arbitrary XS points → YS
void Intrpc(int NCUBIC,
            const double* XCUBES, const double* AS,
            const double* BS, const double* CS, const double* DS,
            int NPTS, const double* XS, double* YS);
void Intrpc(int NCUBIC, const std::vector<double>& XCUBES,
            const std::vector<double>& AS, const std::vector<double>& BS,
            const std::vector<double>& CS, const std::vector<double>& DS,
            int NPTS, const std::vector<double>& XS, std::vector<double>& YS);
