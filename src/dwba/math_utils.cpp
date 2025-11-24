#include "math_utils.h"
#include "JSymbols.h"
#include <cmath>
#include <iostream>

// Wrapper functions for JSymbols.h
double ClebschGordan(double J1, double m1, double J2, double m2, double J,
                     double m) {
  return CGcoeff(J, m, J1, m1, J2, m2);
}

double ThreeJ(double J1, double m1, double J2, double m2, double J3,
              double m3) {
  return ThreeJSymbol(J1, m1, J2, m2, J3, m3);
}

double SixJ(double J1, double J2, double J3, double J4, double J5, double J6) {
  return SixJSymbol(J1, J2, J3, J4, J5, J6);
}

double NineJ(double J1, double J2, double J3, double J4, double J5, double J6,
             double J7, double J8, double J9) {
  return NineJSymbol(J1, J2, J3, J4, J5, J6, J7, J8, J9);
}

// Factorial function (exposed from JSymbols.h if needed, or re-implemented)
double Factorial(double n) { return factorial(n); }
// ...
void GaussLegendre(int n, double x1, double x2, std::vector<double> &x,
                   std::vector<double> &w) {
  x.resize(n);
  w.resize(n);
  double m = (n + 1) / 2.0;
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);
  for (int i = 0; i < m; i++) {
    double z = std::cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5));
    double p1, p2, p3, pp;
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 0; j < n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z = z - p1 / pp;
    } while (std::abs(p1 / pp) > 1e-14);
    x[i] = xm - xl * z;
    x[n - 1 - i] = xm + xl * z;
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
}
