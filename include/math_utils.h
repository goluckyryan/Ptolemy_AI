#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <vector>

double ClebschGordan(double J1, double m1, double J2, double m2, double J,
                     double m);
double ThreeJ(double J1, double m1, double J2, double m2, double J3, double m3);
double SixJ(double J1, double J2, double J3, double J4, double J5, double J6);
double NineJ(double J1, double J2, double J3, double J4, double J5, double J6,
             double J7, double J8, double J9);
double Factorial(double n);
void GaussLegendre(int n, double x1, double x2, std::vector<double> &x,
                   std::vector<double> &w);

#endif
