#include "av18_potential.h"
#include "av18.h"
#include <iostream>

// Global instance of AV18 class
AV18 av18_instance;

void CalculateAV18NN(double r) {
    av18_instance.CalNN(r);
}

void CalculateAV18EM(double r) {
    av18_instance.CalEM(r);
}

void CalculateAV18PW(double r, double L, double S, double J, double T, double T1z, double T2z) {
    av18_instance.CalPW(r, L, S, J, T, T1z, T2z);
}

double GetAV18NN(unsigned short index) {
    return av18_instance.GetNN(index);
}

double GetAV18EM(unsigned short index) {
    return av18_instance.GetEM(index);
}

double* GetAV18PW() {
    return av18_instance.GetPW();
}
