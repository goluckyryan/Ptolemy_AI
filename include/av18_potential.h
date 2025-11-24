#ifndef AV18_POTENTIAL_H
#define AV18_POTENTIAL_H

void CalculateAV18NN(double r);
void CalculateAV18EM(double r);
void CalculateAV18PW(double r, double L, double S, double J, double T, double T1z, double T2z);

double GetAV18NN(unsigned short index);
double GetAV18EM(unsigned short index);
double* GetAV18PW();

#endif
