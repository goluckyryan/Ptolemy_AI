// test_rcwfn_accuracy.cpp
// Tests Rcwfn accuracy for L=11..36, rho=89.39, eta=5.52341
// Checks Wronskian W(F,G) = F*G' - G*F' (should = 1 in rho units, i.e. F*G'-G*F' = 1)
// Also checks two-point determinant det = f2*g1 - f1*g2 (should be nonzero/stable)

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "include/rcwfn.h"

int main() {
    // Physical parameters for 148Sm(a,a) at 50 MeV
    double k = 3.01198;   // fm^-1
    double eta = 5.52341;
    double h = 0.1;       // fm grid spacing
    int N = 300;
    int n_match = N - 3;  // = 297

    double R1 = (n_match - 1) * h;  // r at n_match-1
    double R2 = n_match * h;         // r at n_match

    double rho1 = k * R1;
    double rho2 = k * R2;

    std::cout << "k = " << k << " fm^-1, eta = " << eta << "\n";
    std::cout << "n_match = " << n_match << ", h = " << h << " fm\n";
    std::cout << "R1 = " << R1 << " fm, rho1 = " << rho1 << "\n";
    std::cout << "R2 = " << R2 << " fm, rho2 = " << rho2 << "\n\n";

    std::cout << std::setw(5)  << "L"
              << std::setw(14) << "F(rho2)"
              << std::setw(14) << "G(rho2)"
              << std::setw(14) << "FP(rho2)"
              << std::setw(14) << "GP(rho2)"
              << std::setw(16) << "Wronsk(rho2)"
              << std::setw(16) << "det(f2g1-f1g2)"
              << "\n";
    std::cout << std::string(105, '-') << "\n";

    for (int L = 11; L <= 36; ++L) {
        std::vector<double> FC1(L+2), FCP1(L+2), GC1(L+2), GCP1(L+2);
        std::vector<double> FC2(L+2), FCP2(L+2), GC2(L+2), GCP2(L+2);

        int ret1 = Rcwfn(rho1, eta, L, L, FC1, FCP1, GC1, GCP1);
        int ret2 = Rcwfn(rho2, eta, L, L, FC2, FCP2, GC2, GCP2);

        double f1 = FC1[L], g1 = GC1[L], fp1 = FCP1[L], gp1 = GCP1[L];
        double f2 = FC2[L], g2 = GC2[L], fp2 = FCP2[L], gp2 = GCP2[L];

        // Wronskian in rho units: F*G' - G*F' should = 1
        double wronsk1 = f1 * gp1 - g1 * fp1;
        double wronsk2 = f2 * gp2 - g2 * fp2;

        // Two-point determinant det = f2*g1 - f1*g2
        double det = f2 * g1 - f1 * g2;

        std::cout << std::setw(5) << L
                  << std::scientific << std::setprecision(5)
                  << std::setw(14) << f2
                  << std::setw(14) << g2
                  << std::setw(14) << fp2
                  << std::setw(14) << gp2
                  << std::setw(16) << wronsk2
                  << std::setw(16) << det
                  << "  ret=" << ret2
                  << "\n";
    }

    std::cout << "\nNote: Wronskian should be 1.0 in rho-units (F*G'-G*F'=1)\n";
    std::cout << "Note: det = f2*g1 - f1*g2 is the two-point matching denominator\n";

    return 0;
}
