// smat_o16dp.cpp — Compute elastic S-matrices for 16O(d,p)17O channels
// Channel a: d + 16O at Elab=20 MeV (An-Cai 2006 OMP)
// Channel b: p + 17O at Elab=23.21 MeV (Koning-Delaroche OMP)
//
// Build:
//   g++ -O2 -std=c++17 -Iinclude smat_o16dp.cpp src/elastic/elastic.cpp \
//       src/dwba/rcwfn.cpp src/dwba/math_utils.cpp -o smat_o16dp -lm

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "elastic.h"

int main() {

    // ====================================================================
    // Channel a: d + 16O at Elab = 20 MeV  (An-Cai 2006)
    // ====================================================================
    std::cout << "===== d + 16O at Elab=20 MeV (An-Cai 2006) =====\n";
    {
        ElasticSolver s;
        // Ap=2, Zp=1 (deuteron); At=16, Zt=8
        s.SetSystem(2, 1, 16, 8, 20.0);

        // Volume real WS: -V (attractive)
        s.AddVolumeWS ({-88.9546,  0.0000}, 1.1489, 0.7508);
        // Volume imaginary WS (different R0, a)
        s.AddVolumeWS ({  0.0000, -2.3480}, 1.3446, 0.6030);
        // Surface imaginary WS deriv
        s.AddSurfaceWS({  0.0000,-10.2180}, 1.3943, 0.6872);
        // Spin-orbit: physical VSO = 7.114/2 = 3.557 MeV (C++ uses correct L·S)
        // VSOI = 0 for d channel
        s.AddSpinOrbit({ -3.5570,  0.0000}, 0.9720, 1.0110);
        // Coulomb: RC_phys = RC0 * A_target^(1/3)
        s.AddCoulomb(1.3030);

        s.SetLmax(50);
        s.CalcKinematics();
        s.PrintKinematics();
        s.CalcScatteringMatrix();
        s.PrintSMatrix(12);

        // Save S-matrix: L  J  Re(S)  Im(S)  |S|
        {
            std::ofstream f("Raphael_AI/d16O_smat_cpp.txt");
            f << "# d + 16O elastic S-matrix at Elab=20 MeV (An-Cai 2006, C++)\n";
            f << "# L   J      Re(S)         Im(S)         |S|\n";
            f << std::fixed << std::setprecision(8);

            // spin-1 deuteron: J = L-1, L, L+1 for each L
            // twoJ values: 2(L-1), 2L, 2(L+1)
            int Lmax = 30;
            for (int L = 0; L <= Lmax; L++) {
                // J values: L-S to L+S, S=1
                int twoS = 2;
                for (int di = -twoS; di <= twoS; di += 2) {
                    int twoJ = 2*L + di;
                    if (twoJ < 0) continue;  // J must be >= 0
                    // J must be >= 0, and for L=0 only J=1 exists (L-1=-1 invalid)
                    auto Sval = s.GetSMatrix(L, twoJ);
                    double J = twoJ / 2.0;
                    double absSmat = std::abs(Sval);
                    f << std::setw(3) << L
                      << "  " << std::setw(5) << J
                      << "  " << std::setw(14) << Sval.real()
                      << "  " << std::setw(14) << Sval.imag()
                      << "  " << std::setw(14) << absSmat
                      << "\n";
                }
            }
            std::cout << "Saved Raphael_AI/d16O_smat_cpp.txt\n";
        }
    }

    // ====================================================================
    // Channel b: p + 17O at Elab = 23.21 MeV  (Koning-Delaroche OMP)
    // ====================================================================
    std::cout << "\n===== p + 17O at Elab=23.21 MeV (Koning-Delaroche) =====\n";
    {
        ElasticSolver s;
        // Ap=1, Zp=1 (proton); At=17, Zt=8
        s.SetSystem(1, 1, 17, 8, 23.21);

        // Volume WS: real and imaginary share same R0, a — single AddVolumeWS
        s.AddVolumeWS ({-49.5434, -2.0611}, 1.1462, 0.6753);
        // Surface imaginary WS deriv
        s.AddSurfaceWS({  0.0000, -7.6703}, 1.3016, 0.5275);
        // Spin-orbit (real + imaginary)
        // Note: proton spin S=1/2, C++ has correct L·S
        // Ptolemy VSO=5.2956 is already physical for spin-1/2 proton
        s.AddSpinOrbit({ -5.2956,  0.1059}, 0.9338, 0.5900);
        // Coulomb
        s.AddCoulomb(1.3030);

        s.SetLmax(30);
        s.CalcKinematics();
        s.PrintKinematics();
        s.CalcScatteringMatrix();
        s.PrintSMatrix(12);

        // Save S-matrix: L  J  Re(S)  Im(S)  |S|
        {
            std::ofstream f("Raphael_AI/p17O_smat_cpp.txt");
            f << "# p + 17O elastic S-matrix at Elab=23.21 MeV (KD OMP, C++)\n";
            f << "# L   J      Re(S)         Im(S)         |S|\n";
            f << std::fixed << std::setprecision(8);

            // spin-1/2 proton: J = L-1/2 or L+1/2
            int Lmax = 20;
            for (int L = 0; L <= Lmax; L++) {
                // twoJ = 2L-1 and 2L+1
                for (int twoJ : {2*L - 1, 2*L + 1}) {
                    if (twoJ < 1) continue;  // J must be >= 1/2
                    auto Sval = s.GetSMatrix(L, twoJ);
                    double J = twoJ / 2.0;
                    double absSmat = std::abs(Sval);
                    f << std::setw(3) << L
                      << "  " << std::setw(5) << J
                      << "  " << std::setw(14) << Sval.real()
                      << "  " << std::setw(14) << Sval.imag()
                      << "  " << std::setw(14) << absSmat
                      << "\n";
                }
            }
            std::cout << "Saved Raphael_AI/p17O_smat_cpp.txt\n";
        }
    }

    return 0;
}
