// chi_test.cpp — standalone distorted wave comparison
//
// Prints chi_a (d+16O, Elab=20 MeV) and chi_b (p+17O, Elab_p=20.848 MeV)
// radial wavefunctions u_L(r) = r*chi_L(r) for L=0..12 at selected radii.
//
// Also prints S-matrix |S_L| and arg(S_L) for each L and J sub-channel.
// Compare these against Ptolemy's S-matrix (from print=1 output) and
// against the elastic DCS (already validated to <0.7% for both channels).
//
// Ptolemy normalization: u_L(r) is Numerov-integrated, then normalized so
// that at large r it matches F_L + Re[S_L]*G_L + i*Im[S_L]*G_L.
// Our ElasticSolver uses the same convention.

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include "include/elastic.h"

static void printChannel(const char* label,
                         int Ap, int Zp, int At, int Zt,
                         double Elab_MeV,
                         double V,  double R0,  double A0,
                         double VI, double RI0, double AI0,
                         double VSO, double RSO0, double ASO,
                         double RC0)
{
    std::cout << "\n=== " << label << " ===" << std::endl;

    ElasticSolver e;
    e.SetSystem(Ap, Zp, At, Zt, Elab_MeV);
    e.SetGrid(0.1, 200);  // h=0.1 fm, 200 steps → rmax=20 fm

    // Real central WS
    e.AddVolumeWS({-V,   0.0}, R0, A0);
    // Imaginary volume WS (separate geometry)
    e.AddVolumeWS({ 0.0, -VI}, RI0, AI0);
    // Spin-orbit
    e.AddSpinOrbit({-VSO, 0.0}, RSO0, ASO);
    // Coulomb
    e.AddCoulomb(RC0);

    e.CalcKinematics();
    e.PrintKinematics();
    e.CalcScatteringMatrix();

    double h = e.GetH();
    int twoS = (Ap == 2) ? 2 : 1;  // deuteron S=1, proton S=1/2

    // S-matrix table
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n  L   J      |S_L|     arg(S)/deg\n";
    for (int L = 0; L <= 12; L++) {
        for (int twoJ = std::abs(2*L - twoS); twoJ <= 2*L + twoS; twoJ += 2) {
            auto S = e.GetSMatrix(L, twoJ);
            std::cout << "  " << std::setw(2) << L
                      << " " << std::setw(3) << twoJ << "/2"
                      << "  " << std::setw(9) << std::abs(S)
                      << "  " << std::setw(10) << std::arg(S)*180.0/M_PI
                      << "\n";
        }
    }

    // Wavefunction table: print u_L(r) at key radii for L=0,1,2,3,5,7
    std::cout << "\n  Wavefunction u_L(r) = r*chi_L(r) at selected r (fm):\n";
    std::cout << "  L   J    ";
    for (double r : {1.0, 2.0, 3.0, 5.0, 7.0, 10.0})
        std::cout << std::setw(22) << ("r=" + std::to_string((int)r) + "fm");
    std::cout << "\n";

    for (int L : {0, 1, 2, 3, 5, 7}) {
        // Use J = L + S (dominant sub-channel)
        int twoJ_hi = 2*L + twoS;
        int twoJ_lo = std::abs(2*L - twoS);
        for (int twoJ : {twoJ_lo, twoJ_hi}) {
            if (twoJ == twoJ_lo && twoJ == twoJ_hi) { if (twoJ != twoJ_lo) continue; }
            const auto& wf = e.GetWavefunction(L, twoJ);
            std::cout << "  " << std::setw(2) << L
                      << " " << std::setw(3) << twoJ << "/2  ";
            for (double r : {1.0, 2.0, 3.0, 5.0, 7.0, 10.0}) {
                int idx = (int)std::round(r / h);
                std::complex<double> u = (idx < (int)wf.size()) ? wf[idx] : 0;
                std::cout << std::setw(10) << u.real()
                          << "+" << std::setw(10) << u.imag() << "i  ";
            }
            std::cout << "\n";
            if (twoJ_lo == twoJ_hi) break;  // only one J for L=0
        }
    }
}

int main() {
    // Channel a: d + 16O at Elab=20 MeV
    printChannel("chi_a: d + 16O, Elab=20 MeV",
        2, 1, 16, 8, 20.0,
        88.9546, 1.1489, 0.7508,
        2.3480,  1.3446, 0.6030,
        3.5570,  0.9720, 1.0110,
        1.3);

    // Channel b: p + 17O at Elab_p=20.848 MeV (lab energy of outgoing proton)
    printChannel("chi_b: p + 17O, Elab=20.848 MeV",
        1, 1, 17, 8, 20.848,
        49.5434, 1.1462, 0.6753,
        2.0611,  1.1462, 0.6753,
        5.2956,  0.9338, 0.5900,
        1.3);

    std::cout << "\nDone." << std::endl;
    return 0;
}
