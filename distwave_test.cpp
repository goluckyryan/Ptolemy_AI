// distwave_test.cpp — 16O(d,p)17O distorted waves
// Incoming: d + 16O @ 20 MeV    (Table 5.1)
// Outgoing: p + 17O @ 20.85 MeV (Table 5.2)
// Uses ElasticSolver — validated vs Ptolemy to <0.7%

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include "include/elastic.h"

static const double HBARC = 197.32697;
static const double AMU   = 931.494061;
static const double FINE  = 1.0/137.035999;

int main() {

    // ── Incoming channel: d + 16O @ 20 MeV ─────────────────────────
    double m_d   = 2.01410;
    double m_16O = 15.99491;
    double mu_in_amu = m_d * m_16O / (m_d + m_16O);
    double Elab_in   = 20.0;
    double Ecm_in    = Elab_in * m_16O / (m_d + m_16O);
    double mu_in_MeV = mu_in_amu * AMU;
    double k_in      = std::sqrt(2.0 * mu_in_MeV * Ecm_in) / HBARC;
    double eta_in    = 1.0 * 8.0 * FINE * mu_in_MeV * HBARC / (k_in * HBARC * HBARC);
    // eta = Z1*Z2*e^2/(hbar*v) = Z1*Z2/(137) * mu*c^2 / (hbar*k*c)
    eta_in = 1.0 * 8.0 / 137.035999 * mu_in_MeV / (HBARC * k_in);

    std::cout << "=== Incoming: d + 16O @ " << Elab_in << " MeV ===" << std::endl;
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  mu=" << mu_in_amu << " AMU  Ecm=" << Ecm_in
              << " MeV  k=" << k_in << " fm^-1  eta=" << eta_in << std::endl;

    ElasticSolver solI;
    solI.SetTarget(16, 8);
    solI.SetProjectile(2, 1);    // deuteron S=1
    solI.SetKinematics(k_in, eta_in, mu_in_amu);
    solI.SetLmax(20);
    solI.SetGrid(0.01, 30.0);

    // Table 5.1: d + 16O at 10 MeV/u
    solI.AddVolumeWS( {-88.955,  0.0},   1.149, 0.751);  // real V
    solI.AddVolumeWS( {0.0,     -2.348}, 1.345, 0.603);  // imag VI (separate r0/a0)
    solI.AddSurfaceWS({0.0,    -10.218}, 1.397, 0.687);  // surface
    solI.AddSpinOrbit({-3.557,   0.0},   0.972, 1.011);  // SO
    solI.AddCoulomb(1.303);
    solI.Solve();

    std::cout << "  L      J    Re(S)      Im(S)      |S|" << std::endl;
    std::ofstream of_in("/tmp/cpp_smat_inc.txt");
    of_in << "# L J Re Im |S|\n";
    for (int L = 0; L <= 15; L++) {
        for (int twoJ : {2*L-2, 2*L, 2*L+2}) {
            if (twoJ < 0) continue;
            double J = twoJ / 2.0;
            auto S = solI.GetSMatrix(L, twoJ);
            std::cout << "  " << std::setw(2) << L << "  " << std::setw(5) << J
                      << "  " << std::setw(9) << S.real()
                      << "  " << std::setw(9) << S.imag()
                      << "  " << std::abs(S) << std::endl;
            of_in << L << " " << J << " " << S.real()
                  << " " << S.imag() << " " << std::abs(S) << "\n";
        }
    }
    of_in.close();
    std::cout << "  S-matrix saved: /tmp/cpp_smat_inc.txt" << std::endl;

    // ── Outgoing channel: p + 17O @ 20.85 MeV ──────────────────────
    double m_p   = 1.00728;
    double m_17O = 16.99913;
    double Eout_lab  = 20.8505;   // from Raphael ReactionData
    double mu_out_amu = m_p * m_17O / (m_p + m_17O);
    double Ecm_out    = Eout_lab * m_17O / (m_p + m_17O);
    double mu_out_MeV = mu_out_amu * AMU;
    double k_out      = std::sqrt(2.0 * mu_out_MeV * Ecm_out) / HBARC;
    double eta_out    = 1.0 * 8.0 / 137.035999 * mu_out_MeV / (HBARC * k_out);

    std::cout << "\n=== Outgoing: p + 17O @ " << Eout_lab << " MeV ===" << std::endl;
    std::cout << "  mu=" << mu_out_amu << " AMU  Ecm=" << Ecm_out
              << " MeV  k=" << k_out << " fm^-1  eta=" << eta_out << std::endl;

    ElasticSolver solO;
    solO.SetTarget(17, 8);
    solO.SetProjectile(1, 1);    // proton S=1/2
    solO.SetKinematics(k_out, eta_out, mu_out_amu);
    solO.SetLmax(15);
    solO.SetGrid(0.01, 30.0);

    // Table 5.2: p + 17O at 10 MeV/u (V and VI same r0/a0 -> single complex call)
    solO.AddVolumeWS( {-49.544, -2.061}, 1.146, 0.675);
    solO.AddSurfaceWS({0.0,     -7.670}, 1.302, 0.528);
    solO.AddSpinOrbit({-5.296,  +0.106}, 0.934, 0.590);
    solO.AddCoulomb(1.419);
    solO.Solve();

    std::cout << "  L      J    Re(S)      Im(S)      |S|" << std::endl;
    std::ofstream of_out("/tmp/cpp_smat_out.txt");
    of_out << "# L J Re Im |S|\n";
    for (int L = 0; L <= 14; L++) {
        for (int twoJ : {2*L-1, 2*L+1}) {
            if (twoJ < 0) continue;
            double J = twoJ / 2.0;
            auto S = solO.GetSMatrix(L, twoJ);
            std::cout << "  " << std::setw(2) << L << "  " << std::setw(5) << J
                      << "  " << std::setw(9) << S.real()
                      << "  " << std::setw(9) << S.imag()
                      << "  " << std::abs(S) << std::endl;
            of_out << L << " " << J << " " << S.real()
                   << " " << S.imag() << " " << std::abs(S) << "\n";
        }
    }
    of_out.close();
    std::cout << "  S-matrix saved: /tmp/cpp_smat_out.txt" << std::endl;

    return 0;
}
