// elastic.h — Unified elastic scattering solver
// Handles spin-0 and spin-1/2 projectiles, any combination of potentials.
// Modeled after Raphael's DistortedWave + solveSE Python classes.
#pragma once
#include <array>
#include <vector>
#include <complex>
#include <string>
#include <cmath>

// Potential type IDs (matching Raphael's pot.id)
enum PotType { kVolumeWS = 0, kSurfaceWS = 1, kSpinOrbit = 2, kCoulomb = 3 };

struct PotEntry {
    PotType   type;
    std::complex<double> V;   // depth in MeV (complex: real+imag parts)
    double    r0;             // radius parameter (fm/A^{1/3})
    double    a0;             // diffuseness (fm)
};

class ElasticSolver {
public:
    ElasticSolver();

    // --- System setup ---
    void SetSystem(int Ap, int Zp, int At, int Zt, double Elab_MeV);
    void SetGrid(double h = 0.05, int N = 601);
    void SetLmax(int Lmax);   // -1 = auto

    // --- Add potentials (Raphael-compatible interface) ---
    // V is a complex depth in MeV: real = real part, imag = imaginary part
    void AddVolumeWS (std::complex<double> V,   double r0, double a0);
    void AddSurfaceWS(std::complex<double> V,   double r0, double a0);
    void AddSpinOrbit(std::complex<double> Vso, double r0, double a0);
    void AddCoulomb  (double rc0);
    void ClearPotentials();

    // --- Compute ---
    void CalcKinematics();          // must call before CalcScatteringMatrix
    void CalcScatteringMatrix();

    // --- Output ---
    // Returns DCS in mb/sr (unpolarized)
    double DCSUnpolarized(double theta_deg) const;

    void SaveSMatrix(const std::string& filename) const;
    void SaveDCS(const std::string& filename,
                 double th_min = 1.0, double th_max = 180.0, double dth = 1.0) const;
    void PrintKinematics() const;
    void PrintSMatrix(int Lmax_print = -1) const;

    // Accessors
    double k()   const { return k_;   }
    double eta() const { return eta_; }
    double mu()  const { return mu_;  } // MeV/c^2

private:
    // System
    int Ap_, Zp_, At_, Zt_;
    double Elab_;
    double S_;    // projectile spin (0.0 for alpha, 0.5 for proton)

    // Grid
    double h_;
    int    N_;
    int    Lmax_;

    // Kinematics
    double k_, eta_, mu_, Ecm_, f_conv_;

    // Potentials
    std::vector<PotEntry> pots_;
    double rc0_;           // Coulomb radius parameter
    bool   hasCoulomb_;

    // S-matrix: Smat_[L] is a vector of 2S+1 entries for J = L-S, ..., L+S
    // For spin-0: 1 entry.  spin-1/2: 2 entries (J=L+1/2, J=L-1/2).  spin-1: 3 entries.
    // Index: idx = J - (L - S)  i.e. J=L-S → idx=0, J=L → idx=S, J=L+S → idx=2S
    std::vector<std::vector<std::complex<double>>> Smat_;
    std::vector<double> CoulPhase_;

    // Internal helpers
    double NuclearMass(int A, int Z) const;
    double CoulombPhase(int L) const;
    std::complex<double> RunNumerov(int L, double LS_val,
        const std::vector<double>& Vr, const std::vector<double>& Wi,
        const std::vector<double>& Vc,
        const std::vector<double>& VsoRe, const std::vector<double>& VsoIm,
        const std::vector<double>& FC1, const std::vector<double>& GC1,
        const std::vector<double>& FC2, const std::vector<double>& GC2) const;

    static double LegendreP(int L, double x);
    static double LegendrePprime(int L, double x);

    // DCS helpers — general spin formula (GMatrix / Clebsch-Gordan)
    std::complex<double> CoulombAmp(double theta_deg) const;
    // Returns associated Legendre P_l^m(cos θ)  (unnormalized, Condon-Shortley)
    static double AssocLegendreP(int l, int m, double x);
    // Nuclear scattering amplitude f_N(v, v0, θ) for spin projection v0→v
    std::complex<double> NuclearAmp(double v, double v0, double theta_deg) const;
    // Clebsch-Gordan coefficient (thin wrapper around math_utils)
    static double CG(double j1, double m1, double j2, double m2, double J, double M);
};
