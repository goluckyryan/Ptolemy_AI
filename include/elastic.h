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
    void SetWynn(bool enable) { useWynn_ = enable; }  // Enable Wynn epsilon series acceleration

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
    double Ecm() const { return Ecm_; }
    double S()   const { return S_;   }

    // S-matrix accessor: twoJ = 2*J (integer)
    // idx = J - (L - S) = (twoJ - (2L - 2S)) / 2 = (twoJ - 2L + 2S) / 2... 
    // Smat_[L] index: J=L-S -> idx=0, J=L+S -> idx=2S
    std::complex<double> GetSMatrix(int L, int twoJ) const {
        if (L < 0 || L >= (int)Smat_.size()) return {0,0};
        int twoS = (int)std::round(2*S_);
        int idx = (twoJ - (2*L - twoS)) / 2;
        if (idx < 0 || idx >= (int)Smat_[L].size()) return {0,0};
        return Smat_[L][idx];
    }

    // SetGrid overload: (h, rmax) -> N = rmax/h
    void SetGrid(double h, double rmax) { h_ = h; N_ = (int)(rmax/h + 0.5); }

    // Wavefunction accessor: returns normalized u(r) array for (L, twoJ=2J).
    // Grid: r[i] = i * h, i = 0..N. Call after Solve().
    const std::vector<std::complex<double>>& GetWavefunction(int L, int twoJ) const {
        static const std::vector<std::complex<double>> empty;
        if (L < 0 || L >= (int)wfunc_.size()) return empty;
        int twoS = (int)std::round(2*S_);
        int idx = (twoJ - (2*L - twoS)) / 2;
        if (idx < 0 || idx >= (int)wfunc_[L].size()) return empty;
        return wfunc_[L][idx];
    }
    double GetH()   const { return h_; }
    int    GetN()   const { return N_; }
    double GetK()   const { return k_; }
    double GetEta() const { return eta_; }
    double GetCoulombPhase(int L) const {
        if (L < 0 || L >= (int)CoulPhase_.size()) return 0.0;
        return CoulPhase_[L];
    }

    // SetKinematics: bypass SetSystem/CalcKinematics (for external k/eta/mu)
    void SetProjectile(int Ap, int Zp) { Ap_ = Ap; Zp_ = Zp;
        S_ = (Ap==1)?0.5:(Ap==2)?1.0:0.0; }
    void SetTarget(int At, int Zt) { At_ = At; Zt_ = Zt; }
    void SetKinematics(double k, double eta, double mu_amu) {
        k_ = k; eta_ = eta;
        mu_ = mu_amu * 931.494061;  // store in MeV/c^2
        Ecm_ = k_ * k_ * 197.32697 * 197.32697 / (2.0 * mu_);
        f_conv_ = 2.0 * mu_ / (197.32697 * 197.32697);  // AFAC = 2mu/hbar^2
    }
    void Solve() { CalcScatteringMatrix(); }

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

    // Options
    bool   useWynn_ = false;  // Enable Wynn epsilon series acceleration for nuclear amp

    // S-matrix: Smat_[L] is a vector of 2S+1 entries for J = L-S, ..., L+S
    // For spin-0: 1 entry.  spin-1/2: 2 entries (J=L+1/2, J=L-1/2).  spin-1: 3 entries.
    // Index: idx = J - (L - S)  i.e. J=L-S → idx=0, J=L → idx=S, J=L+S → idx=2S
    std::vector<std::vector<std::complex<double>>> Smat_;
    std::vector<double> CoulPhase_;

    // Wavefunctions: wfunc_[L][idx] = normalized u(r) on the Numerov grid
    // Same indexing as Smat_. Stored after CalcScatteringMatrix().
    std::vector<std::vector<std::vector<std::complex<double>>>> wfunc_;

    // Internal helpers
    double NuclearMass(int A, int Z) const;
    double CoulombPhase(int L) const;
    std::complex<double> RunNumerov(int L, double LS_val,
        const std::vector<double>& Vr, const std::vector<double>& Wi,
        const std::vector<double>& Vc,
        const std::vector<double>& VsoRe, const std::vector<double>& VsoIm,
        const std::vector<double>& FC1, const std::vector<double>& GC1,
        const std::vector<double>& FC2, const std::vector<double>& GC2,
        std::vector<std::complex<double>>* wf_out = nullptr) const;

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
    // Wynn epsilon algorithm: accelerate convergence of partial sum sequence
    // Returns best estimate of the limit of the series {S_0, S_1, S_2, ...}
    static std::complex<double> WynnEpsilon(const std::vector<std::complex<double>>& partial_sums);
};
