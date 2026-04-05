#ifndef COULIN_H
#define COULIN_H

#include <vector>

// RCASYM: Asymptotic series in 1/rho for Coulomb functions
// Ported from Ptolemy Fortran RCASYM subroutine
// Returns 0 on success, negative on error
int Rcasym(int L, double eta, double rho, double sigL,
           double &zeta, double &phi, double &dzeta,
           double &F, double &Fp, double &G, double &Gp,
           std::vector<double> &Z_coeffs,
           double eps, int nmax, int &ntz);

// CLINTS: Coulomb integral for single (LI,LF) pair from RLOWER to infinity
// Computes ∫(RLOWER,∞) X(LF,etaF,kF*r) * Y(LI,etaI,kI*r) / r^N dr
// where X,Y are regular (F) or irregular (G) Coulomb functions
// Returns 0 on success
int Clints(double rlower, double etaI, double etaF, double kI, double kF,
           double sigI, double sigF, double accura, double &asmult,
           double &ffint, double &fgint, double &gfint, double &ggint,
           int N, int LI, int LF, int nterms, int npts);

// COULIN result storage
// FF(id,il) stored as flat vector, column-major matching Fortran
// Access: FF[id + il*ldldim]  (0-based id and il)
struct CoulinResult {
    std::vector<double> FF, FG, GF, GG;
    int ldldim;
    int nils;  // = 2*(LMAX-LMIN)+1
    
    double getFF(int id, int il) const { return FF[id + il*ldldim]; }
    double getFG(int id, int il) const { return FG[id + il*ldldim]; }
    double getGF(int id, int il) const { return GF[id + il*ldldim]; }
    double getGG(int id, int il) const { return GG[id + il*ldldim]; }
};

// COULIN: Coulomb integrals by L-recursion for all (LIN,LOUT) pairs
// Computes ∫(R,∞) F_out(r)*F_in(r)/r^N dr (and FG,GF,GG variants)
// N = power of r (typically LX+1)
// maxdel = max |LIN-LOUT| (typically LX)
// lmin,lmax = range of (LIN+LOUT)/2
// R = lower limit (SUMMAX for tail, 0 for pure Coulomb)
// allSw = true → compute all 4 types; false → FF only
// Returns 0 on success
int Coulin(int N, int maxdel, int lmin, int lmax,
           double etaOut, double kOut, const std::vector<double> &sigOut,
           double etaIn, double kIn, const std::vector<double> &sigIn,
           double R, bool allSw,
           CoulinResult &result,
           double accura, double asmult, int nterms, int npts);

#endif // COULIN_H
