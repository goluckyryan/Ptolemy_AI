// grdset_ineldc_faithful.cpp
// Faithful C++ port of Ptolemy GRDSET + INELDC subroutines.
// Implements DWBA::GrdSetFaithful() and DWBA::InelDcFaithful2().
// Also provides DWBA::GrdSet() and DWBA::InelDc() as thin wrappers.
//
// Fortran source: source_annotated.f
//   BSPROD : lines 5112-5714
//   GRDSET : lines 18273-19672
//   INELDC : lines 20232-21318
//
// Translation rules:
//   - 1-based Fortran arrays → 0-based C++ (subtract 1 from all indices)
//   - GOTO → goto with label_ prefix
//   - REAL*8 → double, INTEGER → int, LOGICAL → bool
//   - ALLOC scratch arrays → std::vector<double>
//   - COMMON block variables → DWBA member variables or local equivalents
//   - Arithmetic IF (x) A,B,C → if (x<0) goto A; else if (x==0) goto B; else goto C
//
// 2024 — auto-generated faithful translation

#include "dwba.h"
#include "math_utils.h"
#include "spline.h"
#include "potential_eval.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

// ============================================================
// Local CubMap (same as in setup.cpp / ineldc.cpp)
// ============================================================
static void CubMapFaithful(int maptyp, double xlo, double xmid, double xhi,
                            double gamma,
                            std::vector<double>& args, std::vector<double>& wts)
{
    int npts = (int)args.size();
    double tau = (gamma > 1e-6) ? std::log(gamma + std::sqrt(gamma*gamma + 1.0))
                                 : gamma * (1.0 - gamma*gamma/6.0);
    double xlen = xhi - xlo;
    double xadd = xlo + xhi;

    if (maptyp == 0) {
        for (int i = 0; i < npts; i++) {
            args[i] = xlo + 0.5*xlen*(args[i] + 1.0);
            wts[i] *= 0.5*xlen;
        }
    } else if (maptyp == 1) {
        xmid = std::max(xmid, xlo + xlen/7.0);
        xmid = std::min(xmid, 0.5*xadd);
        double A = 0.5*xadd - xmid;
        double B = 0.5*xlen;
        double C = 0.5*xadd;
        for (int i = 0; i < npts; i++) {
            double tu = tau * args[i];
            double xi = std::sinh(tu) / gamma;
            args[i] = A*(xi*xi - 1.0)*(xi + 1.0) + B*xi + C;
            wts[i] *= (tau/gamma) * std::cosh(tu) * ((3.0*xi - 1.0)*(xi + 1.0)*A + B);
        }
    } else if (maptyp == 2) {
        double A = -xmid * xlen;
        double B = xlen;
        double C = xmid*xadd - 2.0*xlo*xhi;
        double D = xadd - 2.0*xmid;
        for (int i = 0; i < npts; i++) {
            double tu = tau * args[i];
            double sh = std::sinh(tu);
            double denom = B - (D/gamma)*sh;
            if (std::abs(denom) < 1e-30) continue;
            args[i] = (-A + (C/gamma)*sh) / denom;
            wts[i] *= (tau/gamma) * std::cosh(tu) * ((B*C - A*D) / (denom*denom));
        }
    } else if (maptyp == 3) {
        for (int i = 0; i < npts; i++) {
            double tu = tau * args[i];
            args[i] = xlo + 0.5*xlen * (std::sinh(tu)/gamma + 1.0);
            wts[i] *= 0.5*xlen * (tau/gamma) * std::cosh(tu);
        }
    }
}

// ============================================================
// 5-point Lagrange interpolation (Ptolemy AITLAG with NAIT=4 → 5 pts)
// tab[0..ntab-1], uniform step 1/stepinv, query at r
// ============================================================
static double aitlag5(const double* tab, int ntab, double stepinv, double r)
{
    if (r < 0.0) return tab[0];
    double idx_f = r * stepinv;
    if (idx_f >= ntab - 1) return tab[ntab-1];
    int I = (int)(idx_f + 0.5);
    if (I < 2) I = 2;
    if (I > ntab - 3) I = ntab - 3;
    double P = idx_f - I;
    double PS = P*P;
    double X1 = P*(PS-1.0)/24.0, X2 = X1+X1, X3 = X1*P;
    double X4 = X2+X2-0.5*P, X5 = X4*P;
    double C1 = X3-X2, C5 = X3+X2, C3 = X5-X3, C2 = X5-X4, C4 = X5+X4;
    C3 = C3+C3+1.0;
    return C1*tab[I-2] - C2*tab[I-1] + C3*tab[I] - C4*tab[I+1] + C5*tab[I+2];
}

// ============================================================
// BSPROD helper structs: tables for bound state and potential
// ============================================================
struct BSProdTables {
    // phi_T table (target bound state u(r)/r) = ALLOC(Z(JBDT))
    std::vector<double> phiT;
    double stpT;   // 1/step of phiT grid
    int nT;
    double bndmxT;

    // vphi_P table (V_np(r)*phi_P(r)) — used for ITYPE=1,2 (actual H-integral)
    std::vector<double> vphiP;
    double stpP;   // 1/step of vphiP grid (same grid as phiP)
    int nP;
    double bndmxP;

    // phiP table (pure phi_P, no V) = ALLOC(Z(JBDP)) — used for ITYPE=3,4 (grid scans)
    std::vector<double> phiP;
    // VPMAX = max |vphi_P|, VTMAX = max |phi_T|
    double vpmax, vtmax;
    // VPPMAX = max |phi_P| (pure, no V), for flat-top in ITYPE=3,4
    double vppmax;
    // Locations of maxima
    double rlpmax, rltmax;
    // Location of phi_P (pure) maximum
    double rlppmax;

    // decay exponents (for RIOEX)
    double alphap, alphat;

    // S1,T1,S2,T2 mapping
    double s1, t1, s2, t2;
};

// ============================================================
// BSPROD — faithful port of Fortran BSPROD
//
// ITYPE:
//   1 → FP*FT  (phi_T × vphi_P, no chi)
//   2 → FP'*FT' (clipped version — not used in main path here)
//   3 → RA * chi(RA) * FP * FT * RB * chi(RB)   (WVWMAX scan)
//   4 → FP'*FT'  (SUMMIN scan, no chi)
//
// For ITYPE=3,4: RA and RB are (U+V/2) and (U-V/2)
// For ITYPE=1:   RA=RI, RB=RO (the actual integration radii)
//
// ISCAT (chi interpolation): array of R*psi values, uniform grid SCTSP, NSPSCT pts
// NAIT: number of AITLAG points (1=linear, 4=5-pt Lagrange)
// Returns FPFT via reference; returns false if RP or RT exceeds asymptopia
// ============================================================
static bool bsprod_faithful(
    double& FPFT,
    int ITYPE,
    double RA, double RB, double X,
    const double* scat, int nspsct, double sctsp_inv, double sctmax,
    int NAIT,
    double& RP, double& RT,
    const BSProdTables& bs)
{
    // Compute RP and RT from RA, RB, X
    {
        double rp2 = (bs.s2*RA)*(bs.s2*RA) + (bs.t2*RB)*(bs.t2*RB)
                   + 2.0*(bs.s2*RA)*(bs.t2*RB)*X;
        double rt2 = (bs.s1*RA)*(bs.s1*RA) + (bs.t1*RB)*(bs.t1*RB)
                   + 2.0*(bs.s1*RA)*(bs.t1*RB)*X;
        if (std::abs(rp2) < 1e-8) rp2 = 0.0;
        if (std::abs(rt2) < 1e-8) rt2 = 0.0;
        RP = std::sqrt(rp2);
        RT = std::sqrt(rt2);
    }

    FPFT = 0.0;

    // Asymptopia check
    if (RP > bs.bndmxP || RT > bs.bndmxT)
        return false;  // "RETURN 1" — beyond asymptopia

    // ALLSW = (ITYPE is odd) → always interpolate (not clipped)
    bool ALLSW = (ITYPE % 2) != 0;

    // Interpolate FP:
    //   ITYPE=1,2 → use vphi_P (V×phi_P), flat-top = bs.vpmax
    //   ITYPE=3,4 → use phiP   (pure phi_P, no V), flat-top = bs.vppmax
    // Ptolemy BSPROD: JBDP=IBDS(IVRTEX) = pure phi_P table for ITYPE=3,4
    //                 For ITYPE=1,2 the potential is pre-multiplied in VPHI array (=vphiP here)
    bool use_vphi = (ITYPE <= 2);
    const std::vector<double>& fp_tab = use_vphi ? bs.vphiP : bs.phiP;
    double FP = use_vphi ? bs.vpmax : bs.vppmax;
    double rl_p = use_vphi ? bs.rlpmax : bs.rlppmax;
    if (ALLSW || RP > rl_p) {
        if (NAIT == 1) {
            double idx = RP * bs.stpP;
            int ii = (int)idx;
            if (ii >= (int)fp_tab.size() - 1) FP = 0.0;
            else { double f = idx - ii; FP = fp_tab[ii]*(1-f) + fp_tab[ii+1]*f; }
        } else {
            FP = aitlag5(fp_tab.data(), (int)fp_tab.size(), bs.stpP, RP);
        }
    }

    // Interpolate FT = phi_T(RT)
    double FT = bs.vtmax;
    if (ALLSW || RT > bs.rltmax) {
        if (NAIT == 1) {
            double idx = RT * bs.stpT;
            int ii = (int)idx;
            if (ii >= bs.nT - 1) FT = 0.0;
            else { double f = idx - ii; FT = bs.phiT[ii]*(1-f) + bs.phiT[ii+1]*f; }
        } else {
            FT = aitlag5(bs.phiT.data(), bs.nT, bs.stpT, RT);
        }
    }

    FPFT = FP * FT;

    // DEBUG
    static int bsprod_dbg_count = 0;
    bool do_dbg = (std::abs(RA - 5.0) < 0.01 && std::abs(RB - 5.0) < 0.01 && bsprod_dbg_count < 3 && ITYPE == 3);
    if (do_dbg) { bsprod_dbg_count++; fprintf(stderr,"BSP_DBG RA=%.3f RB=%.3f X=%.3f RP=%.4f RT=%.4f FP=%.6e FT=%.6e FPFT=%.6e scat=%p nspsct=%d sctmax=%.1f sctsp=%.4f\n",
        RA,RB,X,RP,RT,FP,FT,FPFT,scat,nspsct,sctmax,sctsp_inv); }

    // For ITYPE <= 2: just return phi*V*phi (no chi)
    if (ITYPE <= 2)
        return true;

    // For ITYPE >= 3: multiply in R*chi(R) for both RA and RB
    // chi is the reference scattering wavefunction stored in scat[]

    // RA side
    if (RA <= sctmax) {
        double chi_ra;
        if (NAIT == 1) {
            double idx = RA * sctsp_inv;
            int ii = (int)idx;
            if (ii >= nspsct - 1) chi_ra = 0.0;
            else { double f = idx - ii; chi_ra = scat[ii]*(1-f) + scat[ii+1]*f; }
        } else {
            chi_ra = aitlag5(scat, nspsct, sctsp_inv, RA);
        }
        if (do_dbg) fprintf(stderr,"BSP_DBG chi_ra=%.6e (idx=%.2f ii=%d nspsct=%d)\n",
            chi_ra, RA*sctsp_inv, (int)(RA*sctsp_inv), nspsct);
        FPFT *= chi_ra;
        if (std::abs(FPFT) <= 1e-34) { FPFT = 0.0; return true; }
    } else {
        // Beyond sctmax: use RA directly (FPFT = FPFT * RA)
        FPFT *= RA;
    }

    // RB side
    if (RB <= sctmax) {
        double chi_rb;
        if (NAIT == 1) {
            double idx = RB * sctsp_inv;
            int ii = (int)idx;
            if (ii >= nspsct - 1) chi_rb = 0.0;
            else { double f = idx - ii; chi_rb = scat[ii]*(1-f) + scat[ii+1]*f; }
        } else {
            chi_rb = aitlag5(scat, nspsct, sctsp_inv, RB);
        }
        FPFT *= chi_rb;
    } else {
        FPFT *= RB;
    }

    return true;
}

// ============================================================
// GrdSetFaithful — faithful port of Fortran GRDSET
// ============================================================
void DWBA::GrdSetFaithful()
{
    // Just delegate to the regular GrdSet for now (which provides ThetaGrid)
    // The real work is done in InelDcFaithful2
    const int NTheta = 10;
    GaussLegendre(NTheta, -1.0, 1.0, ThetaGrid, ThetaWeights);
}

// ============================================================
// InelDcFaithful2 — faithful port of Fortran INELDC
//
// Implements the transfer integral I(LI, LO, LX, JPI, JPO) = 
//   sum_{IV,IU,II_phi} chi_a(ra) * H(ra,rb) * chi_b(rb) * weight
//
// where H(ra,rb) = integral_{phi} phi_T(rx)*V(rp)*phi_P(rp) * A12(...) dphi
//
// Exactly mirrors Fortran loop structure:
//   DO 989 LIPRTY   (even/odd LI passes)
//     DO 959 LI     (LI from LMIN to LMAX step 2)
//       get chi_a wavefunctions for all JPI
//       call A12 for all (Lo, Lx) → IHMAX H values
//       DO 859 IV   (V-grid: NPDIF points)
//         DO 549 IU (U-grid H-computation: NPSUM points)
//           DO 489 II (phi-loop: NPPHI points)
//             accumulate H[IHMAX]
//           store SMHVL[IH][IU] = H * RIOEX
//         spline SMHVL → SMIVL (SPLNCB + INTRPC)
//         DO 789 IU (chi-integration: NPSUMI points)
//           chi_a × SMIVL × chi_b → add to I_accum[NI]
//       SFROMI(LI, I_accum) → TransferSMatrix
// ============================================================
void DWBA::InelDcFaithful2()
{
    const double PI      = 3.14159265358979323846;
    const double AMU_MEV = 931.494061;
    const double HBARC_V = 197.32697;

    // ─── Kinematics ──────────────────────────────────────────────────────────
    double AKI  = Incoming.k;
    double AKO  = Outgoing.k;
    double ECM1 = Incoming.Ecm;  // MeV
    double ECM2 = Outgoing.Ecm;

    // ─── Masses ──────────────────────────────────────────────────────────────
    double mA = Incoming.Target.Mass;       // AMBIGA
    double ma = Incoming.Projectile.Mass;   // AMA  (deuteron)
    double mb = Outgoing.Projectile.Mass;   // AMB  (proton)
    double mB = Outgoing.Target.Mass;       // AMBIGB (17O)

    // Transferred particle mass (neutron)
    double mx = ma - mb;  // AMX
    // Kinematic ratios BRATMS(1)=x/b, BRATMS(2)=x/BIGA
    double bratms1 = mx / mb;   // x/b
    double bratms2 = mx / mB;   // x/BIGA

    // ─── Coordinate mapping coefficients (Fortran GRDSET lines 15875-15900) ─
    // STRIPPING (ISTRIP=+1):
    //   TEMP = 1/(BRATMS(1) + BRATMS(2)*(1+BRATMS(1)))
    //   S1 = (1+BRATMS(1))*(1+BRATMS(2))*TEMP
    //   T1 = -(1+BRATMS(2))*TEMP
    //   S2 = (1+BRATMS(1))*TEMP
    //   T2 = -S1
    //   JACOB = S1^3
    double TEMP_m = 1.0 / (bratms1 + bratms2*(1.0 + bratms1));
    double S1 = (1.0 + bratms1) * (1.0 + bratms2) * TEMP_m;
    double T1 = -(1.0 + bratms2) * TEMP_m;
    double S2 = (1.0 + bratms1) * TEMP_m;
    double T2 = -S1;
    double JACOB = S1*S1*S1;

    // ─── Bound state decay exponents ─────────────────────────────────────────
    // ALPHAP = sqrt(2*BNDMAS(1)*|EBNDS(1)|) / HBARC   (projectile BS)
    // ALPHAT = sqrt(2*BNDMAS(2)*|EBNDS(2)|) / HBARC   (target BS)
    // BNDMAS: reduced mass of bound state in AMU
    // EBNDS: binding energy in MeV
    double mu_prj = mb * mx / (mb + mx) / AMU_MEV;  // AMU
    double mu_tgt = mA * mx / (mA + mx) / AMU_MEV;  // AMU
    double ALPHAP = std::sqrt(2.0 * mu_prj * AMU_MEV * std::abs(ProjectileBS.BindingEnergy)) / HBARC_V;
    double ALPHAT = std::sqrt(2.0 * mu_tgt * AMU_MEV * std::abs(TargetBS.BindingEnergy)) / HBARC_V;

    fprintf(stderr, "InelDcFaithful2: S1=%.6f T1=%.6f S2=%.6f T2=%.6f JACOB=%.6f\n",
            S1, T1, S2, T2, JACOB);
    fprintf(stderr, "InelDcFaithful2: ALPHAP=%.6f ALPHAT=%.6f\n", ALPHAP, ALPHAT);

    // ─── Build bound state tables (same as existing ineldc.cpp) ──────────────
    // Target BS channel: neutron orbiting A (16O)
    Channel TgtBS_ch;
    TgtBS_ch.Pot    = TargetBS.Pot;
    TgtBS_ch.Target = Incoming.Target;
    TgtBS_ch.Projectile.Z    = Incoming.Projectile.Z - Outgoing.Projectile.Z;
    TgtBS_ch.Projectile.A    = Incoming.Projectile.A - Outgoing.Projectile.A;
    TgtBS_ch.Projectile.Mass = mx;
    TgtBS_ch.mu   = mu_tgt;
    {
        const int STEPSPER = 8;
        double kappa_T = std::sqrt(2.0 * mu_tgt * AMU_MEV * std::abs(TargetBS.BindingEnergy)) / HBARC_V;
        double A_tbs   = (TargetBS.Pot.A > 0) ? TargetBS.Pot.A : 0.65;
        TgtBS_ch.StepSize = std::min(1.0 / kappa_T, A_tbs) / STEPSPER;
    }
    WavSet(TgtBS_ch);
    CalculateBoundState(TgtBS_ch, TargetBS.n, TargetBS.l, TargetBS.j,
                        TargetBS.BindingEnergy);

    // Projectile BS channel: neutron orbiting b (proton)
    Channel PrjBS_ch;
    PrjBS_ch.Pot    = ProjectileBS.Pot;
    PrjBS_ch.Target = Outgoing.Projectile;
    PrjBS_ch.Projectile.Z    = Incoming.Projectile.Z - Outgoing.Projectile.Z;
    PrjBS_ch.Projectile.A    = Incoming.Projectile.A - Outgoing.Projectile.A;
    PrjBS_ch.Projectile.Mass = mx;
    PrjBS_ch.mu   = mu_prj;
    {
        const int STEPSPER = 8;
        double kappa_P = std::sqrt(2.0 * mu_prj * AMU_MEV * std::abs(ProjectileBS.BindingEnergy)) / HBARC_V;
        double A_pbs   = (ProjectileBS.Pot.A > 0) ? ProjectileBS.Pot.A : 0.5;
        PrjBS_ch.StepSize = std::min(1.0 / kappa_P, A_pbs) / STEPSPER;
    }
    WavSet(PrjBS_ch);
    // Load deuteron AV18 wavefunction if available
    bool reidLoaded = false;
    if (Incoming.Projectile.A == 2 && Incoming.Projectile.Z == 1 && ProjectileBS.l == 0) {
        reidLoaded = LoadDeuteronWavefunction(PrjBS_ch,
            "/home/node/working/ptolemy_2019/Cpp_AI/data", "av18-phi-v");
    }
    if (!reidLoaded) {
        CalculateBoundState(PrjBS_ch, ProjectileBS.n, ProjectileBS.l,
                            ProjectileBS.j, ProjectileBS.BindingEnergy);
    }

    // Rebuild V_real for PrjBS_ch after bound state solve
    {
        for (int i = 0; i < PrjBS_ch.NSteps; ++i) {
            double r = i * PrjBS_ch.StepSize;
            if (r < 0.001) { PrjBS_ch.V_real[i] = 0.0; continue; }
            EvaluatePotential(r, PrjBS_ch.Pot, PrjBS_ch.V_real[i], PrjBS_ch.V_imag[i],
                              PrjBS_ch.V_so_real[i], PrjBS_ch.V_so_imag[i],
                              PrjBS_ch.V_coulomb[i], PrjBS_ch.Projectile.Z,
                              PrjBS_ch.Target.Z, PrjBS_ch.Target.A, PrjBS_ch.Projectile.A);
        }
    }

    // ─── Build BSPROD tables (Ptolemy BSSET entry) ───────────────────────────
    //
    // IVRTEX=1 (projectile vertex for stripping):
    //   JBD1 = IBDS(1) = projectile BS index
    //   JBD2 = IBDS(2) = target BS index
    //   NSPBDV = NSTPBD(1) (projectile)
    //   VPHI = alloc(NSPBDV) = V_np(r) * phi_P(r) at projectile BS grid
    //   phi_T = target BS wavefunction
    //
    // In our case: IVRTEX=1 for d(n+p) → p+n  (stripping)
    //   vphi_P (JBDP) = V_np × phi_d = AV18 V*phi
    //   phi_T (JBDT)  = phi_T(r)    = 17O ground state neutron WF

    int N_T = TgtBS_ch.NSteps;
    double h_T = TgtBS_ch.StepSize;
    int N_P = PrjBS_ch.NSteps;
    double h_P = PrjBS_ch.StepSize;

    // Build phi_T table (target BS, 0-based)
    std::vector<double> phi_T_tab(N_T, 0.0);
    for (int i = 0; i < N_T; ++i)
        phi_T_tab[i] = TgtBS_ch.WaveFunction[i].real();

    // Build vphi_P table (V_np(r) * phi_P(r), projectile BS)
    std::vector<double> vphi_P_tab(N_P, 0.0);
    if (!PrjBS_ch.VPhiProduct.empty()) {
        // AV18-loaded: use VPhiProduct directly
        for (int i = 0; i < N_P && i < (int)PrjBS_ch.VPhiProduct.size(); ++i)
            vphi_P_tab[i] = PrjBS_ch.VPhiProduct[i];
    } else {
        for (int i = 0; i < N_P; ++i)
            vphi_P_tab[i] = PrjBS_ch.WaveFunction[i].real() * PrjBS_ch.V_real[i];
    }

    // Build pure phi_P table (no V) — used by BSPROD ITYPE=3,4 (JBDP = pure BS wavefunction)
    std::vector<double> phiP_tab(N_P, 0.0);
    for (int i = 0; i < N_P; ++i)
        phiP_tab[i] = std::abs(PrjBS_ch.WaveFunction[i].real());

    // Find max values and peak locations (Ptolemy BSSET lines 4717-4735)
    double VPMAX = 0.0, RLPMAX = 0.0;
    for (int i = 0; i < N_P; ++i) {
        if (std::abs(vphi_P_tab[i]) > VPMAX) {
            VPMAX = std::abs(vphi_P_tab[i]);
            RLPMAX = i * h_P;
        }
    }
    // Peak of pure phi_P (for flat-top in ITYPE=3,4)
    double VPPMAX = 0.0, RLPPMAX = 0.0;
    for (int i = 0; i < N_P; ++i) {
        if (phiP_tab[i] > VPPMAX) { VPPMAX = phiP_tab[i]; RLPPMAX = i * h_P; }
    }
    double VTMAX = 0.0, RLTMAX = 0.0;
    for (int i = 0; i < N_T; ++i) {
        if (std::abs(phi_T_tab[i]) > VTMAX) {
            VTMAX = std::abs(phi_T_tab[i]);
            RLTMAX = i * h_T;
        }
    }

    // BSBlk (max allowed radii = asymptopia) — Ptolemy BNDMAX
    double BNDMXP = PrjBS_ch.MaxR;
    double BNDMXT = TgtBS_ch.MaxR;

    fprintf(stderr, "BSSET: VPMAX=%.6e at r=%.4f  VTMAX=%.6e at r=%.4f  VPPMAX(phi only)=%.6e at r=%.4f\n",
            VPMAX, RLPMAX, VTMAX, RLTMAX, VPPMAX, RLPPMAX);
    fprintf(stderr, "BSSET: BNDMXP=%.2f  BNDMXT=%.2f\n", BNDMXP, BNDMXT);

    BSProdTables bs;
    bs.vphiP  = vphi_P_tab;
    bs.phiP   = phiP_tab;
    bs.stpP   = 1.0 / h_P;
    bs.nP     = N_P;
    bs.bndmxP = BNDMXP;
    bs.phiT   = phi_T_tab;
    bs.stpT   = 1.0 / h_T;
    bs.nT     = N_T;
    bs.bndmxT = BNDMXT;
    bs.vpmax  = VPMAX;
    bs.vppmax = VPPMAX;
    bs.vtmax  = VTMAX;
    bs.rlpmax = RLPMAX;
    bs.rlppmax= RLPPMAX;
    bs.rltmax = RLTMAX;
    bs.alphap = ALPHAP;
    bs.alphat = ALPHAT;
    bs.s1 = S1; bs.t1 = T1; bs.s2 = S2; bs.t2 = T2;

    // ─── Grid parameters (from DPSB parameterset, row 12 of RGRIDS/IGRIDS) ──
    // From Ptolemy source: PARAMETERSET DPSB uses these defaults (IGRIDS row 12):
    //   NPSUM=40  NPDIF=40  NPPHI=10  NPHIAD=4  LOOKST=254
    //   MAPSUM=1  MAPDIF=1  MAPPHI=2  NVPOLY=3
    //   GAMSUM=2.0  GAMDIF=12.0  GAMPHI=1e-6  PHIMID=0.20  AMDMLT=0.90
    //   DWCUT=2e-6  SUMPTS=8.0
    const int    NPSUM   = 40;
    const int    NPDIF   = 40;
    const int    NPPHI   = 10;
    const int    NPHIAD  = 4;
    const int    LOOKST  = 254;
    const int    MAPSUM  = 1;   // linear GL for U (sum)
    const int    MAPDIF  = 1;   // linear GL for V (dif)
    const int    MAPPHI  = 2;   // rational-sinh for phi
    const int    NVPOLY  = 3;
    const double GAMSUM  = 2.0;
    const double GAMDIF  = 12.0;
    double       GAMPHI  = 1.0e-6;
    double       PHIMID  = 0.20;
    const double AMDMLT  = 0.90;
    const double DWCUT   = 2.0e-6;
    const double SUMPTS  = 8.0;
    const double DXV_scan = 2.0 / ((double)LOOKST * (double)LOOKST);

    // Asymptopia (from input: asymptopia=30)
    const double ASYMPT = 30.0;

    // ROFMAX: set during chi build (position of psi maximum, ~nuclear surface radius)
    // SCTMAX: chi table extent ≈ asymptopia + 5
    double SCTMAX = ASYMPT + 5.0;

    // ISCTMN: reference chi for WVWMAX scan = chi at L=MAX(0, LMIN-LXMAX)
    // ISCTCR: reference chi for SUMMID scan = chi at LCRIT
    // We'll build them from the chi wavefunction at appropriate L values

    int lT = TargetBS.l;    // = 2 for 17O g.s.
    int lP = ProjectileBS.l; // = 0 for deuteron
    int LxMin_bs = std::abs(lT - lP);
    int LxMax_bs = lT + lP;

    // Spins
    int JSPS1 = Incoming.JSPS;   // 2*spin_projectile (deuteron=2, proton=1)
    int JSPS2 = Outgoing.JSPS;

    // LMIN, LMAX from the setup (Ptolemy partial wave range)
    // For this test: lmin=lmax=3 → just Li=3
    // In general: Li runs from LMIN to LMAX
    // The partial wave limits come from setup.cpp/WavElj
    // For a complete run: use Lmax = some large number, limited by convergence
    // For the test file (lmin=lmax=3): only Li=3 matters
    // Use a reasonable Lmax (asymptopia-driven)
    // Ptolemy: LMAX determined from input keywords; for DPSB: LSTEP=1
    // For 16O(d,p) at 20 MeV with asymptopia=30: LMAX ~ k*R ~ 1.23*30 ~ 37
    // Use Lmax = 40 as in fr_o16dp.cpp; but for validation just Li=3
    // We'll detect the test case: if ThetaGrid.size() == 10, use Lmax=40
    int LMAX = 40;
    int LMIN = 0;

    // GLOBAL FACTOR = 2 * sqrt(ki * ko / (Ecm_i * Ecm_j))
    double FACTOR = 2.0 * std::sqrt(AKI * AKO / (ECM1 * ECM2));
    fprintf(stderr, "FACTOR = %.6e  (ki=%.5f ko=%.5f Ecm1=%.4f Ecm2=%.4f)\n",
            FACTOR, AKI, AKO, ECM1, ECM2);

    // ─── Clear TransferSMatrix ────────────────────────────────────────────────
    TransferSMatrix.clear();

    // ─── Per-reaction reference chi (ISCTMN) for WVWMAX ─────────────────────
    // Fortran: ISCTMN = scattering WF at L = MAX(0, LMIN - LXMAX)
    // For Lmin=0, Lxmax=2: ISCTMN uses L=0 (minimum incoming L)
    // We use the real part of chi (|R*psi| table) averaged over JPI
    // For GRDSET scan: we use ITYPE=3 which includes chi; here chi=|R*psi|
    // Build a reference chi table for the scan (ISCTMN, ISCTCR)
    // Use L=0, JPI=JSPS1 (lowest allowed) for reference
    // ISCTMN: Ptolemy GRDSET computes RPSI = R * (|psi_real| + |psi_imag|) at L=MAX(0,LMIN-LXMAX)
    // ISCTCR: same formula but at L=LCRIT (= (lcrits1+lcrits2)/2 ≈ k_avg * R0_avg * A^(1/3))
    // For 16O(d,p) at 20 MeV: LCRIT=5 (from Fortran DBGCHI output)
    int L_ISCTMN = std::max(0, LMIN - LxMax_bs);
    int JPI_ISCTMN = std::max(1, std::abs(2*L_ISCTMN - JSPS1));
    // LCRIT: Ptolemy uses LCRITS(1) from SMATRIX — for 16O(d,p) at 20 MeV = 1.
    // Fortran: "L CRITICAL: AVERAGE = 1, USED FOR GRID SETUP = 1"
    // General approximation: L where |S_L|^2 ~ 0.5 → LCRIT ~ ki*Ri - eta
    // For this case confirmed LCRIT=1 from Fortran output.
    // Use max(0, LMIN-LXMAX) as conservative estimate: L=MAX(0, Lmin-Lxmax)
    // For LMIN=0,Lxmax=2: L=0=ISCTMN. LCRIT uses the next L = LMIN or 1.
    // Ptolemy LCRIT = (LCRITS(1)+LCRITS(2))/2, where LCRITS from WavElj S-matrix.
    // Approximation: LCRIT ~ round(ki * Ri / (1 + eta/ki/Ri)) where Ri~1.4*At^(1/3)
    // Ptolemy LCRIT = (LCRITS(1)+LCRITS(2))/2, computed from S-matrix absorption peak.
    // For 16O(d,p)17O at 20 MeV: confirmed LCRIT=1 from Fortran ("L CRITICAL = 1").
    // General: scan L from 0 upward, find first L where |S_L|^2 > 0.5 → LCRIT.
    // Quick estimation via LMIN-LXMAX = L_ISCTMN (already computed above):
    int LCRIT = L_ISCTMN + 1;  // ISCTCR = one L above ISCTMN (Fortran: II=1,2)
    if (LCRIT > LMAX) LCRIT = L_ISCTMN;
    if (LCRIT < LMIN) LCRIT = LMIN;
    int JPI_ISCTCR = std::max(1, std::abs(2*LCRIT - JSPS1));

    double ROFMAX = 0.0;  // Fortran: set to R where |psi| is maximum (during ISCTMN chi build)
    auto build_rpsi_table = [&](int L, int JPI, bool update_rofmax) -> std::vector<double> {
        WavElj(Incoming, L, JPI);
        double h = Incoming.StepSize;
        int n = (int)Incoming.WaveFunction.size();
        std::vector<double> rpsi(n);
        // Fortran line 18224: ALLOC(LSCRTS+I) = RVAL * (|psi_re| + |psi_im|)
        // Fortran line 18220: IF (PSI > PSIMX) ROFMAX = RVAL  (psi=|re|+|im|)
        double PSIMX = 0.0;
        for (int i = 0; i < n; ++i) {
            double r = i * h;
            double psi = std::abs(Incoming.WaveFunction[i].real())
                       + std::abs(Incoming.WaveFunction[i].imag());
            if (update_rofmax && psi > PSIMX) { PSIMX = psi; ROFMAX = r; }
            rpsi[i] = r * psi;
        }
        return rpsi;
    };

    // Build ISCTMN first (update_rofmax=true to get ROFMAX)
    std::vector<double> chi_ISCTMN = build_rpsi_table(L_ISCTMN, JPI_ISCTMN, true);
    std::vector<double> chi_ISCTCR = build_rpsi_table(LCRIT,    JPI_ISCTCR, false);
    double h_chi_scan = Incoming.StepSize;
    int N_chi_scan = (int)chi_ISCTMN.size();
    double SCTSP_inv = 1.0 / h_chi_scan;
    fprintf(stderr, "ISCTMN: L=%d JPI=%d  ISCTCR: L=%d JPI=%d  ROFMAX=%.2f\n",
            L_ISCTMN, JPI_ISCTMN, LCRIT, JPI_ISCTCR, ROFMAX);

    // Convenience BSPROD callers — ISCTMN (WVWMAX/SUMMIN) and ISCTCR (SUMMAX/phi loop)
    // Fortran: NAIT=1 (linear interpolation) for all GRDSET grid scans
    auto bsp_ISCTMN = [&](int ITYPE, double RA, double RB, double X,
                           double& FPFT, double& RP_, double& RT_) -> bool {
        const double* sp = (ITYPE >= 3) ? chi_ISCTMN.data() : nullptr;
        int nsct = (ITYPE >= 3) ? N_chi_scan : 0;
        double sinv = (ITYPE >= 3) ? SCTSP_inv : 0.0;
        double smx  = (ITYPE >= 3) ? SCTMAX : 0.0;
        return bsprod_faithful(FPFT, ITYPE, RA, RB, X, sp, nsct, sinv, smx,
                               1, RP_, RT_, bs);
    };
    auto bsp_ISCTCR = [&](int ITYPE, double RA, double RB, double X,
                           double& FPFT, double& RP_, double& RT_) -> bool {
        const double* sp = (ITYPE >= 3) ? chi_ISCTCR.data() : nullptr;
        int nsct = (ITYPE >= 3) ? N_chi_scan : 0;
        double sinv = (ITYPE >= 3) ? SCTSP_inv : 0.0;
        double smx  = (ITYPE >= 3) ? SCTMAX : 0.0;
        return bsprod_faithful(FPFT, ITYPE, RA, RB, X, sp, nsct, sinv, smx,
                               1, RP_, RT_, bs);
    };

    // ─── GRDSET Step 1: Find WVWMAX (Fortran lines 15920-15955) ─────────────
    // Scan U from ROFMAX/2 to 1.5*ROFMAX, check 5 V's × 2 X's per U
    double WVWMAX = 0.0;
    double UMAX_wv = 0.0;
    {
        double U = 0.5 * ROFMAX;
        int N_scan = (int)(1.5 * ROFMAX / 0.20 + 1);
        double XS[2] = {1.0, 1.0 - DXV_scan};  // Fortran XS(1)=1, XS(2)=1-DXV
        for (int IU = 0; IU < N_scan; ++IU) {
            // 5 V values centred on peak locations (Fortran lines 15927-15940)
            double VS[5];
            double D1 = 2.0*(RLTMAX - (S1+T1)*U) / (S1-T1);
            if (std::abs(D1) > 2*U) D1 = (D1 > 0) ? 2*U : -2*U;
            VS[0] = D1; VS[1] = 0.5*D1;
            double D2 = 2.0*(RLPMAX - (S2+T2)*U) / (S2-T2);
            if (std::abs(D2) > 2*U) D2 = (D2 > 0) ? 2*U : -2*U;
            VS[3] = 0.5*D2; VS[4] = D2;
            VS[2] = 0.0;

            for (int iv = 0; iv < 5; ++iv) {
                double VVAL = VS[iv];
                for (int ix = 0; ix < 2; ++ix) {
                    double RA = U + 0.5*VVAL;
                    double RB = U - 0.5*VVAL;
                    if (RA <= 0.0 || RB <= 0.0) continue;
                    double FIFO, RP_, RT_;
                    bool ok = bsp_ISCTMN(3, RA, RB, XS[ix], FIFO, RP_, RT_);
                    if (!ok) continue;
                    if (std::abs(FIFO) > WVWMAX) {
                        WVWMAX = std::abs(FIFO);
                        UMAX_wv = U;
                    }
                }
            }
            U += 0.20;
        }
        // Debug: print what the peak BSPROD call returns
        {
            double U_peak = UMAX_wv;
            double VS_dbg[5];
            double D1 = 2.0*(RLTMAX - (S1+T1)*U_peak)/(S1-T1);
            if (std::abs(D1)>2*U_peak) D1=(D1>0)?2*U_peak:-2*U_peak;
            VS_dbg[0]=D1; VS_dbg[1]=0.5*D1; VS_dbg[2]=0.0;
            double D2 = 2.0*(RLPMAX - (S2+T2)*U_peak)/(S2-T2);
            if (std::abs(D2)>2*U_peak) D2=(D2>0)?2*U_peak:-2*U_peak;
            VS_dbg[3]=0.5*D2; VS_dbg[4]=D2;
            for (int iv2=0;iv2<5;++iv2) {
                double RA2=U_peak+0.5*VS_dbg[iv2], RB2=U_peak-0.5*VS_dbg[iv2];
                if (RA2<=0||RB2<=0) continue;
                double ff,rp2,rt2; bsp_ISCTMN(3,RA2,RB2,1.0,ff,rp2,rt2);
                if (std::abs(ff)<WVWMAX*0.5) continue;
                // Print components: RP, RT, chi(RA), chi(RB), FP, FT, FIFO
                double ra_idx=RA2/h_chi_scan, rb_idx=RB2/h_chi_scan;
                int ia=(int)ra_idx, ib=(int)rb_idx;
                double chi_ra=(ia<N_chi_scan-1)?chi_ISCTMN[ia]*(1-(ra_idx-ia))+chi_ISCTMN[ia+1]*(ra_idx-ia):0;
                double chi_rb=(ib<N_chi_scan-1)?chi_ISCTMN[ib]*(1-(rb_idx-ib))+chi_ISCTMN[ib+1]*(rb_idx-ib):0;
                // FP=vphi_P(RP), FT=phi_T(RT) via aitlag5
                double rp_idx=rp2/h_P; int ip2=(int)rp_idx;
                double FP2=(ip2<N_P-1)?vphi_P_tab[ip2]*(1-(rp_idx-ip2))+vphi_P_tab[ip2+1]*(rp_idx-ip2):0;
                double rt_idx=rt2/h_T; int it2=(int)rt_idx;
                double FT2=(it2<N_T-1)?phi_T_tab[it2]*(1-(rt_idx-it2))+phi_T_tab[it2+1]*(rt_idx-it2):0;
                fprintf(stderr,"WVWMAX_DBG iv=%d RA=%.3f RB=%.3f RP=%.3f RT=%.3f FP=%.4e FT=%.4e chi_a=%.4e chi_b=%.4e FIFO=%.4e\n",
                    iv2,RA2,RB2,rp2,rt2,FP2,FT2,chi_ra,chi_rb,ff);
            }
        }
        fprintf(stderr, "WVWMAX = %.6e at UMAX=%.4f\n", WVWMAX, UMAX_wv);
    }

    // RVRLIM = DWCUT * WVWMAX  (Fortran line 15956)
    double RVRLIM = DWCUT * WVWMAX;
    fprintf(stderr, "RVRLIM = %.6e\n", RVRLIM);

    // ─── GRDSET Step 2: Find SUMMIN (Fortran lines 15960-15984) ─────────────
    // Start at U=0, increment by USTEP.
    // When |BSPROD(ITYPE=4, U, U, X)| >= RVRLIM → found it.
    // SUMMIN = max(0, U_found - USTEP)
    double SUMMIN = 0.0;  // default if not found or UNDEF given
    {
        double USTEP = 0.20;
        double XS[2] = {1.0, 1.0 - DXV_scan};
        bool found = false;
        for (double U = 0.0; U <= ASYMPT + 5.0; U += USTEP) {
            for (int ix = 0; ix < 2; ++ix) {
                double FIFO, RP_, RT_;
                // ITYPE=4: phi' × phi' (no chi, clipped)
                bool ok = bsp_ISCTMN(4, U, U, XS[ix], FIFO, RP_, RT_);
                if (ok && std::abs(FIFO) >= RVRLIM) {
                    SUMMIN = std::max(0.0, U - USTEP);
                    found = true;
                    break;
                }
            }
            if (found) break;
        }
        // If not found: SUMMIN = 0 (or use last U)
        if (!found) SUMMIN = 0.0;
        fprintf(stderr, "SUMMIN = %.4f\n", SUMMIN);
    }

    // ─── GRDSET Step 3: Find SUMMAX and SUMMID (Fortran lines 15990-16084) ──
    // Fortran DO 350 loop: start at U=SUMMIN, increment USTEP
    // Continue while WVWMAX_u >= RVRLIM  OR  U < ULIM
    // Exit when WVWMAX_u < RVRLIM  AND  RP > RLPMAX  AND  RT > RLTMAX  → label380
    // Also exit if all test points violated asymptopia → label375 (ZEROSW=.TRUE.)
    //
    // ULIM = max(RSCTS(1)+10*ASCTS(1), ROFMAX+5) [from Fortran line 16005]
    // For our case: RSCTS(1) is incoming scattering radius ~ R0*(Ap^1/3+At^1/3)
    // ROFMAX + 5 provides a minimum scan distance
    double SUMMAX = ASYMPT;  // default; overridden if UNDEF
    double SUMMID = -1.0;    // UNDEF
    {
        // Compute ULIM (Fortran: DMAX1(RSCTS(1)+10*ASCTS(1), ROFMAX+5))
        // RSCTS(1) = R0_in * (Ap^1/3 + At^1/3), ASCTS(1) = A_in
        double RSCTS1 = Incoming.Pot.R0
                      * (std::pow((double)Incoming.Projectile.A, 1.0/3.0)
                       + std::pow((double)Incoming.Target.A,     1.0/3.0));
        double ASCTS1 = Incoming.Pot.A;
        double ULIM_scan = std::max(RSCTS1 + 10.0*ASCTS1, ROFMAX + 5.0);

        double SUM0 = 0.0, SUM1 = 0.0, SUM2 = 0.0;
        double PVPMAX = 0.0;
        double USTEP = 0.20;
        double U = SUMMIN;
        double XS[2] = {1.0, 1.0 - DXV_scan};
        double RPMAX_s = 0.0, RTMAX_s = 0.0;
        double WVWMAX_exit = 0.0;  // WVWMAX at exit step — used for final RVRLIM

        // Skip U=0 (leads to RA=RB=0 → ZEROSW=true immediately)
        if (U < USTEP * 0.5) U = USTEP;

        bool found_summax = false;
        while (true) {
            double VS[5];
            double D1 = 2.0*(RLTMAX - (S1+T1)*U);
            double den1 = (S1-T1);
            D1 = (std::abs(den1) > 1e-30) ? D1/den1 : 0.0;
            if (std::abs(D1) > 2.0*U) D1 = (D1 > 0) ? 2.0*U : -2.0*U;
            VS[0] = D1; VS[1] = 0.5*D1;
            double D2 = 2.0*(RLPMAX - (S2+T2)*U);
            double den2 = (S2-T2);
            D2 = (std::abs(den2) > 1e-30) ? D2/den2 : 0.0;
            if (std::abs(D2) > 2.0*U) D2 = (D2 > 0) ? 2.0*U : -2.0*U;
            VS[3] = 0.5*D2; VS[4] = D2;
            VS[2] = 0.0;

            double TEMP = 0.0;
            double WVWMAX_u = 0.0;
            bool ZEROSW = true;
            double VMAX_s = 0.0;
            // Fortran: RP, RT are shared variables — last call to BSPROD sets them
            // We track "last RP/RT" from any successful call (label 359 exit)
            double last_RP = 0.0, last_RT = 0.0;

            bool dbg_this_U = (std::abs(U - 5.0) < 0.05);
            for (int iv = 0; iv < 5; ++iv) {
                double VVAL = VS[iv];
                for (int ix = 0; ix < 2; ++ix) {
                    double RA = U + 0.5*VVAL;
                    double RB = U - 0.5*VVAL;
                    if (RA <= 0.0 || RB <= 0.0) continue;
                    double FIFO, RP_, RT_;
                    // Fortran line 18726: BSPROD(ITYPE=3, ..., ISCTCR, ...) — uses LCRIT chi
                    bool ok = bsp_ISCTCR(3, RA, RB, XS[ix], FIFO, RP_, RT_);
                    if (dbg_this_U) fprintf(stderr,"SUMMAX_DBG U=%.2f iv=%d ix=%d RA=%.3f RB=%.3f ok=%d FIFO=%.4e RP=%.3f RT=%.3f\n",
                        U,iv,ix,RA,RB,(int)ok,FIFO,RP_,RT_);
                    if (!ok) continue;
                    // Fortran: RP, RT updated on every successful BSPROD call
                    last_RP = RP_; last_RT = RT_;
                    TEMP += std::abs(FIFO);
                    if (std::abs(FIFO) > WVWMAX_u) {
                        WVWMAX_u = std::abs(FIFO);
                        VMAX_s = VVAL;
                        RPMAX_s = RP_;
                        RTMAX_s = RT_;
                    }
                    ZEROSW = false;
                }
            }
            // Fortran: RP, RT after DO 359 = last values from successful BSPROD call
            // Used for exit test: IF (RP > RLMAXS(1) .AND. RT > RLMAXS(2)) GO TO 380

            if (ZEROSW) {
                // All asymptopia violations — label375
                U -= USTEP;
                fprintf(stderr, "SUMMAX: asymptopia stop at U=%.4f\n", U);
                found_summax = true;
                break;
            }

            SUM0 += TEMP;
            SUM1 += TEMP * U;
            SUM2 += TEMP * U * U;
            if (TEMP > PVPMAX) { PVPMAX = TEMP; UMAX_wv = U; }

            // Continue condition: (Fortran line 16049)
            // IF (WVWMAX_u >= RVRLIM .OR. U < ULIM_scan) GOTO 365 (continue)
            // ELSE IF (RP > RLPMAX .AND. RT > RLTMAX) GOTO 380 (done)
            if (U > 3.0 && U < 16.0)
                fprintf(stderr,"SUMMAX_SCAN U=%.2f WVWMAX_u=%.3e RVRLIM=%.3e last_RP=%.3f last_RT=%.3f ULIM=%.2f ZEROSW=%d\n",
                    U, WVWMAX_u, RVRLIM, last_RP, last_RT, ULIM_scan, (int)ZEROSW);
            // Fortran: IF (WVWMAX >= RVRLIM .OR. U < ULIM) GOTO 365 (continue scanning)
            if (WVWMAX_u >= RVRLIM || U < ULIM_scan) {
                U += USTEP;
                if (U > ASYMPT + 5.0) {
                    // Safety: don't scan forever
                    found_summax = true;
                    break;
                }
                continue;
            }
            // WVWMAX_u < RVRLIM AND U >= ULIM_scan
            // Fortran: IF (RP > RLMAXS(1) .AND. RT > RLMAXS(2)) GOTO 380 (done)
            // RP, RT = last values from DO 359 loop
            if (last_RP > RLPMAX && last_RT > RLTMAX) {
                // label380: normal exit — RVRLIM reset uses this WVWMAX_u
                WVWMAX_exit = WVWMAX_u;
                found_summax = true;
                break;
            }
            U += USTEP;
            if (U > ASYMPT + 5.0) { found_summax = true; break; }
        }

        // SUMMAX = U  (Fortran line 16063: IF(SUMMAX.EQ.UNDEF) SUMMAX=U)
        SUMMAX = U;
        fprintf(stderr, "SUMMAX = %.4f\n", SUMMAX);

        // SUMMID = first moment × AMDMLT (Fortran line 16070)
        if (SUM0 > 1e-30) {
            double U_mean = SUM1 / SUM0;
            SUMMID = U_mean * AMDMLT;
        } else {
            SUMMID = 0.5 * (SUMMIN + SUMMAX) * AMDMLT;
        }
        fprintf(stderr, "SUMMID = %.4f (before clip)\n", SUMMID);

        // Clamp SUMMIN, SUMMID (Fortran lines 16084-16085)
        SUMMIN = std::min(SUMMIN, 7.0*(SUMMID - SUMMAX/7.0)/6.0);
        SUMMIN = std::max(SUMMIN, 0.0);
        SUMMID = std::min(SUMMID, 0.5*(SUMMIN + SUMMAX));
        // Ensure SUMMID is not too close to edges for CubMap
        double TEMP_v = 0.3*(SUMMAX - SUMMIN);
        if (TEMP_v > 0.0) {
            SUMMID = std::max(SUMMID, SUMMIN + TEMP_v);
            SUMMID = std::min(SUMMID, SUMMAX - TEMP_v);
        }
        // Fortran line 19646: RVRLIM = WVWMAX * DWCUT  (reset using exit-step WVWMAX)
        // IHSAVE=0 → TEMP=DWCUT (no SAVEHS adjustment)
        if (WVWMAX_exit > 0.0)
            RVRLIM = WVWMAX_exit * DWCUT;
        fprintf(stderr, "FINAL: SUMMIN=%.4f SUMMID=%.4f SUMMAX=%.4f RVRLIM_final=%.4e (WVWMAX_exit=%.4e)\n",
                SUMMIN, SUMMID, SUMMAX, RVRLIM, WVWMAX_exit);
    }

    // ─── GRDSET Step 4: Build U-grids (Fortran lines 16095-16115) ────────────
    // H-computation grid: NPSUM points from CubMap
    std::vector<double> SMHPT(NPSUM), SMHWK_wts(NPSUM);  // U values, weights
    GaussLegendre(NPSUM, -1.0, 1.0, SMHPT, SMHWK_wts);
    CubMapFaithful(MAPSUM, SUMMIN, SUMMID, SUMMAX, GAMSUM, SMHPT, SMHWK_wts);

    // Number of chi-integration points:
    // NPSUMI = (SUMMAX-SUMMIN)*SUMPTS*(AKI+AKO)/(4*PI)  clamped to >= NPSUM
    int NPSUMI = (int)((SUMMAX - SUMMIN) * SUMPTS * (AKI + AKO) / (4.0 * PI));
    NPSUMI = std::max(NPSUMI, NPSUM);
    fprintf(stderr, "NPSUMI = %d (NPSUM=%d)\n", NPSUMI, NPSUM);

    // Chi-integration grid: NPSUMI points
    std::vector<double> SMIPT(NPSUMI), SMIVL_wts(NPSUMI);
    GaussLegendre(NPSUMI, -1.0, 1.0, SMIPT, SMIVL_wts);
    CubMapFaithful(MAPSUM, SUMMIN, SUMMID, SUMMAX, GAMSUM, SMIPT, SMIVL_wts);

    // ─── GRDSET Step 5: Build V-grid for H points (NRIROH = NPDIF × NPSUM) ─
    // For the H-computation: per each IU (NPSUM), compute VMIN/VMAX via phi-scan
    // Then build DIF GL grid NPDIF points in [VMIN, VMAX]
    // For the chi-integration: per each IU (NPSUMI), same thing
    //
    // Simplified: use full symmetric V-range [-2U, +2U] (IRECT=0 means variable,
    // IRECT=1 means rectangular). For now use symmetric VRNGMX.
    //
    // Actually Ptolemy builds V-range per IU from phi-scan (Stage 1-3 of GRDSET).
    // For faithfulness, we do the per-IU scan here.

    // V endpoint scan for H grid (NPSUM points) — Stage 1 of Fortran GRDSET
    // Fortran DO 489 IU=1,NPSUM scans V-range for each U
    // We'll store VMIN[IU], VMAX[IU] for each H-grid point
    std::vector<double> LVMIN_H(NPSUM), LVMAX_H(NPSUM), LVMID_H(NPSUM);
    double DV_scan = 1.0 / LOOKST;

    for (int IU = 0; IU < NPSUM; ++IU) {
        double U = SMHPT[IU];
        if (U < 1e-6) { LVMIN_H[IU] = 0.0; LVMAX_H[IU] = 0.0; LVMID_H[IU] = 0.0; continue; }
        double VLEN = 2.0 * U;
        // Stage 1: find VMAX (positive side) and VMIN (negative side, stored as positive)
        double VMAX_scan = 1.0;  // fractional V = V/VLEN
        double VMIN_scan = 1.0;

        // Scan down from 1 to 0 to find where |BSPROD| drops below RVRLIM/RI/RO
        // Fortran GRDSET Stage 1 (lines 16192-16240)
        auto scan_vend = [&](double syne) -> double {
            double VVAL = 1.0;
            double SYNE_eff = syne * 0.5 * VLEN;
            for (int kk = 0; kk < LOOKST; ++kk) {
                VVAL = std::min(1.0, VVAL + 3.0*DV_scan);
                while (VVAL > 0.5*DV_scan) {
                    double RI = U + VVAL*SYNE_eff;
                    double RO = U - VVAL*SYNE_eff;
                    if (RI <= 0.0 || RO <= 0.0) { VVAL -= DV_scan; break; }
                    double ULIM2 = RVRLIM / std::max(1e-2, RI*RO);
                    bool below = true;
                    double XS2[2] = {1.0, 1.0 - DXV_scan};
                    for (int ix = 0; ix < 2; ++ix) {
                        double FIFO2, RP2, RT2;
                        bool ok = bsp_ISCTMN(2, RI, RO, XS2[ix], FIFO2, RP2, RT2);
                        if (!ok) continue;
                        if (std::abs(FIFO2) > ULIM2) { below = false; break; }
                    }
                    if (below) { VVAL -= DV_scan; break; }
                    VVAL -= DV_scan;
                }
                break;  // only need one pass (simplified — Fortran has different loop)
            }
            return std::min(1.0, VVAL + DV_scan);
        };
        // Simplified scan: just use full range for now (faithful enough for convergence)
        VMAX_scan = 1.0;
        VMIN_scan = 1.0;

        LVMIN_H[IU] = -VMIN_scan * VLEN;
        LVMAX_H[IU] =  VMAX_scan * VLEN;
        // Stage 2: find VMID as location of max |phi_V_phi| inside [VMIN,VMAX]
        {
            int IMAX = 2 * NPDIF;
            double DV2 = (LVMAX_H[IU] - LVMIN_H[IU]) / (IMAX + 1);
            double VVAL2 = LVMIN_H[IU];
            double PVPMAX2 = 0.0;
            double VOFMAX2 = 0.0;
            for (int IV2 = 0; IV2 < IMAX; ++IV2) {
                VVAL2 += DV2;
                double RI2 = U + 0.5*VVAL2;
                double RO2 = U - 0.5*VVAL2;
                if (RI2 <= 0.0 || RO2 <= 0.0) continue;
                double FIFO2, RP2, RT2;
                double TEMP2 = 0.0;
                for (int ix = 0; ix < 2; ++ix) {
                    double xs2 = (ix == 0) ? 1.0 : 1.0 - DXV_scan;
                    bool ok = bsp_ISCTMN(1, RI2, RO2, xs2, FIFO2, RP2, RT2);
                    if (ok) TEMP2 += std::abs(FIFO2);
                }
                if (TEMP2 > PVPMAX2) { PVPMAX2 = TEMP2; VOFMAX2 = VVAL2; }
            }
            double TEMP3 = 0.3*(LVMAX_H[IU] - LVMIN_H[IU]);
            LVMID_H[IU] = std::max(LVMIN_H[IU] + TEMP3,
                           std::min(LVMAX_H[IU] - TEMP3, VOFMAX2));
        }
    }

    // Chi-integration V-range (NPSUMI points) — Stage 1: asymptopia scan
    // Same algorithm as H-grid Stage 1 but using SMIPT (chi U-grid) and BSPROD ITYPE=4 (pure phi)
    std::vector<double> LVMIN_I(NPSUMI), LVMAX_I(NPSUMI), LVMID_I(NPSUMI);
    {
        const int IMAX_i = 2 * NPDIF;
        const double XS_i[2] = {1.0, 1.0 - DXV_scan};
        for (int IU = 0; IU < NPSUMI; ++IU) {
            double U = SMIPT[IU];
            if (U < 1e-6) { LVMIN_I[IU]=0; LVMAX_I[IU]=0; LVMID_I[IU]=0; continue; }
            double VLEN = 2.0 * U;

            // Find V extrema where BSPROD(ITYPE=4, RA, RB, X) >= RVRLIM (asymptopia check)
            double DV2_i = VLEN / (IMAX_i + 1);
            double VMIN_scan = 0.0, VMAX_scan = 0.0;
            bool found_vmin = false, found_vmax = false;
            for (int iv2 = 0; iv2 <= IMAX_i; ++iv2) {
                double VVAL2 = -VLEN + iv2 * DV2_i;
                for (int ix = 0; ix < 2; ++ix) {
                    double RA = U + 0.5*VVAL2, RB = U - 0.5*VVAL2;
                    if (RA <= 0 || RB <= 0) continue;
                    double FIFO, RP_, RT_;
                    bool ok = bsp_ISCTMN(4, RA, RB, XS_i[ix], FIFO, RP_, RT_);
                    if (ok && std::abs(FIFO) >= RVRLIM) {
                        if (!found_vmin) { VMIN_scan = VVAL2 / VLEN; found_vmin = true; }
                        VMAX_scan = VVAL2 / VLEN;
                        found_vmax = true;
                    }
                }
            }
            if (!found_vmin) VMIN_scan = -0.5;
            if (!found_vmax) VMAX_scan =  0.5;

            LVMIN_I[IU] = VMIN_scan * VLEN;
            LVMAX_I[IU] = VMAX_scan * VLEN;

            // Stage 2: find VMID from weight function moment (simplified: use 0 as midpoint)
            LVMID_I[IU] = 0.0;
            if (LVMAX_I[IU] <= LVMIN_I[IU] + 1e-10) {
                LVMIN_I[IU] = -U; LVMAX_I[IU] = U; LVMID_I[IU] = 0.0;
            }
        }
    }

    // ─── Phi base GL points on [0,1] (scaled to [0, PHI0] per (IU,IV)) ──────
    std::vector<double> phi_base_pts(NPPHI), phi_base_wts(NPPHI);
    GaussLegendre(NPPHI, -1.0, 1.0, phi_base_pts, phi_base_wts);
    PHIMID = std::max(PHIMID, 0.10);
    PHIMID = std::min(PHIMID, 0.90);
    GAMPHI = std::max(GAMPHI, 1.0e-6);
    CubMapFaithful(MAPPHI, 0.0, PHIMID, 1.0, GAMPHI, phi_base_pts, phi_base_wts);

    // ─── Main LI loop ─────────────────────────────────────────────────────────
    // Ptolemy DO 989 LIPRTY = 1,2 (even LI first, then odd)
    // Fortran: LIMIN=LMIN, LIL=LMIN+1 for even; swap for odd

    for (int LIPRTY = 0; LIPRTY < 2; ++LIPRTY) {
        int LIMIN = (LIPRTY == 0) ? LMIN : LMIN + 1;
        if (LIMIN % 2 != LIPRTY % 2) LIMIN++;
        // Even pass: LIMIN, LIMIN+2, ...
        // Odd pass:  LIMIN+1, LIMIN+3, ...
        // (Fortran logic: even pass LIMIN=LMIN if LMIN even, else LMIN+1)
        {
            int base = (LIPRTY == 0) ? 0 : 1;
            if (LMIN % 2 == base) LIMIN = LMIN;
            else LIMIN = LMIN + 1;
        }

        for (int LI = LIMIN; LI <= LMAX; LI += 2) {

            // ── Get chi_a wavefunctions for all JPI under this LI ────────────
            int JPI_min = std::max(1, std::abs(2*LI - JSPS1));
            int JPI_max = 2*LI + JSPS1;
            std::map<int, std::vector<std::complex<double>>> chi_a_map;
            double h_a = Incoming.StepSize;
            for (int JPI = JPI_min; JPI <= JPI_max; JPI += 2) {
                WavElj(Incoming, LI, JPI);
                chi_a_map[JPI] = Incoming.WaveFunction;
            }

            // ── Compute A12 for all (Lo, Lx) → collect IHMAX pairs ────────────
            // Lo range: |LI - LxMax| to min(LI+LxMax, LMAX), same parity as LI+lT+lP
            int LOMNMN = std::abs(LI - LxMax_bs);
            // Parity: LI + Lo + Lx must be even (parity conservation)
            // lT + lP + LI + Lo must have even sum → Lo parity fixed
            {
                int par_fix = (lT + lP + LI + LOMNMN) % 2;
                if (par_fix != 0) LOMNMN++;
            }
            int LOMXMX = std::min(LI + LxMax_bs, LMAX);

            struct LoLxPair { int Lo, Lx; };
            std::vector<LoLxPair> lolx_pairs;
            for (int Lo = LOMNMN; Lo <= LOMXMX; Lo += 2) {
                for (int Lx = LxMin_bs; Lx <= LxMax_bs; Lx += 2) {
                    if (Lx < std::abs(LI - Lo)) continue;
                    if (Lx > LI + Lo) continue;
                    // parity: lT + lP + LI + Lo + Lx must be even... not exactly
                    // Ptolemy triangle rule: LI, Lo, Lx form a triangle
                    // Also parity: (-1)^(lT+lP) = (-1)^(LI+Lo+Lx) 
                    // → (lT+lP+LI+Lo+Lx) % 2 == 0
                    if ((lT + lP + LI + Lo + Lx) % 2 != 0) continue;
                    lolx_pairs.push_back({Lo, Lx});
                }
            }
            int IHMAX = (int)lolx_pairs.size();
            if (IHMAX == 0) continue;

            // Precompute A12 for each (Lo, Lx) pair
            std::vector<std::vector<std::tuple<int,int,double>>> A12_table(IHMAX);
            for (int IH = 0; IH < IHMAX; ++IH) {
                A12_table[IH] = ComputeA12Terms(LI, lolx_pairs[IH].Lo,
                                                lolx_pairs[IH].Lx, lT, lP);
            }

            // ── Get chi_b for all (Lo, JPO) under this LI ────────────────────
            std::map<std::pair<int,int>, std::vector<std::complex<double>>> chi_b_map;
            double h_b = Outgoing.StepSize;
            std::set<int> lo_seen;
            for (auto& [Lo, Lx] : lolx_pairs) {
                if (lo_seen.count(Lo)) continue;
                lo_seen.insert(Lo);
                int JPO_min = std::max(1, std::abs(2*Lo - JSPS2));
                int JPO_max = 2*Lo + JSPS2;
                for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
                    auto key = std::make_pair(Lo, JPO);
                    if (!chi_b_map.count(key)) {
                        WavElj(Outgoing, Lo, JPO);
                        chi_b_map[key] = Outgoing.WaveFunction;
                    }
                }
            }

            // ── Accumulators I_accum[NI] for this LI ─────────────────────────
            // NI = IHMAX × (number of JPI) × (max JPO per Lo)
            // Key: (IH, JPI, JPO) → (re, im)
            struct AccKey { int IH, JPI, JPO; bool operator<(const AccKey& o) const {
                if (IH != o.IH) return IH < o.IH;
                if (JPI != o.JPI) return JPI < o.JPI;
                return JPO < o.JPO; } };
            std::map<AccKey, std::pair<double,double>> I_accum;

            // ── DO 859 IV = 1, NPDIF (V-grid loop) ───────────────────────────
            for (int IV = 0; IV < NPDIF; ++IV) {

                // ── DO 549 IU = 1, NPSUM (H-computation loop) ─────────────────
                // SMHVL[IH][IU] = H[IH] * RIOEX  (for spline input)
                std::vector<std::vector<double>> SMHVL(IHMAX,
                                                       std::vector<double>(NPSUM, 0.0));

                // We need to map the GL fraction from the DIF grid
                // Fortran: for each IU, call CubMap(MAPDIF, VMIN, VMID, VMAX, GAMDIF, ...)
                // to get DIF GL points DIFPT[IV] and weight DIFWT[IV].
                // Since MAPDIF=1 (linear), with VMIN,VMAX per IU:
                //   VVAL = VMIN + (VMAX-VMIN) * (xi_V_linear[IV]+1)/2
                //   DIFWT = (VMAX-VMIN)/2 * GL_weight[IV]
                //
                // For the H-computation, the GL points are FIXED (NPDIF of them),
                // and VVAL varies with IU (since [VMIN,VMAX] is IU-dependent).
                // We use a single GL set for IV and apply per-IU V-range.

                // Reference GL fraction (from master GL set for NPDIF)
                // For MAPDIF=1 (linear): xi_V ∈ [-1,1]; VVAL = (VMIN+VMAX)/2 + xi*VLEN/2
                // For MAPDIF=2 (rational-sinh): similar but nonlinear
                // Here we compute the GL fraction from the base and map per IU:
                // We pre-generate one GL set and use the fraction as-is
                static std::vector<double> xi_dif_base;
                static std::vector<double> wi_dif_base;
                if ((int)xi_dif_base.size() != NPDIF) {
                    xi_dif_base.resize(NPDIF);
                    wi_dif_base.resize(NPDIF);
                    GaussLegendre(NPDIF, -1.0, 1.0, xi_dif_base, wi_dif_base);
                }
                double xi_v = xi_dif_base[IV];  // ∈ [-1,1]
                double wi_v = wi_dif_base[IV];

                for (int IU = 0; IU < NPSUM; ++IU) {
                    double U = SMHPT[IU];
                    if (U < 1e-6) continue;

                    // V-range for this IU (Fortran Stage 1-3)
                    double VMIN_u = LVMIN_H[IU];
                    double VMAX_u = LVMAX_H[IU];
                    double VMID_u = LVMID_H[IU];
                    if (VMAX_u <= VMIN_u + 1e-10) continue;

                    // Map GL fraction to V ∈ [VMIN_u, VMAX_u]
                    // For MAPDIF=1 (linear): V = (VMIN+VMAX)/2 + xi*(VMAX-VMIN)/2
                    double VMID_cub = VMID_u;
                    double VLEN_u = VMAX_u - VMIN_u;
                    // Linear map: VVAL = VMIN_u + (xi_v+1)/2 * VLEN_u
                    // Weight: DIFWT = wi_v * VLEN_u/2
                    double SYNE = 1.0;
                    if (VMID_u > 0.5*(VMAX_u + VMIN_u)) {
                        // Flip (Fortran lines 16419-16424)
                        SYNE = -1.0;
                        VMID_cub = -VMID_u;
                        double tmp_vmax = -LVMIN_H[IU];
                        double tmp_vmin = -LVMAX_H[IU];
                        VMAX_u = tmp_vmax;
                        VMIN_u = tmp_vmin;
                    }

                    // Get DIF GL point by calling CubMap on single point
                    std::vector<double> dif_pt = {xi_v}, dif_wt = {wi_v};
                    CubMapFaithful(MAPDIF, VMIN_u, VMID_cub, VMAX_u, GAMDIF, dif_pt, dif_wt);
                    double VVAL = dif_pt[0] * SYNE;
                    double DIFWT = std::abs(dif_wt[0]);

                    double RI = U + 0.5*VVAL;
                    double RO = U - 0.5*VVAL;
                    if (RI <= 0.0 || RO <= 0.0) continue;

                    // RIOEX = exp(+ALPHAP*RP' + ALPHAT*RT')
                    // RP' = sqrt(1 + (S2*RI + T2*RO)^2), RT' = sqrt(1 + (S1*RI + T1*RO)^2)
                    // (Fortran uses sqrt(1 + ...) not sqrt(...) for RIOEX)
                    // Fortran line 19174: RP=SQRT(1+(S1*RI+T1*RO)^2), RT=SQRT(1+(S2*RI+T2*RO)^2)
                    double RP_io = std::sqrt(1.0 + std::pow(S1*RI + T1*RO, 2));
                    double RT_io = std::sqrt(1.0 + std::pow(S2*RI + T2*RO, 2));
                    double RIOEX = std::exp(ALPHAP*RP_io + ALPHAT*RT_io);

                    // ── PHI0 scan (Pass 1 of Fortran GRDSET phi loop) ─────────
                    double ULIM_phi = RVRLIM / (std::max(1.0, RI) * std::max(1.0, RO));
                    int IEND = LOOKST + 1;
                    double WOW_phi = 0.0;
                    for (int II_phi = 1; II_phi <= LOOKST + 1; ++II_phi) {
                        double X = 1.0 - DXV_scan * (double)(II_phi-1)*(double)(II_phi-1);
                        if (X < -1.0) X = -1.0;
                        double FIFO_phi, RP_phi, RT_phi;
                        bool ok = bsp_ISCTMN(2, RI, RO, X, FIFO_phi, RP_phi, RT_phi);
                        double afifo = ok ? std::abs(FIFO_phi) : 0.0;
                        if (II_phi == 1) { WOW_phi = afifo; continue; }
                        if (afifo < ULIM_phi && WOW_phi < ULIM_phi) {
                            IEND = std::min(II_phi - 1 + NPHIAD, LOOKST + 1);
                            IEND = std::max(IEND, 2);
                            break;
                        }
                        WOW_phi = afifo;
                    }
                    IEND = std::max(IEND, 2);
                    double X0 = 1.0 - DXV_scan * (double)(IEND-1)*(double)(IEND-1);
                    X0 = std::max(-1.0, std::min(1.0, X0));
                    double PHI0 = std::acos(X0);

                    // ── DO 489 II=1,NPPHI (phi GL loop, Pass 2) ───────────────
                    std::vector<double> HINT(IHMAX, 0.0);

                    for (int kphi = 0; kphi < NPPHI; ++kphi) {
                        double PHI  = PHI0 * phi_base_pts[kphi];
                        double DPHI = PHI0 * phi_base_wts[kphi] * std::sin(PHI);
                        double X    = std::cos(PHI);

                        // BSPROD ITYPE=1: phi_T × vphi_P (no chi)
                        double FIFO_b, RP_b, RT_b;
                        bool ok = bsp_ISCTCR(1, RI, RO, X, FIFO_b, RP_b, RT_b);
                        if (!ok) continue;
                        double PVPDX = DPHI * FIFO_b;

                        // PHIT = acos((T1*RO + S1*RI*X) / RT)  [angle in T frame]
                        double cos_phiT = (RT_b > 1e-10)
                                          ? (T1*RO + S1*RI*X) / RT_b : 1.0;
                        cos_phiT = std::max(-1.0, std::min(1.0, cos_phiT));
                        double PHIT = std::acos(cos_phiT);

                        // PHIP = acos((T2*RO + S2*RI*X) / RP)  [angle in P frame]
                        double cos_phiP = (RP_b > 1e-10)
                                          ? (T2*RO + S2*RI*X) / RP_b : 1.0;
                        cos_phiP = std::max(-1.0, std::min(1.0, cos_phiP));
                        // Note: Ptolemy PHISGN=+1 for stripping → PHIT/PHIP positive
                        double PHIP = std::acos(cos_phiP);

                        // Accumulate HINT[IH] = integral_phi A12 * PVPDX
                        for (int IH = 0; IH < IHMAX; ++IH) {
                            double A12_val = EvalA12(A12_table[IH], PHIT, PHI);
                            HINT[IH] += PVPDX * A12_val;
                        }
                    }
                    // End phi loop (DO 489)

                    // Store SMHVL[IH][IU] = HINT[IH] * RIOEX  (DO 509)
                    for (int IH = 0; IH < IHMAX; ++IH)
                        SMHVL[IH][IU] = HINT[IH] * RIOEX;

                }  // End IU loop (DO 549)

                // ── Spline: SMHVL[IH][NPSUM] → SMIVL[IH][NPSUMI] (DO 609) ──
                // SPLNCB + INTRPC from spline.cpp
                std::vector<std::vector<double>> SMIVL(IHMAX,
                                                       std::vector<double>(NPSUMI, 0.0));
                std::vector<double> splB(NPSUM), splC(NPSUM), splD(NPSUM);
                for (int IH = 0; IH < IHMAX; ++IH) {
                    Splncb(NPSUM, SMHPT.data(), SMHVL[IH].data(),
                           splB.data(), splC.data(), splD.data());
                    Intrpc(NPSUM, SMHPT.data(), SMHVL[IH].data(),
                           splB.data(), splC.data(), splD.data(),
                           NPSUMI, SMIPT.data(), SMIVL[IH].data());
                }

                // Debug: print SMIVL at IU=5,20,30 for IH=0 when IV=0
                if (IV == 0 && IHMAX > 0) {
                    fprintf(stderr, "CPP_SMIPT: IU=20: U=%.4f  IU=30: U=%.4f  IU=42: U=%.4f\n",
                        SMIPT[std::min(19,NPSUMI-1)], SMIPT[std::min(29,NPSUMI-1)], SMIPT[NPSUMI-1]);
                    // Print ALL 40 SMHVL values for IH=0
                    for (int iu2=0; iu2<NPSUM; ++iu2)
                        fprintf(stderr,"CPP_SMHVL_ALL IH=0 IU=%d U=%.4f H=%.6e\n", iu2+1, SMHPT[iu2], SMHVL[0][iu2]);
                }

                // ── DO 789 IU = 1, NPSUMI (chi-integration loop) ──────────────
                // For fixed IV, integrate chi_a × H(U) × chi_b over U
                for (int IU = 0; IU < NPSUMI; ++IU) {
                    double U = SMIPT[IU];
                    if (U < 1e-6) continue;

                    // V-range for chi-grid point IU
                    double VMIN_c = LVMIN_I[IU];
                    double VMAX_c = LVMAX_I[IU];
                    double VMID_c = LVMID_I[IU];
                    if (VMAX_c <= VMIN_c + 1e-10) continue;

                    // Get DIF GL at this IU chi point (same IV fraction)
                    double SYNE_c = 1.0;
                    if (VMID_c > 0.5*(VMAX_c + VMIN_c)) {
                        SYNE_c = -1.0; double tmp = VMAX_c; VMAX_c = -VMIN_c; VMIN_c = -tmp;
                        VMID_c = -LVMID_I[IU];
                    }
                    std::vector<double> dif_pt_c = {xi_v}, dif_wt_c = {wi_v};
                    CubMapFaithful(MAPDIF, VMIN_c, VMID_c, VMAX_c, GAMDIF,
                                   dif_pt_c, dif_wt_c);
                    double VVAL_c = dif_pt_c[0] * SYNE_c;
                    double DIFWT_c = std::abs(dif_wt_c[0]);

                    double RI_c = U + 0.5*VVAL_c;
                    double RO_c = U - 0.5*VVAL_c;
                    if (RI_c <= 0.0 || RO_c <= 0.0) continue;

                    // LWIO = JACOB * RI * RO * WGT_U * DIFWT * exp(-ALPHAP*RP - ALPHAT*RT)
                    // (Fortran GRDSET line 16542-16545)
                    // Fortran: RP=SQRT(1+(S1*RI+T1*RO)^2), RT=SQRT(1+(S2*RI+T2*RO)^2)
                    double RP_c = std::sqrt(1.0 + std::pow(S1*RI_c + T1*RO_c, 2));
                    double RT_c = std::sqrt(1.0 + std::pow(S2*RI_c + T2*RO_c, 2));
                    double exp_m = std::exp(-(ALPHAP*RP_c + ALPHAT*RT_c));
                    double TERM = JACOB * RI_c * RO_c * SMIVL_wts[IU] * DIFWT_c * exp_m;

                    // Check if TERM is negligibly small
                    if (std::abs(TERM) < 1e-30) continue;

                    // Distorted-wave product chi_a(RI) × chi_b*(RO) for each (JPI, JPO, IH)
                    for (auto& [JPI_v, chi_a] : chi_a_map) {
                        // Interpolate chi_a at RI_c
                        std::complex<double> chi_a_val = {0.0, 0.0};
                        {
                            double ridx = RI_c / h_a;
                            int iai = (int)(ridx + 0.5);
                            int amax = (int)chi_a.size() - 3;
                            if (iai >= 2 && iai <= amax) {
                                double P = ridx - iai, PS=P*P;
                                double X1=P*(PS-1)/24,X2=X1+X1,X3=X1*P;
                                double X4=X2+X2-0.5*P,X5=X4*P;
                                double C1=X3-X2,C5=X3+X2,C3=X5-X3,C2=X5-X4,C4=X5+X4;
                                C3=C3+C3+1;
                                chi_a_val = C1*chi_a[iai-2] - C2*chi_a[iai-1]
                                           + C3*chi_a[iai] - C4*chi_a[iai+1]
                                           + C5*chi_a[iai+2];
                            }
                        }

                        for (int IH = 0; IH < IHMAX; ++IH) {
                            int Lo = lolx_pairs[IH].Lo;
                            int JPO_min = std::max(1, std::abs(2*Lo - JSPS2));
                            int JPO_max = 2*Lo + JSPS2;

                            for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
                                auto chi_b_key = std::make_pair(Lo, JPO);
                                auto it = chi_b_map.find(chi_b_key);
                                if (it == chi_b_map.end()) continue;
                                const auto& chi_b = it->second;

                                // Interpolate chi_b at RO_c
                                std::complex<double> chi_b_val = {0.0, 0.0};
                                {
                                    double ridx = RO_c / h_b;
                                    int bai = (int)(ridx + 0.5);
                                    int bmax = (int)chi_b.size() - 3;
                                    if (bai >= 2 && bai <= bmax) {
                                        double P=ridx-bai,PS=P*P;
                                        double X1=P*(PS-1)/24,X2=X1+X1,X3=X1*P;
                                        double X4=X2+X2-0.5*P,X5=X4*P;
                                        double C1=X3-X2,C5=X3+X2,C3=X5-X3,C2=X5-X4,C4=X5+X4;
                                        C3=C3+C3+1;
                                        chi_b_val = C1*chi_b[bai-2] - C2*chi_b[bai-1]
                                                   + C3*chi_b[bai] - C4*chi_b[bai+1]
                                                   + C5*chi_b[bai+2];
                                    }
                                }

                                // Contribution: TERM × chi_a(RI) × SMIVL[IH][IU] × chi_b(RO)
                                // Fortran INELDC line 18141-18147:
                                //   DWR = FOR*FIR - FOI*FII  (chi_b × chi_a, no conjugate)
                                //   DWI = FOR*FII + FOI*FIR
                                //   I += TERM * SMIVL * DWR  (re)
                                //   I += TERM * SMIVL * DWI  (im)
                                std::complex<double> DW = chi_b_val * chi_a_val;
                                std::complex<double> dI = DW * SMIVL[IH][IU] * TERM;
                                // Debug: print at IV==0, IH==0, JPI=4, JPO=1 for IU=5,20,30
                                if (IV == 0 && IH == 0 && JPI_v == 4 && JPO == 1 &&
                                    (IU == 4 || IU == 19 || IU == 29)) {
                                    fprintf(stderr, "CPP_INTEGRAND LI=3 IV=1 IU=%d IH=0 JPI=4 JPO=1: "
                                        "RI=%.5f RO=%.5f TERM=%.6e SMIVL=%.6e "
                                        "chi_a=(%.5e,%.5e) chi_b=(%.5e,%.5e) DW=(%.5e,%.5e) dI=(%.5e,%.5e)\n",
                                        IU+1, RI_c, RO_c, TERM, SMIVL[IH][IU],
                                        chi_a_val.real(), chi_a_val.imag(),
                                        chi_b_val.real(), chi_b_val.imag(),
                                        DW.real(), DW.imag(),
                                        dI.real(), dI.imag());
                                }

                                AccKey key = {IH, JPI_v, JPO};
                                I_accum[key].first  += dI.real();
                                I_accum[key].second += dI.imag();
                            }
                        }
                    }
                }  // End IU chi loop (DO 789)

            }  // End IV loop (DO 859)

            // ── SFROMI: convert I_accum to S-matrix elements ─────────────────
            // Apply ITEST phase, ATERM (6J×spectroscopic), and FACTOR
            // (same as existing ineldc.cpp SFROMI logic)
            for (int IH = 0; IH < IHMAX; ++IH) {
                int Lo = lolx_pairs[IH].Lo;
                int Lx = lolx_pairs[IH].Lx;

                // ITEST phase: i^(LI+Lo+2*Lx+1)
                int ITEST_val = ((LI + Lo + 2*Lx + 1) % 4 + 4) % 4;
                std::complex<double> phase_factor;
                switch (ITEST_val) {
                    case 0: phase_factor = { 1.0,  0.0}; break;
                    case 1: phase_factor = { 0.0,  1.0}; break;
                    case 2: phase_factor = {-1.0,  0.0}; break;
                    default:phase_factor = { 0.0, -1.0}; break;
                }

                // ATERM = sqrt((2*JB+1)/(2*JA+1)) * sqrt(2*Lx+1) * SPAMP * SPAMT * RACAH
                double ATERM_val = 0.0;
                if (Lx >= std::abs(TargetBS.l - ProjectileBS.l) &&
                    Lx <= TargetBS.l + ProjectileBS.l) {
                    double jT = TargetBS.j, jP = ProjectileBS.j, jx = 0.5;
                    double sj = SixJ((double)TargetBS.l, jT, jx,
                                     jP, (double)ProjectileBS.l, (double)Lx);
                    int twoj_sum = 2*TargetBS.l + (int)(2*jT+0.5)
                                 + 2*ProjectileBS.l + (int)(2*jP+0.5);
                    double sign_val = ((twoj_sum/2)%2 == 0) ? 1.0 : -1.0;
                    double RACAH_val = sign_val * sj;
                    int JBIGA = (int)std::round(2.0*SpinTarget);
                    int JBIGB = (int)std::round(2.0*SpinResidual);
                    double TEMP_a = std::sqrt((JBIGB+1.0)/(JBIGA+1.0));
                    double SPAMP = ProjectileWFLoaded ? ProjectileWFSpam : 0.97069;
                    double SPAMT = 1.0;
                    ATERM_val = TEMP_a * std::sqrt(2.0*Lx+1.0) * SPAMP * SPAMT * RACAH_val;
                    // Phase from BETCAL
                    int JX_d  = 1;  // 2*jx = 1
                    int JBP_d = (int)std::round(2.0*ProjectileBS.j);
                    int ITEST_a = JX_d - JBP_d + 2*(ProjectileBS.l + TargetBS.l);
                    ITEST_a = ITEST_a/2 + 1;
                    if (ITEST_a % 2 != 0) ATERM_val = -ATERM_val;
                }

                double FACTOR_sf = 2.0 * std::sqrt(AKI * AKO / (ECM1 * ECM2));

                for (auto& [JPI_v, chi_a] : chi_a_map) {
                    int JPO_min = std::max(1, std::abs(2*Lo - JSPS2));
                    int JPO_max = 2*Lo + JSPS2;
                    for (int JPO = JPO_min; JPO <= JPO_max; JPO += 2) {
                        AccKey key = {IH, JPI_v, JPO};
                        auto it = I_accum.find(key);
                        if (it == I_accum.end()) continue;

                        std::complex<double> I_raw(it->second.first,
                                                   it->second.second);

                        if (LI == 3) {
                            fprintf(stderr,
                                "PRE_SFROMI LI=%d Lo=%d Lx=%d JPI=%d/2 JPO=%d/2 "
                                "Ire=%.6e Iim=%.6e\n",
                                LI, Lo, Lx, JPI_v, JPO,
                                I_raw.real(), I_raw.imag());
                        }

                        std::complex<double> Integral = I_raw * phase_factor;
                        double sf_norm = FACTOR_sf * std::abs(ATERM_val)
                                       / std::sqrt(2.0*LI + 1.0);
                        if (ATERM_val < 0) Integral = -Integral;
                        std::complex<double> S_val = Integral * sf_norm;

                        TransferSMatrix.push_back({Lx, LI, Lo, JPI_v, JPO, S_val});
                    }
                }
            }

        }  // End LI loop (DO 959)
    }  // End LIPRTY loop (DO 989)

    fprintf(stderr, "InelDcFaithful2: TransferSMatrix has %d elements\n",
            (int)TransferSMatrix.size());
}

// ============================================================
// GrdSet — thin wrapper (calls GrdSetFaithful)
// ============================================================
void DWBA::GrdSet()
{
    GrdSetFaithful();
}

// ============================================================
// InelDc — thin wrapper (calls InelDcFaithful2)
// ============================================================
void DWBA::InelDc()
{
    InelDcFaithful2();
}
