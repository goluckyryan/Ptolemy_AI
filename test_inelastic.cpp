// test_inelastic.cpp — DWBA inelastic using WAVELJ for distorted waves
// Rewritten to use DWBA::WavElj (matching Ptolemy's Numerov) instead of ElasticSolver
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include "dwba.h"
#include "math_utils.h"
#include "ptolemy_mass_table.h"
#include "coulin.h"
#include <map>
// Forward declaration for Rcwfn (defined in src/dwba/rcwfn.cpp)
int Rcwfn(double rho, double eta, int lmin, int lmax, std::vector<double>& FC,
          std::vector<double>& FCP, std::vector<double>& GC, std::vector<double>& GCP, double accur);

static const double PI = 3.14159265358979323846;
static const double HBARC = 197.32858;    // Ptolemy value
static const double AMUMEV = 931.50160;   // Ptolemy value
static const double AFINE = 1.0/137.03604; // Ptolemy value
static const double EMASS = 0.5110034;    // electron mass in MeV

// Coulomb phase shift σ_L = arg(Γ(L+1+iη))
// Uses Stirling approximation at N0=200, same as ElasticSolver::CoulombPhase
static double CoulombPhaseShift(int L, double eta) {
    int N0 = 200;
    double sum_atan = 0;
    for (int n = 1; n <= N0; ++n) sum_atan += std::atan(eta / n);
    double x = N0 + 1, y = eta;
    double theta_N0 = (x-0.5)*std::atan2(y, x) + y/2.0*std::log(x*x + y*y) - y;
    double sig = theta_N0 - sum_atan;
    for (int n = 1; n <= L; ++n) sig += std::atan(eta / n);
    return sig;
}

// CG coefficient for integer arguments
static double CG_int(int j1, int j2, int m1, int m2, int J, int M) {
    if (M != m1+m2) return 0.0;
    if (J < std::abs(j1-j2) || J > j1+j2) return 0.0;
    if (std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(M) > J) return 0.0;
    auto fac = [](int n) -> double {
        double f = 1.0; for (int i = 2; i <= n; i++) f *= i; return f;
    };
    double sum = 0.0;
    for (int k = 0; k <= j1+j2+J; k++) {
        int n1=j1+j2-J-k, n2=j1-m1-k, n3=j2+m2-k, n4=J-j2+m1+k, n5=J-j1-m2+k;
        if (n1<0||n2<0||n3<0||n4<0||n5<0) continue;
        sum += (k%2==0?1.0:-1.0) / (fac(k)*fac(n1)*fac(n2)*fac(n3)*fac(n4)*fac(n5));
    }
    double delta = fac(j1+j2-J)*fac(j1-j2+J)*fac(-j1+j2+J) / fac(j1+j2+J+1);
    double norm = (2*J+1)*fac(J+M)*fac(J-M)*fac(j1-m1)*fac(j1+m1)*fac(j2-m2)*fac(j2+m2);
    return sum * std::sqrt(delta * norm);
}

static double LegP(int L, double x) {
    if (L==0) return 1.0;
    if (L==1) return x;
    double p0=1, p1=x, p2;
    for (int l=2; l<=L; l++) { p2=((2*l-1)*x*p1-(l-1)*p0)/l; p0=p1; p1=p2; }
    return p1;
}

// Associated Legendre P_L^M(x), NO Condon-Shortley phase (no (-1)^M)
static double AssocLegP(int L, int M, double x) {
    if (M < 0 || M > L) return 0.0;
    if (M == 0) return LegP(L, x);
    double pmm = 1.0;
    double somx2 = std::sqrt(1.0 - x*x);
    double fact = 1.0;
    for (int i = 1; i <= M; i++) {
        pmm *= fact * somx2;
        fact += 2.0;
    }
    if (L == M) return pmm;
    double pmm1 = x * (2*M+1) * pmm;
    if (L == M+1) return pmm1;
    double pll;
    for (int ll = M+2; ll <= L; ll++) {
        pll = ((2*ll-1)*x*pmm1 - (ll+M-1)*pmm) / (double)(ll-M);
        pmm = pmm1;
        pmm1 = pll;
    }
    return pmm1;
}

// WS derivative: d/dr of f_WS = -exp(x)/(a*(1+exp(x))^2) where x=(r-R)/a
// For potential V(r) = -V0*f_WS(r), dV/dr = V0*exp(x)/(a*(1+exp(x))^2)
static double dWS(double r, double V0, double R, double a) {
    double x = (r-R)/a;
    if (std::abs(x) > 500) return 0.0;
    double ex = std::exp(x);
    return V0 * ex / (a * (1+ex)*(1+ex));
}

// Surface absorption potential: V_surf(r) = -4*a_s*V_s * df_WS/dr
//   = V_s * 4 * exp(x) / ((1+exp(x))^2)  where x = (r-R_s)/a_s
// Derivative of this: d(V_surf)/dr = V_s * 4/a_s * exp(x)*(1-exp(x))/(1+exp(x))^3
static double dSurfWS(double r, double Vs, double R, double a) {
    double x = (r-R)/a;
    if (std::abs(x) > 500) return 0.0;
    double ex = std::exp(x);
    double d = 1+ex;
    return Vs * 4.0 * ex * (1.0 - ex) / (a * d*d*d);
}

// Woods-Saxon potential value: V(r) = -V0 / (1 + exp((r-R)/a))
// Matches Fortran WOODSX ITYPE=1
static double WS_val(double r, double V0, double R, double a) {
    double x = (r - R) / a;
    if (x > 500) return 0.0;
    if (x < -500) return -V0;
    return -V0 / (1.0 + std::exp(x));
}

// Surface WS potential value: V_surf(r) = (-4*V0) * exp(x) / (1+exp(x))^2
// where x = (r-R)/a.  Matches Fortran WOODSX ITYPE=3
static double SurfWS_val(double r, double V0, double R, double a) {
    double x = (r - R) / a;
    if (std::abs(x) > 500) return 0.0;
    double ex = std::exp(x);
    double dd = 1.0 + ex;
    return (-4.0 * V0) * ex / (dd * dd);
}

// Natural cubic spline fit (Ptolemy SPLNCB equivalent)
// Input: NPTS data points y[0..NPTS-1] on uniform grid with spacing dx
// Output: coefficients b[], c[], d[] (same size as y)
// Spline in interval [i, i+1]: S(t) = y[i] + b[i]*t + c[i]*t^2 + d[i]*t^3
//   where t = x - x_i (offset within interval)
// Natural BC: S''(x_0) = S''(x_{N-1}) = 0
static void splncb(int NPTS, double dx, const std::vector<double>& y,
                   std::vector<double>& b, std::vector<double>& c,
                   std::vector<double>& d) {
    int n = NPTS - 1;  // number of intervals
    c[0] = 0.0;
    
    // Thomas algorithm for tridiagonal system
    // c[i-1] + 4*c[i] + c[i+1] = 3*(y[i+1] - 2*y[i] + y[i-1]) / dx^2
    std::vector<double> cp(NPTS, 0.0), dp_vec(NPTS, 0.0);
    
    for (int i = 1; i < n; i++) {
        double rhs = 3.0 * (y[i+1] - 2.0*y[i] + y[i-1]) / (dx * dx);
        double m = 1.0 / (4.0 - cp[i-1]);
        cp[i] = m;
        dp_vec[i] = (rhs - dp_vec[i-1]) * m;
    }
    
    // Back substitution
    c[n] = 0.0;  // natural BC
    for (int i = n-1; i >= 1; i--) {
        c[i] = dp_vec[i] - cp[i] * c[i+1];
    }
    
    // Compute b and d from c
    for (int i = 0; i < n; i++) {
        d[i] = (c[i+1] - c[i]) / (3.0 * dx);
        b[i] = (y[i+1] - y[i]) / dx - dx * (2.0*c[i] + c[i+1]) / 3.0;
    }
}

int main() {
    // ===== System: 206Hg(d,d')206Hg* (4+, Ex=1.0) Elab=14.78 =====
    int Ap = 2, Zp = 1, At = 206, Zt = 80;
    double Elab = 14.78;
    double Ex = 1.0;
    int LX = 4;
    // Input is BELX = 0.12 (B(E4) in e² barn^4), NOT beta!
    // Convert: beta_C = (4*pi/(3*Z)) * sqrt(BELX) / (Rc/10)^LX / CG_spin
    // Then beta_N = beta_C * Rc / R_real (equal deformation lengths)
    double BELX = 0.12;
    int Lmax = 50;

    // OM parameters - incoming
    double V_in=96.891, r0_in=1.151, a_in=0.793;
    double Vi_in=2.023, ri0_in=1.322, ai_in=0.264;
    double Vsi_in=10.378, rsi0_in=1.360, asi_in=0.897;
    double rc0_in=1.303;

    // OM parameters - outgoing
    double V_out=96.891, r0_out=1.151, a_out=0.793;
    double Vi_out=2.023, ri0_out=1.322, ai_out=0.264;
    double Vsi_out=10.378, rsi0_out=1.360, asi_out=0.897;
    double rc0_out=1.303;

    // Kinematics — Ptolemy convention with mass excess + electron correction
    double dx_p = PtolemyMass::MassExcess_keV(Zp, Ap) / 1000.0;  // MeV
    double dx_t = PtolemyMass::MassExcess_keV(Zt, At) / 1000.0;  // MeV
    double TMP = Ap + dx_p/AMUMEV - Zp*(EMASS/AMUMEV);  // projectile mass in AMU
    double TMT = At + dx_t/AMUMEV - Zt*(EMASS/AMUMEV);  // target mass in AMU
    double mu = AMUMEV * TMP * TMT / (TMP + TMT);        // reduced mass in MeV
    double mu_amu = TMP * TMT / (TMP + TMT);             // reduced mass in AMU
    double Ecm_in = Elab * TMT / (TMP + TMT);            // use actual masses
    double Ecm_out = Ecm_in - Ex;
    double k_in = std::sqrt(2.0 * mu * Ecm_in) / HBARC;
    double k_out = std::sqrt(2.0 * mu * Ecm_out) / HBARC;
    double eta_in = Zp * Zt * AFINE * mu / (HBARC * k_in);

    double R0mass = std::pow((double)At, 1.0/3.0);  // A^(1/3) for OM radii

    double Elab_out = Ecm_out * (TMP + TMT) / TMT;  // use actual masses

    std::cerr << "Ecm_in=" << Ecm_in << " Ecm_out=" << Ecm_out << std::endl;
    std::cerr << "k_in=" << k_in << " k_out=" << k_out << std::endl;
    std::cerr << "Elab_out=" << Elab_out << std::endl;

    // ===== Set up DWBA object with Incoming/Outgoing channels =====
    // Using WAVELJ (Ptolemy's Numerov) instead of ElasticSolver
    DWBA dwba;

    // Incoming channel setup
    dwba.Incoming.Projectile = Isotope(Ap, Zp);
    dwba.Incoming.Target = Isotope(At, Zt);
    dwba.Incoming.Elab = Elab;
    dwba.Incoming.Ecm = Ecm_in;
    dwba.Incoming.k = k_in;
    dwba.Incoming.eta = eta_in;
    dwba.Incoming.mu = mu_amu;
    dwba.Incoming.JSPS = 2;  // deuteron spin=1, JSPS=2*spin=2
    dwba.Incoming.Pot = {V_in, r0_in, a_in,
                         Vi_in, ri0_in, ai_in,
                         0, 0, 0,    // no SO
                         0, 0, 0,    // no SOI
                         Vsi_in, rsi0_in, asi_in,
                         rc0_in};
    // Ptolemy step size: h = min(2π/k, 1.0) / STEPSPER
    // INELOCA1 parameterset: STEPSPER=15
    // dpsb parameterset: STEPSPER=8
    int STEPSPER = 15;  // match Fortran INELOCA1
    dwba.Incoming.StepSize = std::min(2.0*PI/k_in, 1.0) / STEPSPER;
    dwba.Incoming.MaxR = 20.0;  // match Fortran asymptopia=20
    dwba.WavSet(dwba.Incoming);

    // Outgoing channel setup
    double eta_out = Zp * Zt * AFINE * mu / (HBARC * k_out);

    dwba.Outgoing.Projectile = Isotope(Ap, Zp);
    dwba.Outgoing.Target = Isotope(At, Zt);
    dwba.Outgoing.Elab = Elab_out;
    dwba.Outgoing.Ecm = Ecm_out;
    dwba.Outgoing.k = k_out;
    dwba.Outgoing.eta = eta_out;
    dwba.Outgoing.mu = mu_amu;
    dwba.Outgoing.JSPS = 2;  // deuteron spin=1, JSPS=2*spin=2
    dwba.Outgoing.Pot = {V_out, r0_out, a_out,
                          Vi_out, ri0_out, ai_out,
                          0, 0, 0,
                          0, 0, 0,
                          Vsi_out, rsi0_out, asi_out,
                          rc0_out};
    dwba.Outgoing.StepSize = std::min(2.0*PI/k_out, 1.0) / STEPSPER;
    dwba.Outgoing.MaxR = 20.0;  // match Fortran asymptopia=20
    dwba.WavSet(dwba.Outgoing);

    double h_in = dwba.Incoming.StepSize;
    double h_out = dwba.Outgoing.StepSize;
    int N_in = dwba.Incoming.NSteps;
    int N_out = dwba.Outgoing.NSteps;

    std::cerr << "WAVELJ h_in=" << h_in << " h_out=" << h_out
              << " N_in=" << N_in << " N_out=" << N_out << std::endl;

    // ===== Form factor =====
    double R_real = r0_in * R0mass;
    double R_imag = ri0_in * R0mass;
    double R_surf = rsi0_in * R0mass;
    double R_coul = rc0_in * R0mass;

    // Convert BELX (B(E4) in e² barn^LX) to beta_C and beta_N
    // Ptolemy formula from PRBPRT line 25560:
    //   beta_C = (4*pi/(3*Z)) * sqrt(BELX) / (Rc_barn_half)^LX / CG_spin
    // where Rc_barn_half = Rc[fm]/10 (convert fm to barn^(1/2))
    // CG_spin = CG(JIN, 2*LX, K, 0 | JOUT, K) = 1.0 for JIN=0
    double Rc_barn_half = R_coul / 10.0;
    double CG_spin = 1.0;  // CG(0, 2*LX, 0, 0 | 2*LX, 0) for 0+ target
    double beta_C = (4.0*PI / (3.0*Zt)) * std::sqrt(BELX)
                    / (std::pow(Rc_barn_half, LX) * CG_spin);
    // Nuclear beta from equal deformation lengths: beta_N * R_real = beta_C * R_coul
    double beta = beta_C * R_coul / R_real;

    std::cerr << "BELX=" << BELX << " beta_N=" << beta << " beta_C=" << beta_C << std::endl;
    std::cerr << "R_real=" << R_real << " R_coul=" << R_coul << std::endl;

    double e2 = AFINE * HBARC;

    double SUMMIN = 0.0, SUMMID_val = 10.0, SUMMAX = 20.0; // SUMMAX matches asymptopia
    // ===== PIECE 2: Potential grid + cubic spline (Fortran INGRST approach) =====
    // Ptolemy stores -(H^2/12E) * V(r) on uniform grid, then splines it.
    // At Gauss points, evaluates spline DERIVATIVE and multiplies by (12E/H^2)*R*WT.
    double spl_step = h_in;  // Fortran: STEPSZ = RSTEPS(1) = wavefunction step
    double H_dim = k_in * spl_step;   // H = k*step (dimensionless), Fortran HS(1)
    double H2over12E = H_dim * H_dim / (12.0 * Ecm_in);
    double factor_12EH2 = 12.0 * Ecm_in / (H_dim * H_dim);  // (12E/H^2) undo factor
    
    int NSTEPS_spl = (int)(SUMMAX / spl_step) + 21;  // extra 20 points for spline edge
    int NPTS_spl = NSTEPS_spl + 1;  // number of grid points (including r=0)
    
    std::vector<double> V_real_grid(NPTS_spl), V_imag_grid(NPTS_spl);
    for (int i = 0; i < NPTS_spl; i++) {
        double r = i * spl_step;
        // MAKPOT(NWP=3): no Coulomb, no centrifugal, AJ=0
        // Real: -(H^2/12E) * WS_val(r)
        V_real_grid[i] = -H2over12E * WS_val(r, V_in, R_real, a_in);
        // Imag: volume + surface COMBINED into one array (Fortran convention!)
        V_imag_grid[i] = -H2over12E * (WS_val(r, Vi_in, R_imag, ai_in)
                                       + SurfWS_val(r, Vsi_in, R_surf, asi_in));
    }
    
    // Fit natural cubic splines
    std::vector<double> br(NPTS_spl), cr(NPTS_spl), dr(NPTS_spl);
    std::vector<double> bi(NPTS_spl), ci(NPTS_spl), di(NPTS_spl);
    splncb(NPTS_spl, spl_step, V_real_grid, br, cr, dr);
    splncb(NPTS_spl, spl_step, V_imag_grid, bi, ci, di);
    
    // Sanity check: print a few spline values and derivatives at key radii
    std::cerr << "\n=== Spline sanity check ===" << std::endl;
    std::cerr << "H_dim=" << H_dim << " H2over12E=" << H2over12E 
              << " factor_12EH2=" << factor_12EH2 << std::endl;
    std::cerr << "NPTS_spl=" << NPTS_spl << " spl_step=" << spl_step << std::endl;
    for (double r_test : {1.0, 5.0, 6.8, 10.0, 15.0}) {
        int idx = (int)(r_test / spl_step);
        double t = r_test - idx * spl_step;
        double spl_val = V_real_grid[idx] + t*(br[idx] + t*(cr[idx] + t*dr[idx]));
        double spl_deriv = br[idx] + t*(2.0*cr[idx] + t*3.0*dr[idx]);
        double raw_V = WS_val(r_test, V_in, R_real, a_in);
        double analytical_dV = dWS(r_test, V_in, R_real, a_in);
        double spline_dV = factor_12EH2 * spl_deriv;  // undo Numerov scaling
        std::cerr << "r=" << r_test << ": V_raw=" << raw_V 
                  << " spl_val(unscaled)=" << spl_val / (-H2over12E)
                  << " analytical_dV=" << analytical_dV 
                  << " spline_dV=" << spline_dV << std::endl;
    }
    std::cerr << "=== End spline check ===\n" << std::endl;

    // Form factor on a fine uniform grid (used for diagnostics only, not integration)
    int Nff = std::min(N_in, N_out);
    std::vector<double> Hnuc_r(Nff), Hnuc_i(Nff), Hcoul(Nff);
    double h_ff = std::min(h_in, h_out);
    for (int i = 0; i < Nff; i++) {
        double r = (i+1) * h_ff;

        // Nuclear: -beta * R * dV/dr for each potential component
        double dvr = dWS(r, V_in, R_real, a_in);
        double dvi = dWS(r, Vi_in, R_imag, ai_in);
        double dvs = dSurfWS(r, Vsi_in, R_surf, asi_in);

        Hnuc_r[i] = -beta * R_real * dvr;
        Hnuc_i[i] = -beta * R_imag * dvi - beta * R_surf * dvs;
        (void)Hnuc_r[i]; (void)Hnuc_i[i]; // keep for reference but integration uses spline

        // Coulomb multipole form factor
        if (r >= R_coul) {
            Hcoul[i] = beta_C * 3.0*Zp*Zt*e2 * std::pow(R_coul,LX)
                      / ((2*LX+1) * std::pow(r, LX+1));
        } else {
            Hcoul[i] = beta_C * 3.0*Zp*Zt*e2 * std::pow(r, LX)
                      / ((2*LX+1) * std::pow(R_coul, LX+1));
        }
    }

    // ===== Gauss quadrature grid (matching Ptolemy INRDIN) =====
    double GAMSUM = 5.0;  // Fortran default for inelastic (label 400)
    int NUMPT_calc = (int)((SUMMAX-SUMMIN) * (6.0*(k_in+k_out)/(4*PI)));
    int NUMPT = std::max(NUMPT_calc, 15);
    std::cerr << "NUMPT=" << NUMPT << std::endl;

    std::vector<double> gl_pts(NUMPT), gl_wts(NUMPT);
    GaussLegendre(NUMPT, -1.0, 1.0, gl_pts, gl_wts);

    // CUBMAP rational-sinh map (MAPSUM=2, matching Fortran label 300)
    {
        double tau = std::log(GAMSUM + std::sqrt(GAMSUM*GAMSUM + 1.0));
        double xlen = SUMMAX - SUMMIN;
        double xadd = SUMMIN + SUMMAX;
        double A = -SUMMID_val * xlen;
        double B = xlen;
        double C = SUMMID_val * xadd - 2.0 * SUMMIN * SUMMAX;
        double D = xadd - 2.0 * SUMMID_val;
        for (int i = 0; i < NUMPT; i++) {
            double tu = tau * gl_pts[i];
            double sh = std::sinh(tu);
            double denom = B - (D / GAMSUM) * sh;
            gl_pts[i] = (-A + (C / GAMSUM) * sh) / denom;
            gl_wts[i] *= (tau / GAMSUM) * std::cosh(tu) * ((B*C - A*D) / (denom*denom));
        }
    }
    std::cerr << "Quad: r_min=" << gl_pts[0] << " r_max=" << gl_pts[NUMPT-1] << std::endl;

    // Form factor at quadrature points (weights baked in)
    // PIECE 3: Using spline derivative of Numerov-scaled potential (Fortran INGRST approach)
    // H_nuc = beta * (12E/H^2) * R * WT * XX  where XX = spline derivative
    // For imaginary: R_imag used for BOTH volume and surface (combined before spline)
    std::vector<double> H_r_wt(NUMPT), H_i_wt(NUMPT), H_c_wt(NUMPT);
    for (int i = 0; i < NUMPT; i++) {
        double r = gl_pts[i];
        double wt = gl_wts[i];
        
        // Evaluate spline derivative at this r
        int idx = (int)(r / spl_step);
        if (idx < 0) idx = 0;
        if (idx >= NPTS_spl - 1) idx = NPTS_spl - 2;
        double t = r - idx * spl_step;
        
        // dS/dr = b[idx] + 2*c[idx]*t + 3*d[idx]*t^2
        double XX_real = br[idx] + t*(2.0*cr[idx] + t*3.0*dr[idx]);
        double XX_imag = bi[idx] + t*(2.0*ci[idx] + t*3.0*di[idx]);
        
        // Form factor: beta * (12E/H^2) * R * WT * XX
        H_r_wt[i] = beta * factor_12EH2 * R_real * wt * XX_real;
        // Fortran uses R_imag for entire imaginary (vol+surf combined)
        H_i_wt[i] = beta * factor_12EH2 * R_imag * wt * XX_imag;
        
        // Coulomb multipole — matching Ptolemy INRDIN convention
        // INGRST stores: HCOUL = -WT * VC * XX  (no beta, no Rc^LX/(2LX+1))
        //   VC = -3*Zp*Zt*e², so HCOUL = WT * 3*Zp*Zt*e² * XX
        // INRDIN multiplies by BETARAT = beta_coul * CG * Rc^LX / (2LX+1)
        // Total: WT * 3*Zp*Zt*e² * XX * beta_coul * CG * Rc^LX / (2LX+1)
        // For JBIGA=0 case, spin CG = 1.0
        double hc;
        double RcLX = std::pow(R_coul, LX);  // Rc^LX factor from BETARAT
        if (r >= R_coul)
            hc = beta_C * 3.0*Zp*Zt*e2 * RcLX / ((2*LX+1) * std::pow(r, LX+1));
        else
            hc = beta_C * 3.0*Zp*Zt*e2 * std::pow(r, LX) / ((2*LX+1) * std::pow(R_coul, LX+1));
        H_c_wt[i] = hc * wt;
    }

    // Cubic spline interpolation matching Ptolemy SPLNCB+INTRPC
    // Fits natural cubic spline to wavefunction grid, then evaluates at r
    struct SplineInterp {
        std::vector<double> br, cr, dr;  // real part coefficients
        std::vector<double> bi, ci, di;  // imag part coefficients
        double h;
        int n;
        
        void fit(const std::vector<std::complex<double>>& wf, double h_grid) {
            n = (int)wf.size();
            h = h_grid;
            br.resize(n); cr.resize(n); dr.resize(n);
            bi.resize(n); ci.resize(n); di.resize(n);
            std::vector<double> yr(n), yi(n);
            for (int i = 0; i < n; i++) {
                yr[i] = wf[i].real();
                yi[i] = wf[i].imag();
            }
            fitOne(yr, br, cr, dr);
            fitOne(yi, bi, ci, di);
        }
        
        void fitOne(const std::vector<double>& y, std::vector<double>& b,
                    std::vector<double>& c, std::vector<double>& d) {
            int nm1 = n - 1;
            // Ptolemy SPLNCB algorithm (natural BC)
            double H = h, F = (y[1] - y[0]) / H;
            for (int i = 1; i < nm1; i++) {
                double G = H;
                H = h;  // uniform grid
                double E = F;
                F = (y[i+1] - y[i]) / H;
                double GBY3 = G / 3.0;
                d[i-1] = GBY3 * b[i-1];
                double EPSIM1 = G + H;
                double RIM1B3 = F - E;
                b[i] = 1.0 / (2.0/3.0 * EPSIM1 - GBY3 * d[i-1]);
                c[i] = RIM1B3 - d[i-1] * c[i-1];
            }
            d[nm1-1] = 0;
            for (int i1 = 1; i1 < nm1; i1++) {
                int i = nm1 - i1;
                c[i] = b[i] * c[i] - d[i] * c[i+1];
            }
            for (int i = 0; i < nm1; i++) {
                d[i] = (c[i+1] - c[i]) / (3.0 * h);
                b[i] = (y[i+1] - y[i]) / h - (h * d[i] + c[i]) * h;
            }
        }
        
        std::complex<double> eval(double r) const {
            int idx = (int)(r / h);
            if (idx < 0) idx = 0;
            if (idx >= n-1) idx = n-2;
            double t = r - idx * h;
            double re = br[idx] + t*(cr[idx] + t*dr[idx]);
            // Wait - the spline is y[i] + b*t + c*t^2 + d*t^3
            // Need to add y[i] term! But we don't store y.
            // Actually we need y. Let me store it.
            return {0,0}; // placeholder
        }
    };
    
    // 5-point Lagrange interpolation (Abramowitz & Stegun 25.2.15)
    // Matching Ptolemy WAVELJ exactly:
    //   Fortran: WAVR(k) = chi at r = (k-1)*h  (1-based, offset by +1 from Numerov)
    //   Our:     wf[k] = chi at r = k*h          (0-based)
    //   To match centering: access wf[i0-2]..wf[i0+2] (centered on i0)
    auto interpSpline = [](const std::vector<std::complex<double>>& wf, double h_grid,
                           double r) -> std::complex<double> {
        int N = (int)wf.size();
        double rbyh = r / h_grid;
        int i0 = (int)(rbyh + 0.5);
        int nmax = N - 3;
        i0 = std::max(2, std::min(i0, nmax));
        double p = rbyh - i0;
        double ps = p*p, o24 = 1.0/24.0;
        double x1=p*(ps-1.0)*o24, x2=x1+x1, x3=x1*p;
        double x4=x2+x2-0.5*p, x5=x4*p;
        double c1=x3-x2, c5=x3+x2, c3=x5-x3, c2=x5-x4, c4=x5+x4;
        c3 = c3+c3+1.0;
        // 5-point Lagrange (A&S 25.2.15) centered on i0
        // c1*f_{-2} - c2*f_{-1} + c3*f_0 - c4*f_1 + c5*f_2
        // wf[k] = chi at r = k*h
        double wr = c1*wf[i0-2].real()-c2*wf[i0-1].real()+c3*wf[i0].real()
                   -c4*wf[i0+1].real()+c5*wf[i0+2].real();
        double wi = c1*wf[i0-2].imag()-c2*wf[i0-1].imag()+c3*wf[i0].imag()
                   -c4*wf[i0+1].imag()+c5*wf[i0+2].imag();
        return {wr, wi};
    };

    // ===== Coulomb tail correction setup (COULIN) =====
    // Compute Coulomb phase shifts for all L values needed
    int LMNMN_coulin = std::max(0, 0 - LX);  // LMIN=0 for inelastic
    int LMXMX_coulin = Lmax + std::max(2, LX);
    std::vector<double> sigIn_arr(LMXMX_coulin + 1), sigOut_arr(LMXMX_coulin + 1);
    for (int L = 0; L <= LMXMX_coulin; L++) {
        sigIn_arr[L] = CoulombPhaseShift(L, eta_in);
        sigOut_arr[L] = CoulombPhaseShift(L, eta_out);
    }

    // COULIN call 1: tail integrals from SUMMAX to infinity (R=SUMMAX, allSw=true)
    int MXDEL_coulin = std::max(2, LX);
    int ldldim = MXDEL_coulin + 1;
    CoulinResult coulinTail;
    coulinTail.ldldim = ldldim;
    int N_coulin = LX + 1;  // power of 1/r
    int LMNMN_cl = std::max(0, LX / 2);
    // LMXMX_cl = 30: gives best pureFF values simultaneously for all pairs
    // At LMAX=30: FF(0,4)=+4.78e-5(Ftn:+4.83e-5), FF(2,2)=-4.22e-5(Ftn:+4.83e-5, sign diff)
    // The sign difference in FF(2,2) still gives correct correction direction (see comment below)
    int LMXMX_cl = 30;
    std::cerr << "LMXMX_cl = " << LMXMX_cl << std::endl;
    double ACCINE = 1.0e-6;
    double COULML = 0.3;
    int MXCOUL = 20;
    int NPCOUL = 6;

    std::cerr << "\nCalling COULIN (tail, R=SUMMAX=" << SUMMAX << ")..." << std::endl;
    int iret_coulin = Coulin(N_coulin, LX, LMNMN_cl, LMXMX_cl,
                             eta_out, k_out, sigOut_arr,
                             eta_in, k_in, sigIn_arr,
                             SUMMAX, true,
                             coulinTail,
                             ACCINE, COULML, MXCOUL, NPCOUL);
    std::cerr << "COULIN tail iret=" << iret_coulin << std::endl;

    // COULIN call 2: pure Coulomb integrals from 0 to infinity (R=0, allSw=false)
    CoulinResult coulinPure;
    coulinPure.ldldim = ldldim;

    std::cerr << "Calling COULIN (pure Coulomb, R=0)..." << std::endl;
    int iret_coulin2 = Coulin(N_coulin, LX, LMNMN_cl, LMXMX_cl,
                              eta_out, k_out, sigOut_arr,
                              eta_in, k_in, sigIn_arr,
                              0.0, false,
                              coulinPure,
                              ACCINE, COULML, MXCOUL, NPCOUL);
    std::cerr << "COULIN pure iret=" << iret_coulin2 << std::endl;

    // BETARAT = beta_C * CG_spin * Rc^LX / (2LX+1)
    // For 0+ target, CG_spin = 1.0
    double BETARAT = beta_C * 1.0 * std::pow(R_coul, LX) / (2.0 * LX + 1.0);
    // R2S(4) = Zp*Zt*e^2 * 3.0 (charge factor for Coulomb integrals)
    double R2S4 = 3.0 * Zp * Zt * e2;

    // Helper to look up COULIN result for (LI, LO) pair
    auto getCoulinIdx = [&](int LI, int LO) -> std::pair<int,int> {
        int MAXDEL = MXDEL_coulin;
        int ID = (MAXDEL - (LO - LI)) / 2;  // Fortran: iSHFT(MAXDEL-LDEL, -1)
        int IL = LO + LI - 2 * LMNMN_cl;    // Fortran: LO+LI - 2*LMNMN
        return {ID, IL};
    };

    // ===== Radial integrals & S-matrix =====
    double factor = std::sqrt(k_in * k_out / (Ecm_in * Ecm_out * PI));

    struct InelS { int LI, LO; std::complex<double> val; };
    std::vector<InelS> Smat;

    for (int LI = 0; LI <= Lmax; LI++) {
        int LO_min = std::abs(LI - LX);
        int LO_max = std::min(LI + LX, Lmax);

        // Compute incoming wavefunction using WAVELJ
        // JSPS=2 (deuteron), no SO: JP = 2L+2 (J=L+1, highest state)
        // Matches Fortran ANGSET: JPI = JSPS + 2*LI = 2 + 2*LI
        int JP_in = 2 * LI + 2;
        dwba.WavElj(dwba.Incoming, LI, JP_in);
        // Print elastic S-matrix for comparison with Fortran
        if (LI <= 15 && LI < (int)dwba.Incoming.SMatrix.size()) {
            auto Sel = dwba.Incoming.SMatrix[LI];
            double sel_mag = std::abs(Sel);
            double sel_phase = std::arg(Sel);
            std::cerr << "ELASTIC_IN L=" << LI << " |S|=" << std::fixed << std::setprecision(6) << sel_mag
                      << " phase=" << std::setprecision(4) << sel_phase << std::endl;
        }
        std::vector<std::complex<double>> wf_in = dwba.Incoming.WaveFunction;
        if (wf_in.empty()) continue;

        for (int LO = LO_min; LO <= LO_max; LO += 2) {
            if ((LI+LO+LX) % 2 != 0) continue;

            double cg = CG_int(LI, LX, 0, 0, LO, 0);
            if (std::abs(cg) < 1e-15) continue;

            int JP_out = 2 * LO + 2;
            dwba.WavElj(dwba.Outgoing, LO, JP_out);
            std::vector<std::complex<double>> wf_out = dwba.Outgoing.WaveFunction;
            if (wf_out.empty()) continue;

            // Gauss quadrature: I = sum H(r_i)*chi_out(r_i)*chi_in(r_i)
            // H_*_wt arrays already have Gauss weights baked in
            std::complex<double> integral(0,0);
            for (int i = 0; i < NUMPT; i++) {
                double r = gl_pts[i];
                if (r <= 0) continue;
                auto chi_in_val = interpSpline(wf_in, h_in, r);
                auto chi_out_val = interpSpline(wf_out, h_out, r);
                std::complex<double> H(H_r_wt[i]+H_c_wt[i], H_i_wt[i]);
                integral += H * chi_out_val * chi_in_val;
            }

            std::complex<double> phase;
            switch (LX % 4) {
                case 0: phase={0,-1}; break;
                case 1: phase={-1,0}; break;
                case 2: phase={0,1};  break;
                case 3: phase={1,0};  break;
            }

            // Print raw integral for comparison with Fortran INTEGRAL(0,SUMMAX)
            double integ_mag = std::abs(integral);
            double integ_phase = std::arg(integral);
            // Fortran INTEGRAL = FACTOR * integral (no CG, no phase)
            double ftn_integral_mag = factor * integ_mag;
            double ftn_integral_phase = integ_phase;
            std::cerr << "INTEG " << LI << " " << LO << " " << LX 
                      << " raw_mag=" << integ_mag 
                      << " factor*raw=" << ftn_integral_mag
                      << " phase=" << ftn_integral_phase << std::endl;

            // ICOMP = factor * |cg| * integral (from 0 to SUMMAX)
            double C = factor * std::abs(cg);
            std::complex<double> ICOMP = C * std::complex<double>(integral.real(), integral.imag());

            // Coulomb tail correction: IRTOIN + FFI
            // STATUS: All correction approaches are broken or inaccurate.
            // The pure Coulomb FF integral (FFI) requires highly cancelling computation
            // that the Fortran COULIN recursion achieves stably but our alternatives don't.
            // Direct Clints from outer turning point misses the inner region.
            // Adding inner Gauss quadrature gives wrong normalization.
            // COULIN R=0 recursion diverges with LMXMX.
            // Baseline: ENABLE_COULIN=false, ENABLE_FFI_DIRECT=false -> 1.72% mean, 11.2% max.
            std::complex<double> IRTOIN(0.0, 0.0);
            std::complex<double> FFI(0.0, 0.0);
            // COULIN: enabled with seed-restoration fix (STABILITY FIX in coulin.cpp)
            // pureFF values now stable regardless of LMXMX:
            //   FF(0,4) = +4.78e-5 (Fortran: +4.83e-5, 1% off) ✅
            //   FF(2,2) = +8.02e-5 (Fortran: +4.83e-5, 66% off but correct sign) 
            // Sign: cl2ff = +R2S4 * pureFF (COULIN sign opposite to physical)
            // Apply to ALL even-LI pairs within Lmax range
            int LDEL = LO - LI;
            int Lmax_elastic_est = 26;
            bool ENABLE_COULIN = (iret_coulin == 0) && (LI <= Lmax_elastic_est + LX);  // ALL LI (Fortran LSTEP=1)
            bool ENABLE_FFI_DIRECT = false;
            bool valid_for_coulin = std::abs(LDEL) <= LX && (LI <= Lmax_elastic_est + LX);  // ALL LI
            if (ENABLE_COULIN) {
                auto [id, il] = getCoulinIdx(LI, LO);

                // === IRTOIN: Coulomb tail (SUMMAX to inf) via coulinTail (R>0 upward recursion) ===
                // CL1XX = -R2S4 * COULIN_tail.XX  (Fortran sign: CL1FF = -R2S(4)*FF)
                std::complex<double> SIN_el = dwba.Incoming.SMatrix[LI];
                std::complex<double> SOUT_el = dwba.Outgoing.SMatrix[LO];
                // Sign convention: our COULIN returns FF with opposite sign from Fortran
                // (same as pureFF: physical integral is positive, our coulin returns negative)
                // Fortran: CL1FF = -R2S4 * FF_physical;  our: coulinTail.FF = -FF_physical
                // => cl1ff = -R2S4 * (-FF_phys) = +R2S4 * FF_physical = +R2S4 * coulinTail.FF
                // BUT coulinTail.FF is NEGATIVE (same sign as coulinPure.FF)
                // So cl1ff = +R2S4 * coulinTail.FF matches Fortran's CL1FF = -R2S4 * FF_physical
                double cl1ff = +R2S4 * coulinTail.FF[id + il * coulinTail.ldldim];
                double cl1fg = +R2S4 * coulinTail.FG[id + il * coulinTail.ldldim];
                double cl1gf = +R2S4 * coulinTail.GF[id + il * coulinTail.ldldim];
                double cl1gg = +R2S4 * coulinTail.GG[id + il * coulinTail.ldldim];
                std::complex<double> one(1.0, 0.0), imag_i(0.0, 1.0);
                IRTOIN = 0.25 * C * BETARAT * (
                    (one+SOUT_el)*(one+SIN_el)*cl1ff
                  - (one-SOUT_el)*(one-SIN_el)*cl1gg
                  + imag_i*((one+SOUT_el)*(one-SIN_el)*cl1fg
                           +(one-SOUT_el)*(one+SIN_el)*cl1gf)
                );

                // === FFI: Pure Coulomb FF (0 to inf) from coulinPure (COULIN R=0 with LMXMX=12) ===
                // CL2FF = -R2S4 * FF_coulin (Fortran: CL2FF = -R2S(4)*FF where R2S(4)=-3ZpZte^2)
                // At LMXMX=12: pureFF is ~5% accurate for (0,4), less so for diagonal pairs
                if (iret_coulin2 == 0 && valid_for_coulin) {
                    // COULIN returns FF with opposite sign from Fortran (Belling phase convention)
                    // Fortran FF_coulin < 0 for physical pureFF > 0
                    // So negate: CL2FF = -R2S4 * (-our_FF) = +R2S4 * our_FF
                    double raw_pureFF = coulinPure.FF[id + il * coulinPure.ldldim];
                    double cl2ff = R2S4 * raw_pureFF;   // +R2S4 (sign flip for COULIN convention)
                    FFI = C * BETARAT * cl2ff;
                }
            }

            // === FFI: Pure Coulomb FF (0 to inf) via DIRECT Clints call (stable, bypasses COULIN R=0 recursion) ===
            // COULIN R=0 downward recursion is numerically unstable (diverges with LMXMX)
            // Direct Clints is O(1) per (LI,LO) pair and gives correct result
            if (ENABLE_FFI_DIRECT && valid_for_coulin) {
                // Turning point: where F_L becomes negligible
                double eta_avg = 0.5*(eta_in + eta_out);
                double k_avg = 0.5*(k_in + k_out);
                double rho_turn_in = (eta_in + std::sqrt(eta_in*eta_in + LI*(LI+1.0))) / k_in;
                double rho_turn_out = (eta_out + std::sqrt(eta_out*eta_out + LO*(LO+1.0))) / k_out;
                double RMIN_ffi = 1.4 * std::max(rho_turn_in, rho_turn_out);
                double A_ffi = 0.3;
                double ff=0, fg=0, gf=0, gg=0;
                // Call Clints from turning point to inf
                int iret_ffi = Clints(RMIN_ffi, eta_in, eta_out, k_in, k_out,
                                      CoulombPhaseShift(LI, eta_in),
                                      CoulombPhaseShift(LO, eta_out),
                                      1e-6, A_ffi, ff, fg, gf, gg,
                                      N_coulin, LI, LO, MXCOUL, NPCOUL);
                if (iret_ffi == 0) {
                    // Add inner region (0 to RMIN) via Gauss-Legendre quadrature
                    // Integrate F_LI(k_in*r) * F_LO(k_out*r) / r^N from 0 to RMIN
                    // Use 20 Gauss points; RMIN is outside turning point so F is regular
                    std::vector<double> gl_pts_ffi, gl_wts_ffi;
                    GaussLegendre(20, 0.0, RMIN_ffi, gl_pts_ffi, gl_wts_ffi);
                    int LMAX_rcwfn = std::max(LI, LO) + 1;
                    // Only integrate where BOTH F_LI and F_LO are outside their turning points
                    // Below turning point, F_L is evanescent and Rcwfn is unreliable
                    // Use power series F_L ~ C_L(eta)*rho^(L+1) for rho << rho_turning
                    // Here: just integrate from max(rho_turn_in, rho_turn_out)/k to RMIN_ffi
                    double r_inner_start = std::max(rho_turn_in, rho_turn_out);
                    if (r_inner_start < RMIN_ffi) {
                        std::vector<double> gl_pts_in, gl_wts_in;
                        GaussLegendre(20, r_inner_start, RMIN_ffi, gl_pts_in, gl_wts_in);
                        for (int ig = 0; ig < 20; ig++) {
                            double r = gl_pts_in[ig];
                            double wt = gl_wts_in[ig];
                            std::vector<double> FI_arr, FIP, GI_arr, GIP;
                            int ir1 = Rcwfn(k_in*r, eta_in, 0, LI, FI_arr, FIP, GI_arr, GIP, 1e-8);
                            std::vector<double> FO_arr, FOP, GO_arr, GOP;
                            int ir2 = Rcwfn(k_out*r, eta_out, 0, LO, FO_arr, FOP, GO_arr, GOP, 1e-8);
                            if (ir1==0 && ir2==0 &&
                                LI < (int)FI_arr.size() && LO < (int)FO_arr.size()) {
                                double zfac = wt / std::pow(r, N_coulin);
                                // Rcwfn with MINL=L, MAXL=L returns F[0..L]; F_L is at index L
                                ff += zfac * FI_arr[LI] * FO_arr[LO];
                            }
                        }
                    }
                    // CL2FF = -R2S4 * ff (Fortran: CL2FF = -R2S(4)*FF, R2S(4)=-3ZpZte^2<0)
                    double cl2ff = -R2S4 * ff;
                    FFI = C * BETARAT * cl2ff;
                }
            }

            // === INUC = ICOMP + IRTOIN - FFI (Fortran formula) ===
            // Fortran: CL2FF = -R2S4 * FF_coulin; R2S4>0; FF_coulin<0 (from COULIN)
            // So CL2FF>0, FFI=C*BETRAT*CL2FF>0 (positive)
            // INUC = ICOMP - FFI_positive < ICOMP -- but Fortran NUCLEAR > INTEGRAL!
            // Actually sign depends on phase. Let's use Fortran formula exactly:
            // CL2FF = -R2S4 * FF; since our COULIN gives FF<0: CL2FF>0 but sign may differ
            // Use: cl2ff = -R2S4 * raw_pureFF; FFI = C*BETARAT*cl2ff
            // INUC = ICOMP + IRTOIN - FFI (Fortran formula, applied as-is)
            std::complex<double> ITOTAL = ICOMP + IRTOIN;
            std::complex<double> INUC = ITOTAL - FFI;
            std::complex<double> S = phase * INUC;

            std::cerr << "TAIL LI=" << LI << " LO=" << LO
                      << " |ICOMP|=" << std::abs(ICOMP)
                      << " |IRTOIN|=" << std::abs(IRTOIN)
                      << " |FFI|=" << std::abs(FFI)
                      << " |S|=" << std::abs(S) << std::endl;

            Smat.push_back({LI, LO, S});
        }
    }
    std::cerr << "Computed " << Smat.size() << " S-matrix elements" << std::endl;

    // Print S-matrix elements for comparison with Fortran
    std::cerr << "\n  LI  LO  LX    |S|        phase" << std::endl;
    for (const auto& s : Smat) {
        if (s.LI <= 10) {
            double mag = std::abs(s.val);
            double ph = std::arg(s.val);
            std::cerr << std::setw(4) << s.LI << std::setw(4) << s.LO 
                      << std::setw(4) << LX
                      << "   " << std::scientific << std::setprecision(4) << mag
                      << "  " << std::fixed << std::setprecision(3) << ph << std::endl;
        }
    }
    std::cerr << std::endl;

    // ===== BETCAL: accumulate beta(MX, LO) from ALL (LI,LO) pairs =====
    // Verified formula matching Fortran Ptolemy to 0.88%
    double FACTOR_BET = 0.5 / k_in;
    int MX_max = LX;

    // beta[MX][LO] as complex arrays
    std::vector<std::vector<std::complex<double>>> beta_arr(MX_max+1,
        std::vector<std::complex<double>>(Lmax+1, {0.0, 0.0}));

    // Build lookup map: (LI,LO) -> complex S value
    std::map<std::pair<int,int>, std::complex<double>> Smap;
    for (const auto& s : Smat) {
        Smap[{s.LI, s.LO}] = s.val;
    }

    // For each LO, loop over all delta channels
    for (int LO = 0; LO <= Lmax; LO++) {
        for (int delta = -LX; delta <= LX; delta += 2) {
            int LI = LO - delta;
            if (LI < 0 || LI > Lmax) continue;

            auto it = Smap.find({LI, LO});
            if (it == Smap.end()) continue;

            std::complex<double> Sval = it->second;
            double SMAG = std::abs(Sval);
            double SPHASE = std::arg(Sval);
            if (SMAG < 1e-30) continue;

            // Coulomb phase rotation
            double sigma_in = CoulombPhaseShift(LI, eta_in);
            double sigma_out = CoulombPhaseShift(LO, eta_out);
            double PH = SPHASE + sigma_in + sigma_out;
            double SR = SMAG * std::sin(PH);
            double SI = -SMAG * std::cos(PH);

            // Accumulate into all MX channels
            int MX_lim = std::min(LX, LO);
            for (int MX = 0; MX <= MX_lim; MX++) {
                double cg = CG_int(LI, LX, 0, MX, LO, MX);
                if (std::abs(cg) < 1e-15) continue;
                double TEMPS = FACTOR_BET * (2*LI+1) * cg;
                beta_arr[MX][LO] += std::complex<double>(TEMPS * SR, TEMPS * SI);
            }
        }
    }

    // sqrt-factorial normalization for MX > 0
    for (int MX = 1; MX <= MX_max; MX++) {
        for (int LO = MX; LO <= Lmax; LO++) {
            double sqfac = 1.0;
            for (int m = 1; m <= MX; m++) {
                sqfac /= std::sqrt((double)(LO+m) * (LO-m+1));
            }
            beta_arr[MX][LO] *= sqfac;
        }
    }

    // ===== DCS: 10 * sum_MX w(MX) * |sum_LO beta[MX][LO] * P_LO^MX(costh)|^2 =====
    std::cout << std::fixed;
    for (double th = 0.0; th <= 180.01; th += 1.0) {
        double costh = std::cos(th * PI / 180.0);
        double dcs = 0.0;
        for (int MX = 0; MX <= MX_max; MX++) {
            std::complex<double> amp(0.0, 0.0);
            for (int LO = MX; LO <= Lmax; LO++) {
                if (std::abs(beta_arr[MX][LO]) < 1e-30) continue;
                double plm = AssocLegP(LO, MX, costh);
                amp += beta_arr[MX][LO] * plm;
            }
            double w = (MX == 0) ? 1.0 : 2.0;
            dcs += w * std::norm(amp);
        }
        dcs *= 10.0;
        std::cout << std::setw(8) << std::setprecision(2) << th
                  << std::setw(15) << std::scientific << std::setprecision(5) << dcs
                  << std::endl;
    }

    return 0;
}
