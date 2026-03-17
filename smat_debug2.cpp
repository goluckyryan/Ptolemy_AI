// smat_debug2.cpp — trace Numerov integration for d+16O L=0
// Count BIGNUM rescalings and compare f(r) between ES and WavElj

#include "elastic.h"
#include "dwba.h"
#include "potential_eval.h"
#include "rcwfn.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

class DWBATest : public DWBA {
public:
    void testWavSet(Channel& ch) { WavSet(ch); }
    void testWavElj(Channel& ch, int L, int JP) { WavElj(ch, L, JP); }
};

int main() {
    const double AMU_MEV = 931.494061, HBARC = 197.32697;
    double mA=16.0, ma=2.0;
    double Ecm = mA/(mA+ma)*20.0;
    double mu  = mA*ma/(mA+ma);
    double k   = std::sqrt(2*mu*AMU_MEV*Ecm)/HBARC;
    double eta = 1.0*8.0*mu*AMU_MEV/(137.036*HBARC*k);

    printf("d+16O: k=%.6f eta=%.6f\n", k, eta);

    // Build potential arrays the same way ElasticSolver does
    // to compare f(r) at key points
    double h = 0.1;
    int N = 301;
    int L = 0;

    // ---- Manually build f(r) from WavElj potentials (via EvaluatePotential) ----
    Channel ch;
    ch.Target.A=16; ch.Target.Z=8; ch.Target.Mass=16.0;
    ch.Projectile.A=2; ch.Projectile.Z=1; ch.Projectile.Mass=2.0;
    ch.mu=mu; ch.k=k; ch.eta=eta; ch.JSPS=2;
    ch.Pot.V=88.9546; ch.Pot.R0=1.1489; ch.Pot.A=0.7508;
    ch.Pot.VI=2.3480;  ch.Pot.RI0=1.3446; ch.Pot.AI=0.6030;
    ch.Pot.VSI=10.2180; ch.Pot.RSI0=1.3943; ch.Pot.ASI=0.6872;
    ch.Pot.VSO=7.1140;  ch.Pot.RSO0=0.9720; ch.Pot.ASO=1.0110;
    ch.Pot.VSOI=0.0; ch.Pot.RC0=1.3030;

    DWBATest dw;
    dw.testWavSet(ch);  // fills V_real, V_imag, etc.

    double f_conv = 2.0*mu*AMU_MEV/(HBARC*HBARC);
    double k2 = k*k;
    int JSPS = 2;  // 2*S for deuteron
    // JP=2, L=0: SDOTL = 0.25*(2*(2+2) - 2*(2+2)) - 0 = 0.25*(8-8)=0? No:
    // SDOTL_raw = 0.25*(JP*(JP+2) - JSPS*(JSPS+2)) - L*(L+1)
    // JP=2, JSPS=2, L=0: 0.25*(2*4 - 2*4) - 0 = 0
    // spin_dot_L = 0/2 = 0 for L=0!
    double JP = 2;
    double SDOTL_raw = 0.25*(JP*(JP+2) - JSPS*(JSPS+2)) - (double)L*(L+1);
    double spin_dot_L = (JSPS > 0) ? SDOTL_raw / JSPS : 0.0;
    printf("L=0, JP=2, JSPS=2: SDOTL_raw=%.4f, spin_dot_L=%.4f\n\n", SDOTL_raw, spin_dot_L);

    // Build f(r) from WavElj grid
    printf("  r       Vr        Wi        Vc        Vso     f_re_WJ    f_im_WJ\n");
    for (double r : {0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0}) {
        int idx = (int)(r/h);
        if (idx >= N) idx = N-1;
        double Vr  = ch.V_real[idx];
        double Wi  = ch.V_imag[idx];
        double Vc  = ch.V_coulomb[idx];
        double Vso = ch.V_so_real[idx];
        double LL1 = (double)L*(L+1);
        double f_re = k2 - LL1/(r*r) - f_conv*Vc + f_conv*Vr + spin_dot_L*f_conv*Vso;
        double f_im = f_conv*Wi + spin_dot_L*f_conv*ch.V_so_imag[idx];
        printf("  %5.1f   %7.3f   %6.3f   %7.3f   %7.4f   %+.4f   %+.4f\n",
               r, Vr, Wi, Vc, Vso, f_re, f_im);
    }

    printf("\n--- Numerov trace for WavElj (L=0 JP=2) ---\n");
    // Reproduce WavElj's Numerov to count BIGNUM hits and see u at matching
    std::vector<std::complex<double>> f(N+2, 0.0);
    for (int i = 0; i <= N; ++i) {
        double r = (i==0) ? 1e-10 : ((i<N) ? ch.RGrid[i] : i*h);
        int idx = std::min(i, N-1);
        double Vr=ch.V_real[idx], Wi=ch.V_imag[idx], Vc=ch.V_coulomb[idx];
        double Vso=ch.V_so_real[idx], Vsoi=ch.V_so_imag[idx];
        double LL1=(double)L*(L+1);
        double f_re = k2 - LL1/(r*r) - f_conv*Vc + f_conv*Vr + spin_dot_L*f_conv*Vso;
        double f_im = f_conv*Wi + spin_dot_L*f_conv*Vsoi;
        f[i] = {f_re, f_im};
    }

    double h2_12 = h*h/12.0;
    std::vector<std::complex<double>> u(N+2, 0.0);
    u[0]=0.0; u[1]=std::pow(h, L+1);
    const double BIGNUM=1e30, STEPI=1e-10;
    int ISTRT=0, n_rescale=0;

    for (int i=1; i<N; ++i) {
        auto t1 = 2.0*(1.0-5.0*h2_12*f[i])*u[i];
        auto t2 = (1.0+h2_12*f[i-1])*u[i-1];
        auto dn = 1.0+h2_12*f[i+1];
        if (std::abs(dn)<1e-300) dn=1e-300;
        u[i+1]=(t1-t2)/dn;
        double mag=std::abs(u[i+1]);
        if (mag>BIGNUM) {
            n_rescale++;
            double inv_mag=1.0/mag;
            double threshold=STEPI*mag;
            int new_istrt=i+1;
            for (int j=ISTRT; j<=i+1; ++j) {
                if (std::abs(u[j])>=threshold) { new_istrt=j; break; }
            }
            for (int j=ISTRT; j<new_istrt; ++j) u[j]=0.0;
            ISTRT=new_istrt;
            for (int j=ISTRT; j<=i+1; ++j) u[j]*=inv_mag;
        }
    }

    int n_match = N-3;
    printf("BIGNUM rescalings: %d\n", n_rescale);
    printf("ISTRT after integration: %d (r=%.3f fm)\n", ISTRT, ISTRT*h);
    printf("u[n_match=%d] = (%+.6e, %+.6e)\n", n_match, u[n_match].real(), u[n_match].imag());

    // Get F and G at matching
    std::vector<double> FC(L+2), FCP(L+2), GC(L+2), GCP(L+2);
    double rho_m = k*(n_match*h);
    Rcwfn(rho_m, eta, 0, L, FC, FCP, GC, GCP);
    printf("F[L=%d]=%.6e  G[L=%d]=%.6e  at rho=%.4f\n", L, FC[L], L, GC[L], rho_m);

    // Compute S-matrix (WavElj single-point Wronskian)
    double FL=FC[L], GL=GC[L];
    double FLP=FCP[L]*k, GLP=GCP[L]*k;
    // 5-point derivative
    auto u_prime = (u[n_match-2]-8.0*u[n_match-1]+8.0*u[n_match+1]-u[n_match+2])/(12.0*h);
    auto u_m = u[n_match];
    double ur=u_m.real(), ui=u_m.imag(), upr=u_prime.real(), upi=u_prime.imag();
    double Ar=FLP*ur-upr*FL, Ai=FLP*ui-upi*FL;
    double Br=upr*GL-GLP*ur, Bi=upi*GL-GLP*ui;
    double num_r=Br-Ai, num_i=Bi+Ar;
    double den_r=Br+Ai, den_i=Bi-Ar;
    double den=den_r*den_r+den_i*den_i;
    double SJR=(den>1e-60)?(num_r*den_r+num_i*den_i)/den:0.0;
    double SJI=(den>1e-60)?(num_i*den_r-num_r*den_i)/den:0.0;
    printf("WavElj S (manual): (%+.6f, %+.6f)\n", SJR, SJI);

    // --- Now do ElasticSolver's two-point matching ---
    // ElasticSolver rescales WITHOUT zeroing inner region
    std::vector<std::complex<double>> u_es(N+2, 0.0);
    u_es[0]=0.0; u_es[1]=std::pow(h,L+1);
    int n_rescale_es=0;
    for (int i=1; i<N; ++i) {
        auto t1=2.0*(1.0-5.0*h2_12*f[i])*u_es[i];
        auto t2=(1.0+h2_12*f[i-1])*u_es[i-1];
        auto dn=1.0+h2_12*f[i+1];
        if (std::abs(dn)<1e-300) dn=1e-300;
        u_es[i+1]=(t1-t2)/dn;
        double mag=std::abs(u_es[i+1]);
        if (mag>BIGNUM) {
            n_rescale_es++;
            for (int j=0; j<=i+1; ++j) u_es[j]/=mag;
        }
    }

    int n1=N-4, n2=N-3;
    std::vector<double> FC1(L+2),FCP1(L+2),GC1(L+2),GCP1(L+2);
    std::vector<double> FC2(L+2),FCP2(L+2),GC2(L+2),GCP2(L+2);
    Rcwfn(k*(n1*h), eta, 0, L, FC1, FCP1, GC1, GCP1);
    Rcwfn(k*(n2*h), eta, 0, L, FC2, FCP2, GC2, GCP2);
    double f1=FC1[L], g1=GC1[L], f2=FC2[L], g2=GC2[L];
    std::complex<double> u1=u_es[n1], u2=u_es[n2];
    double det=f2*g1-f1*g2;
    using namespace std::complex_literals;
    std::complex<double> A=(f2*u1-u2*f1)/det;
    std::complex<double> B=(u2*g1-g2*u1)/det;
    std::complex<double> S_es=(B+1i*A)/(B-1i*A);
    printf("ElasticSolver S (manual): (%+.6f, %+.6f)\n", S_es.real(), S_es.imag());
    printf("ES rescalings: %d\n", n_rescale_es);

    // Print u at matching for both
    printf("\nu_es[n1=%d] = (%+.6e, %+.6e)\n", n1, u_es[n1].real(), u_es[n1].imag());
    printf("u_es[n2=%d] = (%+.6e, %+.6e)\n", n2, u_es[n2].real(), u_es[n2].imag());
    printf("u_wj[n_m=%d] = (%+.6e, %+.6e)\n", n_match, u[n_match].real(), u[n_match].imag());

    return 0;
}
