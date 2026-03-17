// smat_debug3.cpp — Compare ElasticSolver's internal f(r) vs WavElj's f(r)
// for d+16O L=0

#include "elastic.h"
#include "dwba.h"
#include "potential_eval.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

class DWBATest : public DWBA {
public:
    void testWavSet(Channel& ch) { WavSet(ch); }
};

int main() {
    const double AMU_MEV = 931.494061, HBARC = 197.32697;
    double mA=16.0, ma=2.0;
    double Ecm = mA/(mA+ma)*20.0;
    double mu  = mA*ma/(mA+ma);
    double k   = std::sqrt(2*mu*AMU_MEV*Ecm)/HBARC;
    double eta = 1.0*8.0*mu*AMU_MEV/(137.036*HBARC*k);

    double h = 0.1;
    int N = 301, L = 0;
    double f_conv = 2.0*mu*AMU_MEV/(HBARC*HBARC);
    double k2 = k*k;

    // --- WavElj potential (via EvaluatePotential) ---
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
    dw.testWavSet(ch);

    // --- ElasticSolver's BuildPotentialArrays ---
    ElasticSolver es;
    es.SetSystem(2, 1, 16, 8, 20.0);
    es.AddVolumeWS ({88.9546, 0.0   }, 1.1489, 0.7508);
    es.AddVolumeWS ({0.0,     2.3480}, 1.3446, 0.6030);
    es.AddSurfaceWS({0.0,  -10.2180 }, 1.3943, 0.6872);
    es.AddSpinOrbit({7.1140,  0.0   }, 0.9720, 1.0110);
    es.AddCoulomb(1.3030);
    es.SetGrid(h, 30.0);
    es.SetLmax(2);

    // Access ElasticSolver's potential via its internal arrays
    // We can't directly, so compute them externally using the same WS formulas

    // ElasticSolver::BuildPotentialArrays fills Vr, Wi, Vc, VsoRe, VsoIm
    // Let's call Solve() but intercept: since we can't, just compare
    // the f(r) values by running both through Numerov with debug prints

    // Actually: print WavElj's raw potentials vs what ElasticSolver should get
    printf("Comparing raw potential components at key radii:\n");
    printf("(WavElj via EvaluatePotential)\n\n");

    printf("  r    V_real(WJ)  V_imag(WJ)  V_coul(WJ)  Vso_re(WJ)\n");
    for (double r : {0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0}) {
        int idx = (int)(r/h);
        if (idx >= N) idx = N-1;
        printf("  %4.1f  %10.4f  %10.4f  %10.4f  %10.4f\n",
               r, ch.V_real[idx], ch.V_imag[idx], ch.V_coulomb[idx], ch.V_so_real[idx]);
    }

    // Now compute what ElasticSolver's WS formula gives at the same points
    // ElasticSolver AddVolumeWS({V,W}, r0, a): Vr = V*fWS, Wi = W*fWS
    // AddSurfaceWS({V,W}, r0, a): Wi += W * 4*exp/(1+exp)^2 (note: W<0 for absorptive)
    // AddSpinOrbit({Vso,Vsoi}, r0, a): VsoRe = Vso * 4*(df/dr)/r, VsoIm = Vsoi*same
    // Coulomb: standard

    // Let's compute manually for r=2 fm
    printf("\n--- Manual computation at r=2.0 fm ---\n");
    double r = 2.0;
    double At = 16.0;
    double R0_target = std::pow(At, 1.0/3.0);  // r0target convention: At^1/3

    // Volume WS (real)
    double R_vol = 1.1489 * R0_target;
    double fWS_vol = 1.0/(1.0+std::exp((r-R_vol)/0.7508));
    double Vr_ES = 88.9546 * fWS_vol;
    printf("  ES V_real: R_vol=%.4f, fWS=%.6f, Vr=%.4f\n", R_vol, fWS_vol, Vr_ES);

    // Volume WS (imag)
    double R_vol_i = 1.3446 * R0_target;
    double fWS_vol_i = 1.0/(1.0+std::exp((r-R_vol_i)/0.6030));
    double Wi_vol_ES = 2.3480 * fWS_vol_i;
    printf("  ES V_imag(vol): R_vol_i=%.4f, fWS=%.6f, Wi=%.4f\n", R_vol_i, fWS_vol_i, Wi_vol_ES);

    // Surface WS (imag)
    double R_si = 1.3943 * R0_target;
    double e_si = std::exp((r-R_si)/0.6872);
    double fWS_si = e_si/(1.0+e_si);  // exp/(1+exp)
    double dWS_si = fWS_si*(1.0-fWS_si)/0.6872;  // derivative / (1+exp)^2 * exp term
    // Actually surface = -d/dr[fWS] = exp/(1+exp)^2 / a → standard form
    double surf_si = e_si/((1.0+e_si)*(1.0+e_si)*0.6872);  // 4*a is implicit? No.
    // ElasticSolver AddSurfaceWS uses: Wi_surf = W * 4*exp/(1+exp)^2 ?
    // Check elastic.cpp BuildPotentialArrays
    printf("  ES V_imag(surf) needs to check BuildPotentialArrays implementation\n");

    // WavElj at r=2 fm:
    int idx = (int)(r/h);
    printf("\n  WavElj at r=%.1f:\n", r);
    printf("    V_real=%.6f  V_imag=%.6f  V_coul=%.6f  Vso=%.6f\n",
           ch.V_real[idx], ch.V_imag[idx], ch.V_coulomb[idx], ch.V_so_real[idx]);

    return 0;
}
