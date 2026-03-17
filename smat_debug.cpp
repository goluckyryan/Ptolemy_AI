// smat_debug.cpp — Compare ElasticSolver vs WavElj S-matrix for d+16O and 60Ni(d,d)
// Goal: find exactly why d+16O S-matrix differs between the two implementations

#include "elastic.h"
#include "dwba.h"
#include "potential_eval.h"
#include <cstdio>
#include <cmath>

class DWBATest : public DWBA {
public:
    void testWavSet(Channel& ch) { WavSet(ch); }
    void testWavElj(Channel& ch, int L, int JP) { WavElj(ch, L, JP); }
};

void test_case(const char* label,
               // ElasticSolver params
               int Ap, int Zp, int At, int Zt, double Elab,
               // Potentials: V, RI0=R0, AI=A, W_vol, RI0_vol, AI_vol, VSI, RSI0, ASI, VSO, RSO0, ASO, RC0
               double V, double R0, double A,
               double W_vol, double RI0_vol, double AI_vol,
               double VSI, double RSI0, double ASI,
               double VSO, double RSO0, double ASO,
               double VSOI, double RSOI0, double ASOI,
               double RC0,
               // Ptolemy reference S-matrix for L=0
               double pto_SJR, double pto_SJI, int twoJP_L0)
{
    printf("\n=== %s ===\n", label);

    // --- ElasticSolver ---
    ElasticSolver es;
    es.SetSystem(Ap, Zp, At, Zt, Elab);
    es.AddVolumeWS ({V,     0.0  }, R0,    A);
    es.AddVolumeWS ({0.0,   W_vol}, RI0_vol, AI_vol);
    es.AddSurfaceWS({0.0,  -VSI  }, RSI0,  ASI);
    es.AddSpinOrbit({VSO,  VSOI  }, RSO0,  ASO);
    es.AddCoulomb(RC0);
    es.SetGrid(0.1, 30.0);
    es.SetLmax(4);
    es.Solve();
    auto Ses = es.GetSMatrix(0, twoJP_L0);
    printf("  ElasticSolver  L=0 JP=%d/2: S=(%+.5f, %+.5f)\n", twoJP_L0, Ses.real(), Ses.imag());

    // --- WavElj ---
    const double AMU_MEV = 931.494061, HBARC = 197.32697;
    double mA = At, ma = Ap;
    double Ecm = mA/(mA+ma)*Elab;
    double mu  = mA*ma/(mA+ma);
    double k   = std::sqrt(2*mu*AMU_MEV*Ecm)/HBARC;
    double eta = (double)Zp*(double)Zt*mu*AMU_MEV/(137.036*HBARC*k);

    Channel ch;
    ch.Target.A=At; ch.Target.Z=Zt; ch.Target.Mass=(double)At;
    ch.Projectile.A=Ap; ch.Projectile.Z=Zp; ch.Projectile.Mass=(double)Ap;
    ch.mu=mu; ch.k=k; ch.eta=eta;
    ch.JSPS = (Ap==2) ? 2 : 1;   // deuteron→2, proton→1
    ch.Pot.V=V; ch.Pot.R0=R0; ch.Pot.A=A;
    ch.Pot.VI=W_vol; ch.Pot.RI0=RI0_vol; ch.Pot.AI=AI_vol;
    ch.Pot.VSI=VSI; ch.Pot.RSI0=RSI0; ch.Pot.ASI=ASI;
    ch.Pot.VSO=VSO; ch.Pot.RSO0=RSO0; ch.Pot.ASO=ASO;
    ch.Pot.VSOI=VSOI; ch.Pot.RSOI0=RSOI0; ch.Pot.ASOI=ASOI;
    ch.Pot.RC0=RC0;

    DWBATest dw;
    dw.testWavSet(ch);
    dw.testWavElj(ch, 0, twoJP_L0);
    auto Swj = ch.SMatrix[0];
    printf("  WavElj         L=0 JP=%d/2: S=(%+.5f, %+.5f)\n", twoJP_L0, Swj.real(), Swj.imag());
    printf("  Ptolemy        L=0 JP=%d/2: S=(%+.5f, %+.5f)\n", twoJP_L0, pto_SJR, pto_SJI);

    printf("  ES vs Pto:  dRe=%+.5f dIm=%+.5f\n", Ses.real()-pto_SJR, Ses.imag()-pto_SJI);
    printf("  WJ vs Pto:  dRe=%+.5f dIm=%+.5f\n", Swj.real()-pto_SJR, Swj.imag()-pto_SJI);

    // Key params printed by ElasticSolver
    printf("  k=%.6f  eta=%.6f  Ap=%d  At=%d\n", k, eta, Ap, At);

    // Print what R values the potentials get
    double R_target = std::pow((double)At, 1.0/3.0);
    printf("  At^1/3=%.5f  R_vol=%.5f  RC=%.5f\n",
           R_target, R0*R_target, RC0*R_target);
}

int main() {
    // Case 1: d + 16O at Elab=20 MeV (INCOMING channel)
    // Ptolemy: L=0 JP=2/2: S=(0.14429, 0.11977)
    test_case("d + 16O at Elab=20 MeV",
        2, 1, 16, 8, 20.0,
        88.9546, 1.1489, 0.7508,   // V, R0, A
        2.3480, 1.3446, 0.6030,    // VI, RI0, AI  (volume imag)
        10.2180, 1.3943, 0.6872,   // VSI, RSI0, ASI (surface)
        7.1140, 0.9720, 1.0110,    // VSO, RSO0, ASO
        0.0, 0.0, 0.0,             // VSOI
        1.3030,                    // RC0
        0.14429, 0.11977, 2);      // Ptolemy ref, twoJP=2

    // Case 2: 60Ni + d at Elab=60 MeV (validated)
    // From elastic_test.cpp validated case
    // An-Cai 2006 energy-dependent params at 60 MeV
    test_case("d + 60Ni at Elab=60 MeV (validated)",
        2, 1, 60, 28, 60.0,
        81.919, 1.150, 0.768,      // V, R0, A
        4.836, 1.330, 0.464,       // VI (separate R0/a!)
        8.994, 1.373, 0.774,       // VSI
        3.557, 0.972, 1.011,       // VSO (half-strength for Ptolemy match)
        0.0, 0.0, 0.0,
        1.303,
        // Ptolemy JP-basis L=0 JP=2/2 for 60Ni(d,d) — need to check
        // From our previous validation: we compared to JP-basis S-matrix
        // Using elastic_test reference values
        0.0, 0.0, 2);  // placeholder

    // Case 3: d + 16O — try with HALF VSO (Ptolemy convention)
    test_case("d + 16O at Elab=20 MeV (VSO=3.557, half)",
        2, 1, 16, 8, 20.0,
        88.9546, 1.1489, 0.7508,
        2.3480, 1.3446, 0.6030,
        10.2180, 1.3943, 0.6872,
        3.5570, 0.9720, 1.0110,    // VSO halved
        0.0, 0.0, 0.0,
        1.3030,
        0.14429, 0.11977, 2);

    // Case 4: d + 16O with L=0 JP=0/2 (J=0, lower branch — should be zero for L=0, S=1)
    // Check if twoJP=0 is the issue (invalid for L=0, S=1)
    printf("\n=== Check valid JP values for d+16O L=0 (S=1, JSPS=2) ===\n");
    printf("  Valid JP for L=0, S=1: J=|L-S|..L+S = 1 only (J=0 invalid!)\n");
    printf("  So twoJP must be 2 (J=1). twoJP=0 (J=0) is INVALID.\n");
    printf("  In ElasticSolver for S=1: twoS=2, twoJ_min=|2L-2|=0 for L=0 → J=0 included!\n");
    printf("  But J=0 is INVALID for L=0, S=1 (triangle: |L-S|≤J≤L+S → 1≤J≤1)\n");

    // Case 5: Verify the AddSurfaceWS sign convention
    // In ElasticSolver: AddSurfaceWS({0, -VSI}, RSI0, ASI) stores WS as -VSI (attractive)
    // WavElj: ch.Pot.VSI = 10.218 (positive), and EvaluatePotential uses sign convention?
    printf("\n=== Check potential evaluation sign for surface term ===\n");
    printf("  ElasticSolver: AddSurfaceWS({0,-10.218}, 1.3943, 0.6872)\n");
    printf("  WavElj: ch.Pot.VSI=10.218 — check EvaluatePotential sign convention\n");

    return 0;
}
