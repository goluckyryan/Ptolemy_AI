// wftest_main.cpp — Test WavElj normalization directly
// Builds a Channel manually and calls WavElj for d+16O, L=0 JP=2

#include "dwba.h"
#include "potential_eval.h"
#include <cstdio>
#include <cmath>

// Forward declare a minimal DWBA subclass to expose private WavElj/WavSet
class DWBATest : public DWBA {
public:
    void testWavElj(Channel& ch, int L, int JP) { WavElj(ch, L, JP); }
    void testWavSet(Channel& ch) { WavSet(ch); }
};

int main() {
    const double AMU_MEV = 931.494061;
    const double HBARC   = 197.32697;

    // Incoming: d+16O, Elab=20 MeV
    double mA=16.0, ma=2.0;
    double Ecm = mA/(mA+ma)*20.0;
    double mu  = mA*ma/(mA+ma);
    double k   = std::sqrt(2*mu*AMU_MEV*Ecm)/HBARC;
    double eta = 1.0*8.0*mu*AMU_MEV/(137.036*HBARC*k);

    printf("Incoming d+16O: k=%.6f, eta=%.6f\n", k, eta);
    printf("Ptolemy: k=1.2328, step=0.125 fm\n\n");

    Channel ch;
    ch.Target.A     = 16; ch.Target.Z     = 8;  ch.Target.Mass  = 16.0;
    ch.Projectile.A = 2;  ch.Projectile.Z = 1;  ch.Projectile.Mass = 2.0;
    ch.mu   = mu;
    ch.k    = k;
    ch.eta  = eta;
    ch.JSPS = 2;  // 2*S for deuteron

    ch.Pot.V    = 88.9546; ch.Pot.R0   = 1.1489; ch.Pot.A   = 0.7508;
    ch.Pot.VI   = 2.3480;  ch.Pot.RI0  = 1.3446; ch.Pot.AI  = 0.6030;
    ch.Pot.VSI  = 10.2180; ch.Pot.RSI0 = 1.3943; ch.Pot.ASI = 0.6872;
    ch.Pot.VSO  = 7.1140;  ch.Pot.RSO0 = 0.9720; ch.Pot.ASO = 1.0110;
    ch.Pot.VSOI = 0.0;
    ch.Pot.RC0  = 1.3030;

    DWBATest dw;
    dw.testWavSet(ch);

    printf("h=%.4f fm, N=%d\n\n", ch.StepSize, ch.NSteps);

    // Test L=0, JP=2 (incoming, J=1 for deuteron)
    dw.testWavElj(ch, 0, 2);

    auto S = ch.SMatrix[0];
    printf("WavElj S-matrix L=0 JP=2/2: (%+.5f, %+.5f)\n", S.real(), S.imag());
    printf("Ptolemy:                      (+0.14429, +0.11977)\n\n");

    printf("Wavefunction chi (WavElj normalization):\n");
    printf("Ptolemy STARTS: r=0.125: (0.09815, 0.03753), r=0.250: (0.18354, 0.06980)\n\n");
    printf("  r(fm)   Re(chi)         Im(chi)\n");
    for (int i = 0; i <= (int)(5.0/ch.StepSize); i++) {
        if (i >= (int)ch.WaveFunction.size()) break;
        double r = i * ch.StepSize;
        printf("  %.3f  %+.7e  %+.7e\n", r, ch.WaveFunction[i].real(), ch.WaveFunction[i].imag());
    }

    return 0;
}
