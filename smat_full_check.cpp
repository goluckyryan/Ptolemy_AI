// smat_full_check.cpp — Compare WavElj S-matrix vs Ptolemy for both channels
// Tests all L=0..8 with all valid J for both incoming (d+16O) and outgoing (p+17O)

#include "dwba.h"
#include "potential_eval.h"
#include "rcwfn.h"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

class DWBATest : public DWBA {
public:
    void setup_and_run(Channel& ch, int L, int JP) {
        WavSet(ch);
        WavElj(ch, L, JP);
    }
    Channel& incoming() { return Incoming; }
    Channel& outgoing() { return Outgoing; }
};

static void setup_incoming(Channel& ch) {
    // d + 16O, Elab=20 MeV, An-Cai 2006 OMP
    const double AMU = 931.494061, HBARC = 197.32697;
    ch.Projectile.A = 2; ch.Projectile.Z = 1;
    ch.Target.A = 16;    ch.Target.Z = 8;
    ch.JSPS = 2;
    double Ap=2.0, At=16.0;
    double mu = Ap*At/(Ap+At);
    double Ecm = At/(Ap+At)*20.0;
    ch.mu = mu; ch.Elab = 20.0;
    ch.k   = std::sqrt(2.0*mu*AMU*Ecm)/HBARC;
    ch.eta = 1.0*8.0*mu*AMU/(137.036*HBARC*ch.k);

    double scale = std::cbrt(At);  // r0target: R = r0 * At^(1/3)
    ch.Pot.V    = 88.9546; ch.Pot.R0   = 1.1489*scale; ch.Pot.A    = 0.7508;
    ch.Pot.VI   = 2.3480;  ch.Pot.RI0  = 1.3446*scale; ch.Pot.AI   = 0.6030;
    ch.Pot.VSI  = 10.2180; ch.Pot.RSI0 = 1.3943*scale; ch.Pot.ASI  = 0.6872;
    ch.Pot.VSO  = 7.1140;  ch.Pot.RSO0 = 0.9720*scale; ch.Pot.ASO  = 1.0110;
    ch.Pot.VSOI = 0.0;     ch.Pot.RC0  = 1.3030*scale;
}

static void setup_outgoing(Channel& ch) {
    // p + 17O, Koning-Delaroche OMP
    const double AMU = 931.494061, HBARC = 197.32697;
    ch.Projectile.A = 1; ch.Projectile.Z = 1;
    ch.Target.A = 17;    ch.Target.Z = 8;
    ch.JSPS = 1;
    double Ap=1.0, At=17.0;
    double mu = Ap*At/(Ap+At);
    ch.mu = mu;
    // Use k from Ptolemy output: k_out = 0.94628 fm^-1
    ch.k   = 0.94628;
    ch.eta = 1.0*8.0*mu*AMU/(137.036*HBARC*ch.k);

    double scale = std::cbrt(At);
    ch.Pot.V    = 49.5434; ch.Pot.R0   = 1.1462*scale; ch.Pot.A    = 0.6753;
    ch.Pot.VI   = 2.0611;  ch.Pot.RI0  = 1.1462*scale; ch.Pot.AI   = 0.6753;
    ch.Pot.VSI  = 7.6703;  ch.Pot.RSI0 = 1.3016*scale; ch.Pot.ASI  = 0.5275;
    ch.Pot.VSO  = 5.2956;  ch.Pot.RSO0 = 0.9338*scale; ch.Pot.ASO  = 0.5900;
    ch.Pot.VSOI =-0.1059;  ch.Pot.RSOI0= 0.9338*scale; ch.Pot.ASOI = 0.5900;
    ch.Pot.RC0  = 1.3030*scale;
}

int main() {
    // ── Ptolemy reference (from print=1001 output) ──────────────────────────
    // INCOMING d+16O (JSPS=2, d has 3 J channels per L)
    struct Ref { int L, JP; double Re, Im; };

    std::vector<Ref> ptol_in = {
        {0, 2,  0.14429,  0.11977},
        {1, 0,  0.14914, -0.14607},
        {1, 2,  0.18217, -0.11199},
        {1, 4,  0.21878, -0.02362},
        {2, 2, -0.09686, -0.12363},
        {2, 4, -0.01539, -0.15313},
        {2, 6,  0.10997, -0.12204},
        {3, 4, -0.18101,  0.02907},
        {3, 6, -0.19169, -0.07134},
        {3, 8, -0.09711, -0.17963},
        {4, 6,  0.11158, -0.00815},
        {4, 8,  0.00336,  0.03495},
        {4,10, -0.10824, -0.04910},
    };

    // OUTGOING p+17O (JSPS=1, proton has 2 J channels per L)
    std::vector<Ref> ptol_out = {
        {0, 1,  0.49629, -0.01952},
        {1, 1,  0.02830, -0.40199},
        {1, 3,  0.22806, -0.30539},
        {2, 3, -0.32679, -0.31757},
        {2, 5, -0.04084, -0.54572},
        {3, 5,  0.08283,  0.14530},
        {3, 7, -0.07999, -0.13415},
        {4, 7,  0.71811,  0.20865},
        {4, 9,  0.61227,  0.25849},
        {5, 9,  0.94558,  0.07602},
        {5,11,  0.93736,  0.09693},
        {6,11,  0.99003,  0.02091},
        {6,13,  0.98956,  0.02517},
        {7,13,  0.99814,  0.00542},
        {7,15,  0.99813,  0.00628},
        {8,15,  0.99965,  0.00138},
        {8,17,  0.99965,  0.00156},
    };

    DWBATest dw;
    printf("%-9s %-6s  %-22s  %-22s  %7s  %7s\n",
           "Channel","L JP","C++ S","Ptolemy S","d|Re|","d|Im|");
    printf("%s\n", std::string(90,'-').c_str());

    // ── INCOMING ────────────────────────────────────────────────────────────
    Channel ch_in;
    setup_incoming(ch_in);
    dw.setup_and_run(ch_in, 0, 2);  // first call sets up the grid

    for (auto& ref : ptol_in) {
        setup_incoming(ch_in);
        dw.setup_and_run(ch_in, ref.L, ref.JP);
        std::complex<double> S = (ch_in.SMatrix.size() > (size_t)ref.L)
                                  ? ch_in.SMatrix[ref.L]
                                  : std::complex<double>(0,0);
        double dRe = S.real() - ref.Re;
        double dIm = S.imag() - ref.Im;
        printf("INCOMING  L=%-2d J=%2d/2  cpp=(%+7.4f,%+7.4f)  ptol=(%+7.4f,%+7.4f)  dRe=%+6.3f  dIm=%+6.3f%s\n",
               ref.L, ref.JP, S.real(), S.imag(), ref.Re, ref.Im,
               dRe, dIm,
               (std::abs(dRe)>0.05 || std::abs(dIm)>0.05) ? "  ❌" : "  ✅");
    }

    printf("%s\n", std::string(90,'-').c_str());

    // ── OUTGOING ────────────────────────────────────────────────────────────
    for (auto& ref : ptol_out) {
        Channel ch_out;
        setup_outgoing(ch_out);
        dw.setup_and_run(ch_out, ref.L, ref.JP);
        std::complex<double> S = (ch_out.SMatrix.size() > (size_t)ref.L)
                                  ? ch_out.SMatrix[ref.L]
                                  : std::complex<double>(0,0);
        double dRe = S.real() - ref.Re;
        double dIm = S.imag() - ref.Im;
        printf("OUTGOING  L=%-2d J=%2d/2  cpp=(%+7.4f,%+7.4f)  ptol=(%+7.4f,%+7.4f)  dRe=%+6.3f  dIm=%+6.3f%s\n",
               ref.L, ref.JP, S.real(), S.imag(), ref.Re, ref.Im,
               dRe, dIm,
               (std::abs(dRe)>0.05 || std::abs(dIm)>0.05) ? "  ❌" : "  ✅");
    }

    return 0;
}
