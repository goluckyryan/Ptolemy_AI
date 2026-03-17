// smat_compare.cpp — Compare WavElj S-matrix vs Ptolemy print=1001
// for INCOMING (d+16O) and OUTGOING (p+17O) channels, L=0..8, all valid J.
//
// Uses DWBA::TestWavElj() public wrapper — see dwba.h
//
// Build:
//   g++ -O2 -std=c++17 -Iinclude smat_compare.cpp \
//       src/dwba/wavelj.cpp src/dwba/potential_eval.cpp \
//       src/dwba/rcwfn.cpp src/dwba/math_utils.cpp \
//       src/input/Isotope.cpp src/input/Potentials.cpp \
//       src/dwba/setup.cpp src/dwba/av18_potential.cpp \
//       src/dwba/bound.cpp src/dwba/ineldc_zr.cpp \
//       src/dwba/grdset.cpp src/dwba/ineldc.cpp src/dwba/a12.cpp \
//       src/dwba/xsectn.cpp src/dwba/wavelj.cpp src/elastic/elastic.cpp \
//       -o smat_compare -lm

#include "dwba.h"
#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>
#include <string>

// Ptolemy reference values from print=1001 (lmax=8, both channels)
struct RefVal { int L, JP2; double Re, Im; };

// INCOMING: d + 16O (JSPS=2, deuteron)
// NOTE: fr_o16dp.cpp uses VSO=3.557 (physical, no /2S bug)
// BUT Ptolemy uses VSO=7.114 in the input file.
// For this S-matrix test, we use VSO=7.114 to match Ptolemy's elastic S-matrix.
static const RefVal ref_in[] = {
  {0, 2,  0.14429,   0.11977},
  {1, 0,  0.14914,  -0.14607},
  {1, 2,  0.18217,  -0.11199},
  {1, 4,  0.21878,  -0.023624},
  {2, 2, -0.096864, -0.12363},
  {2, 4, -0.015389, -0.15313},
  {2, 6,  0.10997,  -0.12204},
  {3, 4, -0.18101,   0.029071},
  {3, 6, -0.19169,  -0.071342},
  {3, 8, -0.097105, -0.17963},
  {4, 6,  0.11158,  -0.0081513},
  {4, 8,  0.0033604, 0.034952},
  {4,10, -0.10824,  -0.049103},
  {6,10,  0.52565,   0.097773},
  {6,12,  0.44438,   0.11223},
  {6,14,  0.39504,   0.091358},
  {8,14,  0.94380,   0.030498},
  {8,16,  0.93700,   0.067879},
  {8,18,  0.92588,   0.11200},
};

// OUTGOING: p + 17O (JSPS=1, proton)
static const RefVal ref_out[] = {
  {0, 1,  0.49629,   -0.019521},
  {1, 1,  0.028301,  -0.40199},
  {1, 3,  0.22806,   -0.30539},
  {2, 3, -0.32679,   -0.31757},
  {2, 5, -0.040840,  -0.54572},
  {3, 5,  0.082832,   0.14530},
  {3, 7, -0.079999,  -0.13415},
  {4, 7,  0.71811,    0.20865},
  {4, 9,  0.61227,    0.25849},
  {5, 9,  0.94558,    0.076015},
  {5,11,  0.93736,    0.096934},
  {6,11,  0.99003,    0.020914},
  {6,13,  0.98956,    0.025171},
  {7,13,  0.99814,    0.0054184},
  {7,15,  0.99813,    0.0062810},
  {8,15,  0.99965,    0.0013792},
  {8,17,  0.99965,    0.0015580},
};

int main() {
  DWBA dwba;

  // Setup the full reaction so kinematics are computed correctly
  dwba.SetReaction("16O", "2H", "17O", "1H");
  dwba.SetEnergy(20.0);

  // ── Incoming potential (d+16O) — use Ptolemy's input VSO=7.114 for S-matrix match
  {
    ChannelPotential p = {};
    p.V    =  88.9546; p.R0   = 1.1489; p.A   = 0.7508;
    p.VI   =   2.3480; p.RI0  = 1.3446; p.AI  = 0.6030;
    p.VSI  =  10.2180; p.RSI0 = 1.3943; p.ASI = 0.6872;
    p.VSO  =   7.1140; p.RSO0 = 0.9720; p.ASO = 1.0110;  // Ptolemy input value
    p.VSOI =   0.0;    p.RSOI0= 0.9720; p.ASOI= 1.0110;
    p.RC0  =   1.3030;
    dwba.SetIncomingPotential(p);
  }

  // ── Outgoing potential (p+17O)
  {
    ChannelPotential p = {};
    p.V    =  49.5434; p.R0   = 1.1462; p.A   = 0.6753;
    p.VI   =   2.0611; p.RI0  = 1.1462; p.AI  = 0.6753;
    p.VSI  =   7.6703; p.RSI0 = 1.3016; p.ASI = 0.5275;
    p.VSO  =   5.2956; p.RSO0 = 0.9338; p.ASO = 0.5900;
    p.VSOI =  -0.1059; p.RSOI0= 0.9338; p.ASOI= 0.5900;
    p.RC0  =   1.3030;
    dwba.SetOutgoingPotential(p);
  }

  // Dummy bound states (required by SetReaction)
  { ChannelPotential bs={}; bs.V=60; bs.R0=1.1; bs.A=0.65; bs.RC0=1.3;
    dwba.SetTargetBoundState(0,2,2.5,4.143,bs); }
  { ChannelPotential bs={}; bs.V=72; bs.R0=1.0; bs.A=0.5;
    dwba.SetProjectileBoundState(0,0,0.5,2.224,bs); }
  dwba.SetTargetSpin(0.0);
  dwba.SetResidualSpin(2.5);

  // Kinematics mode: false = NR (Ptolemy default), true = relativistic
  // dwba.SetRelativisticKinematics(true);  // uncomment to use relativistic
  dwba.SetupChannels();

  // ── INCOMING channel ──
  printf("=== INCOMING: d + 16O, Elab=20 MeV (VSO=7.114 as in Ptolemy input) ===\n");
  printf("  (Note: fr_o16dp uses VSO=3.557 to correct Ptolemy SO bug — comparing raw input here)\n");
  printf("%-4s %-5s  %-28s  %-28s  %s\n",
         "L","JP","C++ S","Ptolemy S","Diff");
  for (const auto& r : ref_in) {
    auto S = dwba.TestWavElj(true, r.L, r.JP2);
    double dRe = S.real()-r.Re, dIm = S.imag()-r.Im;
    double dS = std::sqrt(dRe*dRe+dIm*dIm);
    printf("L=%-2d JP=%2d/2  C++=(%+.5f,%+.5f)  Pto=(%+.5f,%+.5f)  |ΔS|=%.4f %s\n",
           r.L, r.JP2, S.real(), S.imag(), r.Re, r.Im, dS, dS>0.02?"❌":"✅");
  }

  // ── OUTGOING channel ──
  printf("\n=== OUTGOING: p + 17O ===\n");
  printf("%-4s %-5s  %-28s  %-28s  %s\n",
         "L","JP","C++ S","Ptolemy S","Diff");
  for (const auto& r : ref_out) {
    auto S = dwba.TestWavElj(false, r.L, r.JP2);
    double dRe = S.real()-r.Re, dIm = S.imag()-r.Im;
    double dS = std::sqrt(dRe*dRe+dIm*dIm);
    printf("L=%-2d JP=%2d/2  C++=(%+.5f,%+.5f)  Pto=(%+.5f,%+.5f)  |ΔS|=%.4f %s\n",
           r.L, r.JP2, S.real(), S.imag(), r.Re, r.Im, dS, dS>0.02?"❌":"✅");
  }

  return 0;
}
