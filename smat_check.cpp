// smat_check.cpp — Compare C++ WavElj S-matrix vs Ptolemy for both channels
// Incoming: d+16O, Outgoing: p+17O
// Tests multiple L and J values

#include "dwba.h"
#include "potential_eval.h"
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
  DWBA dw;

  // ── Incoming channel: d + 16O ──────────────────────────────────────────
  // An-Cai 2006 OMP (same as fr_o16dp.cpp)
  // VSO=7.1140 (Ptolemy convention, includes ×2 for deuteron SO bug)
  {
    Channel &ch = dw.Incoming;
    ch.Projectile.A = 2; ch.Projectile.Z = 1; ch.Projectile.Name = "d";
    ch.Target.A = 16;    ch.Target.Z = 8;     ch.Target.Name = "16O";
    ch.JSPS = 2;  // 2*S for deuteron (S=1)
    ch.Elab = 20.0;

    // CM energy and kinematics
    double Ap = 2.0, At = 16.0;
    double mu = Ap * At / (Ap + At);
    double Ecm = At / (Ap + At) * ch.Elab;
    ch.mu = mu;
    double HBARC = 197.32697, AMU = 931.494;
    ch.k = std::sqrt(2.0 * mu * AMU * Ecm) / HBARC;
    ch.eta = ch.Projectile.Z * ch.Target.Z * mu * AMU / (137.036 * HBARC * ch.k);

    // Potentials (An-Cai 2006)
    ch.Pot.V    = 88.9546; ch.Pot.R0   = 1.1489; ch.Pot.A0   = 0.7508;
    ch.Pot.VI   = 2.3480;  ch.Pot.RI0  = 1.3446; ch.Pot.AI0  = 0.6030;
    ch.Pot.VSI  = 10.2180; ch.Pot.RSI0 = 1.3943; ch.Pot.ASI  = 0.6872;
    ch.Pot.VSO  = 7.1140;  ch.Pot.RSO0 = 0.9720; ch.Pot.ASO  = 1.0110;
    ch.Pot.VSOI = 0.0000;  ch.Pot.RC0  = 1.3030;

    // Scale radii by At^(1/3) — r0target convention
    double scale = std::cbrt(At);
    ch.Pot.R0   *= scale; ch.Pot.RI0  *= scale; ch.Pot.RSI0 *= scale;
    ch.Pot.RSO0 *= scale; ch.Pot.RC0  *= scale;

    dw.WavSet(ch);

    printf("=== INCOMING: d + 16O  Elab=20 MeV  Ecm=%.4f  k=%.5f  eta=%.5f ===\n",
           Ecm, ch.k, ch.eta);
    printf("%-8s %-6s  %12s %12s  %12s\n", "L", "JP", "Re(S)", "Im(S)", "|S|");

    // Test L=0,1,2,3,4 with all valid J
    for (int L = 0; L <= 4; ++L) {
      for (int JP = std::abs(2*L - 2); JP <= 2*L + 2; JP += 2) {
        if (JP < 0) continue;
        dw.WavElj(ch, L, JP);
        if (ch.SMatrix.size() > (size_t)L) {
          auto S = ch.SMatrix[L];
          printf("L=%-2d  JP=%2d/2  Re=%+10.5f  Im=%+10.5f  |S|=%.5f\n",
                 L, JP, S.real(), S.imag(), std::abs(S));
        }
      }
    }
    printf("\n");
  }

  // ── Outgoing channel: p + 17O ──────────────────────────────────────────
  // Koning-Delaroche OMP
  {
    Channel &ch = dw.Outgoing;
    ch.Projectile.A = 1; ch.Projectile.Z = 1; ch.Projectile.Name = "p";
    ch.Target.A = 17;    ch.Target.Z = 8;     ch.Target.Name = "17O";
    ch.JSPS = 1;  // 2*S for proton (S=1/2)

    // Q-value and kinematics
    // BE(17O g.s.) = 4.143 MeV → Q = 4.143 - 2.224 = 1.919 MeV (approx from fr_o16dp)
    // Use same Ecm_out as fr_o16dp.cpp
    double Ap = 1.0, At = 17.0;
    double mu = Ap * At / (Ap + At);
    // From fr_o16dp: Ecm_out computed from Q + Ecm_in
    // Ecm_in = 16/18 * 20 = 17.778; Q=1.919; Ecm_out = Ecm_in + Q = 19.697 (approx)
    // Use fr_o16dp value directly: kout = 0.94628 (from Ptolemy output)
    double HBARC = 197.32697, AMU = 931.494;
    // Reverse-engineer Ecm from k: Ecm = (k*HBARC)^2 / (2*mu*AMU)
    double k_out = 0.94628;  // from Ptolemy output
    ch.k   = k_out;
    double Ecm_out = k_out * k_out * HBARC * HBARC / (2.0 * mu * AMU);
    ch.mu  = mu;
    ch.eta = ch.Projectile.Z * ch.Target.Z * mu * AMU / (137.036 * HBARC * ch.k);

    // Potentials (Koning-Delaroche)
    ch.Pot.V    = 49.5434; ch.Pot.R0   = 1.1462; ch.Pot.A0   = 0.6753;
    ch.Pot.VI   = 2.0611;  ch.Pot.RI0  = 1.1462; ch.Pot.AI0  = 0.6753;
    ch.Pot.VSI  = 7.6703;  ch.Pot.RSI0 = 1.3016; ch.Pot.ASI  = 0.5275;
    ch.Pot.VSO  = 5.2956;  ch.Pot.RSO0 = 0.9338; ch.Pot.ASO  = 0.5900;
    ch.Pot.VSOI =-0.1059;  ch.Pot.RSOI0= 0.9338; ch.Pot.ASOI = 0.5900;
    ch.Pot.RC0  = 1.3030;

    // Scale radii by At^(1/3)
    double scale = std::cbrt(At);
    ch.Pot.R0   *= scale; ch.Pot.RI0  *= scale; ch.Pot.RSI0 *= scale;
    ch.Pot.RSO0 *= scale; ch.Pot.RSOI0*= scale; ch.Pot.RC0  *= scale;

    dw.WavSet(ch);

    printf("=== OUTGOING: p + 17O  Ecm=%.4f  k=%.5f  eta=%.5f ===\n",
           Ecm_out, ch.k, ch.eta);
    printf("%-8s %-6s  %12s %12s  %12s\n", "L", "JP", "Re(S)", "Im(S)", "|S|");

    for (int L = 0; L <= 8; ++L) {
      for (int JP = std::abs(2*L - 1); JP <= 2*L + 1; JP += 2) {
        if (JP < 0) continue;
        dw.WavElj(ch, L, JP);
        if (ch.SMatrix.size() > (size_t)L) {
          auto S = ch.SMatrix[L];
          printf("L=%-2d  JP=%2d/2  Re=%+10.5f  Im=%+10.5f  |S|=%.5f\n",
                 L, JP, S.real(), S.imag(), std::abs(S));
        }
      }
    }
  }

  return 0;
}
