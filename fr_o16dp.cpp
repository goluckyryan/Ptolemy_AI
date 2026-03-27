// fr_o16dp.cpp — Full FR DWBA for 16O(d,p)17O g.s. at Elab=20 MeV
// Compares to Ptolemy: 0°=41.629 mb/sr (no spectroscopic factor)
// Potentials: An-Cai 2006 (d+16O), Koning-Delaroche (p+17O)
// r0target convention throughout

#include "dwba.h"
#include <iostream>
#include <iomanip>

int main() {
    DWBA dwba;

    // ── Incoming: d + 16O at Elab=20 MeV ──
    // SetReaction(target, projectile, ejectile, recoil)
    // Internally: Outgoing.Projectile = Isotope(recoil), Outgoing.Target = Isotope(ejectile)
    // For (d,p): ejectile=proton(1H), recoil=17O
    // We need Outgoing.Projectile = proton → pass "1H" as recoil arg
    dwba.SetReaction("16O", "2H", "17O", "1H");
    dwba.SetEnergy(20.0);
    dwba.SetAngles(0.0, 180.0, 1.0);

    {
        ChannelPotential p = {};
        p.V    =  88.9546; p.R0   = 1.1489; p.A   = 0.7508;
        p.VI   =   2.3480; p.RI0  = 1.3446; p.AI  = 0.6030;
        p.VSI  =  10.2180; p.RSI0 = 1.3943; p.ASI = 0.6872;
        // VSO=7.1140 matches Ptolemy input. Ptolemy WOODSX type2 = 2*VSO/r*dWS/dr.
        // C++ V_so = 2*VSO/r*|dWS/dr|. SDOTL=<L·S> in both. So use same VSO as input.
        p.VSO  =   7.1140; p.RSO0 = 0.9720; p.ASO = 1.0110;
        p.VSOI =   0.0;    p.RSOI0= 0.9720; p.ASOI= 1.0110;
        p.RC0  =   1.3030;
        dwba.SetIncomingPotential(p);
    }

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

    // ── Target bound state: neutron in 17O g.s. (0d5/2, BE=4.143 MeV) ──
    {
        ChannelPotential bs = {};
        bs.V   = 60.0;  // initial guess; bisected to ~71.42 MeV
        bs.R0  =  1.10;
        bs.A   =  0.65;
        bs.RC0 =  1.30;
        dwba.SetTargetBoundState(0, 2, 2.5, 4.143, bs);
    }

    // ── Projectile bound state: neutron in deuteron (1s1/2, BE=2.224 MeV) ──
    {
        ChannelPotential bs = {};
        bs.V   = 72.0;
        bs.R0  =  1.00;
        bs.A   =  0.50;
        bs.RC0 =  0.00;
        dwba.SetProjectileBoundState(0, 0, 0.5, 2.224, bs);
    }

    // ── Nuclear spins ──
    dwba.SetTargetSpin(0.0);    // 16O: 0+
    dwba.SetResidualSpin(2.5);  // 17O g.s.: 5/2+
    dwba.SetLmin(0);            // lmin=0 (full partial wave sum)
    dwba.SetLmax(40);           // lmax=40

    // ── Run FR DWBA ──
    std::cout << "=== FR DWBA: 16O(d,p)17O g.s. at Elab=20 MeV ===\n";
    std::cout << "Ptolemy reference: 0 deg = 41.629 mb/sr\n\n";

    dwba.Calculate();

    // ── Print transfer S-matrix (before 9-J) for comparison with Ptolemy PRINT=2 output ──
    // Ptolemy format: L_IN  JP_IN  L_OUT  JP_OUT  LX   Re(S)   Im(S)   |S|
    std::cout << "\n=== Transfer S-matrix (before 9-J, C++) ===\n";
    std::cout << "   Li  JPI   Lo  JPO  Lx              Re(S)           Im(S)         |S|\n";
    for (auto& e : dwba.GetTransferSMatrix()) {
        double mag = std::abs(e.S);
        std::cout << std::setw(5) << e.Li
                  << std::setw(5) << e.JPI << "/2"
                  << std::setw(5) << e.Lo
                  << std::setw(5) << e.JPO << "/2"
                  << std::setw(4) << e.Lx
                  << "  " << std::setw(14) << std::scientific << std::setprecision(4) << e.S.real()
                  << "  " << std::setw(14) << e.S.imag()
                  << "  " << std::setw(14) << mag << "\n";
    }

    return 0;
}
