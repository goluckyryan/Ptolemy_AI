// a12_validate.cpp — Compare C++ ComputeA12Terms vs Ptolemy print=60000 ANSWER lines
// For 16O(d,p)17O: lT=2, lP=0, Lx=2

#include "dwba.h"
#include "math_utils.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <tuple>

class DWBATest : public DWBA {
public:
    std::vector<std::tuple<int,int,double>> testA12(int Li,int Lo,int Lx,int lT,int lP){
        return ComputeA12Terms(Li,Lo,Lx,lT,lP);
    }
};

// Ptolemy ANSWER format: MP, MT, MU, LX, LO, LOMNMN, XLOTMP, OUTTMP, TTT, TEMP
// TEMP = XLOTMP * OUTTMP * TTT * doubling   (the A12 coefficient stored)
// For lP=0: MP=0 always. 
// Our ComputeA12Terms stores (MT, MU, coeff) where coeff = sum over (MT,MP) for given (MT,MU)
// For lP=0: MP=0 so MT=MX. Our code uses MX=MT (since MP=0).

// Parse Ptolemy ANSWER lines for Li=0 (MUSTRT=0 → MU only 0)
// and compare to our A12 terms

struct PtolemyA12 {
    int MP, MT, MU, LX, LO;
    double XLOTMP, OUTTMP, TTT, TEMP;
};

int main() {
    DWBATest dw;
    int lT=2, lP=0, Lx=2;

    printf("=== A12 Validation: lT=%d, lP=%d, Lx=%d ===\n\n", lT, lP, Lx);
    printf("Ptolemy ANSWER format: MP, MT, MU, LX, LO, LOMNMN, XLOTMP, OUTTMP, TTT, TEMP\n");
    printf("Our A12: (MT, MU, coeff) where coeff should equal TEMP (for matching LX,LO)\n\n");

    // Test cases: (Li, Lo) pairs
    // For lP=0, Lx=2: triangle requires |Li-Lo|<=2<=Li+Lo
    // and (Li+Lo+Lx) even → Li+Lo even → Li and Lo same parity

    struct TestCase { int Li, Lo; };
    std::vector<TestCase> cases = {
        {0, 2},   // main case for DWBA
        {2, 0},
        {2, 2},
        {2, 4},
        {4, 2},
        {4, 4},
    };

    // Ptolemy ANSWER lines from print=60000 for Li=0 (from START line: LBT=2, LBP=0, LO-range, LI=?)
    // First START line: "START 2 0 0 2 2 0 4 6 -2 0 4 15 3"
    // LBT=2, LBP=0, [something]=0, lxmin=2, lxmax=2, lmin=0, lmax=4, lomost=6
    // The first ANSWER block corresponds to LI=0 (from MUSTRT=0, MU starts at 0)
    // Li=0: MU=0 only
    // MT=-2: ANSWER  0 -2 0  2 2 2  2 (XLOTMP=1.3693, OUTTMP=0.30619, TTT=0.44721, TEMP=0.18750)
    // MT=0:  ANSWER  0  0 0  2 2 2  2 (XLOTMP=-1.1180, OUTTMP=-0.25000, TTT=0.44721, TEMP=0.12500)  
    //   Wait: TEMP=XLOTMP*OUTTMP*TTT = (-1.1180)*(-0.25000)*(0.44721) = 0.1250 but that's +0.125
    //   But the doubling for MU=0 in HALFSW=FALSE case: NO doubling when MU=0
    // MT=2: ANSWER  0  2 0  2 2 2  2 (XLOTMP=1.3693, OUTTMP=0.30619, TTT=0.44721, TEMP=0.18750)

    // So for Li=0, Lo=2, Lx=2: A12 terms are:
    //   (MT=-2, MU=0, coeff=0.18750)
    //   (MT=0,  MU=0, coeff=0.12500)  [note: negative TEMP? Let me check]
    //   (MT=2,  MU=0, coeff=0.18750)

    // Actually TEMP = XLOTMP * OUTTMP * TTT (then doubled if MU!=0 or HALFSW)
    // For Li=0: MU=0, HALFSW=FALSE → no doubling → TEMP = XLOTMP*OUTTMP*TTT
    // MT=0: TEMP = (-1.1180)*(-0.25)*0.44721 = 0.1250  ... but it says 0.12500 ✓

    // From the ANSWER lines for Li=0 block (first START block):
    // These correspond to Li being iterated inside A12; the first block is LI from LMIN

    printf("--- Li=0, Lo=2, Lx=2 ---\n");
    auto terms = dw.testA12(0, 2, 2, 2, 0);
    printf("  C++ A12 terms (MT, MU, coeff):\n");
    for (auto& [mt, mu, c] : terms)
        printf("    MT=%+d  MU=%d  coeff=%+.6f\n", mt, mu, c);

    // Ptolemy reference (from ANSWER lines for Li=0, Lo=2, Lx=2):
    printf("  Ptolemy ANSWER (MT, MU → TEMP for LX=2, LO=2):\n");
    printf("    MT=-2  MU=0  TEMP=+0.18750\n");
    printf("    MT= 0  MU=0  TEMP=+0.12500  (note: XLOTMP=-1.1180, OUTTMP=-0.25, TTT=0.44721)\n");
    printf("    MT=+2  MU=0  TEMP=+0.18750\n\n");

    printf("--- Li=2, Lo=0, Lx=2 ---\n");
    auto terms20 = dw.testA12(2, 0, 2, 2, 0);
    printf("  C++ A12 terms:\n");
    for (auto& [mt, mu, c] : terms20)
        printf("    MT=%+d  MU=%d  coeff=%+.6f\n", mt, mu, c);
    printf("  (Ptolemy: need to identify Li=2 block)\n\n");

    printf("--- Li=2, Lo=2, Lx=2 ---\n");
    auto terms22 = dw.testA12(2, 2, 2, 2, 0);
    printf("  C++ A12 terms:\n");
    for (auto& [mt, mu, c] : terms22)
        printf("    MT=%+d  MU=%d  coeff=%+.6f\n", mt, mu, c);

    // From Ptolemy for Li=2, Lo=2, Lx=2 (MU=0 and MU=2 contributions):
    // MU=0: MT=-2 ANSWER: XLOTMP=0.68465, OUTTMP=0.68465, TTT=0.23905, TEMP=-0.11205
    //         MT=0        XLOTMP=0.55902, OUTTMP=-0.55902, TTT=-0.23905, TEMP=+0.074702? 
    //   Hmm need to be careful
    // MU=2: MT=-2 → TEMP=0.18792e-01; MT=0 → TEMP=0.062639; MT=2 → TEMP=0.18792e-01
    printf("  Ptolemy ANSWER for Li=2, Lo=2, Lx=2:\n");
    printf("    MU=0: MT=-2: -0.11205, MT=0: +0.07470, MT=+2: -0.11205\n");
    printf("    MU=2: MT=-2: +0.01879×2, MT=0: +0.06264×2, MT=+2: +0.01879×2  (doubled)\n\n");

    printf("--- Li=0, Lo=2, Lx=2 (detailed) ---\n");
    printf("Check EvalA12 at phi_T=0, phi_ab=0:\n");
    auto t = dw.testA12(0, 2, 2, 2, 0);
    double a12_zero = 0.0;
    for (auto& [mt, mu, c] : t) a12_zero += c;  // cos(0)=1 for all terms
    printf("  C++  A12(0,0) = %.6f\n", a12_zero);
    printf("  Pto  A12(0,0) = 0.18750+0.12500+0.18750 = %.6f\n", 0.18750+0.12500+0.18750);

    return 0;
}
