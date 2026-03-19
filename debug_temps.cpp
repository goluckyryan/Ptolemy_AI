#include <cmath>
#include <cstdio>
#include <cstdint>
#include <complex>
#include <vector>
#include "include/math_utils.h"

int main() {
    // Trace BETCAL TEMPS for LO=0, LI=2, LX=2, MX=0,1,2
    int LMN=0, LMX=20, LXMN=2, LXMX=2, LDELMX=2;
    int NMX=3, NMLX=1;
    double UK=1.23280, FACTOR=0.5/UK;
    
    int lo=0, LIMN=lo-LDELMX; // LIMN=-2
    
    // LI=0 case
    for (int li=LIMN; li<=lo+LDELMX; li+=2) {
        if (li < LMN) continue;
        if (li > LMX) continue;
        int LX1=std::max(std::abs(li-lo), LXMN);
        for (int lx=LX1; lx<=LXMX; lx++) {
            int mxz=(lx+li-lo)%2; if(mxz<0) mxz+=2;
            int base = NMX*(lx-LXMN + NMLX*(li-LIMN)/2);
            for (int mx=mxz; mx<=lx; mx++) {
                double cg = ClebschGordan((double)li,0.0,(double)lx,(double)mx,(double)lo,(double)mx);
                double t = FACTOR*(2*li+1)*cg;
                printf("LI=%d LX=%d MX=%d: base=%d idx=%d CG=%.8f TEMPS=%.8f\n",
                       li,lx,mx,base,base+mx,cg,t);
            }
        }
    }
    return 0;
}
