#include <cstdio>
#include "include/math_utils.h"

int main() {
    // Fortran inner loop: DO JP=0,JPMX=5,2 → JP=0,2,4
    // For each JP, triangle constraint: LXMN2 = max(|JBT-JP|/2, |LO-LI|)
    //                                   LXMX2 = min((JBT+JP)/2, LO+LI)
    // LI=0, LO=2: |LO-LI|=2, LO+LI=2 → LXMN2=2, LXMX2=2 → only LX=2
    // JP=0: LXMN2=max((5-0)/2,2)=max(2,2)=2  LXMX2=min((5+0)/2,2)=min(2,2)=2 ✓
    // JP=2: LXMN2=max((5-2)/2,2)=max(1,2)=2  LXMX2=min((5+2)/2,2)=min(3,2)=2 ✓
    // JP=4: LXMN2=max((5-4)/2,2)=max(0,2)=2  LXMX2=min((5+4)/2,2)=min(4,2)=2 ✓
    // So LX=2 always selected. Now test 2nd 9-J for each JP:
    int JBT=5, JBP=1, JA=2, JB=1, LXP=2, LI=0, LO=2, JPI=2, JPO=3;
    printf("LI=0,LO=2,JPI=2,JPO=3: 2nd 9-J for each JP (LX=2):\n");
    for (int JP=0; JP<=5; JP+=2) {
        double w2 = NineJ(JBT/2.0, 2.0, JP/2.0,
                          JPI/2.0, (double)LI, JA/2.0,
                          JPO/2.0, (double)LO, JB/2.0);
        printf("  JP=%d: w2=%.8f\n", JP, w2);
    }
    // Also try JPO=5
    JPO=5;
    printf("LI=0,LO=2,JPI=2,JPO=5:\n");
    for (int JP=0; JP<=5; JP+=2) {
        double w2 = NineJ(JBT/2.0, 2.0, JP/2.0,
                          JPI/2.0, (double)LI, JA/2.0,
                          JPO/2.0, (double)LO, JB/2.0);
        printf("  JP=%d: w2=%.8f\n", JP, w2);
    }
    return 0;
}
