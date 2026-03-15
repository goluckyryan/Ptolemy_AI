// test_betcal_inject.cpp
// Inject Ptolemy's post-9J SMAG/SPHASE into C++ BETCAL/AMPCAL.
// Uses FLAT KOFFS indexing exactly matching the Fortran test_betcal.f90.
// If correct, should give 1.863 mb/sr at 0° and peak at 15°.

#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include "include/dwba.h"
#include "src/dwba/math_utils.cpp"

// Coulomb phases σ_L for incoming (d+33Si) and outgoing (p+34Si)
static double sigin[31] = {
    -0.293, 0.318, 0.654, 0.883, 1.056, 1.195, 1.311, 1.411,
     1.498, 1.576, 1.646, 1.709, 1.767, 1.821, 1.871, 1.918,
     1.961, 2.003, 2.041, 2.078, 2.113, 2.147, 2.180, 2.212,
     2.243, 2.273, 2.302, 2.331, 2.359, 2.387, 2.414
};
static double sigot[31] = {
    -0.224, 0.193, 0.412, 0.558, 0.669, 0.757, 0.831, 0.895,
     0.950, 0.999, 1.044, 1.084, 1.121, 1.155, 1.187, 1.216,
     1.244, 1.270, 1.295, 1.318, 1.340, 1.361, 1.381, 1.401,
     1.419, 1.437, 1.455, 1.472, 1.489, 1.506, 1.522
};

// JTOCS[koffs][0=dL, 1=LX, 2=2JP, 3=2JT] for 12 KOFFS (0-indexed)
static int JTOCS[12][4] = {
    { 0, 1, 1, 3},  // k=0  KOFFS=1: JP=1/2, LX=1, dL=0
    {-2, 2, 1, 3},  // k=1  KOFFS=2: JP=1/2, LX=2, dL=-2
    { 0, 2, 1, 3},  // k=2  KOFFS=3: JP=1/2, LX=2, dL=0
    { 2, 2, 1, 3},  // k=3  KOFFS=4: JP=1/2, LX=2, dL=+2
    { 0, 0, 3, 3},  // k=4  KOFFS=5: JP=3/2, LX=0, dL=0
    { 0, 1, 3, 3},  // k=5  KOFFS=6: JP=3/2, LX=1, dL=0
    {-2, 2, 3, 3},  // k=6  KOFFS=7: JP=3/2, LX=2, dL=-2
    { 0, 2, 3, 3},  // k=7  KOFFS=8: JP=3/2, LX=2, dL=0
    { 2, 2, 3, 3},  // k=8  KOFFS=9: JP=3/2, LX=2, dL=+2
    {-2, 3, 3, 3},  // k=9  KOFFS=10: JP=3/2, LX=3, dL=-2
    { 0, 3, 3, 3},  // k=10 KOFFS=11: JP=3/2, LX=3, dL=0
    { 2, 3, 3, 3},  // k=11 KOFFS=12: JP=3/2, LX=3, dL=+2
};
static const int NSPL = 12;

// SMAG[koffs][LI], SPHASE[koffs][LI] — all 0 unless set
static float SMAG[12][31];
static float SPHASE_ARR[12][31];

void fill_smag() {
    memset(SMAG, 0, sizeof(SMAG));
    memset(SPHASE_ARR, 0, sizeof(SPHASE_ARR));

    // KOFFS=1 (k=0): JP=1/2, LX=1, dL=0
    SMAG[0][1]=0.14585E-02f; SPHASE_ARR[0][1]=8.568f;
    SMAG[0][2]=0.22091E-02f; SPHASE_ARR[0][2]=7.504f;
    SMAG[0][3]=0.29519E-02f; SPHASE_ARR[0][3]=6.151f;
    SMAG[0][4]=0.25856E-02f; SPHASE_ARR[0][4]=4.715f;
    SMAG[0][5]=0.22830E-02f; SPHASE_ARR[0][5]=2.888f;
    SMAG[0][6]=0.17036E-02f; SPHASE_ARR[0][6]=1.352f;
    SMAG[0][7]=0.29488E-03f; SPHASE_ARR[0][7]=0.549f;
    SMAG[0][8]=0.33918E-03f; SPHASE_ARR[0][8]=0.670f;
    SMAG[0][9]=0.10178E-03f; SPHASE_ARR[0][9]=-0.072f;
    SMAG[0][10]=0.24836E-04f; SPHASE_ARR[0][10]=-0.343f;
    SMAG[0][11]=0.62382E-05f; SPHASE_ARR[0][11]=-0.439f;
    SMAG[0][12]=0.16177E-05f; SPHASE_ARR[0][12]=-0.475f;
    SMAG[0][13]=0.42746E-06f; SPHASE_ARR[0][13]=-0.488f;
    SMAG[0][14]=0.11410E-06f; SPHASE_ARR[0][14]=-0.494f;
    SMAG[0][15]=0.30628E-07f; SPHASE_ARR[0][15]=-0.497f;

    // KOFFS=2 (k=1): JP=1/2, LX=2, dL=-2 → LI=LO+2, LI min=2
    SMAG[1][2]=0.20288E-01f; SPHASE_ARR[1][2]=3.538f;
    SMAG[1][3]=0.28015E-01f; SPHASE_ARR[1][3]=2.172f;
    SMAG[1][4]=0.26046E-01f; SPHASE_ARR[1][4]=1.265f;
    SMAG[1][5]=0.27465E-01f; SPHASE_ARR[1][5]=-0.680f;
    SMAG[1][6]=0.30268E-01f; SPHASE_ARR[1][6]=-1.574f;
    SMAG[1][7]=0.27662E-01f; SPHASE_ARR[1][7]=-0.262f;
    SMAG[1][8]=0.47739E-01f; SPHASE_ARR[1][8]=-0.977f;
    SMAG[1][9]=0.29899E-01f; SPHASE_ARR[1][9]=-1.365f;
    SMAG[1][10]=0.16492E-01f; SPHASE_ARR[1][10]=-1.499f;
    SMAG[1][11]=0.90187E-02f; SPHASE_ARR[1][11]=-1.545f;
    SMAG[1][12]=0.49538E-02f; SPHASE_ARR[1][12]=-1.562f;
    SMAG[1][13]=0.27331E-02f; SPHASE_ARR[1][13]=-1.568f;
    SMAG[1][14]=0.15128E-02f; SPHASE_ARR[1][14]=-1.570f;
    SMAG[1][15]=0.83940E-03f; SPHASE_ARR[1][15]=-1.570f;

    // KOFFS=3 (k=2): JP=1/2, LX=2, dL=0
    SMAG[2][1]=0.28498E-01f; SPHASE_ARR[2][1]=0.325f;
    SMAG[2][2]=0.19032E-01f; SPHASE_ARR[2][2]=-0.455f;
    SMAG[2][3]=0.17730E-01f; SPHASE_ARR[2][3]=-2.011f;
    SMAG[2][4]=0.17251E-01f; SPHASE_ARR[2][4]=-2.645f;
    SMAG[2][5]=0.10069E-01f; SPHASE_ARR[2][5]=-0.754f;
    SMAG[2][6]=0.25254E-01f; SPHASE_ARR[2][6]=-1.007f;
    SMAG[2][7]=0.21121E-01f; SPHASE_ARR[2][7]=-1.198f;
    SMAG[2][8]=0.15263E-01f; SPHASE_ARR[2][8]=-1.316f;
    SMAG[2][9]=0.90015E-02f; SPHASE_ARR[2][9]=-1.463f;
    SMAG[2][10]=0.48676E-02f; SPHASE_ARR[2][10]=-1.529f;
    SMAG[2][11]=0.26025E-02f; SPHASE_ARR[2][11]=-1.555f;
    SMAG[2][12]=0.13955E-02f; SPHASE_ARR[2][12]=-1.565f;
    SMAG[2][13]=0.75197E-03f; SPHASE_ARR[2][13]=-1.568f;
    SMAG[2][14]=0.40714E-03f; SPHASE_ARR[2][14]=-1.570f;
    SMAG[2][15]=0.22140E-03f; SPHASE_ARR[2][15]=-1.570f;

    // KOFFS=4 (k=3): JP=1/2, LX=2, dL=+2
    SMAG[3][0]=0.31541E-01f; SPHASE_ARR[3][0]=-2.236f;
    SMAG[3][1]=0.21680E-01f; SPHASE_ARR[3][1]=-3.835f;
    SMAG[3][2]=0.17769E-01f; SPHASE_ARR[3][2]=-4.403f;
    SMAG[3][3]=0.62396E-02f; SPHASE_ARR[3][3]=-2.169f;
    SMAG[3][4]=0.15631E-01f; SPHASE_ARR[3][4]=-2.182f;
    SMAG[3][5]=0.11666E-01f; SPHASE_ARR[3][5]=-2.110f;
    SMAG[3][6]=0.71561E-02f; SPHASE_ARR[3][6]=-1.778f;
    SMAG[3][7]=0.52177E-02f; SPHASE_ARR[3][7]=-1.544f;
    SMAG[3][8]=0.34708E-02f; SPHASE_ARR[3][8]=-1.468f;
    SMAG[3][9]=0.21234E-02f; SPHASE_ARR[3][9]=-1.505f;
    SMAG[3][10]=0.11820E-02f; SPHASE_ARR[3][10]=-1.542f;
    SMAG[3][11]=0.64164E-03f; SPHASE_ARR[3][11]=-1.559f;
    SMAG[3][12]=0.34699E-03f; SPHASE_ARR[3][12]=-1.566f;
    SMAG[3][13]=0.18802E-03f; SPHASE_ARR[3][13]=-1.569f;

    // KOFFS=5 (k=4): JP=3/2, LX=0, dL=0
    SMAG[4][1]=0.42671E-03f; SPHASE_ARR[4][1]=1.476f;
    SMAG[4][2]=0.29893E-03f; SPHASE_ARR[4][2]=0.403f;
    SMAG[4][3]=0.15009E-03f; SPHASE_ARR[4][3]=-0.887f;
    SMAG[4][4]=0.86380E-04f; SPHASE_ARR[4][4]=-1.561f;
    SMAG[4][5]=0.67459E-04f; SPHASE_ARR[4][5]=-2.038f;

    // KOFFS=6 (k=5): JP=3/2, LX=1, dL=0
    SMAG[5][1]=0.79059E-02f; SPHASE_ARR[5][1]=3.226f;
    SMAG[5][2]=0.50619E-02f; SPHASE_ARR[5][2]=2.142f;
    SMAG[5][3]=0.35960E-02f; SPHASE_ARR[5][3]=1.102f;
    SMAG[5][4]=0.29087E-02f; SPHASE_ARR[5][4]=-0.182f;
    SMAG[5][5]=0.22698E-02f; SPHASE_ARR[5][5]=-1.266f;
    SMAG[5][6]=0.11283E-02f; SPHASE_ARR[5][6]=-0.870f;
    SMAG[5][7]=0.13218E-03f; SPHASE_ARR[5][7]=0.066f;
    SMAG[5][8]=0.14898E-03f; SPHASE_ARR[5][8]=0.742f;

    // KOFFS=7 (k=6): JP=3/2, LX=2, dL=-2
    SMAG[6][2]=0.56289E-02f; SPHASE_ARR[6][2]=3.516f;
    SMAG[6][3]=0.73001E-02f; SPHASE_ARR[6][3]=2.152f;
    SMAG[6][4]=0.66325E-02f; SPHASE_ARR[6][4]=1.237f;
    SMAG[6][5]=0.68413E-02f; SPHASE_ARR[6][5]=-0.731f;
    SMAG[6][6]=0.75254E-02f; SPHASE_ARR[6][6]=-1.615f;
    SMAG[6][7]=0.66609E-02f; SPHASE_ARR[6][7]=-0.315f;
    SMAG[6][8]=0.10891E-01f; SPHASE_ARR[6][8]=-1.001f;
    SMAG[6][9]=0.64736E-02f; SPHASE_ARR[6][9]=-1.382f;
    SMAG[6][10]=0.34165E-02f; SPHASE_ARR[6][10]=-1.510f;
    SMAG[6][11]=0.18032E-02f; SPHASE_ARR[6][11]=-1.552f;

    // KOFFS=8 (k=7): JP=3/2, LX=2, dL=0
    SMAG[7][1]=0.26918E-02f; SPHASE_ARR[7][1]=10.133f;
    SMAG[7][2]=0.22235E-02f; SPHASE_ARR[7][2]=9.175f;
    SMAG[7][3]=0.24132E-02f; SPHASE_ARR[7][3]=8.063f;
    SMAG[7][4]=0.20447E-02f; SPHASE_ARR[7][4]=6.598f;
    SMAG[7][5]=0.10321E-02f; SPHASE_ARR[7][5]=4.709f;
    SMAG[7][6]=0.28501E-03f; SPHASE_ARR[7][6]=2.559f;
    SMAG[7][7]=0.51383E-04f; SPHASE_ARR[7][7]=-0.105f;
    SMAG[7][8]=0.34257E-04f; SPHASE_ARR[7][8]=1.383f;
    SMAG[7][9]=0.91073E-05f; SPHASE_ARR[7][9]=0.166f;
    SMAG[7][10]=0.18984E-05f; SPHASE_ARR[7][10]=-0.243f;

    // KOFFS=9 (k=8): JP=3/2, LX=2, dL=+2
    SMAG[8][0]=0.58798E-02f; SPHASE_ARR[8][0]=1.401f;
    SMAG[8][1]=0.59054E-02f; SPHASE_ARR[8][1]=0.188f;
    SMAG[8][2]=0.46052E-02f; SPHASE_ARR[8][2]=-1.137f;
    SMAG[8][3]=0.22122E-02f; SPHASE_ARR[8][3]=-2.542f;
    SMAG[8][4]=0.61042E-03f; SPHASE_ARR[8][4]=-3.162f;
    SMAG[8][5]=0.15802E-03f; SPHASE_ARR[8][5]=-3.714f;
    SMAG[8][6]=0.89429E-04f; SPHASE_ARR[8][6]=-1.264f;
    SMAG[8][7]=0.58317E-04f; SPHASE_ARR[8][7]=-2.419f;
    SMAG[8][8]=0.26459E-04f; SPHASE_ARR[8][8]=-2.160f;

    // KOFFS=10 (k=9): JP=3/2, LX=3, dL=-2
    SMAG[9][3]=0.32545E-02f; SPHASE_ARR[9][3]=9.647f;
    SMAG[9][4]=0.56884E-02f; SPHASE_ARR[9][4]=8.497f;
    SMAG[9][5]=0.98185E-02f; SPHASE_ARR[9][5]=6.956f;
    SMAG[9][6]=0.10207E-01f; SPHASE_ARR[9][6]=5.053f;
    SMAG[9][7]=0.71335E-02f; SPHASE_ARR[9][7]=2.690f;
    SMAG[9][8]=0.19228E-02f; SPHASE_ARR[9][8]=0.615f;
    SMAG[9][9]=0.35070E-03f; SPHASE_ARR[9][9]=-0.356f;
    SMAG[9][10]=0.70986E-04f; SPHASE_ARR[9][10]=-0.722f;

    // KOFFS=11 (k=10): JP=3/2, LX=3, dL=0
    SMAG[10][2]=0.34132E-02f; SPHASE_ARR[10][2]=6.692f;
    SMAG[10][3]=0.57078E-02f; SPHASE_ARR[10][3]=5.433f;
    SMAG[10][4]=0.56089E-02f; SPHASE_ARR[10][4]=3.915f;
    SMAG[10][5]=0.36845E-02f; SPHASE_ARR[10][5]=2.225f;
    SMAG[10][6]=0.15595E-02f; SPHASE_ARR[10][6]=1.053f;
    SMAG[10][7]=0.30303E-03f; SPHASE_ARR[10][7]=0.457f;
    SMAG[10][8]=0.23996E-03f; SPHASE_ARR[10][8]=0.516f;

    // KOFFS=12 (k=11): JP=3/2, LX=3, dL=+2
    SMAG[11][1]=0.30682E-02f; SPHASE_ARR[11][1]=3.526f;
    SMAG[11][2]=0.29373E-02f; SPHASE_ARR[11][2]=2.158f;
    SMAG[11][3]=0.16742E-02f; SPHASE_ARR[11][3]=0.779f;
    SMAG[11][4]=0.63385E-03f; SPHASE_ARR[11][4]=0.151f;
    SMAG[11][5]=0.16273E-03f; SPHASE_ARR[11][5]=-0.658f;
    SMAG[11][6]=0.11928E-03f; SPHASE_ARR[11][6]=1.964f;
    SMAG[11][7]=0.69306E-04f; SPHASE_ARR[11][7]=0.764f;
    SMAG[11][8]=0.31249E-04f; SPHASE_ARR[11][8]=1.003f;
}

int main() {
    fill_smag();

    const int LMX = 30;
    const double UK_IN = 1.3082;
    const double FACTOR_BET = 0.5 / UK_IN;

    // BETAS[flat_koffs][Lo] — flat_koffs same as KOFFS-1 (0-indexed) from Fortran
    // Matches Fortran: BETAS(1:2, KOFFS, ILO) where ILO = Lo+1
    double BETAS_R[12][31] = {};  // real parts
    double BETAS_I[12][31] = {};  // imaginary parts

    // BETCAL — matching Fortran exactly
    int LIPREV;
    int MXZ, KOFFZ;

    for (int Lo = 0; Lo <= LMX; ++Lo) {
        int LIMN = Lo - 2;  // minimum LI (LDELMX = 2)

        // Compute TEMPS(I+MX) for this Lo — same as Fortran
        // TEMPS indexed by I+MX where I = 1 + NMX*(LX - LX_MN + NMLX*((LI-LIMN)/2))
        // For simplicity, compute directly within the KOFFS loop

        LIPREV = -100000;
        for (int k = 0; k < NSPL; ++k) {
            int dL = JTOCS[k][0];
            int Lx = JTOCS[k][1];
            int Li = Lo - dL;

            // New group when LI does not decrease
            if (Li >= LIPREV) {
                MXZ = Lx + dL;   // = LX + dL
                KOFFZ = k - MXZ; // index of MX=0 slot for this group
            }
            LIPREV = Li;

            if (Li < 0 || Li > LMX) continue;

            float amag = SMAG[k][Li];
            if (amag == 0.0f) continue;

            double phase_val = (double)SPHASE_ARR[k][Li] + sigin[Li] + sigot[Lo];
            double smatr =  (double)amag * std::sin(phase_val);
            double smati = -(double)amag * std::cos(phase_val);

            // CG factor for this (Li, Lx, Lo, Mx) combination
            for (int Mx_loop = std::max(0, MXZ); Mx_loop <= Lx; ++Mx_loop) {
                if (Lo < Mx_loop) continue;
                double cg = ClebschGordan((double)Li, 0.0, (double)Lx, (double)Mx_loop,
                                          (double)Lo, (double)Mx_loop);
                if (!std::isfinite(cg) || std::abs(cg) < 1e-14) continue;
                double TEMPS = FACTOR_BET * (2.0*Li + 1.0) * cg;

                // Store in slot KOFFZ + Mx_loop (matches Fortran BETAS[KOFFZ+MX, ILO])
                int slot = KOFFZ + Mx_loop;
                if (slot < 0 || slot >= NSPL) continue;
                BETAS_R[slot][Lo] += TEMPS * smatr;
                BETAS_I[slot][Lo] += TEMPS * smati;
            }
        }

        // Apply sqrt_factorial for each KOFFS slot
        for (int k = 0; k < NSPL; ++k) {
            int dL = JTOCS[k][0];
            int Lx = JTOCS[k][1];
            int Mx = (dL + Lx + 1) / 2;  // integer division — exact Fortran formula
            // sqrt_fac = prod_{n=1}^{Mx} 1/sqrt((Lo+n)(Lo-n+1))
            if (Mx > Lo) { BETAS_R[k][Lo] = 0.0; BETAS_I[k][Lo] = 0.0; continue; }
            double sf = 1.0;
            for (int n = 1; n <= Mx; ++n)
                sf /= std::sqrt((double)(Lo+n) * (double)(Lo-n+1));
            BETAS_R[k][Lo] *= sf;
            BETAS_I[k][Lo] *= sf;
        }
    }

    // AMPCAL — angular distribution
    printf("\nAngular distribution with injected Ptolemy SMAG (flat-KOFFS BETCAL):\n");
    printf("Angle   dSigma/dOmega  Ptolemy_ref  Ratio\n");
    double ptol_ref[] = {1.863, 1.905, 2.167, 2.535, 2.457, 1.759, 0.905};
    double test_angles[] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};

    for (int iang = 0; iang < 7; ++iang) {
        double theta_deg = test_angles[iang];
        double theta_rad = theta_deg * M_PI / 180.0;
        double cos_t = std::cos(theta_rad);

        double sigma = 0.0;
        for (int k = 0; k < NSPL; ++k) {
            int dL = JTOCS[k][0];
            int Lx = JTOCS[k][1];
            int Mx = (dL + Lx + 1) / 2;
            double FMNEG = (Mx == 0) ? 1.0 : 2.0;

            double Fr = 0.0, Fi = 0.0;
            for (int Lo = std::max(0, Mx); Lo <= LMX; ++Lo) {
                // PLM without Condon-Shortley phase: use assoc_legendre and remove (-1)^Mx
                double plm = std::assoc_legendre(Lo, Mx, cos_t);
                if (Mx % 2 != 0) plm = -plm;  // remove C-S phase from std library
                Fr += BETAS_R[k][Lo] * plm;
                Fi += BETAS_I[k][Lo] * plm;
            }
            sigma += FMNEG * (Fr*Fr + Fi*Fi);
        }
        sigma *= 10.0;
        printf("  %5.1f°  %8.4f mb/sr   (Ptolemy: %.4f)  ratio=%.4f\n",
               theta_deg, sigma, ptol_ref[iang], sigma/ptol_ref[iang]);
    }

    return 0;
}
