// lsqpol_test.cpp — standalone LSQPOL + MATINV C++ implementation
// Compare against Fortran LSQPOL output via a reference data file.
// Also runs with random inputs to check self-consistency.
//
// Usage:
//   ./lsqpol_test                 — run random tests only
//   ./lsqpol_test --ref <file>    — also compare vs Fortran reference

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <algorithm>
#include <string>
#include <random>

// ─── MATINV (faithful port of Fortran MATINV in fortlib.mor) ──────────────
// Solves A*X=B in-place (A → A^-1, B → X).
// A: N×N,  B: N×M  (Fortran column-major → here row-major with NMAX stride).
// NOTE: Fortran MATINV stores rows as the FIRST index (A(I,J) = row I, col J).
// Our C++ uses row-major: A[i][j] = row i, col j. Same layout.
//
// Key quirk: row swap during pivot, THEN column swap in post-loop (lines 600-710).
static void matinv(std::vector<std::vector<double>>& A, int N,
                   std::vector<std::vector<double>>& B, int M_rhs) {
    // PIVOT array: tracks which COLUMN has been used as pivot
    std::vector<double> pivot(N, 0.0);
    // INDEX array: INDEX[I] = IROW*4096 + ICOLUM (packed)
    std::vector<int> index_arr(N, 0);

    for (int I = 0; I < N; ++I) {
        // Search for pivot: max |A[j][k]| where pivot[j]==0 AND pivot[k]==0
        double amax = 0.0;
        int irow = 0, icolum = 0;
        for (int j = 0; j < N; ++j) {
            if (pivot[j] != 0.0) continue;
            for (int k = 0; k < N; ++k) {
                if (pivot[k] != 0.0) continue;
                double tmp = std::abs(A[j][k]);
                if (tmp >= amax) { amax = tmp; irow = j; icolum = k; }
            }
        }
        // Pack: INDEX[I] = IROW*4096 + ICOLUM
        index_arr[I] = irow * 4096 + icolum;
        pivot[icolum] = amax;   // mark column used (non-zero)

        // Interchange rows IROW <-> ICOLUM (puts pivot on diagonal)
        if (irow != icolum) {
            std::swap(A[irow], A[icolum]);
            if (M_rhs > 0) {
                for (int k = 0; k < M_rhs; ++k)
                    std::swap(B[irow][k], B[icolum][k]);
            }
        }

        // Divide pivot row by pivot element
        double pivot_val = A[icolum][icolum];
        A[icolum][icolum] = 1.0;
        for (int k = 0; k < N; ++k)       A[icolum][k] /= pivot_val;
        for (int k = 0; k < M_rhs; ++k)   B[icolum][k] /= pivot_val;

        // Eliminate all other rows
        for (int j = 0; j < N; ++j) {
            if (j == icolum) continue;
            double T = A[j][icolum];
            A[j][icolum] = 0.0;
            for (int k = 0; k < N; ++k)     A[j][k] -= A[icolum][k] * T;
            for (int k = 0; k < M_rhs; ++k) B[j][k] -= B[icolum][k] * T;
        }
    }

    // POST: interchange COLUMNS of A in reverse order (Fortran lines 600-710)
    // Fortran: DO I1=N,1,-1: K=INDEX(I1)/4096; ICOLUM=INDEX(I1)-4096*K
    //          swap A(*,K) <-> A(*,ICOLUM)
    for (int I1 = N - 1; I1 >= 0; --I1) {
        int K      = index_arr[I1] / 4096;
        int ICOLUM = index_arr[I1] % 4096;
        if (K != ICOLUM) {
            for (int j = 0; j < N; ++j)
                std::swap(A[j][K], A[j][ICOLUM]);
        }
    }
    // Note: B is already solved — column-swap on A does NOT apply to B.
    // The Fortran MATINV column swap reconstructs A^-1 properly, but
    // B[i][j] already contains the correct solution X after the row-ops.
}

// ─── LSQPOL (faithful port of Fortran LSQPOL in fortlib.mor) ─────────────
// Fits MSUB-term polynomial to LSUB right-hand-side vectors simultaneously.
//
// X[0..NSUB-1]            — x values
// Y[0..NSUB-1][0..LSUB-1] — y values (NSUB points, LSUB columns)
// W[0..NSUB-1]            — weights
// B_out[0..MSUB-1][0..LSUB-1] — output polynomial coefficients (x^0..x^(MSUB-1))
// RESID_out[0..NSUB-1][0..LSUB-1] — residuals: poly(X[k]) - Y[k][j]
// SUM_out[0..LSUB-1]      — weighted sum of squared residuals
static void lsqpol(const std::vector<double>& X,
                   const std::vector<std::vector<double>>& Y,  // [NSUB][LSUB]
                   const std::vector<double>& W,
                   int NSUB, int LSUB, int MSUB,
                   std::vector<std::vector<double>>& B_out,    // [MSUB][LSUB]
                   std::vector<std::vector<double>>& RESID_out, // [NSUB][LSUB]
                   std::vector<double>& SUM_out)
{
    int N = NSUB, L = LSUB, M = MSUB;

    // Scale X to [-1,1] to prevent overflow (Fortran LSQPOL does this)
    double Xmax = 0.0;
    for (int k = 0; k < N; ++k) Xmax = std::max(Xmax, std::abs(X[k]));
    if (Xmax < 1e-15) Xmax = 1.0;

    std::vector<double> Xs(N);
    for (int k = 0; k < N; ++k) Xs[k] = X[k] / Xmax;

    // XPOWER array (Fortran /F402/ common): holds running power sums
    // Indices used: M..3M-1 for A, 3M..4M-1 for B (0-indexed)
    // Fortran: XPOWER(M1..M31), M1=M+1, M3=2M+M=3M, M31=M3-1, M41=M31+M=4M-1
    // In 0-based C++: indices M..3M-1 for A matrix sums, 3M..4M-1 for B sums
    const int XPOW_SZ = 4 * M + 2;
    std::vector<double> xpow(XPOW_SZ, 0.0);

    // Formation of normal equations (A matrix power sums)
    // Fortran: DO K1=1,N: TERM=W(K1); DO K2=M1,M31: xp(K2)+=TERM; TERM*=X(K1)
    // In 0-based: K2 = M..3M-2 (inclusive), accumulating M+1..3M-1 powers of X
    // But Fortran M1=M+1, M31=3M-1 (all 1-based), so 0-based: M..3M-2
    for (int k = 0; k < N; ++k) {
        double term = W[k];
        for (int kk = M; kk <= 3*M - 2; ++kk) {
            xpow[kk] += term;
            term *= Xs[k];
        }
    }

    // Build A matrix: A(I,J) = XPOWER[I+J+M-1]  (1-based Fortran)
    // 0-based C++: A[i][j] = xpow[i+j+M]  (since i,j are 0-based, equivalent to (i+1)+(j+1)+M-1 = i+j+M+1, but Fortran 1-based I+J+M-1 → 0-based i+j+M-1+2 = i+j+M+1 ? )
    // Let's be careful:
    // Fortran (1-based): A(I,J) = XPOWER(I+J+M-1), I=1..M, J=1..M
    //   when I=1,J=1: XPOWER(1+1+M-1) = XPOWER(M+1) → index M (0-based)
    //   when I=M,J=M: XPOWER(2M+M-1) = XPOWER(3M-1) → index 3M-2 (0-based)
    // 0-based C++: A[i][j] (i=0..M-1, j=0..M-1)
    //   Fortran I=i+1, J=j+1: xpow[(i+1)+(j+1)+M-1 - 1] = xpow[i+j+M]
    std::vector<std::vector<double>> A(M, std::vector<double>(M, 0.0));
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            A[i][j] = xpow[i + j + M];

    // Build B matrix (RHS): B(I,J) = Σ_k W[k]*Y[k][J]*X[k]^I
    // Fortran: DO J=1,L: clear xpow[3M..4M-1]; DO K=1,N: TERM=W(K)*Y(K,J);
    //          DO K2=3M,4M-1: xpow[K2]+=TERM; TERM*=X(K)
    //          B(I,J) = xpow[I+3M-1] (1-based I)
    // 0-based: B[i][j] = xpow[3M + i]
    std::vector<std::vector<double>> B_mat(M, std::vector<double>(L, 0.0));
    for (int j = 0; j < L; ++j) {
        for (int kk = 3*M; kk <= 4*M - 1; ++kk) xpow[kk] = 0.0;
        for (int k = 0; k < N; ++k) {
            double term = W[k] * Y[k][j];
            for (int kk = 3*M; kk <= 4*M - 1; ++kk) {
                xpow[kk] += term;
                term *= Xs[k];
            }
        }
        for (int i = 0; i < M; ++i)
            B_mat[i][j] = xpow[3*M + i];
    }

    // Solve A * coeff = B_mat via MATINV
    matinv(A, M, B_mat, L);
    // B_mat now contains coefficients: B_mat[i][j] = coeff of x^i for column j

    // Un-scale coefficients: undo the Xscale on X
    // Fortran: DO J=1,L; DO I2=2,M; DO I=I2,M: B(I,J) *= XSCALE
    // (applying XSCALE=1/Xmax once per I2-I nesting level)
    // This correctly transforms from scaled-x polynomial to unscaled:
    // if p(xs) = Σ b_i * xs^i and xs = x/Xmax,
    // then b_i → b_i / Xmax^i, achieved by the triangular loop.
    double Xscale = 1.0 / Xmax;
    for (int j = 0; j < L; ++j)
        for (int i2 = 1; i2 < M; ++i2)   // I2=2..M (1-based) → 1..M-1 (0-based)
            for (int i = i2; i < M; ++i)   // I=I2..M (1-based) → i2..M-1 (0-based)
                B_mat[i][j] *= Xscale;

    // Evaluate residuals: RESID[k][j] = poly(Xs[k]) - Y[k][j]
    // Horner (Fortran): POLY=0; DO I2=1,M: I=M+1-I2; POLY=X[k]*POLY+B(I,J)
    // 0-based: POLY=0; for i2=0..M-1: i=M-1-i2; POLY=xs*POLY+B[M-1-i2][j]
    // = B[M-1][j] * xs^0 + B[M-2][j]*xs^1 + ... — wait, unroll:
    // i2=0: i=M-1; POLY=xs*0+B[M-1][j] = B[M-1][j]
    // i2=1: i=M-2; POLY=xs*B[M-1][j]+B[M-2][j]
    // → POLY = B[0][j] + B[1][j]*xs + ... + B[M-1][j]*xs^(M-1)
    // BUT at residual eval we use the SCALED xs, not unscaled x!
    // Fortran evaluates RESID using X[k] (which is STILL SCALED at this point —
    // the un-scale of X happens AFTER residual computation).
    RESID_out.assign(N, std::vector<double>(L, 0.0));
    SUM_out.assign(L, 0.0);
    for (int j = 0; j < L; ++j) {
        for (int k = 0; k < N; ++k) {
            double xs = Xs[k];   // still scaled
            double poly = 0.0;
            for (int i2 = 0; i2 < M; ++i2)
                poly = xs * poly + B_mat[M-1-i2][j];
            // poly is evaluated with SCALED coefficients (before un-scale)
            // But wait: B_mat was already un-scaled above!
            // So we need to re-evaluate with unscaled x:
            double x_unscaled = X[k];
            poly = 0.0;
            for (int i2 = 0; i2 < M; ++i2)
                poly = x_unscaled * poly + B_mat[M-1-i2][j];
            RESID_out[k][j] = poly - Y[k][j];
            SUM_out[j] += W[k] * RESID_out[k][j] * RESID_out[k][j];
        }
    }

    // Copy coefficients to output
    B_out = B_mat;
}

// ─── helpers ──────────────────────────────────────────────────────────────
static double poly_eval(const std::vector<std::vector<double>>& B, int j, double x) {
    int M = (int)B.size();
    double p = 0.0;
    for (int i2 = 0; i2 < M; ++i2)
        p = x * p + B[M-1-i2][j];
    return p;
}

// ─── Fortran reference interface ──────────────────────────────────────────
// We'll write a Fortran driver, compile and run it, then compare.
static bool write_fortran_driver(const std::vector<double>& X,
                                  const std::vector<std::vector<double>>& Y,
                                  const std::vector<double>& W,
                                  int NSUB, int LSUB, int MSUB,
                                  const char* fortran_src, const char* data_file) {
    // Write data file
    FILE* fp = fopen(data_file, "w");
    if (!fp) return false;
    fprintf(fp, "%d %d %d\n", NSUB, LSUB, MSUB);
    for (int k = 0; k < NSUB; ++k) {
        fprintf(fp, "%.15e", X[k]);
        for (int j = 0; j < LSUB; ++j)
            fprintf(fp, " %.15e", Y[k][j]);
        fprintf(fp, " %.15e\n", W[k]);
    }
    fclose(fp);

    // Write Fortran driver
    fp = fopen(fortran_src, "w");
    if (!fp) return false;
    fprintf(fp, R"(
      PROGRAM LSQTEST
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER (NMAX=200, MMAX=25, LMAX=10)
      DIMENSION X(NMAX), Y(NMAX,LMAX), W(NMAX), RESID(NMAX,LMAX)
      DIMENSION SUM(LMAX), A(MMAX,MMAX), B(MMAX,LMAX)
      CHARACTER*80 DFILE
      COMMON /F402/ XPOWER(100)

      OPEN(10, FILE='%s', STATUS='OLD')
      READ(10,*) NSUB, LSUB, MSUB
      DO K=1,NSUB
        READ(10,*) X(K), (Y(K,J), J=1,LSUB), W(K)
      ENDDO
      CLOSE(10)

      CALL LSQPOL(X, Y, W, RESID, NSUB, SUM, LSUB, A, B, MSUB, NMAX, MMAX)

      WRITE(*,'(A)') 'COEFFICIENTS B(I,J):'
      DO I=1,MSUB
        WRITE(*,'(3(1X,E20.12))') (B(I,J), J=1,LSUB)
      ENDDO
      WRITE(*,'(A)') 'RESIDUALS:'
      DO K=1,NSUB
        WRITE(*,'(3(1X,E20.12))') (RESID(K,J), J=1,LSUB)
      ENDDO
      WRITE(*,'(A)') 'SUMS:'
      WRITE(*,'(3(1X,E20.12))') (SUM(J), J=1,LSUB)
      STOP
      END
)", data_file);
    fclose(fp);
    return true;
}

// ─── main ──────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    printf("=== LSQPOL C++ Standalone Test ===\n\n");

    // ── Test cases ──
    // 1. Simple known polynomial: y = 1 + 2x + 3x² + 4x³, fit should recover coeffs
    // 2. Random data

    struct TestCase {
        std::string name;
        std::vector<double> X, W;
        std::vector<std::vector<double>> Y;  // [NSUB][LSUB]
        int MSUB;
    };

    std::vector<TestCase> tests;

    // ── Test 1: known cubic polynomial (no noise) ──
    {
        TestCase tc;
        tc.name = "Known cubic (y=1+2x+3x²+4x³)";
        tc.MSUB = 4;
        int N = 12;
        // coefficients for 3 columns
        double c0[3] = {1.0, 2.0, -1.0};
        double c1[3] = {2.0, -3.0, 0.5};
        double c2[3] = {3.0, 1.0, 2.0};
        double c3[3] = {4.0, -0.5, 1.5};
        for (int k = 0; k < N; ++k) {
            double x = 0.5 + k * 1.0;  // x = 0.5, 1.5, ..., 11.5
            tc.X.push_back(x);
            tc.W.push_back(1.0);
            std::vector<double> yvec(3);
            for (int j = 0; j < 3; ++j)
                yvec[j] = c0[j] + c1[j]*x + c2[j]*x*x + c3[j]*x*x*x;
            tc.Y.push_back(yvec);
        }
        tests.push_back(tc);
    }

    // ── Test 2: random data (NPSUM=40 like Ptolemy, LSUB=3, MSUB=4) ──
    {
        TestCase tc;
        tc.name = "Random (NPSUM=40, LSUB=3, MSUB=4)";
        tc.MSUB = 4;
        int N = 40;
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> dx(0.1, 30.0);
        std::uniform_real_distribution<double> dy(-15.0, 15.0);
        std::uniform_real_distribution<double> dw(0.01, 1.0);
        for (int k = 0; k < N; ++k) {
            tc.X.push_back(0.2 + k * 0.77);  // monotone, like CUBMAP output
            tc.W.push_back(dw(rng));
            std::vector<double> yvec = {dy(rng), dy(rng), dy(rng)};
            tc.Y.push_back(yvec);
        }
        tests.push_back(tc);
    }

    // ── Test 3: simulate Ptolemy's actual VMIN/VMID/VMAX data structure ──
    // VMIN negative, VMID fractional [0,1], VMAX positive
    {
        TestCase tc;
        tc.name = "Ptolemy VMIN/VMID_frac/VMAX (N=40, MSUB=4)";
        tc.MSUB = 4;
        int N = 40;
        // xi_s: CUBMAP(2, 0, 15.2, 30.4, 2.0) gives roughly sinh-spaced U points
        // Approximate: U[k] ~ SUMMAX * (k+0.5)/N for uniform, but let's use linear
        for (int k = 0; k < N; ++k) {
            double U = 0.5 + k * (30.4 - 0.5) / (N - 1);
            tc.X.push_back(U);
            double rng_vmax =  2.0 * U * (0.3 + 0.5 * k / N);
            double rng_vmin = -rng_vmax * 0.9;
            double vmid_frac = 0.5 + 0.1 * std::sin(k * 0.3);
            std::vector<double> yvec = {rng_vmin, vmid_frac, rng_vmax};
            tc.Y.push_back(yvec);
            // Weight = 1/(VMAX-VMIN)²
            double rng = rng_vmax - rng_vmin;
            tc.W.push_back(1.0 / (rng * rng));
        }
        tests.push_back(tc);
    }

    int total_tests = 0, pass = 0;

    for (auto& tc : tests) {
        printf("─── Test: %s ───\n", tc.name.c_str());
        int NSUB = (int)tc.X.size(), LSUB = 3, MSUB = tc.MSUB;

        std::vector<std::vector<double>> B_cpp, RESID_cpp;
        std::vector<double> SUM_cpp;
        lsqpol(tc.X, tc.Y, tc.W, NSUB, LSUB, MSUB, B_cpp, RESID_cpp, SUM_cpp);

        // Print coefficients
        printf("Polynomial coefficients B[i][j] (i=0=const, j=0,1,2):\n");
        for (int i = 0; i < MSUB; ++i)
            printf("  x^%d: %14.6e  %14.6e  %14.6e\n", i,
                   B_cpp[i][0], B_cpp[i][1], B_cpp[i][2]);

        // For Test 1 (known polynomial): check that coefficients match
        if (tc.name.find("Known") != std::string::npos) {
            double c[4][3] = {{1.0,2.0,-1.0},{2.0,-3.0,0.5},{3.0,1.0,2.0},{4.0,-0.5,1.5}};
            double max_err = 0;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 3; ++j)
                    max_err = std::max(max_err, std::abs(B_cpp[i][j] - c[i][j]));
            printf("Known-coeff max error: %.3e  %s\n", max_err,
                   max_err < 1e-8 ? "✅ PASS" : "❌ FAIL");
            total_tests++; if (max_err < 1e-8) pass++;
        }

        // Check residuals: sum should be ~0 when MSUB >= true degree
        if (tc.name.find("Known") != std::string::npos) {
            double max_resid = 0;
            for (int k = 0; k < NSUB; ++k)
                for (int j = 0; j < LSUB; ++j)
                    max_resid = std::max(max_resid, std::abs(RESID_cpp[k][j]));
            printf("Max residual (should be ~0): %.3e  %s\n", max_resid,
                   max_resid < 1e-8 ? "✅ PASS" : "❌ FAIL");
            total_tests++; if (max_resid < 1e-8) pass++;
        }

        // For random tests: verify residuals are consistent with B coefficients
        {
            double max_inconsistency = 0;
            for (int k = 0; k < NSUB; ++k)
                for (int j = 0; j < LSUB; ++j) {
                    double poly = poly_eval(B_cpp, j, tc.X[k]);
                    double expected_resid = poly - tc.Y[k][j];
                    max_inconsistency = std::max(max_inconsistency,
                        std::abs(RESID_cpp[k][j] - expected_resid));
                }
            printf("Residual self-consistency max err: %.3e  %s\n",
                   max_inconsistency,
                   max_inconsistency < 1e-8 ? "✅ PASS" : "❌ FAIL");
            total_tests++; if (max_inconsistency < 1e-8) pass++;
        }

        // Check that poly interpolates correctly at data points (for known case)
        {
            double max_poly_err = 0;
            for (int k = 0; k < NSUB; ++k)
                for (int j = 0; j < LSUB; ++j) {
                    double poly = poly_eval(B_cpp, j, tc.X[k]);
                    // poly(X[k]) = Y[k][j] + RESID[k][j]
                    double expect = tc.Y[k][j] + RESID_cpp[k][j];
                    max_poly_err = std::max(max_poly_err, std::abs(poly - expect));
                }
            printf("poly(X[k]) == Y[k][j]+RESID[k][j] max err: %.3e  %s\n",
                   max_poly_err, max_poly_err < 1e-10 ? "✅ PASS" : "❌ FAIL");
            total_tests++; if (max_poly_err < 1e-10) pass++;
        }

        printf("\n");
    }

    // ── Fortran comparison ──
    printf("=== Fortran comparison ===\n");
    printf("Writing Fortran driver to /tmp/lsqpol_fort_test.f ...\n");

    // Use test case 0 (known cubic, N=12, simple)
    auto& tc = tests[0];
    int NSUB = (int)tc.X.size(), LSUB = 3, MSUB = tc.MSUB;

    bool ok = write_fortran_driver(tc.X, tc.Y, tc.W, NSUB, LSUB, MSUB,
                                    "/tmp/lsqpol_fort_test.f",
                                    "/tmp/lsqpol_test_data.txt");
    if (!ok) {
        printf("Failed to write Fortran driver.\n");
    } else {
        // Compile using the Ptolemy fortlib
        int rc = system(
            "gfortran -O2 -o /tmp/lsqpol_fort "
            "/tmp/lsqpol_fort_test.f "
            "/home/node/working/ptolemy_2019/build/fortlib.mor "
            "2>/tmp/lsqpol_fort_compile.log"
        );
        if (rc != 0) {
            printf("Fortran compile failed. Log:\n");
            system("cat /tmp/lsqpol_fort_compile.log");
        } else {
            // Run Fortran
            system("/tmp/lsqpol_fort > /tmp/lsqpol_fort_out.txt 2>&1");
            printf("Fortran output:\n");
            system("cat /tmp/lsqpol_fort_out.txt");

            // Re-run C++ for test case 0
            std::vector<std::vector<double>> B_cpp, RESID_cpp;
            std::vector<double> SUM_cpp;
            lsqpol(tc.X, tc.Y, tc.W, NSUB, LSUB, MSUB, B_cpp, RESID_cpp, SUM_cpp);

            printf("\nC++ output:\n");
            printf("COEFFICIENTS B(I,J):\n");
            for (int i = 0; i < MSUB; ++i)
                printf(" %20.12e %20.12e %20.12e\n",
                       B_cpp[i][0], B_cpp[i][1], B_cpp[i][2]);
            printf("RESIDUALS:\n");
            for (int k = 0; k < NSUB; ++k)
                printf(" %20.12e %20.12e %20.12e\n",
                       RESID_cpp[k][0], RESID_cpp[k][1], RESID_cpp[k][2]);

            // Now parse Fortran output and compare
            // Parse Fortran B matrix
            FILE* fp = fopen("/tmp/lsqpol_fort_out.txt", "r");
            if (fp) {
                char line[256];
                // Skip "COEFFICIENTS B(I,J):" header
                while (fgets(line, sizeof(line), fp))
                    if (strstr(line, "COEFFICIENTS")) break;

                double max_B_err = 0, max_R_err = 0;
                for (int i = 0; i < MSUB; ++i) {
                    double b0, b1, b2;
                    if (fscanf(fp, "%lf %lf %lf", &b0, &b1, &b2) == 3) {
                        max_B_err = std::max(max_B_err, std::abs(b0 - B_cpp[i][0]));
                        max_B_err = std::max(max_B_err, std::abs(b1 - B_cpp[i][1]));
                        max_B_err = std::max(max_B_err, std::abs(b2 - B_cpp[i][2]));
                    }
                }
                // Skip "RESIDUALS:" header
                fgets(line, sizeof(line), fp);  // newline after last coeff
                while (fgets(line, sizeof(line), fp))
                    if (strstr(line, "RESIDUALS")) break;
                for (int k = 0; k < NSUB; ++k) {
                    double r0, r1, r2;
                    if (fscanf(fp, "%lf %lf %lf", &r0, &r1, &r2) == 3) {
                        max_R_err = std::max(max_R_err, std::abs(r0 - RESID_cpp[k][0]));
                        max_R_err = std::max(max_R_err, std::abs(r1 - RESID_cpp[k][1]));
                        max_R_err = std::max(max_R_err, std::abs(r2 - RESID_cpp[k][2]));
                    }
                }
                fclose(fp);

                printf("\n--- Fortran vs C++ comparison ---\n");
                printf("Max coeff error:   %.3e  %s\n", max_B_err,
                       max_B_err < 1e-8 ? "✅ PASS" : "❌ FAIL");
                printf("Max residual error: %.3e  %s\n", max_R_err,
                       max_R_err < 1e-8 ? "✅ PASS" : "❌ FAIL");
                total_tests += 2;
                if (max_B_err < 1e-8) pass++;
                if (max_R_err < 1e-8) pass++;
            }
        }
    }

    printf("\n=== SUMMARY: %d / %d tests passed ===\n", pass, total_tests);
    return (pass == total_tests) ? 0 : 1;
}
