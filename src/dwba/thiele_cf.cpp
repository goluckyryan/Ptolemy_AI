// thiele_cf.cpp — Thiele continued fraction interpolation/extrapolation
// Ported from Ptolemy Fortran CCNFRC/CCONTF (fortlib.f)
// Fits rational continued fraction to complex-valued tabulated data
// Then evaluates at new points for extrapolation

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cstdio>

// CCNFRC: Setup the continued fraction coefficients
// Input: xs[i], ys[i] = (x, y) pairs where y = f(x)
// Output: xs and ys are rearranged; ys becomes the CF coefficients
// Returns nmax = number of coefficients used (<= numpts-1)
int ThieleCF_Setup(int numpts, 
                   std::vector<std::complex<double>>& xs, 
                   std::vector<std::complex<double>>& ys) {
    const double TINY = 1.0e-14;
    int nmax = numpts - 1;
    if (nmax <= 0) return 0;
    
    // Find the largest |y| and swap to position 0
    double comp = -1.0;
    int k = 0;
    for (int i = 0; i < numpts; i++) {
        double temp = std::norm(ys[i]);  // |y|^2
        if (temp > comp) { comp = temp; k = i; }
    }
    
    // Swap largest to front
    if (k != 0) {
        std::swap(ys[k], ys[0]);
        std::swap(xs[k], xs[0]);
    }
    
    // Calculate candidates for second coefficient
    auto yj = ys[0];
    comp = -1.0;
    k = 1;
    for (int i = 1; i < numpts; i++) {
        auto d = 1.0 - yj / ys[i];
        ys[i] = d;
        double temp = std::norm(d);
        if (temp > comp) { comp = temp; k = i; }
    }
    
    // Loop through rest of coefficients
    for (int j = 1; j < numpts; j++) {
        // Swap largest candidate to position j
        auto xj = xs[k];
        yj = ys[k];
        if (k != j) {
            ys[k] = ys[j];
            xs[k] = xs[j];
            xs[j] = xj;
        }
        // Finish calculating this coefficient
        yj = yj / (xs[j-1] - xj);
        ys[j] = yj;
        
        if (j == numpts - 1) break;
        
        // Quit if coefficient is tiny
        if (comp < TINY) {
            nmax = j;
            fprintf(stderr, "**** CF used only %d of %d points\n", nmax+1, numpts);
            break;
        }
        
        // Calculate candidates for next coefficient
        comp = -1.0;
        xj = xs[j-1];
        for (int i = j + 1; i < numpts; i++) {
            auto d = 1.0 + yj * (xs[i] - xj) / ys[i];
            ys[i] = d;
            double temp = std::norm(d);
            if (temp > comp) { comp = temp; k = i; }
        }
    }
    
    return nmax;
}

// CCONTF: Evaluate the continued fraction at point x
// xs, ys must be the output from ThieleCF_Setup
// nmax = return value from ThieleCF_Setup
std::complex<double> ThieleCF_Eval(int nmax,
                                    const std::vector<std::complex<double>>& xs,
                                    const std::vector<std::complex<double>>& ys,
                                    std::complex<double> x) {
    std::complex<double> y(0.0, 0.0);
    if (nmax < 1) return ys[0] / (1.0 + y);
    
    for (int j = 0; j < nmax; j++) {
        int k = nmax - j;  // k goes nmax, nmax-1, ..., 1
        y = ys[k] * (x - xs[k-1]) / (1.0 + y);
    }
    y = ys[0] / (1.0 + y);
    return y;
}
