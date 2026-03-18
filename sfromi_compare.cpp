// sfromi_compare.cpp
// Compare C++ SFROMI core (pre-9J: TEMP, ITEST, SR/SI) vs Fortran output.
//
// Tests: FACTOR*ATERM/sqrt(2*Li+1), ITEST = Li+Lo+2*Lx+1,
//        sign flip if MOD(ITEST,4)>=2, rotation by i if MOD(ITEST,2)==1.
//
// Build: g++ -O2 -std=c++17 sfromi_compare.cpp -o sfromi_compare
// Run:   ../fortran_testing/sfromi_test > /tmp/fort_sfromi.txt
//        ./sfromi_compare < /tmp/fort_sfromi.txt

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

// ---- 32-bit LCG matching Fortran INTEGER*4 signed arithmetic ----
static int32_t lcg_next(int32_t seed) {
    // Fortran: ISEED4 = 1664525*ISEED4 + 1013904223  (INTEGER*4, wraps at 2^32)
    return (int32_t)((int64_t)1664525 * seed + 1013904223);
}

// ---- SFROMI core (C++) ----
struct SFromiResult {
    double SR, SI;
    int ITEST;
};

static SFromiResult sfromi_core(int Li, int Lo, int Lx,
                                 double XI_R, double XI_I,
                                 double FACTOR, double ATERM) {
    double TEMP = FACTOR * ATERM / std::sqrt(2.0 * Li + 1.0);
    int ITEST = Li + Lo + 2*Lx + 1;
    if ((ITEST % 4 + 4) % 4 >= 2) TEMP = -TEMP;   // MOD(ITEST,4)>=2
    double TEMPR = TEMP * XI_R;
    double TEMPI = TEMP * XI_I;
    double SR = TEMPR, SI = TEMPI;
    if (ITEST % 2 != 0) {   // MOD(ITEST,2) != 0 → rotate by i
        SR = -TEMPI;
        SI =  TEMPR;
    }
    return {SR, SI, ITEST};
}

int main() {
    const double FACTOR = 0.11100;
    const double ATERM  = -1.681284;

    // ---- Generate same (Li,Lo,Lx,XI_R,XI_I) sequence as Fortran ----
    struct TestRow { int Li, Lo, Lx, ITEST; double xi_r, xi_i, sr, si; };
    std::vector<TestRow> cpp_rows;

    int32_t iseed = 12345;
    for (int Li = 0; Li <= 8; Li++) {
        for (int Lo = 0; Lo <= 8; Lo++) {
            int Lx = 2;
            iseed = lcg_next(iseed);
            double xi_r = (double)iseed * 4.656612873e-10 * 0.2;
            iseed = lcg_next(iseed);
            double xi_i = (double)iseed * 4.656612873e-10 * 0.2;

            auto res = sfromi_core(Li, Lo, Lx, xi_r, xi_i, FACTOR, ATERM);
            cpp_rows.push_back({Li, Lo, Lx, res.ITEST, xi_r, xi_i, res.SR, res.SI});
        }
    }

    // ---- Read Fortran output ----
    struct FRow { int Li, Lo, Lx, ITEST; double xi_r, xi_i, sr, si; };
    std::vector<FRow> fort_rows;
    std::string line;
    bool header_done = false;
    while (std::getline(std::cin, line)) {
        if (line.find("Li") != std::string::npos) { header_done = true; continue; }
        if (!header_done) continue;
        std::istringstream ss(line);
        FRow r;
        if (ss >> r.Li >> r.Lo >> r.Lx >> r.ITEST >> r.xi_r >> r.xi_i >> r.sr >> r.si)
            fort_rows.push_back(r);
    }

    // ---- Compare ----
    double max_rel_sr = 0, max_rel_si = 0, max_rel_xi = 0;
    int nbad = 0;

    std::cout << std::setw(3)<<"Li"<<std::setw(3)<<"Lo"<<std::setw(3)<<"Lx"
              <<std::setw(5)<<"ITEST"
              <<std::setw(16)<<"F_SR"<<std::setw(16)<<"C_SR"<<std::setw(11)<<"relΔ_SR"
              <<std::setw(16)<<"F_SI"<<std::setw(16)<<"C_SI"<<std::setw(11)<<"relΔ_SI"
              <<std::setw(11)<<"relΔ_XI_R"<<"\n";

    int n = std::min((int)fort_rows.size(), (int)cpp_rows.size());
    for (int i = 0; i < n; i++) {
        auto& f = fort_rows[i];
        auto& c = cpp_rows[i];

        // Check XI_R matches (LCG must be identical)
        double rel_xi = (std::abs(f.xi_r) > 1e-30)
            ? std::abs(c.xi_r - f.xi_r)/std::abs(f.xi_r) : std::abs(c.xi_r - f.xi_r);
        max_rel_xi = std::max(max_rel_xi, rel_xi);

        double ref_sr = std::abs(f.sr), ref_si = std::abs(f.si);
        double rel_sr = (ref_sr > 1e-30) ? std::abs(c.sr - f.sr)/ref_sr : std::abs(c.sr - f.sr);
        double rel_si = (ref_si > 1e-30) ? std::abs(c.si - f.si)/ref_si : std::abs(c.si - f.si);
        max_rel_sr = std::max(max_rel_sr, rel_sr);
        max_rel_si = std::max(max_rel_si, rel_si);
        if (rel_sr > 1e-10 || rel_si > 1e-10) nbad++;

        // Flag ITEST mismatch
        std::string flag = (f.ITEST != c.ITEST) ? " <<<ITEST!" : "";
        flag += (rel_sr > 1e-10 || rel_si > 1e-10) ? " <<<VAL!" : "";

        std::cout << std::setw(3) << f.Li << std::setw(3) << f.Lo << std::setw(3) << f.Lx
                  << std::setw(5) << f.ITEST
                  << std::setw(16) << std::setprecision(8) << std::scientific << f.sr
                  << std::setw(16) << c.sr
                  << std::setw(11) << std::setprecision(2) << rel_sr
                  << std::setw(16) << std::setprecision(8) << f.si
                  << std::setw(16) << c.si
                  << std::setw(11) << std::setprecision(2) << rel_si
                  << std::setw(11) << std::setprecision(2) << rel_xi
                  << flag << "\n";
    }

    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "N compared: " << n << "\n";
    std::cout << "Max rel error SR:    " << std::scientific << max_rel_sr << "\n";
    std::cout << "Max rel error SI:    " << max_rel_si << "\n";
    std::cout << "Max rel error XI_R:  " << max_rel_xi << " (LCG seed check)\n";
    std::cout << "Points with |Δ|>1e-10: " << nbad << "\n";
    // Threshold: 1e-9 (FP rounding in sqrt+multiply chain, same as BSPROD/CUBMAP)
    if (max_rel_sr < 1e-9 && max_rel_si < 1e-9 && max_rel_xi < 1e-9)
        std::cout << "RESULT: PASS ✓ — C++ SFROMI core matches Fortran to FP precision (thr 1e-9)\n";
    else
        std::cout << "RESULT: FAIL ✗ — discrepancies exceed 1e-9!\n";

    return 0;
}
