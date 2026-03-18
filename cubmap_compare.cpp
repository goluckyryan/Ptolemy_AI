// cubmap_compare.cpp
// Compare C++ CubMap + GaussLegendre vs Fortran CUBMAP + GAUSSL output.
// Reads fortran output from stdin (or file), prints comparison table.
//
// Build: g++ -O2 -std=c++17 cubmap_compare.cpp -o cubmap_compare
// Run:   ./fortran_testing/cubmap_test > /tmp/fort_cubmap.txt
//        ./cubmap_compare < /tmp/fort_cubmap.txt

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ---- Gauss-Legendre points on [-1,1] (same algorithm as Ptolemy GAUSSL) ----
static void GaussLegendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n); w.resize(n);
    if (n == 1) { x[0]=0; w[0]=2; return; }
    int m = (n+1)/2;
    double e1 = n*(n+1);
    double e2 = M_PI/(4*n+2);
    double e3 = 1.0-(1.0-1.0/n)/(8.0*n*n);
    for (int i=1; i<=m; i++) {
        double t = (4*i-1)*e2;
        double x0 = e3*std::cos(t);
        double pkm1=1, pk=x0, fk=1;
        for (int k=2; k<=n; k++) {
            fk++;
            double t1=x0*pk;
            double pkp1=t1-pkm1-(t1-pkm1)/fk+t1;
            pkm1=pk; pk=pkp1;
        }
        double den=1-x0*x0;
        double d1=n*(pkm1-x0*pk);
        double dpn=d1/den;
        double d2pn=(2*x0*dpn-e1*pk)/den;
        double d3pn=(4*x0*d2pn+(2-e1)*dpn)/den;
        double d4pn=(6*x0*d3pn+(6-e1)*d2pn)/den;
        double u=pk/dpn, v=d2pn/dpn;
        double h1=1+0.5*u*(v+u*(v*v-u*d3pn/(3*dpn)));
        double h=-u*h1;
        double p=pk+h*(dpn+0.5*h*(d2pn+h/3*(d3pn+0.25*h*d4pn)));
        double dp=dpn+h*(d2pn+0.5*h*(d3pn+h*d4pn/3));
        h=h-p/dp;
        x[n-i]=x0+h; x[i-1]=-(x0+h);
        double fx=d1-h*e1*(pk+0.5*h*(dpn+h/3*(d2pn+0.25*h*(d3pn+0.2*h*d4pn))));
        w[n-i]=2*(1-x[i-1]*x[i-1])/(fx*fx);
        w[i-1]=w[n-i];
    }
    if ((m+m)>n) x[m-1]=0;
}

// ---- CubMap (matches Ptolemy CUBMAP, modifies args/wts in place) ----
static void CubMap(int maptyp, double xlo, double xmid, double xhi, double gamma,
                   std::vector<double>& args, std::vector<double>& wts) {
    int npts = (int)args.size();
    double tau = (gamma > 1e-6) ? std::log(gamma + std::sqrt(gamma*gamma + 1.0))
                                : gamma*(1.0 - gamma*gamma/6.0);
    double xlen = xhi - xlo;
    double xadd = xlo + xhi;

    if (maptyp == 0) {
        for (int i=0; i<npts; i++) {
            args[i] = xlo + 0.5*xlen*(args[i]+1.0);
            wts[i] *= 0.5*xlen;
        }
    } else if (maptyp == 1) {
        xmid = std::max(xmid, xlo + xlen/7.0);
        xmid = std::min(xmid, 0.5*xadd);
        double A = 0.5*xadd - xmid;
        double B = 0.5*xlen;
        double C = 0.5*xadd;
        for (int i=0; i<npts; i++) {
            double tu = tau*args[i];
            double xi = std::sinh(tu)/gamma;
            args[i] = A*(xi*xi-1.0)*(xi+1.0) + B*xi + C;
            wts[i] *= (tau/gamma)*std::cosh(tu)*((3.0*xi-1.0)*(xi+1.0)*A + B);
        }
    } else if (maptyp == 2) {
        double A = -xmid*xlen;
        double B = xlen;
        double C = xmid*xadd - 2.0*xlo*xhi;
        double D = xadd - 2.0*xmid;
        for (int i=0; i<npts; i++) {
            double tu = tau*args[i];
            double sh = std::sinh(tu);
            double denom = B - (D/gamma)*sh;
            args[i] = (-A + (C/gamma)*sh)/denom;
            wts[i] *= (tau/gamma)*std::cosh(tu)*((B*C - A*D)/(denom*denom));
        }
    } else if (maptyp == 3) {
        for (int i=0; i<npts; i++) {
            double tu = tau*args[i];
            args[i] = xlo + 0.5*xlen*(std::sinh(tu)/gamma + 1.0);
            wts[i] *= 0.5*xlen*(tau/gamma)*std::cosh(tu);
        }
    }
}

struct TestCase {
    int maptyp, npts;
    double xlo, xmid, xhi, gamma;
    std::string label;
};

int main() {
    // Same 7 test cases as Fortran
    std::vector<TestCase> tests = {
        {0, 10, 1.0, 3.0, 5.0, 2.0, "T1"},
        {1, 10, 0.0, 3.0, 8.0, 12.0, "T2"},
        {2, 10, 0.0, 4.0, 15.0, 2.0, "T3"},
        {3, 10, 0.0, 2.0, 10.0, 5.0, "T4"},
        {2, 40, 1.0, 4.0, 15.0, 2.0, "T5"},
        {1, 40, -3.0, 0.0, 3.0, 12.0, "T6"},
        {2, 20, 0.0, 0.20, 1.0, 1e-6, "T7"},
    };

    // Read Fortran output
    std::vector<std::vector<std::pair<double,double>>> fort_data(7);
    std::string line;
    int test_idx = -1;
    while (std::getline(std::cin, line)) {
        if (line.size() >= 2 && line[0]=='T' && line[1]>='1' && line[1]<='7') {
            test_idx = line[1]-'1';
            continue;
        }
        if (test_idx < 0) continue;
        std::istringstream ss(line);
        int idx; double a, wt;
        if (ss >> idx >> a >> wt) {
            fort_data[test_idx].emplace_back(a, wt);
        }
    }

    // Compare
    double max_rel_arg = 0, max_rel_wt = 0;
    int total_bad_arg = 0, total_bad_wt = 0;

    for (int t=0; t<(int)tests.size(); t++) {
        auto& tc = tests[t];
        std::vector<double> args(tc.npts), wts(tc.npts);
        GaussLegendre(tc.npts, args, wts);
        CubMap(tc.maptyp, tc.xlo, tc.xmid, tc.xhi, tc.gamma, args, wts);

        std::cout << "\n=== " << tc.label
                  << " MAPTYP=" << tc.maptyp
                  << " N=" << tc.npts
                  << " XLO=" << tc.xlo << " XMID=" << tc.xmid
                  << " XHI=" << tc.xhi << " GAM=" << tc.gamma << " ===\n";
        std::cout << std::setw(4) << "I"
                  << std::setw(22) << "F_ARG"
                  << std::setw(22) << "C_ARG"
                  << std::setw(12) << "relΔARG"
                  << std::setw(22) << "F_WT"
                  << std::setw(22) << "C_WT"
                  << std::setw(12) << "relΔWT" << "\n";

        auto& fd = fort_data[t];
        for (int i=0; i<tc.npts && i<(int)fd.size(); i++) {
            double fa = fd[i].first, fw = fd[i].second;
            double ca = args[i], cw = wts[i];
            double rel_a = (fa != 0) ? std::abs(ca-fa)/std::abs(fa) : std::abs(ca);
            double rel_w = (fw != 0) ? std::abs(cw-fw)/std::abs(fw) : std::abs(cw);
            max_rel_arg = std::max(max_rel_arg, rel_a);
            max_rel_wt  = std::max(max_rel_wt,  rel_w);
            if (rel_a > 1e-10) total_bad_arg++;
            if (rel_w > 1e-10) total_bad_wt++;
            std::cout << std::setw(4) << i+1
                      << std::setw(22) << std::setprecision(14) << fa
                      << std::setw(22) << ca
                      << std::setw(12) << std::setprecision(3) << std::scientific << rel_a
                      << std::setw(22) << std::setprecision(14) << std::fixed << fw
                      << std::setw(22) << cw
                      << std::setw(12) << std::setprecision(3) << std::scientific << rel_w
                      << "\n";
        }
    }

    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "Max rel error ARG: " << std::scientific << max_rel_arg << "\n";
    std::cout << "Max rel error WT:  " << std::scientific << max_rel_wt  << "\n";
    std::cout << "Points with |ΔARG|>1e-10: " << total_bad_arg << "\n";
    std::cout << "Points with |ΔWT|>1e-10:  " << total_bad_wt  << "\n";
    if (max_rel_arg < 1e-10 && max_rel_wt < 1e-10)
        std::cout << "RESULT: PASS ✓ — C++ CubMap matches Fortran to floating-point precision\n";
    else
        std::cout << "RESULT: FAIL ✗ — discrepancies found!\n";

    return 0;
}
