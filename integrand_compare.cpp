// integrand_compare.cpp
// Compare C++ H-function integrand vs Fortran integrand_test output.
//
// H(ra,rb) = sum_k DPHI_k * phi_T(rx_k) * IVPHI_P(rp_k) * A12(PHIT_k, PHI_k)
//
// Uses same analytic tables, same CUBMAP phi grid, same kinematics as integrand_test.f
//
// Build: g++ -O2 -std=c++17 integrand_compare.cpp -o integrand_compare
// Run:   ../fortran_testing/integrand_test > /tmp/fort_integrand.txt
//        ./integrand_compare < /tmp/fort_integrand.txt

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// ---- Gauss-Legendre (same as bsprod/cubmap tests) ----
static void GaussLegendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n); w.resize(n);
    if (n == 1) { x[0]=0; w[0]=2; return; }
    int m=(n+1)/2;
    double e1=n*(n+1), e2=M_PI/(4*n+2), e3=1.0-(1.0-1.0/n)/(8.0*n*n);
    for (int i=1; i<=m; i++) {
        double t=(4*i-1)*e2, x0=e3*std::cos(t);
        double pkm1=1,pk=x0,fk=1;
        for (int k=2;k<=n;k++){fk++;double t1=x0*pk;double pkp1=t1-pkm1-(t1-pkm1)/fk+t1;pkm1=pk;pk=pkp1;}
        double den=1-x0*x0,d1=n*(pkm1-x0*pk),dpn=d1/den;
        double d2pn=(2*x0*dpn-e1*pk)/den,d3pn=(4*x0*d2pn+(2-e1)*dpn)/den,d4pn=(6*x0*d3pn+(6-e1)*d2pn)/den;
        double u=pk/dpn,v=d2pn/dpn,h1=1+0.5*u*(v+u*(v*v-u*d3pn/(3*dpn))),h=-u*h1;
        double p=pk+h*(dpn+0.5*h*(d2pn+h/3*(d3pn+0.25*h*d4pn)));
        double dp=dpn+h*(d2pn+0.5*h*(d3pn+h*d4pn/3)); h=h-p/dp;
        x[n-i]=x0+h; x[i-1]=-(x0+h);
        double fx=d1-h*e1*(pk+0.5*h*(dpn+h/3*(d2pn+0.25*h*(d3pn+0.2*h*d4pn))));
        w[n-i]=2*(1-x[i-1]*x[i-1])/(fx*fx); w[i-1]=w[n-i];
    }
    if ((m+m)>n) x[m-1]=0;
}

// ---- CubMap ----
static void CubMap(int maptyp, double xlo, double xmid, double xhi, double gamma,
                   std::vector<double>& args, std::vector<double>& wts) {
    int npts=(int)args.size();
    double tau=(gamma>1e-6)?std::log(gamma+std::sqrt(gamma*gamma+1.0)):gamma*(1.0-gamma*gamma/6.0);
    double xlen=xhi-xlo, xadd=xlo+xhi;
    if (maptyp==0){for(int i=0;i<npts;i++){args[i]=xlo+0.5*xlen*(args[i]+1.0);wts[i]*=0.5*xlen;}}
    else if (maptyp==1){xmid=std::max(xmid,xlo+xlen/7.0);xmid=std::min(xmid,0.5*xadd);double A=0.5*xadd-xmid,B=0.5*xlen,C=0.5*xadd;for(int i=0;i<npts;i++){double tu=tau*args[i],xi=std::sinh(tu)/gamma;args[i]=A*(xi*xi-1.0)*(xi+1.0)+B*xi+C;wts[i]*=(tau/gamma)*std::cosh(tu)*((3.0*xi-1.0)*(xi+1.0)*A+B);}}
    else if (maptyp==2){double A=-xmid*xlen,B=xlen,C=xmid*xadd-2.0*xlo*xhi,D=xadd-2.0*xmid;for(int i=0;i<npts;i++){double tu=tau*args[i],sh=std::sinh(tu),denom=B-(D/gamma)*sh;args[i]=(-A+(C/gamma)*sh)/denom;wts[i]*=(tau/gamma)*std::cosh(tu)*((B*C-A*D)/(denom*denom));}}
    else if (maptyp==3){for(int i=0;i<npts;i++){double tu=tau*args[i];args[i]=xlo+0.5*xlen*(std::sinh(tu)/gamma+1.0);wts[i]*=0.5*xlen*(tau/gamma)*std::cosh(tu);}}
}

// ---- AITLAG ----
static double Aitlag(double x, double stpinv, const std::vector<double>& table, int nait) {
    int ntable=(int)table.size();
    if (x<0) return table[0];
    double xbyh=x*stpinv; int itable=(int)xbyh;
    if (itable>=ntable) return table[ntable-1];
    static const double inv[17]={0,1.0,0.5,1.0/3,0.25,0.2,1.0/6,1.0/7,0.125,1.0/9,0.1,1.0/11,1.0/12,1.0/13,1.0/14,1.0/15,1.0/16};
    int nait2=nait/2,nstart=itable-nait2; if(nstart<0)nstart=0;
    int nend=nstart+nait; if(nend>=ntable)nend=ntable-1; nstart=nend-nait;
    std::vector<double> fs(nait+2),dels(nait+2);
    fs[1]=table[nstart]; double del1=xbyh-nstart; dels[1]=del1; double f=0;
    for(int i=1;i<=nait;i++){f=table[nstart+i];del1-=1.0;for(int j=1;j<=i;j++){f=(f*dels[j]-fs[j]*del1)*inv[i+1-j];}dels[i+1]=del1;fs[i+1]=f;}
    return f;
}

// ---- A12 kernel for Li=0,Lo=2,Lx=2,lT=2,lP=0 (same as Fortran EVAL_A12) ----
static double EvalA12(double PHIT, double /*PHI*/) {
    // From validated a12.cpp: terms (MT=0,MU=0,C0) and (MT=2,MU=0,C2)
    const double C0 =  0.2820947917738781;
    const double C2 = -0.2308356256733685;
    return C0 + C2 * std::cos(2.0*PHIT);
}

int main() {
    // ---- Build same tables as Fortran ----
    const int NTAB_T=641, NTAB_P=641;
    const double H_T=0.05, H_P=0.05;
    std::vector<double> phi_T(NTAB_T), ivphi_P(NTAB_P);
    for(int i=0;i<NTAB_T;i++){double r=i*H_T;phi_T[i]=(i==0)?0.0:r*std::exp(-0.433*r);}
    for(int i=0;i<NTAB_P;i++){double r=i*H_P;ivphi_P[i]=-0.8*std::exp(-0.3*r*r);}
    const double bndmxp=(NTAB_P-1)*H_P, bndmxt=(NTAB_T-1)*H_T;
    const int NAIT=4, NPPHI=10;

    // Kinematics
    const double BRATMS1=1.0, BRATMS2=1.0/16.0;
    double temp=1.0/(BRATMS1+BRATMS2*(1.0+BRATMS1));
    double S1=(1.0+BRATMS1)*(1.0+BRATMS2)*temp;
    double T1=-(1.0+BRATMS2)*temp;
    double S2=(1.0+BRATMS1)*temp;
    double T2=-S1;

    // PHI grid: CUBMAP(MAPPHI=2, XLO=0, XMID=0.20, XHI=1, GAMPHI=1e-6)
    std::vector<double> phi_pts(NPPHI), phi_wts(NPPHI);
    GaussLegendre(NPPHI, phi_pts, phi_wts);
    double xmid_phi = 0.20;
    CubMap(2, 0.0, xmid_phi, 1.0, 1e-6, phi_pts, phi_wts);

    // ---- Read Fortran output ----
    struct FRow { double ra, rb, phi0, H; };
    std::vector<FRow> fort_rows;
    std::string line;
    bool header_done = false;
    while (std::getline(std::cin, line)) {
        if (line.find("RA") != std::string::npos) { header_done=true; continue; }
        if (!header_done) continue;
        if (line.find("phi") != std::string::npos) continue;
        std::istringstream ss(line);
        FRow r;
        if (ss >> r.ra >> r.rb >> r.phi0 >> r.H) fort_rows.push_back(r);
    }

    // ---- Compare ----
    double max_rel = 0;
    int nbad = 0;
    std::cout << std::setw(5)<<"#"<<std::setw(6)<<"RA"<<std::setw(6)<<"RB"
              <<std::setw(20)<<"F_H"<<std::setw(20)<<"C_H"<<std::setw(12)<<"relΔ"<<"\n";

    int idx=0;
    for (int ira=1; ira<=5; ira++) {
        double ra=ira*1.0;
        for (int irb=1; irb<=5; irb++) {
            double rb=irb*1.0;
            const double PHI0 = M_PI/2.0;

            // C++ H computation (same as ineldc.cpp phi loop)
            double H_cpp = 0.0;
            for (int k=0; k<NPPHI; k++) {
                double phi_frac = phi_pts[k];
                double PHI      = PHI0 * phi_frac;
                double DPHI     = PHI0 * phi_wts[k] * std::sin(PHI);
                double X        = std::cos(PHI);

                double rx2 = S1*S1*ra*ra + T1*T1*rb*rb + 2.0*S1*T1*ra*rb*X;
                double rp2 = S2*S2*ra*ra + T2*T2*rb*rb + 2.0*S2*T2*ra*rb*X;
                if (rx2<0) rx2=0; if (rp2<0) rp2=0;
                double rx=std::sqrt(rx2), rp=std::sqrt(rp2);

                if (rx>bndmxt || rp>bndmxp) continue;

                double FT = Aitlag(rx, 1.0/H_T, phi_T,  NAIT);
                double FP = Aitlag(rp, 1.0/H_P, ivphi_P, NAIT);
                double PVPDX = DPHI * FT * FP;

                double cos_phiT = (rx2>1e-30) ? (T1*rb + S1*ra*X)/rx : 1.0;
                cos_phiT = std::max(-1.0, std::min(1.0, cos_phiT));
                double PHIT = std::acos(cos_phiT);

                double A12k = EvalA12(PHIT, PHI);
                H_cpp += PVPDX * A12k;
            }

            if (idx < (int)fort_rows.size()) {
                double fH = fort_rows[idx].H;
                double ref = std::abs(fH);
                double rel = (ref > 1e-30) ? std::abs(H_cpp - fH)/ref : std::abs(H_cpp - fH);
                max_rel = std::max(max_rel, rel);
                if (rel > 1e-9) { nbad++; }

                std::cout << std::setw(5)<<idx+1<<std::setw(6)<<(int)ra<<std::setw(6)<<(int)rb
                          <<std::setw(20)<<std::setprecision(10)<<std::scientific<<fH
                          <<std::setw(20)<<H_cpp
                          <<std::setw(12)<<std::setprecision(3)<<rel
                          <<(rel>1e-9?" <<<":"")<<"\n";
                idx++;
            }
        }
    }

    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "N compared: " << idx << "\n";
    std::cout << "Max rel error H: " << std::scientific << max_rel << "\n";
    std::cout << "Points |Δ|>1e-9: " << nbad << "\n";
    if (max_rel < 1e-9)
        std::cout << "RESULT: PASS ✓ — C++ H-integrand matches Fortran to FP precision\n";
    else
        std::cout << "RESULT: FAIL ✗ — discrepancies found!\n";

    return 0;
}
