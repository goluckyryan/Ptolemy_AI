// spline_compare.cpp
// Validate C++ SPLNCB+INTRPC port against Fortran output.
// Reads Fortran spline_test output from stdin, recomputes in C++, compares.
//
// Build: g++ -O2 -std=c++17 spline_compare.cpp -o spline_compare
// Run:   ../fortran_testing/spline_test > /tmp/fort_spline.txt
//        ./spline_compare < /tmp/fort_spline.txt

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>

// =============================================================
// Ptolemy SPLNCB — natural cubic spline coefficients
// Exact port of fortlib.mor lines 4946-5120.
// Inputs:  N, X[N], Y[N]
// Outputs: B[N], C[N], D[N]  (same indexing as Fortran, 1-based stored as 0-based here)
// The cubic on interval [X[i], X[i+1]] is:
//   Y(x) = Y[i] + (x-X[i])*(B[i] + (x-X[i])*(C[i] + (x-X[i])*D[i]))
// =============================================================
// Ptolemy SPLNCB — faithful C++ port of fortlib.mor lines 4946-5120.
// Key indexing: Fortran arrays are 1-based. IBASE is the Fortran 0-offset
// (IBASE_fortran = I_1based - 1 = I_0based). J is 1-based in Fortran.
// NUSE = J_1based - IBASE_fortran = (J_0based+1) - IBASE_0based.
void Splncb(int N, const std::vector<double>& X, const std::vector<double>& Y,
            std::vector<double>& B, std::vector<double>& C, std::vector<double>& D) {
    B.assign(N, 0.0); C.assign(N, 0.0); D.assign(N, 0.0);
    if (N <= 1) return;

    int NM1 = N - 1;

    // ---- Forward scan: find IBASE (0-based offset, = Fortran I_1based - 1) ----
    // Fortran: DO 74 I=1,NM1; IF 1e15*|Y(I+1)-Y(I)|>|Y(I)|: GO TO 80
    // IBASE = I - 1  (Fortran assigns after loop, using last I=NM1 as default)
    int IBASE = NM1 - 1;  // default: Fortran falls through to I=NM1, IBASE=NM1-1
    for (int ii = 0; ii < NM1; ii++) {
        // B,C,D already zeroed by assign above
        if (1.0e15 * std::abs(Y[ii+1] - Y[ii]) > std::abs(Y[ii])) {
            IBASE = ii;  // Fortran: IBASE = I-1 = (ii+1)-1 = ii
            goto found_ibase;
        }
    }
    IBASE = NM1 - 1;
found_ibase:;

    // ---- Backward scan: find J (1-based in Fortran) ----
    // Fortran: DO 84 I1=I,NM1: J(1based)=N+I_1based-I1_1based
    // In 0-based: I1_0based from IBASE to NM1-1
    //             J_1based = N + (IBASE+1) - (I1_0based+1) = N + IBASE - I1_0based
    //             J_0based = J_1based - 1 = N + IBASE - I1_0based - 1
    // NUSE = J_1based - IBASE_fortran = J_1based - IBASE
    //      = (N + IBASE - I1_0based) - IBASE = N - I1_0based
    int J_0based = N - 1;   // default: falls through → J_1based=IBASE+1, but Fortran default is different
    // Actually Fortran: if loop completes without firing, J is left at last value = IBASE+1 (1based)
    // Let's just faithfully step:
    {
        int I_1based = IBASE + 1;  // the Fortran I value
        for (int I1 = I_1based; I1 <= NM1; I1++) {
            int Jf = N + I_1based - I1;  // Fortran J (1-based)
            int J0 = Jf - 1;             // 0-based
            // B[J0]=C[J0]=D[J0]=0;  already zeroed
            if (1.0e15 * std::abs(Y[J0-1] - Y[J0]) > std::abs(Y[J0])) {
                J_0based = J0;
                goto found_j;
            }
        }
        // Loop completed: J_1based = I_1based+1 (last iteration I1=NM1: Jf=N+I_1based-NM1)
        J_0based = N + I_1based - NM1 - 1;  // 0-based
    }
found_j:;

    // NUSE = J_1based - IBASE_fortran = (J_0based+1) - IBASE
    int NUSE = (J_0based + 1) - IBASE;
    if (NUSE <= 1) return;
    NM1 = NUSE - 1;

    // ---- Forward sweep (builds intermediate B,C,D for tridiagonal) ----
    double H = X[IBASE+1] - X[IBASE];
    double F = (Y[IBASE+1] - Y[IBASE]) / H;
    if (NUSE == 2) {
        D[IBASE] = 0.0;
        B[IBASE] = F;
        return;
    }

    // Fortran: DO 399 I=2,NM1 (1-based) → 0-based offset i=1..NM1-1
    for (int i = 1; i < NM1; i++) {
        double G = H;
        H = X[IBASE+i+1] - X[IBASE+i];
        double E = F;
        F = (Y[IBASE+i+1] - Y[IBASE+i]) / H;
        double GBY3 = G / 3.0;
        D[IBASE+i-1] = GBY3 * B[IBASE+i-1];
        double EPSIM1 = G + H;
        double RIM1B3 = F - E;
        B[IBASE+i] = 1.0 / ((2.0/3.0)*EPSIM1 - GBY3*D[IBASE+i-1]);
        C[IBASE+i] = RIM1B3 - D[IBASE+i-1]*C[IBASE+i-1];
    }
    D[IBASE+NM1] = 0.0;

    // ---- Backward sweep ----
    // Fortran: DO 500 I1=2,NM1 (1-based): I=NM1+2-I1 (1-based)
    // At I1=2: I=NM1 (1-based) = NM1-1 (0-based)
    // At I1=NM1: I=2 (1-based) = 1 (0-based)
    // 0-based i = I_1based - 1 = (NM1+2-I1) - 1 = NM1+1-I1
    // I1 goes 2..NM1 → i goes NM1-1..1
    for (int I1 = 2; I1 <= NM1; I1++) {
        int i = NM1 + 1 - I1;   // 0-based: NM1-1 down to 1
        C[IBASE+i] = B[IBASE+i]*C[IBASE+i] - D[IBASE+i]*C[IBASE+i+1];
    }

    // ---- Final B,C,D coefficients ----
    // Fortran: DO 699 I=1,NM1 (1-based) → 0-based i=0..NM1-1
    for (int i = 0; i < NM1; i++) {
        H = X[IBASE+i+1] - X[IBASE+i];
        D[IBASE+i] = (C[IBASE+i+1] - C[IBASE+i]) / (3.0*H);
        B[IBASE+i] = (Y[IBASE+i+1]-Y[IBASE+i])/H - (H*D[IBASE+i] + C[IBASE+i])*H;
    }

    // ---- Extend N-th cubic (= N-1-th cubic at Fortran IBASE+NUSE index) ----
    // Fortran: D(IBASE+NUSE) = D(IBASE+NM1); B(IBASE+NUSE) = B(IBASE+NM1) + ...
    // 0-based: index IBASE+NM1 (= IBASE+NUSE-1)
    // Fortran IBASE+NUSE (1-based) = IBASE+NUSE-1 (0-based) = IBASE+NM1
    // Fortran IBASE+NM1  (1-based) = IBASE+NM1-1 (0-based) = IBASE+NM1-1
    D[IBASE+NM1]     = D[IBASE+NM1-1];
    B[IBASE+NM1]     = B[IBASE+NM1-1] + (2.0*C[IBASE+NM1-1] + 3.0*D[IBASE+NM1-1]*H)*H;
    // H is from the last iteration of DO 699: interval i=NM1-1, i.e. X[IBASE+NM1]-X[IBASE+NM1-1]
}

// =============================================================
// Ptolemy INTRPC — evaluate spline at new points
// Exact port of fortlib.mor lines 2313-2420.
// Inputs: NCUBIC, XCUBES[NCUBIC], AS/BS/CS/DS[NCUBIC], NPTS, XS[NPTS]
// Output: YS[NPTS]
// =============================================================
void Intrpc(int NCUBIC, const std::vector<double>& XCUBES,
            const std::vector<double>& AS, const std::vector<double>& BS,
            const std::vector<double>& CS, const std::vector<double>& DS,
            int NPTS, const std::vector<double>& XS, std::vector<double>& YS) {
    YS.resize(NPTS, 0.0);
    if (NCUBIC < 2) return;

    const double INF = 1.0e300;
    double XSIGN = (XCUBES[1] > XCUBES[0]) ? 1.0 : -1.0;

    double XPREV = INF;
    double XNEXT = INF;
    double XBASE = 0.0;
    int    N = 0;
    double A = 0, B = 0, CC = 0, D = 0;

    for (int ii = 0; ii < NPTS; ii++) {
        double x = XS[ii];
        double xc = XSIGN * x;

loop100:
        if (xc < XNEXT) goto loop200;

        // Advance to next spline
        XPREV = XNEXT;
        XBASE = XNEXT * XSIGN;
        N = N + 1;
        XNEXT = XCUBES[N+1] * XSIGN;
        if (N == NCUBIC - 1) XNEXT = INF;

loop150:
        A  = AS[N];
        B  = BS[N];
        CC = CS[N];
        D  = DS[N];
        goto loop100;

loop200:
        if (xc >= XPREV) goto loop300;

        // Re-initialize (non-monotonic XS)
        N = 0;
        XPREV = -INF;
        XNEXT = XCUBES[1] * XSIGN;
        XBASE = XCUBES[0];
        goto loop150;

loop300:
        double DEL = x - XBASE;
        YS[ii] = A + DEL*(B + DEL*(CC + DEL*D));
    }
}

// ---- GaussLegendre ----
static void GaussLegendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n); w.resize(n);
    if (n==1){x[0]=0;w[0]=2;return;}
    int m=(n+1)/2;
    double e1=n*(n+1),e2=M_PI/(4*n+2),e3=1.0-(1.0-1.0/n)/(8.0*n*n);
    for(int i=1;i<=m;i++){
        double t=(4*i-1)*e2,x0=e3*std::cos(t),pkm1=1,pk=x0,fk=1;
        for(int k=2;k<=n;k++){fk++;double t1=x0*pk,pkp1=t1-pkm1-(t1-pkm1)/fk+t1;pkm1=pk;pk=pkp1;}
        double den=1-x0*x0,d1=n*(pkm1-x0*pk),dpn=d1/den;
        double d2pn=(2*x0*dpn-e1*pk)/den,d3pn=(4*x0*d2pn+(2-e1)*dpn)/den,d4pn=(6*x0*d3pn+(6-e1)*d2pn)/den;
        double u=pk/dpn,v=d2pn/dpn,h1=1+0.5*u*(v+u*(v*v-u*d3pn/(3*dpn))),h=-u*h1;
        double p=pk+h*(dpn+0.5*h*(d2pn+h/3*(d3pn+0.25*h*d4pn)));
        double dp=dpn+h*(d2pn+0.5*h*(d3pn+h*d4pn/3));h=h-p/dp;
        x[n-i]=x0+h;x[i-1]=-(x0+h);
        double fx=d1-h*e1*(pk+0.5*h*(dpn+h/3*(d2pn+0.25*h*(d3pn+0.2*h*d4pn))));
        w[n-i]=2*(1-x[i-1]*x[i-1])/(fx*fx);w[i-1]=w[n-i];
    }
    if((m+m)>n) x[m-1]=0;
}

// ---- CubMap ----
static void CubMap(int maptyp, double xlo, double xmid, double xhi, double gamma,
                   std::vector<double>& args, std::vector<double>& wts) {
    int npts=(int)args.size();
    double tau=(gamma>1e-6)?std::log(gamma+std::sqrt(gamma*gamma+1.0)):gamma*(1.0-gamma*gamma/6.0);
    double xlen=xhi-xlo,xadd=xlo+xhi;
    if(maptyp==0){for(int i=0;i<npts;i++){args[i]=xlo+0.5*xlen*(args[i]+1.0);wts[i]*=0.5*xlen;}}
    else if(maptyp==1){xmid=std::max(xmid,xlo+xlen/7.0);xmid=std::min(xmid,0.5*xadd);double A=0.5*xadd-xmid,B=0.5*xlen,C=0.5*xadd;for(int i=0;i<npts;i++){double tu=tau*args[i],xi=std::sinh(tu)/gamma;args[i]=A*(xi*xi-1.0)*(xi+1.0)+B*xi+C;wts[i]*=(tau/gamma)*std::cosh(tu)*((3.0*xi-1.0)*(xi+1.0)*A+B);}}
    else if(maptyp==2){double A=-xmid*xlen,B=xlen,C=xmid*xadd-2.0*xlo*xhi,D=xadd-2.0*xmid;for(int i=0;i<npts;i++){double tu=tau*args[i],sh=std::sinh(tu),denom=B-(D/gamma)*sh;args[i]=(-A+(C/gamma)*sh)/denom;wts[i]*=(tau/gamma)*std::cosh(tu)*((B*C-A*D)/(denom*denom));}}
    else if(maptyp==3){for(int i=0;i<npts;i++){double tu=tau*args[i];args[i]=xlo+0.5*xlen*(std::sinh(tu)/gamma+1.0);wts[i]*=0.5*xlen*(tau/gamma)*std::cosh(tu);}}
}

int main() {
    // -------------------------------------------------------
    // Test 1: cubic y=x^3-2x^2+1, N=7 uniform, Nout=11
    // -------------------------------------------------------
    {
        int N=7, NOUT=11;
        std::vector<double> X(N),Y(N),B(N),C(N),D(N);
        for(int i=0;i<N;i++){X[i]=i*0.5;Y[i]=X[i]*X[i]*X[i]-2*X[i]*X[i]+1;}
        Splncb(N,X,Y,B,C,D);
        std::vector<double> XS(NOUT),YS;
        for(int i=0;i<NOUT;i++) XS[i]=i*0.25;
        Intrpc(N,X,Y,B,C,D,NOUT,XS,YS);
        std::cout<<"=C1= Cubic N=7 Nout=11\n";
        for(int i=0;i<NOUT;i++){
            double yex=XS[i]*XS[i]*XS[i]-2*XS[i]*XS[i]+1;
            std::cout<<std::fixed<<std::setprecision(6)<<XS[i]<<"  "
                     <<std::setprecision(9)<<YS[i]<<"  "<<yex<<"  "
                     <<std::scientific<<std::setprecision(3)<<YS[i]-yex<<"\n";
        }
    }

    // -------------------------------------------------------
    // Test 2: Gaussian exp(-x^2), non-uniform N=10, Nout=19
    // -------------------------------------------------------
    {
        int N=10, NOUT=19;
        std::vector<double> X={0.0,0.2,0.5,0.8,1.1,1.5,2.0,2.4,2.8,3.0};
        std::vector<double> Y(N),B(N),C(N),D(N);
        for(int i=0;i<N;i++) Y[i]=std::exp(-X[i]*X[i]);
        Splncb(N,X,Y,B,C,D);
        std::vector<double> XS(NOUT),YS;
        for(int i=0;i<NOUT;i++) XS[i]=i*3.0/18.0;
        Intrpc(N,X,Y,B,C,D,NOUT,XS,YS);
        std::cout<<"=C2= Gaussian N=10 Nout=19\n";
        for(int i=0;i<NOUT;i++){
            double yex=std::exp(-XS[i]*XS[i]);
            std::cout<<std::fixed<<std::setprecision(6)<<XS[i]<<"  "
                     <<std::setprecision(10)<<YS[i]<<"  "<<yex<<"  "
                     <<std::scientific<<std::setprecision(3)<<YS[i]-yex<<"\n";
        }
    }

    // -------------------------------------------------------
    // Test 3: INELDC realistic spline NPSUM=40 → NPSUMI=42
    // y = U^2*exp(-0.3*U^2) on CUBMAP(rational-sinh) grid
    // -------------------------------------------------------
    {
        int NPSUM=40, NPSUMI=42;
        double SUMMIN=0.1, SUMMAX=15.0, SUMMID=2.5, GAMSUM=2.0;
        std::vector<double> X(NPSUM), W(NPSUM), Y(NPSUM), B(NPSUM), C(NPSUM), D(NPSUM);
        GaussLegendre(NPSUM,X,W);
        double xm=SUMMID; CubMap(2,SUMMIN,xm,SUMMAX,GAMSUM,X,W);
        for(int i=0;i<NPSUM;i++) Y[i]=X[i]*X[i]*std::exp(-0.3*X[i]*X[i]);
        Splncb(NPSUM,X,Y,B,C,D);

        std::vector<double> XS(NPSUMI),WS(NPSUMI),YS;
        GaussLegendre(NPSUMI,XS,WS);
        double xm2=SUMMID; CubMap(2,SUMMIN,xm2,SUMMAX,GAMSUM,XS,WS);
        Intrpc(NPSUM,X,Y,B,C,D,NPSUMI,XS,YS);
        std::cout<<"=C3= INELDC spline NPSUM=40 NPSUMI=42\n";
        for(int i=0;i<NPSUMI;i++){
            double yex=XS[i]*XS[i]*std::exp(-0.3*XS[i]*XS[i]);
            double rel=(std::abs(yex)>1e-30)?(YS[i]-yex)/yex:0.0;
            std::cout<<std::fixed<<std::setprecision(6)<<XS[i]<<"  "
                     <<std::scientific<<std::setprecision(9)<<YS[i]<<"  "<<yex<<"  "
                     <<std::setprecision(3)<<rel<<"\n";
        }
    }

    {
        int N=10;
        std::vector<double> X2={0.0,0.2,0.5,0.8,1.1,1.5,2.0,2.4,2.8,3.0};
        std::vector<double> Y2(N),B2(N),C2(N),D2(N);
        for(int i=0;i<N;i++) Y2[i]=std::exp(-X2[i]*X2[i]);
        Splncb(N,X2,Y2,B2,C2,D2);
        std::cout<<"=DEBUG_BCD=\n";
        for(int i=0;i<N;i++) std::cout<<"i="<<i<<" B="<<std::setprecision(12)<<B2[i]<<" C="<<C2[i]<<" D="<<D2[i]<<"\n";
    }
    return 0;
}

// Quick debug entry point: print B,C,D coefficients for T2
void debug_t2() {
    int N=10;
    std::vector<double> X={0.0,0.2,0.5,0.8,1.1,1.5,2.0,2.4,2.8,3.0};
    std::vector<double> Y(N),B(N),C(N),D(N);
    for(int i=0;i<N;i++) Y[i]=std::exp(-X[i]*X[i]);
    Splncb(N,X,Y,B,C,D);
    std::cout<<"=DEBUG T2 B,C,D=\n";
    for(int i=0;i<N;i++)
        std::cout<<"  i="<<i<<" X="<<X[i]<<" Y="<<std::setprecision(12)<<Y[i]
                 <<" B="<<B[i]<<" C="<<C[i]<<" D="<<D[i]<<"\n";
}
