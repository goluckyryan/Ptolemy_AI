#ifndef AV18_C
#define AV18_C

#include "constant.h"

#include <cstdio>

class AV18 {

public:
  AV18(): b(4.27), 
          gamma(0.577216), 
          beta(0.0189), 
          b3(pow(b,3)), 
          alphaHbarC3(FINE_STRCUTURE_CONSTANT * pow(HBARC, 3)),
          fsq(0.075),
          cpi(2.1),
          rws(0.5),
          aiws(0.2),
          mass_pi( (MASS_PI0 + 2*MASS_PIC)/3.),
          mu0(MASS_PI0/HBARC),
          muC(MASS_PIC/HBARC){
    mu = mass_pi/HBARC;
  }

  void CalNN(double r); // strong NN interaction
  void CalEM(double r); // EM interaction
  void CalPW(double r, double L, double S, double J, double T, double T1z, double T2z); // partial wave

  double GetNN(unsigned short index){return NN[index];}
  double GetEM(unsigned short index){return EM[index];}
  double * GetPW() {return *PW;}

private:

  const double b;
  const double gamma;
  const double beta;

  const double b3;
  const double alphaHbarC3;

  const double fsq; // square of coupling constant
  const double cpi;
  const double rws, aiws; // Woods-Saxon parameter
  const double mass_pi; // average of pion mass
  const double mu0, muC;
  double mu;

  double NN[18]; // store the 18 values for strong force;
  double EM[14]; //store the 14 values for EM components;

  double PW[2][2];

};


inline void AV18::CalEM(double r){ 
 double br = b * r;
 double br2 = pow(br, 2);
 double br3 = pow(br, 3);
 double br4 = pow(br, 4);
 double br5 = pow(br, 5);
 double expbr = exp(-br);

 double b3 = pow(b,3);
 
 double fCoul = b * 5. /16.;
 double ft = b3 * br2 / 720.;
 double fLS = b3 / 48.;
 double kr = MASS_ELECTRON * 1e-5 / HBARC;
 if (r > 1e-5 ){
   fCoul = (1 - (1 + 11./16. * br + 3./16. * br2 + br3 / 48.) * expbr )/r;
   ft = (1 - (1 + br + br2/2 + br3/6 + br4/24 + br5/144) * expbr ) / r /r / r;
   fLS = (1 - (1 + br + br2/2 + br3*7./48 + br4/48) * expbr ) / r /r / r;
   kr = MASS_ELECTRON * r / HBARC;
 }
 double fIVP = - gamma - 5./6. + abs(log(kr)) + 6./8 * M_PI * kr;
 double fDelta = b3 * (3 + 3 * br + br2) * expbr / 48.;
 double fnpr = b3 * (15 + 15 * br + 6 * br2 + br3) * expbr / 384.;

 double v1 = FINE_STRCUTURE_CONSTANT * HBARC * fCoul;

 /*  C1(pp)               */ EM[0] = v1;
 /*  DF(pp), Darwin-Foldy */ EM[1] = - alphaHbarC3 * fDelta / 4. / pow(MASS_PROTON, 2);
 /*  C2(pp)               */ EM[2] = - pow(v1, 2)/MASS_PROTON;
 /*  VP(pp), Vaccum pol.  */ EM[3] = 2 * FINE_STRCUTURE_CONSTANT * v1 * fIVP / 3. / M_PI;
 /*  C1(np)               */ EM[4] = FINE_STRCUTURE_CONSTANT * HBARC * beta * fnpr;
 /*  s1.s2(pp)            */ EM[5] = - alphaHbarC3 * pow(MU_PROTON/MASS_PROTON, 2) * fDelta / 6;
 /*  s1.s2(nn)            */ EM[6] = - alphaHbarC3 * pow(MU_NEUTRON/MASS_NEUTRON, 2) * fDelta / 6;
 /*  s1.s2(np)            */ EM[7] = - alphaHbarC3 * MU_NEUTRON * MU_PROTON * fDelta / 6. / MASS_NEUTRON / MASS_PROTON;
 /*  S12(pp)              */ EM[8] = - alphaHbarC3 * pow(MU_PROTON/MASS_PROTON,2) * ft / 4.;
 /*  S12(nn)              */ EM[9] = - alphaHbarC3 * pow(MU_NEUTRON/MASS_NEUTRON,2) * ft / 4.;
 /*  S12(np)              */ EM[10] = - alphaHbarC3 *  MU_NEUTRON * MU_PROTON * ft / 4. / MASS_NEUTRON / MASS_PROTON;
 /*  L.S(pp)              */ EM[11] = - alphaHbarC3 * ( 4 * MU_PROTON - 1) * fLS / 2 / MASS_PROTON / MASS_PROTON;
 /*  L.S(nn)              */ EM[12] = 0;
 /*  L.S(np)              */ EM[13] = - alphaHbarC3 * MU_NEUTRON * fLS / 2 / MASS_NEUTRON / MASS_NUCLEON;

}

inline void AV18::CalNN(double r){
  double mur = mu * r;
  double expur = exp(- mur );
  double haha = 1 - exp( - cpi * r * r);
  double mu0r = mu0 * r;
  double muCr = muC * r;

  double ypi = 0 ;
  if( r > 0 ) ypi = expur * haha / mur;

  double tpi  = 3 * pow(cpi,2) * r / pow(mu,3);
  double ypi0 = pow( MASS_PI0 / MASS_PIC, 2) * (MASS_PI0 / 3.) * cpi / mu0 * r;
  double tpi0 = 3 * cpi * ypi0 / pow(mu0, 2);
  double ypiC = (MASS_PIC / 3.) * ( cpi / muC ) * r;
  double tpiC = 3 * cpi * ypiC / pow(muC, 2);
  if( r > 1e-4 ){
    tpi  = (1 +  3. / mur + 3. / pow(mur,2) ) * ypi * haha;
    ypi0 = pow( MASS_PI0 / MASS_PIC, 2) * (MASS_PI0 / 3.) * haha * exp(- mu0r ) / mu0r;
    tpi0 = (1 + 3. / mu0r +  3. / pow(mu0r, 2) ) * ypi0 * haha;
    ypiC = (MASS_PIC / 3.) * exp( - muCr )  * haha / muCr;
    tpiC = (1 + 3. / muCr +  3. / pow(muCr, 2) ) * ypiC * haha;
  }

  // printf("%.10f, %.10f, %.10f, %.10f, %.10f, %.10f\n", ypi, ypi0, ypiC, tpi,  tpi0, tpiC);

  ypi0 = fsq * ypi0;
  ypiC = fsq * ypiC;
  tpi0 = fsq * tpi0;
  tpiC = fsq * tpiC;
  double tpi2 = pow( tpi, 2); // 2 pion exchange

  double ws = 1. / (1 + exp( (r - rws)/aiws ) );
  double ws0 = 1. / (1 + exp( - rws/aiws ) );
  double wsp = ws * ( 1 + ws0 * r / aiws * exp( - rws / aiws) );
  double wsx = ws * mur;
  double wsx2 = ws * pow( mur, 2);

  double dypi00 = pow( MASS_PI0/MASS_PIC ,2) * (MASS_PI0 / 3.) * cpi / mu0;
  double dypiC0 = (MASS_PIC / 3.) * cpi / muC;
  double ypi0p = ypi0 - fsq * dypi00 * ws / ws0 * r;
  double ypiCp = ypiC - fsq * dypiC0 * ws / ws0 * r;

  // printf("%.10f, %.10f\n", ypi0, dypi00);

  // ypi = (ypi0 + 2 * ypiC) / 3.;
  // tpi = (tpi0 + 2 * tpiC) / 3.;

  // S = 1 , T = 1
  double p11pp = -7.62701 * tpi2 + 1815.4920 * wsp + 1847.8059 * wsx2 + ypi0p; // central pp, v_{11,pp}^c 
  double p11np = -7.62701 * tpi2 + 1813.5315 * wsp + 1847.8059 * wsx2 - ypi0p + 2 * ypiCp; // central np, v_{11,np}^c 
  double p11nn = -7.62701 * tpi2 + 1811.5710 * wsp + 1847.8059 * wsx2 + ypi0p; // central  nn, v_{11,nn}^c 
  double pt1pp =  1.07985 * tpi2 -  190.0949 * wsx -  811.2040 * wsx2 + tpi0; //  tensor,  pp, v_{11,pn}^(S12} or v_{11,pn}^(t} 
  double pt1np =  1.07985 * tpi2 -  190.0949 * wsx -  811.2040 * wsx2 - tpi0 + 2 * tpiC; // tensor, np 
  double pt1nn =  1.07985 * tpi2 -  190.0949 * wsx -  811.2040 * wsx2 + tpi0; // tensor, nn 
  double pls1  = -0.62697 * tpi2 -  570.5571 * wsp +  819.1222 * wsx2;  // ls, v_{11,NN}^(ls} 
  double pl211 =  0.06709 * tpi2 +  342.0669 * wsp -  615.2339 * wsx2; // l^2, v_{11,NN}^(l^2} 
  double pls21 =  0.74129 * tpi2 +    9.3418 * wsp -  376.4384 * wsx2; // (ls)^2, v_{11,NN}^((ls)^2}

  // S=1, T=0 
  double p10   = -8.62770  * tpi2 + 2605.2682 * wsp + 441.9733 * wsx2 - ypi0p - 2 * ypiCp; // Central 
  double pt0   =  1.485601 * tpi2 - 1126.8359 * wsx + 370.1324 * wsx2 - tpi0 - 2 * tpiC;//  tensor 
  double pls0  =  0.10180  * tpi2 +   86.0658 * wsp - 356.5175 * wsx2; //  ls 
  double pl210 = -0.13201  * tpi2 +  253.4350 * wsp -   1.0076 * wsx2; //  l^2 
  double pls20 =  0.07357  * tpi2 -  217.5791 * wsp +  18.3935 * wsx2; // (ls)^2 

  // printf("%15.12f \n", ypi0p);
  // printf("%15.12f \n", ypiCp);

  // S=0, T=1 
  double p01pp = -11.27028 * tpi2 + 3346.6874 * wsp - 3 * ypi0p;// central pp 
  double p01np = -10.66788 * tpi2 + 3126.5542 * wsp - 3 * (-ypi0p + 2 * ypiCp); // central np 
  double p01nn = -11.27028 * tpi2 + 3342.7664 * wsp - 3 * ypi0p;// central nn 
  double pl201 =   0.12472 * tpi2 +   16.7780 * wsp ; // l^2 
 
  // S=0, T=0 
  double p00   = -2.09971 * tpi2 + 1204.4301 * wsp - 3 * (- ypi0p - 2 * ypiCp); // central 
  double pl200 = -0.31452 * tpi2 +  217.4559 * wsp; // l^2 

  // Projection 
  // S=1, T=1, central 
  double p11   = 1. / 3. * (p11pp + p11nn +  p11np); // charge independent, v_{11}^{CI} 
  double p11cd = 1. / 6. * ((p11pp + p11nn)/2 -  p11np); // charge dependent tensor, v_{11}^{CD} 
  double p11cs = 1. / 4. * (p11pp - p11nn); // charge asymmetric, v_{11}^{CA} 

  // S=1, T = 1, tensor
  double pt1   = 1. / 3. * (pt1pp + pt1nn + pt1np); // v_{11}^{t} 
  double pt1cd = 1. / 6. * ((pt1pp + pt1nn)/2. - pt1np);
  double pt1cs = 1. / 4. * (pt1pp - pt1nn);

  // S=0, T=1, cental
  double p01   = 1. / 3. * (p01pp + p01nn + p01np);
  double p01cd = 1. / 6. * ((p01pp + p01nn)/2. - p01np);
  double p01cs = 1. / 4. * (p01pp - p01nn);

  // printf("%15.12f \n", p11);
  // printf("%15.12f \n", p10);
  // printf("%15.12f \n", p01);
  // printf("%15.12f \n", p00);
  
  // printf("%15.12f \n", p01pp);
  // printf("%15.12f \n", p01nn);
  // printf("%15.12f \n", p01np);


  NN[0]  = 1. / 16. * (9 * p11 + 3 * p10 + 3 * p01 + p00); //1, 1 
  NN[1]  = 1. / 16. * (3 * p11 - 3 * p10 +     p01 - p00); //2, t1.t2 
  NN[2]  = 1. / 16. * (3 * p11 +     p10 - 3 * p01 - p00); //3, s1.s2 
  NN[3]  = 1. / 16. * (    p11 -     p10 -     p01 + p00); //4, (s1.s2)(t1.t2) 
  NN[4]  = 1. /  4. * (3 * pt1 + pt0); //5, S12 = 3(s1.r)(s2.r)-s1.s2 
  NN[5]  = 1. /  4. * (    pt1 - pt0); //6, S12(t1.t2) 
  NN[6]  = 1. /  4. * (3 * pls1 + pls0); //7, L.S 
  NN[7]  = 1. /  4. * (    pls1 - pls0); //8, L.S(t1.t2) 
  NN[8]  = 1. / 16. * (9 * pl211 + 3 * pl210 + 3 * pl201 + pl200); //9, L^2 
  NN[9]  = 1. / 16. * (3 * pl211 - 3 * pl210 +     pl201 - pl200); //10, L^2(t1.t2) 
  NN[10] = 1. / 16. * (3 * pl211 +     pl210 - 3 * pl201 - pl200); //11, L^2(s1.s2) 
  NN[11] = 1. / 16. * (    pl211 -     pl210 -     pl201 + pl200); //12, L^2(s1.s2)(t1.t2) 
  NN[12] = 1. /  4. * (3 * pls21 + pls20); //13, (L.S)^2 
  NN[13] = 1. /  4. * (    pls21 - pls20); //14, (L.S)^2 (t1.t2) 
  NN[14] = 1. /  4. * (3 * p11cd + p01cd); //15, T12 
  NN[15] = 1. /  4. * (    p11cd - p01cd); //16, T12(s1.s2) 
  NN[16] = pt1cd; //17, S12 T12 
  NN[17] = p01cs; //18, t1z + t2z 
 
}

inline void AV18::CalPW(double r, double L, double S, double J, double T, double T1z, double T2z){
  
  CalNN(r);
  double s1ds2 = 4 * S - 3; // s1.s2
  double t1dt2 = 4 * T - 3; // t1.t2

  double t12 = 3 * T1z * T2z - t1dt2;
  double vc = NN[0] + t1dt2 * NN[1] + s1ds2 * NN[2] + s1ds2 * t1dt2 * NN[3] + t12 * NN[14] + s1ds2 * t12 * NN[15] + (T1z + T2z) * NN[16];
  double vt = NN[4] + t1dt2 * NN[5] + t12 * NN[16];
  double vls = NN[6] + t1dt2 * NN[7];
  double vl2 = NN[8] + t1dt2 * NN[9] + s1ds2 * NN[10] + s1ds2 * t1dt2 * NN[11];
  double vls2 = NN[12] + t1dt2 * NN[13];

  CalEM(r);

  if( T1z + T2z == - 1) { // nn 
    vc += s1ds2 * EM[6];
    vt += EM[9];
  }

  if( T1z + T2z == 0) { // np 
    vc  +=  EM[4] + s1ds2 * EM[7];
    vt  +=  EM[10];
    vls +=  EM[13];
  };
  
  if( T1z + T2z == 1) { // pp 
    vc  += EM[0] + EM[1] + EM[2] + EM[3] +  s1ds2 * EM[5];
    vt  += EM[8];
    vls += EM[11];
  };

  if(S == 1 && J == L + 1) {// coupled channel 
    double s12m = -2 * (J - 1)/(2 * J + 1); 
    double s12  = 6 * sqrt(J * (J + 1))/(2 * J + 1); 
    double s12p = -2 * (J + 2)/(2 * J + 1);
    double lsm = J - 1; 
    double lsp = - J - 2;

    PW[0][0] = vc + s12m * vt + lsm * vls + L * (L + 1) * vl2 + lsm * lsm * vls2;
    PW[0][1] = s12 * vt;
    PW[1][0] = s12 * vt; 
    PW[1][1] = vc + s12p * vt + lsp * vls + (L + 2) * (L + 3) * vl2 + lsp * lsp * vls2;

  }else{ //(* single channel *)
    double s12 = 0;
    if(S == 1 && L == J) s12 = 2;
    if(L == J + 1) s12 = -2. * (J + 2)/(2 * J + 1);
 
    double ls = 0.5 * (J * (J + 1.) - L * (L + 1.) - S * (S + 1.));

    PW[0][0] = vc + s12 * vt + ls * vls + L * (L + 1) * vl2 + ls * ls * vls2;
    PW[0][1] = 0 ;
    PW[1][0] = 0 ;
    PW[1][1] = 0 ;
  };

}

#endif