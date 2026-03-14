#ifndef JSYMBOLS_C
#define JSYMBOLS_C

#include <stdlib.h>
#include <cmath>
#include <algorithm>
 
double factorial(double n){
  if( n < 0 ) return -100.;
  return (n == 1. || n == 0.) ? 1. : factorial(n-1) * n ;
}
 
double CGcoeff(double J, double m, double J1, double m1, double J2, double m2){
  // (J1,m1) + (J2, m2) = (J, m)
 
  if( m != m1 + m2 ) return 0;
 
  double Jmin = abs(J1 - J2);
  double Jmax = J1+J2;
 
  if( J < Jmin || Jmax < J ) return 0;
 
  double s0 = (2*J+1.) * factorial(J+J1-J2) * factorial(J-J1+J2) * factorial(J1+J2-J) / factorial(J+J1+J2 + 1.);
  s0 = sqrt(s0);
 
  double s = factorial(J +m ) * factorial(J -m);
  double s1 = factorial(J1+m1) * factorial(J1-m1);
  double s2 = factorial(J2+m2) * factorial(J2-m2);
  s = sqrt(s * s1 * s2);
 
  //printf(" s0, s = %f , %f \n", s0, s);
 
  int kMax = std::min( std::min( J1+J2-J, J1 - m1), J2 + m2);
 
  double CG = 0.;
  for( int k = 0; k <= kMax; k++){
    double k1 = factorial(J1+J2-J-k);
    double k2 = factorial(J1-m1-k);
    double k3 = factorial(J2+m2-k);
    double k4 = factorial(J - J2 + m1 +k);
    double k5 = factorial(J - J1 - m2 +k);
    double temp = pow(-1, k) / (factorial(k) * k1 * k2 * k3 * k4 * k5);
    if( k1 == -100. || k2 == -100. || k3 == -100. || k4 == -100. || k5 == -100. ) temp = 0.;
 
    //printf(" %d | %f \n", k, temp);
    CG += temp;
  }
 
  return s0*s*CG;
 
}
 
double ThreeJSymbol(double J1, double m1, double J2, double m2, double J3, double m3){
 
  // ( J1 J2 J3 ) = (-1)^(J1-J2 - m3)/ sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2) 
  // ( m1 m2 m3 )
  
  // The phase (-1)^(J1-J2-m3) is only well-defined when J1-J2-m3 is an integer.
  // If not (e.g., J1 half-int, J2,m3 integer), the 3J is zero by parity.
  double phase_exp = J1 - J2 - m3;
  double phase_exp_rounded = std::round(phase_exp);
  if( std::abs(phase_exp - phase_exp_rounded) > 1e-6 ) return 0.0;  // non-integer exponent → 3J = 0
  
  int phase_int = (int)phase_exp_rounded;
  double sign = (phase_int % 2 == 0) ? 1.0 : -1.0;
  
  return sign / sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2);
 
}
 
double SixJSymbol(double J1, double J2, double J3, double J4, double J5, double J6){
 
  // coupling of j1, j2, j3 to J-J1
  // J1 = j1
  // J2 = j2
  // J3 = j12 = j1 + j2
  // J4 = j3
  // J5 = J = j1 + j2 + j3
  // J6 = j23 = j2 + j3
 
  //check triangle condition
  if( J3 < abs(J1 - J2 ) || J1 + J2 < J3 ) return 0; 
  if( J6 < abs(J2 - J4 ) || J2 + J4 < J6 ) return 0;
  if( J5 < abs(J1 - J6 ) || J1 + J6 < J5 ) return 0;
  if( J5 < abs(J3 - J4 ) || J3 + J4 < J5 ) return 0;
 
  // Use half-integer index loops to avoid float precision issues.
  // Convert J values to doubled integers: J -> 2J (exact integers)
  int tJ1 = (int)std::round(2*J1), tJ2 = (int)std::round(2*J2), tJ3 = (int)std::round(2*J3);
  int tJ4 = (int)std::round(2*J4), tJ5 = (int)std::round(2*J5), tJ6 = (int)std::round(2*J6);
  // f = sum of (Ji - mi), each term is integer since Ji and mi are both half-integers
  // Total sum f is always integer when all m-conditions are satisfied.
  // Sum of all 2*Ji - sum of 2*mi_in = 2*(J1+J2+J3+J4+J5+J6) when all m-sums vanish.
  // f_total_doubled = tJ1+tJ2+tJ3+tJ4+tJ5+tJ6 (all doubled, integer)
  int f_total_doubled = tJ1 + tJ2 + tJ3 + tJ4 + tJ5 + tJ6;
  // sign = (-1)^(f_total) = (-1)^((tJ1+...+tJ6)/2) -- but this is pre-multiplication sign
  // We'll compute the sign per-term using integer arithmetic.

  double sixJ = 0;
  // Loop with integer 2*m values to avoid float precision
  for( int tm1 = -tJ1; tm1 <= tJ1; tm1 += 2){
    double m1 = tm1 * 0.5;
    for( int tm2 = -tJ2; tm2 <= tJ2; tm2 += 2){
      double m2 = tm2 * 0.5;
      for( int tm3 = -tJ3; tm3 <= tJ3; tm3 += 2){
        double m3 = tm3 * 0.5;
        for( int tm4 = -tJ4; tm4 <= tJ4; tm4 += 2){
          double m4 = tm4 * 0.5;
          for( int tm5 = -tJ5; tm5 <= tJ5; tm5 += 2){
            double m5 = tm5 * 0.5;
            for( int tm6 = -tJ6; tm6 <= tJ6; tm6 += 2){
              double m6 = tm6 * 0.5;

              // f = (J1-m1)+(J2-m2)+(J3-m3)+(J4-m4)+(J5-m5)+(J6-m6)
              // In doubled ints: f_doubled = (tJ1-tm1+tJ2-tm2+...+tJ6-tm6)/2*2 = tJ1-tm1+...
              // f is always an integer since Ji+mi is integer or 2*(half-int)
              int f2 = (tJ1 - tm1) + (tJ2 - tm2) + (tJ3 - tm3) + (tJ4 - tm4) + (tJ5 - tm5) + (tJ6 - tm6);
              // f2 = 2*f, so (-1)^f = (-1)^(f2/2). f2 is always even (since J-m is integer for any J,m pair)
              // Actually f2/2 = f. Check: for integer J, J-m = integer; for half-int J, J-m = integer too.
              // f2 is the sum of (tJ-tm), each of which equals 2*(Ji-mi) = even integer, so f2 is divisible by 2? No.
              // Wait: for half-int J=1.5, tJ=3, m=0.5, tm=1: tJ-tm=2 (even). For m=-0.5, tm=-1: tJ-tm=4 (even).
              // For integer J=2, tJ=4, m=1, tm=2: tJ-tm=2 (even). For m=0, tm=0: tJ-tm=4 (even).
              // So ALL (tJ-tm) values are even! Therefore f2 is divisible by 2.
              // => f = f2/2 is always an integer. Use (-1)^f = 1 if f%2==0, -1 if f%2==1.
              int f_int = f2 / 2;
              double sign = (f_int % 2 == 0) ? 1.0 : -1.0;

              double a1 = ThreeJSymbol( J1, -m1, J2, -m2, J3, -m3);
              if( a1 == 0 ) continue;
              double a2 = ThreeJSymbol( J1, m1, J5, -m5, J6, m6);
              if( a2 == 0 ) continue;
              double a3 = ThreeJSymbol( J4, m4, J2, m2, J6, -m6);
              if( a3 == 0 ) continue;
              double a4 = ThreeJSymbol( J4, -m4, J5, m5, J3, m3);
              if( a4 == 0 ) continue;

              sixJ += sign * a1 * a2 * a3 * a4;
            }
          }
        }
      }
    }
  }
 
  return sixJ;
}
 
double NineJSymbol( double J1, double J2, double J3, double J4, double J5, double J6, double J7, double J8, double J9){
 
  // Standard formula (Varshalovich et al.):
  // {j1 j2 j3; j4 j5 j6; j7 j8 j9} = sum_x (-1)^(2x)(2x+1)
  //    * {j1 j2 j3; j6 j9 x} * {j4 j5 j6; j2 x j8} * {j7 j8 j9; x j1 j4}
  // x range: max(|j1-j9|, |j2-j6|, |j4-j8|) <= x <= min(j1+j9, j2+j6, j4+j8)
  double xMin = std::max( std::max( std::abs(J1-J9), std::abs(J2-J6) ), std::abs(J4-J8) );
  double xMax = std::min( std::min( J1+J9, J2+J6 ), J4+J8 );
 
  double nineJ = 0;
  for( double x = xMin; x <= xMax + 1e-9; x = x + 1.0){
    double s1 = SixJSymbol(J1, J2, J3, J6, J9, x);
    if( s1 == 0 || std::isnan(s1) ) continue;
    double s2 = SixJSymbol(J4, J5, J6, J2, x, J8);
    if( s2 == 0 || std::isnan(s2) ) continue;
    double s3 = SixJSymbol(J7, J8, J9, x, J1, J4);
    if( s3 == 0 || std::isnan(s3) ) continue;
    double f = pow(-1.0, 2.0*x) * (2.0*x+1.0);
    nineJ += f * s1 * s2 * s3;
  }
 
  return nineJ;
}

#endif