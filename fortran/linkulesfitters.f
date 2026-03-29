c*id* aitken     linkulesfitters  collection
      SUBROUTINE AITKEN ( NDEGRE, DELTA, EPSLON, NIN, XIN, FIN,
     1  NOUT, XOUT, FOUT, NFAIL, WORST, WORK )
C
C     AITKEN-LAGRANGE INTERPOLATION ON A NONUNIFORM GRID
C
C
C     The nearest points that surround each X are used in
C     interpolating to X.  The first point is always the one that
C     is closest to X and the second is always the nearest point on the    10   
C     other side of X.  The third point is choosen to be the next
C     nearest point that is closer to X( (NIN+1)/2 ) (i.e., the center
C     of the array.  Thereafter points are choosen on alternating sides
C     of X.  Thus if the input X's (XIN) were
C       1, 1/2, 1/4, 1/8, 1/16, ..., 1/1024
C     and it was desired to interpolate to 1/5, the points used
C     (in order of increasing degree of polynomial) would be
C       1/4, 1/8, 1/16, 1/2, 1/32, 1, 1/64, 1/128, 1/256, ...
C     This order of choosing points guarantees that (except for
C     NDEGRE = 0) there will be no discontinuities in the                  20   
C     interpolated function.  However, discontinuities of the first
C     derivative can occur no matter what order is used.
C        X NEED NOT BE IN THE DOMAIN SPANNED BY THE DEFINITION OF F;
C     BUT THE USE OF POLYNOMIALS TO EXTRAPOLATE IS USUALLY NOT
C     SUCCESSFUL.
C
C
C     NDEGRE - THE ORDER ( 0 =< NDEGRE ) OF THE LAGRANGE INTERPOLATION
C          POLYNOMIAL TO USE OR, IF DELTA > 0, THE MAXIMUM ORDER TO USE.
C                                                                          30   
C     DELTA, EPSLON - SPECIFY THE ACCURACY WITH WHICH THE INTERPOLATION
C          IS TO BE DONE.  THE ITERATIONS CONTINUE UNTIL ONE OF THE
C          FOLLOWING OCCUR:
C               1)  |F(I)-F(I-1)|/|F(I)|  =<  DELTA
C               2)  |F(I)|  =<  EPSLON
C               3)  I = NDEGRE
C          WHERE I IS THE DEGREE OF THE LAGRANGE POLYNOMIAL.
C          IF  DELTA =< 0, THEN NO CONVERGENCE TESTS ARE MADE AND
C          NDEGRE POLYNOMIALS ARE ALWAYS USED.
C                                                                          40   
C     NIN - NUMBER OF INPUT POINTS.  NIN >= NDEGRE+1 IS REQUIRED.
C     XIN, FIN - ARRAYS OF DIMENSION NIN THAT SPECIFY THE INPUT.
C          XIN MUST BE MONOTONIC (EITHER INCREASING OR DECREASING) AND
C          AND DOES NOT NEED TO BE EQUAL SPACED.
C          FIN MUST CONTAIN
C             FIN(I) = F(XIN(I))    1 =< I =< NIN.
C
C     NOUT - NUMBER OF DESIRED INTERPOLATIONS.
C     XOUT - ARRAY OF DIMENSION NOUT THAT CONTAINS THE X'S TOWHICH
C          INTERPOLATION IS TO BE DONE.  IT CAN HAVE ANY ORDER OR          50   
C          DISORDER.
C     FOUT - ARRAY OF DIMENSION NOUT THAT WILL BE SET TO
C          FOUT(I) = F(XOUT(I))  1 =< I =< NOUT.
C
C     NFAIL - WILL BE SET TO THE NUMBER OF FOUT'S THAT DID NOT CONVERGE
C          TO ACCURACY DELTA OR EPSLON.  FOR DELTA =< 0, NFAIL = 0.
C     WORST - WILL BE SET TO THE LARGEST RELATIVE DIFERENCE IN THE
C          LAST TWO ITERATIONS OF THE INTERPOLATION; I.E. WORST WILL
C          BE THE MAXIMUM FINAL NUMBER USED IN THE CONVERGENCE TEST.
C          IF DELTA =< 0, WORST IS SET TO 0.                               60   
C     BOTH NFAIL AND WORST ARE SET TO 0 IF NDEGRE IS 0 (EVEN IF
C     DELTA > 0.
C
C     WORK - A WORK ARRAY OF DIMENSION ATLEAST 2*(NDEGRE+1).
C
C     MAY 15, 1976 - FIRST VERSION - S. PIEPER
C     NOV 5, 1976 - ALWAYS USE SURROUNDING POINTS
C
C
      IMPLICIT  REAL*8 ( A-H, O-Z )                                        70   
C
      DIMENSION  XIN(NIN), FIN(NIN), XOUT(NOUT), FOUT(NOUT), WORK(2,1)
C
      LOGICAL DELTSW, UPSW
C
      DELTSW = DELTA .GT. 0
C
      NI = 1
      NFAIL = 0
      WORST = 0                                                            80   
      XMID = XIN( (NIN+1)/2 )
      NDEGR = MIN0( NDEGRE, NIN-1 )
C
      DO 799  NO = 1, NOUT
         X = XOUT(NO)
C
C     FIND NEAREST POINT
C
         UPSW = .FALSE.
         D = DABS( X - XIN(NI) )                                           90   
 120     IF ( NI .EQ. NIN )  GO TO 150
            D2 = DABS( X - XIN(NI+1) )
            IF ( D2 .GE. D )  GO TO 150
            UPSW = .TRUE.
            NI = NI + 1
            D = D2
         GO TO 120
C
C     GETTING LARGER MAKES IT WORSE; HAVE WE PASSED A MINIMUM.
C                                                                         100   
 150     IF ( UPSW )  GO TO 200
 160     IF ( NI .EQ. 1 )  GO TO 200
            D2 = DABS( X - XIN(NI-1) )
            IF ( D2 .GE. D )  GO TO 200
            NI = NI - 1
            D = D2
         GO TO 160
C
C     HAVE NOW FOUND NEAREST XIN
C                                                                         110   
C     NL IS THE INDEX OF THE NEXT LOWER POINT TO USE.
C     NU IS THE INDEX OF THE NEXT UPPER POINT TO USE.
C
 200     NL = NI - 1
         NU = NI + 1
         F = FIN(NI)
         WORK(1,1) = X - XIN(NI)
         UPSW = WORK(1,1) .LT. 0
         WORK(2,1) = F
         ERROR = 0                                                        120   
         IF ( NDEGR .LE. 0 )  GO TO 700
C
C     DO UP TO NDEGR ITERATIONS
C
         DO 599  ND = 1, NDEGR
            FLAST = F
C
C     FIND NEXT POINT
C
C     BIAS CHOICE OF 3RD POINT TOWARDS CENTER                             130   
C
            IF ( ND .EQ. 2 )  UPSW = X .LT. XMID
            IF ( UPSW )  GO TO 350
C
C     LAST POINT WAS A LOWER - NOW USE UPPER
C
 320        I = NU
            NU = NU + 1
            UPSW = .TRUE.
            IF ( I .LE. NIN )  GO TO 400                                  140   
C
C     LAST POINT WAS AN UPPER - NOW USE LOWER
C
 350        I = NL
            NL = NL - 1
            UPSW = .FALSE.
            IF ( I .LT. 1 )  GO TO 320
C
C     NOW USE NEXT POINT
C                                                                         150   
 400        DEL = X - XIN(I)
            F = FIN(I)
            DO 489  J = 1, ND
               F = (WORK(1,J)*F - DEL*WORK(2,J))/(WORK(1,J)-DEL)
 489        CONTINUE
C
C     HAVE WE CONVERGED
C
            IF ( .NOT. DELTSW )  GO TO 550
            ERROR = DABS(F)                                               160   
            IF ( ERROR .LE. EPSLON )  GO TO 700
            ERROR = DABS( (F-FLAST)/F )
            IF ( ERROR .LE. DELTA )  GO TO 700
 550        WORK(1,ND+1) = DEL
            WORK(2,ND+1) = F
 599     CONTINUE
C
C     STOPPED BY THE DEGREE
C
         IF (DELTSW)  NFAIL = NFAIL + 1                                   170   
C
C     CONVERGED
C
 700     FOUT(NO) = F
         WORST = DMAX1( WORST, ERROR )
 799  CONTINUE
C
      RETURN
      END
!!!      subroutine AITKEN ( i, x, y, LNSTRT, R, V,                       180   
!!!     1  NSTRT2, ALLOC, ARRAY1,  NFAIL, WORST, A )
!!!c
!!!c    fake aitken for unused linkules
!!!c
!!!        write (6,*) '!!!AITKEN was called.'
!!!        stop 8888
!!!        end
C*ID* BKGPTELP C109.PTOLEMY.LINKULES                        PER704  14:
cccCNI%'LINKUL' = 'BKGPTElp'
      SUBROUTINE BKGPTElp ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,        190   
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC,
     5  IPARAM )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),               200   
     1  FLOAT(100), TEMPVS(6), PARAM(20), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1),
     4  IPARAM(5)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     BALTZ, KAUFFMANN, GLENDENNING AND PRUESS  2+ TELP
C
C     COMPUTES THE IMAGINARY POTENTIAL AS                                 210   
C
C        VI, RI, AI, IMAPOWER   WOODS-SAXON
C      + VSI, RSI, ASI          SURFACE TERM
C      + BKGP TELP FOR A 2+ STATE.
C
C     THE TELP PART IS CONTROLLED BY
C
C     IPARAM5 - 1 = PROJECTILE EXCITATION
C               2 = TARGET EXCITATION
C               IF NOT DEFINED, THEN  JA @= JB  OR  JBIGA @= JBIGB        220   
C          IS USED TO DETERMIN WHICH NUCLEUS IS EXCITED.
C     PARAM20 - CONTAINS THE G2(ETA,XI) FUDGE FACTOR.  IF NOT DEFINED,
C               THE VALUE 1 IS USED.
C     BETACOULOMB
C
C
C     WE STORE SOME CONSTANTS IN THE BEGINNING OF THE WORK AREA
C     ALLOC(LOC(IWRK)+I):
C     I   STORED
C     0   R2FAC   R2 = R2FAC*RC                                           230   
C     1   BETACOULOMB
C     2   ETA
C     3   WPP = (3/8) WP
C     4   AC  = A PARAMETER FOR L = L(CRITICAL)
C     5   BC
C     6   CC
C     7   N1 = BOUNDARY AT RC.
C     8   AK (K-VALUE)
C     9   RC
C                                                                         240   
C     3/13/79 - MADE FROM LTSTELP - S.P.
C     7/30/79 - IPARAM ADDED TO ARGUMENT LIST - S.P.
C     1/16/80 - STANDARD MORTRAN - RPG
C
C
C
      CHARACTER*4  ID
      CHARACTER*4  WORKNM(2,2) / '*000', 'SHAP', '*000', 'WORK' /
      CHARACTER*8  WORDS(2) / 'PROJCTIL', '  TARGET' /
C                                                                         250   
C
C     DEFINITIONS OF THE COEFFICENTS IN THE 1/R EXPANSION
C
      ATERM(DL) = ETA*AK**2/DL**2 *
     1  ( ETA*(3*DL**2+ETASQ)/(DL**2+ETASQ)**2
     2    - DATAN(DL/ETA)/DL )
      BTERM(DL) = 4*ETA*AK*DL**2/(DL**2+ETASQ)**2
      CTERM(DL) = 2*DL**4/(DL**2+ETASQ)**2
C
C                                                                         260   
      NUMOUT = 0
      IRETUR = 0
C
      IF ( IREQUE .EQ. 4 )  GO TO 700
C
      UNDEF = FINTRN(1)
      NOTDEF = INTRNL(3)
C
      VI = TEMPVS(2)
      AI = TEMPVS(6)                                                      270   
      RI = FLOAT(47)
      POWER = FLOAT(103)
      VSI = FLOAT(71)
      RSI = FLOAT(69)
      ASI = FLOAT(68)
      G2 = PARAM(20)
      IF ( G2 .EQ. UNDEF )  G2 = 1
      ICHANW = INTRNL(5)
      AK  = FKANDM(ICHANW)
      RC = FLOAT(45)                                                      280   
      ZP = INTGER(24)
      ZT = INTGER(25)
      ETA = FKANDM(ICHANW+16)
      ALPHA = 1/CNSTNT(8)
      PI = CNSTNT(1)
      HBARC = CNSTNT(6)
      E = FLOAT(12)
      AM = FLOAT(24)
C
      IF ( IREQUE .GT. 1 )  GO TO 200                                     290   
C
C     IS THIS PROJECTILE OR TARGET EXCITATION.
C     WE USE THIS INFORMATION TO CONSTRUCT R2FAC WHICH MAKES THE
C     SINGLE NUCLEON RADIUS OUT OF THE SUM (RC).
C
      IPORT = IPARAM(5)
      IF ( IPORT .NE. NOTDEF )  GO TO 50
      IF ( JBLOCK(2) .NE. JBLOCK(3) )  IPORT = 1
      IF ( JBLOCK(4) .NE. JBLOCK(5) )  IPORT = 2
 50   IF ( IPORT .NE. 1  .AND.  IPORT .NE. 2 )  GO TO 930                 300   
C      AMP = FLOAT(40)
C      AMT = FLOAT(41)
      R0MASS = FINTRN(10)
      R2FAC = FLOAT(39+IPORT)**(1.D0/3.D0) / R0MASS
C
      IB = NAMLOC('BETACOUL')
      IF ( IB .EQ. 0 )  GO TO 940
C
      IF ( VI .EQ. 0 )  GO TO 65
      IF ( RI .EQ. UNDEF  .OR.  AI .EQ. UNDEF )  GO TO 950                310   
  65  IF ( VSI .EQ. 0 )  GO TO 70
      IF ( RSI .EQ. UNDEF  .OR.  ASI .EQ. UNDEF )  GO TO 945
C
C     GENERATE WORK ARRAY TO SAVE THE W.S.
C
  70  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
      IWRK = NALLOC( NUMPTS+10, WORKNM(1,1) )
      LWRK = LOC(IWRK)                                                    320   
      ALLOC(LWRK) = R2FAC
      ALLOC(LWRK+1) = ALLOC(LOC(IB))
      GO TO 300
C
C     COMPUTE THE CONSTANTS
C     THIS PART IS DONE FOR IREQUE = 2 (PRINTING) AND 3 (NEW POTENTIAL)
C
 200  IWRK = MYINTS(1)
      LWRK = LOC(IWRK)
      IF ( VI .NE. 0  .AND.  ( AI .LE. 0  .OR.  RI .LE. 0 ) )             330   
     1  GO TO 960
      IF ( VSI .NE. 0  .AND.  ( ASI .LE. 0  .OR.  RSI .LE. 0 ) )
     1  GO TO 960
      RC2 = RC*ALLOC(LWRK)
      BETA = ALLOC(LWRK+1)
      RD = ZP*ZT*ALPHA*HBARC/E
      WP = 3 * (ZP*ZT*ALPHA*BETA)**2 * AM * G2
     1  * RC2**4
     2  / (50*PI*AK)
      WPP = 3*WP/8.D0                                                     340   
      ALLOC(LWRK+2) = ETA
      ALLOC(LWRK+3) = WPP
      ALLOC(LWRK+8) = AK
      ETASQ = ETA**2
C
C     ESTIMATE CRITICAL L (IN VERY SIMPLE FASHION) AND GET CONSTANT
C     (IN L) PARAT OF THE POTENTIAL PARAMETERS.
C
      DLCSQ = AK*RC*(AK*RC - 2*ETA)
      IF ( DLCSQ .LT. .25 )  DLCSQ = .25                                  350   
      DLC = DSQRT(DLCSQ)
      AC = ATERM(DLC)
      BC = BTERM(DLC)
      CC = CTERM(DLC)
      ALLOC(LWRK+4) = AC
      ALLOC(LWRK+5) = BC
      ALLOC(LWRK+6) = CC
      FATRC = ( (CC/RC + BC)/RC + AC )/RC**3
C
C                                                                         360   
 300  GO TO ( 400, 500, 600, 700 ), IREQUE
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = IWRK
      MYINTS(2) = II
C
C     HERE WE INDICATE THAT THIS IS A L-DEPENDENT POTENTIAL
C
      IRETUR = +1                                                         370   
C
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT
C
 500  IF ( VI .EQ. 0 )  GO TO 510
      R0 = RI/R0MASS
      WRITE (IUNIT, 503)  VI, RI, AI, R0, POWER
 503  FORMAT ( ' W.S. WELL:', T20, G14.3, T37, F7.4,                      380   
     1  T49, F7.4, T67, F7.4, T80, F7.3 )
 510  IF ( VSI .EQ. 0 )  GO TO 520
      R0 = RSI/R0MASS
      WRITE (IUNIT, 513)  VSI, RSI, ASI, R0
 513  FORMAT ( ' SURFACE ABSORPTION:', T20, G14.3, T37, F7.4,
     1  T49, F7.4, T67, F7.4, T80, F7.3 )
C
 520  WPRC5 = WP/RC**5
      WRITE (IUNIT, 523)  WORDS(IPORT), BETA, G2, RC, RD, WPRC5
 523  FORMAT ( 1X, A8, ' 2+ BKGP TELP WITH BETACOULOMB =', F7.4,          390   
     1  5X, 'G2 =', G12.5,  5X, 'RC =', F7.4,
     2  5X, 'RD =', F7.4,  5X, 'WP/RC**5 =', G12.5, ' MEV' )
      FFC = FATRC*RC**5
      AAC = AC*RC**2
      BBC = BC*RC
      WRITE (IUNIT, 527)  DLC, FFC, AAC, BBC, CC
 527  FORMAT ( ' SIMPLE L-HAT CRITICAL =', F8.1,
     1  5X, 'F(LCRIT, RC)*RC**5 =', G13.5,
     2  5X, 'TERMS =', 23G14.5 )
      GO TO 900                                                           400   
C
C     SET UP THE TWO PIECES
C
 600  LL = LWRK+10
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   1, VI, RI, AI, POWER )
      DO 619  I = 1, NUMPTS
         ARRAY1(I) = ALLOC(LL-1+I)
 619  CONTINUE
C                                                                         410   
      IF ( VSI .EQ. 0 )  GO TO 650
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   3, VSI, RSI, ASI, 1.D0 )
      DO 629  I = N1, N2
         ARRAY1(I) = ARRAY1(I) + ALLOC(LL-1+I)
 629  CONTINUE
C
C     NOW THE BKGP TELP FOR L = LCRIT
C
 650  N1 = (RC-RSTART)/STEPSZ + 1.5                                       420   
      N1 = MIN0(N1, NUMPTS)
      ALLOC(LWRK+7) = N1
      ALLOC(LWRK+9) = RC
      DO 679  I = 1, N1
         ARRAY1(I) = ARRAY1(I) - WPP*FATRC
 679  CONTINUE
C
      N1 = N1+1
      IF ( N1 .GT. NUMPTS )  RETURN
      R = RSTART + (N1-1)*STEPSZ                                          430   
      DO 689  I = N1, NUMPTS
         RINV = 1/R
         ARRAY1(I) = ARRAY1(I)
     1      - WPP * ( (CC*RINV +BC)*RINV + AC )*RINV**3
         R = R+STEPSZ
 689  CONTINUE
      RETURN
C
C
C     ADD IN THE L-DEPENDANT PART                                         440   
C
C     NOTE THAT WE ARE REQUIRED TO ADD IN TO ARRAY1 THE QUANTITY
C
C        - H**2/(12*E) * DELTA V(L)
C
C     WHERE  DELTA V(L)  HAS ITS TRUE SIGN (NEGATIVE FOR MORE ATTRACTION
C     -H**2/(12*E) IS TO BE FOUND IN ARRAY2(1)
C
C
 700  LWRK = LOC(MYINTS(1))                                               450   
      LL = LWRK+10
      ETA = ALLOC(LWRK+2)
      WPP = ALLOC(LWRK+3) * ARRAY2(1)
      AC = ALLOC(LWRK+4)
      BC = ALLOC(LWRK+5)
      CC = ALLOC(LWRK+6)
      N1 = ALLOC(LWRK+7)
      AK = ALLOC(LWRK+8)
      RC = ALLOC(LWRK+9)
      ETASQ = ETA**2                                                      460   
      DL = L+.50D0
      A = ATERM(DL)
      B = BTERM(DL)
      C = CTERM(DL)
      FATRC = ( ( (C-CC)/RC + (B-BC) )/RC + (A-AC) ) / RC**3
      DO 779  I = 1, N1
         ARRAY1(I) = ARRAY1(I) - WPP*FATRC
 779  CONTINUE
      R = RSTART + N1*STEPSZ
      N1 = N1+1                                                           470   
      IF ( N1 .GT. NUMPTS )  RETURN
      DO 789  I = N1, NUMPTS
         RINV = 1/R
         ARRAY1(I) = ARRAY1(I)
     1      - WPP * ( ( (C-CC)*RINV + (B-BC) )*RINV + (A-AC) ) * RINV**3
         R = R + STEPSZ
 789  CONTINUE
      RETURN
 900  NUMOUT = 1
cib      REWIND IUNIT                                                     480   
      RETURN
C
 930  WRITE (IUNIT, 933) IPORT
 933  FORMAT ( '0**** IPARAM5 =', I3, ' MUST BE 1 OR 2 TO INDICATE',
     1  ' PROJECTILE OR TARGET EXCITATION FOR BKGPTELP.' )
      GO TO 955
 940  WRITE (IUNIT, 943)
 943  FORMAT ( '0**** BETACOULOMB MUST BE DEFINED FOR BKGPTELP.' )
      GO TO 955
 945  WRITE (IUNIT, 947)                                                  490   
 947  FORMAT ( '0**** RSI (OR RSI0) AND ASI MUST BE DEFINED.' )
      GO TO 955
 950  WRITE (IUNIT, 953)
 953  FORMAT ('0**** RI (OR RI0) AND AI MUST BE DEFINED.' )
 955  IRETUR = -1
      GO TO 900
C
 960  IRETUR = -1
      RETURN
C                                                                         500   
      END
C*ID* DEFORMED C109.PTOLEMY.LINKULES                        PER704  13:
cccCNI%'LINKUL' = 'DEFORMEd'
      SUBROUTINE DEFORMEd ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC,
     5  IPARAM )                                                          510   
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1),
     4  IPARAM(5)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4                  520   
C
C     DEFORMED - DEFORMED WOODS-SAXON OPTICAL POTENTIAL
C
C     THE POTENTIAL IS DEFINED BY  V, R, A, AND POWER
C     THE DEFORMATION IS DEFINED BY THE ARRAY BETAP OR BETAT, WHICH
C     CONTAIN BETAS FOR  LXMIN, LXMIN+1, LXMIN+2, ...
C     IPARAM4 = FIRST LX (LXMIN) FOR WHICH THE BETA ARRAYS ARE DEFINED.
C               THE DEFAULT IS 2. IF BOTH BETAP AND BETAT ARE
C               DEFINED, BOTH MUST START FROM THE SAME LXMIN.
C     IF EXCITATION IS OCCURING (I.E. NOT OPTICAL POTENTIAL)              530   
C     THEN JBLOCK(1) IS USED TO DETERMIN WHICH NUCLEUS IS BEING
C     EXCITED.  JJ IS MOD(JBLOCK(1), 100) - 3 AND HAS THE MEANINGS
C          JJ = 1 FOR PROJECTILE DEFORMATION
C               2 FOR TARGET
C               3 FOR MUTUAL DEFORMATION
C     FURTHERMORE, FOR MUTUAL EXCITATION THE INDIVIDUAL MULTIPOLES
C     ARE FOUND AS
C          LXP = MOD( JBLOCK(1)/100, 100 )   AND
C          LXT = MOD( JBLOCK(1)/10000, 100 ) ,
C     IN THIS CASE LX IS IGNORED                                          540   
C     NCOSINE = NUMBER OF GAUSS PTS ON (0,1) (10=DEFAULT)
C     LX IS USED AS THE DESIRED PROJECTION
C     LX = -1 MEANS THAT THE DEFORMED OPTICAL POT IS TO BE FOUND
C
C     7/27/81 - BASED ON FIXEDWOO - S.P.
C     9/22/81 - SET LXOUT=-1 TO 0 - S.P.
C     10/16/81 - ALLOW BOTH NUCLEI TO BE DEFORMED - S.P.
C     5/14/82 - CONSIDER LX=NOTDEF TO BE LX=-1 - S.P.
C     12/19/84 - ADD DEFORMED SURFACE POT - S.P.
C     8/19/85 - REWIND IUNIT ONLY ON IBM - S.P.                           550   
C
C
      DIMENSION  BETAS(2,20), RBETAS(2,20), IBETA(2), LDEFMX(2),
     1   RDEFS(2), LXS(2)
      CHARACTER*8  DEFNAM(3) / 'PROJECTI', 'TARGET', 'SUMMED' /
C
      COMMON /CNSTNT/  PI, RT4PI, PIINV, RADIAN, DEGREE,                /cnstnt/
     &  HBARC, AMUMEV, AFINE, BIGEST, SMLEST, BIGNUM, SMLNUM, BIGLOG
 
C
      NUMOUT = 0
      IRETUR = 0                                                          560   
      UNDEF = FINTRN(1)
      NOTDEF = INTRNL(3)
      PI = CNSTNT(1)
      BIGLOG = CNSTNT(13)
C
      IF ( IPOTYP .EQ. 13 )  GO TO 20
      V = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)                                 570   
      POWER = FLOAT(101+IPOTYP)
      GO TO 30
C
C     GET PARAMETERS FOR DEFORMED SURFACE POTENTIAL
C
 20   V = FLOAT(71)
      R = FLOAT(69)
      A = FLOAT(68)
      POWER = 1
C                                                                         580   
 30   VV = 0
      IF ( V .EQ. 0 )  GO TO 100
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  V .EQ. UNDEF
     1   )  GO TO 950
      IF ( R .LE. 0  .OR.  A .LE. 0 )  GO TO 960
      LXOUT = INTGER(5)
      IF ( LXOUT .EQ. -1  .OR.  LXOUT .EQ. NOTDEF )  LXOUT = 0
      LXS(1) = 0
      LXS(2) = 0
      LXMIN = IPARAM(4)                                                   590   
      IF (LXMIN .EQ. NOTDEF )  LXMIN = 2
      IBETA(1) = MYINTS(1)
      IBETA(2) = MYINTS(2)
      IF ( IREQUE .NE. 1 )  GO TO 40
C
C     FIND BETAP, BETAT
C
      IBETA(1) = NAMLOC( 'BETAP   ' )
      IBETA(2) = NAMLOC( 'BETAT   ' )
      IF ( IBETA(1)+IBETA(2) .EQ. 0 )  GO TO 970                          600   
C
C     GET THE BETAS AND COMPUTE R*BETA
C
  40  R0MASS = FINTRN(10)
      DO 49  I = 1, 2
         LDEFMX(I) = 0
         IF ( IBETA(I) .EQ. 0 )  GO TO 49
         LDEFMX(I) = LXMIN + LENGTH(IBETA(I))-1
         LXMAX = LDEFMX(I)
         RDEF = R*FLOAT(39+I)**(1.D0/3.D0) / R0MASS                       610   
         RDEFS(I) = RDEF
         LBETA = LOC(IBETA(I))
         DO 44  LX = LXMIN, LXMAX
            BETA = ALLOC(LBETA)
            BETAS(I,LX) = BETA
            RBETAS(I,LX) = RDEF*BETA
            LBETA = LBETA+1
 44      CONTINUE
 49   CONTINUE
      NCOSIN = INTGER(11)                                                 620   
      IF ( NCOSIN .EQ . NOTDEF )  NCOSIN = 10
      IF ( NCOSIN .GT. 100 )  NCOSIN = 10
      VV = V
      IPORT = 0
      IF ( INTGER(5) .EQ. -1  .OR.  INTGER(5) .EQ. NOTDEF )  GO TO 100
C
C     DETERMIN WHAT PROJECTION IS NEEDED FOR EXCITATION
C
      IF ( IREQUE .NE. 3 )  GO TO 100
      IPORT = MOD( JBLOCK(1), 100 ) - 2                                   630   
      IF ( IPORT .GE. 3 )  GO TO 80
      IF ( IPORT .LT. 1  .OR.  IBETA(IPORT) .EQ. 0 )  GO TO 965
      LXS(IPORT) = LXOUT
      VV = -V/RDEFS(IPORT)
      GO TO 100
  80  LXS(1) = MOD( JBLOCK(1)/100, 100 )
      LXS(2) = MOD( JBLOCK(1)/10000, 100 )
      VV = -V / ( RDEFS(1)*RDEFS(2) )
C
C                                                                         640   
 100  GO TO ( 200, 500, 600 ), IREQUE
 200  MYINTS(1) = IBETA(1)
      MYINTS(2) = IBETA(2)
      RETURN
C
C     PRINT IT OUT
C
 500  IF ( V .EQ. 0 )  GO TO 550
      IF ( IPOTYP .NE. 13 )  WRITE (IUNIT, 503)
 503  FORMAT ( ' THE POTENTIAL IS A DEFORMED WOODS-SAXON:' )              650   
      IF ( IPOTYP .EQ. 13 )  WRITE (IUNIT, 504)
 504  FORMAT ( ' THE POTENTIAL IS A DEFORMED DERIVATIVE OF A ',
     1  'WOODS-SAXON:' )
      WRITE (IUNIT, 508)  V, R, A, POWER
 508  FORMAT (
     1  ' V =', F7.2, '   R =', F7.4, '   A = ', F7.4,
     1  5X, 'POWER =', F8.4 )
      DO 519  I = 1, 2
         IF ( IBETA(I) .EQ. 0 )  GO TO 519
         LXMAX = LDEFMX(I)                                                660   
      WRITE (IUNIT, 513)   DEFNAM(I), RDEFS(I),
     1  ( LX, BETAS(I,LX), RBETAS(I,LX), LX = LXMIN, LXMAX )
  513 FORMAT ( '0THE ', A8, ' RADIUS OF', F8.4,
     2     ' FM. IS BEING DEFORMED.' /
     3  ' LX   BETA   DEFORMATION LENGTH' /
     4  ( I3, F8.4, F15.4 ) )
  519 CONTINUE
      WRITE ( IUNIT, 523 )  NCOSIN
 523  FORMAT ( I5, '-POINT GAUSS INTEGRATION IS BEING USED.' )
      GO TO 900                                                           670   
 550  WRITE ( IUNIT, 553 )
 553  FORMAT ( ' THE POTENTIAL IS IDENTICALLY ZERO.' )
      GO TO 900
C
C     DO IT
C
 600  ITYPE = 1
      IF ( IPOTYP .EQ. 13 )  ITYPE = 3
      CALL DEFWOO ( NUMPTS, RSTART, STEPSZ, ARRAY1, N1, N2,
     1   ITYPE, VV, R, A, POWER, IPORT, LXMIN, LDEFMX,                    680   
     2   RBETAS, LXS, NCOSIN )
C
      IF ( VV .EQ. 0 )  GO TO 900
      VV = VV/V
      IF ( IPORT .GT. 0 )  GO TO 610
      WRITE ( IUNIT, 603 )  VV, LDEFMX
 603  FORMAT ( ' DEFORMED OPTICAL POTENTIAL WEIGHTED BY', G15.5,
     1  ' WAS COMPUTED WITH LXP AND LXT MAXIMUM =', 2I3 )
      GO TO 680
 610  IF ( IPORT .GE. 3 )  GO TO 620                                      690   
      WRITE ( IUNIT, 613 )  LXOUT, VV, DEFNAM(IPORT)
 613  FORMAT ( ' THE LAMBDA =', I2,
     1 ' PROJECTION OF THE DEFORMED POTENTIAL WEIGHTED BY',
     2 G15.5, ' WAS COMPUTED WITH A DEFORMED ', A8,
     3  ' RADIUS.' )
      GO TO 680
 620  WRITE ( IUNIT, 623 )  VV, LXS
 623  FORMAT ( ' THE MUTUAL DEFORMED POTENTIAL WEIGHTED BY',
     1  G15.5, ' WAS COMPUTED FOR LXP, LXT =', 2I3 )
C                                                                         700   
C     INDICATE THAT NO DERIVATIVE IS NEEDED IN THE FORM FACTOR
C
  680 IF ( IPORT .GT. 0 )  ISWTCH(30) = 1
C
C
 900  NUMOUT = 1
CIB      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)                                                  710   
 953  FORMAT ( '0**** R, V, AND A MUST BE DEFINED.' )
C
 955  IRETUR = -1
      GO TO 900
C
 960  WRITE ( IUNIT, 963 )  R, A
 963  FORMAT ( '0*** R OR A IS INVALID:', 2G20.10 )
      GO TO 955
C
 965  WRITE ( IUNIT, 968 ) JBLOCK(1), IPORT, IREQUE, LXS, IBETA           720   
 968  FORMAT ( '0*** WRONG NUCLEUS DEFORMED:', I15, 8I8 )
      GO TO 955
C
 970  WRITE ( IUNIT, 973 )
 973  FORMAT ( '0*** BETAP OR BETAT MUST BE DEFINED' )
      GO TO 955
      END
C*ID* FIXEDWOO C109.PTOLEMY.LINKULES                        PER704  14:
CDC%'LINKUL' = 'FIXEDWOo'
      SUBROUTINE FIXEDWOo ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,        730   
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),          740   
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     FIXEDWOODS - WOODS SAXON FIXED AT A POINT
C
C     THE POTENTIAL IS DEFINED BY  R, A, AND POWER AND HAS DEPTH
C     V AT RADIUS PARAM1
C                                                                         750   
C     4/6/78 - BASED ON SHAPE - S.P.
C     7/7/78 - CDC VERSION
C
C
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /
C
      NUMOUT = 0
      IRETUR = 0                                                          760   
      UNDEF = FINTRN(1)
C
      RFIX = PARAM(IPOTYP)
      VFIX = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)
      POWER = FLOAT(101+IPOTYP)
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  VFIX .EQ. UNDEF
     1  .OR.  RFIX .EQ. UNDEF )  GO TO 950                                770   
      IF ( R .LE. 0  .OR.  A .LE. 0  .OR.  RFIX .LE. 0 )  GO TO 960
      V = VFIX * ( 1 + DEXP((RFIX-R)/A) )**POWER
C
C
      GO TO ( 100, 500, 600 ), IREQUE
C
C     GENERATE RADIUS ARRAY FOR INTRPC
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID                                                    780   
C
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT                                                        790   
C
 500  WRITE (IUNIT, 503)  V, R, A, POWER, VFIX, RFIX
 503  FORMAT ( ' THE POTENTIAL IS A WOODS-SAXON:', 3G15.5 ,
     1  5X, 'POWER =', G14.5 /
     1  ' WITH V =', G15.5, 5X, 'AT R =', G15.5 )
      GO TO 900
C
C     DO IT
C
 600  CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ARRAY1, N1, N2,               800   
     1   1, V, R, A, POWER )
C
      RETURN
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** R, V, A, AND PARAM1 MUST BE DEFINED.' )             810   
      IRETUR = -1
      GO TO 900
C
 960  IRETUR = -1
      RETURN
C
      END
C*ID* GAUSSIAN C109.PTOLEMY.LINKULES                        PER704  14:
CDC%'LINKUL' = 'GAUSSIAn'
      SUBROUTINE GAUSSIAn ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,        820   
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),          830   
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     GAUSSIAN - GENERALIZED GAUSSIAN POTENTIAL
C
C     THE POTENTIAL IS DEFINED BY  V, R, A, AND REALPOWER,  OR BY
C     W, RI, AI, AND IMAGPOWER.  THE FORM IS
C                                                                         840   
C     POT(X) = - V EXP( ( R**(2P) - X**(2P) )/A**(2P) )
C
C     WHERE P IS ____POWER.  THUS POT(X=R) IS V.  NOTE THAT V AND R
C     ARE NOT INDEPENDANT PARAMETERS SINCE THE SAME POTENTIAL OBTAINS
C     FOR  V EXP((R/A)**(2P))  A CONSTANT.  THE PARAMETER R IS PROVIDED
C     TO AVOID OVERFLOW BY ALLOWING THE DEFINITION OF THE DEPTH
C     OF THE POTENTIAL AT A PHYSICALLY MEANINGFUL POINT.
C
C     5/30/78 - MADED FROM FIXEDWOOD - S.P.
C     7/7/78 - CDC VERSION                                                850   
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /
C
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
C
      V = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)                                                860   
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)
      POWER = FLOAT(101+IPOTYP)
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  V .EQ. UNDEF
     1  .OR.  POWER .EQ. UNDEF )  GO TO 950
      IF ( R .LE. 0  .OR.  A .LE. 0  .OR.  POWER .LE. 0 )  GO TO 960
C
C
      GO TO ( 100, 500, 600 ), IREQUE
C                                                                         870   
C     GENERATE RADIUS ARRAY FOR INTRPC
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II                                                      880   
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  V, R, A, POWER
 503  FORMAT ( ' THE POTENTIAL IS A GENERALIZED GAUSSIAN:' /
     1  ' V =', G13.5, ' AT R =', F6.2, ' FM;   A =', F8.3,
     2  ' FM,   POWER =', F8.4 )
      GO TO 900                                                           890   
C
C     DO IT
C
 600  TWOP = POWER+POWER
C
C     IN THE FOLLOWING  46 IS LN(1E+20)
C     AND 23 IS LN(1E+10)
C
      NMAX = ( R**TWOP + A**TWOP*(46+DLOG(DABS(V))) )**(1/TWOP) / STEPSZ
      NMAX = MIN0( NUMPTS, NMAX )                                         900   
      NMIN = 1
      TEMP = R**TWOP - A**TWOP*(23-DLOG(DABS(V)))
      IF ( TEMP .GT. 0 )  NMIN = TEMP**(1/TWOP) / STEPSZ
      NMIN = MAX0( NMIN, 1 )
C
      IF ( NMIN .EQ. 1 )  GO TO 650
      DO 639  I = 1, NMIN
         ARRAY1(I) = -1.D+10
 639  CONTINUE
C                                                                         910   
 650  RVAL = (NMIN-1)*STEPSZ/A
      DO 679  I = NMIN, NMAX
         ARRAY1(I) = - V * DEXP( (R/A)**TWOP - RVAL**TWOP )
         RVAL = RVAL + (STEPSZ/A)
 679  CONTINUE
C
      IF ( NMAX .GE. NUMPTS )  RETURN
      NMAX = NMAX + 1
      DO 699  I = NMAX, NUMPTS
         ARRAY1(I) = 0                                                    920   
 699  CONTINUE
C
      RETURN
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** R, V, AND A MUST BE DEFINED.' )                     930   
      IRETUR = -1
      GO TO 900
C
 960  IRETUR = -1
      RETURN
C
      END
C*ID* JDEPEN   C109.PTOLEMY.LINKULES                        PER704  17:
ccc%'LINKUL' = 'JDEPEN'
      SUBROUTINE JDEPEN ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR, IUNIT,   940   
     1  NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),          950   
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     JDEPEN - J-DEPENDANT DEPTH FOR POTENTIAL
C
C     THE POTENTIAL IS A WOODS SAXON WITH J DEPENDANT DEPTH:
C
C       V = (1+PARAM1*J+PARAM2*J**2) * WOODS SAXON( V R A )               960   
C
C     6/13/78 - BASED ON FIXEDWOOD - S.P.
C     7/11/78 - CDC VERSION AND IRETUR CORRECTION
C     11/16/81 - JDEPEN BASED ON PARITWOOD
C
C
      COMMON /CNSTNT/  PI, RT4PI, PIINV, RADIAN, DEGREE,                /cnstnt/
     &  HBARC, AMUMEV, AFINE, BIGEST, SMLEST, BIGNUM, SMLNUM, BIGLOG
 
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /         970   
C
C
C
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
      BIGLOG = CNSTNT(13)
C
      PARAM1 = PARAM(1)
      PARAM2 = PARAM(2)                                                   980   
      IF ( PARAM2 .EQ. UNDEF )  PARAM2 = 0
      V = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)
      POWER = FLOAT(101+IPOTYP)
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  V .EQ. UNDEF
     1  .OR.  PARAM1 .EQ. UNDEF )  GO TO 950
      IF ( R .LE. 0  .OR.  A .LE. 0  )  GO TO 960
C                                                                         990   
C
      GO TO ( 100, 500, 600, 700 ), IREQUE
C
C     GENERATE WORK ARRAY TO SAVE THE W.S.
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
      I = NALLOC( NUMPTS, WORKNM(1,1) )
C                                                                        1000   
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
C
C     HERE WE INDICATE THAT THIS IS A L-DEPENDENT POTENTIAL
C
      IRETUR = +1
C
      IF ( NUMOUT .NE. 0 )  GO TO 900                                    1010   
      RETURN
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  V, R, A, POWER, PARAM1, PARAM2
 503  FORMAT ( ' THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:', 3G15.5 ,
     1  5X, 'POWER =', G14.5 /
     1  ' WITH THE J-DEPENDANT DEPTH FACTOR',
     2  '  1 + J *', G13.5, ' + J**2 *', G13.5 )
      GO TO 900                                                          1020   
C
C     SET UP THE TWO PIECES
C
 600  LL = LOC(MYINTS(1))
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   1, V, R, A, POWER )
      DO 619  I = 1, NUMPTS
         ARRAY1(I) = ALLOC(LL-1+I)
 619  CONTINUE
C                                                                        1030   
      RETURN
C
C
C     ADD IN THE J-DEPENDANT PART
C
 700  LL = LOC(MYINTS(1))
      DJ = .5*JBLOCK(1)
      FAC = ARRAY2(1)*DJ*(PARAM1 + DJ*PARAM2)
      N1 = RSTART/STEPSZ + .5
      LL = LL - 1 + N1                                                   1040   
      DO 719  I = 1, NUMPTS
         ARRAY1(I) = ARRAY1(I) + FAC*ALLOC(LL+I)
 719  CONTINUE
CCC      WRITE ( IUNIT, 723 )  DJ, FAC, RSTART, N1, NUMPTS
CCC  723 FORMAT ( ' J-DEPENDANT DEPTH:  J, FAC =', 2G15.8,
CCC     1  5X, 'RSTART, N1 =', G15.8, 2I6 )
C
      RETURN
C
C                                                                        1050   
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** R, V, A, AND PARAM1 MUST BE DEFINED.' )
      IRETUR = -1
      GO TO 900
C
 960  IRETUR = -1                                                        1060   
      RETURN
C
      END
C*ID* JDEPENWS C109.PTOLEMY.LINKULES                        PER704  12:
cccCDC%'LINKUL' = 'JDEPENWS'
      SUBROUTINE JDEPENWS ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,                                          1070   
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4                 1080   
C
C     JDEPENWS - J-DEPENDANT W.S. DEPTH FOR POTENTIAL
C
C     THE POTENTIAL IS A WOODS SAXON WITH J-DEPENDANT W.S. DEPTH:
C
C       V = F(J) * WOODS SAXON( V R A )
C       F(J) = 1 / ( 1 + EXP( (J-PARM100)/PARM101 ) )
C
C     6/13/78 - BASED ON FIXEDWOOD - S.P.
C     7/11/78 - CDC VERSION AND IRETUR CORRECTION                        1090   
C     11/16/81 - JDEPEN BASED ON PARITWOOD
C     5/11/82 - JDEPENWS MADE FROM JDEPEN - S.P.
C
C
      COMMON /CNSTNT/  PI, RT4PI, PIINV, RADIAN, DEGREE,                /cnstnt/
     &  HBARC, AMUMEV, AFINE, BIGEST, SMLEST, BIGNUM, SMLNUM, BIGLOG
 
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /
C
C                                                                        1100   
C
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
      BIGLOG = CNSTNT(13)
C
      PARM10 = PARAM(10)
      PARM11 = PARAM(11)
      V = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)                                               1110   
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)
      POWER = FLOAT(101+IPOTYP)
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  V .EQ. UNDEF
     1  .OR.  PARM10 .EQ. UNDEF  .OR.  PARM11 .EQ. UNDEF )  GO TO 950
      IF ( R .LE. 0  .OR.  A .LE. 0  )  GO TO 960
C
C
      GO TO ( 100, 500, 600, 700 ), IREQUE
C                                                                        1120   
C     GENERATE WORK ARRAY TO SAVE THE W.S.
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
      I = NALLOC( NUMPTS, WORKNM(1,1) )
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I                                                      1130   
      MYINTS(2) = II
C
C     HERE WE INDICATE THAT THIS IS A L-DEPENDENT POTENTIAL
C
      IRETUR = +1
C
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT                                                       1140   
C
 500  WRITE (IUNIT, 503)  V, R, A, POWER, PARM10, PARM11
 503  FORMAT ( ' THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:', 3G15.5 ,
     1  5X, 'POWER =', G14.5 /
     1  ' WITH THE J-DEPENDANT DEPTH FACTOR',
     2  '  1 / ( 1 + EXP( ( J - ', G12.5, ') / ', G12.5, ') )' )
      GO TO 900
C
C     SET UP THE TWO PIECES
C                                                                        1150   
 600  LL = LOC(MYINTS(1))
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   1, V, R, A, POWER )
      DO 619  I = 1, NUMPTS
         ARRAY1(I) = ALLOC(LL-1+I)
 619  CONTINUE
C
      RETURN
C
C                                                                        1160   
C     ADD IN THE J-DEPENDANT PART
C
 700  LL = LOC(MYINTS(1))
      DJ = .5*JBLOCK(1)
      FAC = 1
      X = (DJ-PARM10)/PARM11
      IF ( X .LT. -30 )  GO TO 710
      FAC = 0
      IF ( X .GT. 30 )  GO TO 710
      FAC = 1/(1+DEXP(X))                                                1170   
 710  FAC = ARRAY2(1) * ( FAC - 1 )
      N1 = RSTART/STEPSZ + .5
      LL = LL - 1 + N1
      DO 719  I = 1, NUMPTS
         ARRAY1(I) = ARRAY1(I) + FAC*ALLOC(LL+I)
 719  CONTINUE
CCC      WRITE ( IUNIT, 723 )  DJ, FAC, RSTART, N1, NUMPTS
CCC  723 FORMAT ( ' J-DEPENDANT DEPTH:  J, FAC =', 2G15.8,
CCC     1  5X, 'RSTART, N1 =', G15.8, 2I6 )
C                                                                        1180   
      RETURN
C
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** R, V, A, PARAM10 & PARAM11 MUST BE DEFINED.' )
      IRETUR = -1                                                        1190   
      GO TO 900
C
 960  IRETUR = -1
      RETURN
C
      END
C*ID* LAGRANGE C109.PTOLEMY.LINKULES                        PER704  14:
CDC%'LINKUL' = 'LAGRANGe'
      SUBROUTINE LAGRANGE ( ALIAS, MYINTS, IPOTYP, IREQUE,
     1  IRETUR, IUNIT, NUMOUT,                                           1200   
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(2,10), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),                   1210   
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL  ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     LAGRANGE - CONSTRUCT POTENTIAL BY LAGRANGE INTERPOLATION
C
C     THIS LINKULE ASSUMES THAT THE PAIRS OF NUMBERS:
C       (PARAM1, PARAM2), (PARAM3, PARAM4), ..., (PARAM19, PARAM20)
C     ARE TO BE INTERPERTED AS
C       (R1, V1), (R2, V2), ..., (R10, V10)                              1220   
C     AND DEFINE A POTENTIAL V(R).  FEWER THAN 10 POINTS MAY
C     BE USED BY LEAVING THE REMAINING V'S OR R'S UNDEFINED.
C     THE R'S DO NOT HAVE TO BE ORDERED.
C     THE POTENTIAL IS DEFINED AS (NOTE SIGN CHANGE)
C
C       V(R) =  V1    R < R1 ,
C               V'(R)   R1 < R < RLAST
C               -A*EXP(-B*R) R > RLAST
C            WHERE A AND B ARE FIT TO THE LAST POINTS.
C                                                                        1230   
C     WHERE V'(R) IS THE RESULT OF LAGRANGE INTERPOLATION OR ORDER IPARA
C
C     IF THE POTENTIAL DOES NOT CHANGE SIGN, THEN THE LAGRANGE INTERPOLA
C     IS APPLIED TO LOG|V|.  IF IT DOES CHANGE SIGN, THEN LET N BE
C     THE LARGEST INDEX SUCH THAT  SIGN(V(N-1)) ^= SIGN(V(N)).
C     LAGRANGE INTERPOLATION IS APPLIED TO V FOR  R < R(N)  AND TO
C     LOG|V| FOR R > R(N).
C
C     4/10/78 - BASED ON SPLINE - S.P.
C     2/28/79 - ALLOW CHANGES OF SIGN - S.P.                             1240   
C     3/14/79 - CDC VERSION - S.P.
C
C
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /
      DIMENSION  R(20), V(20), A(30)
C
C
      NUMOUT = 0                                                         1250   
      IRETUR = 0
      UNDEF = FINTRN(1)
C
C     INTERPOLATION ORDER
C
      IORDER = INTGER(38)
CIB      IF ( IORDER .EQ. INTRNL(3) )  GO TO 960
CDC      IF ( IORDER .EQ. INTRNL(2) )  GO TO 960
C
C     GET THE DEFINED POINTS AND VALUES AND ORDER THEM                   1260   
C
      NV = 1
      IV = 1
  40  R(1) = PARAM(1,IV)
      V(1) = PARAM(2,IV)
      IF ( R(1) .NE. UNDEF  .AND.  V(1) .NE. UNDEF )  GO TO 50
      IV = IV+1
      IF ( IV .LT. 10 )  GO TO 40
      GO TO 950
C                                                                        1270   
 50   IV = IV + 1
      IF ( IV .GT. 10 )  GO TO 90
      IF ( PARAM(1,IV) .EQ. UNDEF  .OR.  PARAM(2,IV) .EQ. UNDEF )
     1  GO TO 50
      DO 59  I = 1, NV
         IF ( PARAM(1,IV) .LT. R(I) )  GO TO 70
 59   CONTINUE
C
C     ADD TO END
C                                                                        1280   
      NV = NV+1
      R(NV) = PARAM(1,IV)
      V(NV) = PARAM(2,IV)
      GO TO 50
C
C     INSERT POINT; FIRST MOVE REST RIGHT ONE
C
 70   DO 74  II = I, NV
         R(NV+I-II+1) = R(NV+I-II)
         V(NV+I-II+1) = V(NV+I-II)                                       1290   
 74   CONTINUE
      NV = NV+1
      R(I) = PARAM(1,IV)
      V(I) = PARAM(2,IV)
      GO TO 50
C
C     ALL DEFINITIONS FOUND AND ORDERED
C
 90   IF ( NV .LE. 2 )  GO TO 950
C                                                                        1300   
      GO TO ( 100, 500, 600 ), IREQUE
C
C     GENERATE RADIUS ARRAY FOR INTRPC
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
      I = NALLOC( NUMPTS, WORKNM(1,1), LOC )
C
      RR = RSTART
      DO 159  II = 1, NUMPTS                                             1310   
         ALLOC(LOC(I)-1+II) = RR
         RR = RR + STEPSZ
 159  CONTINUE
C
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
      IF ( NUMOUT .NE. 0 )  GO TO 900                                    1320   
      RETURN
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  IORDER,  ( R(I), V(I), I = 1, NV )
 503  FORMAT ( ' THE POTENTIAL IS DEFINED BY ',
     1    'LAGRANGE INTERPOLATION OF ORDER', I2,
     2    ' BASED ON THE POINTS:' /
     3   ( 5( F10.3, G13.5 ) ) )
      GO TO 900                                                          1330   
C
C     DO IT
C
 600  LL = LOC(MYINTS(1)) - 1
C
C     SETUP OUTPUT X'S, DO NOT EXTRAPOLATE; RATHER USE
C     F(1) FOR R < R1  AND  EXPONENTIAL DECAY FOR R > R2 .
C
 270  RR = RSTART
      NSTRT = 1                                                          1340   
      DO 299  N = 1, NUMPTS
         IF ( RR .GE. R(1) )  GO TO 280
         ARRAY1(N) = -V(1)
         NSTRT = N+1
         GO TO 290
 280     IF ( RR .LE. R(NV) )  GO TO 285
         ARRAY1(N) = 0
         GO TO 290
 285     NLAST = N
 290     RR = RR + STEPSZ                                                1350   
 299  CONTINUE
      NOUT = NLAST+1-NSTRT
C
C     DOES THE POTENTIAL CHANGE SIGN
C
      LNORDR = IORDER
      VSIGN = -DSIGN(1.0D0, V(NV) )
      LNSTRT = 1
      DO 659  II = 2, NV
         I = NV+2-II                                                     1360   
         IF ( VSIGN*V(I-1) .GE. 0 )  GO TO 660
 659  CONTINUE
      GO TO 700
C
C     YES, USE STRAIGHT INTERPOLATION UP TO THE LAST CHANGE.
C     FIND THIS TRANSISTION POINT IN THE OUTPUT ARRAYS.
C
 660  LNSTRT = I
      NSTRT2 = ( R(LNSTRT)-RSTART )/STEPSZ + 2
      CALL AITKEN ( MIN0(IORDER,LNSTRT-1), 0.D0, 0.D0, LNSTRT, R, V,     1370   
     1  NSTRT2-NSTRT, ALLOC(LL+NSTRT), ARRAY1(NSTRT),
     2  NFAIL, WORST, A )
      DO 679  I = NSTRT, NSTRT2
         ARRAY1(I) = -ARRAY1(I)
 679  CONTINUE
C
      LNORDR = MIN0( IORDER, NV-LNSTRT )
      IF ( LNORDR .LE. 0 )  RETURN
      NSTRT = NSTRT2
      NOUT = NLAST+1-NSTRT2                                              1380   
C
C     INTERPOLATE THE LOGS
C
 700  DO 719  I = LNSTRT, NV
         V(I) = DLOG( DMAX1( DABS(V(I)), 1.D-15 ) )
 719  CONTINUE
C
      CALL AITKEN ( LNORDR, 0.D0, 0.D0, NV-LNSTRT+1,
     1  R(LNSTRT), V(LNSTRT),
     1  NOUT, ALLOC(LL+NSTRT), ARRAY1(NSTRT),                            1390   
     2  NFAIL, WORST, A )
C
      DO 749  I = NSTRT, NLAST
         X = ARRAY1(I)
         X = DMAX1( X, -69.D0 )
         X = DMIN1( X, 69.D0 )
         ARRAY1(I) = VSIGN*DEXP(X)
 749  CONTINUE
C
C                                                                        1400   
C     FILL IN THE END WITH AN EXPONENTIAL DECAY
C
      TERM = ARRAY1(NLAST)
      NLAST = NLAST+1
      IF ( NLAST .GT. NUMPTS )  RETURN
      STEP = TERM/ARRAY1(NLAST-2)
      TERM = TERM*STEP
      IF ( 1-STEP .LT. 1.D-5 )  RETURN
      DO 769  I = NLAST, NUMPTS
         ARRAY1(I) = TERM                                                1410   
         IF ( TERM .GT. -1.E-30 )  RETURN
         TERM = TERM*STEP
 769  CONTINUE
C
      RETURN
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C                                                                        1420   
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** POTENTIAL MUST BE DEFINED ON AT LEAST',
     1  ' TWO POINTS FOR LAGRANGE.' )
      IRETUR = -1
      GO TO 900
C
 960  WRITE (IUNIT, 963)
 963  FORMAT ( ' IPARAM1 MUST BE DEFINED AS ORDER OF INTERPOLATION',
     1  ' FOR LAGRANGE.' )
      IRETUR = -1                                                        1430   
      GO TO 900
C
      END
C*ID* LTSTELP  C109.PTOLEMY.LINKULES                        PER704  14:
cccCNI%'LINKUL' = 'LTSTELp'
      SUBROUTINE LTSTELp ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,                                          1440   
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC,
     5  IPARAM )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(20), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1),
     4  IPARAM(5)                                                        1450   
      INTEGER  FACFR4
      CHARACTER*4 ID
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4
C
C     LOVE, TERASAWA AND SATCHLER  2+ TELP
C
C     COMPUTES THE IMAGINARY POTENTIAL AS
C
C        VI, RI, AI, IMAPOWER   WOODS-SAXON
C      + VSI, RSI, ASI          SURFACE TERM                             1460   
C      + LTS TELP FOR A 2+ STATE.
C
C     THE TELP PART IS CONTROLLED BY
C
C     IPARAM5 - 1 = PROJECTILE EXCITATION
C               2 = TARGET EXCITATION
C               IF NOT DEFINED, THEN  JA @= JB  OR  JBIGA @= JBIGB
C          IS USED TO DETERMIN WHICH NUCLEUS IS EXCITED.
C     PARAM20 - CONTAINS THE G2(ETA,XI) FUDGE FACTOR.  IF NOT DEFINED,
C               THE VALUE 1 IS USED.                                     1470   
C     BETACOULOMB
C
C
C     WE STORE SOME CONSTANTS IN THE BEGINNING OF THE WORK AREA
C     ALLOC(LOC(IWRK)+I):
C     I   STORED
C     0   R2FAC   R2 = R2FAC*RC
C     1   BETACOULOMB
C
C     3/11/79 - MADE FROM PARITWOOD - S.P.                               1480   
C     8/1/79 - IPARAM IS A SEPARATE ARRAY - S.P.
C     1/18/80 - STANDARD MORTRAN - RPG
C
C
C
      CHARACTER*4  WORKNM(2,2) / '*000', 'SHAP', '*000', 'WORK' /
      CHARACTER*8  WORDS(2) / 'PROJCTIL', '  TARGET' /
C
C
C                                                                        1490   
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
      NOTDEF = INTRNL(3)
C
      VI = TEMPVS(2)
      AI = TEMPVS(6)
      RI = FLOAT(47)
      POWER = FLOAT(103)
      VSI = FLOAT(71)                                                    1500   
      RSI = FLOAT(69)
      ASI = FLOAT(68)
      G2 = PARAM(20)
      IF ( G2 .EQ. UNDEF )  G2 = 1
      ICHANW = INTRNL(5)
      AK  = FKANDM(ICHANW)
      RC = FLOAT(45)
      ZP = INTGER(24)
      ZT = INTGER(25)
      ETA = FKANDM(ICHANW+16)                                            1510   
      ALPHA = 1/CNSTNT(8)
      PI = CNSTNT(1)
      HBARC = CNSTNT(6)
      E = FLOAT(12)
      AM = FLOAT(24)
C
      IF ( IREQUE .GT. 1 )  GO TO 200
C
C     IS THIS PROJECTILE OR TARGET EXCITATION.
C     WE USE THIS INFORMATION TO CONSTRUCT R2FAC WHICH MAKES THE         1520   
C     SINGLE NUCLEON RADIUS OUT OF THE SUM (RC).
C
      IPORT = IPARAM(5)
      IF ( IPORT .NE. NOTDEF )  GO TO 50
      IF ( JBLOCK(2) .NE. JBLOCK(3) )  IPORT = 1
      IF ( JBLOCK(4) .NE. JBLOCK(5) )  IPORT = 2
 50   IF ( IPORT .NE. 1  .AND.  IPORT .NE. 2 )  GO TO 930
C      AMP = FLOAT(40)
C      AMT = FLOAT(41)
      R0MASS = FINTRN(10)                                                1530   
      R2FAC = FLOAT(39+IPORT)**(1.D0/3.D0) / R0MASS
C
      IB = NAMLOC('BETACOUL')
      IF ( IB .EQ. 0 )  GO TO 940
C
      IF ( VI .EQ. 0 )  GO TO 65
      IF ( RI .EQ. UNDEF  .OR.  AI .EQ. UNDEF )  GO TO 950
  65  IF ( VSI .EQ. 0 )  GO TO 70
      IF ( RSI .EQ. UNDEF  .OR.  ASI .EQ. UNDEF )  GO TO 945
C                                                                        1540   
C     GENERATE WORK ARRAY TO SAVE THE W.S.
C
  70  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
      IWRK = NALLOC( NUMPTS+10, WORKNM(1,1) )
      LWRK = LOC(IWRK)
      ALLOC(LWRK) = R2FAC
      ALLOC(LWRK+1) = ALLOC(LOC(IB))
      GO TO 300                                                          1550   
C
C     COMPUTE THE CONSTANTS
C     THIS PART IS DONE FOR IREQUE = 2 (PRINTING) AND 3 (NEW POTENTIAL)
C
 200  IWRK = MYINTS(1)
      LWRK = LOC(IWRK)
      IF ( VI .NE. 0  .AND.  ( AI .LE. 0  .OR.  RI .LE. 0 ) )
     1  GO TO 960
      IF ( VSI .NE. 0  .AND.  ( ASI .LE. 0  .OR.  RSI .LE. 0 ) )
     1  GO TO 960                                                        1560   
      RC2 = RC*ALLOC(LWRK)
      BETA = ALLOC(LWRK+1)
      RD = ZP*ZT*ALPHA*HBARC/E
      WP = 3 * (ZP*ZT*ALPHA*BETA)**2 * AM * G2
     1  * RC2**4
     2  / (50*PI*AK)
C
C
 300  GO TO ( 400, 500, 600, 700 ), IREQUE
C                                                                        1570   
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = IWRK
      MYINTS(2) = II
CCCC
CCCC     HERE WE INDICATE THAT THIS IS A L-DEPENDENT POTENTIAL
CCCC
CCC      IRETUR = +1
      IRETUR = 0
C                                                                        1580   
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT
C
 500  IF ( VI .EQ. 0 )  GO TO 510
      R0 = RI/R0MASS
      WRITE (IUNIT, 503)  VI, RI, AI, R0, POWER
 503  FORMAT ( ' W.S. WELL:', T20, G14.3, T37, F7.4,
     1  T49, F7.4, T67, F7.4, T80, F7.3 )                                1590   
 510  IF ( VSI .EQ. 0 )  GO TO 520
      R0 = RSI/R0MASS
      WRITE (IUNIT, 513)  VSI, RSI, ASI, R0
 513  FORMAT ( ' SURFACE ABSORPTION:', T20, G14.3, T37, F7.4,
     1  T49, F7.4, T67, F7.4, T80, F7.3 )
C
 520  WPRC5 = WP/RC**5
      WRITE (IUNIT, 523)  WORDS(IPORT), BETA, G2, RC, RD, WPRC5
 523  FORMAT ( 1X, A8, ' 2+ LTS TELP WITH BETACOULOMB =', F7.4,
     1  5X, 'G2 =', G12.5,  5X, 'RC =', F7.4,                            1600   
     2  5X, 'RD =', F7.4,  5X, 'WP/RC**5 =', G12.5, ' MEV' )
      GO TO 900
C
C     SET UP THE TWO PIECES
C
 600  LL = LWRK+10
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   1, VI, RI, AI, POWER )
      DO 619  I = 1, NUMPTS
         ARRAY1(I) = ALLOC(LL-1+I)                                       1610   
 619  CONTINUE
C
      IF ( VSI .EQ. 0 )  GO TO 650
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   3, VSI, RSI, ASI, 1.D0 )
      DO 629  I = N1, N2
         ARRAY1(I) = ARRAY1(I) + ALLOC(LL-1+I)
 629  CONTINUE
C
C     NOW THE LTS TELP                                                   1620   
C      FIRST THE "COULOMB BRAKING TERM" WHICH IS ARBITRARILY TRUNCATED
C
 650  X = 1/DSQRT(.10D0)
      N1 = (RD-RSTART)/(.90D0*STEPSZ) + 1.5
      N1 = MIN0(N1, NUMPTS)
      DO 654  I = 1, N1
         ALLOC(LL-1+I) = X
 654  CONTINUE
      R = N1*STEPSZ + RSTART
      DO 659  I = N1, NUMPTS                                             1630   
         ALLOC(LL-1+I) = 1/DSQRT(1-RD/R)
         R = R+STEPSZ
 659  CONTINUE
C
C     NOW THE LTS TELP TERM WHICH HAS TWO REGIONS OF DEFINITION
C
      N1 = (RC-RSTART)/STEPSZ + 1.5
      N1 = MIN0(N1, NUMPTS)
      R = RSTART
      DO 679  I = 1, N1                                                  1640   
         ARRAY1(I) = ARRAY1(I)
     1      - (2*WP/(3*RC**9)) * ALLOC(LL-1+I) * R**4
         R = R+STEPSZ
 679  CONTINUE
C
      N1 = N1+1
      IF ( N1 .GT. NUMPTS )  RETURN
      DO 689  I = N1, NUMPTS
         ARRAY1(I) = ARRAY1(I)
     1      - WP*ALLOC(LL-1+I) * ( 1 - (RC/R)**2 * ( 2.D0/7.D0           1650   
     2         + (RC/R)**2 * (1.D0/21.D0) )) / R**5
         R = R+STEPSZ
 689  CONTINUE
      RETURN
C
C
C
C     ADD IN THE L-DEPENDANT PART
C
 700  LWRK = LOC(MYINTS(1))                                              1660   
      LL = LWRK+10
CCC      FAC = ARRAY2(1)*PCOEF
CCC      IF ( TBIT(L,31) )  FAC = -FAC
CCC      DO 719  I = 1, NUMPTS
CCC         ARRAY1(I) = ARRAY1(I) + FAC*ALLOC(LL-1+I)
CCC 719  CONTINUE
      RETURN
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN                                                             1670   
C
 930  WRITE (IUNIT, 933) IPORT
 933  FORMAT ( '0**** IPARAM5 =', I3, ' MUST BE 1 OR 2 TO INDICATE',
     1  ' PROJECTILE OR TARGET EXCITATION FOR LTSTELP.' )
      GO TO 955
 940  WRITE (IUNIT, 943)
 943  FORMAT ( '0**** BETACOULOMB MUST BE DEFINED FOR LTSTELP.' )
      GO TO 955
 945  WRITE (IUNIT, 947)
 947  FORMAT ( '0**** RSI (OR RSI0) AND ASI MUST BE DEFINED.' )          1680   
      GO TO 955
 950  WRITE (IUNIT, 953)
 953  FORMAT ('0**** RI (OR RI0) AND AI MUST BE DEFINED.' )
 955  IRETUR = -1
      GO TO 900
C
 960  IRETUR = -1
      RETURN
C
      END                                                                1690   
C*ID* OHTA     C109.PTOLEMY.LINKULES                        PER704  13:
cccCNI%'LINKUL' = 'OHTA'
      SUBROUTINE OHTA ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR, IUNIT,
     1  NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, WAVR, WAVI,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC,
     5  IPARAM )
C                                                                        1700   
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), WAVR(NUMPTS), WAVI(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(13),
     3  IWAVBL(200), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1),
     4  IPARAM(5)
      INTEGER  FACFR4
C
C     LINKULE FOR INCOMING WAVES WITH FIXED LN DERIV (TIMES K)           1710   
C
C     2/15/78 - FIRST VERSION
C     6/1/78 - CHANGES FOR NEW WAVLJ STRUCTURE - S.P.
C     7/7/78 - CDC VERSION
C     1/16/80 - STANDARD MORTRAN, USE SMLNUM - RPG
C     8/4/83 - ADD IPARAM; CDC MACRO - S.P.
C
      LOGICAL PBUGSW
      EQUIVALENCE ( IBUGSW, PBUGSW )
C                                                                        1720   
c     wavefunction linkules use id to pass the channel number as an
c     integer which may be a problem!
c
ccc      character*4 id
      character*4  worknm(2) /  '*000', 'WORK' /
      DIMENSION UDERIV(2, 4),  FDERIV(2, 6),  BC(5),
     1  FTERMS(2, 6, 2),  FSUMS(2, 2)
      DIMENSION XS(22), BS(22), CS(22), DS(22)
CDC      COMMON /OHTABL/ XS, BS, CS, DS
CDC      ALLOCLEVEL XS, BS, CS, DS                                       1730   
C
      DATA TWELTH / .0833333333333333333333 0D0 /
C
      NUMOUT = 0
      IRETUR = 0
      IBUGSW = IWAVBL(114)
CDC      IBUGSW = IWAVBL(76)
      IBUGSW = IWAVBL(114)                                              cnd
      BIGNUM = CNSTNT(11)
      SMLNUM = CNSTNT(12)                                                1740   
C
      GO TO ( 100, 500, 600 ), IREQUE
C
C     LOCATE THE BOUNDRY RADIUS (PARAM1) AND DERIVATIVE (PARAM2)
C
 100  IF ( PARAM(1) .NE. FINTRN(1)  .AND.  PARAM(2) .NE. FINTRN(1) )
     1    GO TO 120
      WRITE ( IUNIT, 113 )
 113  FORMAT ( '0**** BOTH PARAM1 AND PARAM2 MUST BE DEFINED',
     1  ' AS THE BOUNDARY RADIUS AND DERIVATIVE.' )                      1750   
      IRETUR = -1
      GO TO 900
 120  continue
cccc      WORKNM(1) = ??? char(ID)
CCCC      IWRK = NALLOC( 100,WORKNM)
      UROLD = -1.E-30
      UIOLD = -1.E-30
C
C     SAVE ADDRESS OF WORK ARRAY
C                                                                        1760   
      MYINTS(1) = IWRK
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT
C
 500  DERIV2 = 0
      IF ( PARAM(3) .NE. FINTRN(1) )  DERIV2 = PARAM(3)
      WRITE (IUNIT, 503)  PARAM(1), PARAM(2), DERIV2, FKANDM(ID)
 503  FORMAT ( ' THE WAVEFUNCTION WILL HAVE AN INCOMING',                1770   
     1  ' BOUNDARY CONDITION AT R =', F6.2,
     2  ' FM'  /  '  WITH A LOG DERIVATIVE OF -I * (',
     2  G13.6, ' +', G13.6, 'I ) * K',
     3  '  (K =', F8.3, ' FM**-1) ',
     3  ' FOR ALL VALUES OF L.' )
      GO TO 900
C
C     DO IT
C
C     GET SOME ARGUMENTS INTO THEIR CONVENTIONAL NAMES                   1780   
C
 600  NFIRST = IPOTYP
      H = STEPSZ
      ISTRT = NFIRST-2
      NSTEP = NUMPTS-1
      NWP = ID
CCCC      IWRK = MYINTS(1)
CCCC      LWRK = LOC(IWRK)
      STPSZ = H/FKANDM(NWP)
      RB = PARAM(1)                                                      1790   
      DERIV1 = PARAM(2)
      DERIV2 = 0
      IF ( PARAM(3) .NE. FINTRN(1) )  DERIV2 = PARAM(3)
      DERIV1 = FKANDM(NWP) * DERIV1
      DERIV2 = FKANDM(NWP) * DERIV2
C
C     NOW GET STARTING VALUES OF THE POTENTIAL AND SAVE STARTING
C     VALUES OF THE WAVEFUNCTION.
C
C                                                                        1800   
      N1 = RB/STPSZ + .2
      R1 = N1*STPSZ
      X = RB - R1
      IF ( N1 .GT. ISTRT )  GO TO 250
C
C     THE BOUNDRY IS INSIDE THE FIRST POINT SO IT DOESNT COUNT
C
      GO TO 350
C
C     WE WILL START AT N1                                                1810   
C
 250  ABACK = WAVR(N1+4)
      BBACK = WAVI(N1+4)
C
C     IT MAY BE POSSIBLE TO JUST CORRECT THE PREVIOUS DERIVATIVES
C     AND THUS SAVE THE CUBIC SPLINE TIME
C
      IF ( BBACK .NE. UIOLD )  GO TO 255
      DEL = L*(L+1) - LOLD*(LOLD+1)
      IF ( DABS( UROLD+DEL/(12.D0*N1*N1) - ABACK )  .GT. 1.D-12 )        1820   
     1   GO TO 255
C
C     YES, JUST ADJUST IT
C
      DEL = -DEL/RB
      DO 252  I = 1, 4
         DEL = -I*DEL/RB
         UDERIV(1,I) = UDERIV(1,I) + DEL
 252  CONTINUE
      GO TO 270                                                          1830   
C
C
C     GET THE DERIVATIVES OF THE POTENTIAL FROM CUBIC SPLINES
C
 255  I1 = MAX0( N1-10, ISTRT+1 )
      DO 254  I = 1, 21
         XS(I) = I*STPSZ
 254  CONTINUE
      I = N1-I1+1
      FACTOR = 12/STPSZ**2                                               1840   
      IF ( PBUGSW )  WRITE ( IUNIT, 263 )  ISTRT, N1, R1, X
 263  FORMAT ( ' BARRIER ISTRT, N1, R1, X =', 2I6, 2G15.5 )
      IF (PBUGSW)  NUMOUT = 1
C
C     FIRST THE REAL POTENTIAL DERIVATIVES
C
      CALL SPLNCB ( 21, XS, WAVR(I1+4), BS, CS, DS )
      A = WAVR(N1+4)
      B = BS(I)
      C = CS(I)                                                          1850   
      D = DS(I)
      UDERIV(1,1) = FACTOR*( A-1 + X*(B + X*(C + X*D )))
      UDERIV(1,2) = FACTOR*( B + X*(2*C + 3*D*X))
      UDERIV(1,3) = FACTOR* 2*( C + 3*D*X )
      UDERIV(1,4) = FACTOR * 6*D
C
C     NOW THE IMAG POTENTIAL DERIVATIVES
C
      CALL SPLNCB ( 21, XS, WAVI(I1+4), BS, CS, DS )
      A = WAVI(N1+4)                                                     1860   
      B = BS(I)
      C = CS(I)
      D = DS(I)
      UDERIV(2,1) = FACTOR*( A + X*(B + X*(C + X*D )))
      UDERIV(2,2) = FACTOR*( B + X*(2*C + 3*D*X))
      UDERIV(2,3) = FACTOR* 2*( C + 3*D*X )
      UDERIV(2,4) = FACTOR * 6*D
C
C     SAVE STUFF TO CHECK NEXT TIME FOR CHANGES
C                                                                        1870   
 270  UROLD = ABACK
      UIOLD = BBACK
      LOLD = L
      IF (PBUGSW)  WRITE (IUNIT, 273)  UDERIV
 273  FORMAT ( ' POT DERIV', ( T12, 2( 2G11.3, 3X ) ) )
C
C     NOW GENERATE THE DERIVATIVES OF THE WAVEFUNCTION
C
      FDERIV(1,1) = 1
      FDERIV(2,1) = 0                                                    1880   
      FDERIV(1,2) = +DERIV2 * FDERIV(1,1)
      FDERIV(2,2) = -DERIV1 * FDERIV(1,1)
C
C     WE GENERATE THE HIGHER DERIVATIVES BY RECURSION ON THE
C     SCHROEDINGER EQUATION.
C
      BC(1) = 1
      DO 349  N = 2, 5
         BC(N) = 1
         NM = N-1                                                        1890   
         FDERIV(1,N+1) = 0
         FDERIV(2,N+1) = 0
         DO 329  I = 1, NM
            FDERIV(1,N+1) = FDERIV(1,N+1) - BC(N-I) * (
     1         FDERIV(1,I)*UDERIV(1,N-I) - FDERIV(2,I)*UDERIV(2,N-I) )
            FDERIV(2,N+1) = FDERIV(2,N+1) - BC(N-I) * (
     1         FDERIV(2,I)*UDERIV(1,N-I) + FDERIV(1,I)*UDERIV(2,N-I) )
            IF ( I .NE. NM )  BC(N-I) = BC(N-I) + BC(N-I-1)
 329     CONTINUE
 349  CONTINUE                                                           1900   
C
C     NOW GENERATE THE TERMS AND SUM THEM TO GET F(R1) AND F(R2)
C
      X = -X
      DO 389  II = 1, 2
         DO 379  I = 1, 2
            FSUMS(I,II) = 0
            FAC = 1
            DO 379  N = 1, 6
               FTERMS(I,N,II) = FAC*FDERIV(I,N)                          1910   
               FSUMS(I,II) = FSUMS(I,II) + FTERMS(I,N,II)
               IF ( DABS(FAC) .LT. SMLNUM )  FAC = 0
               FAC = X*FAC/N
 379        CONTINUE
            X = X + STPSZ
 389  CONTINUE
C
      IF (PBUGSW)  WRITE (IUNIT, 393)  FTERMS, FSUMS
 393  FORMAT (    ' FTERMS', 3( T9, 2( 2G11.3, 3X ) / ),
     1            ' FTERMS', 3( T9, 2( 2G11.3, 3X ) / ),                 1920   
     1  ' FSUMS', ( T9, 2( 2G11.3, 3X ) ) )
C
C     NOW START THE INTEGRATION FROM N1
C
      BACKR = 1.D-15 * FSUMS(1,1)
      BACKI = 1.D-15 * FSUMS(2,1)
      THISR = 1.D-15 * FSUMS(1,2)
      THISI = 1.D-15 * FSUMS(2,2)
C
C     FILL IN THE INTERIOR WAVEFUNCTION WITH ZEROS                       1930   
C
      DO 399  I = ISTRT, N1
         WAVR(I+1) = 0
         WAVI(I+1) = 0
 399  CONTINUE
C
      ISTRT = N1
      NFIRST = ISTRT + 2
C
C                                                                        1940   
      WAVR(ISTRT+1) = BACKR
      WAVI(ISTRT+1) = BACKI
      WAVR(ISTRT+2)=THISR
      WAVI(ISTRT+2)=THISI
C
C     COMPUTE FIRST TWO XSI'S
C
      II = ISTRT+1
      DO 309  I = ISTRT, II
         WAVR(I+3) =                                                     1950   
     1      TWELTH*( WAVR(I+1)*WAVR(I+4) - WAVI(I+1)*WAVI(I+4) )
         WAVI(I+3) =
     1      TWELTH*( WAVR(I+1)*WAVI(I+4) + WAVI(I+1)*WAVR(I+4))
 309  CONTINUE
C
C
C     WILL DO INTEGRATION HERE
C
 350  CONTINUE
C                                                                        1960   
C     MODIFIED NUMEROV (A LA RAYNAL) IS USED IN THE FOLLOWING LOOP
C
C     WE DEFINE
C       W(I) = 1 - (H**2/12)*(V(I)/E - 1)
C       XSI(I) = W(I)*U(I)
C     WHERE
C       V(I) = V(I*STEPSZ) = V(I*H/K)
C     CONTAINS THE COMPLETE NUCLEAR, COULOMB AND CENTRIFIGAL POTENTIALS.
C     MOST OF W(I) IS COMPUTED AND STORED IN MAKPOT, HOWEVER
C     THE L-DEPENDENT PIECES ARE ADDED IN JUST ABOVE.                    1970   
C
C     THE DESIRED WAVEFUNCTION IS
C       U(I) = U(I*STEPSZ)
C     EXCEPT THAT AT THIS STAGE THE NORMALIZATION IS ARBITRARY:
C       U(R) ---> CONS*R**(L+1) .
C
C     AT STEP I IN THE LOOP, THE WAVR AND WAVI ARRAYS CONTAIN:
C
C     WAV_(...)   CONTENTS
C          I-1    U(I-2)     NOT USED                                    1980   
C          I      U(I-1)
C          I+1    XSI(I-2),   CHANGED TO U(I) DURING STEP
C          I+2    XSI(I-1)
C          I+3    W(I-1), NOT USED, CHANGED TO XSI(I) DURING STEP
C          I+4    W(I)
C
C     THE SMALLEST POSSIBLE NFIRST IS 2 WHICH RESULTS IN ISTRT=0.
C     THE U AND XSI FOR NFIRST-2, NFIRST-1 HAVE ALREADY BEEN STORED.
C
C                                                                        1990   
C     ON THE FIRST PARTIAL WAVE, WE MUST FROM TIME TO TIME CHECK
C     FOR IMMENENT OVERFLOW AND RENORMALIZE DOWN.  AFTER THE FIRST
C     PARTIAL WAVE WE KNOW WHERE TO START AND NEED NOT CHECK.
C
      N1 = NFIRST
      N2 = NSTEP
      IF ( L .GE. IWAVBL(114+NWP) )  GO TO 420
C
CCCCC     GET INTERVAL INWHICH WE MAY REACH BIGNUM
CCCCC                                                                    2000   
CCCC 410  N2 = N1*(BIGNUM/THISR)**(1./(DL+1))
CCCC      N2 = MAX0( N2, N1+5 )
 410  N2 = N1 + 20
      N2 = MIN0(N2, NSTEP)
C
 420  D = ( 12. / ( WAVR(I+4)**2 + WAVI(I+4)**2 ) )
C
      DO 449  I = N1, N2
C
         WAVR(I+3) = ((-10.)*WAVR(I+2)) + (WAVR(I)-WAVR(I+1))            2010   
         WAVI(I+3) = ((-10.)*WAVI(I+2)) + (WAVI(I)-WAVI(I+1))
C
         WAVR(I+1) = ( WAVR(I+4)*WAVR(I+3) + WAVI(I+4)*WAVI(I+3) ) * D
         WAVI(I+1) = ( WAVR(I+4)*WAVI(I+3) - WAVI(I+4)*WAVR(I+3) ) * D
C
         D = ( 12. / ( WAVR(I+5)**2 + WAVI(I+5)**2 ) )
C
 449  CONTINUE
      IF ( N2 .EQ. NSTEP )  GO TO 800
C                                                                        2020   
C     CHECK FOR NEARNESS OF OVERFLOW
C
      THISR = DABS(WAVR(N2+1)) + DABS(WAVI(N2+1))
      IF ( PBUGSW )  WRITE (IUNIT, 463) NWP, L, N1, N2, THISR
 463  FORMAT ( ' CHECKING SIZE: NWP, L, N1, N2, MAG =',
     1   4I5, G13.3 )
      IF ( PBUGSW )  NUMOUT = 1
C
      IF ( THISR .LT. BIGNUM )  GO TO 490
C                                                                        2030   
C     MUST SCALE DOWN.  WE WILL MULTIPLY EVERYTHING BY SMLNUM.
C     FIRST WE FIND WHERE THIS WILL RESULT IN \U\ = SMLNUM.
C
      THISR = THISR*SMLNUM
      N1 = ISTRT
      DO 479  ISTRT = N1, N2
         IF ( DABS(WAVR(ISTRT+1)) .GE. 1 )  GO TO 480
         WAVR(ISTRT+1) = 0
         WAVI(ISTRT+1) = 0.
 479  CONTINUE                                                           2040   
C
 480  IF ( PBUGSW )  WRITE (IUNIT, 483) ISTRT
 483  FORMAT ( 5X, 'RESCALING: NEW ISTRT =', I6 )
C
      N1 = N2+2
      DO 489  I = ISTRT, N1
         WAVR(I+1) = SMLNUM*WAVR(I+1)
         WAVI(I+1) = SMLNUM*WAVI(I+1)
 489  CONTINUE
C                                                                        2050   
 490  N1 = N2+1
      GO TO 410
C
C     END OF THE INTEGRATION
C
C
 800  IF ( NUMOUT .GT. 0 )  GO TO 900
      RETURN
C
 900  NUMOUT = 1                                                         2060   
cib      REWIND IUNIT
      RETURN
      END
C*ID* PARITWOO C109.PTOLEMY.LINKULES                        PER704  14:
ccc%'LINKUL' = 'PARITWOO'
      SUBROUTINE PARITWOO ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,                                          2070   
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
CDC      ALLOCLEVEL ARRAY1, ARRAY2, ALLOC, ILLOC, FACFR4                 2080   
C
C     PARITWOODSAXON - PARITY DEPENDANT WOODS SAXON
C
C     THE POTENTIAL IS A WOODS SAXON WITH PARITY DEPENDANT DEPTH:
C
C       V = (1+PARAM1*(-1)**L) * WOODS SAXON( V R A )
C
C     6/13/78 - BASED ON FIXEDWOOD - S.P.
C     7/11/78 - CDC VERSION AND IRETUR CORRECTION
C                                                                        2090   
C
C
      character*4 id
      character*4  worknm(2,2) / '*000', 'SHAP', '*000', 'WORK' /
C
C
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
C                                                                        2100   
      PCOEF = PARAM(IPOTYP)
      V = TEMPVS(IPOTYP)
      A = TEMPVS(4+IPOTYP)
      R = FLOAT(43)
      IF ( IPOTYP .EQ. 2 )  R = FLOAT(47)
      POWER = FLOAT(101+IPOTYP)
      IF ( R .EQ. UNDEF  .OR.  A .EQ. UNDEF  .OR.  V .EQ. UNDEF
     1  .OR.  PCOEF .EQ. UNDEF )  GO TO 950
      IF ( R .LE. 0  .OR.  A .LE. 0  )  GO TO 960
C                                                                        2110   
C
      GO TO ( 100, 500, 600, 700 ), IREQUE
C
C     GENERATE WORK ARRAY TO SAVE THE W.S.
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
C
      I = NALLOC( NUMPTS, WORKNM(1,1) )
C                                                                        2120   
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
C
C     HERE WE INDICATE THAT THIS IS A L-DEPENDENT POTENTIAL
C
      IRETUR = +1
C
      IF ( NUMOUT .NE. 0 )  GO TO 900                                    2130   
      RETURN
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  V, R, A, POWER, PCOEF
 503  FORMAT ( ' THE POTENTIAL IS A PARITY WOODS-SAXON:', 3G15.5 ,
     1  5X, 'POWER =', G14.5 /
     1  ' WITH THE PARITY-DEPENDANT DEPTH FACTOR',
     2  '  1 + (-PARAM1)**L:  PARAM1 =', G15.5 )
      GO TO 900                                                          2140   
C
C     SET UP THE TWO PIECES
C
 600  LL = LOC(MYINTS(1))
      CALL WOODSX ( NUMPTS, RSTART, STEPSZ, ALLOC(LL), N1, N2,
     1   1, V, R, A, POWER )
      DO 619  I = 1, NUMPTS
         ARRAY1(I) = ALLOC(LL-1+I)
 619  CONTINUE
C                                                                        2150   
      RETURN
C
C
C     ADD IN THE L-DEPENDANT PART
C
 700  LL = LOC(MYINTS(1))
      FAC = ARRAY2(1)*PCOEF
      IF ( btest(L, 0) )  FAC = -FAC
      DO 719  I = 1, NUMPTS
         ARRAY1(I) = ARRAY1(I) + FAC*ALLOC(LL-1+I)                       2160   
 719  CONTINUE
      RETURN
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** R, V, A, AND PARAM1 MUST BE DEFINED.' )
      IRETUR = -1
      GO TO 900                                                          2170   
C
 960  IRETUR = -1
      RETURN
C
      END
C*ID* RAWITSCH C109.PTOLEMY.LINKULES                        PER704  13:
cccCNI%'LINKUL' = 'RAWITSch'
      SUBROUTINE RAWITSch ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, WAVR, WAVI,                        2180   
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC,
     5  IPARAM )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), WAVR(NUMPTS), WAVI(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(13),                  2190   
     3  IWAVBL(200), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1),
     4  IPARAM(5)
      INTEGER  FACFR4
C
C     LINKULE FOR INCOMING WAVES WITH RAWITSCH PRESCRIPTION
C
C     2/15/78 - FIRST VERSION
C     6/1/78 - CHANGES FOR NEW WAVLJ STRUCTURE - S.P.
C     7/7/78 - CDC VERSION
C     1/16/80 - STANDARD MORTRAN, USE SMLNUM - RPG                       2200   
C     8/4/83 - ADD IPARAM; CDC MACRO - S.P.
C     8/6/83 - RAWITSCH REMADE FROM OHTA - S.P.
C     8/17/85 - REWIND IUNIT ONLY FOR IBM
C
      LOGICAL PBUGSW
      EQUIVALENCE ( IBUGSW, PBUGSW )
C
c     wavefunction linkules use id to pass the channel number as an
c     integer which may be a problem!
c                                                                        2210   
ccc      character*4 id
      character*4  worknm(2) /  '*000', 'WORK' /
      DIMENSION   FSUMS(2, 2)
      COMPLEX*16  CSUM, CK, CF(2)
      EQUIVALENCE ( CF(1), FSUMS(1,1) )
      DIMENSION XS(22), BS(22), CS(22), DS(22)
CDC      COMMON /OHTABL/ XS, BS, CS, DS
CDC      ALLOCLEVEL XS, BS, CS, DS
C
      DATA TWELTH / .0833333333333333333333 0D0 /                        2220   
C
      NUMOUT = 0
      IRETUR = 0
      IBUGSW = IWAVBL(114)
CDC      IBUGSW = IWAVBL(76)
      BIGNUM = CNSTNT(11)
      SMLNUM = CNSTNT(12)
      UNDEF = FINTRN(1)
C
      GO TO ( 100, 500, 600 ), IREQUE                                    2230   
C
C     LOCATE THE BOUNDRY RADIUS (PARAM1)
C
 100  IF ( PARAM(1) .NE. UNDEF  )
     1    GO TO 120
      WRITE ( IUNIT, 113 )
 113  FORMAT ( '0**** PARAM1 MUST BE DEFINED',
     1  ' AS THE BOUNDARY RADIUS.' )
      IRETUR = -1
      GO TO 900                                                          2240   
 120  continue
cccc      WORKNM(1) = ??? char(ID)
CCCC      IWRK = NALLOC( 100,WORKNM)
      UROLD = -1.E-30
      UIOLD = -1.E-30
C
C     SAVE ADDRESS OF WORK ARRAY
C
      MYINTS(1) = IWRK
      IF ( NUMOUT .NE. 0 )  GO TO 900                                    2250   
      RETURN
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  PARAM(1)
 503  FORMAT ( ' THE WAVEFUNCTION WILL HAVE A',
     1  ' RAWITSCHER-TYPE BOUNDARY CONDITION AT R =', F6.2 )
      GO TO 900
C
C                                                                        2260   
C     DO IT
C
C     GET SOME ARGUMENTS INTO THEIR CONVENTIONAL NAMES
C
 600  NFIRST = IPOTYP
      H = STEPSZ
      ISTRT = NFIRST-2
      NSTEP = NUMPTS-1
      NWP = ID
      write (6,*) ' !!!!! channel indicator in linkule:', i20            2270   
CCCC      IWRK = MYINTS(1)
CCCC      LWRK = LOC(IWRK)
      STPSZ = H/FKANDM(NWP)
      RB = PARAM(1)
C
C     NOW GET STARTING VALUES OF THE POTENTIAL AND SAVE STARTING
C     VALUES OF THE WAVEFUNCTION.
C
C
      N1 = RB/STPSZ + .2                                                 2280   
      R1 = N1*STPSZ
      X = RB - R1
      IF ( N1 .GT. ISTRT )  GO TO 250
C
C     THE BOUNDRY IS INSIDE THE FIRST POINT SO IT DOESNT COUNT
C
      GO TO 350
C
C     WE WILL START AT N1
C                                                                        2290   
 250  ABACK = WAVR(N1+4)
      BBACK = WAVI(N1+4)
C
C
C     GET THE DERIVATIVES OF THE POTENTIAL FROM CUBIC SPLINES
C
 255  I1 = MAX0( N1-10, ISTRT+1 )
      DO 254  I = 1, 21
         XS(I) = I*STPSZ
 254  CONTINUE                                                           2300   
      I = N1-I1+1
      FACTOR = 12/STPSZ**2
      IF ( PBUGSW )  WRITE ( IUNIT, 263 )  ISTRT, N1, R1, X
 263  FORMAT ( ' BARRIER ISTRT, N1, R1, X =', 2I6, 2G15.5 )
      IF (PBUGSW)  NUMOUT = 1
C
C     FIRST THE REAL POTENTIAL DERIVATIVES
C
      CALL SPLNCB ( 21, XS, WAVR(I1+4), BS, CS, DS )
      A1 = WAVR(N1+4)                                                    2310   
      B1 = BS(I)
      C1 = CS(I)
      D1 = DS(I)
C
C     NOW THE IMAG POTENTIAL DERIVATIVES
C
      CALL SPLNCB ( 21, XS, WAVI(I1+4), BS, CS, DS )
      A2 = WAVI(N1+4)
      B2 = BS(I)
      C2 = CS(I)                                                         2320   
      D2 = DS(I)
C
C     INITIALIZE COMPUTATION OF EFFECTIVE K
C
      CALL CKSET ( R1, STPSZ, A1, B1, C1, D1, A2, B2, C2, D2 )
C
C     DO RAWITSCHER INTEGRALS BY SIMPLE SIMPSON'S RULE
C
      R = R1
      DO 359  I = 1, 2                                                   2330   
         DEL = (R-RB)/4
         CSUM = CK(RB) + 4*CK(RB+DEL) + 2*CK(RB+2*DEL)
     1       + 4*CK(RB+3*DEL) + CK(R)
         CSUM = CSUM * DEL/3
         CSUM = CDEXP( (0.D0,-1.D0)*CSUM )                              cnv
CVA         CSUM = CEXP( (0.,-1.)*CSUM )
CVA%M
         CF(I) = CSUM / CDSQRT(CK(R))
CVA%F
         R = R + STPSZ                                                   2340   
 359  CONTINUE
C
C
      IF (PBUGSW)  WRITE (IUNIT, 393)  FSUMS
 393  FORMAT (
     1  ' FSUMS', ( T9, 2( 2G11.3, 3X ) ) )
C
C     NOW START THE INTEGRATION FROM N1
C
      BACKR = 1.D-15 * FSUMS(1,1)                                        2350   
      BACKI = 1.D-15 * FSUMS(2,1)
      THISR = 1.D-15 * FSUMS(1,2)
      THISI = 1.D-15 * FSUMS(2,2)
C
C     FILL IN THE INTERIOR WAVEFUNCTION WITH ZEROS
C
      DO 399  I = ISTRT, N1
         WAVR(I+1) = 0
         WAVI(I+1) = 0
 399  CONTINUE                                                           2360   
C
      ISTRT = N1
      NFIRST = ISTRT + 2
C
C
      WAVR(ISTRT+1) = BACKR
      WAVI(ISTRT+1) = BACKI
      WAVR(ISTRT+2)=THISR
      WAVI(ISTRT+2)=THISI
C                                                                        2370   
C     COMPUTE FIRST TWO XSI'S
C
      II = ISTRT+1
      DO 309  I = ISTRT, II
         WAVR(I+3) =
     1      TWELTH*( WAVR(I+1)*WAVR(I+4) - WAVI(I+1)*WAVI(I+4) )
         WAVI(I+3) =
     1      TWELTH*( WAVR(I+1)*WAVI(I+4) + WAVI(I+1)*WAVR(I+4))
 309  CONTINUE
C                                                                        2380   
C
C     WILL DO INTEGRATION HERE
C
 350  CONTINUE
C
C     MODIFIED NUMEROV (A LA RAYNAL) IS USED IN THE FOLLOWING LOOP
C
C     WE DEFINE
C       W(I) = 1 - (H**2/12)*(V(I)/E - 1)
C       XSI(I) = W(I)*U(I)                                               2390   
C     WHERE
C       V(I) = V(I*STEPSZ) = V(I*H/K)
C     CONTAINS THE COMPLETE NUCLEAR, COULOMB AND CENTRIFIGAL POTENTIALS.
C     MOST OF W(I) IS COMPUTED AND STORED IN MAKPOT, HOWEVER
C     THE L-DEPENDENT PIECES ARE ADDED IN JUST ABOVE.
C
C     THE DESIRED WAVEFUNCTION IS
C       U(I) = U(I*STEPSZ)
C     EXCEPT THAT AT THIS STAGE THE NORMALIZATION IS ARBITRARY:
C       U(R) ---> CONS*R**(L+1) .                                        2400   
C
C     AT STEP I IN THE LOOP, THE WAVR AND WAVI ARRAYS CONTAIN:
C
C     WAV_(...)   CONTENTS
C          I-1    U(I-2)     NOT USED
C          I      U(I-1)
C          I+1    XSI(I-2),   CHANGED TO U(I) DURING STEP
C          I+2    XSI(I-1)
C          I+3    W(I-1), NOT USED, CHANGED TO XSI(I) DURING STEP
C          I+4    W(I)                                                   2410   
C
C     THE SMALLEST POSSIBLE NFIRST IS 2 WHICH RESULTS IN ISTRT=0.
C     THE U AND XSI FOR NFIRST-2, NFIRST-1 HAVE ALREADY BEEN STORED.
C
C
C     ON THE FIRST PARTIAL WAVE, WE MUST FROM TIME TO TIME CHECK
C     FOR IMMENENT OVERFLOW AND RENORMALIZE DOWN.  AFTER THE FIRST
C     PARTIAL WAVE WE KNOW WHERE TO START AND NEED NOT CHECK.
C
      N1 = NFIRST                                                        2420   
      N2 = NSTEP
      IF ( L .GE. IWAVBL(114+NWP) )  GO TO 420
C
CCCCC     GET INTERVAL INWHICH WE MAY REACH BIGNUM
CCCCC
CCCC 410  N2 = N1*(BIGNUM/THISR)**(1./(DL+1))
CCCC      N2 = MAX0( N2, N1+5 )
 410  N2 = N1 + 20
      N2 = MIN0(N2, NSTEP)
C                                                                        2430   
 420  D = ( 12. / ( WAVR(I+4)**2 + WAVI(I+4)**2 ) )
C
      DO 449  I = N1, N2
C
         WAVR(I+3) = ((-10.)*WAVR(I+2)) + (WAVR(I)-WAVR(I+1))
         WAVI(I+3) = ((-10.)*WAVI(I+2)) + (WAVI(I)-WAVI(I+1))
C
         WAVR(I+1) = ( WAVR(I+4)*WAVR(I+3) + WAVI(I+4)*WAVI(I+3) ) * D
         WAVI(I+1) = ( WAVR(I+4)*WAVI(I+3) - WAVI(I+4)*WAVR(I+3) ) * D
C                                                                        2440   
         D = ( 12. / ( WAVR(I+5)**2 + WAVI(I+5)**2 ) )
C
 449  CONTINUE
      IF ( N2 .EQ. NSTEP )  GO TO 800
C
C     CHECK FOR NEARNESS OF OVERFLOW
C
      THISR = DABS(WAVR(N2+1)) + DABS(WAVI(N2+1))
      IF ( PBUGSW )  WRITE (IUNIT, 463) NWP, L, N1, N2, THISR
 463  FORMAT ( ' CHECKING SIZE: NWP, L, N1, N2, MAG =',                  2450   
     1   4I5, G13.3 )
      IF ( PBUGSW )  NUMOUT = 1
C
      IF ( THISR .LT. BIGNUM )  GO TO 490
C
C     MUST SCALE DOWN.  WE WILL MULTIPLY EVERYTHING BY SMLNUM.
C     FIRST WE FIND WHERE THIS WILL RESULT IN \U\ = SMLNUM.
C
      THISR = THISR*SMLNUM
      N1 = ISTRT                                                         2460   
      DO 479  ISTRT = N1, N2
         IF ( DABS(WAVR(ISTRT+1)) .GE. 1 )  GO TO 480
         WAVR(ISTRT+1) = 0
         WAVI(ISTRT+1) = 0.
 479  CONTINUE
C
 480  IF ( PBUGSW )  WRITE (IUNIT, 483) ISTRT
 483  FORMAT ( 5X, 'RESCALING: NEW ISTRT =', I6 )
C
      N1 = N2+2                                                          2470   
      DO 489  I = ISTRT, N1
         WAVR(I+1) = SMLNUM*WAVR(I+1)
         WAVI(I+1) = SMLNUM*WAVI(I+1)
 489  CONTINUE
C
 490  N1 = N2+1
      GO TO 410
C
C     END OF THE INTEGRATION
C                                                                        2480   
C
 800  IF ( NUMOUT .GT. 0 )  GO TO 900
      RETURN
C
 900  NUMOUT = 1
CIB      REWIND IUNIT
      RETURN
      END
      FUNCTION CK( R )
C                                                                        2490   
C     COMPUTES THE EFFECTIVE COMPLEX K
C
      IMPLICIT REAL*8 ( A-B, D-H, O-Z ), COMPLEX*16 ( C )
CIB      GENERIC
      COMMON /CKBLOK/  CPARAM(4), R1, FACTOR
C
      X = R-R1
      CSUM = CPARAM(1) + X*( CPARAM(2) + X*( CPARAM(3)
     1   + X*CPARAM(4) ) )
CVA%M                                                                    2500   
      CK = CDSQRT( FACTOR*CSUM )
CVA%F
      RETURN
      END
C
      SUBROUTINE CKSET ( RR1, STEPSZ, A1, B1, DC1, D1, A2, B2,
     1   DC2, D2 )
C
      IMPLICIT REAL*8 ( A-B, D-H, O-Z ), COMPLEX*16 ( C )
CIB      GENERIC                                                         2510   
      COMMON /CKBLOK/  CPARAM(4), R1, FACTOR
CVA%M
      R1 = RR1
      CPARAM(1) = DCMPLX( A1-1, A2 )
      CPARAM(2) = DCMPLX( B1, B2 )
      CPARAM(3) = DCMPLX( DC1, DC2 )
      CPARAM(4) = DCMPLX( D1, D2 )
CVA%F
      FACTOR = 12/STEPSZ**2
      RETURN                                                             2520   
      END
C*ID* REID     C109.PTOLEMY.LINKULES                        PER704  14:
cccCNI%'LINKUL' = 'REID'
      SUBROUTINE REID ( ALIAS, SAVINT, IPOTYP, IREQUE, IRETUR,
     1   IUNIT, NUMOUT, L, JP, RSTART, STEPSZ, NUMPTS,
     2   ARRAY1, ARRAY2, FLT, TEMPVS, PARAM, INTGER, JBLOCK,
     3   ISWTCH, INTRNL, FINTRN, CNSTNT, WAVBL, FKANDM, ID,
     4   ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC, IPARAM )
C
C     LINKULE FOR REID SOFT-CORE DEUTERON WAVE FUNCTION                  2530   
C
C
C     RELEVANT ARGUMENTS ARE:
C
C     ALIAS = LINKULE NAME, USED IN ERROR MESSAGES (INPUT).
C     SAVINT = A TWO-ELEMENT REAL*4 ARRAY TO STORE THE INTEGRALS
C        OF WF**2 AND WF*R**(L+1).
C     IPOTYP = 6 (BOTH WF AND POTENTIAL) OR 1 (POTENTIAL ONLY) (INPUT).
C     IREQUE = 1 FOR INITIALIZATION, 2 FOR PRINTING, 3 FOR CALCULATION
C        (INPUT).                                                        2540   
C     IRETUR = ERROR RETURN CODE:  <0 = ERROR (OUTPUT ).
C     IUNIT = UNIT NO. FOR PRINTED OUTPUT (INPUT).
C     NUMOUT = NO. OF LINES TO BE PRINTED (OUTPUT).
C     L = ORBITAL ANGULAR MOMENTUM, MUST = 0 OR 2 (INPUT).
C     JP = 2*PROJECTILE "TOTAL" ANGULAR MOMENTUM; SET TO
C        1 OR 3 FOR L = 0 OR 2 RESPECTIVELY (OUTPUT).
C     RSTART = STARTING R VALUE (FERMIS) (INPUT).
C     STEPSZ = GRID SPACING (FERMIS) (INPUT).
C     NUMPTS = NO. OF GRID POINTS TO BE CALCULATED (INPUT).
C     ARRAY1 (LENGTH NUMPTS) = PRIMARY OUTPUT ARRAY, SET TO              2550   
C        ARRAY2(1)*THE POTENTIAL IF IPOTYP=1, OR TO THE
C        WAVE FUNCTION IF IPOTYP=6.
C     ARRAY2: IF IPOTYP=1, ARRAY2(1) IS THE FACTOR BY WHICH THE LINKULE
C        MUST MULTIPLY THE POTENTIAL.
C        IF IPOTYP=2, ARRAY2 (LENGTH NUMPTS) IS SET TO
C        (-1/V) * THE POTENTIAL (OUTPUT).
C     FLT(62) = V, SET TO 1.0 (OUTPUT).
C     FLT(43) = R, SET TO 1.0 IF UNDEFINED.
C     FLT(1) = A, SET TO 0.4 IF UNDEFINED.
C     FLT(12) = E, SET TO -2.224644 IF UNDEFINED.                        2560   
C     FLT(24) = AM = REDUCED MASS IN MEV/C**2 (INPUT).
C     FLT(53) = SPAM = SPECTROSCOPIC AMPLITUDE.
C        SET TO SQRT(1-PD) OR SQRT(PD) IF L = 0 OR 2, RESPECTIVELY,
C        IF INITIALLY UNDEFINED, WHERE PD = 0.064696.
C     FLT(54) = SPAMP = PROJECTILE SPAM.  SET TO SPAM IF NEXT = 1.
C     FLT(55) = SPAMT = TARGET SPAM.  SET TO SPAM IF NEXT = 2.
C     PARAM(1) = R SCALE COMPRESSION FACTOR.
C     INTGER(12) = NODES, SET TO 0 (OUTPUT).
C     INTGER(17) = IPRINT = PRINT CONTROL.  DUMPS IF THE 1 DIGIT
C        IS >= 6.                                                        2570   
C     JBLOCK(1) = 2*TOTAL ANGULAR MOMENTUM, SET TO 2 (OUTPUT).
C     JBLOCK(8) = 2*PROJECTILE SPIN, SET TO 1 (OUTPUT).
C     JBLOCK(9) = 2*TARGET SPIN, SET TO 1 (OUTPUT).
C     ISWTCH(11) = NEXT = 1 IF THIS IS THE PROJECTILE, 2 IF TARGET.
C     INTRNL(3) = NOTDEF = CONSTANT FOR UNDEFINED INTEGERS.
C     FINTRN(1) = UNDEF = CONSTANT FOR UNDEFINED REAL NUMBERS (INPUT).
C     CNSTNT(1) = PI (INPUT).
C     CNSTNT(2) = SQRT(4*PI) (INPUT).
C     CNSTNT(6) = HBARC = HBAR*C (MEV FM) (INPUT).
C     CNSTNT(13) = BIGLOG = NATURAL LOG OF A BIG NUMBER WHOSE SQUARE     2580   
C        DOESN'S QUITE OVERFLOW (INPUT).
C
C
C     REFERENCE:  R. V. REID, ANN. OF PHYSICS 50 (1968) 411-448.
C
C
C           THIS LINKULE MAY BE USED IN EITHER OF TWO WAYS:
C
C           TO CALCULATE BOTH THE WAVEFUNCTION AND THE POTENTIAL, THE
C     MINIMUM INPUT IS                                                   2590   
C        "CHANNEL:  N + P = D"
C        (OR AN APPROPRIATE "REACTION" AND "PROJECTILE" OR "TARGET")
C        "R0 = 1   A = 0.5"      (TO AVOID COMPLAINTS)
C        "WAVEFUNCTION = REID"
C        "L = 0;"  (OR "L = 2;")
C     THEN THE LINKULE RETURNS THE WAVEFUNCTION IN ARRAY1, AND
C     (-1/V)*THE POTENTIAL IN ARRAY2.
C           TO CALCULATE ONLY THE EFFECTIVE POTENTIAL, AND LET BOUND
C     CALCULATE THE WAVE FUNCTION, USE "REALPOTENTIAL" INSTEAD OF
C     "WAVEFUNCTION" IN THE INPUT ABOVE, AND ADD "NODES = 0".            2600   
C     THEN THE LINKULE RETURNS ARRAY2(1)*THE POTENTIAL IN ARRAY1.
C     (NOTE:  THE PRESENT VERSION OF  BOUND  (10/79) WILL NOT CALCULATE
C     THE WAVE FUCTION FROM THIS POTENTIAL BECAUSE IT CHOKES ON THE
C     CORE.)
C
C
C           NORMALIZATION:  THE INTEGRAL OF THE SQUARE OF THE WAVE
C     FUNCTION RETURNED IS 1.0.  SPECTROSCOPIC INFORMATION IS ENTERED
C     USING THE KEYWORD "SPAM", WHICH IS NOT ACTUALLY USED UNTIL MUCH
C     LATER IN THE DWBA CALCULATION.  IF "SPAM" IS NOT DEFINED, THIS     2610   
C     LINKULE SUPPLIES DEFAULT VALUES FROM THE ACTUAL WAVE FUNCTIONS:
C     SPAM(L=0) = SQRT(1-PD), SPAM(L=2) = SQRT(PD),
C     WHERE PD = 0.064696, N = 0.8776, RHO = 0.02622,
C     D0 = 1.251*D00 = -125.2 MEV FM**(3/2), D2 = 0.4839 FM**2.
C
C           NOTE:  IN DWBA CALCULATIONS, THE DEFAULT SPECTROSCOPIC
C     AMPLITUDE IS ALWAYS USED UNLESS THE KEYWORD "SPAM" IS USED IN
C     THE DESCRIPTION OF THIS CHANNEL.  SETTING "SPAMP" OR "SPAMT" WILL
C     NOT WORK.  IN STAND-ALONE CALCULATIONS, YOU MUST "UNDEFINE SPAM"
C     TO FORCE USE OF THE DEFAULT IN THE SECOND AND SUBSEQUENT           2620   
C     CALCULATIONS.  F "SPAM" IS USED, REMEMBER THAT A PROLATE
C     DEUTERON REQUIRES A NEGATIVE D-STATE AMPLITUDE.
C
C           R COMPRESSION:  FOR DWBA TESTS, IT IS CONVENIENT TO CON-
C     STRUCT A FAKE DEUTERON WHICH IS SMALLER THAN A REAL DEUTERON,
C     WHILE PRESERVING ITS LEADING LOW-MOMENTUM TERMS.  THIS IS
C     ACCOMPLISHED BY ENTERING THE COMPRESSION FACTOR USING THE
C     "PARAM1" KEYWORD.  FOR EXAMPLE, THE FOLLOWING INPUT (IN ADDITION
C     TO THE INPUT ABOVE) WILL CREATE A "DEUTERON" ONLY HALF AS BIG AS
C     A REAL DEUTERON, BUT WITH THE SAME D(L) AND PROBABILITY INTEGRAL:  2630   
C        "PARAM1 = 2
C        "STEPSIZE = 0.05"  (OR "A = 0.25" IF "STEPSPER" IS DEFINED)
C        "BOUNDASY = 10   LOOKSTEP = 500"
C     THIS SHOULD GIVE DWBA RESULTS SIMILAR TO A NORMAL CALCULATION
C     WITH STEPSIZE=0.1, BOUNDASY=20, AND LOOKSTEP=250, EXCEPT THE
C     FINITE-RANGE EFFECTS SHOULD BE REDUCED BY A FACTOR OF 4.
C     FOR COMPRESSED CALCULATIONS, ALWAYS USE "WAVEFUCNTION",
C     NOT 'REALPOTENTIAL" (EVEN IF THE LATTER WORKS).
C
C                                                                        2640   
C     METHOD:
C           THE WAVE FUNCTIONS (BOTH S AND D) ARE CALCULATED BY
C     INTERPOLATING BETWEEN THE VALUES TABULATED IN THE APPENDIX
C     OF REID'S PAPER.  CUBIC INTERPOLATION IS USED, AS REID SUGGESTS.
C     THE CENTRAL, TENSOR, AND SPIN-ORBIT POTENTIALS ARE CALCULATED
C     USING REID'S EQUATION 30; THEY ARE SUMS OF YUKAWA TERMS.
C     THE "EFFECTIVE POTENTIAL" FOR THE S (D) WAVE IS JUST THE
C     DIAGONAL (CENTRAL, SPIN-ORBIT, AND TENSOR) POTENTIAL, PLUS
C     THE OFF-DIAGONAL TENSOR POTENTIAL MULTIPLIED BY THE D/S (S/D)
C     WAVE FUNCTION RATIO.                                               2650   
C
C
C     7/10/79 - NEW ROUTINE - RPG
C     9/27/79 - REID CREATED FROM MCGEE - RPG
C     11/19/79 - ADD R COMPRESSION FACTOR - RPG
C     12/12/79 - PRINT COMPRESSION FACTOR. - RPG
C     1/18/80 - CNI%'LINKUL' = 'REID' - RPG
C     8/8/80 - FIX UNDEFINED PRINTOUT BUG - S.P.
C
      implicit real*8 ( a-h, o-z )                                      implicit
      CHARACTER*8 ALIAS
      CHARACTER*1 SD(2) / 'S', 'D' /
      REAL*4  SAVINT(2)
c
      DIMENSION  ARRAY1(NUMPTS), ARRAY2(NUMPTS), FLT(152),
     1   TEMPVS(6), PARAM(20), INTGER(50), JBLOCK(12), ISWTCH(23),
     2   INTRNL(74), FINTRN(37), CNSTNT(13), WAVBL(240), FKANDM(26),
     3   ALLOC(1), ILLOC(1), LOC(1), LENGTH(1), IPARAM(5)
CCCC      EXTERNAL NALLOC, NAMLOC
C                                                                        2670   
      LOGICAL L2SW, WFSW, PBUGSW
      DIMENSION  XX(34), UI(34,2), UPI(34,2), YNRMS(2), SPAMS(2)
C
C     THESE WAVE FUNCTION VALUES ARE FROM THE APPENDIX OF REID'S
C     PAPER.  NOTE THAT X = 0.7*R, AND UP = DU/DX.
C
      DATA  XX /
     1  .0100D0, .041250D0, .07250D0, .1350D0, .19750D0,
     2  .2600D0, .32250D0, .3850D0, .44750D0, .5100D0,
     3  .57250D0, .6350D0, .69750D0, .7600D0, .8850D0,                   2680   
     4  1.0100D0, 1.1350D0, 1.2600D0, 1.3850D0, 1.5100D0,
     5  1.7600D0, 2.0100D0, 2.5100D0, 3.0100D0, 3.5100D0,
     6  4.0100D0, 4.5100D0, 5.0100D0, 5.5100D0, 6.0100D0,
     7  7.0100D0, 8.0100D0, 9.0100D0, 10.0100D0   /
      DATA  UI /
     1  .000000D0, .333730D-4, .239010D-3, .276210D-2, .127370D-1,
     2  .360620D-1, .753590D-1, .128470D0, .189930D0, .253490D0,
     3  .313900D0, .367700D0, .413170D0, .449920D0, .499530D0,
     4  .524060D0, .531660D0, .528640D0, .519260D0, .506210D0,
     5  .475050D0, .442000D0, .378640D0, .322490D0, .273990D0,           2690   
     6  .232510D0, .197190D0, .167180D0, .141720D0, .120120D0,
     7  .862900D-1, .619830D-1, .445230D-1, .319820D-1,
     1  .000000D0, .108500D-4, .840730D-4, .103690D-2, .496420D-2,
     2  .144460D-1, .307950D-1, .531570D-1, .789950D-1, .105250D0,
     3  .129330D0, .149580D0, .165290D0, .176450D0, .187100D0,
     4  .186540D0, .179460D0, .169100D0, .157420D0, .145530D0,
     5  .123140D0, .103730D0, .738590D-1, .532930D-1, .390770D-1,
     6  .291150D-1, .220160D-1, .168710D-1, .130790D-1, .102430D-1,
     7  .644120D-2, .415750D-2, .273630D-2, .182770D-2  /
      DATA  UPI /                                                        2700   
     1  .297510D-3, .261270D-2, .123350D-1, .829510D-1, .253260D0,
     2  .500520D0, .750720D0, .933490D0, .101620D1, .100340D1,
     3  .920420D0, .796740D0, .657440D0, .519850D0, .284940D0,
     4  .118650D0, .113970D-1, -.540900D-1, -.925230D-1, -.114180D0,
     5  -.130970D0, -.131930D0, -.120040D0, -.104530D0, -.897090D-1,
     6  -.765100D-1, -.650540D-1, -.552290D-1, -.468500D-1, -.397260D-1,
     7  -.285470D-1, -.205080D-1, -.147310D-1, -.105820D-1,
     1  .742020D-4, .893210D-3, .447550D-2, .318920D-1, .101210D0,
     2  .205960D0, .314770D0, .393710D0, .424660D0, .408420D0,
     3  .357720D0, .288480D0, .214280D0, .144270D0, .332940D-1,          2710   
     4  -.358990D-1, -.731510D-1, -.900910D-1, -.952990D-1, -.941580D-1,
     5  -.839730D-1, -.712820D-1, -.493270D-1, -.339400D-1, -.236120D-1,
     6  -.166910D-1, -.120020D-1, -.877680D-2, -.652010D-2, -.491360D-2,
     7  -.289560D-2, -.177580D-2, -.112270D-2, -.726440D-3  /
      DATA  YNRMS / 0.877580D0, 0.023012780D0 /
C     DEFAULT SPEC. AMPLITUDES (NOTE SIGNS):
      DATA  SPAMS / 0.9671110D0, -0.254354 /
      DATA  D00 / -100.05650D0 /, ROOT8 / 2.8284271250D0 /
C
C                                                                        2720   
C
C     SKIP TO SETUP, PRINTING, OR CALCULATION CODE.
C
      IRETUR = 0
      NUMOUT = 0
      IL = L/2 + 1
      UNDEF = FINTRN(1)
      IF ( IREQUE - 2 ) 100, 200, 300
C
C     SETUP:  CHECK FOR ERROR IN L (ONLY LIKELY ERROR).                  2730   
C
 100  IF ( L .EQ. 0  .OR.  L .EQ. 2 ) GO TO 120
      WRITE ( IUNIT, 103 ) ALIAS, L
 103  FORMAT ( '0****** ERROR IN LINKULE ', A8 /
     1   ' ******  L =', I5, '   L MUST = 0 OR 2' )
      NUMOUT = 2
cib      REWIND IUNIT
      IRETUR = -1
      RETURN
C                                                                        2740   
C     SET JP, NODES, V, R, A, E, J, JSP, JST, SPAM
C
 120  JP = L + 1
      INTGER(12) = 0
C     THE POTENTIAL DEPTH IS ALWAYS SET TO 1.0
      FLT(62) = 1.
      IF ( FLT(43) .EQ. UNDEF ) FLT(43) = 1.
      IF ( FLT(1) .EQ. UNDEF ) FLT(1) = 0.4
      IF ( FLT(12) .EQ. UNDEF ) FLT(12) = -2.224644
      JBLOCK(1) = 2                                                      2750   
      JBLOCK(8) = 1
      JBLOCK(9) = 1
C     THE DEFAULT SPEC. AMPLITUDES GIVE PD = 0.064696
      IF ( FLT(53) .EQ. UNDEF ) FLT(53) = SPAMS(IL)
C     "SPAMP" OR "SPAMT", AS APPROPRIATE, IS ALWAYS SET TO "SPAM".
      I = ISWTCH(11)
      IF ( I .EQ. 1  .OR.  I .EQ. 2 )  FLT(53+I) = FLT(53)
      RETURN
C
C                                                                        2760   
C     IREQUE = 2:  PRINT STUFF ABOUT THE LINKULE.
C
 200  IF ( IPOTYP .NE. 6 ) GO TO 250
C
C     PRINTOUT FOR WAVE FUNCTION CALCULATION.
C
      WRITE ( IUNIT, 203 ) SD(IL), SAVINT(1)
 203  FORMAT ( 10X, 'REID SOFT-CORE DEUTERON WAVE FUNCTION' /
     1   10X, A1, '-STATE PROBABILITY =', F7.4 )
      NUMOUT = 1                                                         2770   
      IF ( PARAM(1) .EQ. UNDEF )  GO TO 220
      WRITE ( IUNIT, 213 )  PARAM(1)
 213  FORMAT ( 10X, 'R SCALE COMPRESSION FACTOR =', F7.4 )
      NUMOUT = 2
 220  IF ( IL .EQ. 2 )  GO TO 230
      T = SAVINT(2)/D00
      WRITE ( IUNIT, 223 )  SAVINT(2), T
 223  FORMAT ( 10X, 'D0 =', F8.2, ' MEV FM**(3/2) =', F7.4, '*D00' )
      GO TO 270
 230  WRITE ( IUNIT, 233 )  SAVINT(2)                                    2780   
 233  FORMAT ( 10X, 'D2*D0 =', F9.4, ' MEV FM**(7/2)' )
      GO TO 270
C
C     PRINTOUT FOR POTENTIAL CALCULATION
C
 250  WRITE ( IUNIT, 253 ) FLT(62)
 253  FORMAT ( 10X, 'EFFECTIVE POTENTIAL DERIVED FROM REID SOFT-CORE',
     1   ' DEUTERON WAVEFUNCTION.' /
     2   10X, 'MULTIPLIER =', F9.6 )
C                                                                        2790   
C     PRINTOUT FOR BOTH WF AND POTENTIAL CALCULATIONS.
C
 270  WRITE (IUNIT, 273 ) FLT(53), YNRMS(IL)
 273  FORMAT ( 10X, 'SPECTROSCOPIC AMPLITUDE ("SPAM") =', F8.5 /
     1   10X, 'ASYMPTOTIC NORMALIZATION =', F8.5, 'FM**(-1/2)' /
     2   10X, 'REFERENCE:  R. V. REID, ANNALS OF PHYSICS 50 (1968)',
     3   ' 411-448' )
      NUMOUT = NUMOUT + 5
cib      REWIND IUNIT
      RETURN                                                             2800   
C
C
C     SET UP FOR THE MAIN CALCULATION LOOP.
C
 300  L2SW = L .NE. 0
      AM = FLT(24)
      E = FLT(12)
      PBUGSW = MOD(INTGER(17), 10) .GE. 6
      RT4PI = CNSTNT(2)
      HBARC = CNSTNT(6)                                                  2810   
      ALSQ = -2.*AM*E/HBARC**2
      ALPHA = DSQRT( ALSQ )
C
C     THE "PARAM1" KEYWORD GIVES THE RADIAL COMPRESSION FACTOR
C     (2 MEANS FACTOR-OF-2 COMPRESSION, I.E. SMALL DEUTERON)
C
      SCALX = 1.0
      IF ( PARAM(1) .NE. UNDEF )  SCALX = PARAM(1)
      SCALW = DSQRT( SCALX )
      SCALV = SCALW*SCALX                                                2820   
      IF ( L2SW )  SCALV = SCALV*SCALX**2
      WFSW = IPOTYP .EQ. 6
      ASSIGN 332 TO JUMPL
      IF ( L2SW ) ASSIGN 334 TO JUMPL
      ASSIGN 370 TO JUMPW
      IF ( WFSW ) ASSIGN 360 TO JUMPW
      ASSIGN 350 TO JUMPP
      IF ( PBUGSW )  ASSIGN 340 TO JUMPP
      IF ( PBUGSW )  WRITE (IUNIT, 313 )
 313  FORMAT ( '0', 5X, 'R', 7X, 'WF', 11X, 'V', 8X, 'V*PHI', 8X,        2830   
     1   'DINT' )
      VMULT = ARRAY2(1)
      IF ( WFSW ) VMULT = -1./FLT(62)
      DMULT = FLT(53)*RT4PI*STEPSZ/(VMULT*SPAMS(IL))
      IF ( L2SW )  DMULT = DMULT/15.
      VMULT = VMULT*SCALV
      WINT = 0
      DINT = 0
      RT = RSTART
      IF ( RT .EQ. 0 )  RT = 1.D-10                                      2840   
      X = 0.7*SCALX*RT
      XSTEP = 0.7*SCALX*STEPSZ
      Y = DEXP( -X )
      YM = DEXP( -XSTEP )
      IEND = (XX(34)-X)/XSTEP
      IEND = MIN0( IEND, NUMPTS )
      IEND = MAX0( IEND, 0 )
      IF ( IEND .EQ. 0 )  GO TO 400
      WF = 0
      POT = 0                                                            2850   
      VPHI = 0
      XXP = XX(1)
      ASSIGN 338 TO JUMPZ
      II = 1
C
C     LOOP THROUGH RADII, UP TO THE HIGHEST TABULATED VALUE.
C
      DO 399  I = 1, IEND
C
C        IF NECESSARY, FIND THE CORRECT INTERVAL AND CALCULATE THE       2860   
C        COEFFICIENTS.
C
         IF ( X .LE. XXP )  GO TO 325
 320        II = II + 1
            XXP = XX(II)
            IF ( X .GT. XXP )  GO TO 320
            IM = II - 1
            XXM = XX(IM)
            H = XXP - XXM
            A1 = UI(IM,1)                                                2870   
            A2 = H*UPI(IM,1)
            UDIF = UI(II,1) - UI(IM,1)
            A3 = 3.*UDIF - H*( 2.*UPI(IM,1) + UPI(II,1) )
            A4 = -2.*UDIF + H*( UPI(IM,1) + UPI(II,1) )
            B1 = UI(IM,2)
            B2 = H*UPI(IM,2)
            UDIF = UI(II,2) - UI(IM,2)
            B3 = 3.*UDIF - H*( 2.*UPI(IM,2) + UPI(II,2) )
            B4 = -2.*UDIF + H*( UPI(IM,2) + UPI(II,2) )
            ASSIGN 330 TO JUMPZ                                          2880   
C
 325     GO TO JUMPZ, ( 330, 338 )
C
C        COEFFICIENTS ARE CALCULATED.  CALCULATE THE DESIRED WAVE
C        FUNCTION WF = EITHER U OR W, AND CALCULATE THE EFFECTIVE
C        POTENTIAL POT FROM THE CENTRAL, TENSOR, AND SPIN-ORBIT
C        PARTS (VC, VT, VLS) AND THE WAVE FUNCTION RATIO.
C
 330     P = (X-XXM)/H
         U = A1 + P*( A2 + P*( A3 + P*A4 ) )                             2890   
         W = B1 + P*( B2 + P*( B3 + P*B4 ) )
         Y2 = Y**2
         Y4 = Y2**2
         Y6 = Y4*Y2
         VC = ( -10.463*Y + 105.468*Y2 - 3187.8*Y4 + 9924.3*Y6 )/X
         VT = ( -10.463*( Y + 3.*( (Y-4.*Y4) + (Y-Y4)/X )/X )
     1      + 351.77*Y4 - 1673.5*Y6 )/X
         GO TO JUMPL, ( 332, 334 )
C
C        S-WAVE EFFECTIVE POTENTIAL AND WAVEFUNCTION                     2900   
C
 332        WF = U
            POT = VC + ROOT8*VT*W/U
            GO TO 336
C
C        D-WAVE EFFECTIVE POTENTIAL AND WAVEFUNCTION
C
 334        WF = W
            VLS = ( 708.91*Y4 - 2713.1*Y6 )/X
            POT = VC - 2.*VT - 3.*VLS + ROOT8*VT*U/W                     2910   
C
 336     POT = POT*VMULT
         WF = WF*SCALW
C
C        ACCUMULATE THE INTEGRALS.  WINT IS THE TOTAL
C        PROBABILITY INTEGRAL, AND DINT IS THE
C        INTEGRAL USED TO CALCULATE D0 OR D2.
C
         WINT = WINT + STEPSZ*WF*WF
         VPHI = WF*POT                                                   2920   
         DTERM = DMULT*VPHI*RT
         IF ( L2SW )  DTERM = DTERM*RT*RT
         DINT = DINT + DTERM
C
C        PRINT IF REQUIRED.
C
 338     GO TO JUMPP, ( 340, 350 )
 340     VPHI = VPHI/RT
         WRITE ( IUNIT, 343 )  RT, WF, POT, VPHI, DINT
 343     FORMAT ( 1X, F8.3, 4G12.4 )                                     2930   
C
C        STORE THE WAVE FUNCTION AND POTENTIAL (IF IPOTYP = 6)
C        OR JUST THE POTENTIAL.
C
 350     GO TO JUMPW, ( 360, 370 )
 360     ARRAY1(I) = WF
         ARRAY2(I) = POT
         GO TO 380
 370     ARRAY1(I) = POT
C                                                                        2940   
C        INCREMENT R AND END LOOP.
C
 380     RT = RT + STEPSZ
         X = X + XSTEP
         Y = Y*YM
 399  CONTINUE
C
 400  IEND = IEND + 1
C
C     IF WE ARE CALCULATING ONLY THE POTENTIAL, ZERO OUT THE REST        2950   
C     OF IT AND QUIT.
C
      IF ( WFSW )  GO TO 430
      IF ( IEND .GT. NUMPTS )  GO TO 600
      DO 419  I = IEND, NUMPTS
         ARRAY1(I) = 0
 419  CONTINUE
      GO TO 600
C
C     WE ARE CALCULATING WAVE FUNCTIONS.                                 2960   
C     FINISH THE INTEGRAL ANALYTICALLY, SINCE THE REST IS JUST AN
C     EXPONENTIAL (TIMES A POLYNOMIAL IN 1/R FOR THE D STATE).
C
 430  XA = SCALX*ALPHA*RT
      YA = SCALW*YNRMS(IL)*DEXP( -XA )
      G = YA
      U = 0.5*YA**2/(SCALX*ALPHA)
      IF ( L2SW )  GO TO 440
C     S-WAVE
      ASSIGN 470 TO JUMPL                                                2970   
      GO TO 450
C     D-WAVE
 440  T = 1./XA
      G = G*( 1. + T*( 3. + 3.*T ) )
      U = U*( 1. + T*( 6. + T*( 12. + T*6. ) ) )
      ASSIGN 460 TO JUMPL
C     BOTH S AND D
 450  WINT = WINT + 0.5*G*G*STEPSZ + U
      SAVINT(1) = WINT
      SAVINT(2) = DINT                                                   2980   
C
C     IN THE NEXT LOOP WE CALCULATE FROM THE END OF THE
C     TABLES ONWARD, WHERE THE NUCLEAR POTENTIALS ARE NEGLIGIBLE.
C
      IF ( IEND .GT. NUMPTS )  GO TO 500
      XP = SCALX*ALPHA*STEPSZ
      YM = DEXP( -XP )
      ASSIGN 490 TO JUMPP
      IF ( PBUGSW )  ASSIGN 480 TO JUMPP
C                                                                        2990   
      DO 499  I = IEND, NUMPTS
         WF = YA
         GO TO JUMPL, ( 460, 470 )
C
C        D-STATE MODIFICATIONS
C
 460     WF = WF*( 1. + ( 3. + 3./XA )/XA )
         XA = XA + XP
C
 470     GO TO JUMPP, ( 480, 490 )                                       3000   
 480     WRITE ( IUNIT, 343 )  RT, WF
         RT = RT + STEPSZ
C
 490     ARRAY1(I) = WF
         ARRAY2(I) = 0
         YA = YA*YM
C
 499  CONTINUE
C
C     THE WAVE FUNCTION CALCULATION IS FINISHED.                         3010   
C     NOW GO BACK AND RENORMALIZE IT SO ITS INTEGRAL IS
C     UNITY, AS PTOLEMY EXPECTS.
C
 500  T = 1./DSQRT( WINT )
      DO 509  I = 1, NUMPTS
         ARRAY1(I) = T*ARRAY1(I)
 509  CONTINUE
C
C     FINISH THE OUTPUT IF NECESSARY, AND QUIT.
C                                                                        3020   
 600  IF ( .NOT. PBUGSW )  RETURN
      WRITE ( IUNIT, 603 )  WINT, DINT
 603  FORMAT ( '0WINT =', F10.7, '   DINT =', F11.6 )
      NUMOUT = NUMPTS + 2
cib      REWIND IUNIT
      RETURN
      END
C*ID* SHAPE    C109.PTOLEMY.LINKULES                        PER704  16:
cccCNI%'LINKUL' = 'SHAPE'
      SUBROUTINE SHAPE ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR, IUNIT,   3030   
     1  NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),         3040   
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
      CHARACTER*4 ID
C
C     SIMPLE LINKULE TO USE A FIXED SHAPE TIMES A DEPTH PARAMETER
C
C     DEPENDING ON THE POTENTIAL TYPE, THE NAME OF THE SHAPE IS
C     ONE OF REALSHAP, IMAGSHAP, REALSOSH, OR IMAGSOSH.  IF NONE
C     OF THESE OBJECTS EXISTS, THEN SHAPE IS USED.                       3050   
C     IF THE CORRESPONDING ARRAY REALSCAL, IMAGSCAL, REALSOSC, IMAGSOSC,
C     OR SHAPESCA IS DEFINED THEN IT MUST HAVE TWO ELEMENTS THAT ARE
C     THE INITIAL AND FINAL COORDINATE VALUES FORWHICH THE POTENTIAL
C     IS GIVEN.  IN THIS CASE THE ARRAY WILL BE INTERPOLATED IF NECESSAR
C     NECESSARY TO THE REQUIRED GRID.
C     ALTERNATIVELY THE SCALE ARRAY MAY HAVE AS MANY ELEMENTS
C     AS THE SHAPE ARRAY.  IN THIS CASE IT GIVE THE R VALUES
C     (POSSIBLY UNEQUALLY SPACED) OF THE ELEMENTS IN SHAPE.
C     IF THERE IS NO SCALE ARRAY TO BE DEFINED ON THE CORRECT GRID.
C                                                                        3060   
C     12/8/77 - FIRST VERSION - S.P.
C     3/15/78 - ALLOW R GRID TO BE SPECIFIED - S.P.
C     1/16/80 - STANDARD MORTRAN - RPG
C     4/14/82 - EXTRAPOLATE FOR LARGE R - S.P.
C
      CHARACTER*8 NAMES(5)  / 'REALSHAP', 'IMAGSHAP', 'REALSOSH',
     1  'IMAGSOSH', 'SHAPE' /
      CHARACTER*8 SCALNM(5) / 'REALSCAL', 'IMAGSCAL', 'REALSOSC',
     1  'IMAGSOSC',  'SHAPESCA' /
C                                                                        3070   
      CHARACTER*4  WORKNM(2,2) / '*000', 'SHAP', '*000', 'WORK' /
      REAL*8 WORK(20)
      LOGICAL RSW
C
C
      SMLNUM = CNSTNT(12)
      NUMOUT = 0
      IRETUR = 0
      CONS = TEMPVS(IPOTYP)
      IF ( CONS .NE. FINTRN(1) )  GO TO 45                               3080   
      WRITE ( IUNIT, 43 )
 43   FORMAT ( '0**** THE WELL DEPTH PARAMETER (V, VI, VSO, VSOI)',
     1   ' MUST BE DEFINED FOR THE SHAPE LINKULE.' )
      IRETUR = -1
      GO TO 900
 45   CONS = -CONS
C
      GO TO ( 100, 500, 600 ), IREQUE
C
C     LOCATE THE SHAPE                                                   3090   
C
 100  II = IPOTYP
      I = NAMLOC(NAMES(II))
      IF ( I .NE. 0 )  GO TO 150
      II = 5
      I = NAMLOC(NAMES(II))
      IF ( I .NE. 0 )  GO TO 150
      WRITE ( IUNIT, 123 )  NAMES(IPOTYP)
 123  FORMAT ('0**** AN OBJECT WITH THE NAME ', A8,
     1  ' OR SHAPE MUST BE DEFINED TO USE THE SHAPE LINKULE.' )          3100   
      IRETUR = -1
      GO TO 900
C
C     IS THE SCALING ARRAY DEFINED
C
 150  ISCAL = NAMLOC( SCALNM(II) )
      RMAX = RSTART + (NUMPTS-1)*STEPSZ
      NIN = LENGTH(I)
      IF ( ISCAL .NE. 0 )  GO TO 200
C                                                                        3110   
C     NO, ALL WE CAN DO IS CHECK THE NUMBER OF POINTS
C
      IF ( NUMPTS .EQ. NIN )  GO TO 180
      WRITE ( IUNIT, 163 )  RMAX, STEPSZ, NUMPTS, NAMES(II), NIN,
     1  SCALNM(II)
 163  FORMAT ( '0**** ASYMPTOPIA AND STEPSIZE =', 2G15.5,
     1  ' AND REQUIRE', I6, ' POINTS.' /
     2  6X, 'HOWEVER, ', A8, ' HAS', I6, ' POINTS AND ',
     3    A8, ' IS NOT DEFINED.' )
      IRETUR = -1                                                        3120   
      GO TO 900
C
 180  WRITE ( IUNIT, 183 )  NAMES(II), RSTART, RMAX, STEPSZ
 183  FORMAT ( '0', A8, ' IS BEING ASSUMED TO BE DEFINED FOR',
     1    G15.5, ' =< R =<', G15.5, ' WITH STEPSIZE =', G15.5 )
      NUMOUT = 1
      GO TO 400
C
C     SCALLING ARRAY IS DEFINED, WILL INTERPOLATION BE REQUIRED
C                                                                        3130   
 200  IF ( LENGTH(ISCAL) .EQ. 2 )  GO TO 210
      IF ( LENGTH(ISCAL) .EQ. NIN )  GO TO 220
      WRITE ( IUNIT, 203 ) SCALNM(II), NIN, LENGTH(ISCAL)
 203  FORMAT ( '0**** SCALLING ARRAY ', A8, ' SHOULD HAVE 2 OR',
     1  I4, ' ELEMENTS,',
     1    ' BUT IT HAS', I6, ' ELEMENTS.' )
      IRETUR = -1
      GO TO 900
C
C     SCALLING ARRAY SPECIFIES START, STOP                               3140   
C
 210  R1 = ALLOC(LOC(ISCAL))
      R2 = ALLOC(LOC(ISCAL)+1)
      RSW = .FALSE.
      IF ( R1 .EQ. RSTART  .AND.  NIN .EQ. NUMPTS  .AND.
     1  DABS(R2-RMAX)/STEPSZ .LT. 1.E-5 )  GO TO 400
      GO TO 240
C
C     SCALLING SPECIFIES COMPLETE R-GRID
C                                                                        3150   
 220  R1 = ALLOC(LOC(ISCAL))
      R2 = ALLOC(LOC(ISCAL)+NIN-1)
      RSW = .TRUE.
C
C     MUST RESCALE
C
 240  IOLD = I
      WORKNM(1,1) = ID
      WORKNM(1,2) = ID
      ISIZE = NIN+NUMPTS                                                 3160   
      I = NALLOC( NUMPTS, WORKNM(1,1), LOC )
      IWORK = NALLOC( ISIZE, WORKNM(1,2), LOC )
      LXIN = LOC(IWORK)-1
      LXOUT = LXIN+NIN
      LYIN = LOC(IOLD)-1
      LYOUT = LOC(I)-1
      IF ( RSW )  LXIN = LOC(ISCAL) - 1
      IF ( RSW )  GO TO 260
C
C     SETUP INPUT X'S                                                    3170   
C
      R = R1
      STEPIN = (R2-R1)/(NIN-1)
      DO 259  N = 1, NIN
         ALLOC(LXIN+N) = R
         R = R+STEPIN
 259  CONTINUE
C
C     GET EXTRAPOLATING FUNCTION
C                                                                        3180   
 260  VSTEP = 0
      AVAL = 0
      YNM1 = ALLOC(LYIN+NIN-1)
      YN = ALLOC(LYIN+NIN)
      IF ( YN*YNM1 .LE. 0 )  GO TO 270
      AVAL = YN
      BVAL = DLOG(YNM1/AVAL) / (R2-ALLOC(LXIN+NIN-1))
      VSTEP = DEXP(-BVAL*STEPSZ)
C
C     SETUP OUTPUT X'S,  USE                                             3190   
C     F(1) FOR R < R1  AND  A*EXP(-B*R) FOR R > R2 .
C
 270  R = RSTART
      NSTRT = 1
      VAL0 = ALLOC(LYIN+1)
      DO 299  N = 1, NUMPTS
         IF ( R .GE. R1 )  GO TO 280
         ALLOC(LYOUT+N) = VAL0
         NSTRT = N+1
         GO TO 290                                                       3200   
 280     IF ( R .LE. R2 )  GO TO 285
         AVAL = AVAL*VSTEP
         IF ( DABS(AVAL) .LT. SMLNUM )  AVAL = 0
         ALLOC(LYOUT+N) = AVAL
         GO TO 290
 285     ALLOC(LXOUT+N) = R
         NLAST = N
 290     R = R + STEPSZ
 299  CONTINUE
      NOUT = NLAST+1-NSTRT                                               3210   
C
C     DO THE INTERPOLATION
C
      CALL AITKEN ( 5, 0.0D0, 0.0D0, NIN, ALLOC(LXIN+1), ALLOC(LYIN+1),
     1  NOUT, ALLOC(LXOUT+NSTRT), ALLOC(LYOUT+NSTRT),
     2  NFAIL, WORST, WORK )
C
C     FREE WORK AREA
C
      LOC(IWORK) = -LOC(IWORK)                                           3220   
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT
C                                                                        3230   
 500  WRITE (IUNIT, 503)  NAMES(MYINTS(2)), CONS
 503  FORMAT ( ' THE POTENTIAL SHAPE STORED IN OBJECT ', A8,
     1  ' IS BEING MULTIPLIED BY', G15.5 )
      GO TO 900
C
C     DO IT
C
 600  LL = LOC(MYINTS(1)) - 1
      DO 629  I = 1, NUMPTS
         ARRAY1(I) = CONS*ALLOC(LL+I)                                    3240   
 629  CONTINUE
      RETURN
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
      END
C*ID* SPLINE   C109.PTOLEMY.LINKULES                        PER704  14:
cccCNI%'LINKUL' = 'SPLINE'
      SUBROUTINE SPLINE ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR, IUNIT,  3250   
     1  NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(2,10), INTGER(42), JBLOCK(12),      3260   
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
      CHARACTER*4  ID
C
C     SPLINE - CONTRUCT POTENTIAL BY SPLINES
C
C     THIS LINKULE ASSUMES THAT THE PAIRS OF NUMBERS:
C       (PARAM1, PARAM2), (PARAM3, PARAM4), ..., (PARAM19, PARAM20)
C     ARE TO BE INTERPERTED AS                                           3270   
C       (R1, V1), (R2, V2), ..., (R10, V10)
C     AND DEFINE A POTENTIAL V(R).  FEWER THAN 10 POINTS MAY
C     BE USED BY LEAVING THE REMAINING V'S OR R'S UNDEFINED.
C     THE R'S DO NOT HAVE TO BE ORDERED.
C     THE POTENTIAL IS DEFINED AS
C
C       V(R) =  V1    R < R1 ,
C               V'(R)   R1 < R < RLAST
C               0       R > RLAST
C                                                                        3280   
C     WHERE V'(R) IS THE RESULT OF CUBIC SPLINE INTERPOLATION
C     OF THE LOG OF THE POTENTIAL.
C
C     4/6/78 - BASED ON SHAPE - S.P.
C     1/16/80 - STANDARD MORTRAN - RPG
C
C
C
      CHARACTER*4  WORKNM(2,2) / '*000', 'SHAP', '*000', 'WORK' /
      DIMENSION  R(20), V(20), B(20), C(20), D(20)                       3290   
C
C
      NUMOUT = 0
      IRETUR = 0
      UNDEF = FINTRN(1)
C
C     GET THE DEFINED POINTS AND VALUES AND ORDER THEM
C
      NV = 1
      IV = 1                                                             3300   
  40  R(1) = PARAM(1,IV)
      V(1) = PARAM(2,IV)
      IF ( R(1) .NE. UNDEF  .AND.  V(1) .NE. UNDEF )  GO TO 50
      IV = IV+1
      IF ( IV .LT. 10 )  GO TO 40
      GO TO 950
C
 50   IV = IV + 1
      IF ( IV .GT. 10 )  GO TO 90
      IF ( PARAM(1,IV) .EQ. UNDEF  .OR.  PARAM(2,IV) .EQ. UNDEF )        3310   
     1  GO TO 50
      DO 59  I = 1, NV
         IF ( PARAM(1,IV) .LT. R(I) )  GO TO 70
 59   CONTINUE
C
C     ADD TO END
C
      NV = NV+1
      R(NV) = PARAM(1,IV)
      V(NV) = PARAM(2,IV)                                                3320   
      GO TO 50
C
C     INSERT POINT; FIRST MOVE REST RIGHT ONE
C
 70   DO 74  II = I, NV
         R(NV+I-II+1) = R(NV+I-II)
         V(NV+I-II+1) = V(NV+I-II)
 74   CONTINUE
      NV = NV+1
      R(I) = PARAM(1,IV)                                                 3330   
      V(I) = PARAM(2,IV)
      GO TO 50
C
C     ALL DEFINITIONS FOUND AND ORDERED
C
 90   IF ( NV .LE. 2 )  GO TO 950
C
      GO TO ( 100, 500, 600 ), IREQUE
C
C     GENERATE RADIUS ARRAY FOR INTRPC                                   3340   
C
 100  WORKNM(1,1) = ID
      WORKNM(1,2) = ID
      I = NALLOC( NUMPTS, WORKNM(1,1), LOC )
C
      RR = RSTART
      DO 159  II = 1, NUMPTS
         ALLOC(LOC(I)-1+II) = RR
         RR = RR + STEPSZ
 159  CONTINUE                                                           3350   
C
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = I
      MYINTS(2) = II
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN
C
C     PRINT IT OUT                                                       3360   
C
 500  WRITE (IUNIT, 503)  ( R(I), V(I), I = 1, NV )
 503  FORMAT ( ' THE POTENTIAL IS DEFINED BY THE POINTS:' /
     1   ( 5( F10.3, G13.5 ) ) )
      GO TO 900
C
C     DO IT
C
 600  LL = LOC(MYINTS(1)) - 1
C                                                                        3370   
C     SETUP OUTPUT X'S, DO NOT EXTRAPOLATE; RATHER USE
C     F(1) FOR R < R1  AND  EXPONENTIAL DECAY FOR R > R2 .
C
 270  RR = RSTART
      NSTRT = 1
      DO 299  N = 1, NUMPTS
         IF ( RR .GE. R(1) )  GO TO 280
         ARRAY1(N) = -V(1)
         NSTRT = N+1
         GO TO 290                                                       3380   
 280     IF ( RR .LE. R(NV) )  GO TO 285
         ARRAY1(N) = 0
         GO TO 290
 285     NLAST = N
 290     RR = RR + STEPSZ
 299  CONTINUE
      NOUT = NLAST+1-NSTRT
C
C     INTERPOLATE THE LOGS
C                                                                        3390   
      DO 719  I = 1, NV
         V(I) = DLOG(V(I))
 719  CONTINUE
C
      CALL SPLNCB ( NV, R, V, B, C, D )
C
      CALL INTRPC ( NV, R, V, B, C, D,
     1  NOUT, ALLOC(LL+NSTRT), ARRAY1(NSTRT) )
C
      DO 749  I = NSTRT, NLAST                                           3400   
         ARRAY1(I) = -DEXP(ARRAY1(I))
 749  CONTINUE
C
C
C     FILL IN THE END WITH AN EXPONENTIAL DECAY
C
      TERM = ARRAY1(NLAST)
      NLAST = NLAST+1
      IF ( NLAST .GT. NUMPTS )  RETURN
      STEP = DEXP( B(NV)*STEPSZ )                                        3410   
      TERM = TERM*STEP
      DO 769  I = NLAST, NUMPTS
         ARRAY1(I) = TERM
         IF ( TERM .GT. -1.E-30 )  RETURN
         TERM = TERM*STEP
 769  CONTINUE
C
      RETURN
C
 900  NUMOUT = 1                                                         3420   
cib      REWIND IUNIT
      RETURN
C
 950  WRITE (IUNIT, 953)
 953  FORMAT ( '0**** POTENTIAL MUST BE DEFINED ON AT LEAST',
     1  ' TWO POINTS FOR SPLINE.' )
      GO TO 900
C
      END
C*ID* TWOSHAPE C109.PTOLEMY.LINKULES                        PER704  14:  3430   
cccCNI%'LINKUL' = 'TWOSHApe'
      SUBROUTINE TWOSHApe ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     2  FLOAT, TEMPVS, PARAM, INTGER, JBLOCK, ISWTCH, INTRNL, FINTRN,
     3  CNSTNT, IWAVBL, FKANDM,
     4  ID, ALLOC, ILLOC, FACFR4, LOC, LENGTH, NALLOC, NAMLOC )
C
      IMPLICIT  REAL*8  ( A-H, O-Z )
C                                                                        3440   
      DIMENSION  MYINTS(2), ARRAY1(NUMPTS), ARRAY2(NUMPTS),
     1  FLOAT(100), TEMPVS(6), PARAM(5), INTGER(42), JBLOCK(12),
     2  ISWTCH(19), INTRNL(22), FINTRN(10), CNSTNT(8),
     3  IWAVBL(50), FKANDM(20), ALLOC(1), ILLOC(1), LOC(1), LENGTH(1)
      INTEGER  FACFR4
      CHARACTER*4  ID
C
C     TWOSHAPE - USE TWO FIXED SHAPES TIMES DEPTH PARAMETERS
C
C     DEPENDING ON THE POTENTIAL TYPE, THE NAMES OF THE SHAPES ARE       3450   
C     ONE OF (REALSHAP, REALADD), (IMAGSHAP, IMAGADD),
C     (REALSOSH, REALSOAD), OR (IMAGSOSH, IMAGSOAD).
C     IF NONE OF THESE OBJECTS EXISTS, THEN (SHAPE, SHAPEADD) IS USED.
C     IF THE CORRESPONDING ARRAY REALSCAL, IMAGSCAL, REALSOSC, IMAGSOSC,
C     OR SHAPESCA IS DEFINED THEN IT MUST HAVE TWO ELEMENTS THAT ARE
C     THE INITIAL AND FINAL COORDINATE VALUES FORWHICH THE SHAPES ARE
C     GIVEN.  IN THIS CASE THE ARRAYS WILL BE INTERPOLATED IF NECESSARY
C     NECESSARY TO THE REQUIRED GRID.
C     NOTE THAT BOTH SHAPES MUST BE GIVEN ON THE SAME GRID.
C     ALTERNATIVELY THE SCALE ARRAY MAY HAVE AS MANY ELEMENTS            3460   
C     AS THE SHAPE ARRAY.  IN THIS CASE IT GIVE THE R VALUES
C     (POSSIBLY UNEQUALLY SPACED) OF THE ELEMENTS IN SHAPE.
C     IF THERE IS NO SCALE ARRAY TO BE DEFINED ON THE CORRECT GRID.
C     THE POTENTIAL USED IS
C
C   -V*REALSHAP - PARAM1*REALADD ,
C       -VI*IMAGSHAP - PARAM2*IMAGADD ,
C      -VSO*REALSOSH - PARAM3*REALSOAD , OR
C     -VSOI*IMAGSOSH - PARAM4*IMAGSOAD .
C                                                                        3470   
C     12/8/77 - FIRST VERSION - S.P.
C     3/15/78 - ALLOW R GRID TO BE SPECIFIED - S.P.
C     3/25/78 - TSOSHAPE MADE FROM SHAPE - S.P.
C     1/16/80 - STANDARD MORTRAN - RPG
C
      CHARACTER*8 NAMES(5)  / 'REALSHAP', 'IMAGSHAP', 'REALSOSH',
     1  'IMAGSOSH', 'SHAPE' /
      CHARACTER*8 NAMES2(5) / 'REALADD', 'IMAGADD', 'REALSOAD',
     1  'IMAGSOAD', 'SHAPEADD' /
      CHARACTER*8 SCALNM(5) / 'REALSCAL', 'IMAGSCAL', 'REALSOSC',        3480   
     1  'IMAGSOSC', 'SHAPESCA' /
C
      CHARACTER*4  WORKNM(2,3) / '*000', 'SHAP', '*000', 'SHP2' ,
     1    '*000', 'WORK' /
      REAL*8 WORK(20)
      LOGICAL RSW
C
c
      NUMOUT = 0
      IRETUR = 0                                                         3490   
      CONS = TEMPVS(IPOTYP)
      CONS2 = PARAM(IPOTYP)
      IF ( CONS .NE. FINTRN(1)  .AND.  CONS2 .NE. FINTRN(1) )  GO TO 45
      WRITE ( IUNIT, 43 )
 43   FORMAT ( '0**** THE WELL DEPTH PARAMETER (V, VI, VSO, VSOI)',
     1   ' AND CORRESPONDINGLY PARAM1 - PARAM4',
     1   ' MUST BE DEFINED FOR THE TWOSHAPE LINKULE.' )
      IRETUR = -1
      GO TO 900
 45   CONS = -CONS                                                       3500   
      CONS2 = -CONS2
C
      GO TO ( 100, 500, 600 ), IREQUE
C
C     LOCATE THE SHAPE
C
 100  II = IPOTYP
      I = NAMLOC(NAMES(II))
      IF ( I .NE. 0 )  GO TO 130
      II = 5                                                             3510   
      I = NAMLOC(NAMES(II))
      IF ( I .NE. 0 )  GO TO 130
      WRITE ( IUNIT, 123 )  NAMES(IPOTYP)
 123  FORMAT ('0**** AN OBJECT WITH THE NAME ', A8,
     1  ' OR TWOSHAPE MUST BE DEFINED TO USE THE SHAPE LINKULE.' )
      IRETUR = -1
      GO TO 900
C
C     NOW LOCATE SECOND OBJECT
C                                                                        3520   
 130  I2 = NAMLOC(NAMES2(II))
      IF (I2 .NE. 0 )  GO TO 140
      WRITE ( IUNIT, 133 ) NAMES2(II)
 133  FORMAT ( '0**** AN OBJECT WITH THE NAME ', A8,
     1  ' MUST BE DEFINED FOR THE TWOSHAPE LINKULE.' )
      IRETUR = -1
      GO TO 900
C
C     CHECK THE TWO ARRAYS FOR COMPATABILITY
C                                                                        3530   
 140  NIN = LENGTH(I)
      IF ( NIN .EQ. LENGTH(I2) )  GO TO 150
      WRITE ( IUNIT, 143 )  NAMES(II), NAMES2(II)
 143  FORMAT ( '0**** ', A8, ' AND ', A8,
     1  ' MUST BE DEFINED ON THE SAME GRID.' )
      IRETUR = -1
      GO TO 900
C
C     IS THE SCALING ARRAY DEFINED
C                                                                        3540   
 150  ISCAL = NAMLOC( SCALNM(II) )
      RMAX = RSTART + (NUMPTS-1)*STEPSZ
      IF ( ISCAL .NE. 0 )  GO TO 200
C
C     NO, ALL WE CAN DO IS CHECK THE NUMBER OF POINTS
C
      IF ( NUMPTS .EQ. NIN )  GO TO 180
      WRITE ( IUNIT, 163 )  RMAX, STEPSZ, NUMPTS, NAMES(II), NIN,
     1  SCALNM(II)
 163  FORMAT ( '0**** ASYMPTOPIA AND STEPSIZE =', 2G15.5,                3550   
     1  ' AND REQUIRE', I6, ' POINTS.' /
     2  6X, 'HOWEVER, ', A8, ' HAS', I6, ' POINTS AND ',
     3    A8, ' IS NOT DEFINED.' )
      IRETUR = -1
      GO TO 900
C
 180  WRITE ( IUNIT, 183 )  NAMES(II), NAMES2(II), RSTART, RMAX, STEPSZ
 183  FORMAT ( '0', A8, ' AND ', A8, ' ARE BEING ASSUMED TO BE DEFINED',
     1  ' FOR',  G15.5, ' =< R =<', G15.5, ' WITH STEPSIZE =', G15.5 )
      NUMOUT = 1                                                         3560   
      GO TO 400
C
C     SCALLING ARRAY IS DEFINED, WILL INTERPOLATION BE REQUIRED
C
 200  IF ( LENGTH(ISCAL) .EQ. 2 )  GO TO 210
      IF ( LENGTH(ISCAL) .EQ. NIN )  GO TO 220
      WRITE ( IUNIT, 203 ) SCALNM(II), NIN, LENGTH(ISCAL)
 203  FORMAT ( '0**** SCALLING ARRAY ', A8, ' SHOULD HAVE 2 OR',
     1  I4, ' ELEMENTS,',
     1    ' BUT IT HAS', I6, ' ELEMENTS.' )                              3570   
      IRETUR = -1
      GO TO 900
C
C     SCALLING ARRAY SPECIFIES START, STOP
C
 210  R1 = ALLOC(LOC(ISCAL))
      R2 = ALLOC(LOC(ISCAL)+1)
      RSW = .FALSE.
      IF ( R1 .EQ. RSTART  .AND.  NIN .EQ. NUMPTS  .AND.
     1  DABS(R2-RMAX)/STEPSZ .LT. 1.E-5 )  GO TO 400                     3580   
      GO TO 240
C
C     SCALLING SPECIFIES COMPLETE R-GRID
C
 220  R1 = ALLOC(LOC(ISCAL))
      R2 = ALLOC(LOC(ISCAL)+NIN-1)
      RSW = .TRUE.
C
C     MUST RESCALE
C                                                                        3590   
 240  IOLD = I
      IOLD2 = I2
      WORKNM(1,1) = ID
      WORKNM(1,2) = ID
      WORKNM(1,3) = ID
      ISIZE = NIN+NUMPTS
      I = NALLOC( NUMPTS, WORKNM(1,1), LOC )
      I2 = NALLOC( NUMPTS, WORKNM(1,2), LOC )
      IWORK = NALLOC( ISIZE, WORKNM(1,3), LOC )
      LXIN = LOC(IWORK)-1                                                3600   
      LXOUT = LXIN+NIN
      LYIN = LOC(IOLD)-1
      LYOUT = LOC(I)-1
      LYIN2 = LOC(IOLD2)-1
      LYOUT2 = LOC(I2)-1
      IF ( RSW )  LXIN = LOC(ISCAL) - 1
      IF ( RSW )  GO TO 270
C
C     SETUP INPUT X'S
C                                                                        3610   
      R = R1
      STEPIN = (R2-R1)/(NIN-1)
      DO 269  N = 1, NIN
         ALLOC(LXIN+N) = R
         R = R+STEPIN
 269  CONTINUE
C
C     SETUP OUTPUT X'S, DO NOT EXTRAPOLATE; RATHER USE
C     F(1) FOR R < R1  AND  0 FOR R > R2 .
C                                                                        3620   
 270  R = RSTART
      NSTRT = 1
      VAL0 = ALLOC(LYIN+1)
      VAL02 = ALLOC(LYIN2+1)
      DO 299  N = 1, NUMPTS
         IF ( R .GE. R1 )  GO TO 280
         ALLOC(LYOUT+N) = VAL0
         ALLOC(LYOUT2+N) = VAL02
         NSTRT = N+1
         GO TO 290                                                       3630   
 280     IF ( R .LE. R2 )  GO TO 285
         ALLOC(LYOUT+N) = 0
         ALLOC(LYOUT2+N) = 0
         GO TO 290
 285     ALLOC(LXOUT+N) = R
         NLAST = N
 290     R = R + STEPSZ
 299  CONTINUE
      NOUT = NLAST+1-NSTRT
C                                                                        3640   
C     DO THE INTERPOLATION
C
      CALL AITKEN ( 5, 0.0D0, 0.0D0, NIN, ALLOC(LXIN+1), ALLOC(LYIN+1),
     1  NOUT, ALLOC(LXOUT+NSTRT), ALLOC(LYOUT+NSTRT),
     2  NFAIL, WORST, WORK )
C
      CALL AITKEN ( 5, 0.0D0, 0.0D0, NIN, ALLOC(LXIN+1), ALLOC(LYIN2+1),
     1  NOUT, ALLOC(LXOUT+NSTRT), ALLOC(LYOUT2+NSTRT),
     2  NFAIL, WORST, WORK )
C                                                                        3650   
C     FREE WORK AREA
C
      LOC(IWORK) = -LOC(IWORK)
C
C     SAVE ADDRESS OF SHAPE ARRAY
C
 400  MYINTS(1) = 256*I + I2
      MYINTS(2) = II
      IF ( NUMOUT .NE. 0 )  GO TO 900
      RETURN                                                             3660   
C
C     PRINT IT OUT
C
 500  WRITE (IUNIT, 503)  CONS, NAMES(MYINTS(2)), CONS2,
     1  NAMES2(MYINTS(2))
 503  FORMAT ( ' THE POTENTIAL IS', T20, G15.5,
     1    ' TIMES THE SHAPE IN ', A8 /
     2  T16, 'PLUS', G15.5, ' TIMES THE SHAPE IN ', A8 )
      GO TO 900
C                                                                        3670   
C     DO IT
C
 600  I = MYINTS(1)/256
      I2 = MYINTS(1) - 256*I
      LL = LOC(I) - 1
      LL2 = LOC(I2) - 1
      DO 629  I = 1, NUMPTS
         ARRAY1(I) = CONS*ALLOC(LL+I) + CONS2*ALLOC(LL2+I)
 629  CONTINUE
      RETURN                                                             3680   
C
 900  NUMOUT = 1
cib      REWIND IUNIT
      RETURN
      END
c=======*******  fitters
C*ID* CMTRIX   C109.PTOLEMY.FITTERS                         PER705  21:
C
C     ----------------------------------------------------------------
C                                                                        3690   
      SUBROUTINE CMTRIX(N,M,RJ,F,RJTJ,RJTF)
C
CDC      ALLOCLEVEL RJ,F,RJTJ,RJTF
      INTEGER I,J,K,M,N
      REAL*8  RJ(M,N),F(M),RJTJ(N,N),RJTF(N)
C
C       THIS ROUTINE FORMS THE MATRIX
C     (JACOBIAN TRANSPOSE) * JACOBIAN  AND THE VECTOR
C     (JACOBIAN TRANSPOSE) * F.
C                                                                        3700   
C     ON INPUT
C
C        N  IS THE DIMENSION OF THE VECTOR  RJTF, THE COLUMN
C           DIMENSION OF THE ARRAY  RJ  AND THE DIMENSION OF THE SQUARE
C           ARRAY  RJTJ.
C
C        M  IS THE DIMENSION OF THE VECTOR  F  AND THE ROW DIMENSION
C           OF THE ARRAY  RJ.
C
C        RJ  CONTAINS THE ELEMENTS OF THE JACOBIAN MATRIX.               3710   
C
C        F  CONTAINS THE ELEMENTS OF THE FUNCTION VECTOR.
C
C     ON OUTPUT
C
C        RJTJ  CONTAINS THE LOWER TRIANGLE OF THE MATRIX
C           (JACOBIAN TRANSPOSE) * JACOBIAN  STORED BY ROWS.
C
C        RJTF  CONTAINS THE ELEMENTS OF THE VECTOR
C           (JACOBIAN TRANSPOSE) * F.                                    3720   
C
  100 DO 200 I = 1, N
      DO 200 J = 1, I
         RJTJ(I,J) = 0.0D0
C
         DO 200 K = 1, M
            RJTJ(I,J) = RJTJ(I,J) + RJ(K,I) * RJ(K,J)
  200 CONTINUE
C
      DO 300 I = 1, N                                                    3730   
         RJTF(I) = 0.0D0
C
         DO 300 K = 1, M
            RJTF(I) = RJTF(I) + F(K) * RJ(K,I)
  300 CONTINUE
C
      RETURN
C     ----- LAST CARD OF SUBROUTINE  CMTRIX  -----
      END
C*ID* DLSMIN   C109.PTOLEMY.FITTERS                         PER705  22:  3740   
C
C     ---------------------------------------------------------------
C
      SUBROUTINE DLSMIN(N,M,X,DISMAX,ACC,MAXFCN,FCN,IERR,F,W,MONIT,IDEB,
     1                  STEPH,BOUND)
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      DOUBLE PRECISION ACC,AD,ANGLE,ANMULT,AP,DELTA,DELTSQ,DISMAX,DM,
     1                 DMULT,DPARM,DS,DTEST,DW,F,FERR,FMIN,FNP,FSQ,      3750   
     2                 GDOTMQ,PJ,PRED,PTM,QPARM,SP,SPARM,SPP,SQFERR,
     3                 SQGRAD,SQMQ,SQPARM,SQSTEP,SS,ST,TEMP,TINC,W,X
      INTEGER I,IC,IERR,IS,J,K,KK,KS,L,M,MAXFCN,MPN,N,NC,ND,NF,NFCALL,
     1        NI,IDEB,NTEST,NT,NTPARM,NV,NW,NX
      DIMENSION X(N),F(M),W(1)
      EXTERNAL FCN,MONIT
      DATA DELTA/Z'3A10000000000000'/
C
C     LEAST SQUARES MINIMIZATION WITHOUT ANALYTIC DERIVATIVES
C                                                                        3760   
C
C       THIS SUBROUTINE IS A MODIFICATION OF  VA05A  CONTAINED IN THE
C     HARWELL  SUBROUTINE LIBRARY. THE ROUTINE WAS AUTHORED BY
C     M.J.D. POWELL.
C
C
C     THIS ROUTINE HAS RECEIVED MINOR MODIFICATIONS BY LARRY NAZARITH
C     AND STEVEN PIEPER AS FOLLOWS:
C       1)  MONIT AND IDEB WERE ADDED;
C       2)  STEPH AND BOUND WERE ADDED;                                  3770   
C       3)  IMPLICIT STATEMENT WAS ADDED.
C
C     10/31/75 - ABOVE MINOR MODIFICATIONS - S.P. AND L.M.
C
C
C       THIS SUBROUTINE MINIMIZES A SUM OF SQUARES  FSQ(X)  OF FUNCTIONS
C     F(J)(X(1),X(2),...,X(N)), J = 1,2,...,M  USING A COMBINATION OF
C     THE LEVENBERG-MARQUARDT METHOD AND THE METHOD OF STEEPEST
C     DESCENTS APPLIED TO THE SUM OF SQUARES  FSQ(X). THE ALGORITHM
C     ASSUMES THE  F(J)  HAVE CONTINUOUS FIRST DERIVITIVES, ESTIMATING   3780   
C     THEM WITH FINITE DIFFERENCES.
C
C     ON INPUT
C
C         N  IS THE NUMBER OF VARIABLES. IT IS ALSO THE DIMENSION OF
C            THE VECTOR  X. A DIAGNOSTIC IS RETURNED FOR  N  LESS
C            THAN 2.
C
C         M  IS THE NUMBER OF FUNCTIONS  F(J). IT IS ALSO THE DIMENSION
C            OF THE VECTOR  F. A DIAGNOSTIC IS RETURNED FOR  M  LESS     3790   
C            THAN N.
C
C         X  CONTAINS AN ESTIMATE OF THE SOLUTION VECTOR
C            (X(1),X(2),...,X(N)).
C
C         DISMAX  IS A GENEROUS ESTIMATE OF THE EUCLIDEAN METRIC
C            DISTANCE OF THE SOLUTION VECTOR FROM THE INITIAL GUESS
C            IN  X.  A DIAGNOSTIC IS RETURNED FOR  DISMAX  .LE.  0.0D0.
C
C         ACC  IS THE ACCURACY REQUESTED IN THE SOLUTION, I.E.‡, A       3800   
C            NORMAL RETURN FROM THIS SUBROUTINE OCCURS WHEN IT IS
C            PREDICTED THAT THE BEST CALCULATED VALUE OF  FSQ  IS
C            NOT MORE THAN  ACC  GREATER THAN THE MINIMUM VALUE.
C
C         MAXFCN  IS A LIMIT ON THE NUMBER OF CALLS TO THE FUNCTION
C            EVALUATION ROUTINE  FCN. NOTE THAT  MAXFCN  IS ALTERED
C            BY THE ROUTINE. THUS THE CORRESPONDING ACTUAL PARAMETER
C            SHOULD NOT BE AN EXPRESSION OR A CONSTANT.
C
C         FCN  IS A USER SUPPLIED SUBROUTINE TO EVALUATE THE FUNCTIONS   3810   
C            F(J), J = 1,2,...,M  AT THE ESTIMATE  X(I),
C            I = 1,2,...,N. THE USER MAY FORCE AN EXIT FROM  DLSMIN
C            BY SETTING  IERR  TO A NEGATIVE INTEGER IN  FCN. FCN,
C            MUST BE DECLARED IN AN  EXTERNAL  STATEMENT IN THE USERS
C            MAIN PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS
C
C            SUBROUTINE FCN(N,M,X,F,G,IERR)
C            DOUBLE PRECISION X(N),F(M),G
C            INTEGER IERR,N,M
C            -----------------------------------                         3820   
C            STATEMENTS TO EVALUATE THE  F(J)  AND TO SET  IERR  IF
C            DESIRED
C            -----------------------------------
C            RETURN
C            END
C
C            THE ARGUMENT  G  IS A DUMMY ARGUMENT INSERTED TO MAKE
C            THE  DLSMIN  USE OF  FCN  COMPATIBLE WITH THAT OF OTHER
C            MINIMIZATION PROGRAMS IN THE  MINPACK  PACKAGE
C        MONIT  IS A USER SUPPLIED MONITOR ROUTINE                       3830   
C        IDEB   = 0  THEN NO MONITORING
C               = 1  THEN MONITOR AFTER EVERY ITERATION
C        STEPH  IS THE STEP LENGTH FOR ESTIMATING THE INITIAL JACOBIAN
C        BOUND  IS A LOWER BOUND ON THIS STEP LENGTH.
C
C     ON OUTPUT
C
C         X  CONTAINS THE BEST AVAILABLE ESTIMATE OF THE SOLUTION VECTOR
C
C         MAXFCN  IS THE NUMBER OF CALLS TO  FCN.                        3840   
C
C         IERR = K GIVES THE FOLLOWING TERMINATION INDICATIONS
C            NORMAL TERMINATION,
C              K = 0,
C            INTERMEDIATE TERMINATION,
C              K = -N(N ANY INTEGER) - USER TERMINATION
C              K = 1 - FAILURE TO CONVERGE IN  MAXFCN  CALLS OF  FCN,
C              K = 2 - N+4 CALLS OF  FCN  FAILED TO IMPROVE RESIDUALS,
C            INITIAL TERMINATION,
C              DUE TO ILLEGAL INPUT PARAMETER VALUES                     3850   
C              K = 10 * I1 + 100 * I2 + 1000 * I3
C              WHERE I1, I2 AND I3 CAN BE ZERO OR 1,
C              BUT AT LEAST ONE OF THEM IS NON-ZERO.
C              ERROR INDICATIONS ARE
C              I1 = 1 - N .LT. 2,
C              I2 = 1 - M .LT N,
C              I3 = 1 - DISMAX .LE. 0.0D0.
C
C         F  CONTAINS THE FUNCTION VALUES CORRESPONDING TO THE  X
C            VECTOR WITH THE EXCEPTION THAT, THE CONTENTS ARE            3860   
C            MEANINGLESS WHEN  IERR  .GE. 10, OR HAVE BEEN DETERMINED
C            BY THE USER WHEN  IERR  .LT. 0.
C
C         W  IS AN ARRAY OF DIMENSION AT LEAST 2*M*(N+1)+N*(2*N+5)
C            USED BY THE SUBROUTINE FOR WORKING SPACE. ON NORMAL OR
C            INTERMEDIATE TERMINATION, THE MOST RECENT APPROXIMATION
C            TO THE JACOBIAN MATRIX IS STORED ROW-WISE AT THE
C            BEGINNING OF THE ARRAY.
C
C                                                                        3870   
C    --------------------------------------------------------------
C
C     TEST THE INPUT PARAMETERS FOR CONSISTANCY
C
  100 IERR = 0
      IF (N .LT. 2) IERR = 10
      IF (M .LT. N) IERR = IERR + 100
      IF (DISMAX .LE. 0.0D0) IERR = IERR + 1000
      IF (IERR .NE. 0) GO TO 730
C                                                                        3880   
C     INITIALIZE THE FOLLOWING PARAMETERS
C
C        SQSTEP -  THE SQUARE OF THE CURRENT STEPSIZE CONSTRAINT
C        DM     -  THE SQUARE OF THE UPPER BOUND ON SQSTEP
C        DELTSQ -  THE SQUARE OF DELTA. THE LOWER BOUND ON SQSTEP
C        DTEST  -  USED IN THE TEST FOR LINEAR INDEPENDENCE OF THE
C                  ELEMENTS OF THE DIRECTION MATRIX
C        IS     -  A PARAMETER CONTROLLING THE LOGICAL FLOW WITHIN
C                  FCNCHK  AND  DLSMIN. IN PARTICULAR
C                  IS = 0  CALLS FOR AN EXIT                             3890   
C                  IS = 1  CALLS FOR A REVISION OF  SQSTEP
C                  IS = 2  CALLS FOR AN UPDATE OF THE JACOBIAN AND THE
C                          GENERALIZED INVERSE
C                  IS = 3  CALLS FOR THE INITIAL FUNCTION AND
C                          JACOBIAN CALCULATIONS
C                  IS = 4  CALLS FOR BYPASSING THE UPDATING OF THE
C                          JACOBIAN, GENERALIZED INVERSE AND DIRECTION
C                          MATRICES BECAUSE THE CURRENT STEP IS TOO
C                          SMALL
C        NFCALL -  THE NUMBER OF CALLS TO THE SUBROUTINE  FCN            3900   
C        NTEST  -  IF  F(X)  DOES NOT DECREASE FOR  NTEST  CONSECUTIVE
C                  TRIES, AN ERROR RETURN WILL BE MADE WITH  IERR = 2
C        TINC   -  USED IN THE CRITERION TO INCREASE STEP SIZE
C        SPARM  -  THE LEAST VALUE OF THE MARQUARDT PARAMETER
C        DPARM  -  USED IN THE REGULATION OF THE MARQUARDT PARAMETER
C
      SQSTEP = 0.0D0
      DM = DISMAX * DISMAX
      DELTSQ = DELTA * DELTA
      DTEST = DFLOAT(N+N) - 0.5D0                                        3910   
      IS = 3
      NFCALL = 0
      NTEST = N + 2
      TINC = 1.0D0
      MPN = M + N
      SPARM = DSQRT(ACC) / DISMAX
      DPARM = 10.0D0 * DM
      ITCNT = 0
C
C     PARTITION THE SCRATCH ARRAY  W                                     3920   
C
      NI = M * N
      NX = NI + MPN * N
      NF = NX + N
      NC = NF + M
      ND = NC + N
      NW = ND + N * N
      NV = NW + N
      NT = NV + M
C                                                                        3930   
C     THE SCRATCH ARRAY  W  HAS BEEN PARTITIONED AS FOLLOWS
C
C        W(1),...,W(M*N)         - CONTAINS AN APPROXIMATION TO THE
C                                  JACOBIAN
C        W(NI+1),...,W(NI+MPN*N) - CONTAINS AN APPROXIMATION TO THE
C                                  GENERALIZED INVERSE OF THE APPENDED
C                                  JACOBIAN
C        W(NX+1),...,W(NX+N)     - CONTAINS THE CURRENT BEST ESTIMATE OF
C                                  THE SOLUTION VECTOR  X
C        W(NF+1),...,W(NF+M)     - CONTAINS THE FUNCTION VALUES AT THE   3940   
C                                  CURRENT BEST  X
C        W(NC+1),...,W(NC+N)     - CONTAINS THE VECTOR OMEGA, WHERE
C                                  OMEGA(N+1-J)=I(J), AND WHERE THE I(J)
C                                  MOST RECENT CORRECTION VECTORS SPAN J
C                                  DIMENSIONS IN X-SPACE
C        W(ND+1),...,W(ND+N*N)   - CONTAINS A DIRECTION MATRIX OF N
C                                  ORTHONORMAL COLUMN VECTORS, THE LAST
C                                  J OF WHICH SPAN THE SPACE CONTAINING
C                                  THE I(J) MOST RECENT CORRECTION
C                                  VECTORS                               3950   
C        W(NW+1),...,W(NW+N)     - CONTAINS A SCRATCH VECTOR
C        W(NV+1),...,W(NV+M)     - CONTAINS A SCRATCH VECTOR
C        W(NT+1),...,W(NT+N)     - CONTAINS A SCRATCH VECTOR
C
C
C     PREPARE FOR THE FIRST ITERATION BY INITIALIZING THE FUNCTION
C       VALUES
C
      CALL FCNCHK(N,M,X,F,FCN,IS,IERR,ACC,FSQ,FMIN,SQSTEP,DELTSQ,
     X            MAXFCN,NFCALL,NTEST,NX,NF,W)                           3960   
      IF (IS .EQ. 0) GO TO 720
      FMIN = FSQ
C
      DO 110 I = 1, N
         W(NX+I) = X(I)
  110 CONTINUE
C
      DO 120 I = 1, M
         W(NF+I) = F(I)
  120 CONTINUE                                                           3970   
C
C     CALCULATE THE INITIAL JACOBIAN APPROXIMATION
C
      DO 220 IC = 1, N
         SAVE1 = DMAX1(DABS(W(NX+IC)*STEPH),BOUND)
         X(IC) = W(NX+IC) + SAVE1
         CALL FCNCHK(N,M,X,F,FCN,IS,IERR,ACC,FSQ,FMIN,SQSTEP,DELTSQ,
     X            MAXFCN,NFCALL,NTEST,NX,NF,W)
         IF (IS .EQ. 0) GO TO 720
         K = IC                                                          3980   
C
         DO 210 I = 1, M
            W(K) = (F(I) - W(NF+I)) / SAVE1
            K = K + N
  210    CONTINUE
C
C     TEST WHETHER THE MOST RECENT X IS BEST
C
         IF (FMIN .GT. FSQ) GO TO 212
         X(IC) = W(NX+IC)                                                3990   
         GO TO 220
  212    W(NX+IC) = X(IC)
C
         DO 214 I = 1, M
            W(NF+I) = F(I)
  214    CONTINUE
C
         FMIN = FSQ
  220 CONTINUE
C                                                                        4000   
C     SET THE DIRECTION MATRIX TO UNITY
C
      K = ND
C
      DO 230 I = 1, N
C
         DO 225 J = 1, N
            K = K + 1
            W(K) = 0.0D0
  225    CONTINUE                                                        4010   
C
         W(K+I-N) = 1.0D0
         W(NC+I) = 1.0D0 + DFLOAT(N-I)
  230 CONTINUE
C
C     SET THE MARQUARDT PARAMETER  QPARM  TO ITS LEAST VALUE
C
  235 QPARM = SPARM
  240 SQPARM = QPARM * QPARM
      NTPARM = 0                                                         4020   
C
C     COPY THE JACOBIAN  J  AND APPEND THE MARQUARDT MATRIX
C      QPARM * I  WHERE  I  IS THE IDENTITY MATRIX. THEN THE
C      GENERALIZED INVERSE OF THE APPENDED JACOBIAN  JA  IS
C      (((JA TRANSPOSE) * JA) INVERSE) * (JA TRANSPOSE) =
C      (((J TRANSPOSE) * J + (QPARM ** 2) * I) INVERSE) * (JA TRANSPOSE)
C      WHICH IS SIMILAR TO THE LEVENBERG-MARQUARDT FORMULATION FOR
C      THE SOLUTION OF THE SYSTEM OF NORMAL EQUATIONS
C
  245 KK = 0                                                             4030   
      K = NI + NI
C
      DO 260 I = 1, N
C
         DO 250 J = 1, M
            KK = KK + 1
            W(KK+NI) = W(KK)
  250    CONTINUE
C
         DO 255 J = 1, N                                                 4040   
            K = K + 1
            W(K) = 0.0D0
  255    CONTINUE
C
      W(K+I-N) = QPARM
  260 CONTINUE
C
C     CALCULATE THE GENERALIZED INVERSE OF THE APPENDED JACOBIAN
C
      CALL GENINV (N,MPN,W(NI+1),N,W(NW+1))                              4050   
C
C     START THE ITERATION BY PREDICTING THE STEEPEST DESCENT AND
C       MARQUARDT STEPS. IN THIS SECTION
C
C       X      - CONTAINS THE PREDICTED STEEPEST DESCENT DIRECTION  G
C                VECTOR(G) = - (JACOBIAN TRANSPOSE) * VECTOR(F(X))
C       F      - CONTAINS THE PREDICTED MARQUARDT CORRECTION VECTOR  MQ
C                VECTOR(MQ) = - (GENERALIZED INVERSE) * VECTOR(F(X))
C       SQGRAD - CONTAINS THE SQUARE OF THE 2-NORM OF  G
C       SQMQ   - CONTAINS THE SQUARE OF THE 2-NORM OF  MQ                4060   
C       GDOTMQ - CONTAINS THE INNER PRODUCT OF  G  AND  MQ
C
      OLDFUN = FMIN
  300 CONTINUE
      TEMFSQ = 0.0D0
      DO 301 I = 1,M
         TEMFSQ = TEMFSQ +  W(NF+I)**2
  301 CONTINUE
      IF ( TEMFSQ .GE. OLDFUN ) GOTO 302
      OLDFUN = TEMFSQ                                                    4070   
      IF ( IDEB .EQ. 1 ) CALL MONIT(ITCNT,NFCALL,0,N,M,W(NX+1),
     1   TEMFSQ,W(NF+1))
  302 ITCNT = ITCNT + 1
      SQGRAD = 0.0D0
      SQMQ = 0.0D0
      GDOTMQ = 0.0D0
C
      DO 310 I = 1, N
         X(I) = 0.0D0
         F(I) = 0.0D0                                                    4080   
         K = I
C
         DO 305 J = 1, M
            X(I) = X(I) - W(K) * W(NF+J)
            F(I) = F(I) - W(NI+K) * W(NF+J)
            K = K + N
  305    CONTINUE
C
         SQGRAD = SQGRAD + X(I) * X(I)
         SQMQ = SQMQ + F(I) * F(I)                                       4090   
         GDOTMQ = GDOTMQ + X(I) * F(I)
  310 CONTINUE
C
C     PREDICT THE REDUCTION IN  F(X)  DUE TO THE MARQUARDT STEP AND
C       CALCULATE THE LENGTH OF THE OPTIMAL STEEPEST DESCENT STEP.
C       IN THIS SECTION
C
C       AD     - CONTAINS THE SUCCESSIVE ELEMENTS OF
C                JACOBIAN * VECTOR(G)
C       AP     - CONTAINS THE SUCCESSIVE ELEMENTS OF                     4100   
C                JACOBIAN * VECTOR(MQ)
C       SQGRAD - CONTAINS THE SQUARE OF THE 2-NORM OF THE OPTIMAL
C                STEEPEST DESCENT STEP
C       PRED   - CONTAINS THE PREDICTED REDUCTION IN  F(X)
C
      K = 0
      PRED = GDOTMQ + GDOTMQ
      DMULT = 0.0D0
C
      DO 325 I = 1, M                                                    4110   
         AP = 0.0D0
         AD = 0.0D0
C
         DO 320 J = 1, N
            K = K + 1
            AP = AP + W(K) * F(J)
            AD = AD + W(K) * X(J)
  320    CONTINUE
C
         PRED = PRED - AP * AP                                           4120   
         DMULT = DMULT + AD * AD
  325 CONTINUE
C
C     TEST WHETHER CONVERGENCE IS PREDICTED
C
      IF (SQMQ .GT. DM) GO TO 326
      AP = DSQRT(SQMQ)
      IF ((PRED + 2.0D0 * SQPARM * AP * (DISMAX - AP)) .LE. ACC) GO TO
     X    700
      GO TO 328                                                          4130   
  326 IF ((PRED + SQPARM * (DM - SQMQ)) .LE. ACC) GO TO 700
C
C     CONVERGENCE IS NOT PREDICTED
C
  328 DMULT = SQGRAD / DMULT
      SQGRAD = SQGRAD * DMULT * DMULT
C
C     BEGIN THE SELECTION OF THE CORRECTION STEP BY TESTING WHETHER
C       THE MARQUARDT STEP IS ACCEPTABLE
C                                                                        4140   
  330 IS = 2
      IF (SQMQ .GT. SQSTEP) GO TO 335
C
C     MARQUARDT STEP IS ACCEPTABLE. TEST THAT THE MARQUARDT PARAMETER
C       HAS ITS LEAST VALUE
C
      IF (QPARM .GT. SPARM) GO TO 235
C
C     MARQUARDT PARAMETER HAS ITS LEAST VALUE. IF MARQUARDT STEP IS
C       NOT TOO SMALL, UPDATE THE DIRECTION MATRIX. OTHERWISE, PREPARE   4150   
C       TO CALCULATE A NEW STEP FOR THIS PURPOSE
C
      SQSTEP = DMAX1(SQMQ,DELTSQ)
      SQGRAD = 0.25D0 * SQMQ
      TINC = 1.0D0
      IF (SQMQ .GE. DELTSQ) GO TO 400
      IS = 4
      GO TO 470
C
C     MARQUARDT STEP IS TOO LARGE. NOW TEST WHETHER TO INCREASE THE      4160   
C       MARQUARDT PARAMETER
C
  335 IF (SQMQ .GT. DPARM) GO TO 340
      NTPARM = 0
      GO TO 350
  340 IF (NTPARM .GT. 0) GO TO 345
      NTPARM = 1
      PTM = SQMQ
      GO TO 350
  345 NTPARM = NTPARM + 1                                                4170   
      PTM = DMIN1(PTM,SQMQ)
      IF (NTPARM .LT. N + 2) GO TO 350
C
C     SET MARQUARDT PARAMETER TO LARGER VALUE
C
      QPARM = QPARM * (PTM / DM) ** 0.25D0
      IF (6.0D0 * SQSTEP .GE. DM) GO TO 240
      AP = DSQRT(PRED / SQMQ)
      IF (AP .LE. QPARM) GO TO 240
      QPARM = DMIN1(AP,QPARM * (DM / (6.0D0 * SQSTEP)) ** 0.25D0)        4180   
      GO TO 240
C
C     TEST WHETHER THE INITIAL VALUE OF SQSTEP HAS BEEN SET
C
  350 IF (SQSTEP .GT. 0.0D0) GO TO 355
      SQSTEP = DMIN1(DM,SQGRAD)
      IF (SQSTEP .GE. DELTSQ) GO TO 355
      SQSTEP = DELTSQ
      GO TO 330
C                                                                        4190   
C     SQSTEP HAS BEEN INITIALIZED. NOW IF SQGRAD .LT. SQSTEP TAKE
C       THE LARGEST ALLOWABLE STEP WHICH IS A LINEAR COMBINATION OF
C       THE OPTIMAL STEEPEST DESCENT AND MARQUARDT STEPS. OTHERWISE TAKE
C       THE LARGEST STEP IN THE STEEPEST DESCENT DIRECTION
C
  355 IF (SQGRAD .LT. SQSTEP) GO TO 360
      ANMULT = 0.0D0
      DMULT = DMULT * DSQRT(SQSTEP / SQGRAD)
      GO TO 365
  360 GDOTMQ = GDOTMQ * DMULT                                            4200   
      ANMULT = (SQSTEP - SQGRAD) / ((GDOTMQ - SQGRAD) + DSQRT((GDOTMQ
     X         - SQSTEP) ** 2 + (SQMQ - SQSTEP) * (SQSTEP - SQGRAD)))
      DMULT = DMULT * (1.0D0 - ANMULT)
C
C     NOW CALCULATE THE PROPOSED STEP AND ITS ANGLE WITH THE FIRST
C       DIRECTION. WHEN FINISHED
C
C       F    - CONTAINS THE PROPOSED STEP
C       SQMQ - CONTAINS THE SQUARE OF THE 2-NORM OF THE STEP
C                                                                        4210   
  365 SQMQ = 0.0D0
      ANGLE = 0.0D0
C
      DO 370 I = 1, N
         F(I) = DMULT * X(I) + ANMULT * F(I)
         SQMQ = SQMQ + F(I) * F(I)
         ANGLE = ANGLE + F(I) * W(ND + I)
  370 CONTINUE
C
      SQGRAD = 0.25D0 * SQMQ                                             4220   
C
C     IF LINEAR INDEPENDENCE OF SUCCESSIVE STEPS IS MAINTAINED, ACCEPT
C       THIS STEP. OTHERWISE, REJECT THIS STEP IN FAVOR OF ONE
C       PROPORTIONAL TO THE FIRST VECTOR IN THE DIRECTION MATRIX AND
C       UPDATE THE DIRECTION MATRIX
C
      IF ((W(NC+1) .LE. DTEST) .OR.
     1    (ANGLE * ANGLE .GE. SQGRAD)) GO TO 400
C
C     CALCULATE A NEW STEP AND UPDATE THE DIRECTION MATRIX               4230   
C
  372 DO 375 I = 1, N
         L = NC + I
         X(I) = W(NX+I) + DELTA * W(ND+I)
         W(L) = W(L+1) + 1.0D0
  375 CONTINUE
C
      W(ND) = 1.0D0
C
      DO 385 I = 1, N                                                    4240   
         K = ND + I
         TEMP = W(K)
C
         DO 380 J = 2, N
            W(K) = W(K+N)
            K = K + N
  380    CONTINUE
C
         W(K) = TEMP
  385 CONTINUE                                                           4250   
C
      CALL FCNCHK(N,M,X,F,FCN,IS,IERR,ACC,FSQ,FMIN,SQSTEP,DELTSQ,
     X            MAXFCN,NFCALL,NTEST,NX,NF,W)
      IF (IS .EQ. 0) GO TO 720
      GO TO 500
C
C     THE CURRENT STEP IS ACCEPTED. EXPRESS THE NEW DIRECTION IN TERMS
C       OF THOSE OF THE DIRECTION MATRIX, UPDATE THE COUNTS IN
C       W(NC+1) ETC.
C                                                                        4260   
  400 TEMP = 0.0D0
      K = ND
C
      DO 420 I = 1, N
         L = NC + I
         X(I) = DW
         DW = 0.0D0
C
         DO 405 J = 1, N
            K = K + 1                                                    4270   
            DW = DW + F(J) * W(K)
  405    CONTINUE
C
         IF (IS .EQ. 1) GO TO 410
         W(L) = W(L) + 1.0D0
         TEMP = TEMP + DW * DW
         IF (TEMP .LE. SQGRAD) GO TO 420
         IS = 1
         KK = I
         X(1) = DW                                                       4280   
         GO TO 415
  410    X(I) = DW
  415    W(L) = W(L+1) + 1.0D0
  420 CONTINUE
C
      W(ND) = 1.0D0
C
C     REORDER THE DIRECTIONS SO THAT  KK  IS FIRST
C
      IF (KK .LE. 1) GO TO 435                                           4290   
      KS = NC + KK * N
C
      DO 430 I = 1, N
         K = KS + I
         TEMP = W(K)
C
         DO 425 J = 2, KK
            W(K) = W(K-N)
            K = K - N
  425    CONTINUE                                                        4300   
C
         W(K) = TEMP
  430 CONTINUE
C
C     GENERATE THE NEW ORTHOGONAL DIRECTION MATRIX
C
  435 DO 440 I = 1, N
         W(NW+I) = 0.0D0
  440 CONTINUE
C                                                                        4310   
      TEMP = X(1) * X(1)
      K = ND
C
      DO 450 I = 2, N
         SQGRAD = DSQRT(TEMP * (TEMP + X(I) * X(I)))
         DW = TEMP / SQGRAD
         SQGRAD = X(I) / SQGRAD
         TEMP = TEMP + X(I) * X(I)
C
         DO 445 J = 1, N                                                 4320   
            L = NW + J
            K = K + 1
            W(L) = W(L) + X(I-1) * W(K)
            W(K) = DW * W(K+N) - SQGRAD * W(L)
  445    CONTINUE
C
  450 CONTINUE
C
      TEMP = 1.0D0 / DSQRT(SQMQ)
C                                                                        4330   
      DO 455 I = 1, N
         K = K + 1
         W(K) = TEMP * F(I)
  455 CONTINUE
C
C
C     CALCULATE THE NEW  X  AND RELATED QUANTITIES. IN THIS SECTION
C
C       X       - CONTAINS THE OLD  X  PLUS THE CORRECTION STEP
C       W(NW+I) - CONTAINS PREDICTED FUNCTION VALUES FOR THE NEW  X      4340   
C                 VECTOR(W) = VECTOR(F(X)) +
C                 JACOBIAN * VECTOR(CORRECTION STEP)
C       FNP     - CONTAINS VECTOR(W) * VECTOR(W)
C
      FNP = 0.0D0
      K = 0
C
      DO 465 I = 1, M
         L = NW + I
         W(L) = W(NF+I)                                                  4350   
C
         DO 460 J = 1, N
            K = K + 1
            W(L) = W(L) + W(K) * F(J)
  460    CONTINUE
C
         FNP = FNP + W(L) * W(L)
  465 CONTINUE
C
  470 DO 475 I = 1, N                                                    4360   
         X(I) = W(NX+I) + F(I)
  475 CONTINUE
C
      CALL FCNCHK(N,M,X,F,FCN,IS,IERR,ACC,FSQ,FMIN,SQSTEP,DELTSQ,
     X            MAXFCN,NFCALL,NTEST,NX,NF,W)
      IF (IS .EQ. 0) GO TO 720
      IF ((IS .EQ. 2) .OR.
     1    (IS .EQ. 4)) GO TO 500
C
C     TEST WHETHER TO DECREASE THE STEPSIZE CONSTRAINT                   4370   
C
      DMULT = 0.9D0 * FMIN + 0.1D0 * FNP - FSQ
      IF (DMULT .GE. 0.0D0) GO TO 480
C
C     DECREASE THE STEPSIZE CONSTRAINT, BUT NOT BELOW THE MINIMUM.
C       ALSO SET TINC = 1.0
C
      SQSTEP = DMAX1(DELTSQ,0.25D0 * SQSTEP)
      TINC = 1.0D0
      IF (FSQ .GE. FMIN) GO TO 520                                       4380   
      GO TO 502
C
C     STEPSIZE CONSTRAINT IS NOT TO BE DECREASED. PREPARE FOR POSSIBLE
C       INCREASE IN THE CONSTRAINT BY CALCULATING A TENTATIVE
C       MULTIPLIER (.LE. 2) FOR THE STEPSIZE
C
  480 TEMP = 0.0D0
      SQFERR = 0.0D0
C
      DO 485 I = 1, M                                                    4390   
         FERR = F(I) - W(NW+I)
         TEMP = TEMP + DABS(F(I) * FERR)
         SQFERR = SQFERR + FERR * FERR
  485 CONTINUE
C
      PJ = 1.0D0 + DMULT / (TEMP + DSQRT(TEMP*TEMP+DMULT*SQFERR))
      TEMP = DMIN1(4.0D0,TINC,PJ)
C
C     INCREASE THE STEPSIZE, BUT NOT ABOVE THE MAXIMUM
C                                                                        4400   
      TINC = PJ / TEMP
      SQSTEP = DMIN1(DM,TEMP * SQSTEP)
      GO TO 502
C
C     IF  F(X)  IMPROVES STORE THE NEW VALUE OF  X. OTHERWISE COMPUTE
C       ANOTHER  X
C
  500 IF (FSQ .GE. FMIN) GO TO 515
C
C     STORE THE NEW  X  AND  F(X)  BY INTERCHANGING WITH THE OLD         4410   
C
  502 FMIN = FSQ
C
      DO 505 I = 1, N
         L = NX + I
         TEMP = X(I)
         X(I) = W(L)
         W(L) = TEMP
  505 CONTINUE
C                                                                        4420   
      DO 510 I = 1, M
         L = NF + I
         TEMP = F(I)
         F(I) = W(L)
         W(L) = TEMP
  510 CONTINUE
C
  515 IF ((IS .EQ. 1) .OR.
     1    (IS .EQ. 2)) GO TO 520
      IS = 2                                                             4430   
      GO TO 372
C
C     CALCULATE THE INCREMENTS IN  X  AND IN  F
C
  520 DS = 0.0D0
C
      DO 525 I = 1, N
         X(I) = X(I) - W(NX+I)
         DS = DS + X(I) * X(I)
  525 CONTINUE                                                           4440   
C
      DO 530 I = 1, M
         F(I) = F(I) - W(NF+I)
  530 CONTINUE
C
C     CALCULATE THE GENERALIZED INVERSE TIMES THE CHANGE IN X
C
      K = NI
      SS = 0.0D0
C                                                                        4450   
      DO 540 I = 1, MPN
         TEMP = 0.0D0
C
         DO 535 J = 1, N
            K = K + 1
            TEMP = TEMP + W(K) * X(J)
  535    CONTINUE
C
         W(NV+I) = TEMP
         SS = SS + TEMP * TEMP                                           4460   
  540 CONTINUE
C
C     CALCULATE THE JACOBIAN TIMES THE CHANGE IN  F. ALSO APPLY
C       PROJECTION TO THE GENERALIZED INVERSE
C
      DO 560 I = 1, N
         ST = 0.0D0
         K = NI + I
C
         DO 545 J = 1, MPN                                               4470   
            ST = ST + W(K) * W(NV+J)
            K = K + N
  545    CONTINUE
C
         ST = ST / SS
         K = NI + I
C
         DO 550 J = 1, MPN
            W(K) = W(K) - ST * W(NV+J)
            K = K + N                                                    4480   
  550    CONTINUE
C
         ST = SQPARM * X(I)
         K = I
C
         DO 555 J = 1, M
            ST = ST + W(K) * F(J)
            K = K + N
  555    CONTINUE
C                                                                        4490   
         W(NW+I) = ST
  560 CONTINUE
C
C     UPDATE THE APPROXIMATION TO THE JACOBIAN AND CALCULATE A ROW
C       VECTOR USED IN THE CORRECTION TO THE GENERALIZED INVERSE
C
  600 IC = 0
      K = 0
      KK = NI
      SP = 0.0D0                                                         4500   
      SPP = 0.0D0
C
      DO 615 I = 1, M
         SS = F(I)
         ST = F(I)
C
         DO 605 J = 1, N
            IC = IC + 1
            KK = KK + 1
            SS = SS - W(IC) * X(J)                                       4510   
            ST = ST - W(KK) * W(NW+J)
  605    CONTINUE
C
         SS = SS / DS
         W(NV+I) = ST
         SP = SP + F(I) * ST
         SPP = SPP + ST * ST
C
         DO 610 J = 1, N
            K = K + 1                                                    4520   
            W(K) = W(K) + SS * X(J)
  610    CONTINUE
C
  615 CONTINUE
C
      DO 625 I = 1, N
         ST = QPARM * X(I)
C
         DO 620 J = 1, N
            KK = KK + 1                                                  4530   
            ST = ST - W(KK) * W(NW+J)
  620    CONTINUE
C
         W(NT+I) = ST
         SP = SP + QPARM * X(I) * ST
         SPP = SPP + ST * ST
  625 CONTINUE
C
C     IF THE SCALAR PRODUCT  SPP  IS SUFFICIENTLY ACCURATE UPDATE THE
C       GENERALIZED INVERSE. OTHERWISE CALCULATE A NEW GENERALIZED       4540   
C       INVERSE AND BEGIN A NEW ITERATION
C
      IF ((0.01D0 + SPP) .LE. (DABS(SP - SPP))) GO TO 245
C
C     UPDATE THE APPROXIMATION TO THE GENERALIZED INVERSE
C
      DO 645 I = 1, N
         K = NI + I
         ST = X(I)
C                                                                        4550   
         DO 630 J = 1, M
            ST = ST - W(K) * F(J)
            K = K + N
  630    CONTINUE
C
         SS = 0.0D0
C
         DO 635 J = 1, N
            SS = SS + W(K) * X(J)
            K = K + N                                                    4560   
  635    CONTINUE
C
         ST = (ST - QPARM * SS) / SP
         K = NI + I
C
         DO 640 J = 1, MPN
            W(K) = W(K) + ST * W(NV+J)
            K = K + N
  640    CONTINUE
C                                                                        4570   
  645 CONTINUE
C
C     NOW START A NEW ITERATION
C
      GO TO 300
C
C     CONVERGENCE HAS BEEN PREDICTED
C
  700 DO 705 I = 1, N
         X(I) = W(NX+I)                                                  4580   
  705 CONTINUE
C
      DO 710 I = 1, M
         F(I) = W(NF+I)
  710 CONTINUE
C
  720 MAXFCN = NFCALL
  730 RETURN
C     ----- LAST CARD OF SUBROUTINE DLSMIN ----
      END                                                                4590   
C*ID* FCNCHK   C109.PTOLEMY.FITTERS                         PER705  22:
C
C     ------------------------------------------------------------------
C
      SUBROUTINE FCNCHK(N,M,X,F,FCN,IS,IERR,ACC,FSQ,FMIN,SQSTEP,
     X                  DELTSQ,MAXFCN,NFCALL,NTEST,NX,NF,W)
C
      DOUBLE PRECISION ACC,DELTSQ,F,FMIN,FSQ,G,SQSTEP,W,X
      INTEGER I,IERR,IS,M,MAXFCN,N,NF,NFCALL,NTEST,NX
      DIMENSION X(N),F(N),W(1)                                           4600   
C
C       THIS SUBROUTINE INVOKES THE USER SUPPLIED SUBROUTINE  FCN  TO
C     EVALUATE THE FUNCTIONS  F(J)(X(1),X(2),...,X(N))  J = 1,2,...,M,
C     AND CALCULATES THE SUM OF THE SQUARES OF THE  F(J).  VARIOUS TESTS
C     FOR CONVERGENCE AND FOR MALFUNCTION ARE ALSO MADE.
C
C     ------------------------------------------------------------------
C
C     CALL THE SUBROUTINE  FCN, ITERATE  FCN  COUNTER AND EXIT FROM
C       THE CALLING ROUTINE IF  IERR  HAS BEEN SET .LT. 0 BY USER        4610   
C
      CALL FCN (N,M,X,F,G,IERR)
      NFCALL = NFCALL + 1
      IF (IERR .LT. 0) GO TO 90
C
C     COMPUTE SUM OF SQUARES  FSQ  AND TEST FOR CONVERGENCE
C
      FSQ = 0.0D0
C
      DO 20 I = 1, M                                                     4620   
         FSQ = FSQ + F(I) * F(I)
   20 CONTINUE
C
      IF (FSQ .LE. ACC) GO TO 90
C
C     TEST 'IS' SWITCH
C
      IF ((IS .NE. 1) .AND.
     1    (IS .NE. 4)) GO TO 60
C                                                                        4630   
C     TEST FOR ERROR RETURN BECAUSE F(X) DOES NOT DECREASE
C
      IF (FSQ .LT. FMIN) GO TO 50
      IF (SQSTEP .GT. DELTSQ) GO TO 60
      NTEST = NTEST - 1
      IF (NTEST .GT. 0) GO TO 60
C
C     IERR = 2 IF N + 4 CALLS OF FCN FAILED TO IMPROVE RESIDUALS
C     IERR = 3 IF A NEW JACOBIAN IS UNSUCCESSFUL
C                                                                        4640   
      IERR = 2
      IF (NTEST .LT. 0) IERR = 3
      GO TO 70
   50 NTEST = N + 4
C
C     TEST WHETHER THERE HAVE BEEN  MAXFCN  CALLS OF  FCN
C
   60 IF (NFCALL .LT. MAXFCN) GO TO 100
      IERR = 1
      IF (FSQ .LT. FMIN) GO TO 90                                        4650   
C
C     OBTAIN LAST APPROXIMATION
C
   70 DO 80 I = 1, N
         X(I) = W(NX+I)
   80 CONTINUE
C
      DO 85 I = 1, M
         F(I) = W(NF+I)
   85 CONTINUE                                                           4660   
C
   90 IS = 0
  100 RETURN
C     ----- LAST CARD OF SUBROUTINE FCNCHK -----
      END
C*ID* FCUBIC   C109.PTOLEMY.FITTERS                         PER705  21:
      FUNCTION FCUBIC(F0,G0,F1,G1,XLAM)
C     CUBIC INTERPOLATION (ACTON VARIANT) FOR MINIMIZERS
C     GIVEN FUNCTION (F0,F1) AND DERIVATIVE (G0,G1) AT TWO POINTS
C     WITH XLAM=X2-X1 (X2.GT.X1) FITS A CUBIC TO THES 4 ITEMS OF         4670   
C     INFORMATION AND RETURNS MINIMUM PREDICTED BY THIS CUBIC
C     IF WSQ.LT.0 CUBIC HAS NO MINIMUM AND QUADRATIC INTERPOLATION
C     IS USED INSTEAD. IN EITHER CASE FCUBIC=X(MIN)-X1
C     NO TESTS ARE MADE TO ENSURE THAT VALUE RETURNED IS SENSIBLE
C     I.E. 0.LT.FCUBIC.LT.XLAM
C     THIS SHOULD BE DONE IN CALLING PROGRAM
C     MHM JANUARY 4 1976
C
      implicit real*8 ( a-h, o-z )                                      implicit
C                                                                        4680   
      F=(F1-F0)/XLAM
      G=G1-G0
      C=G-2*(F-G0)
      D=G-3*C
      WSQ=D*D-12*C*G0
      IF(WSQ.LT.0) GO TO 1
      W=DSQRT(WSQ)
      Y=-2*G0/(D+W)
      GO TO 2
    1 B=G1-F                                                             4690   
      Y=.5*(1-F/B)
    2 FCUBIC=Y*XLAM
      RETURN
      END
C*ID* GENINV   C109.PTOLEMY.FITTERS                         PER705  22:
C
C     -------------------------------------------------------------
C
      SUBROUTINE GENINV (M,N,A,IA,W)
C                                                                        4700   
      DOUBLE PRECISION A,BSQ,RMAX,SIGMA,SUM,W
      INTEGER I,IA,IR,J,K,KK,KP,L,M,N,NC,NR
      DIMENSION A(IA,1),W(1)
C
C     -------------------------------------------------------------
C
C     PARTITION THE SCRATCH ARRAY  W
C
      NR = M
      NC = M + M                                                         4710   
C
C     THE SCRATCH ARRAY  W  HAS BEEN PARTITIONED AS FOLLOWS
C
C        W(1),...,W(M)       - CONTAINS THE FIRST COMPONENTS OF THE
C                              VECTORS OF THE ELEMENTARY TRANSFORMATIONS
C        W(NR+1),...,W(NR+M) - CONTAINS THE RECORD OF ROW INTERCHANGES
C        W(NC+1),...,W(NC+N) - CONTAINS THE RECORD OF COLUMN
C                              INTERCHANGES
C
C                                                                        4720   
C     SET THE INITIAL RECORDS OF ROW AND COLUMN INTERCHANGES
C
      DO 105 I = 1, M
         W(NR+I) = 0.5D0 + DFLOAT(I)
  105 CONTINUE
C
      DO 110 I = 1, N
         W(NC+I) = 0.5D0 + DFLOAT(I)
  110 CONTINUE
C                                                                        4730   
C     'KK' COUNTS THE SEPARATE ELEMENTARY TRANSFORMATIONS
C
      KK = 1
C
C     FIND LARGEST ROW AND MAKE ROW INTERCHANGES
C
  115 RMAX = 0.0D0
C
      DO 130 I = KK, M
         SUM = 0.0D0                                                     4740   
C
         DO 120 J = KK, N
            SUM = SUM + A(I,J) ** 2
  120    CONTINUE
C
         IF (RMAX .GE. SUM) GO TO 130
         RMAX = SUM
         IR = I
  130 CONTINUE
C                                                                        4750   
      IF (IR .LE. KK) GO TO 140
      K = NR + KK
      L = NR + IR
      SUM = W(K)
      W(K) = W(L)
      W(L) = SUM
C
      DO 135 J = 1, N
         SUM = A(KK,J)
         A(KK,J) = A(IR,J)                                               4760   
         A(IR,J) = SUM
  135 CONTINUE
C
C     FIND LARGEST ELEMENT OF PIVOTAL ROW, AND MAKE COLUMN INTERCHANGES
C
  140 RMAX = 0.0D0
      SUM = 0.0D0
C
      DO 150 J = KK, N
         SUM = SUM + A(KK,J) ** 2                                        4770   
         IF (RMAX .GE. DABS(A(KK,J))) GO TO 150
         RMAX = DABS(A(KK,J))
         IR = J
  150 CONTINUE
C
      IF (IR .LE. KK) GO TO 160
      K = NC + KK
      L = NC + IR
      RMAX = W(K)
      W(K) = W(L)                                                        4780   
      W(L) = RMAX
C
      DO 155 I = 1, M
         RMAX = A(I,KK)
         A(I,KK) = A(I,IR)
         A(I,IR) = RMAX
  155 CONTINUE
C
C     BEGIN THE ORTHOGONAL TRIANGULARIZATION PROCESS BY REPLACING THE
C       THE PIVOTAL ROW BY THE VECTOR OF THE TRANSFORMATION              4790   
C
  160 SIGMA = DSQRT(SUM)
      BSQ = DSQRT(SUM + SIGMA * DABS(A(KK,KK)))
      W(KK) = DSIGN(SIGMA + DABS(A(KK,KK)),A(KK,KK)) / BSQ
      A(KK,KK) = - DSIGN(SIGMA,A(KK,KK))
      KP = KK + 1
      IF (KP .GT. N) GO TO 200
C
      DO 165 J = KP, N
         A(KK,J) = A(KK,J) / BSQ                                         4800   
  165 CONTINUE
C
C     APPLY THE TRANSFORMATION TO THE REMAINING ROWS OF A
C
      IF (KP .GT. M) GO TO 200
C
      DO 180 I = KP, M
         SUM = W(KK) * A(I,KK)
C
         DO 170 J = KP, N                                                4810   
            SUM = SUM + A(KK,J) * A(I,J)
  170    CONTINUE
C
         A(I,KK) = A(I,KK) - SUM * W(KK)
C
         DO 175 J = KP, N
            A(I,J) = A(I,J) - SUM * A(KK,J)
  175    CONTINUE
C
  180 CONTINUE                                                           4820   
C
      KK = KP
      GO TO 115
C
C     THE ORTHOGONAL TRIANGULARIZATION OF  A  IS COMPLETE. NOW BEGIN
C       THE CONSTRUCTION OF THE GENERALIZED INVERSE BY APPLYING THE
C       FIRST ELEMENTARY TRANSFORMATION
C
  200 KK = M
      KP = M + 1                                                         4830   
      SUM = W(M) / A(M,M)
      IF (N .LE. M) GO TO 210
C
      DO 205 J = KP, N
         A(M,J) = - SUM * A(M,J)
  205 CONTINUE
C
  210 A(M,M) = 1.0D0 / A(M,M) - SUM * W(M)
C
C     NOW APPLY THE OTHER (M - 1) TRANSFORMATIONS                        4840   
C
  215 KP = KK
      KK = KP - 1
      IF (KK .LE. 0) GO TO 250
C
C     FIRST TRANSFORM THE LAST (M - KK) ROWS
C
      DO 230 I = KP, M
         SUM = 0.0D0
C                                                                        4850   
         DO 220 J = KP, N
            SUM = SUM + A(KK,J) * A(I,J)
  220    CONTINUE
C
         DO 225 J = KP, N
            A(I,J) = A(I,J) - SUM * A(KK,J)
  225    CONTINUE
C
         W(I) = - SUM * W(KK)
  230 CONTINUE                                                           4860   
C
C     THEN CALCULATE THE NEW ROW IN POSITION KK
C
      DO 240 J = KP, N
         SUM = - W(KK) * A(KK,J)
C
         DO 235 I = KP, M
            SUM = SUM - A(I,KK) * A(I,J)
  235    CONTINUE
C                                                                        4870   
         A(KK,J) = SUM / A(KK,KK)
  240 CONTINUE
C
C     AND REVISE THE COLUMN IN POSITION KK
C
      SUM = 1.0D0 - W(KK) ** 2
C
      DO 245 I = KP, M
         SUM = SUM - A(I,KK) * W(I)
         A(I,KK) = W(I)                                                  4880   
  245 CONTINUE
C
      A(KK,KK) = SUM / A(KK,KK)
      GO TO 215
C
C     COMPLETE THE CONSTRUCTION OF THE GENERALIZED INVERSE BY
C       RESTORING THE ROW AND COLUMN INTERCHANGES
C
  250 DO 265 I = 1, M
  255    K = NR + I                                                      4890   
         IR = W(K)
         IF (I .GE. IR) GO TO 265
         L = NR + IR
         SUM = W(K)
         W(K) = W(L)
         W(L) = SUM
C
      DO 260 J = 1, N
            SUM = A(I,J)
            A(I,J) = A(IR,J)                                             4900   
            A(IR,J) = SUM
  260    CONTINUE
C
      GO TO 255
  265 CONTINUE
      DO 280 J = 1, N
  270    K = NC + J
         IR = W(K)
         IF (J .GE. IR) GO TO 280
         L = NC + IR                                                     4910   
         SUM = W(K)
         W(K) = W(L)
         W(L) = SUM
C
         DO 275 I = 1, M
            SUM = A(I,J)
            A(I,J) = A(I,IR)
            A(I,IR) = SUM
  275    CONTINUE
C                                                                        4920   
         GO TO 270
  280 CONTINUE
C
      RETURN
C     ----- LAST CARD OF SUBROUTINE  GENINV ----
      END
C*ID* LMCHOL   C109.PTOLEMY.FITTERS                         PER705  11:
C
C     ------------------------------------------------------------------
C                                                                        4930   
      SUBROUTINE LMCHOL(NVAR,NFUNS,X,MODE,MAXF,S,EVLFUN, TOL, TOL2,
     1                  FVECT,RJACOB,RJINV,IERR,NFCALL,NJCALL,SCRV1,
     1                  SQSCAL,DJTJ,VJTF,SCRV2,HISTRY, MONIT, IDEB )
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,S,FVECT,RJACOB,RJINV,SCRV1,SQSCAL,DJTJ,VJTF,SCRV2
C
      INTEGER I,ICONV,IERR,MAXF,MODE,NFCALL,NFUNS,NJCALL,NVAR
      REAL*8  X(NVAR),S(NVAR),FVECT(NFUNS),RJACOB(NFUNS,NVAR),
     1        RJINV(NVAR,NVAR),SCRV1(NVAR),SQSCAL(NVAR),DJTJ(NVAR),      4940   
     2        VJTF(NVAR),SCRV2(NFUNS)
      REAL*8  TOL,HISTRY
      REAL*8  ACTRED,DABS,DIRDRV,DSQRT,FSQ,FSQP,PRERED,QCPARM,QPARM,Z
C
C     LEAST SQUARES MINIMIZATION WITH ANALYTIC DERIVATIVES
C
C     THIS ROUTINE HAS BEEN SLIGHTLY MODIFIED BY S. PIEPER:
C       1)  EVLFUN AND EVLJAC HAVE BEEN COMBINED
C       2)  TOL2 PARAMETER ADDED
C       3)  MONIT AND IDEB ADDED                                         4950   
C       4)  IMPLICIT STATEMENT ADDED.
C       5)  A STOP SEARCH BASED ON SMALL CHANGES IN THE SUM OF SQUARES
C           WAS ADDED.
C       6)  RJINV IS SET TO AN ESTIMATE OF THE COVARIANCE MATRIX AT
C           THE END OF THE SEARCH.
C       7)  THE TEST FOR CONVERGENCE DUE TO SMALL STEPS IS MADE
C           BEFORE EACH EVLFUN CALL.  THIS SAVES AN UNNECESSARY
C           EVLFUN CALL AFTER CONVERGENCE HAS BEEN ACHIEVED.
C           ALSO THE FINAL JACOBIAN IS ALWAYS REQUESTED AT THE END
C           OF AN ITERATION.                                             4960   
C       8)  IF EVLFUN RETURNS  IERR = 1,  A STEP OUT-OF-BOUNDS IS
C           IS ASSUMED; THE MARQUARDT PARAMETER IS MULTIPLIED BY 4
C           AND A NEW STEP IS TRIED.
C
C
C       THIS ROUTINE IS A MODULARIZED MODIFICATION OF THE COMPUTER
C     ALGORITHM  VA07A  BY R. FLETCHER IN HARWELL REPORT  AERE - R.6799.
C
C       THIS SUBROUTINE MINIMIZES A SUM OF SQUARES  FSQ(X)  OF
C     FUNCTIONS  F(J)(X(1),X(2),...,X(NVAR)), J = 1,2,...,NFUNS, USING   4970   
C     A MODIFICATION OF THE LEVENBERG-MARQUARDT ALGORITHM. THE USER
C     MUST PROVIDE A SUBROUTINE TO CALCULATE THE FUNCTIONS
C     F(J)  AND THE ELEMENTS OF THE JACOBIAN
C     (FIRST PARTIAL DERIVATIVE) MATRIX.
C
C     TERMINOLOGY:
C
C         QUANTITY TO BE MINIMIZED:
C            FSQ  =  SUM( 1 =< J =< NFUNS )  F(J)**2
C                                                                        4980   
C         JACOBIAN:
C            JACOBIAN(J,N)  =  D(F(J)) / D(X(N)) ,
C              1 =< J =< NFUNS ,
C              1 =< N =< NVAR .
C
C
C     ******************************************************************
C       THIS ROUTINE INVOKES THE  MINPACK  PACKAGE MODULES  CMTRIX,
C     QMTRIX, UPQPAR  AND  SIP  AND THE USER SUPPLIED SUBROUTINES
C     EVLFUN  AND  MONIT. THE  ANL  LIBRARY ROUTINE  ANL F152S, SYMINV   4990   
C     IS ALSO INVOKED.
C     ******************************************************************
C
C     ON INPUT
C
C        NVAR  IS THE NUMBER OF VARIABLES. IT IS ALSO THE DIMENSION OF
C           THE VECTORS  X, S, SCRV1, SQSCQL, DJTJ  AND  VJTF  THE
C           SQUARE ARRAY  RJINV  AND IS THE COLUMN DIMENSION OF THE
C           ARRAY  RJACOB.
C                                                                        5000   
C        NFUNS  IS THE NUMBER OF FUNCTIONS  F(J). IT IS ALSO THE
C           DIMENSION OF THE VECTOR  FVECT  AND  SCRV2  AND IS THE ROW
C           DIMENSION OF THE ARRAY  RJACOB. A DIAGNOSTIC IS RETURNED FOR
C           NFUNS .LT. NVAR.
C
C        X  CONTAINS AN ESTIMATE OF THE SOLUTION VECTOR
C           (X(1),X(2),...,X(NVAR)).
C
C        MODE  IS A PARAMETER SET EQUAL TO  K  WHICH GOVERNS THE
C           METHOD OF SCALING THE VARIABLES THROUGH THE INPUT PARAMETER  5010   
C           S. THE OPTIONS ARE
C
C           K = 1 - THE  S(I)  ARE SET BY  LMCHOL  TO THE DIAGONAL
C                   ELEMENTS OF  (JACOBIAN TRANSPOSE) * JACOBIAN, IF
C                   POSITIVE, AND TO ONE OTHERWISE,
C           K = 2 - THE  S(I)  ARE SET BY  LMCHOL  TO ONE FOR ALL  I,
C           K = 3 - THE  S(I)  ARE SET BY THE USER THROUGH THE
C                   PARAMETER  S,
C
C           IF  MODE  IS SET EQUAL TO ANY OTHER VALUE IT IS RESET TO     5020   
C           ONE BY  LMCHOL.
C
C        MAXF  IS THE LIMIT ON THE NUMBER OF CALLS TO  EVLFUN.
C
C        S  CONTAINS THE MULTIPLIERS WHICH CONTROL SCALING OF THE
C           VARIABLES IF  MODE  IS SET TO 3. OTHERWISE IT IS
C           INITIALIZED INTERNALLY.
C
C        EVLFUN  IS A USER SUPPLIED SUBROUTINE TO CALCULATE FUNCTIONS
C           F(J)  IN  FVECT  AT THE ESTIMATE  X(I),I = 1,2,...,NVAR. THE 5030   
C           USER MAY FORCE AN EXIT FROM  LMCHOL  BY SETTING  IERR  TO A
C           NEGATIVE INTEGER IN  EVLFUN. EVLFUN, MUST BE DECLARED IN AN
C           EXTERNAL  STATEMENT IN THE USERS CALLING PROGRAM, AND
C           SHOULD BE WRITTEN AS FOLLOWS
C
C           SUBROUTINE EVLFUN ( ITYPE, NVAR,NFUNS,X,FVECT, RJACOB,IERR)
C           REAL*8  X(NVAR),FVECT(NFUNS)
C           REAL*8  RJACOB(NFUNS,NVAR)
C           INTEGER  NVAR,NFUNS,IERR
C           -----------------------------                                5040   
C
C       STATEMENTS TO EVALUATE
C         1)  FVECT = F(J)       IF ITYPE = 1 OR 3.
C         2)  RJACOB = JACOBIAN   IF ITYPE = 2 OR 3.
C         3)  SET IERR IF DESIRED TO STOP SEARCH.
C     OR
C         1)  SET IERR = +1 TO INDICATE A STEP OUT OF BOUNDS.
C
C     IF THE USER FINDS IT MORE CONVIENENT (OR OF LITTLE EXTRA COST),
C     TO COMPUTE BOTH QUANTITIES WHEN ITYPE IS 1 OR 2, HE MAY DO SO AND  5050   
C     SHOULD THEN SET ITYPE TO 3.  IN MANY CASES, COMPUTING THE GRADIENT
C     WHEN ONLY THE FUNCTION IS REQUESTED WILL SAVE A SUBSEQUENT REQUEST
C     FOR THE GRADIENT AT THE SAME PLACE ( HOWEVER, SOMETIMES THE
C     GRADIENT WILL NEVER BE USED).
C
C           -----------------------------
C           RETURN
C           END
C
C                                                                        5060   
C        TOL & TOL2  IS THE ACCURACY REQUIRED IN THE SOLUTION, I.E. A
C           NORMAL RETURN FROM THE ROUTINE OCCURS WHEN THE DIFFERENCES
C           BETWEEN THE COMPONENTS OF TWO SUCCESSIVE ESTIMATES OF THE
C           SOLUTION ARE NOT GREATER THAN  MAX(TOL*ABS(X(I)),TOL2)  FOR
C           ALL  I. NOTE THAT A NORMAL RETURN WILL ALSO OCCUR IF THE
C           PREDICTED REDUCTION IN THE SUM OF SQUARES IS .LE. 0.
C          NORMAL RETURN IS ALSO MADE IF
C             FSQ(I-2) - FSQ(I)  <  TOL*FSQ(I)
C          OR IF
C             FSQ(I) < TOL2 .                                            5070   
C          THIS CRITERION IS DESINGED TO PREVENT MINIMIZATION BASED ON
C          INSIGNIFICANT (OR POSSIBLY NOISY) CHANGES IN THE SUM OF THE
C          SQUARES.  NOTE THE COMPARASON OF STEP I WITH THAT FOR
C          I-2 TO ALLOW A MOMENTARY PAUSE IN THE FUNCTION IMPROVEMENT.
C
C        HISTRY  IS A DUMMY PARAMETER.
C
C     MONIT IS A SUBROUTINE THAT WILL BE CALLED AT THE END OF
C        EACH ITERATION IF IDEB IS 1.  SEE BELOW FOR ITS FORM.
C                                                                        5080   
C     IDEB = 0  MEANS DON'T EVER CALL MONIT.
C          = 1  MEANS CALL MONIT AT THE END OF EACH ITERATION.
C
C
C     ON OUTPUT
C
C        X  CONTAINS THE BEST AVAILABLE ESTIMATE OF THE SOLUTION VECTOR.
C
C        FVECT  CONTAINS THE FUNCTION VALUES  F(J)  CORRESPONDING TO
C           THE  X  VECTOR WITH THE EXCEPTION THAT THE CONTENTS ARE      5090   
C           MEANINGLESS WHEN  IERR  .EQ. 10.
C
C        RJACOB  CONTAINS THE ELEMENTS OF THE JACOBIAN MATRIX
C           CORRESPONDING TO THE  X  VECTOR WITH THE EXCEPTION THAT THE
C           CONTENTS ARE MEANINGLESS WHEN  IERR .EQ. 10.
C
C        RJINV  CONTAINS AN ESTIMATE OF THE NVAR X NVAR  COVARIANCE
C            MATRIX.  THE ESTIMATE IS JUST
C               INVERSE( (JACOBIAN TRANSPOSE) * JACOBIAN )
C                                                                        5100   
C        IERR  IS A PARAMETER SET EQUAL TO  K  WHICH GIVES THE FOLLOWING
C           TERMINATION INDICATIONS
C
C           NORMAL TERMINATION,
C             K = 0,
C           INTERMEDIATE TERMINATION,
C             K = -N(N ANY INTEGER) - USER TERMINATION,
C             K =  1 - FAILURE TO CONVERGE IN  MAXF  CALLS OF  EVLFUN,
C             K =  2 - DIRECTIONAL DERIVATIVE IS .GE. 0,
C           INITIAL TERMINATION,                                         5110   
C             K = 10 - NFUNS .LT. NVAR.
C
C        NFCALL IS  THE NUMBER OF CALLS TO  EVLFUN.
C
C        NJCALL IS  THE NUMBER OF CALLS TO  EVLJAC.
C
C        SCRV1  IS A SCRATCH VECTOR.
C
C        SQSCAL  CONTAINS THE SQUARE ROOTS OF THE ELEMENTS OF THE
C           SCALE CONTROL VECTOR  S.                                     5120   
C
C        DJTJ  CONTAINS THE DIAGONAL ELEMENTS OF
C           (JACOBIAN TRANSPOSE) * JACOBIAN.
C
C        VJTF  CONTAINS THE VECTOR  (JACOBIAN TRANSPOSE) * FVECT.
C
C        SCRV2  IS A SCRATCH VECTOR.
C
C     WRITTEN BY  K. E. HILLSTROM, SEPTEMBER, 1975.
C     10/31/75 - SLIGHT MODIFICATIONS - S.P.                             5130   
C     1/1/76   - COMPUTE THE COVARIANCE MATRIX IN RJINV - S.P.
C     2/6/76 - CHECK FOR STEP < TOL BEFORE EVLFUN - S. P.
C     11/25/76 - ALLOW STEPS OUT OF BOUNDS - S.P.
C     2/7/77 - DO NOT STOP IMMEADIATELY AFTER AN OUT-OF-BOUNDS
C
C     ------------------------------------------------------------------
C
C     TEST THE INPUT PARAMETERS  NFUNS  AND  NVAR  FOR CONSISTENCY
C
  100 IERR = 0                                                           5140   
      IF (NFUNS .LT. NVAR) IERR = 10
      IF (IERR .NE. 0) GO TO 325
      IF (MODE .LT. 1 .OR. MODE .GT. 3) MODE = 1
C
C     INITIALIZE THE FOLLOWING PARAMETERS
C
C     QPARM  - THE MARQUARDT ITERATION CONTROL PARAMETER
C     QCPARM - THE MARQUARDT PARAMETER CUT-OFF VALUE
C     NFCALL - THE NUMBER OF CALLS TO  EVLFUN
C     NJCALL - THE NUMBER OF CALLS TO  EVLJAC                            5150   
C
      QPARM = 0.0D0
      QCPARM = 1.0D0
      NFCALL = 1
      NJCALL = 1
      ITCNT = -1
C
C     PREPARE FOR THE FIRST ITERATION BY INITIALIZING THE FUNCTION
C       VALUES , THE SUM OF SQUARES  FSQ  AND THE JACOBIAN MATRIX
C                                                                        5160   
C     ******************************************************************
      ITYPE = 3
      CALL EVLFUN ( ITYPE, NVAR,NFUNS,X,FVECT, RJACOB, IERR )
C     ******************************************************************
C
      IF (IERR .LT. 0) GO TO 325
C
CXXXX C     ************************************************************
CXXXX       CALL EVLJAC(NVAR,NFUNS,X,RJACOB,IERR)
CXXXX C     ************************************************************ 5170   
CXXXX C
CXXXX       IF (IERR .LT. 0) GO TO 325
C
C     ******************************************************************
      CALL SIP(FVECT,1,2,FVECT,1,2,0.0D0,FSQ,NFUNS)
C     ******************************************************************
C
C     CALL  CMTRIX  TO FORM THE MATRIX
C
C       JTJ  = (JACOBIAN TRANSPOSE) * JACOBIAN                           5180   
C       AND THE VECTOR
C       JTF  = (JACOBIAN TRANSPOSE) * FVECT.
C
C       IN THIS SECTION
C
C       RJINV  - CONTAINS THE LOWER TRIANGLE OF THE SYMMETRIC MATRIX
C                JTJ
C
C       ALSO, FROM THIS SECTION ON
C                                                                        5190   
C       VJTF(I)  - CONTAINS THE VECTOR  JTF
C
C     ******************************************************************
      CALL CMTRIX(NVAR,NFUNS,RJACOB,FVECT,RJINV,VJTF)
C     ******************************************************************
C
C     INITIALIZE THE SCALE CONTROL VECTOR  S  ACCORDING TO THE VALUE
C       OF  MODE. FROM THIS SECTION ON
C
C     SQSCAL(I) - CONTAINS THE SQUARE ROOTS OF THE ELEMENTS OF  S        5200   
C
      IF (MODE .EQ. 3) GO TO 110
      IF (MODE .EQ. 2) GO TO 120
C
C     MODE = 1, SET  S(I) = JTJ(I,I) OR 1
C
      DO 105 I = 1, NVAR
         S(I) = RJINV(I,I)
         IF (S(I) .LE. 0.0D0) S(I) = 1.0D0
  105 CONTINUE                                                           5210   
C
C     MODE = 3, S(I) HAS BEEN SET BY THE USER THROUGH THE PARAMETER  S
C
  110 DO 115 I = 1, NVAR
         SQSCAL(I) = DSQRT(S(I))
  115 CONTINUE
C
      GO TO 200
C
C     MODE = 2 SET S(I) = 1                                              5220   
C
  120 DO 125 I = 1, NVAR
         SQSCAL(I) = 1.0D0
         S(I) = 1.0D0
  125 CONTINUE
C
C     BEGIN AN ITERATION BY PLACING DIAGONAL ELEMENTS OF  JTJ  IN
C       DJTJ(I)
C
C     ALLOW USER TO PRINT STUFF EACH ITERATION                           5230   
C
 200  ITCNT = ITCNT + 1
      IF ( IDEB .EQ. 1 )  CALL MONIT ( ITCNT, NFCALL, NJCALL, NVAR,
     1  NFUNS, X, FSQ, FVECT )
      NFLAST = NFCALL
C
      DO 205 I = 1, NVAR
         DJTJ(I) = RJINV(I,I)
  205 CONTINUE
C                                                                        5240   
C     CALL  QMTRIX  TO SOLVE THE MARQUARDT MATRIX EQUATION
C
C       (JTJ + QPARM * S) * DELTA = - JTF.
C
C       NOTE THAT  QMTRIX  INVOKES THE  ANL  LIBRARY ROUTINE
C       ANL F152S, SYMINV  TO SOLVE THE MATRIX EQUATION.
C
C       WHEN FINISHED
C
C       RJINV    - CONTAINS THE UPPER TRIANGLE OF                        5250   
C                  ((JTJ + QPARM * S) INVERSE)  AND THE STRICT LOWER
C                  TRIANGLE OF  JTJ
C       SCRV2(I) - CONTAINS THE NEGATIVE OF THE CORRECTION VECTOR
C                  DELTA
C
C
C     ******************************************************************
  210 CALL QMTRIX(NVAR,RJINV,DJTJ,VJTF,S,SCRV2,QPARM)
C     ******************************************************************
C                                                                        5260   
C     ERROR RETURN IF THE DIRECTIONAL DERIVATIVE ALONG THE CORRECTION,
C       (JTF TRANSPOSE) * DELTA, IS .GE. 0.
C
C     ******************************************************************
      CALL SIP(VJTF,1,2,SCRV2,1,2,0.0D0,DIRDRV,NVAR)
C     ******************************************************************
C
      IF (DIRDRV .LE. 0.0D0) GO TO 290
C
C     COMPUTE  PRERED  THE PREDICTED REDUCTION IN THE SUM OF SQUARES     5270   
C
      DO 230 I = 1, NVAR
         Z = DJTJ(I) * SCRV2(I)
C
C     ******************************************************************
         IF (I .GT. 1) CALL SIP(RJINV(I,1),1,NVAR+1,SCRV2,1,2,Z,Z,
     1   I-1)
C     ******************************************************************
         IF (I .LT. NVAR) CALL SIP(RJINV(I,I),2,3,SCRV2(I),2,3,Z,Z,
     1   NVAR-I)                                                         5280   
C     ******************************************************************
C
         SCRV1(I) = 2.0D0 * VJTF(I) - Z
  230 CONTINUE
C
C     ******************************************************************
      CALL SIP(SCRV1,1,2,SCRV2,1,2,0.0D0,PRERED,NVAR)
C     ******************************************************************
C
C     ICONV  = 0 IF  DABS(DELTA) .LE. MAX(TOL*DABS(X(I)),TOL)  FOR ALL   5290   
C       I. FROM THIS SECTION ON
C
C       SCRV1(I)  -  CONTAINS THE ESTIMATE (XK + DELTA)
C
      ICONV = 0
C
      DO 235 I = 1, NVAR
         Z = TOL * DABS(X(I))
         IF (Z .LE. TOL2) Z = TOL2
         IF (DABS(SCRV2(I)) .GT. Z) ICONV = 1                            5300   
         SCRV1(I) = X(I) - SCRV2(I)
  235 CONTINUE
C
C     NORMAL EXIT IF  ICONV = 0  OR THE PREDICTED REDUCTION IS .LE. 0.
C
C     DO NOT STOP IF THE LAST STEP WAS OUT OF BOUNDS
C
      IF ( IERR .EQ. 1 )  GO TO 240
      IF ((ICONV .EQ. 0) .OR.
     1    (PRERED .LE. 0.0D0)) GO TO 310                                 5310   
C
C     TEST WHETHER THERE HAVE BEEN  MAXFCN  CALLS OF  EVLFUN
C
 240  IF (NFCALL .GE. MAXF) GO TO 295
C
C     COMPUTE  FSQP  THE SUM OF SQUARES AT THE TRIAL ESTIMATE. FROM
C       THIS SECTION ON
C
C       SCRV2(I)  -  CONTAINS THE FUNCTIONS AT THE TRIAL ESTIMATE
C                                                                        5320   
C     ******************************************************************
      ITYPE = 1
      CALL EVLFUN ( ITYPE, NVAR,NFUNS,SCRV1,SCRV2, RJACOB, IERR )
C     ******************************************************************
C
      NFCALL = NFCALL + 1
      IF (IERR .LT. 0) GO TO 310
      IF ( IERR .NE. 1 )  GO TO 260
C
C     STEP WAS OUT OF BOUNDS                                             5330   
C
      QPARM = 4*QPARM
      IF ( QPARM .LT. .01 )  QPARM = .01
      QCPARM = .01*QPARM
      GO TO 210
C
C     ******************************************************************
 260  CALL SIP(SCRV2,1,2,SCRV2,1,2,0.0D0,FSQP,NFUNS)
C     ******************************************************************
C                                                                        5340   
C     COMPUTE  ACTRED  THE ACTUAL REDUCTION IN THE SUM OF SQUARES AND
C       TEST WHETHER TRIAL SUM OF SQUARES IS LEAST
C
      ACTRED = FSQ - FSQP
C
C     CALL  UPQPAR  TO UPDATE THE MARQUARDT PARAMETER
C
C     ******************************************************************
      CALL UPQPAR(NVAR,SQSCAL,RJINV,ACTRED,PRERED,DIRDRV,QPARM,
     1            QCPARM)                                                5350   
C     ******************************************************************
C
C     TEST WHETHER TRIAL SUM OF SQUARES IS LEAST
C
      IF (ACTRED .LE. 0.0D0) GO TO 210
C
C     THE TRIAL SUM OF SQUARES IS LEAST. SET FSQ,  X  AND FVECT  TO
C       BEST CALCULATED VALUES.
C
      DO 275 I = 1, NVAR                                                 5360   
         X(I) = SCRV1(I)
  275 CONTINUE
C
      DO 280 I = 1, NFUNS
         FVECT(I) = SCRV2(I)
  280 CONTINUE
C
C     THE USER MAY HAVE ELECTED TO COMPUTE THE JACOBIAN ON THE PREVIOUS
C     FUNCTION REQUEST.
C                                                                        5370   
      IF ( ITYPE .EQ. 3 ) GO TO 285
C
C     ******************************************************************
      ITYPE = 2
      CALL EVLFUN ( ITYPE, NVAR,NFUNS,X, FVECT, RJACOB, IERR )
C     ******************************************************************
C
 285  NJCALL = NJCALL + 1
      IF (IERR .LT. 0) GO TO 310
C                                                                        5380   
C     SHOULD WE STOP DUE TO INSIGNIFICANT CHANGES IN FSQ.
C
      IF ( ITCNT .LE. 2 )  GO TO 270
      IF ( FSQP .LT. TOL2 )  GO TO 300
      IF ( FSQM1-FSQP .LT. TOL*FSQP )  GO TO 300
C
C     KEEP GOING
C
 270  FSQM1 = FSQ
      FSQ = FSQP                                                         5390   
C
C     CALL  CMTRIX  TO FORM NEW  JTJ  AND  JTF  AND START A NEW
C       ITERATION
C
C     ******************************************************************
      CALL CMTRIX(NVAR,NFUNS,RJACOB,FVECT,RJINV,VJTF)
C     ******************************************************************
C
      GO TO 200
C                                                                        5400   
C     ERROR RETURN BECAUSE DIRECTIONAL DERIVATIVE IS .GE. 0
C
  290 IERR = 2
      GO TO 325
C
C     ERROR RETURN BECAUSE THERE HAVE BEEN  MAXF  CALLS OF  EVLFUN
C
  295 IERR = 1
      GO TO 310
C                                                                        5410   
C     CONVERGENCE CRITERION HAS BEEN SATISFIED OR THE PREDICTED
C       REDUCTION IN THE SUM OF SQUARES IS .LE. 0 OR THE FUNCTION CALL
C       LIMIT HAS BEEN REACHED. X  AND  FVECT ARE UPDATED ONLY IF THE
C       FINAL REDUCTION IN  FSQ  IS POSITIVE
C
  300 ITCNT = ITCNT + 1
C
C     ALLOW USER TO PRINT STUFF EACH ITERATION
C
 310  IF ( IDEB .EQ. 1  .AND.  NFCALL .GT. NFLAST )                      5420   
     1  CALL MONIT ( ITCNT, NFCALL, NJCALL, NVAR,
     2    NFUNS, X, FSQP, FVECT )
C
C     COMPUTE  INVERSE(JTJ)  AS AN ESTIMATE OF THE COVARIANCE MATRIX
C     AT THE FINAL VALUE.
C
 325  CALL QMTRIX ( NVAR, RJINV, DJTJ, VJTF, S, SCRV2, 0.D0 )
C
      DO 339  I = 1, NVAR
         DO 339  J = 1, I                                                5430   
            RJINV(I,J) = RJINV(J,I)
 339  CONTINUE
C
      RETURN
C
C     ----- LAST CARD OF SUBROUTINE  LMCHOL  -----
      END
C*ID* MINFIT   C109.PTOLEMY.FITTERS                         PER705  10:
C
C     ------------------------------------------------------------------ 5440   
C
      SUBROUTINE MINFIT(NM,M,N,A,W,IP,B,IERR,RV1)
C
CDC      ALLOCLEVEL A,W,B,RV1
      INTEGER I,J,K,L,M,N,II,IP,I1,KK,K1,LL,L1,M1,NM,ITS,IERR
      REAL*8 A(NM,N),W(N),B(NM,IP),RV1(N)
      REAL*8 C,F,G,H,S,X,Y,Z,EPS,SCALE,MACHEP
      REAL*8 DSQRT,DMAX1,DABS,DSIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE MINFIT,    5450   
C     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     THIS SUBROUTINE DETERMINES, TOWARDS THE SOLUTION OF THE LINEAR
C                                                        T
C     SYSTEM AX=B, THE SINGULAR VALUE DECOMPOSITION A=USV  OF A REAL
C                                         T
C     M BY N RECTANGULAR MATRIX, FORMING U B RATHER THAN U.  HOUSEHOLDER
C     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
C                                                                        5460   
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
C          AS LARGE AS THE MAXIMUM OF M AND N;
C
C        M IS THE NUMBER OF ROWS OF A AND B;
C
C        N IS THE NUMBER OF COLUMNS OF A AND THE ORDER OF V;             5470   
C
C        A CONTAINS THE RECTANGULAR COEFFICIENT MATRIX OF THE SYSTEM;
C
C        IP IS THE NUMBER OF COLUMNS OF B.  IP CAN BE ZERO;
C
C        B CONTAINS THE CONSTANT COLUMN MATRIX OF THE SYSTEM
C          IF IP IS NOT ZERO.  OTHERWISE B IS NOT REFERENCED.
C
C     ON OUTPUT:
C                                                                        5480   
C        A HAS BEEN OVERWRITTEN BY THE MATRIX V (ORTHOGONAL) OF THE
C          DECOMPOSITION IN ITS FIRST N ROWS AND COLUMNS.  IF AN
C          ERROR EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO
C          INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT;
C
C        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
C          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
C          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,IERR+2,...,N;
C                                                                        5490   
C                                   T
C        B HAS BEEN OVERWRITTEN BY U B.  IF AN ERROR EXIT IS MADE,
C                       T
C          THE ROWS OF U B CORRESPONDING TO INDICES OF CORRECT
C          SINGULAR VALUES SHOULD BE CORRECT;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS;                    5500   
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C
C     12/28/79 - STANDARD MORTRAN - RPG
C
C     ------------------------------------------------------------------ 5510   
C
C
C     ;;;;;; MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360
CIB      DATA MACHEP/Z3410000000000000/
C
      DATA MACHEP / 1.E-14 /                                            cni
C                                                                        5520   
C
C
      IERR = 0
C
C           ;;;;;; HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM
C
      G = 0.0D0
      SCALE = 0.0D0
      X = 0.0D0
C                                                                        5530   
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0D0
         S = 0.0D0
         SCALE = 0.0D0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + DABS(A(K,I))                                    5540   
C
         IF (SCALE .EQ. 0.0D0) GO TO 210
C
         DO 130 K = I, M
            A(K,I) = A(K,I) / SCALE
            S = S + A(K,I)**2
  130    CONTINUE
C
         F = A(I,I)
         G = -DSIGN(DSQRT(S),F)                                          5550   
         H = F * G - S
         A(I,I) = F - G
         IF (I .EQ. N) GO TO 160
C
         DO 150 J = L, N
            S = 0.0D0
C
            DO 140 K = I, M
  140       S = S + A(K,I) * A(K,J)
C                                                                        5560   
            F = S / H
C
            DO 150 K = I, M
               A(K,J) = A(K,J) + F * A(K,I)
  150    CONTINUE
C
  160    IF (IP .EQ. 0) GO TO 190
C
         DO 180 J = 1, IP
            S = 0.0D0                                                    5570   
C
            DO 170 K = I, M
  170       S = S + A(K,I) * B(K,J)
C
            F = S / H
C
            DO 180 K = I, M
               B(K,J) = B(K,J) + F * A(K,I)
  180    CONTINUE
C                                                                        5580   
  190    DO 200 K = I, M
  200    A(K,I) = SCALE * A(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0D0
         S = 0.0D0
         SCALE = 0.0D0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N                                                 5590   
  220    SCALE = SCALE + DABS(A(I,K))
C
         IF (SCALE .EQ. 0.0D0) GO TO 290
C
         DO 230 K = L, N
            A(I,K) = A(I,K) / SCALE
            S = S + A(I,K)**2
  230    CONTINUE
C
         F = A(I,L)                                                      5600   
         G = -DSIGN(DSQRT(S),F)
         H = F * G - S
         A(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = A(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M                                                 5610   
            S = 0.0D0
C
            DO 250 K = L, N
  250       S = S + A(J,K) * A(I,K)
C
            DO 260 K = L, N
               A(J,K) = A(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N                                                 5620   
  280    A(I,K) = SCALE * A(I,K)
C
  290    X = DMAX1(X,DABS(W(I))+DABS(RV1(I)))
  300 CONTINUE
C
C           ;;;;;; ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
C                FOR I=N STEP -1 UNTIL 1 DO --
C
      DO 400 II = 1, N
         I = N + 1 - II                                                  5630   
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0D0) GO TO 360
C
         DO 320 J = L, N
C
C
  320    A(J,I) = (A(I,J) / A(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0D0                                                    5640   
C
            DO 340 K = L, N
  340       S = S + A(I,K) * A(K,J)
C
            DO 350 K = L, N
               A(K,J) = A(K,J) + S * A(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            A(I,J) = 0.0D0                                               5650   
            A(J,I) = 0.0D0
  380    CONTINUE
C
  390    A(I,I) = 1.0D0
         G = RV1(I)
         L = I
  400 CONTINUE
C
      IF (M .GE. N .OR. IP .EQ. 0) GO TO 510
      M1 = M + 1                                                         5660   
C
      DO 500 I = M1, N
C
         DO 500 J = 1, IP
            B(I,J) = 0.0D0
  500 CONTINUE
C
C        ;;;;;; DIAGONALIZATION OF THE BIDIAGONAL FORM
C
  510 EPS = MACHEP * X                                                   5670   
C
C
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C
C       ;;;; TEST FOR SPLITTING.                                    0023
C                FOR L=K STEP -1 UNTIL 1 DO --
  520    DO 530 LL = 1, K                                                5680   
            L1 = K - LL
            L = L1 + 1
            IF (DABS(RV1(L)) .LE. EPS) GO TO 565
C
C       ;;;;; RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT             002
C                THROUGH THE BOTTOM OF THE LOOP
            IF (DABS(W(L1)) .LE. EPS) GO TO 540
  530    CONTINUE
C
C      ;;;;; CANCELLATION OF RV1(L) IF L GREATER THAN 1                  5690   
  540    C = 0.0D0
         S = 1.0D0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (DABS(F) .LE. EPS) GO TO 565
            G = W(I)
            H = DSQRT(F*F+G*G)
            W(I) = H                                                     5700   
            C = G / H
            S = -F / H
            IF (IP .EQ. 0) GO TO 560
C
            DO 550 J = 1, IP
               Y = B(L1,J)
               Z = B(I,J)
               B(L1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C
  550       CONTINUE                                                     5710   
C
  560    CONTINUE
C
C      ;;;;; TEST FOR CONVERGENCE
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C
C    ;;;;; SHIFT FROM BOTTOM 2 BY 2 MINOR
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1                                                   5720   
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = ((Y - Z) * (Y + Z) + (G - H) * (G + H)) / (2.0D0 * H * Y)
         G = DSQRT(F*F+1.0D0)
         F = ((X - Z) * (X + Z) + H * (Y / (F + DSIGN(G,F)) - H)) / X
C
C          ;;;;; NEXT QR TRANSFORMATION
         C = 1.0D0                                                       5730   
         S = 1.0D0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = DSQRT(F*F+H*H)
            RV1(I1) = Z                                                  5740   
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
C
            DO 570 J = 1, N
               X = A(J,I1)
               Z = A(J,I)                                                5750   
               A(J,I1) = X * C + Z * S
               A(J,I) = -X * S + Z * C
  570       CONTINUE
C
            Z = DSQRT(F*F+H*H)
            W(I1) = Z
C
C         ;;;;; ROTATION CAN BE ARBITRARY IF Z IS ZERO
            IF (Z .EQ. 0.0D0) GO TO 580
            C = F / Z                                                    5760   
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (IP .EQ. 0) GO TO 600
C
            DO 590 J = 1, IP
               Y = B(I1,J)
               Z = B(I,J)
               B(I1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C                                   5770   
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0D0
         RV1(K) = F
         W(K) = X
         GO TO 520
C
C      ;;;;; CONVERGENCE ;;;;;;;;;:                                 0023 5780   
  650    IF (Z .GE. 0.0D0) GO TO 700
C
C       ;;;;; W(K) IS MADE NON-NEGATIVE
         W(K) = -Z
C
         DO 690 J = 1, N
  690    A(J,K) = -A(J,K)
C
  700 CONTINUE
C                                                                        5790   
      GO TO 1001
C
C       ;;;;; SET ERROR -- NO CONVERGENCE TO A                       002
C                SINGULAR VALUE AFTER 30 ITERATIONS
 1000 IERR = K
 1001 RETURN
C
C         ;;;;; LAST CARD OF MINFIT
      END
C*ID* MINIM    C109.PTOLEMY.FITTERS                         PER705  21:  5800   
      SUBROUTINE MINIM(DIAG,EPS,X,XP,G,GP,S,SIG,Y,ETA,RHO,H,FS,
     1FJACOB,NA,N,IPRINT,MAXFUN,ISIG)
C     VARIABLE-METRIC MINIMIZATION ; DX = CONST * (-H * G)
C     WHERE G IS LOCAL DERIVATIVE ,H IS INVERSE SECOND-DERIVATIVE MATRIX
C     H UPDATED BY DAVIDON OR FLETCHER PRESCRIPTION
C
C     MHM DECEMBER 29 1975
C
C     DIAG : IF GE.0 MINIMUM NOT ALWAYS LOCATED ALONG -H*G
C            IF LT.0 MINIMUM LOCATED ALONG -H*G AT EACH STEP             5810   
C            /DIAG/ CONTROLS METRIC MODIFICATION AND STEP VATIATION
C            IF /DIAG/.GT.1 FLETCHER WITH FACTOR=1/(/DIAG/)
C            IF /DIAG/.LT.1 DAVIDON WITH FACTOR=/DIAG/
C            IF /DIAG/.EQ.0  OR 1 FACTOR=1/3
C     EPS  : RELATIVE ACCURACY REQUESTED
C     X    : VECTOR OF VARIABLE-VALUES  ( DIMENSION N)
C            INITIALLY CONTAINS GUESS AT POSITION OF MINIMUM
C            RETURNS VALUES AT LAST ITERATION ; THESE WILL BE EITHER
C            VALUES AT MIN TO ACCURACY EPS, OR A BEST APPROACH TO THESE
C            IN MAXFUN CALLS OF SUBROUTINE CALGRA                        5820   
C     XP,G,GP,S,SIG,Y,ETA,RHO : 8 WORK ARRAYS OF LENGTH N
C     H    : N BY N MATRIX - RETURNS VARIANCE-COVARIANCE MATRIX
C     FS   : VECTOR (DIMENSION NA) OF FUNCTIONS IN SUM-OF-SQUARES
C     FJACOB: WORK ARRAY OF DIMENSION NA BY N
C     N    : NO OF PARAMETERS
C     NA   : NO OF TERMS IN CHISQ SUM
C     IPRINT: PRINT CONTROL (USES MINMON)
C     MAXFUN : UPPER LIMIT ON NUMBER OF FUNCTION CALLS
C     ISIG : AT RETURN ISIG = 0 MEANS NORMAL COMPLETION
C                      ISIG = 1 - STOPPED AFTER MAXFUN FUNCTION CALLS    5830   
C                      ISIG =-5 - ABNORMAL CONDITIONS ENCOUNTERED
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,SIG,Y,ETA,RHO,H,FS,FJACOB
C
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),SIG(N),Y(N),ETA(N),RHO(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
C
C     INITIAL FUNCTION-CALL
C                                                                        5840   
      IF(DIAG.GE.0) ND=0
      IF(DIAG.LT.0) ND=1
      DG=DABS(DIAG)
      IF(DG.LT.1) W=DG
      IF((DG.EQ.0).OR.(DG.EQ.1)) W=1./3.
      IF(DG.GT.1) W=1/DG
      ITYP=0
      IF(DG.LT.1) ITYP=1
      write (6, 3985)
 3985 FORMAT(1H0,2X,'++ VARIABLE-METRIC MINIMIZATION ++')                5850   
      IF((ND.EQ.0).AND.(ITYP.EQ.0)) write (6, 4050) W
 4050 FORMAT(1H ,' FLETCHER METRIC,INCOMPLETE LINEAR SEARCH, W=',G12.6)
      IF((ND.EQ.0).AND.(ITYP.EQ.1)) write (6, 4450) W
 4450 FORMAT(1H ,' DAVIDON METRIC,INCOMPLETE LINEAR SEARCH, W=',G12.6)
      IF((ND.EQ.1).AND.(ITYP.EQ.1)) write (6, 4650) W
 4650 FORMAT(1H ,' DAVIDON METRIC,COMPLETE LINEAR SEARCH, W=',G12.6)
      IF((ND.EQ.1).AND.(ITYP.EQ.0)) write (6, 4850) W
 4850 FORMAT(1H ,' FLETCHER METRIC,COMPLETE LINEAR SEARCH, W=',G12.6)
      ISCREW=1
      ISIG=10                                                            5860   
      NIT=1
      AA=1
      UPRLIM=1.D+18
      NCALL=1
      DO 5 I=1,N
    5 XP(I)=X(I)
      chisq = CALGRA(NA,N,FS,XP,FJACOB,GP)
      FXP=0
      DO 10 I=1,NA
   10 FXP=FXP+FS(I)*FS(I)                                                5870   
      FX1=FXP
      FX2=0
      FX3=0
C
C     INITIAL ESTIMATE OF METRIC MATRIX
C
   30 DO 50 I=1,N
      DO 50 J=I,N
      H(I,J)=0
      DO 40 IA=1,NA                                                      5880   
   40 H(I,J)=H(I,J)+FJACOB(IA,I)*FJACOB(IA,J)
      H(I,J)=2*H(I,J)
   50 H(J,I)=H(I,J)
      CALL SYMINV(H,N,Y,0,DET,N)
      DO 60 I=2,N
      JTOP=I-1
      DO 60 J=1,JTOP
   60 H(I,J)=H(J,I)
      GO TO 90
C                                                                        5890   
C     PRINT-OUT OF INFORMATION AT START OF ITERATION
C
   70 ISIG=0
      IF(IPRINT.EQ.1) CALL MINMON(NIT,NCALL,NCALL,N,NA,X,FX,FJACOB)
C
C     MAIN ITERATIVE SECTION  : INCOMPLETE LINEAR SEARCH
C
  132 ALPHA=AA
      IF((AA.LT.1).AND.(ND.EQ.0)) AA=W
      NTURN=1                                                            5900   
  120 SIGX=ALPHA*SXNEW
      IF(SIGX.GE.EPS/10) GO TO 1120
      ISIG=-5
      write (6, 7000) SIGX,NTURN,NIT
 7000 FORMAT(1H0,2X,'-- MINIM TRIED TO TAKE STEP OF LENGTH ',G12.6/
     1 3X,'AT THE ',I2,'TH TURN OF ITERATION ',I3,'--')
      GO TO 250
 1120 DO 130 I=1,N
      SIG(I)=ALPHA*S(I)
  130 XP(I)=X(I)+SIG(I)                                                  5910   
      IF ( NCALL .GE. MAXFUN )  GO TO 215
      chisq = CALGRA(NA,N,FS,XP,FJACOB,GP)
      NCALL=NCALL+1
  360 G1=0
      DO 282 I=1,N
  282 G1=G1+S(I)*GP(I)
  292 FXP=0
      DO 140 I=1,NA
  140 FXP=FXP+FS(I)*FS(I)
      IF(FXP.GT.UPRLIM) GO TO 157                                        5920   
C
C     FIND SATISFACTORY LENGTH FOR STEP ALONG VECTOR S
C
  138 IF(ND.NE.0) GO TO 4138
      IF(FXP.LT.FX) GO TO 159
      IF((G1.GT.0).AND.(NTURN.EQ.1)) GO TO 147
  157 ALPHA=ALPHA*W
      AA=AA*W
      GO TO 148
  147 ALPHA=FCUBIC(FX,G0,FXP,G1,ALPHA)                                   5930   
  148 NTURN=NTURN+1
      GO TO 120
  159 IF(G1.GE.0) GO TO 150
      AA=AA/W
      AA=DMIN1(1.D0,AA)
      GO TO 150
C
C     MAIN ITERATIVE SECTION : LINEAR SEARCH TO MINIMUM
C
 4138 IF(FXP.LT.FX) GO TO 4159                                           5940   
      IF(G1.LT.0) GO TO 4869
      IKEY=1
      GO TO 4140
 4150 IF(FZ.GE.FX) GO TO 4869
      GO TO 4160
 4869 ALPHA=ALPHA*W
      AA=AA*W
      NTURN=NTURN+1
      GO TO 120
 4159 IF(G1.GE.0) GO TO 4871                                             5950   
      AA=AA/W
      AA=DMIN1(1.D0,AA)
      GO TO 150
 4871 IKEY=2
      GO TO 4140
 4170 IF(FZ.GE.FXP) GO TO 4877
 4160 FXP=FZ
      DO 4867 I=1,N
      XP(I)=ETA(I)
 4867 GP(I)=RHO(I)                                                       5960   
      GO TO 150
 4877 AA=AA*W
      GO TO 150
 4140 ALPHA=FCUBIC(FX,G0,FXP,G1,ALPHA)
      DO 4873 I=1,N
 4873 ETA(I)=X(I)+ALPHA*X(I)
      chisq = CALGRA(NA,N,FS,ETA,FJACOB,RHO)
      NCALL=NCALL+1
      FZ=0
      DO 4875 I=1,NA                                                     5970   
 4875 FZ=FZ+FS(I)*FS(I)
      IF(FZ.LT.UPRLIM) GO TO 4220
      AA=AA*W
      GO TO 4869
 4220 GO TO (4150,4170),IKEY
C
C     ITERATION COMPLETE
C     GO ON TO CALCULATE S,UPDATE METRIC AND TEST FOR CONVERGENCE
C
  150 SX=SXNEW                                                           5980   
      T=TNEW
      DO 160 I=1,N
  160 Y(I)=GP(I)-G(I)
      SSX=SX*ALPHA
      FX3=FX2
      FX2=FX1
      FX1=FXP
C
C     UPDATE METRIC MATRIX
C                                                                        5990   
      SY=0
      EY=0
      DO 180 I=1,N
      ETA(I)=0
      DO 170 J=1,N
  170 ETA(I)=ETA(I)+H(I,J)*Y(J)
      SY=SY+SIG(I)*Y(I)
  180 EY=EY+ETA(I)*Y(I)
      IF(ITYP.EQ.0) GO TO 2280
      DO 280 I=1,N                                                       6000   
      DO 280 J=1,N
  280 H(I,J)=H(I,J)+SIG(I)*SIG(J)/SY-ETA(I)*ETA(J)/EY
      GO TO 90
 2280 Z=1+EY/SY
      DO 430 I=1,N
      DO 430 J=1,N
  430 H(I,J)=H(I,J)-(SIG(I)*ETA(J)+ETA(I)*SIG(J)-Z*SIG(I)*SIG(J))/SY
C
C     CALCULATES DERIVATIVE IN DIRECTION S=-H*G
C     CHECKS THAT THIS IS NEGATIVE (H POSITIVE DEFINITE)                 6010   
C     IF NOT RESETS H TO UNIT MATRIX OR QUITS(SECOND TIME)
C
   90 SXNEW=0
      TNEW=0
      G0=0
      DO 136 I=1,N
      S(I)=0
      DO 100 J=1,N
  100 S(I)=S(I)-H(I,J)*GP(J)
      TES=DABS(S(I))                                                     6020   
      IF(TES.LE.SXNEW) GO TO 136
      SXNEW=TES
      TNEW=DABS(XP(I))
  136 G0=G0+S(I)*GP(I)
      IF(G0.LE.0) GO TO 300
  824 IF(ISCREW.GT.2) GO TO 850
      DO 855 I=1,N
      DO 855 J=1,N
      H(I,J)=0
      IF(I.NE.J) GO TO 855                                               6030   
      H(I,I)=1
  855 CONTINUE
      ISCREW=ISCREW+1
      GO TO 90
  850 ISIG=-5
      write (6, 2900) NIT,G0
 2900 FORMAT(1H0,2X,'--- IT ',I4,' G0 = ',G12.6/
     13X,'    METRIC CANNOT BE POSITIVE DEFINITE ___')
      GO TO 250
C                                                                        6040   
C      CONVERGENCE TEST
C      APPLIED ONLY AT 3 RD AND LATER ITERATIONS
C
  300 DO 210 I=1,N
      G(I)=GP(I)
  210 X(I)=XP(I)
      FX=FXP
      IF(G0.EQ.0) GO TO 250
      IF(ISIG.EQ.10) GO TO 70
      IF(NIT.LE.3) GO TO 220                                             6050   
      T1=DABS(FX3-FX1)
      T2=EPS*DABS(FX1)
      IF(((SX.LE.EPS*T).AND.(SXNEW.LE.EPS*TNEW)).OR.(T1.LE.T2))GO TO 250
      IF(NCALL.LT.MAXFUN) GO TO 220
 215  ISIG=1
      GO TO 250
  220 NIT=NIT+1
      GO TO 70
C
  250 DO 444 I=1,N                                                       6060   
      DO 444 J=1,N
  444 H(I,J)=2*H(I,J)
      RETURN
      END
C*ID* MINMON   C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE MINMON ( ITCNT, NFCALL, NJCALL, NPARAM, NPTS, PARAM, F,
     1  FS )
C
C     PRINT ITERATION INFORMATION FOR DLSMIN AND OTHER MINIMIZERS
C                                                                        6070   
C     10/23/75 - FIRST VERSION - S. P.
C
      IMPLICIT REAL*8  ( A-H, O-Z )
C
      DIMENSION PARAM(NPARAM), FS(NPTS)
C
      FF = F/NPTS
      write (6, 53) ITCNT, NFCALL, NJCALL, FF, PARAM
 53   FORMAT ( '0ITERATION', I5, I15, ' FUNCTION CALLS.',
     1  I10, ' JACOBIAN CALLS',                                          6080   
     2  5X, '(CHI. SQ.)/PT =', G16.8, 5X, 'PARAMETERS:' /
     3  ( 8G16.6 ) )
      RETURN
      END
C*ID* OCOPTR   C109.PTOLEMY.FITTERS                         PER705  22:
      SUBROUTINE OCOPTR(IDIM,M,N,X,F,X0,F0,FCN,NFCALL,FJ,MAXFCN,ICONV,
     1                  EPS,H,FL0,FMIN,IDEB)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),X0(N),FJ(IDIM,M)                                    6090   
      EXTERNAL FCN
      DIMENSION FK(50),WORK(50),FN(50),FM(50),P(50),Q(50)
      character*1 ARRAY(8), TRACE(20)
      DATA ARRAY/'T','D','E','M','N','B','1',' '/
      DATA TRACE/20*' '/
C
C     INITIALIZE
C
      NFCALL = 1
      IERR = 0                                                           6100   
      IT = 0
      ICONV = 0
C
C     INITIALIZE J TO UNIT MATRIX
C
      DO 9  I = 1, M
         DO 8  J = 1, N
            FJ(J,I) = 0
  8      CONTINUE
         FJ(I,I) = 1                                                     6110   
  9   CONTINUE
C
C     START
C
      CALL FCN(N,X,F)
      CALL MOVE(X0,X,N)
      F0 = F
      FA = F
      FB = 2*(F - FMIN)
      DO 10 I = 1,M                                                      6120   
C
         DO 5 J = 1,N
            X(J) = X0(J) + FJ(J,I)*H
    5    CONTINUE
         CALL FCN(N,X,F)
         NFCALL = NFCALL + 1
         IF (NFCALL .GT. MAXFCN ) RETURN
         FK(I) = (F - F0)/H - H/2
C
   10 CONTINUE                                                           6130   
C
C     STEP 0
C
    1 CONTINUE
      FNS = 0
C
C     STEP 1
C
  100 CONTINUE
C                                                                        6140   
      IF ( IDEB .EQ. 1) write (6, 3000) (TRACE(I),I=1,10)
 3000 FORMAT(/' TRACE IS',5X,10A1)
      IF(IDEB.EQ.1)write (6, 3010) IT,NFCALL,F0
 3010 FORMAT(2I6,D24.14)
      IT = IT + 1
      DO 3020 I = 1,20
         TRACE(I) = ARRAY(8)
 3020 CONTINUE
      IPT = 1
  101 DO 102 I = 1,N                                                     6150   
         X(I) = 0
         DO 103 J = 1,M
            X(I) = X(I) + FJ(I,J)*FK(J)
  103    CONTINUE
  102 CONTINUE
      A = VTMV(FK,FK,M)
      IF ( A .LT. EPS ) GOTO 1000
C
      T = DSQRT(A)
      IF ( FK(1) .LT. 0 ) T = -T                                         6160   
      FK(1) = FK(1) + T
      DO 105 I = 1,N
         WORK(I) = X(I)/T + FJ(I,1)
  105 CONTINUE
      DO 120 I = 1,N
         DO 110 J = 1,M
            FJ(I,J) = FJ(I,J) - WORK(I)*FK(J)/FK(1)
  110    CONTINUE
  120 CONTINUE
      IF ( FNS .EQ. 0 ) GOTO 200                                         6170   
      TEMP1 = VTMV(FK,FN,M)
      TEMP2 = FK(1)*T
      DO 130 I = 1,M
         FN(I) = FN(I) - FK(I)*TEMP1/TEMP2
  130 CONTINUE
C
C     STEP 2
C
  200 CONTINUE
      S = 0                                                              6180   
      FL = FL0
      FK(1) = -T
      TEMP1 = DMIN1(2*(FB-F0)/A,2*(F0-FMIN)/A)
      IF ( TEMP1 .GT. 1.0D0 ) GOTO 210
      TRACE(IPT) = ARRAY(1)
      IPT = IPT + 1
      T = T * TEMP1
  210 FB = FA
      FA = F0
C                                                                        6190   
C     STEP 3
C
  300 CONTINUE
      IF ( -FK(1)*T .LT. EPS ) GOTO 1100
      DO 310 I = 1,N
         X(I) = X0(I) + FJ(I,1)*T
  310 CONTINUE
      CALL FCN(N,X,F)
      NFCALL = NFCALL + 1
      IF ( F .LT. F0 ) GOTO 350                                          6200   
      T = T/2
      FL = 0.5D0
      IF ( F .LT. FA ) GOTO 400
      TRACE(IPT) = ARRAY(2)
      IPT = IPT + 1
      GOTO 300
C
C     STEP 3B
C
  350 A = (F-F0)/T                                                       6210   
      B = (A-FK(1))/T
      C = T
      S = S + T
      T = T*FL
      CALL MOVE(X0,X,N)
      F0 = F
      FK(1) = 2*A - FK(1)
      IF ( FL .LT. 1.0D0 ) GOTO 400
      IF ( FK(1)/T + 2*B .GE. 0 ) GOTO 400
      TRACE(IPT) = ARRAY(3)                                              6220   
      IPT = IPT + 1
      GOTO 300
C
C     STEP 4
C
  400 DO 410 I = 1,N
         X(I) = X0(I) + FJ(I,1)*H
  410 CONTINUE
      CALL FCN(N,X,F)
      NFCALL = NFCALL + 1                                                6230   
      IF ( NFCALL .GT. MAXFCN ) RETURN
      FK(1) = (F-F0)/H
      B = (FK(1)-A)/(H+C)
      FK(1) = FK(1) - B*H
      IF (FK(1)/T + 2*B .GE. 0 ) GOTO 500
      IF ( FL .LT. 1.0D0 ) GOTO 415
      TRACE(IPT) = ARRAY(4)
      IPT = IPT + 1
      GOTO 300
  415 B = -FK(1)/C                                                       6240   
C
C     STEP 5
C
  500 FM(1) = 1-2*B
      FMS = FM(1)*FM(1)
      IF ( M .LT. 2 ) GOTO 511
      DO 510 I = 2,M
         DO 520 J = 1,N
            X(J) = X0(J) + FJ(J,I)*H
  520    CONTINUE                                                        6250   
         CALL FCN(N,X,F)
         NFCALL = NFCALL + 1
      IF ( NFCALL .GT. MAXFCN ) RETURN
         FK(I) = (F-F0)/H - H/2
         FM(I) = -FK(I)/S
         FMS = FMS + FM(I)*FM(I)
  510 CONTINUE
  511 CONTINUE
      IF ( FMS .LT. EPS ) GOTO 1
      FNU = FM(1)                                                        6260   
      FMU = FNU - FMS
      IF ( FNS .EQ. 0 ) GOTO 700
      A = VTMV(FM,FN,M)
      B = A/FMS
      FNS = FNS - A*B
      IF ( FNS .GE. EPS ) GOTO 530
      A = 0
      TRACE(IPT) = ARRAY(5)
      IPT = IPT + 1
      GOTO 600                                                           6270   
  530 DO 540 J = 1,M
         FN(J) = FN(J) - FM(J)*B
  540 CONTINUE
      A = FN(1)/FNS
C
C     STEP 6
C
  600 CONTINUE
      DO 610 J = 1,M
         FN(J) = FN(J)*A                                                 6280   
  610 CONTINUE
      FNS = FN(1)
      B = FMU*FNU/FMS + FNS
      IF ( EPS .LE. B ) GOTO 800
C
C     STEP 7
C
  700 CONTINUE
      DO 710 J = 1,M
         FN(J) = - FM(J)*FNU/FMS                                         6290   
  710 CONTINUE
      FN(1) = FN(1) + 1
      FNS = FN(1)
      B = 1 - FNU
      TRACE(IPT) = ARRAY(6)
      IPT = IPT + 1
C
C     STEP 8
C
  800 CONTINUE                                                           6300   
      IF ( FMS*FNS .GT. FMU*FNU ) GOTO 810
      TRACE(IPT) = ARRAY(7)
      IPT = IPT + 1
      DEL = DSQRT(FNU/FMU)
      GA = 0
      GOTO 900
  810 A = B - FMU
      C = B + FNU
      DEL = DSQRT(C/A)
      GA = DSQRT((1-FMU*FNU/(FMS*FNS))/(A*B))                            6310   
      IF ( C .LT. A ) GA = -GA
C
C     STEP 9
C
  900 CONTINUE
      FL = FNU + FMU*DEL + FMS*FNS*GA
      TEMP1 = DEL - FNS*GA
      TEMP2 = GA*FNU
      TEMP3 = (1+FNS*GA)/FL
      TEMP4 = GA*FMU/FL                                                  6320   
      DO 910 J = 1,M
         P(J) = FM(J)*TEMP1 + FN(J)*TEMP2
         Q(J) = FM(J)*TEMP3 - FN(J)*TEMP4
  910 CONTINUE
      TEMP1 = VTMV(Q,FK,M)
      DO 920 J = 1,M
         FK(J) = FK(J) + P(J)*TEMP1
  920 CONTINUE
C
      DO 930 I = 1,N                                                     6330   
         TEMP1 = 0
         DO 935 K = 1,M
            TEMP1 = TEMP1 + FJ(I,K)*Q(K)
  935    CONTINUE
         DO 940 J = 1,M
            FJ(I,J) = FJ(I,J) + TEMP1*P(J)
  940    CONTINUE
  930 CONTINUE
      IF ( FNS .EQ. 0 ) GOTO 100
      TEMP1 = FNS*(1+GA*FMU*FNU/FL)                                      6340   
      TEMP2 = (1+DEL)*FMU*FNU/FL
      DO 950 J = 1,M
         FN(J) = FM(J)*TEMP1 - FN(J)*TEMP2
  950 CONTINUE
      FNS = VTMV(FN,FN,M)
      GOTO 100
C
C     STEP 10
C
 1000 DO 1010 J = 1,N                                                    6350   
         X(J) = X0(J) - X(J)
 1010 CONTINUE
      CALL FCN(N,X,F)
      NFCALL = NFCALL+1
      IF ( F0 .LT. F ) GOTO 1100
      CALL MOVE(X0,X,N)
      F0 = F
C
 1100 ICONV = 1
      RETURN                                                             6360   
      END
      DOUBLE PRECISION FUNCTION VTMV(ARRAY1,ARRAY2,IDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY1(IDIM),ARRAY2(IDIM)
      TEMP = 0
      DO 10 I = 1,IDIM
         TEMP = TEMP + ARRAY1(I)*ARRAY2(I)
   10 CONTINUE
      VTMV = TEMP
      RETURN                                                             6370   
      END
C
      SUBROUTINE MTMM(ARRAY1,ARRAY2,RESULT,IDIM1,IDIM2,IDIM3)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY1(IDIM1,IDIM2),ARRAY2(IDIM2,IDIM3),
     1 RESULT(IDIM1,IDIM3)
      DO 10 I = 1,IDIM1
         DO 20 J = 1,IDIM3
            TEMP = 0
            DO 30 K = 1,IDIM2                                            6380   
               TEMP = TEMP + ARRAY1(I,K)*ARRAY2(K,J)
   30       CONTINUE
            RESULT(I,J) = TEMP
   20    CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE MOVE(X0,X,N)
      IMPLICIT REAL*8(A-H,O-Z)                                           6390   
      DIMENSION X0(N),X(N)
      DO 10 I = 1,N
         X0(I) = X(I)
   10 CONTINUE
      RETURN
      END
C*ID* QMTRIX   C109.PTOLEMY.FITTERS                         PER705  21:
C
C     -------------------------------------------------------------
C                                                                        6400   
      SUBROUTINE QMTRIX(N,RJTJ,D,RJTF,S,DELTA,QPARM)
C
CDC      ALLOCLEVEL RJTJ,D,RJTF,S,DELTA
      INTEGER I,J,N
      REAL*8  RJTJ(N,N),D(N),RJTF(N),S(N),DELTA(N)
      REAL*8  QPARM
      REAL*8  DETERM
C
C       THIS ROUTINE SOLVES THE MARQUARDT MATRIX EQUATION
C     ((J TRANSPOSE) * J + QPARM * S) * DELTA = -(J  TRANSPOSE) * F      6410   
C     WHERE  J  IS THE JACOBIAN MATRIX, S  IS THE SCALE CONTROL VECTOR,
C     DELTA  IS THE CORRECTION VECTOR AND  F  IS THE FUNCTION VECTOR.
C     THE MARQUARDT PARAMETER  QPARM  IS INCREASED IN VALUE IF THE
C     ORIGINAL MARQUARDT MATRIX IS SINGULAR.
C
C       THIS ROUTINE REQUIRES THE SUBROUTINE  SYMINV  OR ITS EQUIVALENT.
C     SYMINV  SOLVES THE MATRIX EQUATION  A * X = B  BY GAUSSIAN
C     ELIMINATION. THE MATRIX  A  IS SYMMETRIC. THE DETERMINANT AND THE
C     INVERSE OF  A  ARE ALSO OBTAINED. ONLY THE UPPER TRIANGULAR
C     PORTION OF  A  IS USED AND ONLY THE UPPER TRIANGULAR PORTION OF    6420   
C     (A INVERSE)  IS STORED.
C
C     ON INPUT
C
C        N  IS THE DIMENSION OF THE VECTORS  D, DELTA, RJTF  AND  S
C           AND THE SQUARE MATRIX IN  RJTJ.
C
C        RJTJ  CONTAINS THE LOWER TRIANGLE OF  (J TRANSPOSE) * J
C           STORED BY ROWS
C                                                                        6430   
C        D  CONTAINS THE DIAGONAL ELEMENTS OF  RJTJ.
C
C        RJTF  CONTAINS THE ELEMENTS OF  (J TRANSPOSE) * F.
C
C        S  CONTAINS THE MARQUARDT MULTIPLIERS.
C
C        QPARM  IS THE MARQUARDT PARAMETER.
C
C     ON OUTPUT
C                                                                        6440   
C        RJTJ  CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF THE
C           MARQUARDT MATRIX STORED BY ROWS. THE ELEMENTS IN  RJTJ
C           BELOW THE MAIN DIAGONAL HAVE NOT BEEN ALTERED BY THE
C           ROUTINE.
C
C        DELTA  CONTAINS THE NEGATIVE CORRECTION VECTOR.
C
C        QPARM  IS THE UPDATED MARQUARDT PARAMETER.
C
C     -------------------------------------------------------------      6450   
C
C     ASSEMBLE THE MARQUARDT MATRIX. IN THIS SECTION THE DIAGONAL
C       ELEMENTS OF  (J  TRANSPOSE) * J  ARE AUGMENTED BY  QPARM * S
C       SO THAT WHEN FINISHED
C
C       RJTJ  -  CONTAINS THE LOWER TRIANGLE OF THE MARQUARDT MATRIX
C
C       ALSO
C
C       DELTA(I)  -  CONTAINS THE VECTOR  RJTF                           6460   
C
  100 DO 110 I = 1, N
         DELTA(I) = RJTF(I)
         RJTJ(I,I) = D(I) + QPARM * S(I)
  110 CONTINUE
C
C     SOLVE THE MARQUARDT MATRIX EQUATION. WHEN FINISHED
C
C       RJTJ     - CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF THE
C                  MARQUARDT MATRIX. THE ELEMENTS BELOW THE MAIN         6470   
C                  DIAGONAL HAVE NOT BEEN ALTERED
C       DELTA(I) - CONTAINS THE NEGATIVE CORRECTION VECTOR
C
      DO 120 I = 1, N
      DO 120 J = I, N
         RJTJ(I,J) = RJTJ(J,I)
  120 CONTINUE
C
      CALL SYMINV(RJTJ,N,DELTA,1,DETERM,N)
C                                                                        6480   
C     IF THE DETERMINANT OF THE MARQUARDT MATRIX IS ZERO, QPARM  IS
C       INCREASED IN VALUE AND THE MARQUARDT MATRIX IS REEVALUATED
C
      IF (DETERM .NE. 0.0D0) GO TO 300
      QPARM = 2.0D0 * QPARM
      IF (QPARM .EQ. 0.0D0) QPARM = 1.0D0
      GO TO 100
C
  300 RETURN
C     ----- LAST CARD OF SUBROUTINE  QMTRIX  ----                        6490   
      END
C*ID* QUAVER   C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE QUAVER(DIAG,EPS,TOL,X,XP,G,GP,S,SIG,Y,ETA,RHO,H,FS,
     1FJACOB,NA,N,IPRINT,MAXFUN,ISIG)
C     QUASI-NEWTON MINIMIZATION ; DX = CONST * (-H * G)
C     WHERE G IS LOCAL DERIVATIVE ,H IS INVERSE SECOND-DERIVATIVE MATRIX
C     H IS APPROXIMATED BY JT*J WITH J JACOBIAN OF RESIDUALS
C
C     MHM DECEMBER 29 1975
C                                                                        6500   
C     MODIFIED TO WORK IN TERMS OF PSEUDO-INVERSES
C     MHM JANUARY 29 1976
C
C     DIAG : STEPSIZE CONTROL
C            W = /DIAG/ OR /DIAG/**-1 WHICHEVER IS LESS THAN 1
C            IF DIAG=0 OR /DIAG/=1 W=1/3
C            WHEN CHISQ INCREASES REPLACE STEP BY W*STEP
C     EPS  : RELATIVE ACCURACY REQUESTED
C     TOL  : SINGULAR VALUES OF J SET 0 IF LT. TOL*(MAXIMUM VAL)
C     X    : VECTOR OF VARIABLE-VALUES  ( DIMENSION N)                   6510   
C            INITIALLY CONTAINS GUESS AT POSITION OF MINIMUM
C            RETURNS VALUES AT LAST ITERATION ; THESE WILL BE EITHER
C            VALUES AT MIN TO ACCURACY EPS, OR A BEST APPROACH TO THESE
C            IN MAXFUN CALLS OF SUBROUTINE CALGRA
C     XP,G,GP,S,SIG,Y,RHO : 7 WORK ARRAYS OF LENGTH N
C     ETA  : VECTOR (DIMENSION N) OF SINGULAR VALUES OF JACOBIAN
C     H    : N BY N MATRIX - RETURNS VARIANCE-COVARIANCE MATRIX
C     FS   : VECTOR (DIMENSION NA) OF WEIGHTED RESIDUALS
C     FJACOB: NA BY N JACOBIAN MATRIX DF(I)/DX(J) OF RESIDUALS
C     N    : NO OF PARAMETERS                                            6520   
C     NA   : NO OF TERMS IN CHISQ SUM
C     IPRINT: PRINT CONTROL (USES MINMON)
C     MAXFUN : UPPER LIMIT ON NUMBER OF FUNCTION CALLS
C     ISIG : AT RETURN ISIG = 0 MEANS NORMAL COMPLETION
C                      ISIG = 1 - STOPPED AFTER MAXFUN FUNCTION CALLS
C                      ISIG =-5 - ABNORMAL CONDITIONS ENCOUNTERED
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,SIG,Y,ETA,RHO,H,FS,FJACOB
C                                                                        6530   
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),SIG(N),Y(N),ETA(N),RHO(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
C
C     INITIAL FUNCTION-CALL
C
      DG=DABS(DIAG)
      IF(DG.LT.1) W=DG
      IF((DG.EQ.0).OR.(DG.EQ.1)) W=1./3.
      IF(DG.GT.1) W=1/DG
      write (6, 4050)                                                    6540   
 4050 FORMAT(1H0,2X,'++ QUASI-NEWTON MINIMIZATION (PSEUDO-INVERSE) ++')
      write (6, 4450) W,TOL
 4450 FORMAT(1H ,5X,'STEP-REDUCTION FACTOR = ',G12.6/
     16X,'MAXIMUM RATIO OF LARGEST TO SMALLEST SING.VAL. = ',G12.6)
      ISIG=10
      NIT=0
      AA=1
      UPRLIM=1.D+18
      NCALL=1
      DO 5 I=1,N                                                         6550   
    5 XP(I)=X(I)
      FXP=CALGRA(NA,N,FS,XP,FJACOB,GP)
      FX1=0
      FX2=0
      FX3=0
C
C     DETERMINE DIRECTION FOR NEXT STEP
C     BY SOLUTION OF LINEAR LSQ PROBLEM WITH JACOBIAN AND RESIDUALS
C
   30 DO 50 I=1,N                                                        6560   
      DO 50 IA=1,NA
   50 FJACOB(IA,I)=FJACOB(IA,I)*XP(I)
      CALL MINFIT(NA,NA,N,FJACOB,ETA,1,FS,IERR,RHO)
      IF(IERR.EQ.0) GO TO 32
      write (6, 4500) NIT,IERR
 4500 FORMAT(1H0,2X,'-- AT IT NO.',I3,' TROUBLE WITH ',I2,
     1'TH SINGULAR VAL OF J --')
      ISIG=-5
      GO TO 250
   32 TOP=0                                                              6570   
      DO 55 J=1,N
   55 IF(ETA(J).GT.TOP) TOP=ETA(J)
      THRESH=TOL*TOP
      NCHOP=0
      DO 34 I=1,N
      S(I)=0
      Y(I)=ETA(I)
      IF(ETA(I).GE.THRESH) GO TO 34
      Y(I)=0
      NCHOP=NCHOP+1                                                      6580   
   34 CONTINUE
      DO 60 K=1,N
      IF(Y(K).EQ.0) GO TO 60
      DO 36 I=1,N
   36 S(I)=S(I)-FJACOB(I,K)*FS(K)/Y(K)
   60 CONTINUE
      DO 80 J=1,N
   80 S(J)=XP(J)*S(J)
      IF(NCHOP.GT.0) write (6, 4570) NCHOP
 4570 FORMAT(1H ,2X,'THE ',I2,' SMALLEST SINGULAR VALUES ',              6590   
     1'SET TO ZERO')
      GO TO 90
C
C     PRINT-OUT OF INFORMATION AT START OF ITERATION
C
   70 ISIG=0
      IF(IPRINT.GE.1) CALL MINMON(NIT,NCALL,NCALL,N,NA,X,FX,FJACOB)
      IF(IPRINT.EQ.2) write (6, 4550) (ETA(I),I=1,N)
 4550 FORMAT(1H0,2X,'SINGULAR VALUES OF J '/ ( 3X, 10G13.6 ) )
C                                                                        6600   
C     MAIN ITERATIVE SECTION  : INCOMPLETE LINEAR SEARCH
C
  132 IF(AA.LT.1) AA=W
      ALPHA=AA
      NTURN=1
  120 DO 130 I=1,N
      SIG(I)=ALPHA*S(I)
  130 XP(I)=X(I)+SIG(I)
      FXP=CALGRA(NA,N,FS,XP,FJACOB,GP)
      IF(FXP.GT.UPRLIM) GO TO 157                                        6610   
      NCALL=NCALL+1
  360 G1=0
      DO 282 I=1,N
  282 G1=G1+S(I)*GP(I)
C
C     FIND SATISFACTORY LENGTH FOR STEP ALONG VECTOR S
C
  138 IF(FXP.LT.FX) GO TO 159
      IF((G1.GE.0).AND.(NTURN.EQ.1)) GO TO 147
  157 ALPHA=ALPHA*W                                                      6620   
      AA=AA*W
      GO TO 148
  147 ALPHA=FCUBIC(FX,G0,FXP,G1,ALPHA)
      AA=AA*W
  148 NTURN=NTURN+1
      IF(NTURN.GT.5) GO TO 149
      GO TO 120
  159 IF((G1.GE.0).OR.(NTURN.GT.1)) GO TO 30
      AA=AA/W
      AA=DMIN1(1.D0,AA)                                                  6630   
      GO TO 30
  149 ISIG=-5
      write (6, 3200)
 3200 FORMAT(1H0,2X,'--- F DOES NOT DECREASE ENOUGH IN 5 CHOPS ---')
      GO TO 250
C
C     ITERATION COMPLETE
C     GO ON TO TEST FOR CONVERGENCE
C
   90 FX3=FX2                                                            6640   
      FX2=FX1
      FX1=FXP
      SX=0
      T=0
      G0=0
      DO 136 I=1,N
      TES=DABS(S(I))
      IF(TES.LE.SX) GO TO 136
      SX=TES
      T=DABS(XP(I))                                                      6650   
  136 G0=G0+S(I)*GP(I)
      IF(G0.LE.0) GO TO 300
  850 ISIG=-5
      write (6, 2900) NIT,G0
 2900 FORMAT(1H0,2X,'--- IT ',I4,' G0 = ',G12.6/
     13X,'    METRIC CANNOT BE POSITIVE DEFINITE ___')
      GO TO 250
C
C     CONVERGENCE TEST
C     NOT APPLIED ON FIRST TWO PASSES                                    6660   
C
  300 DO 210 I=1,N
      G(I)=GP(I)
  210 X(I)=XP(I)
      FX=FXP
      IF(G0.EQ.0) GO TO 250
      IF(ISIG.EQ.10) GO TO 70
      IF(NIT.LT.1) GO TO 220
      IF((SX.LE.EPS*T).OR.(DABS(FX3-FX1).LE.EPS*DABS(FX1)))GO TO 250
      IF(NCALL.LT.MAXFUN) GO TO 220                                      6670   
      ISIG=1
      GO TO 250
  220 NIT=NIT+1
      GO TO 70
C
C     CALCULATE FINAL VARIANCE-COVARIANCE MATRIX AND QUIT
C
  250 DO 65 I=1,N
      DO 65 J=1,N
      H(I,J)=0                                                           6680   
      DO 40 K=1,N
      IF(ETA(K).LT.1.D-14) GO TO 40
      H(I,J)=H(I,J)+FJACOB(I,K)*FJACOB(J,K)/(ETA(K)*ETA(K))
   40 CONTINUE
      H(I,J)=H(I,J)*X(I)*X(J)
      IF(I.EQ.J) GO TO 65
      H(J,I)=H(I,J)
   65 CONTINUE
      NIT=NIT+1
      IF(IPRINT.GE.1) CALL MINMON(NIT,NCALL,NCALL,N,NA,X,FX,FJACOB)      6690   
      IF(IPRINT.EQ.2) write (6, 4550) (ETA(I),I=1,N)
      RETURN
      END
C*ID* SDIAG    C109.PTOLEMY.FITTERS                         PER705  10:
      SUBROUTINE SDIAG (NM, N, EPS, ILIMIT, A, D, X, E)
C
C     DIAGONALIZATION OF SYMMETRIC MATRICES
C
C     A - (NM X NM) - SYMMETRIC INPUT MATRIX IN EXPANDED FORM
C     D - (N) - WILL GET THE ARRAY OF N EIGENVALUES                      6700   
C     X - (NM X NM) - WILL GET THE N X N MATRIX OF EIGEVECTORS.
C         X MAY BE THE SAME AS A.
C     E - (N) - WORK ARRAY.
C
C     2/10/77 - SPECIAL CASE FOR 1 X 1 MATRICES - S. P.
C     12/28/79 - CFR%F - RPG
C
      implicit real*8 ( a-h, o-z )                                      implicit
      DIMENSION  A(NM,NM),D(NM),X(NM,NM),E(NM)
C                                                                        6710   
      IF ( N .LE. 1 )  GO TO 500
      DO 100 I = 1,N
      DO 100 J = 1,I
  100 X(I,J) = A(I,J)
C ----HOUSEHOLDER REDUCTION
      NM1 = N - 1
      DO 101 NN = 1,NM1
      I = N +1 - NN
      M = I - 2
      F = X(I,I -1)                                                      6720   
      G = F
      S = 0
      IF (M.LT.1) GOTO 99
      DO 102 K = 1,M
  102 S = S + X(I,K)*X(I,K)
   99 D(I) = S
      IF (S.EQ. 0) GOTO 103
      M = M + 1
      S = S + F * F
      G = DSQRT (S)                                                      6730   
      IF (F.GE.0.0) G=-G
      RH = S - F * G
      X(I,I-1) = F - G
      F = 0
      IF (M .LT. 1) GOTO 103
      DO 107 J = 1,M
      X(J,I) = X(I,J)/ RH
      S = 0
      DO 108 K = 1,J
  108 S = S + X(J,K) * X(I,K)                                            6740   
      JP1 = J + 1
      IF (M.LT.JP1) GOTO 98
      DO 109 K = JP1, M
  109 S = S + X(K,J) * X(I,K)
   98 E(J) = S/RH
  107 F = F + S * X(J,I)
   89 F = F / (RH + RH)
      DO 110 J = 1,M
  110 E(J)= E(J) - F * X(I,J)
      DO 111 J = 1,M                                                     6750   
      DO 111 K = 1,J
  111 X(J,K) = X(J,K) - X(I,J) * E(K) -X(I,K) * E(J)
C-----END  S NOT EQUAL TO ZERO
  103 E(I-1) = G
  101 CONTINUE
C-----END I LOOP
C-----ACCUMULATION  OF TRANSFORMATION  MATRICES
      D(1) = 0
      E(N)= 0
      B = 0                                                              6760   
      F = 0
      DO 200 I = 1,N
      M = I - 1
      IF (D(I).EQ. 0)GOTO 201
      DO 202 J = 1,M
      S = 0
      DO 203 K = 1,M
  203 S = S + X(I,K) * X(K,J)
      DO 202 K = 1,M
C-----THE FOLLOWING NONSENSE STATEMENT IS NEEDED TO COMBAT COMPILER BUG  6770   
      L=K+I-K
  202 X(K,J) = X(K,J) - S * X(K,L)
C-----END  J LOOP
  201 D(I) = X(I,I)
      X(I,I) = 1.0
      IF (M.LT. 1 ) GOTO 200
      DO 204 J = 1,M
      X(I,J) = 0
  204 X(J,I) = 0
  200 CONTINUE                                                           6780   
C-----END I
C-----DIAGONALIZATION  OF THE TRIDIAGONAL MATRIX
      DO 205 L = 1,N
      J = 0
      RH = EPS * (DABS (D(L)) +  DABS (E(L)))
      IF (RH .LE. B) GOTO 206
      B = RH
C-----TEST FOR SPLITTING
  206 DO 207 MM = L,N
      M = MM                                                             6790   
      IF (DABS (E(M)) .LT. B) GOTO 208
  207 CONTINUE
C-----LABEL 208 IS TEST FOR CONVERGENCE
  208 IF (M .EQ. L) GOTO 209
C-----LABEL 209 IS CONVERGENCE
C-----LABEL 210 IS SHIFT
  210 P =.5 *(D(L+1)- D(L))
      R = DSQRT (P * P + E(L)* E(L))
      RH = D(L) + P - R
      IF (P .LT. 0) RH=RH+R+R                                            6800   
      DO 213 I = L,N
  213 D(I) = D(I) - RH
      F = F + RH
C-----QR TRANSFORMATION
      J = J + 1
      P = D(M)
      C = 1.0
      S = 0
      MML = M - L
      DO 214 II = 1,MML                                                  6810   
      I = M - II
      R = DSQRT (P*P + E(I) * E(I))
      G = C * E(I)
      RH = C * P
      E(I+1) = S*R
      C = P/R
      S = E(I)/ R
      D(I+1) = RH + S * (C * G + S * D(I))
      DO 214 K = 1,N
      RH = X(K,I+1)                                                      6820   
      P = C * D(I) - S*G
      X(K,I+1) = X(K,I) * S + RH * C
  214 X(K,I) = X(K,I)  * C - RH * S
C-----END K LOOP AND END I LOOP
  215 E(L) = S * P
      D(L)= C * P
      IF (DABS(E(L)) .LT. B) GOTO 209
      IF (J .LT. ILIMIT) GOTO 210
      GOTO 1000
C-----1000 IS ERROR RETURN INDICATING FAILURE OF QR TO CONVERGE          6830   
  209 D(L) = D(L)+ F
  205 CONTINUE
C-----END L LOOP
C-----ORDERING OF EIGENVALUES
      DO 300 I = 1,N
      K = I
      P = D(I)
      IP1 = I + 1
      IF (N .LT. IP1) GOTO 301
      DO 302 J = IP1 , N                                                 6840   
      IF (D(J) .GE. P) GOTO 302
      K = J
      P = D(J)
  302 CONTINUE
C-----END J LOOP
  301 IF (K.EQ.I) GOTO 300
      D(K) = D(I)
      D(I) = P
      DO 304 J = 1,N
      P = X(J,I)                                                         6850   
      X(J,I) = X(J,K)
  304 X(J,K) = P
  300 CONTINUE
C-----END J,END K NOT EQUAL I, I, SDIAG
      GOTO  1001
C
C     SPECIAL CASE FOR 1 X 1
C
 500  D(1) = A(1,1)
      X(1,1) = 1                                                         6860   
      RETURN
C
 1000 WRITE (6,1002)
 1002 FORMAT  (1H0,25HFAILURE OF QR TO CONVERGE)
 1001 RETURN
      END
C*ID* SIP      C109.PTOLEMY.FITTERS                         PER705  21:
C
C     -----------------------------------------------------------------
C                                                                        6870   
      SUBROUTINE SIP(A,I,J,B,K,L,ALPHA,SUM,M)
C
CDC      ALLOCLEVEL A,B
      INTEGER I,J,K,L,M,N
      REAL*8  A(1),B(1)
      REAL*8  ALPHA,SUM
C
C       THIS ROUTINE SUMS THE PRODUCTS OF ELEMENTS  A(I), A(J),
C     A(I+(N-1)*(J-I))  IN  A  AND THE CORRESPONDING ELEMENTS  B(K),
C     B(L),  B(K+(N-1)*(L-K))  IN  B, FOR N = 3,4,...,M, WITH  ALPHA     6880   
C
C       THE ARRAYS  A  AND  B  MUST BE OF DIMENSION AT LEAST
C     I+(M-1)*(J-I)  AND  K+(M-1)*(L-K)  RESPECTIVELY.
C
C     ----------------------------------------------------------------
C
      SUM = ALPHA + A(I) * B(K)
      IF (M .LE. 1) GO TO 110
      SUM = SUM + A(J) * B(L)
      IF (M .LE. 2) GO TO 110                                            6890   
C
      DO 100 N = 3, M
         SUM = SUM + A(I+(N-1)*(J-I))*B(K+(N-1)*(L-K))
  100 CONTINUE
C
  110 RETURN
C     ------ LAST CARD OF SUBROUTINE  SIP  ------
      END
C*ID* SV01A    C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE SV01A (M,N,VAR,IVDIM,IVPAR, W)                          6900   
C
C     10/11/75 - W ADDED TO ARG. LIST - S.P.
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL VAR,W
      DIMENSION VAR(IVDIM,IVDIM), W(1)
      MN = M + N
      IMIN = MN + N
      GO TO (1,5), IVPAR
  1   IMAX = N*(MN+2)                                                    6910   
      DO 2 I = 1,N
      IMIN = IMIN + 1
      DO 3 J = I,N
      KI = J - I
      SUM = 0.0
      DO 4 K = IMIN,IMAX,MN
      KK = K + KI
      SUM = SUM + W(K)*W(KK)
  4   CONTINUE
      VAR(I,J) = SUM                                                     6920   
      VAR(J,I) = SUM
  3   CONTINUE
  2   CONTINUE
      RETURN
  5   IBAS = (N+1) * (MN + 1)
      DO 6 I = 1,N
      JCUR = IMIN +  I
      DO 7 J = 1,N
      W(J) = 0.0
  7   CONTINUE                                                           6930   
      IDEL = IBAS
      DO 8 J = 1,N
      W(J) = W(J) + W(IDEL)*W(JCUR)
      IF (J.EQ.N) GO TO 8
      KMIN = J + 1
      KCUR = JCUR
      DO 9 K = KMIN,N
      IDEL = IDEL + 1
      KCUR = KCUR + MN
      W(J) = W(J) + W(IDEL) * W(KCUR)                                    6940   
      W(K) = W(K) + W(IDEL) * W(JCUR)
 9    CONTINUE
      IDEL = IDEL + 1
      JCUR = JCUR + MN
  8   CONTINUE
      DO 10 J = I,N
      KCUR  = IMIN + J
      SUM = 0.0
      DO 11 K = 1,N
      SUM = SUM + W(K) * W(KCUR)                                         6950   
      KCUR = KCUR + MN
 11   CONTINUE
      VAR(I,J) = SUM
      VAR(J,I) = SUM
  10  CONTINUE
  6   CONTINUE
      RETURN
      END
C*ID* SYMINV   C109.PTOLEMY.FITTERS                         PER705
      SUBROUTINE SYMINV(A,N,B,M,DETERM,NMAX)                             6960   
C
C     SYMMETRIC MATRIX INVERSION (OPERATES ON UPPER TRIANGLE)
C
C     PTOLEMY VERSION FOR USE BY LMCHOL
C
C       1)  ONLY ONE INHOMOGENEOUS TERM.
C       2)  PIVOT MADE A LOGICAL*1 ARRAY.
C       3)  DETERMINANT CALCULATION REMOVED; LAST PIVOT ELEMENT IS
C           RETURNED AS DETERMINANT.
C       4)  LIMITED TO  25 X 25  MATRICES.                               6970   
C
C     11/5/75 - ABOVE CHANGES FOR PTOLEMY
C     12/28/79 - CDC TO CNI - RPG
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL A,B
C
CIB      LOGICAL*1 PIVOT
      LOGICAL PIVOT                                                     cni
      DIMENSION A(NMAX,N), B(NMAX)                                       6980   
      COMMON /F152/ P(25), Q(25),  PIVOT(25)
C
C     INITIALIZATION
C
      DO 20 I=1,N
      PIVOT(I)=.FALSE.
   20 CONTINUE
C
C     GRAND LOOP
C                                                                        6990   
      DO 550 I=1,N
C
C     SEARCH FOR PIVOT ELEMENT
C
      DETERM=0.0
      DO 100 J=1,N
      IF (PIVOT(J)) GO TO 100
      TEMP=DABS(A(J,J))
      IF (TEMP.LT.DETERM) GO TO 100
      L=J                                                                7000   
      DETERM=TEMP
  100 CONTINUE
      DETERM=A(L,L)
C
C     IF PIVOT ELEMENT IS ZERO, PROGRAM FAILS
C
      IF (DETERM.EQ.0.0) RETURN
C
C     PREPARATION OF ELIMINATION STEP
C                                                                        7010   
      PIVOT(L)= .TRUE.
      P(L)=1.0
      Q(L)=1./DETERM
      A(L,L)=0.0
      QB = B(L)
      B(L)=0.0
  120 CONTINUE
C
  130 IF (L.EQ.1) GO TO 230
      L1=L-1                                                             7020   
      DO 200 J=1,L1
      P(J)=A(J,L)
      Q(J)=P(J)*Q(L)
      IF ( .NOT. PIVOT(J) ) Q(J)=-Q(J)
  200 A(J,L)=0.0
C
  230 IF (L.EQ.N) GO TO 400
      L1=L+1
      DO 300 K=L1,N
      P(K)=-A(L,K)                                                       7030   
      Q(K)=P(K)*Q(L)
      IF ( .NOT. PIVOT(K) ) P(K)=-P(K)
  300 A(L,K)=0.0
C
C     ELIMINATION PROPER
C
  400 DO 550 J=1,N
      DO 450 K=J,N
  450 A(J,K)=A(J,K)+P(J)*Q(K)
  500 B(J)=B(J)+Q(J)*QB                                                  7040   
C
  550 CONTINUE
C
      RETURN
      END
C*ID* UPQPAR   C109.PTOLEMY.FITTERS                         PER705  21:
C
C     -------------------------------------------------------------
C
      SUBROUTINE UPQPAR(N,SQS,RJTJ,ACTRED,PRERED,DIRDRV,QPARM,QCPARM)    7050   
C
CDC      ALLOCLEVEL SQS,RJTJ
      INTEGER I,II,J,N
      REAL*8  SQS(N),RJTJ(N,N)
      REAL*8  ACTRED,PRERED,DIRDRV,QPARM,QCPARM
      REAL*8  DABS,FNORM,TRACE,Y,Z
C
C       THIS ROUTINE UPDATES THE MARQUARDT PARAMETER  QPARM  BASED ON
C     THE DEGREE OF AGREEMENT BETWEEN THE ACTUAL AND PREDICTED REDUCTION
C     IN THE SUM OF SQUARES  FSQ.                                        7060   
C
C     ON INPUT
C
C        N  IS THE DIMENSION OF THE VECTOR  SQS  AND THE SQUARE  MATRIX
C           IN  RJTJ.
C
C        SQS  CONTAINS  SQRT(S(I)), I = 1,2,...,N, WHERE THE  S(I)  ARE
C           ELEMENTS OF THE MARQUARDT SCALE CONTROL VECTOR.
C
C        RJTJ  CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF THE         7070   
C           MARQUARDT MATRIX  ((J TRANSPOSE) * J + QPARM * S).
C
C        ACTRED  IS THE ACTUAL REDUCTION IN THE SUM OF SQUARES  FSQ.
C
C        PRERED  IS THE PREDICTED REDUCTION IN THE SUM OF SQUARES  FSQ.
C
C        DIRDRV  IS THE DIRECTIONAL DERIVATIVE ALONG THE CURRENT
C           CORRECTION VECTOR.
C
C        QPARM  IS THE MARQUARDT PARAMETER.                              7080   
C
C        QCPARM  IS THE LOWER CUT-OFF LIMIT ON THE MARQUARDT
C           PARAMETER.
C
C     ON OUTPUT
C
C        QPARM  IS THE UPDATED MARQUARDT PARAMETER.
C
C        QCPARM  IS THE UPDATED LOWER CUT-OFF LIMIT ON THE
C           MARQUARDT PARAMETER.                                         7090   
C
C     -------------------------------------------------------------
C
C     TEST WHETHER THE RATIO OF THE ACTUAL REDUCTION TO THE PREDICTED
C       REDUCTION IS .GE. 0.25
C
  100 IF (ACTRED .GE. 0.25 * PRERED) GO TO 160
C
C     ACTRED / PRERED  .LT. 0.25. CALCULATE  NU  SUCH THAT
C       2 .LE.  NU  .LE. 10. IN THIS SECTION                             7100   
C
C       Y  -  CONTAINS 1 / NU
C
      Y = 0.5
      Z = 2.0D0 * DIRDRV - ACTRED
      IF (Z .GT. 0.0D0) Y = DIRDRV / Z
      IF (Y .GT. 0.5) Y = 0.5
      IF (Y .LT. 0.1D0) Y = 0.1D0
C
C     TEST WHETHER  QPARM = 0                                            7110   
C
      IF (QPARM .NE. 0.0D0) GO TO 150
C
C     QPARM = 0. CALCULATE  QCPARM = MAX(1/FNORM),(1/TRACE))  AND SET
C       QPARM = QCPARM,  Y = 2 * Y. IN THIS SECTION
C
C       FNORM  -  CONTAINS A SCALED INFINITY NORM OF  RJTJ  INVERSE
C       TRACE  -  CONTAINS A SCALED TRACE OF  RJTJ  INVERSE
C
      Y = 2.0D0 * Y                                                      7120   
      QPARM = 0.0D0
      TRACE = 0.0D0
C
      DO 140 I = 1, N
         TRACE = TRACE + RJTJ(I,I) * SQS(I) * SQS(I)
         FNORM = 0.0D0
C
         DO 110 J = 1, I
            FNORM = FNORM + DABS(RJTJ(J,I)) * SQS(J)
  110    CONTINUE                                                        7130   
C
         IF (I .EQ. N) GO TO 130
         II = I + 1
C
         DO 120 J = II, N
            FNORM = FNORM + DABS(RJTJ(I,J)) * SQS(J)
  120    CONTINUE
C
  130    FNORM = FNORM * SQS(I)
         IF (FNORM .GT. QPARM) QPARM = FNORM                             7140   
  140 CONTINUE
C
      IF (TRACE .LT. QPARM) QPARM = TRACE
      QPARM = 1.0D0 / QPARM
      QCPARM = QPARM
C
C     QPARM = NU * QPARM
C
  150 QPARM = QPARM / Y
      GO TO 200                                                          7150   
C
C     ACTRED / PRERED  .GE.  0.25. TEST WHETHER  ACTRED / PRERED
C       .LE.  0.75
C
  160 IF(ACTRED .LE. 0.75 * PRERED) GO TO 200
C
C     ACTRED / PRERED  .GT.  0.75. SET  QPARM = QPARM / 2. THEN SET
C       QPARM = 0 IF  QPARM  .LT.  QCPARM
C
      QPARM = QPARM * 0.5                                                7160   
      IF (QPARM .LT. QCPARM) QPARM = 0.0D0
  200 RETURN
C     ----- LAST CARD OF SUBROUTINE  UPQPAR  ----
      END
C*ID* VA02A    C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE VA02A ( M, N, F, X, TOL, TOL2, ESCALE, IPRINT, MAXFUN,
     1  IRET, W )
C
C     POWELL 1965 LEAST SQUARES FITTER FROM  AMD ARGONNE
C                                                                        7170   
C     M = NUMBER OF DATA POINTS TO BE FIT TO.
C     N = NUMBER OF PARAMETERS TO VARY.
C     F = COMPUTED ARRAY OF FUNCTION VALUES FOR EACH DATA POINT.
C         DIMENSION IS M.
C     X = PARAMETER ARRAY (DIMENSION N).
C     TOL AND TOL2 DETERMIN THE ACCURACY WITH WHICH THE MINIMUM IS
C     FOUND:  THE FITTER STOPS WHENEVER A STEP IS SUCH THAT
C     1)  A STEP IS SUCH THAT  DELTA(X(I)) < MAX( TOL*]X(I)], TOL2 )
C         FOR EVERY I;  OR
C     2)  THE SUM OF SQUARES FOR 4 CONSEQUTIVE ITERATIONS HAS            7180   
C         CHANGED BY LESS THAN TOL:
C              FSQ(K-4) - FSQ(K) < TOL*FSQ(K)
C     3)  THE SUM OF SQUARES BECOMES LESS THAN TOL2.
C
C     ESCALE - THE MAXIMUM SINGLE STEP IN X(I) WILL BE ESCALE*DELTA.
C           WHERE   DELTA = MAX( TOL*]X(I)], TOL2 )
C     IPRINT -  0 MEANS DON'T CALL MINMON FOR EACH ITERATION
C               1 MEANS DO CALL MINMON FOR EACH ITERATION
C     MAXFUN - MAXIMUM ALLOWED NUMBER OF CALLS TO CALFUN.
C     IRET - 0  FIT HAS COMPLETED O.K.                                   7190   
C          - 1  RAN OUT OF FUNCTION REFS.
C     W - A WORK ARRAY WHOSE DIMENSION MUST BE ATLEAST
C           (N+1)*(M+N+1) + (N*(N+1))/2 .
C
C      8/75  IRET ADDED, A FEW PRINT STATEMSNTS CHANGED. D.G.
C     10/11/75 - W ADDED TO ARG. LIST - S.P.
C     10/17/75 - MORE PRINTOUT AT ST. NO. 5; MORE COMMENTS - S.P.
C     11/6/75 - USE MINMON FOR PRINTING
C     1/22/76 - CHANGE IRET; ADD TOL AND TOL2; CALFUN IS A FUNCTION
C                                                                        7200   
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,F,W
      DIMENSION F(M), X(N), W(1)
C
      ITC = 0
      IRET=0
      MPLUSN=M+N
      KST=N+MPLUSN
      NPLUS=N+1
      KINV=NPLUS*(MPLUSN+1)                                              7210   
      KSTORE=KINV-MPLUSN-1
      FF = CALFUN (M,N,F,X)
      MC = 1
      NN=N+N
      K=NN
C
      DO 1 I=1,M
      K=K+1
      W(K)=F(I)
    1 CONTINUE                                                           7220   
C
      IINV=2
      K=KST
C
C     THIS IS THE START OF A LOOP THAT IS TERMINATED AT ST. NO. 24
C
      I=1
  2   DELTA = DMAX1( TOL*DABS(X(I)), TOL2 )
      X(I)=X(I)+DELTA
      FC = CALFUN (M,N,F,X)                                              7230   
      MC = MC+1
      X(I)=X(I)-DELTA
      DO 3 J=1,N
         K=K+1
         W(K)=0.
         W(J)=0.
    3 CONTINUE
      SUM=0.
      KK=NN
      DO 4 J=1,M                                                         7240   
         KK=KK+1
         F(J)=F(J)-W(KK)
         SUM=SUM+F(J)*F(J)
    4 CONTINUE
      IF (SUM) 5,5,6
5001  FORMAT('0**** ACCURACY FOR THE',I3, 'TH VARIABLE IS UNREASONABLY',
     1  ' SMALL'  /
     2  ' X(I), DELTA =', G20.10, G15.5, 10X, 'SUM =', G15.5 )
5     WRITE(6,5001) I, X(I), DELTA, SUM
      DO 8 J=1,M                                                         7250   
         NN=NN+1
         F(J)=W(NN)
    8 CONTINUE
      GO TO 10
    6 SUM=1./DSQRT(SUM)
      J=K-N+I
      W(J)=DELTA*SUM
      DO 9 J=1,M
         K=K+1
         W(K)=F(J)*SUM                                                   7260   
         KK=NN+J
         DO 11 II=1,I
            KK=KK+MPLUSN
            W(II)=W(II)+W(KK)*W(K)
   11    CONTINUE
    9 CONTINUE
      ILESS=I-1
      IGAMAX=N+I-1
      INCINV=N-ILESS
      INCINP=INCINV+1                                                    7270   
   12 IF (ILESS.NE.0) GO TO 14
      W(KINV)=1.
      GO TO 15
   14 B=1.
      DO 16 J=NPLUS,IGAMAX
         W(J)=0.
   16 CONTINUE
      KK=KINV
      DO 17 II=1,ILESS
         IIP=II+N                                                        7280   
         W(IIP)=W(IIP)+W(KK)*W(II)
         JL=II+1
         IF (JL-ILESS) 18,18,19
   18    DO 20 JJ=JL,ILESS
            KK=KK+1
            JJP=JJ+N
            W(IIP)=W(IIP)+W(KK)*W(JJ)
            W(JJP)=W(JJP)+W(KK)*W(II)
   20    CONTINUE
   19    B=B-W(II)*W(IIP)                                                7290   
         KK=KK+INCINP
   17 CONTINUE
      B=1./B
      KK=KINV
      DO 21 II=NPLUS,IGAMAX
         BB=-B*W(II)
         DO 22 JJ=II,IGAMAX
            W(KK)=W(KK)-BB*W(JJ)
            KK=KK+1
   22    CONTINUE                                                        7300   
         W(KK)=BB
         KK=KK+INCINV
   21 CONTINUE
      W(KK)=B
   15 GO TO (27,24),IINV
   24 I=I+1
      IF (I-N) 2,2,25
C
C     END OF THE LOOP THAT STARTED AT STNO 2
C                                                                        7310   
   25 IINV=1
      FF=0.
      KL=NN
      DO 26 I=1,M
         KL=KL+1
         F(I)=W(KL)
         FF=FF+F(I)*F(I)
   26 CONTINUE
      ISS=1
C                                                                        7320   
C     ALLOW THE USER TO PRINT STUFF FOR EACH ITERATION
C
 27   IF ( IPRINT .EQ. 1 )  CALL MINMON ( ITC, MC, 0, N, M, X, FF, F )
C
C     ARE WE DONE
C     SEE THE INITIAL COMMENTS FOR THE STOPPING CRITERIA
C
      IF ( ITC .LT. 4 )  GO TO 34
      IF ( CHANGE .LE. 1 )  GO TO 33
      IF ( F4-FF .GT. TOL*FF  .AND.  FF .GT. TOL2 )  GO TO 34            7330   
   10 CONTINUE
   33 RETURN
C
C     SAVE FUNCTION VALUES FROM FOUR PAST ITERATIONS
C
 34   F4 = F3
      F3 = F2
      F2 = F1
      F1 = FF
C                                                                        7340   
C
C     START OF A NEW ITERATION
C
      ITC=ITC+1
      K=N
      KK=KST
      DO 39 I=1,N
      K=K+1
      W(K)=0.
      KK=KK+N                                                            7350   
      W(I)=0.
      DO 40 J=1,M
      KK=KK+1
      W(I)=W(I)+W(KK)*F(J)
   40 CONTINUE
   39 CONTINUE
      DM=0.
      K=KINV
      DO 41 II=1,N
      IIP=II+N                                                           7360   
      W(IIP)=W(IIP)+W(K)*W(II)
      JL=II+1
      IF (JL-N) 42,42,43
   42 DO 44 JJ=JL,N
      JJP=JJ+N
      K=K+1
      W(IIP)=W(IIP)+W(K)*W(JJ)
      W(JJP)=W(JJP)+W(K)*W(II)
   44 CONTINUE
      K=K+1                                                              7370   
   43 IF (DM-DABS(W(II)*W(IIP))) 45,41,41
   45 DM=DABS(W(II)*W(IIP))
      KL=II
   41 CONTINUE
C
C     ACCUMULATE THE CHANGE VECTOR IN  W(I)  1 =< I =< N
C     AND THE MAXIMUM (OVER I) CHANGE/TOLERANCE
C
      II=N+MPLUSN*KL
      CHANGE=0.                                                          7380   
      DO 46 I=1,N
         JL=N+I
         W(I)=0.
         DO 47 J=NPLUS,NN
            JL=JL+MPLUSN
            W(I)=W(I)+W(J)*W(JL)
   47    CONTINUE
         II=II+1
         W(II)=W(JL)
         W(JL)=X(I)                                                      7390   
         DELTA = DMAX1( TOL*DABS(X(I)), TOL2 )
         CHANGE = DMAX1( CHANGE, DABS(W(I))/DELTA )
   46 CONTINUE
      DO 49 I=1,M
      II=II+1
      JL=JL+1
      W(II)=W(JL)
      W(JL)=F(I)
   49 CONTINUE
      FC=FF                                                              7400   
      ACC=0.1/CHANGE
      IT=3
      XC=0.
      XL=0.
      IS=3
      XSTEP=-DMIN1(0.50D0,ESCALE/CHANGE)
C
   51 CALL VD01A (IT,XC,FC,6,ACC,0.10D0,XSTEP)
      GO TO (52,53,53,53),IT
C                                                                        7410   
   52 MC=MC+1
      IF (MC-MAXFUN) 54,54,55
55    IRET=1
      ISS=2
      GO TO 53
   54 XL=XC-XL
      DO 57 J=1,N
      X(J)=X(J)+XL*W(J)
   57 CONTINUE
      XL=XC                                                              7420   
      FC = CALFUN (M,N,F,X)
      GO TO (59,59,60),IS
   60 K=N
      IF (FC-FF) 61,51,62
   61 IS=2
      FMIN=FC
      FSEC=FF
      GO TO 63
   62 IS=1
      FMIN=FF                                                            7430   
      FSEC=FC
      GO TO 63
   59 IF (FC-FSEC) 64,51,51
   64 K=KSTORE
      GO TO (75,74),IS
   75 K=N
   74 IF (FC-FMIN) 65,51,66
   66 FSEC=FC
      GO TO 63
   65 IS=3-IS                                                            7440   
      FSEC=FMIN
      FMIN=FC
   63 DO 67 J=1,N
      K=K+1
      W(K)=X(J)
   67 CONTINUE
      DO 68 J=1,M
      K=K+1
      W(K)=F(J)
   68 CONTINUE                                                           7450   
      GO TO 51
   53 K=KSTORE
      KK=N
      GO TO (69,70,69),IS
   70 K=N
      KK=KSTORE
   69 SUM=0.
      DM=0.
      JJ=KSTORE
      DO 71 J=1,N                                                        7460   
      K=K+1
      KK=KK+1
      JJ=JJ+1
      X(J)=W(K)
      W(JJ)=W(K)-W(KK)
   71 CONTINUE
      DO 72 J=1,M
      K=K+1
      KK=KK+1
      JJ=JJ+1                                                            7470   
      F(J)=W(K)
      W(JJ)=W(K)-W(KK)
      SUM=SUM+W(JJ)*W(JJ)
      DM=DM+F(J)*W(JJ)
   72 CONTINUE
      GO TO (73,10),ISS
   73 IF (ILESS.EQ.0) GO TO 83
      J=KINV
      KK=NPLUS-KL
      DO 76 I=1,KL                                                       7480   
      K=J+KL-I
      J=K+KK
      W(I)=W(K)
      W(K)=W(J-1)
   76 CONTINUE
      IF (KL-N) 77,78,78
   77 KL=KL+1
      JJ=K
      DO 79 I=KL,N
      K=K+1                                                              7490   
      J=J+NPLUS-I
      W(I)=W(K)
      W(K)=W(J-1)
   79 CONTINUE
      W(JJ)=W(K)
      B=1./W(KL-1)
      W(KL-1)=W(N)
      GO TO 88
   78 B=1./W(N)
   88 K=KINV                                                             7500   
      DO 80 I=1,ILESS
      BB=B*W(I)
      DO 81 J=I,ILESS
      W(K)=W(K)-BB*W(J)
      K=K+1
   81 CONTINUE
      K=K+1
   80 CONTINUE
   83 IF (FMIN.LT.FF) GO TO 82
      CHANGE=0.                                                          7510   
      GO TO 84
   82 FF=FMIN
      CHANGE=DABS(XC)*CHANGE
   84 XL=-DM/FMIN
      SUM=1./DSQRT(DMAX1(1.0D-16,SUM+DM*XL))
      K=KSTORE
      DO 85 I=1,N
      K=K+1
      W(K)=SUM*W(K)
      W(I)=0.                                                            7520   
   85 CONTINUE
      DO 86 I=1,M
      K=K+1
      W(K)=SUM*(W(K)+XL*F(I))
      KK=NN+I
      DO 87 J=1,N
      KK=KK+MPLUSN
      W(J)=W(J)+W(KK)*W(K)
   87 CONTINUE
   86 CONTINUE                                                           7530   
      GO TO 12
      END
C*ID* VD01A    C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE VD01A (ITEST,X,F,MAXFUN,ABSACC,RELACC,XSTEP)
C
C     THIS ROUTINE IS USED BY VA02A
C
C     11/6/75 - STOLLEN FROM AMDLIB
C
      implicit real*8 ( a-h, o-z )                                      implicit
      GO TO (1,2,2),ITEST
    2 IS=6-ITEST
      ITEST=1
      IINC=1
      XINC=XSTEP+XSTEP
      MC=IS-3
      IF (MC) 4,4,15
    3 MC=MC+1
      IF (MAXFUN-MC) 12,15,15
   12 ITEST=4                                                            7550   
   43 X=DB
      F=FB
      IF (FB-FC) 15,15,44
   44 X=DC
      F=FC
   15 RETURN
    1 GO TO (5,6,7,8),IS
    8 IS=3
    4 DC=X
      FC=F                                                               7560   
      X=X+XSTEP
      GO TO 3
    7 IF (FC-F) 9,10,11
   10 X=X+XINC
      XINC=XINC+XINC
      GO TO 3
    9 DB=X
      FB=F
      XINC=-XINC
      GO TO 13                                                           7570   
   11 DB=DC
      FB=FC
      DC=X
      FC=F
   13 X=DC+DC-DB
      IS=2
      GO TO 3
    6 DA=DB
      DB=DC
      FA=FB                                                              7580   
      FB=FC
   32 DC=X
      FC=F
      GO TO 14
    5 IF (FB-FC) 16,17,17
   17 IF (F-FB) 18,32,32
   18 FA=FB
      DA=DB
   19 FB=F
      DB=X                                                               7590   
      GO TO 14
   16 IF (FA-FC) 21,21,20
   20 XINC=FA
      FA=FC
      FC=XINC
      XINC=DA
      DA=DC
      DC=XINC
   21 XINC=DC
      IF ((D-DB)*(D-DC)) 32,22,22                                        7600   
   22 IF (F-FA) 23,24,24
   23 FC=FB
      DC=DB
      GO TO 19
   24 FA=F
      DA=X
   14 IF (FB-FC) 25,25,29
   25 IINC=2
      XINC=DC
      IF (FB-FC) 29,45,29                                                7610   
   29 D=(FA-FB)/(DA-DB)-(FA-FC)/(DA-DC)
      IF (D*(DB-DC)) 33,37,37
   37 D=0.5*(DB+DC-(FB-FC)/D)
      IF (DABS(D-X)-DABS(ABSACC)) 34,34,35
   35 IF (DABS(D-X)-DABS(D*RELACC)) 34,34,36
   34 ITEST=2
      GO TO 43
   36 IS=1
      X=D
      IF ((DA-DC)*(DC-D)) 3,26,38                                        7620   
   38 IS=2
      GO TO (39,40),IINC
   39 IF (DABS(XINC)-DABS(DC-D)) 41,3,3
   33 IS=2
      GO TO (41,42),IINC
   41 X=DC
      GO TO 10
   40 IF (DABS(XINC-X)-DABS(X-DC)) 42,42,3
   42 X=0.5*(XINC+DC)
      IF ((XINC-X)*(X-DC)) 26,26,3                                       7630   
   45 X=0.5*(DB+DC)
      IF ((DB-X)*(X-DC)) 26,26,3
   26 ITEST=3
      GO TO 43
      END
C*ID* VMMIN    C109.PTOLEMY.FITTERS                         PER705  21:
      SUBROUTINE VMMIN(DIAG,EPS,X,XP,G,GP,S,T,GB,Q1,Q2,H,FS,FJACOB,
     1  NA,N,IPRINT,NRAN,MAXFUN,ISIG)
C
C     AMD-DAVIDON VARIABLE-METRIC MINIMIZATION ADAPTED FOR LEAST-SQUARES 7640   
C
C     THIS IS A VERSION OF THE MAIN PROGRAM "MATMIN" TO BE FOUND IN
C     ARGONNE'S AMDLIB.  IT IS THE ORIGINAL DAVIDON VARIABLE METRIC
C     FITTER CONVERTED TO SUBROUTINE FORM AND PROVIDED WITH THE
C     ABILITY TO BASE THE INITIAL METRIC ON THE (JACOBIAN TRANSPOSE)*
C     JACOBIAN APPROXIMATION TO THE SECOND DERIVATIVE MATRIX THAT IS
C     APPROPRIATE TO SUM-OF-SQUARES FUNCTIONS.  THE STOPPING CRITERION
C     HAS ALSO BEEN CHANGED, AND THE ABILITY TO IMPOSE LINEAR
C     CONSTRAINTS ON THE PARAMETERS HAS BEEN REMOVED.  OTHERWISE
C     THE LOGIC OF THE PROGRAM IS UNCHANGED FROM                         7650   
C     THAT OF DAVIDON -- THE SUM-OF-SQUARES PROPERTY IS USED ONLY FOR
C     THE INITIAL METRIC.
C
C     IT IS POSSIBLE TO USE THIS SUBROUTINE TO FIND THE MINIMA
C     OF GENERAL FUNCTIONS BUT IN SUCH CASES DIAG MUST BE SET
C     TO A NON-ZERO NUMBER.
C
C     PARAMETER:
C
C     DIAG - CONTROLS THE INITIAL CHOICE OF METRIC:                      7660   
C        0 - METRIC = ( 2 JTJ )**(-1)  WHERE J IS THE JACOBIAN.
C       >0 - METRIC =  DIAG * UNIT MATRIX.
C       <0 - METRIC =  ]DIAG] * DELTA(I,J) * X(I)**2; I.E. A SCALED
C                      DIAGONAL MATRIX.
C
C     EPS - REQUIRED RELATIVE ACCURACY.  THE FITTER WILL STOP IF EITHER
C        OF THE FOLLOWING TWO CONDITIONS OBTAINS:
C        A)  ] X(I-1) - X(I) ] < EPS*]X(I)]  FOR ALL PARAMETERS
C        OR
C        B)  FUNC(I-2) - FUNC(I) < EPS*FUNC(I)                           7670   
C        HERE I LABELS THE ITERATION.
C
C     X - INITIAL AND FINAL PARAMETER VECTOR.  THE FIT STARTS AT THE
C        INPUT VALUE OF X AND X CONTAINS THE LOCATION OF THE MINIMUM
C        UPON RETURN.
C
C     XP, G, GP, S, T, GB, Q1, Q2 - WORK ARRAYS OF DIMENSION N.
C
C     H - AN N X N WORK ARRAY THAT WILL CONTAIN THE VARIANCE MATRIX
C            ( .5  D**2 F / DX DX ) ** (-1)                              7680   
C         UPON RETURN.
C
C     FS - WORK ARRAY OF DIMENSION NA
C     FJACOB - WORK ARRAY OF DIMENSION  NA X N
C
C     NA - NUMBER OF TERMS IN THE SUM-OF-SQUARES
C     N - NUMBER OF PARAMETERS THAT ARE TO BE VARIED IN THE FIT.
C
C     IPRINT - CONTROLS PRINTING DURING THE FIT:
C        0 - ONLY PRINT ERROR MESSAGES                                   7690   
C        1 - CALL MINMON AT EACH ITERATION TO ALLOW IT TO PRINT
C
C     NRAN - NUMBER OF RANDOM STEPS TO USE TO CONFIRM A MINIMUM
C
C     MAXFUN - MAXIMUM ALLOWED NUMBER OF FUNCTION REFERENCES.
C
C     ISIG - A RETURN CODE IS PLACED HERE:
C        0 - ALL WENT O.K.
C        1 - EXCEEDED MAXFUN, FIT HAS STOPPED WITH BEST VALUE FOR THE
C            MINIMUM IN X.                                               7700   
C
C
C     SUBROUTINES THAT MUST BE SUPPLIED BY THE USER:
C
C
C     F = CALGRA ( NA, N, FS, X, FJACOB, G )
C
C     THIS IS CALLED EACH TIME A FUNCTION AND ITS GRADIENT IS DESIRED.
C
C     NA AND N ARE AS ABOVE.                                             7710   
C     FS IS AS ABOVE AND SHOULD BE DEFINED TO HAVE THE INDIVIDUAL
C        TERMS IN THE SUM-OF-SQUARES.
C     F (THE FUNCTIONAL VALUE) IS THE SUM OF THE SQUARES OF THE FS.
C     X IS AS ABOVE AND IS THE LOCATION AT WHICH F, FS, FJACOB AND G
C       ARE TO BE EVALUATED.
C     FJACOB IS AS ABOVE AND MUST BE SET TO THE JACOBIAN:
C          FJACOB(I,J) = D(FS(I))/D(X(J))
C     G IS AS ABOVE AND MUST BE SET TO THE GRADIENT OF F:
C          G(J) = SUM(I) 2*FS(I)*FJACOB(I,J)
C                                                                        7720   
C     NOTE:  IT IS NECESSARY TO COMPUTE FS AND FJACOB ONLY ON THE
C            FIRST FUNCTION CALL.  IF DIAG IS NON-ZERO IT IS NEVER
C            NECESSARY TO COMPUTE FS AND FJACOB.  IN THE LATTER CASE
C            IT IS NOT NECESSARY THAT THE FUNCTION BE A SUM-OF-SQUARES
C            AND THE ARRAYS FS AND FJACOB AND THE PARAMETER NA WILL
C            NOT BE USED.
C
C
C     CALL MINMON ( IT, NCALL, NCALL, N, NA, X, F, FJACOB )
C                                                                        7730   
C     IF IPRINT = 1, THIS CALL IS MADE AT EACH ITERATION TO ALLOW
C     PRINTING TO SHOW THE PROGRESS OF THE FIT.  THE PARAMETERS ARE
C     AS ABOVE; THE FIRST THREE PARAMETERS ARE:
C
C     IT - THE ITERATION COUNTER - INCREMENTED BY 1 FOR EACH ITERATION.
C     NCALL - THE NUMBER OF FUNCITON REFERENCES THAT HAVE BEEN MADE.
C
C
C     1/23/76 - CONVERTED FROM MATMIN - M. H. MACFARLANE
C     2/19/76 - ADD COMMENTS, DIAG < 0 FEATURE - S. C. PIEPER            7740   
C
C
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,Q1,Q2,H,FS,FJACOB
C
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L
      COMMON/VMM2/IY
C
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),Q1(N),Q2(N),       7750   
     1H(N,N),FS(NA),FJACOB(NA,N)
C
      ISIG=0
      FAC= DABS(DIAG)
      E=EPS
      NSSW1=IPRINT
      NSSW2=0
      K=IABS(NRAN)
      NC=0
      DELTA=1                                                            7760   
      NCALL=0
      IY=1
C
      IF(FAC.EQ.0) GO TO 1070
C
C     HERE WE USE A DIAGONAL MATRIX AS THE INITIAL METRIC
C
      DO 35 II=1,N
      DO 35 JJ=1,N
      H(II,JJ)=0                                                         7770   
      IF (II-JJ) 35,30,35
   30 H(II,JJ)= FAC
      IF ( DIAG .LT. 0 )  H(II,JJ) = FAC*X(II)**2
   35 CONTINUE
 1070 M1=4
      GO TO 840
   65 M1=1
      IF(FAC.NE.0) GO TO 810
   70 MS=0
C                                                                        7780   
  100 IT=0
      F=CALGRA(NA,N,FS,X,FJACOB,G)
      NCALL=NCALL+1
      F1=F
      F2=0
      F3=0
      FFIRST=F
      IF(NSSW1.EQ.1) CALL MINMON(IT,NCALL,NCALL,N,NA,X,F,FJACOB)
      IF(FAC.NE.0) GO TO 840
      DO 1050 II=1,N                                                     7790   
      DO 1050 JJ=II,N
      H(II,JJ)=0
      DO 1040 IA=1,NA
 1040 H(II,JJ)=H(II,JJ)+FJACOB(IA,II)*FJACOB(IA,JJ)
      H(II,JJ)=2*H(II,JJ)
 1050 H(JJ,II)=H(II,JJ)
      CALL SYMINV(H,N,Q1,0,DET,N)
      DO 1060 II=2,N
      JTOP=II-1
      DO 1060 JJ=1,JTOP                                                  7800   
 1060 H(II,JJ)=H(JJ,II)
      GO TO 840
C
  120 CONTINUE
      M1=2
      IF(NCALL.LE.MAXFUN) GO TO 250
      ISIG=1
      GO TO 800
  250 CALL READY(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      IF(L-2)570,300,500                                                 7810   
  300 CALL AIM(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      IF(L.NE.1)GO TO 500
      CALL FIRE(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      IF(L.GT.2)GO TO 300
  500 CALL DRESS(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      GO TO 120
  570 IF (NSSW1) 600,580,600
  580 M1=5
      GO TO 840
  590 M1=2                                                               7820   
  600 CALL STUFF(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      IF(L.EQ.1)GO TO 100
C
 800  F = CALGRA( NA, N, FS, X, FJACOB, G )
      NCALL=NCALL+1
      M1=3
 810  IF ( M1 .EQ. 1 )  GO TO 70
  840 GO TO (120, 120, 880, 65, 590), M1
  880 CONTINUE
C                                                                        7830   
      DO 1200 I=1,N
      DO 1200 J=1,N
 1200 H(I,J)=2*H(I,J)
      RETURN
C
      END
C
C
      SUBROUTINE STUFF(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,FS,H,FJACOB
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L
      COMMON/VMM2/IY
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
C
      L=2
      K=K-1
      IF (K) 640,610,610                                                 7850   
C
  610 L=1
      PSQ=2*DSQRT(FP)
      P=DSQRT(PSQ)
      MS=MS+1
      write (6, 1) MS,DELTA,GS
      DO 620 I=1,N
      CALL RANDU (IY,YFL)
  620 T(I)=YFL-.5
      CALL MATMPY(N,N,H,T,S)                                             7860   
      CALL MATMPY(1,N,S,T,EL)
      EL=P/DSQRT(EL)
      DO 630 I=1,N
  630 X(I)=X(I)+EL*S(I)
C
  640 RETURN
  1   FORMAT ( '0RANDOM STEP', I4, '  DELTA =', D14.5,
     1  '  GS =', D14.5 )
      END
C                                                                        7870   
C
      SUBROUTINE RANDU(IX,YFL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 YFL
      IY=IX*65539
      IF (IY)5,6,6
    5 IY=IY+2147483647+1
    6 YFL=IY
      YFL=YFL*.4656613E-9
      IX=IY                                                              7880   
      RETURN
      END
C
C
      SUBROUTINE READY(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,FS,H,FJACOB
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),                   7890   
     1H(N,N),FS(NA),FJACOB(NA,N)
C
      L=1
      CALL MATMPY(N,N,H,G,S)
      DO 205 I=1,N
  205 S(I)=-S(I)
      CALL MATMPY(1,N,S,G,GS)
      DO 250 I=1,N
      IF(DABS(S(I)).GT.E*DABS(X(I))) GO TO 210
  250 CONTINUE                                                           7900   
      GO TO 240
C
  210 L=2
      EL=2.0
      T0=EL*F/GS
      IF (T0+EL) 213, 213, 212
  212 EL=-T0
  213 SL=-GS
      DO 215 I=1,N
  215 XP(I)=X(I)+EL*S(I)                                                 7910   
      FP=CALGRA(NA,N,FS,XP,FJACOB,GP)
      NCALL=NCALL+1
      IF(FP.GT.100*FFIRST) FP=10*FFIRST
      F3=F2
      F2=F1
      F1=FP
      IF((DABS(F1-F3).GT.E*DABS(F1)).OR.(FP.EQ.10*FFIRST)) GO TO 260
      L=1
      GO TO 240
  260 CALL MATMPY(1,N,S,GP,GSP)                                          7920   
      IF (-GSP) 240,240,220
  220 IF (F-FP) 240,240,225
C
  225 L=3
      FB=FP
      DO 230 I=1,N
      GB(I)=GP(I)
      T(I)=XP(I)
  230 CONTINUE
      IF (EL-2.0) 240,235,240                                            7930   
C
  235 L=4
      DELTA=DELTA+DELTA
      T0=1.0/SL
C
  240 RETURN
      END
C
C
      SUBROUTINE AIM(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)                  7940   
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,FS,H,FJACOB
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
C
      L=1
      Z=3.0/EL*(F-FP)+GS+GSP
      Q=DABS(Z*DSQRT(1.0-(GS/Z)*(GSP/Z)))                                7950   
      A=(Q-Z+GSP)/(Q+Q-GS+GSP)
      T0=EL/3.0*(Q+Q+Z+GSP)*A*A
      F0=FP-T0
      CALL MATMPY(N,N,H,GP,T)
      DO 305 I=1,N
  305 T(I)=(GSP/SL)*S(I)-T(I)
      CALL MATMPY(1,N,T,GP,GTP)
      IF (T0+T0+GTP) 315,310,310
  310 DO 312 I=1,N
  312 T(I)=XP(I)+A*(X(I)-XP(I))                                          7960   
      GO TO 340
  315 IF (F+F+GTP) 310,320,320
  320 DO 322 I=1,N
  322 T(I)=T(I)+XP(I)
      FB=CALGRA(NA,N,FS,T,FJACOB,GB)
      NCALL=NCALL+1
      IF(FB.GT.100*FFIRST) FB=10*FFIRST
      IF (F0-FB) 310,325,325
C
  325 L=3                                                                7970   
      DO 327 I=1,N
  327 S(I)=T(I)-XP(I)
      CALL MATMPY(1,N,S,GB,GTT)
      GTT=GTT-GTP
      IF (GTT) 340,330,330
C
  330 L=2
      GSS=GTT
      SL=-GTP
      EL=1.0                                                             7980   
C
  340 RETURN
      END
C
C
      SUBROUTINE FIRE(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,FS,H,FJACOB
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L       7990   
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
      EQUIVALENCE (TEMP,GTT)
C
      L=1
      TEMP=A/(1.0-A)
      FB=CALGRA(NA,N,FS,T,FJACOB,GB)
      NCALL=NCALL+1
      IF(FB.GT.100*FFIRST) FB=10*FFIRST
      CALL MATMPY(1,N,S,GB,GSB)                                          8000   
      T0=F
      IF (T0-FP) 403, 403, 402
  402 T0=FP
  403 IF (T0-FB+E*DABS(FB)) 415, 405, 405
  405 GSS=Q+Q
      T0=GSB*(TEMP-1.0/TEMP)
      IF (DABS(T0)-Q) 430,410,410
C
  410 L=2
      GO TO 440                                                          8010   
C
  415 L=3
      IF (FP-F) 425,420,420
C
 420  EL=(1.0-A)*EL
      FP=FB
      GSP=GSB
      DO 422 I=1,N
      XP(I)=T(I)
      GP(I)=GB(I)                                                        8020   
  422 CONTINUE
      GO TO 440
C
 425  EL=EL*A
      F=FB
      GS=GSB
      DO 427 I=1,N
      X(I)=T(I)
      G(I)=GB(I)
  427 CONTINUE                                                           8030   
      GO TO 440
C
  430 GSS=GSS+T0
      DO 435 I=1,N
  435 G(I)=(GB(I)-G(I))*TEMP+(GP(I)-GB(I))/TEMP
C
  440 RETURN
      END
C
C                                                                        8040   
      SUBROUTINE DRESS(NA,N,X,XP,G,GP,S,T,GB,FS,H,FJACOB)
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL X,XP,G,GP,S,T,GB,H,FS,FJACOB
      COMMON/VMM/F,GS,EL,SL,FP,GSP,T0,Z,Q,A,GSS,F0,GTP,FB,GTT,GSB,F1,
     1F2,F3,FFIRST,FAC,DELTA,E,P,NCALL,K,NC,NSSW1,NSSW2,MS,IT,M1,L
      DIMENSION X(N),XP(N),G(N),GP(N),S(N),T(N),GB(N),
     1H(N,N),FS(NA),FJACOB(NA,N)
C
      GO TO (505,520,530,525), L
C                                                                        8050   
  505 CALL MATMPY(N,N,H,G,X)
      CALL MATMPY(1,N,X,G,T0)
      IF (T0-GSS**2/SL-E*DABS(T0)) 515,510,510
  510 DO 512 II=1,N
      DO 512 JJ=1,N
  512 H(II,JJ)=H(II,JJ)-X(II)*X(JJ)/T0
      DELTA=DELTA*(EL*GSS/T0)
      T0=EL/GSS
      GO TO 525
C                                                                        8060   
  515 CONTINUE
C
  520 DELTA=DELTA*(EL*SL/GSS)
      T0=EL/GSS-1.0/SL
C
  525 DO 527 II=1,N
      DO 527 JJ=1,N
  527 H(II,JJ)=H(II,JJ)+T0*S(II)*S(JJ)
C
  530 IT=IT+1                                                            8070   
      F=FB
      DO 532 I=1,N
      G(I)=GB(I)
      X(I)=T(I)
  532 CONTINUE
      IF(NSSW1.EQ.1) CALL MINMON(IT,NCALL,NCALL,N,NA,X,F,FJACOB)
      RETURN
C
      END
C                                                                        8080   
C
      SUBROUTINE MATMPY(M,N,H,G,S)
      implicit real*8 ( a-h, o-z )                                      implicit
CDC      ALLOCLEVEL H,G,S
      REAL*8  H,G,S
      DIMENSION H(N,N),G(N),S(N)
C
      DO 720 II=1,M
      S(II)=0.0
      DO 720 JJ=1,N                                                      8090   
  720 S(II)=S(II)+H(JJ,II)*G(JJ)
C
      RETURN
      END
