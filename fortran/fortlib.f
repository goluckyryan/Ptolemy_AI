C*ID* CCNFRC   C109.PTOLEMY.FORTLIB                         PER704  09:
      SUBROUTINE CCNFRC ( NUMPTS, XS, YS )
C
C     THIS ROUTINE WILL INTERPOLATE OR EXTRAPOLATE A COMPLEX FUNCTION
C     THAT IS TABULATED ON A NONUNIFORM COMPLEX GRID.  A SPECIAL FORM OF
C     CONTINUED FRACTION IS USED.  NO TWO INPUT VALUES OF THE
C     FUNCTION MAY BE EQUAL.  THIS ALGORITHM IS QUITE SUCCESSFUL AT
C     FUNCTIONS AROUND POLES AND MODERATELY SUCCESSFUL AT CONTINUING
C     AROUND CUTS AND OTHER ESSENTIAL SINGULARITIES.
C                                                                          10   
C     ENTRY  CCNFRC  IS USED TO SETUP THE COEFFICIENTS OF THE COMPLEX
C     CONTINUED FRACTION.  THEN THE ENTRY  CCONTF  MAY BE USED AS
C     MANY TIMES AS IS DESIRED TO GET VALUES OF THE FUNCTION AT
C     NEW VALUES OF THE ORDINATE.  IN ADDITION THE CONVERGENCE
C     TO THE INTERPOLATED/EXTRAPOLATED VALUE IS AVAILABLE FOR
C     INSPECTION.
C
C
C     NUMPTS IS THE NUMBER OF INPUT VALUES OF XS AND YS.
C     XS IS THE ARRAY OF INPUT X'S AT WHICH THE FUNCTION VALUES ARE        20   
C        TABULATED.  THIS MUST BE A COMPLEX*16 ARRAY.
C        IT IS REARRANGED DURING THE CCNFRC CALCULATION, AND MUST BE
C        LEFT THAT WAY THROUGH ALL CCONTF CALLS.
C     YS, ON INPUT, IS THE ARRAY OF INPUT FUNCTION VALUES:
C            YS(I) = Y(XS(I))
C        IT MUST BE A COMPLEX*16 ARRAY.
C     ON OUTPUT, YS GETS THE COEFFICENTS.
C
C
C                                                                          30   
C     VERSION OF 4/18/70 = 70.108
C     SAVED ON DISK 7/17/70
C     9/26/72  ARGONNE VERSION FOR CONTIN
C     10/02/72   COMPLEX VERSION
C     5/1/74 - MADE BACK TO REAL VERSION
C     5/23/74 - AND BACK TO COMPLEX AGAIN.
C     8/11/75 - PTOLEMY VERSION BASED ON SPEAKEASY LINKULE VERSION.
C     7/23/79 - CHANGE ORDER, USE BIGGEST FOR BETTER ACCURACY, PASS
C        CCNFRC ARGUMENTS TO CCONTF - RPG
C     1/31/80 - END STATEMENT, COMPLEX*8 VERSION FOR VAX - RPG             40   
C     11/15/01 - SAVE NMAX added
C
      implicit real*8 ( a-h, o-z )                                      implicit
      COMPLEX*16  XS, YS, X, Y
C
      DIMENSION YS(NUMPTS), XS(NUMPTS)
C
      COMPLEX*16 D, XJ, YJ
C
      save nmax                                                            50   
C
CVA      REAL*4  DR
      DIMENSION  DR(2)
      EQUIVALENCE ( DR(1), D )
C     THIS IS THE SQUARE OF THE RELATIVE ACCURACY NEEDED IN THE
C     INTERPOLATED FUNCTION.
      DATA  TINY / 1.D-14 /
C
C
CVAC     REDUCE THE INPUT ARRAYS TO SINGLE PRECISION FOR THE VAX.          60   
CVA      CALL TOSING ( 2*NUMPTS, XS, XS )
CVA      CALL TOSING ( 2*NUMPTS, YS, YS )
      NMAX = NUMPTS - 1
C
C     SEARCH FOR THE LARGEST Y.
C
      COMP = -1
      DO 109 I = 1, NUMPTS
         D = YS(I)
         TEMP = DR(1)**2 + DR(2)**2                                        70   
         IF ( TEMP .LE. COMP )  GO TO 109
         COMP = TEMP
         K = I
 109  CONTINUE
      IF ( NMAX .LE. 0 )  RETURN
C
C     THE FIRST COEFFICIENT IS THE BIGGEST Y.
C     EXCHANGE THE FIRST WITH THE BIGGEST.
C
      YJ = YS(K)                                                           80   
      IF ( K .EQ. 1 )  GO TO 120
      YS(K) = YS(1)
      YS(1) = YJ
      XJ = XS(K)
      XS(K) = XS(1)
      XS(1) = XJ
C
C     CALCULATE THE CANDIDATES FOR THE SECOND COEFFICIENT,
C     AND SEARCH FOR THE LARGEST.
C                                                                          90   
 120  COMP = -1
      DO 139 I = 2, NUMPTS
         D = 1. - YJ/YS(I)
         YS(I) = D
         TEMP = DR(1)**2 + DR(2)**2
         IF ( TEMP .LE. COMP )  GO TO 139
         COMP = TEMP
         K = I
 139  CONTINUE
C                                                                         100   
C     LOOP THROUGH THE REST OF THE COEFFICIENTS.
C
      DO 400 J = 2, NUMPTS
C
C        EXCHANGE THE LARGEST CANDIDATE WITH THE J-TH, AND
C        FINISH CALCULATING THIS COEFFICIENT.
C
         XJ = XS(K)
         YJ = YS(K)
         IF ( K .EQ. J )  GO TO 209                                       110   
         YS(K) = YS(J)
         XS(K) = XS(J)
         XS(J) = XJ
 209     YJ = YJ/(XS(J-1)-XJ)
         YS(J) = YJ
         IF ( J .EQ. NUMPTS )  GO TO 450
C
C        QUIT IF THIS COEFFICIENT IS TINY.
C        NOTE:  THIS COMPARISON IS MEANINGFUL ONLY IF THE SIZE OF THE
C        DIFFERENCE XS(J-1)-XS(J) IS MORE OR LESS TYPICAL OF THE RANGE    120   
C        OF X.  IF THIS DIFFERENCE IS MANY ORDERS OF MAGNITUDE SMALLER
C        THAN TYPICAL X DIFFERENCES, IT MAY QUIT PREMATURELY, BUT THEN
C        THE WHOLE PROBLEM IS PATHOLOGICAL ANYWAY.
C
         IF ( COMP .LT. TINY )  GO TO 430
C
C        CALCULATE THE CANDIDATES FOR THE NEXT COEFFICIENT,
C        AND SEARCH FOR THE LARGEST.
C
         JP1 = J + 1                                                      130   
         COMP = -1
         XJ = XS(J-1)
         DO 309 I = JP1, NUMPTS
            D = 1. + YJ*(XS(I)-XJ)/YS(I)
            YS(I) = D
            TEMP = DR(1)**2 + DR(2)**2
            IF ( TEMP .LE. COMP )  GO TO 309
            COMP = TEMP
            K = I
 309     CONTINUE                                                         140   
C
 400  CONTINUE
C
C     END OF LOOP THROUGH COEFFICIENTS.
C     (IT NEVER GETS HERE ANYWAY.)
C
C     THE INPUT FUNCTION APPEARS TO BE A CONTINUED FRACTION EXACTLY
 430  NMAX = I-1
      write (6, 433)  NMAX, NUMPTS
 433  FORMAT ( '0**** WARNING:  CONTINUED FRACTION USED ONLY', I3,        150   
     1   ' OUT OF', I3, ' POINTS.' )
C
 450  CONTINUE
      RETURN
C
C    HERE WE COMPUTE ONE VALUE OF THE CONTINUED FRACTION
C
C      X IS THE  LOCATION OF THE DESIRED FUNCTION VALUE
C      Y WILL BE THE FUNCTION
C     BOTH X AND Y MUST BE COMPLEX*16.                                    160   
C
C
C
      ENTRY CCONTF ( NUMPTS, XS, YS, X, Y )
C
 500  Y = 0
      IF ( NMAX .LT. 1 )  GO TO 520
      DO 509  J=1, NMAX
         K = NMAX - J + 1
         Y = YS(1+K)*(X-XS(K))/(1.+Y)                                     170   
 509  CONTINUE
 520  Y = YS(1)/(1.+Y)
      RETURN
      END
CVAC
CVAC     THIS SUBROUTINE COMPRESSES A REAL*8 ARRAY INTO A REAL*4 ARRAY.
CVAC     THE ARGUMENTS MAY BE (AND ARE) THE SAME ARRAY.
CVAC
CVA      SUBROUTINE TOSING ( N, D, R )
CVA      REAL*8 D(N)                                                      180   
CVA      REAL*4 R(N)
CVA      DO 9  I = 1, N
CVA         R(I) = D(I)
CVA 9    CONTINUE
CVA      END
C*ID* CLEBSH        CLEBSCH collection for Ptolemy
c
c     This Ptolemy collection depends on the log(n!) table
c     being constructed and pointed to.
c                                                                         190   
      FUNCTION CLEBSH (A, B, X, Y, CIN, ZIN)
C
C     COMPUTES CLEBSCH-GORDAN AND 3-J COEFICENTS.
C
C       THIS ROUTINE COMPUTES THE CLEBSH GORDAN OR 3-J COEFICENT FOR
C     A GIVEN SET OF J'S AND M'S.  THE INPUT CONSISTS OF 6 INTEGERS THAT
C     ARE TWICE THE VALUES OF THE J'S AND M'S.  IN THIS MANNER
C     HALF-INTEGER ANGULAR MOMENTA CAN BE ACCOMMODATED.
C
C       A FACTORIZATION OF THE RACAH SUM SUGGESTED BY                     200   
C         J. G. WILLS, COMPUTER PHYSICS COMM. 2, 381 (1971)
C     IS USED.  THE REQUIRED FACTORIALS THAT MULTIPLY THE SUM ARE
C     FOUND IN ONE OF 3 WAYS:
C       FOR J1+J2+J3 < 96,        A  TABLE OF SQRT(FACTORIAL) IS USED;
C       FOR 95 < J1+J2+J3 < 1000, A TABLE OF LOG(FACTORIAL) IS USED;
C       FOR 1000 <= J1+J2+J3,     THE LOG GAMMA FUNCTION IS USED.
C
C       THESE CHOICES REPRESENT REASONABLE TRADE OFFS BETWEEN SIZE AND
C     SPEED CONSIDERATIONS.  BY RECOMPILING THE BLOCK DATA PROGRAMS
C     FOR THE /FACTRL/ OR /LOGFAC/ COMMON BLOCKS, ONE CAN CHANGE          210   
C     THE ABOVE BOUNDRIES TO EITHER REDUCE CORE (IF ONE IS NOT
C     INTERESTED IN LARGE J'S) OR INCREASE THE SPEED FOR LARGE J'S.
C
C       THE MAGNITUDE OF THE J'S AND M'S THAT CAN BE ACCOMMODATED BY
C     THIS ROUTINE IS DETERMINED BY CANCELLATIONS IN THE
C     RACAH SUM.  AS CAN BE SEEN IN THE TABLE BELOW, IF ALL THREE J'S
C     GET LARGE, SEVERE CANCELLATION SETS IN FOR EACH J OF THE ORDER OF
C     75, WHILE OVERFLOWS BEGIN TO OCCUR FOR THE J'S OF THE ORDER OF
C     150.  HOWEVER IF ONE OF THE J'S IS HELD SMALL, THE PRECISION
C     REMAINS GOOD FOR VERY LARGE VALUES OF THE OTHER TWO J'S AND THERE   220   
C     ARE NO OVERFLOW PROBLEMS (THE LARGEST VALUES TRIED WERE 2000 AS
C     IS SHOWN IN THE TABLE).
C
C       A SPECIAL CASE IS USED FOR M1 = M2 = M3 = 0.  ALL THREE J'S
C     MAY THEN BE VERY LARGE.
C
C       EACH LINE IN THE FOLLOWING TABLE IS THE RESULT OF SEVERAL
C     THOUSAND EVALUATIONS OF CLEBSH GORDAN COEFFICENTS FOR J'S AND
C     M'S IN THE RANGES INDICATED.  PRECISIONS WERE DETERMINED BY
C     COMPARING THE RESULTS TO A QUADRUPOLE PRECISION CLEBSCH GORDAN      230   
C     ROUTINE.  IT IS NOT SIGNIFICANT THAT J3 WAS HELD SMALL; ANY
C     OF THE J'S MAY BE KEPT SMALL WITH THE SAME TIMING AND RELATIVE
C     ERROR RESULTS.  THE EXECUTION TIMES ARE GIVEN IN MICROSECONDS.
C
C       J1, J2    J3   M1, M2     TIMES    RMS ABSOLUTE
C        RANGE   RANGE  LIMIT   /195  /75     ERROR
C        0   5   0   5    5      40   340     1E-16
C        5  10   0  10   10      48   425     8E-17
C       10  15   0  25   13      51   465     8E-17
C       15  30   0  47   30      59   560     7E-17                       240   
C       33  42   0  75   42      84   715     9E-13
C       50  52   0 100   52      95   835     2E-10
C       74  75   0 150   75     103   935     2E-5
C        100     0 200  100     124  1175     .6
C        150     0 300  150     205  2000     4E+8   SOME OVERFLOWS
C        100     0  10   10      63   455     2E-15
C        400     0  10   10      63   455     1E-14
C        500     0  10   10     465  3500     6E-14
C        500     0  10  500     470           6E-14  NOTE LARGE M1, M2
C       2000     0  10   10     460           2E-13                       250   
C        0  10   0  10    0      37   185     8E-17  M1=M2=M3=0
C        100      100     0      37   210     5E-15
C        300      300     0      38   215     3E-12
C        500      500     0     172  1320
C       1000     1000     0     174           3E-12
C
C     NO EVALUATION OF CLEBSCH-GORDAN COEFICENTS THAT WERE TRIVIALLY
C     ZERO WERE INCLUDED IN THE ABOVE TIMINGS.  THE LONG EXECUTION
C     TIMES FOR THE  J1, J2 > 500  ENTRIES ARE DUE TO THE USE OF
C     THE LOG GAMMA FUNCTION.  THEY MAY BE REDUCED BY INCREASING THE      260   
C     SIZE OF THE  /LOGFAC/  TABLE.
C
C
C     THERE ARE THREE ENTRIES:
C
C     CLEBSH (2*J1, 2*J2, 2*M1, 2*M2, 2*J3, 2*M3)   RETURNS A
C                  DOUBLE PRECISION CLEBSCH-GORDAN COEFICENT.
C     THREEJ (2*J1, 2*J2, 2*J3, 2*M1, 2*M2, 2*M3)   RETURNS A
C                  DOUBLE PRECISION 3-J SYMBOL (M3 = -(M1+M2)).
Cnor4 C     COFCG (2*J1, 2*J2, 2*M1, 2*M2, 2*J3, 2*M3)    RETURNS A       270   
Cnor4 C                  SINGLE PRECISION CLEBSCH-GORDAN COEFICENT.
C
C     IN ALL CASES THE TRIANGLE INEQUALITIES, THE REQUIREMENTS THAT
C     EACH (J, M) SATISFY   |M| <= J,   AND THE REQUIREMENT THAT
C     M1+M2 = M3 (FOR 3-J M1+M2+M3 = 0)  ARE CHECKED.  IF THEY ARE
C     VIOLATED, 0 IS RETURNED.  HOWEVER NO CHECK IS MADE TO VERIFY
C     THAT EACH (J, M) PAIR IS EITHER 1/2 INTEGER OR INTEGER NOR
C     THAT J1+J2+J3 IS INTEGER.  THE USER MUST INSURE THAT
C        2*J+2*M IS EVEN   AND THAT  2*J1+2*J2+2*J3  IS EVEN.
C     IF THIS IS NOT SO, MEANINGLESS RESULTS WILL BE RETURNED.            280   
C
C       CORE REQUIREMENTS (BASE 10):
C
C          CLEBSH    4088 BYTES
C         /FACTRL/    798     MAY BE REDUCED
C         /LOGFAC/   8028     MAY BE REDUCED OR INCREASED
C          DLGAMA    1008     (NOT KNOWN FOR VAX/VMS VERSION)
C          TOTAL   13,922
C
C     ALSO "DSQRT", "DLOG", "DEXP" AND THE FORTRAN I/O PACKAGE ARE USED   290   
C     SINCE THESE ARE LIKELY TO BE USED BY OTHER PARTS OF YOUR
C     PROGRAM, THEY ARE NOT INCLUDED IN THE ABOVE (THE FORTRAN I/O
C     PACKAGE IS USED BY DLOG, DEXP, AND DSQRT).  THE  /FACTRL/
C     AND  /LOGFAC/  COMMON BLOCKS AND THE "DLGAMA" FUNCTION WILL BE
C     SHARED WITH "SIXJ" AND "RACAH" SHOULD YOU ALSO USE THEM.
C
C     IN THE VAX/VMS VERSION, THE NAG ROUTINE S14ABF IS USED WITHIN
C     THE DUMMY ROUTINE DLGAMA.  THIS MAY BE EASILY CHANGED TO USE
C     ANOTHER LIBRARY ROUTINE.  IF TRANSFERING THIS PROGRAM TO
C     ANOTHER SITE, IT IS IMPORTANT TO CHECK THAT S14ABF IS THE           300   
C     DOUBLE PRECISION VERSION (RATHER THAN S14ABE) IN THE LOCAL
C     NAG LIBRARY.  THE SOURCE OF DLGAMA IS AT THE END OF THIS FILE.
C                                                  R.OSBORN 4/10/88
C     DLGAMA is standard on RS6000 so we do not use the above
C
C     S. PIEPER
C
C     4/5/71
C     8/21/72  ARGONNE VERSION IN DOUBLE PRECISION
C     8/25/72  INTEGER INPUT VERSION                                      310   
C     9/3/72  SPEEDED UP VERSION, ALSO QUITE ACCURATE
C     9/7/72  USE SQRT OF FACTORIALS AND BINOMIAL COFS
C     9/8/72  SINGLE AND DOUBLE PRECISION ENTRY POINTS
C     10/22/72  ALMOST COMPLETE REWRITE FOR RACAH FORMULA
C     10/23/72  ADD SPECIAL CASE FOR M1 = M2 = M3 = 0 AND 3-J ENTRY
C     5/17/73  CORRECT SIGN ERROR IN THREEJ
C     6/18/74 - TOTALLY NEW METHOD
C     7/8/74 - USE "XL" FUNCTIONS FOR SPEED, REDUCED SIZE.
C     9/11/74 - CORRECT THREEJ DOUBLE/SINGLE PRECISION ERROR
C     4/8/88 - MODIFIED FOR USE ON VAX/VMS - R.OSBORN AND G.L.GOODMAN     320   
C     7/3/91 - made to RS6000 - s.p.
C
C
      IMPLICIT INTEGER*4 (A-C, I-Z),  REAL*8 (D-H)
C
      REAL*8  SUM,ANSWER
      REAL*8  CLEBSH, THREEJ, THRJ
      REAL*8  DLGAMA
C
      COMMON /FACTRL/ MAXFAC, ISPACE, FACTBL(1)                           330   
C
      REAL*8 LF
      COMMON /LOGFAC/  MAXLF, lfbias, LF(1)
C
C     FOLLOWING ARE THE VAX FORTRAN FUNCTIONS THAT WE USE.
C     BTEST CHECKS TO SEE IF A BIT IS ON.
C     BTEST(I,0) CHECKS FOR I BEING ODD, BTEST(I,1) CHECKS FOR I/2 BEING
C     ISHFT(I,-1) SHIFTS THE BITS IN I ONE TO THE RIGHT  = /2
C
C     BTEST & ISHFT are the same on the RS6000!                           340   
C
      LOGICAL DBLSW, THRESW, BTEST
Cnor4       REAL*4 COFCG
C
C
C     FOLLOWING ARRAYS ARE DEFINED TO ALLOW USE OF LOOPS TO SHORTEN CODE
C
      DIMENSION    DRAY(9)
C
      EQUIVALENCE  (DRAY(1), DA1), (DRAY(2), DA2), (DRAY(3), DA3),        350   
     1  (DRAY(4), DY1), (DRAY(5), DX2), (DRAY(6), DY2), (DRAY(7), DX1),
     2  (DRAY(8), DZ1), (DRAY(9), DZ2)
C
C     FOLLOWING DUMMY EXTERNAL CAUSES THE FACTORIAL TABLES TO BE DEFINED
C
      EXTERNAL  FACDUM, LOGDUM
C
C
C      CLEBSH  IS THE AMD COMPATIBLE ENTRY POINT
C                                                                         360   
C
C
      THRESW = .FALSE.
      DBLSW = .TRUE.
      GO TO 10
C
Cnor4 C     COFCG IS SINGLE PRECISION ENTRY
Cnor4 C
Cnor4       ENTRY COFCG (A, B, X, Y, CIN, ZIN)
Cnor4       THRESW = .FALSE.                                              370   
Cnor4       DBLSW = .FALSE.
C
 10   Z = ZIN
C
 20   C = CIN
C
C
      IF (Z .EQ. X+Y)  GO TO 120
C
  110 CLEBSH=0.0                                                          380   
      RETURN
C
C
C
C     X1 = J1 + M1 ;  ETC.
C
C     THE FOLLOWING ISHFT FUNCTIONS EFFECT A DIVIDE BY 2
C     TO GET FROM THE 2*J STUFF TO J STUFF SINCE THE 3
C     (J+M)'S MUST BE INTEGERS.
C                                                                         390   
 120  X1 = ISHFT(A+X, -1)
      Y1 = ISHFT(B+Y, -1)
      Z1 = ISHFT(C+Z, -1)
C
C     THE  ISHFT  INSTRUCTION IS A LOGICAL SHIFT.  THUS IT CONVERTS
C     A NEGATIVE NUMBER INTO A LARGE POSITIVE NUMBER WHICH WILL LOUSE
C     UP THE TRIANGLE CONDITION CHECKS.  HERE WE CHECK TO SEE IF ANY OF
C     THE THREE NUMBERS HAD SIGN BITS THAT WERE SHIFTED INTO THE FIRST
C     POSITION.  IF SO THEN THE CORRESPONDING M WAS GREATER THAN J SO
C     WE JUST RETURN 0.                                                   400   
C     ALL OF THIS IS TO AVOID DIVIDES BY 2 WHICH ARE VERY SLOW ON
C     THE IBM/195.
C
C     T30 IS THE LOCATION OF THE SHIFTED SIGN BIT AND IS 30 FOR THE
C     VAX FORTRAN USED HERE
      T30 = 30
C
      IF (BTEST(X1,T30).OR.BTEST(Y1,T30).OR.BTEST(Z1,T30)) GO TO 110
C
C     X2 = J1 - M1 ;  ETC                                                 410   
C
      X2 = X1 - X
      Y2 = Y1 - Y
      Z2 = Z1 - Z
C
C     R = J1+J2+J3 + 1
C
      R = X1 + Y1 + Z2 + 1
C
C     A(I) = J1+J2+J3 - 2J(I) ;  A1 = -J1 + J2 + J3 , ETC.                420   
C
      A1 = R-1 - A
      A2 = R-1 - B
      A3 = R-1 - C
C
C     B1 = J3 - J1 - M2 ;   B2 = J3 - J2 + M1
C
      B1 = A1 - Y1
      B2 = A2 - X2
C                                                                         430   
C     WE WILL NEED THESE FLOATING POINT NUMBERS AT VARIOUS PLACES
C
      DA1 = A1
      DA2 = A2
      DA3 = A3
      DY1 = Y1
      DX2 = X2
C
      DB1 = DA1 - DY1
      DB2 = DA2 - DX2                                                     440   
C
      SUM = 1
C
      NMIN = MAX0 (0, -B1, -B2)
      NMAX = MIN0 (A3, X2, Y1)
C
C     THIS ONE CHECK MAKES SURE OF THE FOLLOWING
C     1)  ALL 3 TRIANGLE RULES
C     2)  |M|  <=  J  FOR ALL 3 CASES
C     3)  J  >=  0    FOR ALL 3 CASES                                     450   
C
      IF (NMAX - NMIN)  110, 300, 200
C
C     IS IT THE ALL M = 0 CASE
C
 200  IF (IABS(X)+IABS(Y) .EQ. 0)  GO TO 500
C
C     PREPARE TO DO THE SUM
      DX = NMAX
      E1 = -DX                                                            460   
      E2 = -DB1-DX
      E3 = -DB2 - DX
      E4 = DA3 - (DX-1)
      E5 = DX2 - (DX-1)
      E6 = DY1 - (DX-1)
C
      F1 = E1*E2*E3
      F4 = E4*E5*E6
C
      NMAXM2 = NMAX - 2                                                   470   
      SUM = 1 + F4/F1
      IF (NMAXM2 .LT. NMIN)  GO TO 300
C
      G1 = (E1+1)*(E2+1)*(E3+1) - F1
      G4 = (E4+1)*(E5+1)*(E6+1) - F4
      H1 = (E1+2)*(E2+2)*(E3+2) - G1 - G1 - F1
      H4 = (E4+2)*(E5+2)*(E6+2) - G4 - G4 - F4
C
C     COMPUTE THE SUM IN THE RACAH FORMULA
C                                                                         480   
      DO  259  N = NMIN, NMAXM2
         F1 = F1 + G1
         F4 = F4 + G4
         G1 = G1 + H1
         G4 = G4 + H4
         H1 = H1 + 6
         H4 = H4 + 6
         SUM = 1 + SUM*(F4/F1)
 259  CONTINUE
C                                                                         490   
C     NOW WE NEED THE FACTORIALS OUTSIDE THE SUM
C
 300  IF (R .LE. MAXFAC)  GO TO 400
C
C
C    SOME OF THE FACTORIALS WILL BE TOO LARGE TO USE THE TABLE.
C     WE USE ONE OF TWO METHODS*
C
C     1)  FOR MODERATELY SIZED J'S WE USE A TABLE OF
C         LOG FACTORIALS.                                                 500   
C
C     2)  FOR TRUELY LARGE J'S WE USED THE LOG GAMMA FUNCTION
C         TO GET THE LOGS OF THE FACTORIALS.
C
C
      IF (R .GT. MAXLF)  GO TO 350
C
C     THIS IS THE MODERATELY LARGE REGION - USE TABLE
C
C                                                                         510   
      E = lf(lfbias+1+A3) - lf(lfbias+1+NMIN) - lf(lfbias+1+A3-NMIN)
     1  + lf(lfbias+1+A1) - lf(lfbias+1+B1+NMIN) - lf(lfbias+1+Y1-NMIN)
     2  + lf(lfbias+1+A2) - lf(lfbias+1+B2+NMIN) - lf(lfbias+1+X2-NMIN)
      E2 = .5*(-lf(lfbias+1+A1) - lf(lfbias+1+A2) - lf(lfbias+1+A3)
     &   + lf(lfbias+1+X1) + lf(lfbias+1+X2)
     1    + lf(lfbias+1+Y1) + lf(lfbias+1+Y2) + lf(lfbias+1+Z1)
     &   + lf(lfbias+1+Z2) - lf(lfbias+1+R)
     2    + lf(lfbias+2+C) - lf(lfbias+1+C) )
C
      IF (DABS(E) .GT. 40)  GO TO 330                                     520   
C
C     SUM*E**2  MUST BE AN INTEGER.
C     SINCE E IS NOT VERY LARGE IS IS LIKELY TO BE MEANINGFUL TO
C     ATTEMPT TO FORCE THE INTEGER RESULT.
C
      SUM = ANINT(SUM*DEXP(E))
C
      ANSWER = DEXP(E2)*SUM
      GO TO 450
C                                                                         530   
C     HERE  E  IS LARGE AND AN ATTEMPT TO INTERGERIZE  SUM*EXP(E)  WOUL
C     NOT BE SIGNIFICANT ON A 16 PLACE MACHINE.  INSTEAD WE SPREAD
C     THE E FACTOR AMOUNGST THE OTHER STUFF IN AN ATTEMPT TO
C     AVOID OVERFLOWS.
C
 330  E = E + E2
      ANSWER = DEXP(E)*SUM
      GO TO 450
C
C     FOR LARGE VALUES OF J, WE USE THE LOG GAMMA FUNTION TO DIRECTLY     540   
C     FIND THE LOG OF THE COMPLETE MESS OF OUTSIDE FACTORIALS.
C
C     HERE WE AVOID INTEGER TO DOUBLE CONVERSIONS
C
 350  DX1 = DA2 + DA3 - DX2
      DC = C
      DN = NMIN
      DZ1 = DX1 + DY1 - DA3
      DY2 = DA1 + DA3 - DY1
      DZ2 = DC - DZ1                                                      550   
C
      E = DLOG(DC+1) - DLGAMA(2+DX1+DY1+DZ2)
      DO 369  I = 1, 9
         E = E + DLGAMA(1+DRAY(I))
 369  CONTINUE
      E = .5*E
      DA1 = DB1 + (DN+1)
      DA2 = DB2 + (DN+1)
      DA3 = DA3 - (DN-1)
      DY1 = DY1 - (DN-1)                                                  560   
      DX2 = DX2 - (DN-1)
      DY2 = (DN+1)
      DO  379  I = 1, 6
         E = E - DLGAMA(DRAY(I))
 379  CONTINUE
C
      ANSWER = DEXP(E) * SUM
C
      GO TO 450
C                                                                         570   
C     WE CAN USE THE FACTORIAL TABLES AND NOT WORRY
C
 400  F1 = FACTBL(1+A1)
      F2 = FACTBL(1+A2)
      F3 = FACTBL(1+A3)
C
C     HERE WE FORCE THE SUM OF PRODUCTS OF 3 BINOMIAL
C     COEFFICENTS TO BE AN INTEGER
C
      E = (F3/(FACTBL(1+NMIN)*FACTBL(1+A3-NMIN)))                         580   
     1   * (F2/(FACTBL(1+X2-NMIN)*FACTBL(1+B2+NMIN)))
     2   * (F1/(FACTBL(1+Y1-NMIN)*FACTBL(1+B1+NMIN)))
      SUM = ANINT(E*SUM*E)
C
C    FOR THE FOLLOWING NOTE THAT  R = J1+J2+J3 + 1
      ANSWER = SUM * (FACTBL(1+X1)*FACTBL(1+X2)/(F3*F2))
     1  * (FACTBL(1+Y1)*FACTBL(1+Y2)/F1)
     2  * (FACTBL(2+C)/FACTBL(1+C))
     3  * (FACTBL(1+Z1)*FACTBL(1+Z2)/FACTBL(1+R))
C                                                                         590   
 450  IF (BTEST(NMIN,0))  ANSWER = -ANSWER
      GO TO 800
C
C
C
C     USE SPECIAL CASE FOR M1 = M2 = M3 = 0 SINCE RACAH IS
C     NOT TRIVIAL FOR THAT CASE
C
C
C     RESULT IS 0 IF J1+J2+J3 IS ODD                                      600   
 500  IF (.NOT.BTEST(R,0))  GO TO 110
      JSUM = R-1
      IF (R .GT. MAXFAC)  GO TO 550
C
C     USE THE FACTORIAL TABLES
C
      ANSWER = ((FACTBL(1+A1)*FACTBL(1+A2))/FACTBL(1+R))
     1  * ((FACTBL(2+C)*FACTBL(1+A3))/FACTBL(1+C))
     2  * (FACTBL(1+JSUM/2)/
     3     (FACTBL(1+A1/2)*FACTBL(1+A2/2)*FACTBL(1+A3/2)) )**2            610   
      GO TO 570
C
C     USE LOG(FACTORIAL)
C
C     IN THE FOLLOWING WE USE  ISHFT( ,-1)  TO MAKE A /2.
C
 550  IF (R .GT. MAXLF)  GO TO 560
      E = .5 * (lf(lfbias+1+A1) + lf(lfbias+1+A2) + lf(lfbias+1+A3)
     1    + lf(lfbias+2+C) - lf(lfbias+1+C) - lf(lfbias+1+R) )
     2  - lf(lfbias+1+ISHFT(A1,-1)) - lf(lfbias+1+ISHFT(A2,-1))           620   
     3  - lf(lfbias+1+ISHFT(A3,-1)) + lf(lfbias+1+ISHFT(JSUM,-1))
C
      GO TO 565
C
 560  DJSUM = DA1 + DA2 + DA3
      E = DLOG(DFLOAT(C+1)) - DLGAMA(DJSUM+2)
      E2 = -DLGAMA(.5*DJSUM+1)
      DO 564  I = 1, 3
         E = E + DLGAMA(1+DRAY(I))
         E2 = E2 + DLGAMA(1 + .5*DRAY(I))                                 630   
 564  CONTINUE
      E = .5*E - E2
C
C
 565  ANSWER = DEXP(E)
C
C     GET PHASE
 570  IF (BTEST(A3,1))  ANSWER = -ANSWER
C
C     DETERMINE THE REQUESTED OUTPUT MODE                                 640   
C
C     WAS A 3-J REQUESTED
 800  IF (THRESW) GO TO 950
C
 810  IF (DBLSW)  GO TO 850
C
C   ROUND THE ANSWER TO SINGLE PRECISION
C
Cnor4       COFCG = REAL(ANSWER)
Cnor4       RETURN                                                        650   
C
C     DOUBLE PRECISION ENTRY USED
C
 850  CLEBSH = ANSWER
      RETURN
C
C     THIS IS AVAILABLE ONLY IN DOUBLE PRECISION
C
      ENTRY THRJ(A,B,CIN,X,Y,ZIN)
      ENTRY THREEJ (A, B, CIN, X, Y, ZIN)                                 660   
C
      Z = -ZIN
      THRESW = .TRUE.
      DBLSW = .TRUE.
      GO TO 20
C
C     NOW PUT IN THE 3-J FACTOR
C
 950  TEMP = A2-Z2
      IF (BTEST(TEMP, 0))  ANSWER = -ANSWER                               670   
      IF (C .GT. MAXFAC-1)  GO TO 970
      ANSWER = ANSWER * (FACTBL(1+C)/FACTBL(2+C))
      GO TO 810
 970  ANSWER = ANSWER / DSQRT(C+1.D0)
      GO TO 810
C
      END
      FUNCTION SIXJ (A, B, C, X, Y, Z)
C
C                                                                         680   
C     COMPUTES 6-J AND RACAH COEFFICIENTS.
C
C       THIS ROUTINE COMPUTES THE 6-J OR RACAH COEFFICIENT FOR
C     A GIVEN SET OF J'S.  THE INPUT CONSISTS OF 6 INTEGERS THAT
C     ARE TWICE THE VALUES OF THE J'S.  IN THIS MANNER
C     HALF-INTEGER ANGULAR MOMENTA CAN BE ACCOMODATED.
C
C       FOR THE FOLLOWING DISCUSSION CONSIDER THE 6-J COEFFICENT:
C
C           ( J1  J2  J3 )                                                690   
C           ( J4  J5  J6 )
C
C
C       A FACTORIZATION OF THE RACAH SUM SUGGESTED BY
C         J. G. WILLS, COMPUTER PHYSICS COMM. 2, 381 (1971)
C     IS USED.  THE REQUIRED FACTORIALS THAT MULTIPLY THE SUM ARE
C     FOUND IN ONE OF 3 WAYS:
C       FOR J1+J2+J4+J5 < 96,        A  TABLE OF SQRT(FACTORIAL) IS USE
C       FOR 95 < J1+J2+J4+J5 < 1000, A TABLE OF LOG(FACTORIAL) IS USED;
C       FOR 1000 <= J1+J2+J4+J5,     LOG GAMMA IS USED.                   700   
C
C       THESE CHOICES REPRESENT REASONABLE TRADE OFFS BETWEEN SIZE AND
C     SPEED CONSIDERATIONS.  BY RECOMPILING THE BLOCK DATA PROGRAMS
C     FOR THE /FACTRL/ OR /LOGFAC/ COMMON BLOCKS, ONE CAN CHANGE
C     THE ABOVE BOUNDRIES TO EITHER REDUCE CORE (IF ONE IS NOT
C     INTERESTED IN LARGE J'S) OR INCREASE THE SPEED FOR LARGE J'S.
C
C       THE MAGNITUDE OF THE J'S AND M'S THAT CAN BE ACCOMODATED BY
C     THIS ROUTINE IS DETERMINED BY CANCELLATIONS IN THE
C     RACAH SUM.  AS CAN BE SEEN IN THE TABLE BELOW, IF ALL SIX J'S       710   
C     GET LARGE, SEVERE CANCELLATION SETS IN FOR EACH J OF THE ORDER OF
C     100.  HOWEVER IF THREE OF THE J'S ARE HELD SMALL, THE PRECISION
C     REMAINS GOOD FOR VERY LARGE VALUES OF THE OTHER THREE J'S AND
C     THERE ARE NO OVERFLOW PROBLEMS (THE LARGEST VALUES TRIED WERE
C     500 AS SHOWN IN THE TABLE).
C
C
C       EACH LINE IN THE FOLLOWING TABLE IS THE RESULT OF SEVERAL
C     THOUSAND EVALUATIONS OF 6-J COEFFICENTS FOR J'S
C     IN THE RANGES INDICATED.  PRECISIONS WERE DETERMINED BY             720   
C     COMPARING THE RESULTS TO A QUADRUPOLE PRECISION 6-J
C     ROUTINE.  IT IS NOT SIGNIFICANT THAT J3, J4 AND J5 WERE HELD
C     SMALL; ANY THREE OF THE J'S MAY BE KEPT SMALL WITH THE SAME TIMIN
C     AND RELATIVE ERROR RESULTS.  THE TIMING INFORMATION
C     IS IN MICROSECONDS.
C
C
C       J1, J2, J6    J3, J4, J5        TIMES      RMS ABSOLUTE
C         RANGE         RANGE        /195   /75       ERROR
C         0   5         0   5         49    400       4E-17               730   
C         5  10         0  10         64    615       7E-18
C        10  15         0  25         59    550       5E-18
C        15  30         0  47         73    650       2E-16
C        33  42         0  75         95    790       2E-15
C        50  52         0 100        112   1020       3E-13
C          75           0 150        130   1140       3E-10
C         100           0 200        146   1315       5E-6
C         150            150         575   5520       3E+2
C
C          50           0  10         68    485       9E-16               740   
C         100           0  10         68    485       7E-16
C         400           0  10         70    478       6E-15
C         500           0  10        670   5080       9E-15
C
C
C     NO EVALUATION OF 6-J COEFICENTS THAT WERE TRIVIALLY
C     ZERO WERE INCLUDED IN THE ABOVE RESULTS.  THE LONG EXECUTION
C     TIMES IN THE LAST LINE OF THE TABLE ARE DUE TO THE USE OF THE
C     LOG GAMMA FUNCTION.  THEY MAY BE MADE COMPARABLE TO THE TIMES IN
C     THE PREVIOUS LINE BY INCREASING THE SIZE OF THE  /LOGFAC/           750   
C     TABLE.
C
C
C     THERE ARE THREE ENTRIES:
C
C     SIXJ (2*J1, 2*J2, 2*J3, 2*J4, 2*J5, 2*J6)   RETURNS THE
C                  DOUBLE PRECISION 6-J COEFICENT.
C                        ( J1  J2  J3 )
C                        ( J4  J5  J6 )
C                                                                         760   
C     RACAH (2*J1, 2*J2, 2*J5, 2*J4, 2*J3, 2*J6)   RETURNS THE
C                  DOUBLE PRECISION RACAH SYMBOL
C                     W (J1, J2, J5, J4; J3, J6)
C
Cnor4 C     COF6J (2*J1, 2*J2, 2*J3, 2*J4, 2*J5, 2*J6)    RETURNS A
Cnor4 C                  PROPERLY ROUNDED SINGLE PRECISION 6-J COEFICENT
C
C     IN ALL CASES THE TRIANGLE INEQUALITIES ARE CHECKED.  IF THEY ARE
C     VIOLATED, 0 IS RETURNED.  HOWEVER NO CHECK IS MADE TO VERIFY
C     THAT  J1+J2+J3,  J4+J2+J6,  J1+J5+J6,  AND  J4+J5+J3  ARE           770   
C     INTEGERS.  THE USER MUST INSURE THAT  2*J1 + 2*J2 + 2*J3,  ETC,
C     ARE EVEN.  IF THIS IS NOT SO MEANINGLESS RESULTS WILL BE RETURNED
C
C       CORE REQUIREMENTS (BASE 10):
C
C          SIXJ      4000 BYTES
C         /FACTRL/    798     MAY BE MADE SMALLER
C         /LOGFAC/   8028     MAY BE MADE SMALLER OR LARGER
C          DLGAMA    1008
C          TOTAL   13,834                                                 780   
C
C     IN ADDITION "DEXP" AND THE FORTRAN I/O PACKAGE WILL BE USED.
C     SINCE THESE ARE LIKELY TO BE USED BY OTHER PARTS OF YOUR
C     PROGRAM, THEY ARE NOT INCLUDED IN THE ABOVE (THE FORTRAN I/O
C     PACKAGE IS CALLED ONLY BY DEXP).  THE  /FACTRL/
C     AND  /LOGFAC/  COMMON BLOCKS AND THE "DLGAMA" FUNCTION WILL BE
C     SHARED WITH "CLEBSH" AND "THREEJ" SHOULD YOU ALSO USE THEM.
C
C
C     S. PIEPER                                                           790   
C
C      4/1/71  MODIVIED FROM COFJU TO GIVE 6-J
C     8/21/72  ARGONNE DOUBLE PRECISION VERSION
C     9/5/72  INTEGER VERSION -- S. PIEPER
C     9/6/72  USE SQRT OF BINOMIAL COFS AND FACTORIALS
C     9/8/72  SINGLE AND DOUBLE PRECISION ENTRY POINTS
C     9/23/73  SIXJ AND RACAH ENTRIES.
C     7/6/74 - MAJOR REWRITE - TOTALLY NEW METHOD
C     APRIL 1988 - MOVE TO VAX/VMS R.OSBORN AND G.L.GOODMAN
C     7/3/91 - RS6000 version - s.p.                                      800   
C
      IMPLICIT INTEGER*4 (A-C, I-Z),  REAL*8 (D-H)
C
      REAL*8 SUM, ANSWER, SIXJ, RACAH, DLGAMA
Cnor4       REAL*4  COF6J
      LOGICAL  DBLSW, RACHSW, BTEST
C
      COMMON /FACTRL/ MAXFAC, ISPACE, FACTBL(1)
C
      REAL*8 LF                                                           810   
      COMMON /LOGFAC/  MAXLF, lfbias, LF(1)
C
C
C     THE FOLLOWING FAKE EXTERNALS DRAG IN THE REQUIRED BLOCK DATAS
C     FOR THE ABOVE COMMON BLOCKS
C
      EXTERNAL  FACDUM, LOGDUM
C
CCCCC**********
CCCC      COMMON /NUMBLK/  NUMIT, IBIG, INOTSO, ILF, E, SUM               820   
CCCCC**********
C
C     WE USE THE FOLLOWING VAX FORTRAN FUNCTIONS:
C     BTEST - TESTS TO SEE IF A BIT IS ON.  BTEST(I, 0) IS TRUE
C            IF  I  IS ODD.
C     ISHFT(N,-1) SHIFTS THE BITS IN N ONE TO THE RIGHT  = /2
C
C
      DIMENSION DRAY(8)
      EQUIVALENCE  (DA1, DRAY(1)), (DA2, DRAY(2)),                        830   
     1  (DA3, DRAY(3)), (DA4, DRAY(4)),
     2  (DC, DRAY(5)),  (DZ, DRAY(7))
C
C
C    SIXJ IS THE DOUBLE PRECISON ENTRY POINT
C
      DBLSW = .TRUE.
      RACHSW = .FALSE.
      GO TO 100
C                                                                         840   
Cnor4 C    COF6J IS THE SINGLE PRECISION ENTRY POINT
Cnor4 C
Cnor4       ENTRY COF6J (A, B, C, X, Y, Z)
Cnor4       DBLSW = .FALSE.
Cnor4       RACHSW = .FALSE.
Cnor4       GO TO 100
C
C     RACAH ENTRY IS DOUBLE PRECISION
C
      ENTRY  RACAH (A, B, Y, X, C, Z)                                     850   
C
      RACHSW = .TRUE.
      DBLSW = .TRUE.
C
C
C
C
 100  SUM=1
C
CCCCC****************                                                     860   
CCCC      NUMIT = 0
CCCCC****************
C
C
      A2 = ISHFT(X+Y-C, -1)
      A3 = ISHFT(A+Y-Z, -1)
      A4 = ISHFT(B+X-Z, -1)
      A1 = ISHFT(A+B-C, -1)
C
C     THE ISHFT IS A LOGICAL SHIFT AND HENCE CONVERTS NEGATIVE            870   
C     QUANTITIES INTO VERY LARGE POSITIVE QUANTITIES.  IF ANY OF
C     THE AI'S ARE NEGATIVE, THEN THE CORRESPONDING TRIANGLE
C     IS VIOLATED.  WE TEST FOR THAT CASE HERE.
C
      IF (BTEST(A1,30) .OR. BTEST(A2,30) .OR.
     1    BTEST(A3,30) .OR. BTEST(A4,30)) GO TO 900
C
C     B1 = J3 +J6 - J1 - J4
C     B2 = J3 + J6 - J2 - J5
C     R  = J1 + J2 + J4 + J5 + 1                                          880   
      B1 = Y - A2 - A3
      B2 = A3 - A1 + Z - Y
C
      R = A1 + A2 + C + 1
C
      NMAX = MIN0 (A1, A2, A3, A4)
      NMIN = MAX0 (0, -B1, -B2)
C
C
CCCCC***************                                                      890   
CCCC      NUMIT = NMAX - NMIN + 1
CCCCC***************
C
C     CONVERT VARIOUS INTEGERS TO DOUBLE PRECISON FOR NOW AND
C     POSSIBLY LATER.
C
      DA1 = A1
      DA2 = A2
      DA3 = A3
      DA4 = A4                                                            900   
      DY = Y
      DZ = Z
      DC = C
      DB1 = DY - DA2 - DA3
      DB2 = DA3 - DA1 + DZ - DY
      DR = DA1 + DA2 + DC + 1
      DN = NMIN
      DX = NMAX
C
C     THIS ONE TEST CHECKS ALL 12 INEQUALITIES IMPLIED BY THE             910   
C     4 TRIANGLE CONDITIONS.
C
      IF (NMAX - NMIN)  900, 300, 200
C
 200  D1 = DA1 - (DX-1)
      D2 = DA2 - (DX-1)
      D3 = DA3 - (DX-1)
      D4 = DA4 - (DX-1)
      D5 = -DX
      D6 = -DB1 - DX                                                      920   
      D7 = -DB2 - DX
      D8 = DR - (DX-1)
C
      E1 = D1*D2*D3*D4
      E5 = D5*D6*D7*D8
C
      NMAXM2 = NMAX-2
      SUM = 1 + E1/E5
      IF (NMAXM2 .LT. NMIN)  GO TO 300
C                                                                         930   
      F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1
      F5 = (D5+1)*(D6+1)*(D7+1)*(D8+1) - E5
C
      G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1
      G5 = (D5+2)*(D6+2)*(D7+2)*(D8+2) - F5 - F5 - E5
C
      H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1
      H5 = (D5+3)*(D6+3)*(D7+3)*(D8+3) - 3*(G5+F5) - E5
C
C     COMPUTE THE SUM IN THE RACAH FORMULA                                940   
C
      DO 259  N = NMIN, NMAXM2
         E1 = E1 + F1
         E5 = E5 + F5
         F1 = F1 + G1
         F5 = F5 + G5
         G1 = G1 + H1
         G5 = G5 + H5
         H1 = H1 + 24
         H5 = H5 + 24                                                     950   
         SUM = 1 + SUM*(E1/E5)
 259  CONTINUE
C
C     NOW WE NEED THE OUTSIDE FACTORIALS
C     WE USE ONE OF 3 METHODS DEPENDING ON THE SIZE OF J1+J2+J4+J5
C        1)  A TABLE OF SQRT(FACTORIALS)
C        2)  A TABLE OF LN(FACTORIALS)
C        3)  LN(GAMMA) FUNCTION
C
 300  IF (R .GT. MAXFAC)  GO TO 400                                       960   
C
C     FIRST THE STUFF WE FACTORED OUT OF THE SUM
C
      E = FACTBL(1+R-NMIN) / (FACTBL(1+NMIN)
     1  * FACTBL(1+B1+NMIN) * FACTBL(1+B2+NMIN) * FACTBL(1+A1-NMIN)
     2  *  FACTBL(1+A2-NMIN) * FACTBL(1+A3-NMIN) * FACTBL(1+A4-NMIN) )
C
C     SUM IS AN INTEGER NOW
C
      SUM = ANINT(E*SUM*E)                                                970   
C
C     AND NOW THE OUTSIDE FACTORS
C
C
      DENOM = (FACTBL(2+A1+C)/(FACTBL(1+A1)*FACTBL(1+A4+B1)
     1                        *FACTBL(1+A3+B2)))
     2  * (FACTBL(2+A2+C)/(FACTBL(1+A2)*FACTBL(1+A3+B1)
     3                        *FACTBL(1+A4+B2)))
     4  * (FACTBL(2+A3+Z)/(FACTBL(1+A3)*FACTBL(1+A2+B1)
     5                        *FACTBL(1+A1+B2)))                          980   
     6  * (FACTBL(2+A4+Z)/(FACTBL(1+A4)*FACTBL(1+A1+B1)
     7                        *FACTBL(1+A2+B2)))
      ANSWER = SUM/(DENOM)
      GO TO 600
C
C     TRY FOR LOG FACTORIAL TABLE
C
 400  IF (R .GT. MAXLF)  GO TO 500
C
CCCCC**********                                                           990   
CCCC      IBIG = IBIG + 1
CCCC      ILF = ILF + 1
CCCCC************
C
      E = lf(lfbias+1+R-NMIN) - lf(lfbias+1+NMIN)
     1  - lf(lfbias+1+B1+NMIN) - lf(lfbias+1+B2+NMIN)
     &   - lf(lfbias+1+A1-NMIN)
     2  - lf(lfbias+1+A2-NMIN) - lf(lfbias+1+A3-NMIN)
     &   - lf(lfbias+1+A4-NMIN)
      E2 = .5 * ( lf(lfbias+1+A1) + lf(lfbias+1+A4+B1)                   1000   
     &   + lf(lfbias+1+A3+B2)
     1    - lf(lfbias+2+A1+C)
     2  + lf(lfbias+1+A2) + lf(lfbias+1+A3+B1) + lf(lfbias+1+A4+B2)
     &   - lf(lfbias+2+A2+C)
     3  + lf(lfbias+1+A3) + lf(lfbias+1+A2+B1) + lf(lfbias+1+A1+B2)
     &   - lf(lfbias+2+A3+Z)
     4  + lf(lfbias+1+A4) + lf(lfbias+1+A1+B1) + lf(lfbias+1+A2+B2)
     &   - lf(lfbias+2+A4+Z) )
C
C     AT THIS POINT IT WOULD BE POSSIBLE TO TRY TO INTERGERIZE           1010   
C     SUM*EXP(E).  HOWEVER, TESTS SHOW THAT IN PRACTICE THIS WOULD
C     BE USEFUL (RESULT < 10**18) LESS THAN 1% OF THE TIME SO
C     WE DON'T BOTHER WITH IT.
C
 450  ANSWER = DEXP(E+E2) * SUM
      GO TO 600
C
C     HERE THE J'S ARE VERY LARGE - WE USE LN(GAMMA)
C
 500  DRAY(6) = DC                                                       1020   
      DRAY(8) = DZ
      E2 = 0
      E = DLGAMA(1-DN+DR) - DLGAMA(1+DN)
     1  - DLGAMA(1+DN+DB1) - DLGAMA(1+DN+DB2)
      DO 549  I = 1, 4
         E2 = E2 + DLGAMA(1+DRAY(I)) + DLGAMA(1+DB1+DRAY(I))
     1      + DLGAMA(1+DB2+DRAY(I)) - DLGAMA(2+DRAY(I)+DRAY(4+I))
         E = E - DLGAMA(1-DN+DRAY(I))
 549  CONTINUE
      E = E + .5*E2                                                      1030   
C
CCCC      IBIG = IBIG + 1
CCCCC***************
CCCCC
      ANSWER = DEXP(E) * SUM
C
C     NOW GET THE SIGN CORRECT
C
 600  TEMP = NMIN + R
      IF (.NOT. BTEST(TEMP, 0))  ANSWER = -ANSWER                        1040   
C
C
C     CALCULATION OF THE 6-J COMPLETED; WAS RACAH REQUESTED
C
 700  IF (.NOT. RACHSW)  GO TO 800
C
C     YES, CONVERT THE SIGN
      IF (.NOT. BTEST(R, 0))  ANSWER = -ANSWER
C
C     WHAT PRECISION WAS REQUESTED FOR THE RESULT                        1050   
C
 800  IF (DBLSW)  GO TO 850
C
C    ROUND THE ANSWER TO SINGLE PRECISION
C
Cnor4       COF6J = REAL(ANSWER)
Cnor4       RETURN
C
C     DOUBLE PRECISION ENTRY POINT WAS USED
C                                                                        1060   
 850  SIXJ = ANSWER
      RETURN
C
C     THE TRIANGLE RULES ARE NOT SATISFIED
C
 900  COF6J = 0
      IF (DBLSW) SIXJ = 0
      RETURN
      END
      FUNCTION WIG9J ( J1, J2, J3, J4, J5, J6, J7,                       1070   
     1  J8, J9 )
C
C
C   COMPUTES 9-J COEFICENTS.
C
C
C     THIS ROUTNE COMPUTES THE WIGNER 9-J COEFICENT
C
C         (  J1  J2  J3  )
C         (  J4  J5  J6  )                                               1080   
C         (  J7  J8  J9  )
C
C   THE INPUT CONSISTS OF 9 INTEGERS THAT ARE TWICE THE CORRESPONDING
C   J'S.  IN THIS MANNER HALF-INTEGER ANGULAR MOMENTA CAN BE
C   ACCOMODATED.
C
C     THE STANDARD METHOD (COMPUTING A SUM OF PRODUCTS OF THREE 6-J'S)
C   IS USED.  ALL OF THE FACTORIALS COMMON TO ALL TERMS IN THE SUM ARE
C   FACTORED OUT AND EVALUATED ONLY ONCE.  SIMILARLY ALL THE INTEGER TO
C   REAL CONVERSIONS AND OTHER SETUP REQUIRED FOR TH 6-J CALCULATIONS    1090   
C   ARE CARRIED OUT ONLY ONCE.  TO EVALUATE THE 6-J'S A FACTORIZATION O
C   THE RACAH SUM SUGGESTED BY
C       J. G. WILLS, COMPUTER PHYSICS COMM. 2, 38 (1971)
C   IS USED.  THE LOGARITHMS OF THE REQUIRED FACTORIALS ARE FOUND
C   EITHER FROM A TABLE OF  LOG(FACTORIAL)  OR BY USING THE  LOG GAMMA
C   FUNCTION.  IN THE PRESENT IMPLIMENTATION, THE TABLE CONTAINS THE
C   FIRST 1001 FACTORIALS.  THIS SIZE MAY BE REDUCED (TO SAVE CORE) OR
C   INCREASED (TO SAVE TIME FOR VERY LARGE J'S) BY RECOMPILING THE
C   /LOGFAC/ COMMON BLOCK.
C                                                                        1100   
C     THE MAGNITUDE OF THE J'S THAT CAN BE ACCOMODATED BY THIS ROUTINE
C   IS DETERMINED BY CANCELLATIONS IN THE RACAH AND 9-J SUMS.  AS CAN B
C   SEEN IN THE TABLE BELOW, FAIRLY SEVERE CANCELLATIONS OCCUR WHEN ALL
C   NINE J'S ARE OF THE ORDER OF 50.  COMPLETE LOSS OF SIGNIFICANCE AND
C   EXTENSIVE NUMERIC OVERFLOWS OCCUR WHEN ALL THE J'S ARE 75.  HOWEVER
C   IF FIVE OF THE J'S REMAIN SMALL WHILE THE OTHER FOUR GET LARGE, THE
C   THE LARGE J'S MAY GET VERY LARGE.
C
C     THE NUMBER OF TERMS IN THE 9-J SUM IS
C            N = MIN ( J1+J9 , J2+J6 , J4+J8 ) -                         1110   
C                  MAX ( |J1-J9| , |J2-J6| , |J4-J8| )  +  1
C   NO ATTEMPT IS MADE TO MINIMIZE N BY INTERCHANGING ROWS OR COLUMNS.
C   IT MAY BE ADVANTAGEOUS TO THE USER TO ARRANGE THE J'S SUCH THAT N
C   IS MINIMIZED (AT LEAST APPROXIMATELY).  IN PARTICULAR IT SHOULD BE
C   NOTICED THAT  J3, J5, AND J7  ARE NOT USED IN DETERMINING  N  SO
C   THAT A ZERO VALUE OF ONE OF THESE IS "WASTED".
C
C     EACH LINE IN THE FOLLOWING TABLE IS THE RESULT OF SEVERAL THOUSAN
C   COMPUTATIONS WITH J'S IN THE INDICATED RANGES.  ERRORS WERE
C   DETERMINED BY COMPARASON WITH QUADRUPOLE PRECISION 9-J VALUES.       1120   
C   EXECUTION TIMES FOR THE 370/195 ARE GIVEN IN MICROSECONDS.  NONE OF
C   THE 9-J'S EVALUATED WERE TRIVIALLY ZERO.
C
C    J1 J2 J5     J3 J6 J9    J4   J7   J8       TIME     RMS RELATIVE
C     RANGE        RANGE      FIXED VALUES                    ERROR
C
C   1/2 - 5/2    1/2 - 5/2    3/2  3/2  3/2       300         1E-14
C     5 - 10       5 - 10    15/2   8  17/2      1450         1E-12
C    15 - 30       0 - 47     25  45/2   20      3000         3E-12
C                                                                        1130   
C    50 - 54      50 - 54     51   52  107/2    34000         4E-8
C    25 - 30     1/2 -  5     55/2  5   5/2       450         2E-14
C    25 - 30       5 - 10     55/2  7  15/2      1600         2E-13
C      250       1/2 -  5    505/2  2   5/2       440         2E-13
C      250         5 - 10    505/2  7  15/2      1680         5E-13
C      500       1/2 -  5   1005/2  2   5/2      4850         2E-12
C      500         5 - 10   1005/2  7  15/2     13800         2E-12
C     2500       1/2 -  5   5005/2  2   5/2      4700         2E-12
C     2500         5 - 10   5005/2  7  15/2     13700         2E-12
C                                                                        1140   
C
C   USAGE:
C
C     WIG9J (2*J1, 2*J2, 2*J3, 2*J4, 2*J5, 2*J6, 2*J7, 2*J8, 2*J9)
C
C   RETURNS THE DOUBLE PRECISION 9-J COEFICENT.  THE TRIANGLE RULES ARE
C   CHECKED AND 0 IS RETURNED IF THEY ARE VIOLATED.  HOWEVER NO CHECK I
C   MADE FOR INVALID COMBINATIONS OF HALF-INTEGER AND INTEGERS J'S.  TH
C   USER MUST INSURE THAT TWICE THE SUM OF THE J'S IN EACH ROW AND
C   COLUMN IS EVEN.  IF THIS IS NOT SO, MEANINGLESS RESULTS WILL BE      1150   
C   RETURNED.
C
C   CORE REQUIREMENTS (BASE 10):
C
C       WIG9J    4878
C      /LOGFAC/  8028    (MAY BE REDUCED OR INCREASED)
C       DLGAMA   1008
C       TOTAL   13914
C
C     IN ADDITION "DEXP" AND THE FORTRAN I/O PACKAGE ARE REQUIRED.       1160   
C   THESE ARE NOT INCLUDED IN THE ABOVE SINCE YOU ARE LIKELY TO USE THE
C   ANYWAY.  /LOGFAC/ WILL BE SHARED WITH THE CLEBSCH AND
C   6-J ROUTINES SHOULD YOU USE THEM.
C
C
C   S. PIEPER
C
C   9/6/74 - FIRST VERSION
C   4/11/88 - VAX/VMS VERSION - R.OSBORN AND G.L.GOODMAN
C   7/3/91 - RS6000 - s.p. (no problems)                                 1170   
C
C
      IMPLICIT INTEGER*4 (A-C, I-R, T-Z),  REAL*8 (D-H, S)
C
      REAL*8  WIG9J, DLGAMA
C
      REAL*8 LF
      COMMON /LOGFAC/  MAXLF, lfbias, LF(1)
C
C                                                                        1180   
C     THE FOLLOWING FAKE EXTERNALS DRAG IN THE REQUIRED BLOCK DATAS
C     FOR THE ABOVE COMMON BLOCK
C
      EXTERNAL   LOGDUM
C
CCCCC**********
CCCC      COMMON /NUMBLK/  NUMIT, NUM6J, BIGOUT, BIGIN, SUMLOG, SUMOUT
CCCCC**********
C
C     WE USE THE FOLLOWING VAX FORTRAN FUNCTIONS:                        1190   
C     BTEST - TESTS TO SEE IF A BIT IS ON.  BTEST(I, 0) IS TRUE
C            IF  I  IS ODD.
C     ISHFT(N,-1) SHIFTS THE BITS IN N ONE TO THE RIGHT  = /2
C
      LOGICAL BTEST
C
      REAL*8  ZERONE(2) / 0.D0, 1.D0 /
C
      DIMENSION  J1S(5), J2S(5), J3S(5)
      DIMENSION  AS(4, 3),  BS(2, 3),  DAS(4, 3),  DBS(2, 3),            1200   
     1  MX2S(3),  MN2S(3),  MN1S(3),  RS(3), DRS(3),
     2  DMN1S(3),  DMN2S(3), DMX2S(3),  IDS(2, 3)
      DIMENSION SUMS(3)
C
C     FOLLOWING IS FACTORIAL(N-1)
C
      DLF(N) = DLGAMA(DFLOAT(N))
C
CCCCC****************
CCCC      NUMIT = 0                                                      1210   
CCCC      NUM6J = 0
CCCC      BIGOUT = 0
CCCC      BIGIN = 0
CCCCC****************
C
C
C
C     PICK UP THE J'S IN THE APPROPRIATE ORDER
C     NOTE THAT THE JIS ARRAYS REALLY CONTAIN 2*J
C                                                                        1220   
      J1S(1) = J1
      J1S(2) = J6
      J1S(3) = J8
      J2S(1) = J2
      J2S(2) = J4
      J2S(3) = J9
      J3S(1) = J3
      J3S(2) = J5
      J3S(3) = J7
C                                                                        1230   
C     EXPAND THE ARRAYS TO ALLOW CYCLIC INDICES
C     AFTER THIS IS DONE THE 3 6-J'S ARE
C       (  J1S(I)    J2S(I)    J3S(I)   )
C       ( J1S(I+1)  J2S(I+2)     X      )
C     FOR I = 1, 2, 3.
C
      DO 59  I = 1, 2
         J1S(I+3) = J1S(I)
         J2S(I+3) = J2S(I)
         J3S(I+3) = J3S(I)                                               1240   
 59   CONTINUE
C
C     DETERMIN THE RANGE OF 2*X
C
      XSTART = IABS(J2S(1) - J1S(2))
      XEND = J2S(1) + J1S(2)
      XSUM = XEND
      DO  129  I = 2, 3
         XSTART = MAX0 ( XSTART, IABS(J2S(I)-J1S(I+1)) )
         XPIECE = J2S(I) + J1S(I+1)                                      1250   
         XEND = MIN0 ( XEND, XPIECE)
         XSUM = XSUM + XPIECE
 129  CONTINUE
C
      IF (XEND .LT. XSTART)  GO TO 900
C
      XSUM = XSUM - XEND
C
C     EXTRACT 1/2 INTEGER PART OF X AND CONVERT 2*X TO INTPART(X)
C                                                                        1260   
      XHALF = 0
      IF ( BTEST(XSTART, 0) ) XHALF = 1
      XSTART = ISHFT(XSTART, -1)
      XEND = ISHFT(XEND, -1)
      XSUM = ISHFT(XSUM, -1)
      DXHALF = ZERONE(XHALF+1)
      DX = XSTART - 1
      XSUM = XSUM + 1
C
CCCCC******************                                                  1270   
CCCC      NUMIT = XEND - XSTART + 1
CCCCC******************
C
C     SETUP THE J1+J2-J3, ETC
C
      DO 159  I = 1, 3
C
C     THE AS, FOR EACH 6-J THEY ARE
C       A1 = J1+J2-J3
C       A2 = J4+J5-J3                                                    1280   
C       A3 = J1+J5-J6
C       A4 = J4+J2-J6
C     IN A3 AND A4 THE  J6 = X  IS LEFT OUT.
C
         AS(1,I) = ISHFT (J1S(I)+J2S(I)-J3S(I), -1)
         AS(2,I) = ISHFT (J1S(I+1)+J2S(I+2)-J3S(I), -1)
         AS(3,I) = ISHFT (J1S(I)+J2S(I+2), -1)
         AS(4,I) = ISHFT (J2S(I)+J1S(I+1), -1)
C
C     NOW MAKE SURE NONE OF THEM ARE NEGATIVE WHICH WOULD INDICATE       1290   
C     A VIOLATION OF THE TRIANGLE RULES.
C
         IF (BTEST(AS(1,I),30) .OR. BTEST(AS(2,I),30)) GO TO 900
C
C     MAKE SURE NO SUPER LARGE INTEGERS SNEAK THROUGH
C        AS(3,I) = IBCLR(AS(3,I),30)
C        AS(4,I) = IBCLR(AS(4,I),30)
C
C     THE BS, FOR EACH 6-J THEY ARE
C        B1 = J3 + J6 - J1 - J4                                          1300   
C        B2 = J3 + J6 - J2 - J5
C     THE  J6 = X  IS LEFT OUT
C
         BS(1,I) = J2S(I+2) - AS(2,I) - AS(3,I)
         BS(2,I) = J1S(I) - AS(1,I) - AS(3,I)
C
 159  CONTINUE
C
C     SETUP THE QUANTITIES FOR EACH OF THE 3 6-J SUMS.
C     THESE INCLUDE THE STARTING MIN AND MAX VALUES, ETC.                1310   
C
      DO 299  I = 1, 3
         MN1S(I) = MIN0(AS(1,I), AS(2,I))
         MN2S(I) = MIN0(AS(3,I), AS(4,I))
         MX2S(I) = -MIN0(BS(1,I), BS(2,I))
         RS(I) = 1 + XHALF + AS(3,I) + AS(4,I)
         DO 239  J = 1, 4
            DAS(J,I) = AS(J,I)
 239     CONTINUE
         DBS(1,I) = BS(1,I)                                              1320   
         DBS(2,I) = BS(2,I)
         DRS(I) = 1 + DXHALF + DAS(3,I) + DAS(4,I)
         DMN1S(I) = MN1S(I)
         DMN2S(I) = MN2S(I)
         DMX2S(I) = MX2S(I)
C
C     THESE ARE INDICES FOR THE 4 FACTORIALTS THAT DEPEND ON X
C
         IDS(1,I) = AS(2,I) + BS(2,I)
         IDS(2,I) = AS(1,I) + BS(1,I)                                    1330   
 299  CONTINUE
C
C     FIND THE LOG OF THE OUTSIDE FACTORS NOW
C     IT WILL BE USED IN EACH TERM OF THE 9-J SUM TO AVOID UNDER/OVER
C     FLOW PROBLEMS.
C
C
C     THESE ARE ALL THE TRIANGLE FACTORIALS THAT COMPLTELY FACTOR OUT
C     OF THE SUM
C                                                                        1340   
      SUMOUT = 0
C
C     XSUM IS >= ANY OF THE REQUIRED FACTORIAL ARGUMENTS
C
      IF (XSUM .GT. MAXLF)  GO TO 350
C
C     WE CAN USE THE LOG(FACTORIAL) TABLE
C
      DO 339  I = 1, 3
         DO 339  J = 1, 2                                                1350   
            SUMOUT = SUMOUT +
     1        lf(lfbias+1+AS(J,I)) - lf(lfbias+2+J3S(I)+AS(J,I))
     2        + lf(lfbias+1-AS(J,I)+J1S(I-1+J))
     3        + lf(lfbias+1-AS(J,I)+J2S(I-2+2*J))
 339  CONTINUE
      GO TO 400
C
C     WE MUST USE LOG GAMMA
C
 350  DO 359  I = 1, 3                                                   1360   
         DO 359  J = 1, 2
            SUMOUT = SUMOUT +
     1        DLGAMA(1+DAS(J,I)) - DLF(2+J3S(I)+AS(J,I))
     2        + DLF(1-AS(J,I)+J1S(I-1+J))
     3        + DLF(1-AS(J,I)+J2S(I-2+2*J))
 359  CONTINUE
C
CCCCC************
CCCC      BIGOUT = 1
CCCCC***********                                                         1370   
C
 400  SUMOUT = .5*SUMOUT
C
C     WE ARE NOW READY TO DO THE 9-J LOOP
C
      WIG9J = 0
C
      DO 799  X = XSTART, XEND
         DX = DX + 1
C                                                                        1380   
C     DO THE PARTS OF EACH 6-J THAT DEPEND ON X
C
         SUMLOG = SUMOUT
         ISIGN = 0
C
         DO 759  I = 1, 3
C
C     GET THE COEFFICENTS FOR THIS 6-J
C
            MX2 = MX2S(I) - X                                            1390   
            MN2 = MN2S(I) - X
            NMAX = MN1S(I)
            DNMAX = DMN1S(I)
            DNMIN = 0
            NMIN = 0
C
C     GET THE RANGE OF THE 6-J SUM
C
            IF (MN2 .GE. NMAX)  GO TO 520
            NMAX = MN2                                                   1400   
            DNMAX = DMN2S(I) - DX
C
 520        IF (MX2 .LE. 0)  GO TO 540
            NMIN = MX2
            DNMIN = DMX2S(I) - DX
C
 540        SUM = 1
            ISIGN = ISIGN + NMIN
C
C                                                                        1410   
C     THIS ONE TEST CHECKS ALL 12 INEQUALITIES IMPLIED BY THE
C     4 TRIANGLE CONDITIONS.FOR THE GIVEN 6-J.
C
            IF (NMAX - NMIN)  900, 700, 600
C
 600        D1 = DAS(1,I) - (DNMAX-1)
            D2 = DAS(2,I) - (DNMAX-1)
            D8 = DRS(I) - (DNMAX-1)
            D3 = DAS(3,I) - (DNMAX-1+DX)
            D4 = DAS(4,I) - (DNMAX-1+DX)                                 1420   
            D5 = -DNMAX
            D6 = -DBS(1,I) - (DX + DNMAX)
            D7 = -DBS(2,I) - (DX + DNMAX)
C
            E1 = D1*D2*D3*D4
            E5 = D5*D6*D7*D8
C
            NMAXM2 = NMAX-2
            SUM = 1 + E1/E5
            IF (NMAXM2 .LT. NMIN)  GO TO 700                             1430   
C
            F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1
            F5 = (D5+1)*(D6+1)*(D7+1)*(D8+1) - E5
C
            G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1
            G5 = (D5+2)*(D6+2)*(D7+2)*(D8+2) - F5 - F5 - E5
C
            H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1
            H5 = (D5+3)*(D6+3)*(D7+3)*(D8+3) - 3*(G5+F5) - E5
C                                                                        1440   
C           COMPUTE THE SUM IN THE RACAH FORMULA
C
            DO 659  N = NMIN, NMAXM2
               E1 = E1 + F1
               E5 = E5 + F5
               F1 = F1 + G1
               F5 = F5 + G5
               G1 = G1 + H1
               G5 = G5 + H5
               H1 = H1 + 24                                              1450   
               H5 = H5 + 24
               SUM = 1 + SUM*(E1/E5)
 659        CONTINUE
C
C
700         SUMS(I) = SUM
C
CCCCC***********
CCCC            NUM6J = NUM6J + NMAX - NMIN + 1
CCCCC************                                                        1460   
C
C     HERE ARE THE PARTS OF THE FACTORIALS FOR NMIN AND ALSO THOSE
C     FACTORIALS THAT DEPEND ON X
C
      IF (RS(I) .GT. MAXLF)  GO TO 730
C
C     WE CAN USE THE LOG FACTORIAL TABLE
C
            SUMLOG = SUMLOG +
     1        lf(lfbias+1-X+AS(4,I)) + lf(lfbias+1+X+IDS(1,I))           1470   
     2         + lf(lfbias+1+X+IDS(2,I)) - lf(lfbias+2+XHALF+X+AS(4,I))
     3         + lf(lfbias+1-NMIN+RS(I)) - lf(lfbias+1+NMIN)
     4         - lf(lfbias+1+NMIN+X+BS(1,I))
     &         - lf(lfbias+1+NMIN+X+BS(2,I))
     5         - lf(lfbias+1-NMIN+AS(1,I)) - lf(lfbias+1-NMIN+AS(2,I))
     6         - lf(lfbias+1-NMIN-X+AS(3,I))
     &         - lf(lfbias+1-NMIN-X+AS(4,I))
C
      GO TO 759
C                                                                        1480   
C     WE MUST USE LOG GAMMA
C
 730        SUMLOG = SUMLOG +
     1        DLGAMA(1-DX+DAS(4,I)) + DLGAMA(1+DX+DAS(2,I)+DBS(2,I))
     1         + DLGAMA(1+DX+DAS(1,I)+DBS(1,I))
     2         - DLGAMA(2+DXHALF+DX+DAS(4,I))
     3         + DLGAMA(1-DNMIN+DRS(I)) - DLGAMA(1+DNMIN)
     4         - DLGAMA(1+DNMIN+DX+DBS(1,I))
     4          - DLGAMA(1+DNMIN+DX+DBS(2,I))
     5         - DLGAMA(1-DNMIN+DAS(1,I)) - DLGAMA(1-DNMIN+DAS(2,I))     1490   
     6         - DLGAMA(1-DNMIN-DX+DAS(3,I))
     7         - DLGAMA(1-DNMIN-DX+DAS(4,I))
C
CCCCC***************
CCCC      BIGIN = 1
CCCCC**************
C
 759     CONTINUE
C
C     NOW MULTIPLY THE 3 6-J PARTS TOGETHER AND SUM                      1500   
C
         SUMP = SUMS(1) * SUMS(2) * DEXP(SUMLOG) * SUMS(3)
     1      * (1+DXHALF+(DX+DX))
C
         IF (BTEST(ISIGN, 0))   SUMP = -SUMP
         WIG9J =  WIG9J + SUMP
C
 799  CONTINUE
C
C                                                                        1510   
C
      RETURN
C
C     THE TRIANGLE RULES ARE NOT SATISFIED
C
 900  WIG9J = 0
      RETURN
      END
      BLOCKDATA FACDUM
C                                                                        1520   
C     BLOCK DATA TO DEFINE THE SQRT(FACTORIAL) TABLE
C
C     THE TABLE IS USED BY THE CLEBSH, THREEJ, SIXJ AND RACAH
C     ROUTINES.
C
C     IT MAY BE INCLUDED IN YOUR LOAD MODULE FOR OTHER PURPOSES
C     BY PLACING A DUMMY EXTERNAL CARD OF THE FORM
C          EXTERNAL  FACDUM
C     IN YOUR SOURCE.
C                                                                        1530   
C     THE FORM OF THE TABLE IS
C
C     COMMON / FACTRL /  MAX, ISPACE, TABLE
C
C
C     MAX IS THE MAXIMUM NUMBER WHOSE SQRT(FACTORIAL) IS
C         IN TABLE.  TABLE HAS MAX+1 ELEMENTS.
C         IN THE PRESENT VERSION:
C            MAX = 96    (DETERMINED BY EXPONENT RANGE)
C     ISPACE IS NOT USED.                                                1540   
C     TABLE IS THE TABLE.  THE FIRST ELEMENT CORRESPONDS TO 0.
C
C     6/30/74 - FIRST VERSION - S. PIEPER
C     4/10/88 - VAX/VMS VERSION - R. OSBORN AND G.L.GOODMAN
C
C
      INTEGER*4 MAXFAC, SP1
      REAL*8  FACTBL
      COMMON  /FACTRL/  MAXFAC, SP1, FACTBL(97)
C                                                                        1550   
C
      DATA  MAXFAC /96/
C
C     THESE ARE THE SQUARE ROOTS OF FACTORIALS.
C
C     THEY ARE STORED UP TO THE EXPONENT RANGE OF THE MACHINE.
C
      DATA FACTBL( 1) /  0.1000000000000000D+01 /
      DATA FACTBL( 2) /  0.1000000000000000D+01 /
      DATA FACTBL( 3) /  0.1414213562373095D+01 /                        1560   
      DATA FACTBL( 4) /  0.2449489742783178D+01 /
      DATA FACTBL( 5) /  0.4898979485566356D+01 /
      DATA FACTBL( 6) /  0.1095445115010332D+02 /
      DATA FACTBL( 7) /  0.2683281572999748D+02 /
      DATA FACTBL( 8) /  0.7099295739719539D+02 /
      DATA FACTBL( 9) /  0.2007984063681781D+03 /
      DATA FACTBL(10) /  0.6023952191045344D+03 /
      DATA FACTBL(11) /  0.1904940943966505D+04 /
      DATA FACTBL(12) /  0.6317974358922328D+04 /
      DATA FACTBL(13) /  0.2188610518114176D+05 /                        1570   
      DATA FACTBL(14) /  0.7891147445080469D+05 /
      DATA FACTBL(15) /  0.2952597012800765D+06 /
      DATA FACTBL(16) /  0.1143535905863913D+07 /
      DATA FACTBL(17) /  0.4574143623455652D+07 /
      DATA FACTBL(18) /  0.1885967730625315D+08 /
      DATA FACTBL(19) /  0.8001483428544984D+08 /
      DATA FACTBL(20) /  0.3487765766344294D+09 /
      DATA FACTBL(21) /  0.1559776268628498D+10 /
      DATA FACTBL(22) /  0.7147792818185865D+10 /
      DATA FACTBL(23) /  0.3352612008237171D+11 /                        1580   
      DATA FACTBL(24) /  0.1607856235454059D+12 /
      DATA FACTBL(25) /  0.7876854713229383D+12 /
      DATA FACTBL(26) /  0.3938427356614691D+13 /
      DATA FACTBL(27) /  0.2008211794424596D+14 /
      DATA FACTBL(28) /  0.1043497458090740D+15 /
      DATA FACTBL(29) /  0.5521669535672285D+15 /
      DATA FACTBL(30) /  0.2973510046012911D+16 /
      DATA FACTBL(31) /  0.1628658527169496D+17 /
      DATA FACTBL(32) /  0.9067986906793549D+17 /
      DATA FACTBL(33) /  0.5129628026803635D+18 /                        1590   
      DATA FACTBL(34) /  0.2946746955341073D+19 /
      DATA FACTBL(35) /  0.1718233974287565D+20 /
      DATA FACTBL(36) /  0.1016520927791757D+21 /
      DATA FACTBL(37) /  0.6099125566750542D+21 /
      DATA FACTBL(38) /  0.3709953246501409D+22 /
      DATA FACTBL(39) /  0.2286968774309350D+23 /
      DATA FACTBL(40) /  0.1428211541796153D+24 /
      DATA FACTBL(41) /  0.9032802905233224D+24 /
      DATA FACTBL(42) /  0.5783815921445271D+25 /
      DATA FACTBL(43) /  0.3748341123420972D+26 /                        1600   
      DATA FACTBL(44) /  0.2457951648494613D+27 /
      DATA FACTBL(45) /  0.1630420674178431D+28 /
      DATA FACTBL(46) /  0.1093719437815202D+29 /
      DATA FACTBL(47) /  0.7417966136220958D+29 /
      DATA FACTBL(48) /  0.5085501366740237D+30 /
      DATA FACTBL(49) /  0.3523338699662023D+31 /
      DATA FACTBL(50) /  0.2466337089763416D+32 /
      DATA FACTBL(51) /  0.1743963680863606D+33 /
      DATA FACTBL(52) /  0.1245439180886559D+34 /
      DATA FACTBL(53) /  0.8980989654316716D+34 /                        1610   
      DATA FACTBL(54) /  0.6538259159791714D+35 /
      DATA FACTBL(55) /  0.4804619624270389D+36 /
      DATA FACTBL(56) /  0.3563201278858420D+37 /
      DATA FACTBL(57) /  0.2666455677120592D+38 /
      DATA FACTBL(58) /  0.2013129889124823D+39 /
      DATA FACTBL(59) /  0.1533154046820762D+40 /
      DATA FACTBL(60) /  0.1177637968756484D+41 /
      DATA FACTBL(61) /  0.9121944481710788D+41 /
      DATA FACTBL(62) /  0.7124466393192018D+42 /
      DATA FACTBL(63) /  0.5609810447812647D+43 /                        1620   
      DATA FACTBL(64) /  0.4452649004137245D+44 /
      DATA FACTBL(65) /  0.3562119203309796D+45 /
      DATA FACTBL(66) /  0.2871872314724746D+46 /
      DATA FACTBL(67) /  0.2333120097803461D+47 /
      DATA FACTBL(68) /  0.1909741105966688D+48 /
      DATA FACTBL(69) /  0.1574812859496909D+49 /
      DATA FACTBL(70) /  0.1308137807832727D+50 /
      DATA FACTBL(71) /  0.1094466613011557D+51 /
      DATA FACTBL(72) /  0.9222139602976428D+51 /
      DATA FACTBL(73) /  0.7825244940376377D+52 /                        1630   
      DATA FACTBL(74) /  0.6685892207860282D+53 /
      DATA FACTBL(75) /  0.5751421947239992D+54 /
      DATA FACTBL(76) /  0.4980877514193197D+55 /
      DATA FACTBL(77) /  0.4342228346904444D+56 /
      DATA FACTBL(78) /  0.3810289910601106D+57 /
      DATA FACTBL(79) /  0.3365156932181068D+58 /
      DATA FACTBL(80) /  0.2991016905800262D+59 /
      DATA FACTBL(81) /  0.2675246849288189D+60 /
      DATA FACTBL(82) /  0.2407722164359370D+61 /
      DATA FACTBL(83) /  0.2180285150390389D+62 /                        1640   
      DATA FACTBL(84) /  0.1986334304622628D+63 /
      DATA FACTBL(85) /  0.1820505461284133D+64 /
      DATA FACTBL(86) /  0.1678423103505356D+65 /
      DATA FACTBL(87) /  0.1556505553593457D+66 /
      DATA FACTBL(88) /  0.1451811729660402D+67 /
      DATA FACTBL(89) /  0.1361920123419132D+68 /
      DATA FACTBL(90) /  0.1284832874770429D+69 /
      DATA FACTBL(91) /  0.1218899489080934D+70 /
      DATA FACTBL(92) /  0.1162756005221389D+71 /
      DATA FACTBL(93) /  0.1115276380752381D+72 /                        1650   
      DATA FACTBL(94) /  0.1075533591796017D+73 /
      DATA FACTBL(95) /  0.1042768505784838D+74 /
      DATA FACTBL(96) /  0.1016365017512855D+75 /
      DATA FACTBL(97) /  0.9958302741285533D+75 /
      END
C
C
      BLOCKDATA LOGDUM
C
C     BLOCK DATA TO DEFINE THE LOG(FACTORIAL) TABLE.                     1660   
C
C     THIS TABLE IS USED BY THE CLEBSH, THREEJ, SIXJ AND RACAH
C     ROUTINES.
C
C     IT MAY BE INCLUDED IN YOUR LOAD MODULE FOR OTHER PURPOSES
C     BY PLACING A DUMMY EXTERNAL CARD OF THE FORM
C          EXTERNAL  LOGDUM
C     IN YOUR SOURCE.
C
C     THE FORM OF THE TABLE IS                                           1670   
C
C     COMMON / LOGFAC /  MAX, ISPACE, TABLE
C
C
C     MAX IS THE MAXIMUM NUMBER WHOSE LOG(FACTORIAL) IS
C         IN TABLE.  TABLE HAS MAX+1 ELEMENTS.
C         IN THE PRESENT VERSION:
C            MAX = 1000
C     ISPACE IS NOT USED.
C     TABLE IS THE TABLE.  THE FIRST ELEMENT CORRESPONDS TO 0.           1680   
C
C     6/30/74 - FIRST VERSION - S. PIEPER
C     4/10/88 - VAX/VMS VERSION - R. OSBORN AND G.L.GOODMAN
C
C
      INTEGER*4 MAXLF
      REAL*8  LF
      COMMON /LOGFAC/  MAXLF, lfbias, LF(1001)
C
      DATA    MAXLF /0/,  lfbias /0/                                     1690   
C
      DATA LF(   1) /  0.0                    /
      END
 
cvaxCS    REAL FUNCTION ALGAMA(X)
cvax      DOUBLE PRECISION FUNCTION DLGAMA(X)
cvaxC-------------------------------------------------------------------
cvaxC
cvaxC THIS ROUTINE CALCULATES THE LOG(GAMMA) FUNCTION FOR A REAL ARGUMEN
cvaxC     X.  COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1700   
cvaxC     1 AND 2.  THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE
cvaxC     LOG(GAMMA) TO AT LEAST 18 SIGNIFICANT DECIMAL DIGITS.  THE
cvaxC     APPROXIMATION FOR X .GE. 12 IS FROM REFERENCE 3.  APPROXIMATIO
cvaxC     FOR X .LT. 12.0 ARE UNPUBLISHED.  LOWER ORDER APPROXIMATIONS C
cvaxC     BE SUBSTITUTED ON MACHINES WITH LESS PRECISE ARITHMETIC.
cvaxC
cvaxC
cvaxC EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
cvaxC
cvaxC XBIG   - THE LARGEST ARGUMENT FOR WHICH LN(GAMMA(X)) IS REPRESENTA 1710   
cvaxC          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
cvaxC                  LN(GAMMA(XBIG)) = XINF.
cvaxC XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER.
cvaxC EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
cvaxC          1.0+EPS .GT. 1.0
cvaxC FRTBIG - ROUGH ESTIMATE OF THE FOURTH ROOT OF XBIG
cvaxC
cvaxC     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
cvaxC
cvaxC          IBM/370   UNIVAC/110X      VAX 11/780                     1720   
cvaxC           (D.P.)     (D.P.)      (S.P.)     (D.P.)    (G.F.)
cvaxC
cvaxC XBIG    4.293D+73  1.280D+305  2.057E+36  2.057D+36  1.28D305
cvaxC XINF    7.230D+75  8.980D+307  1.701E+38  1.701D+38  8.98D307
cvaxC EPS     2.220D-16  1.735D-018  5.960E-08  1.388D-17  1.70D-15
cvaxC FRTBIG  2.500D+18  1.800D+076  1.100E+09  1.100D+09  1.80D+76
cvaxC
cvaxC
cvaxC ERROR RETURNS
cvaxC                                                                    1730   
cvaxC  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
cvaxC     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
cvaxC     TO BE FREE OF UNDERFLOW AND OVERFLOW.
cvaxC
cvaxC
cvaxC
cvaxC OTHER SUBPROGRAMS REQUIRED (SINGLE PRECISION VERSION)
cvaxC
cvaxC      ALOG,EXP,FLOAT,IFIX,SIN
cvaxC                                                                    1740   
cvaxC OTHER SUBPROGRAMS REQUIRED (DOUBLE PRECISION VERSION)
cvaxC
cvaxC     DBLE,DEXP,DLOG,DSIN,FLOAT,IFIX,SNGL
cvaxC
cvaxC
cvaxC REFERENCES:
cvaxC
cvaxC  1) W. J. CODY AND K. E. HILLSTROM, 'CHEBYSHEV APPROXIMATIONS FOR
cvaxC     THE NATURAL LOGARITHM OF THE GAMMA FUNCTION,' MATH. COMP. 21,
cvaxC     1967, PP. 198-203.                                             1750   
cvaxC
cvaxC  2) K. E. HILLSTROM, ANL/AMD PROGRAM ANLC366S, DGAMMA/DLGAMA, MAY,
cvaxC     1969.
cvaxC
cvaxC  3) HART, ET. AL., COMPUTER APPROXIMATIONS, WILEY AND SONS, NEW
cvaxC     YORK, 1968.
cvaxC
cvaxC
cvaxC  AUTHOR: W. J. CODY
cvaxC          ARGONNE NATIONAL  LABORATORY                              1760   
cvaxC
cvaxC  LATEST MODIFICATION: JULY 14, 1983
cvaxC
cvaxC-------------------------------------------------------------------
cvax      INTEGER I
cvaxCS    REAL             C,CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT6
cvaxCS   1     P1,P2,P4,Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDE
cvaxCS   2     XINF,XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
cvax      DOUBLE PRECISION C,CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT6
cvax     1     P1,P2,P4,Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDE 1770   
cvax     2     XINF,XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
cvax      DIMENSION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)
cvaxC-------------------------------------------------------------------
cvaxC  MATHEMATICAL CONSTANTS
cvaxC-------------------------------------------------------------------
cvaxCS    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/
cvaxCS    DATA FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/
cvaxCS    DATA SQRTPI/0.9189385332046727417803297E0/
cvax      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/
cvax      DATA FOUR,THRHAL,TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/      1780   
cvax      DATA SQRTPI/0.9189385332046727417803297D0/
cvaxC-------------------------------------------------------------------
cvaxC  MACHINE DEPENDENT PARAMETERS
cvaxC-------------------------------------------------------------------
cvaxCS    DATA XBIG,XINF,EPS,FRTBIG/2.057E36,1.701E38,5.96E-8,1.1E9/
cvaxCD    DATA XBIG,XINF,EPS,FRTBIG/2.057D36,1.701D38,1.388D-17,1.1D9/
cvaxC VAX G_FLOAT PARAMETERS
cvax      DATA XBIG,XINF,EPS,FRTBIG/1.28D305,8.98D307,1.70D-16,1.80D+76/
cvaxC-------------------------------------------------------------------
cvaxC  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX       1790   
cvaxC     APPROXIMATION OVER (0.5,1.5).
cvaxC-------------------------------------------------------------------
cvaxCS    DATA D1/-5.772156649015328605195174E-1/
cvaxCS    DATA P1/4.945235359296727046734888E0,2.01811262085677508391556
cvaxCS   1        2.290838373831346393026739E3,1.13196720590338082868504
cvaxCS   2        2.855724635671635335736389E4,3.84849622844379335999026
cvaxCS   3        2.637748787624195437963534E4,7.22581397970028819769896
cvaxCS    DATA Q1/6.748212550303777196073036E1,1.11333239385719932351300
cvaxCS   1        7.738757056935398733233834E3,2.76398707440334070889858
cvaxCS   2        5.499310206226157329794414E4,6.16112218006600212783335 1800   
cvaxCS   3        3.635127591501940507276287E4,8.78553630243101317087083
cvax      DATA D1/-5.772156649015328605195174D-1/
cvax      DATA P1/4.945235359296727046734888D0,2.01811262085677508391556
cvax     1        2.290838373831346393026739D3,1.13196720590338082868504
cvax     2        2.855724635671635335736389D4,3.84849622844379335999026
cvax     3        2.637748787624195437963534D4,7.22581397970028819769896
cvax      DATA Q1/6.748212550303777196073036D1,1.11333239385719932351300
cvax     1        7.738757056935398733233834D3,2.76398707440334070889858
cvax     2        5.499310206226157329794414D4,6.16112218006600212783335
cvax     3        3.635127591501940507276287D4,8.78553630243101317087083 1810   
cvaxC-------------------------------------------------------------------
cvaxC  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
cvaxC     APPROXIMATION OVER (1.5,4.0).
cvaxC-------------------------------------------------------------------
cvaxCS    DATA D2/4.227843350984671393993777E-1/
cvaxCS    DATA P2/4.974607845568932035012064E0,5.42413859989107049410198
cvaxCS   1        1.550693864978364947665077E4,1.84793290444563242541722
cvaxCS   2        1.088204769468828767498470E6,3.33815296798702973591722
cvaxCS   3        5.106661678927352456275255E6,3.07410905485053955625092
cvaxCS    DATA Q2/1.830328399370592604055942E2,7.76504932144500587132304 1820   
cvaxCS   1        1.331903827966074194402448E5,1.13670582132196960893875
cvaxCS   2        5.267964117437946917577538E6,1.34670145431110169229005
cvaxCS   3        1.782736530353274213975932E7,9.53309559184435361339574
cvax      DATA D2/4.227843350984671393993777D-1/
cvax      DATA P2/4.974607845568932035012064D0,5.42413859989107049410198
cvax     1        1.550693864978364947665077D4,1.84793290444563242541722
cvax     2        1.088204769468828767498470D6,3.33815296798702973591722
cvax     3        5.106661678927352456275255D6,3.07410905485053955625092
cvax      DATA Q2/1.830328399370592604055942D2,7.76504932144500587132304
cvax     1        1.331903827966074194402448D5,1.13670582132196960893875 1830   
cvax     2        5.267964117437946917577538D6,1.34670145431110169229005
cvax     3        1.782736530353274213975932D7,9.53309559184435361339574
cvaxC-------------------------------------------------------------------
cvaxC  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
cvaxC     APPROXIMATION OVER (4.0,12.0).
cvaxC-------------------------------------------------------------------
cvaxCS    DATA D4/1.791759469228055000094023E0/
cvaxCS    DATA P4/1.474502166059939948905062E4,2.42681336948670450283631
cvaxCS   1        1.214755574045093227939592E8,2.66343244963097694989807
cvaxCS   2      2.940378956634553899906876E10,1.702665737765398868392998 1840   
cvaxCS   3      4.926125793377430887588120E11,5.606251856223951465078242
cvaxCS    DATA Q4/2.690530175870899333379843E3,6.39388565430009239898423
cvaxCS   2        4.135599930241388052042842E7,1.12087210961614794137657
cvaxCS   3      1.488613728678813811542398E10,1.016803586272438228077304
cvaxCS   4      3.417476345507377132798597E11,4.463158187419713286462081
cvax      DATA D4/1.791759469228055000094023D0/
cvax      DATA P4/1.474502166059939948905062D4,2.42681336948670450283631
cvax     1        1.214755574045093227939592D8,2.66343244963097694989807
cvax     2      2.940378956634553899906876D10,1.702665737765398868392998
cvax     3      4.926125793377430887588120D11,5.606251856223951465078242 1850   
cvax      DATA Q4/2.690530175870899333379843D3,6.39388565430009239898423
cvax     2        4.135599930241388052042842D7,1.12087210961614794137657
cvax     3      1.488613728678813811542398D10,1.016803586272438228077304
cvax     4      3.417476345507377132798597D11,4.463158187419713286462081
cvaxC-------------------------------------------------------------------
cvaxC  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
cvaxC-------------------------------------------------------------------
cvaxCS    DATA C/-1.910444077728E-03,8.4171387781295E-04,
cvaxCS   1     -5.952379913043012E-04,7.93650793500350248E-04,
cvaxCS   2     -2.777777777777681622553E-03,8.333333333333333331554247E- 1860   
cvaxCS   3      5.7083835261E-03/
cvax      DATA C/-1.910444077728D-03,8.4171387781295D-04,
cvax     1     -5.952379913043012D-04,7.93650793500350248D-04,
cvax     2     -2.777777777777681622553D-03,8.333333333333333331554247D-
cvax     3      5.7083835261D-03/
cvaxC-------------------------------------------------------------------
cvax      Y = X
cvax      IF ((Y .LE. ZERO) .OR. (Y .GT. XBIG)) GO TO 700
cvax      IF (Y .GT. TWELVE) GO TO 400
cvax      IF (Y .GT. FOUR) GO TO 300                                     1870   
cvax      IF (Y .GT. THRHAL) GO TO 200
cvax      IF (Y .GE. PNT68) GO TO 100
cvaxCS    CORR = -ALOG(Y)
cvax      CORR = -DLOG(Y)
cvax      XM1 = Y
cvax      IF (Y .GT. EPS) GO TO 120
cvax      RES = CORR
cvax      GO TO 900
cvaxC-------------------------------------------------------------------
cvaxC  0.5 .LT. X .LE. 1.5                                               1880   
cvaxC-------------------------------------------------------------------
cvax  100 CORR = ZERO
cvax      XM1 = (Y - HALF) - HALF
cvax  120 XDEN = ONE
cvax      XNUM = ZERO
cvax      DO 140 I = 1, 8
cvax         XNUM = XNUM*XM1 + P1(I)
cvax         XDEN = XDEN*XM1 + Q1(I)
cvax  140 CONTINUE
cvax      RES = CORR + (XM1 * (D1 + XM1*(XNUM/XDEN)))                    1890   
cvax      GO TO 900
cvaxC-------------------------------------------------------------------
cvaxC  1.5 .LT. X .LE. 4.0
cvaxC-------------------------------------------------------------------
cvax  200 XM2 = Y - TWO
cvax      XDEN = ONE
cvax      XNUM = ZERO
cvax      DO 240 I = 1, 8
cvax         XNUM = XNUM*XM2 + P2(I)
cvax         XDEN = XDEN*XM2 + Q2(I)                                     1900   
cvax  240 CONTINUE
cvax      RES = XM2 * (D2 + XM2*(XNUM/XDEN))
cvax      GO TO 900
cvaxC-------------------------------------------------------------------
cvaxC  4.0 .LT. X .LE. 12.0
cvaxC-------------------------------------------------------------------
cvax  300 XM4 = Y - FOUR
cvax      XDEN = -ONE
cvax      XNUM = ZERO
cvax      DO 340 I = 1, 8                                                1910   
cvax         XNUM = XNUM*XM4 + P4(I)
cvax         XDEN = XDEN*XM4 + Q4(I)
cvax  340 CONTINUE
cvax      RES = D4 + XM4*(XNUM/XDEN)
cvax      GO TO 900
cvaxC-------------------------------------------------------------------
cvaxC  EVALUATE FOR ARGUMENT .GE. 12.0,
cvaxC-------------------------------------------------------------------
cvax  400 RES = ZERO
cvax      IF (Y .GT. FRTBIG) GO TO 460                                   1920   
cvax      RES = C(7)
cvax      YSQ = Y * Y
cvax      DO 450 I = 1, 6
cvax         RES = RES / YSQ + C(I)
cvax  450 CONTINUE
cvax  460 RES = RES/Y
cvaxCS    CORR = ALOG(Y)
cvax      CORR = DLOG(Y)
cvax      RES = RES + SQRTPI - HALF*CORR
cvax      RES = RES + Y*(CORR-ONE)                                       1930   
cvax      GO TO 900
cvaxC-------------------------------------------------------------------
cvaxC  RETURN FOR BAD ARGUMENTS
cvaxC-------------------------------------------------------------------
cvax  700 RES = XINF
cvaxC-------------------------------------------------------------------
cvaxC  FINAL ADJUSTMENTS AND RETURN
cvaxC-------------------------------------------------------------------
cvaxCS900 ALGAMA = RES
cvax  900 DLGAMA = RES                                                   1940   
cvax      RETURN
cvaxC ---------- LAST CARD OF DLGAMA ----------
cvax      END
c*id* dsgmal
      SUBROUTINE DSGMAL(ETA,NO,DSG)
C
C
C     COULOMB PHASES
C
C     ?/?/? - ORIGINAL VERSION                                           1950   
C     12/28/79 - CDC TO CNI, TITLE COMMENT - RPG
C     7/9/80 - USE CDM PREFIX - S.P.
C     7/3/91 - RS6000 version - s.p.
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      DIMENSION P(26),Q(23),DSG(1)
c
      DATA P( 1) /    0.310165781012994887 0D-08 /
CIB      DATA P( 1) /Z39D524F0F53B0885/                                  1960   
      DATA P( 2) /    0.301594747899910910 0D-05 /
CIB      DATA P( 2) /Z3C32996552696E14/
      DATA P( 3) /    0.205644287153878958 0D-03 /
CIB      DATA P( 3) /Z3DD7A237CE95149F/
      DATA P( 4) /    0.444558861472056296 0D-02 /
CIB      DATA P( 4) /Z3F1235899B631A0C/
      DATA P( 5) /    0.424137251918122321 0D-01 /
CIB      DATA P( 5) /Z3FADBA03A99B1830/
      DATA P( 6) /    0.206378868197296478 0D+00 /
CIB      DATA P( 6) /Z4034D53ED97E2D92/                                  1970   
      DATA P( 7) /    0.546892269520120738 0D+00 /
CIB      DATA P( 7) /Z408C0121BC062DCA/
      DATA P( 8) /    0.794948065779220198 0D+00 /
CIB      DATA P( 8) /Z40CB81B7688A4B0A/
      DATA P( 9) /    0.593131560870811980 0D+00 /
CIB      DATA P( 9) /Z4097D778502A62D8/
      DATA P(10) /    0.177060007536960065 0D+00 /
CIB      DATA P(10) /Z402D53CDFDCCFE69/
      DATA P(11) /    0.241087856973594357 0D-08 /
CIB      DATA P(11) /Z39A5AC9FD40B5BA8/                                  1980   
      DATA P(12) /    0.932656994995955490 0D-06 /
CIB      DATA P(12) /Z3BFA5BB35F20D66A/
      DATA P(13) /    0.758581379314270606 0D-04 /
CIB      DATA P(13) /Z3D4F8B038B78FE9E/
      DATA P(14) /    0.200706718860656989 0D-02 /
CIB      DATA P(14) /Z3E8388FFEF99C18E/
      DATA P(15) /    0.188924467702797066 0D-01 /
CIB      DATA P(15) /Z3F4D622A9050FC9E/
      DATA P(16) /    0.444315921800891838 0D-01 /
CIB      DATA P(16) /Z3FB5FDE6B529BF7E/                                  1990   
      DATA P(17) /   -0.159601618385062274 0D+00 /
CIB      DATA P(17) /ZC028DBA6D35A3E69/
      DATA P(18) /   -0.685141773275969618 0D+00 /
CIB      DATA P(18) /ZC0AF65738557FF1A/
      DATA P(19) /   -0.577201920703612828 0D+00 /
CIB      DATA P(19) /ZC093C3814C9C4355/
      DATA P(20) /    0.139541051788112899 0D+03 /
CIB      DATA P(20) /Z428B8A825EB76321/
      DATA P(21) /    0.109052269358000365 0D+04 /
CIB      DATA P(21) /Z4344285CF3F17F1E/                                  2000   
      DATA P(22) /    0.485358524121951234 0D+03 /
CIB      DATA P(22) /Z431E55BC83CA29B9/
      DATA P(23) /   -0.658758053903885070 0D+03 /
CIB      DATA P(23) /ZC3292C20FD215CAA/
      DATA P(24) /   -0.413023138385032667 0D+02 /
CIB      DATA P(24) /ZC2294D6470917FF0/
      DATA P(25) /    0.922933578234238228 0D+01 /
CIB      DATA P(25) /Z4193AB5BFF4F68E2/
      DATA P(26) /   -0.999999999999999944 0D+00 /
CIB      DATA P(26) /ZC0FFFFFFFFFFFFFC/                                  2010   
      DATA Q( 1) /    0.162632093394676580 0D-05 /
CIB      DATA Q( 1) /Z3C1B48FEC701B849/
      DATA Q( 2) /    0.180858599479503871 0D-03 /
CIB      DATA Q( 2) /Z3DBDA4DC51C6A4B7/
      DATA Q( 3) /    0.577143957513920856 0D-02 /
CIB      DATA Q( 3) /Z3F17A3CB039DF446/
      DATA Q( 4) /    0.787659072077549621 0D-01 /
CIB      DATA Q( 4) /Z40142A00A37F4085/
      DATA Q( 5) /    0.545157981064667799 0D+00 /
CIB      DATA Q( 5) /Z408B8F7933D37D9D/                                  2020   
      DATA Q( 6) /    0.208553610110275223 0D+01 /
CIB      DATA Q( 6) /Z41215E5B1A4DD198/
      DATA Q( 7) /    0.457313299223146097 0D+01 /
CIB      DATA Q( 7) /Z41492B8D801E48B1/
      DATA Q( 8) /    0.569717859194618481 0D+01 /
CIB      DATA Q( 8) /Z415B27A4BD3E1545/
      DATA Q( 9) /    0.373731134057316638 0D+01 /
CIB      DATA Q( 9) /Z413BCC06F9EBB4E9/
      DATA Q(10) /    0.841871279655082615 0D-09 /
CIB      DATA Q(10) /Z3939DA5B2E6CDBD7/                                  2030   
      DATA Q(11) /    0.525199068334112846 0D-06 /
CIB      DATA Q(11) /Z3B8CFB67B871E12F/
      DATA Q(12) /    0.638793895356789098 0D-04 /
CIB      DATA Q(12) /Z3D42FB7E391D89CE/
      DATA Q(13) /    0.262822889969226652 0D-02 /
CIB      DATA Q(13) /Z3EAC3E5D2BAB275D/
      DATA Q(14) /    0.448275937331078468 0D-01 /
CIB      DATA Q(14) /Z3FB79D2390AC547B/
      DATA Q(15) /    0.345389220533413407 0D+00 /
CIB      DATA Q(15) /Z40586B6D8E94FA2B/                                  2040   
      DATA Q(16) /    0.122326029540586134 0D+01 /
CIB      DATA Q(16) /Z4113927963343730/
      DATA Q(17) /    0.188103263084231598 0D+01 /
CIB      DATA Q(17) /Z411E18B5AC02D34E/
      DATA Q(18) /   -0.141115652526220753 0D+03 /
CIB      DATA Q(18) /ZC28D1D9B6769D163/
      DATA Q(19) /   -0.996536251962568770 0D+03 /
CIB      DATA Q(19) /ZC33E48947CF01A61/
      DATA Q(20) /   -0.535816773446634613 0D+03 /
CIB      DATA Q(20) /ZC3217D1181089897/                                  2050   
      DATA Q(21) /    0.654427162588470765 0D+03 /
CIB      DATA Q(21) /Z4328E6D5A87038E8/
      DATA Q(22) /    0.421589251537016345 0D+02 /
CIB      DATA Q(22) /Z422A28AF51A1A909/
      DATA Q(23) /   -0.931266911567572087 0D+01 /
CIB      DATA Q(23) /ZC19500B154A4BE4F/
      DATA EZ    /    0.180554707160510697 0D+01 /
CIB      DATA EZ    /Z411CE385537EEB8A/
      DATA EZ1   /    0.180554771423339844 0D+01 /
CIB      DATA EZ1   /Z411CE38600000000/                                  2060   
      DATA EZ2   /    0.642628291517623633 0D-06 /
CIB      DATA EZ2   /Z3BAC811476376D16/
C
      X=DABS(ETA)
      XSQ=X*X
      IF(X.GT.4.0) GO TO 200
      IF(X.GT.2.0) GO TO 100
C
C     ABS(ETA) IN (0,2) RANGE SO=ETA(ETA**2-EZSQ)R(ETA**2)
C                                                                        2070   
      I=7
      J=1
      K=1
      M=1
      RSQ=XSQ
      GO TO 300
C
C     ABS(ETA) IN (2,4) RANGE SO=ETA*R(ETA**2)
C
  100 I=6                                                                2080   
      J=11
      K=10
      M=2
      RSQ=XSQ
      GO TO 300
C
C     ABS(ETA) IN (4,-) RANGE SO=ATAN(ETA)/2+ETA(LN(1+ETA**2)/2
C                                +R(1/ETA**2))
C
  200 I=4                                                                2090   
      J=20
      K=18
      M=3
      RSQ=1.0D0/XSQ
  300 QSUM=Q(K)
      N=K+I
      DO 305 L=K,N
  305 QSUM=QSUM*RSQ+Q(L+1)
      PSUM=P(J)
      N=J+I+1                                                            2100   
      DO 310 L=J,N
  310 PSUM=PSUM*RSQ+P(L+1)
      R=PSUM/(QSUM*RSQ+1.0D0)
      GO TO (400,500,600),M
  400 DSG(1)=X*R*(X+EZ)*((X-EZ1)+EZ2)
      GO TO 610
  500 DSG(1)=X*R
      GO TO 610
  600 DSG(1)=DATAN(X)/2.0D0+X*(DLOG(1.0D0+XSQ)/2.0D0+R)
C                                                                        2110   
C     SL(ETA)=SO(ETA)+SUM(I=1TOL)DATAN(ETA/I)
C
  610 IF(NO.EQ.0) GO TO 650
      R=1.0D0
      DO 620 I=1,NO
      RSQ=X/R
      DSG(I+1)=DSG(I)+DATAN(RSQ)
  620 R=R+1.0D0
      IF  (ETA.GE.0.0D0) GOTO 660
      DO 640 I=1,NO                                                      2120   
  640 DSG(I+1)=-DSG(I+1)
  650 IF (ETA.LT.0.0D0) DSG(1)=-DSG(1)
 660  RETURN
      END
C*ID* get_date
       SUBROUTINE get_date ( date )
C
C         RETURNS THE DATE IN THE FORM
C          19 Jan 01
C                                                                        2130   
C   DATE - a character*9
C
C  11/16/01 - all new
C
      character*9 date
      character*8 date8
      character*10 time
      character*5 zone
      integer values(8)
      character*3 months(12)                                             2140   
      data months / 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     &   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
!
      call date_and_time ( date8, time, zone, values )
C
      date = " "
      date(1:2) = date8(7:8)
      date(4:6) = months(values(2))
      date(8:9) = date8(3:4)
C                                                                        2150   
      RETURN
          END
c*id* gaussl
      SUBROUTINE GAUSSL(N,X,W)
C     POINTS AND WEIGHTS FOR GAUSS-LEGENDRE INTEGRATION OVER (-1 1)
C     COMPUTED BY METHOD GIVEN ON PP.88 AND 89 OF
C     'METHODS OF NUMERICAL INTEGRATION' DAVIS AND RABINOWITZ
C     ACADEMIC PRESS (1975)
C
C     THIS SUBROUTINE IS AN ADAPTATION OF SUBROUTINE GRULE LISTED ON     2160   
C     P.369 (APP 2) OF DAVIS AND RABINOWITZ
C     MHM  JUNE 29, 1975
C     7/5/91 - RS6000 - s.p.
C
C     METHOD: - - - - - - - - - - -
C     BASED ON A VERY CLOSE ESTIMATE FOR ZEROS OF PN(X)
C     THIS ESTIMATE USES THE ASYMPTOTIC RELATION BETWEEN PN(X)
C     AND THE COSINE FUNCTION (WHICH IS REMARKABLY ACCURATE)
C     THE INITIAL ESTIMATE IS IMPROVED BY A FOURTH-ORDER
C     NEWTON-RAPHSON PROCEDURE                                           2170   
C     TO SIMPLIFY MATTERS THIS IS CARRIED OUT IN TWO STEPS
C     1. A THIRD-ORDER NEWTON-RAPHSON FROM X0 TO X0+H
C     2. A FIRST-ORDER NEWTON-RAPHSON FROM X0+H TO X0+H+H1
C     FOR ALL ZEROS AND N UP TO ABOUT 2000,THIS
C     SUFFICES TO GIVE ZEROS TO AT LEAST 13 SIGNIFICANT FIGURES
C     GAUSS POINTS AND WEIGHTS ARE ALSO GOOD TO 13 OR MORE FIGURES
C
C     MORE DETAILS (WITH FORMULAS) ARE GIVEN BELOW CLOSE TO THE
C     APPROPRIATE PARTS OF THE PROGRAM LISTING
C                                                                        2180   
C     ACCURACY AND TIMING - - - - - -
C
C     GAUSSL EXECUTION TIMES : IN MILLISECONDS
C                            : UNCERTAINTY ABOUT .1 MS
C     /195   T = .035  +  .0133 * N  +  .00103 * N**2
C     /75    T = 12.1*T(195)
C
C     NUMBER OF SIGNIFICANT DIGITS : (TO .5 DIGIT)
C     DIGITS(BASE 10)=16.1-.91*LOG10(N)
C                                                                        2190   
C     NUMBER OF BINARY DIGITS LOST : (TO 1.5 BINARY DIGITS)
C     LOST DIGITS(BASE 2) = 2.5 + 3.2 * LOG10(N)
C
C     SUMMARY OF TIMES AND ACCURACIES FOR VARIOUS VALUES OF N
C
C       N   DIGITS(10)  DIGITS(2)    TIME(MS)
C           ACCURACY      LOST      /195    /75
C     ------------------------------------------------
C       10    15.7        3.7        .28     3.5
C       20    15.2        5.5        .75     9.2                         2200   
C       50    14.4        8.1       3.35    40.8
C      100    14.0        9.5      11.7    144.0
C      200    13.8       10.1      44.5    540.
C      500    13.7       10.5     270.0
C     1000    13.5       11.1    1060.0   12700.0
C     2000    13.2       12.1    4110.0
C     ______________________________________________
C
C     ABSCISSAE COME OUT IN ASCENDING ORDER FROM NEAR -1 TO NEAR +1
C -  -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      2210   
      IMPLICIT REAL*8 ( A-H, O-Z )
      DIMENSION X(N),W(N)
C
C     TAKE CARE OF SPECIAL CASES (N.LE.1)
C
CCC        write (0,*) 'gaussl', n
      IF(N.GT.0) GO TO 200
      write (6, 100)N
  100 FORMAT(1H0,2X,'0000 NUMBER OF GAUSS POINTS = ',I5,' 0000'/3X,
     1'0000             ABSURD             0000')                        2220   
      GO TO 600
  200 IF(N.NE.1) GO TO 300
      X(1)=0
      W(1)=2
      GO TO 600
C
C     BEGIN COMPUTATION
C     INITIAL ESTIMATE X0(I,N) OF I'TH ZERO OF PN(X) IS
C     -  -  -
C     X0(I,N)=(1-.125*(N**-2)+.125*(N**-3))*COS(T)  + O(N**-4)           2230   
C     T=((4*I-1)/(4*N+2))*PI
C     -  -  -
C     THE ZEROS ARE ORDERED SO THAT I=1 IS NEAREST TO X=1
C     1 > X(1) > X(2) > X(3) > ............> X(M) GE 0
C
  300 M=(N+1)/2
      E1=N*(N+1)
      E2=3.1415 92653 58979 320D0/(4*N+2)
      E3=1.-(1.-1.0D0/N)/(8.0D0*N*N)
      DO 500 I=1,M                                                       2240   
      T=(4*I-1)*E2
      X0=E3*DCOS(T)
C
C     USE RECURSION FORMULAS AND THE 2ND ORDER DIFFERENTIAL EQUATION
C     SATISFIED BY PN TO COMPUTE PK AND ITS FIRST FOUR DERIVATIVES
C     AT X0
C
      PKM1=1
      PK=X0
      FK=1                                                               2250   
      DO 400 K=2,N
      FK=FK+1
      T1=X0*PK
      PKP1=T1-PKM1-(T1-PKM1)/FK+T1
      PKM1=PK
  400 PK=PKP1
      DEN=1.-X0*X0
      D1=N*(PKM1-X0*PK)
      DPN=D1/DEN
      D2PN=(2.*X0*DPN-E1*PK)/DEN                                         2260   
      D3PN=(4.*X0*D2PN+(2.-E1)*DPN)/DEN
      D4PN=(6.*X0*D3PN+(6.-E1)*D2PN)/DEN
C
C     USE THIRD-ORDER NEWTON-RAPHSON METHOD TO IMPROVE ON THE INITIAL
C     ESTIMATE X0 OF THE ZERO OF PN
C     I.E H IS DETERMINED SUCH THAT PN(X0+H)=0 UP TO AND INCLUDING
C     THE TERM IN D3(PN)/DX IN THE TAYLOR EXPANSION OF PN ABOUT X0
C
      U=PK/DPN
      V=D2PN/DPN                                                         2270   
CCC
CCC  there is some strange Intel compiler bug that results in the
Ccc  simple 1-line computation of H failing on the 2nd call to gaussl
CCC
CCC        write (0,*) i,m,n
CCC        write (0,*) 'pk', pk,dpn,d2pn,d3pn,d4pn
CCC         write(0,*) 'u', u, v
CCC         write (0,*) 1.+.5*U*(V+U*(V*V-U*D3PN/(3.*DPN)))
CCC       h=0
CCC       H=-U*(1.+.5*U*(V+U*(V*V-U*D3PN/(3.*DPN))))                     2280   
CCC        write (0,*) 'h computed'
        h1 = (1.+.5*U*(V+U*(V*V-U*D3PN/(3.*DPN))))
        h2 = -u*h1
CCC         write (0,*) 'h1 h2', h1, h2
        h = h2
CCC        write (0,*) i,m,n,h
C
C     EVALUE PN AT AND D(PN)/DX AT X0+H USING TAYLOR EXPANSION OF PN
C     AT XO UP TO TERMS OF ORDER H**3. FOR D(PN)/DX THIS REQUIRES
C     D4(PN)/DX  AT X0.                                                  2290   
C     THESE VALUES OF PN AND DPN AT X0+H ARE THEN USED TO
C     COMPUTE A FURTHER FIRST-ORDER NEWTON-RAPHSON CORRECTION
C     I.E. FIND H1 SUCH THAT PN(X0+H+H1)=0 UP TO THE FIRST-DERIVATIVE
C     TERM (TERM OF ORDER H1) IN THE TAYLOR EXPANSION OF PN ABOUT X0+H
C
      P=PK+H*(DPN+.5*H*(D2PN+H/3.*(D3PN+.25*H*D4PN)))
      DP=DPN+H*(D2PN+.5*H*(D3PN+H*D4PN/3.))
      H=H-P/DP
      X(N+1-I)=X0+H
      X(I)=-X(N+1-I)                                                     2300   
C
C     CALCULATE WEIGHT W(I) CORRESPONDING TO I'TH GAUSS POINT X0I
C     W(I)=2*(1-X0I**2)/(N*PN-1(X0I))**2
C     PN-1(X) IS THE LEGENDRE POLYNOMIAL PN OF ORDER N-1
C
      FX=D1-H*E1*(PK+.5*H*(DPN+H/3.*(D2PN+.25*H*(D3PN+.2*H*D4PN))))
      W(N+1-I)=2.*(1-X(I)*X(I))/(FX*FX)
  500 W(I)=W(N+1-I)
      IF(M+M.GT.N) X(M)=0
 600  RETURN                                                             2310   
      END
C*ID* INTRPC   C109.PTOLEMY.FORTLIB                         PER704
      SUBROUTINE INTRPC (NCUBIC, XCUBES, AS, BS, CS, DS, NPTS, XS, YS)
C
C     INTERPOLATES TO A GRID OF XS USING INPUT CUBICS.
C
C       THIS ROUTINE MAY BE USED TO FIND INTERPOLATED VALUES OF A
C     FUNCTION AFTER EITHER THE B88888 ROUTINE "SPLNCB" OR THE
C     AMD ROUTINE "SMOOTH" HAS FOUND THE COEFICENTS OF THE INTERPOLATING
C     CUBICS.  FOR EACH DESIRE POINT, THE APPROPRIATE CUBIC IS FOUND     2320   
C     AND IS USED TO DO THE INTERPOLATIONS.  IF NECESSARY THE END
C     CUBICS WILL BE USED FOR EXTRAPOLATION, BUT THIS IS NOT RECOMMENDED
C     NOTE THAT THE X'S ATWHICH THE OUTPUT IS DESIRED NEED NOT BE
C     ORDERED.
C
C
C     THE ARGUMENTS ARE:
C
C     NCUBIC - THE NUMBER OF POINTS USED IN DETERMINING THE CUBICS.
C     XCUBES - THE X VALUES USED IN FINDING THE CUBICS.  THIS ARRAY      2330   
C              HAS NCUBIC ELEMENTS.  IT MUST BE MONOTONIC (INCREASING
C              OR DCREASING).
C     AS, BS, CS, AND DS - THE COEFICENTS OF THE CUBICS.  EACH ARRAY
C              IS NCUBIC-1 ELEMENTS LONG.  THE COEFICENTS ARE
C              DEFINED AS:
C        YS(N) = F(XS(N)) = AS(I) + DEL*(BS(I) + DEL*(CS(I) +
C                 DEL*DS(I)))
C        WHERE
C          DEL = XS(N)-XCUBES(I);  XCUBES(I) < XS(N) < XCUBES(I+1)
C        THE AS, ..., DS MAY BE FOUND USING EITHER SMOOTH OR SPLNCB.     2340   
C     NPTS - THE NUMBER OF OUTPUT VALUES OF YS TO BE GENERATED.
C     XS  - THE ARRAY OF X-VALUES AT WHICH THE INTERPOLATED VALUES
C           ARE DESIRED.  THIS ARRAY HAS NPTS ELEMENTS.  IF ANY
C           XS(I) IS NOT IN THE RANGE DEFINED BY XCUBES,
C           THE EXTREME CUBICS WILL BE USED FOR EXTRAPOLATION.
C     YS  - THIS ARRAY WILL GET F(XS(I)).  IT WILL HAVE NPTS ELEMENTS.
C
C     S. PIEPER
C
C     PDP  2/21/72                                                       2350   
C     12/10/72  /360 VERSION
C     6/26/73  REMOVE PRINTING OF ERRORS
C     7/15/73  GENERAL XS
C     8/29/74 - FIX ERRORS FOR NON-INCREASING XS.
C     9/4/74 - CHANGE NAME TO INTRPC FROM INTRP1 - S. P.
C     12/28/79 - INF FOR VAX, UNIVAC - RPG
C
      implicit real*8 ( a-h, o-z )                                      implicit
C
      DIMENSION XCUBES(*), AS(*), BS(*), CS(*), DS(*), XS(*), YS(*)      2360   
      REAL*8 INF /1.E300/                                               ci3e
CIB      REAL*8 INF / 1.D77 /
CVA      REAL*8 INF / 1.D38 /
CUN      REAL*8 INF / 1.D300 /
C
      XSIGN = DSIGN (1.D0, XCUBES(2)-XCUBES(1))
C
C     THESE ARE FAKE STARTING CONDITIONS THAT WILL CAUSE DIRECT
C     BRANCH TO STMNT 200 INSIDE THE LOOP.
      XPREV = INF                                                        2370   
      XNEXT = INF
C
C
      DO 399  I = 1, NPTS
      X = XS(I)
      XC = XSIGN*X
C
C     IS X STILL INSIDE THIS SPLINE
C
 100  IF (XC .LT. XNEXT)  GO TO 200                                      2380   
C
C     GO TO THE NEXT SPLINE.
C
      XPREV = XNEXT
      XBASE = XNEXT*XSIGN
      N = N+1
      XNEXT = XCUBES(N+1)*XSIGN
C     WHEN WE ARE AT THE LAST SPLINE, ALLOW IT TO GO TO INFINITY
      IF (N .EQ. NCUBIC-1)  XNEXT = INF
C                                                                        2390   
 150  A = AS(N)
      B = BS(N)
      C = CS(N)
      D = DS(N)
C     SEE IF X IS INSIDE THIS NEW SPLINE
      GO TO 100
C
C     WE HAVE NOT PASSED THE END OF THIS SPLINE, ARE WE STILL
C     BEYOND ITS BEGINNING.
C                                                                        2400   
 200  IF (XC .GE. XPREV)  GO TO 300
C
C     THE XS ARRAY IS NOT MONOTONIC INCREASING -- GO BACK TO START
C
C     ALSO INITIALIZE FOR FIRST TIME
      N = 1
      XPREV = -INF
      XNEXT = XCUBES(2)*XSIGN
      XBASE = XCUBES(1)
      GO TO 150                                                          2410   
C
C     THE SET OF SPLINES TO USE HAS BEEN FOUND, USE IT
C
 300  DEL = X-XBASE
      YS(I) = A + DEL*(B + DEL*(C + DEL*D))
 399  CONTINUE
C
      RETURN
      END
c*id* laguer                                                             2420   
      SUBROUTINE LAGUER ( NN, X, W, ALF, EPS, CSX, CSW, TSX, TSW,
     1   A, B, C )
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C     POINTS AND WEIGHTS FOR GAUSS-LAGUERRE INTEGRATION
C
C     THE ARGUMENTS ARE:
C
C     NN - THE NUMBER OF POINTS DESIRED.                                 2430   
C     X - WILL BE SET TO THE POINTS WITH THE MINIMUM POINT FIRST.
C     W - WILL BE SET TO THE WEIGHTS.
C         X AND W MUST BE REAL*8 ARRAYS OF DIMENSION NN.
C     ALF - THE POWER OF X TO BE INCLUDED IN THE WEIGHT FUNCTION.
C           MUST HAVE -1 < ALF.  LARGE ALF CAN RESULT IN OVERFLOWS.
C           NOTE THAT ALF DOES NOT NEED TO BE AN INTEGER.
C     EPS - THE DESIRED ACCURACY OF THE ROOTS.  EPS = 1.D-20 IS
C           REASONABLE.
C     CSX, CSW, TSX, TSW - WILL BE SET TO AN INDICATION OF THE
C          ACCURACY OF X AND W AS FOLLOWS:                               2440   
C             CSX = CALC SUM X(I)        TSX = TRUE SUM X(I)
C             CSW = CALC SUM W(I)        TSW = TRUE SUM W(I)
C     A, B, C - REAL*8 SCRATCH ARRAYS OF DIMENSION NN.
C
C     THE POINTS AND WEIGHTS THAT ARE RETURNED ARE USED AS:
C
C     INTEGRAL( 0 TO INFINITY ) DX EXP(-X) X**ALF F(X)
C              =  SUM  W(I) F( X(I) )
C
C                   THIS PROGRAM IS TAKEN FROM:                          2450   
C                   STROUD AND SECREST, GAUSSIAN
C                   QUADRATURE FORMULAS,PRENTICE
C                   HALL (1966)  ..  PAGES 32-33
C
C
C     THIS ROUTINE CALLS  LGROOT, LGRECR , LAGBC  AND  DGAMMA
C
C     4/22/75 - F. SERDUKE
C     1/27/80 - USE STANDARD (A&S) LAGUERRE POLYNOMIALS TO PREVENT
C        OVERFLOW; PASS A, B, C AS ARGUMENTS - RPG                       2460   
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      DIMENSION X(NN) , W(NN) , A(NN) , B(NN) , C(NN)
C
C     THE COEFFICIENTS
C
C                     A(N)  =  1 / N                                     2470   
C
C                     B(N)  =  (ALF + 2N - 1)
C
C                     C(N)  =  (ALF + N - 1) / N
C
C      IN THE RECURRENCE RELATION
C
C             L(N)  =  A(N) * (B(N) - X) * L(N-1)  -  C(N)*L(N-2)
C
C     THESE ARE THE LAGUERRE POLYNOMIALS DEFINED IN                      2480   
C     ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS.
C
      CALL LAGBC ( NN , ALF , A , B , C )
      FN = NN
      CSX = 0.D0
      CSW = 0.D0
      CC = DGAMMA(ALF+1.D0)
      TSX = FN*(FN+ALF)
      TSW = CC
      DO 1 J = 2,NN                                                      2490   
    1 CC = CC*C(J)
      DO 7 I = 1,NN
      IF ( I - 2 ) 2, 4, 5
C         SMALLEST ZERO
    2 XT = (1.D0+ALF)*(3.D0+0.920D0*ALF) / (1.D0+2.40D0*FN + 1.80D0*ALF)
      GO TO 6
C        SECOND ZERO
    4 XT = XT + (15.D0 + 6.250D0*ALF) / (1.D0 + 0.90D0*ALF + 2.50D0*FN)
      GO TO 6
C        ALL OTHER ZEROS                                                 2500   
    5 FI = I-2
      R1 = (1.D0+2.550D0*FI) / (1.90D0*FI)
      R2 = 1.260D0 * FI*ALF / (1.D0+3.50D0*FI)
      RATIO = (R1+R2) / (1.D0 + 0.30D0*ALF)
      XT = XT + RATIO*(XT-X(I-2))
C
    6 CONTINUE
      CALL LGROOT(XT,NN,ALF,DPN,PN1,A,B,C,EPS)
      X(I) = XT
      W(I) = -CC / DPN / PN1                                             2510   
      CSX = CSX + XT
    7 CSW = CSW + W(I)
      RETURN
      END
      SUBROUTINE LGROOT ( X , NN , ALF , DPN , PN1 , A , B , C , EPS )
      IMPLICIT REAL*8 ( A-H, O-Z )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C        THIS SUBROUITNE IMPROVES THE APPROXIMATE ROOT S OF THE
C        NN-TH ORDER LAGUERRE POLYNOMIAL.  IN ADDITION, IT ALSO
C        RETURNS                                                         2520   
C                DPN  =  DERIVATIVE OF P(N) AT X
C                PN1  =  VALUE OF P(N-1) AT X
C        WHICH ARE NEEDED IN THE CALCULATION OF THE WEIGHTS FOR
C        NN-TH ORDER LAGUERRE QUADRATURE FORMULA.
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DIMENSION A(NN), B(NN) , C(NN)
C
      DO 2 ITER = 1, 10                                                  2530   
      CALL LGRECR ( P , DP , PN1 , X , NN , ALF , A , B , C )
      D = P/DP
      X = X - D
      IF ( DABS(D/X) .LT. EPS ) GO TO 3
    2 CONTINUE
    3 DPN = DP
      RETURN
      END
      SUBROUTINE LGRECR ( PN , DPN , PN1 , X , NN , ALF , A , B , C )
      IMPLICIT REAL*8 ( A-H, O-Z )                                       2540   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C        THIS SUBROUTINE IS THE RECURRENCE RELATION FOR THE LAGUERRE
C        POLYNOMIALS
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DIMENSION A(NN) , B(NN) , C(NN)
C
C
      P1 = 1.D0
      P = ALF + 1.D0 - X
      DP1 = 0.D0                                                         2550   
      DP = -1.D0
      DO 1 J = 2,NN
      Q = A(J) * (B(J) - X) * P - C(J) * P1
      DQ = A(J) * (B(J) - X) * DP - C(J) * DP1 - A(J)*P
      P1 = P
      P = Q
      DP1 = DP
    1 DP = DQ
      PN = P
      DPN = DP                                                           2560   
      PN1 = P1
      RETURN
      END
      SUBROUTINE LAGBC ( NN , ALF , A , B , C )
      IMPLICIT REAL*8 ( A-H, O-Z )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS A(N), B(N), AND C(N)
C     REQUIRED BY THE LAGUERRE QUADRATURE W AND X GENERATION
C     SUBROUTINE.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  2570   
      DIMENSION A(NN) , B(NN) , C(NN)
C
      DO 10 N = 1,NN
      EN = N
      A(N) = 1.D0 / EN
      B(N) = ALF + EN + EN - 1.D0
      C(N) = (ALF + EN - 1.D0) / EN
   10 CONTINUE
      RETURN
      END                                                                2580   
C*ID* LINSAV   C109.PTOLEMY.FORTLIB                         PER704  21:
      SUBROUTINE LINSAV (ASV, NMAX, NDIM, BSV, MMAX, A, B, POWER, DET)
C
C     UNIVAC VERSION OF 7/21/70
C     7/14/72  ARGONNE VERSION
C     11/06/72  REVERSE ORDER OF SUBSCRIPTS FOR /195
C     1/16/79 - CORRECT FORMAT ERROR - RPG
C
C     SOLUTION OF SEVERAL SETS OF SIMULTANEOUS LINEAR EQUATIONS.
C     THE MATRIX MUST BE SYMETRIC.  ITS MODE MAY BE CHANGED BY CHANGING  2590   
C      THE MARKED CARDS BELOW.  NOTE THAT EVEN COMPLEX MATRICES ARE
C      SYMETRIC... NOT HERMEATIAN.
C     SOLVE AX = B
C     THERE IS NO PIVOT SEARCHING.  THE DETERMINANT IS COMPUTED BUT
C     NO CHECK FOR SINGULARITY IS MADE.
C
C     'LINSAV'... THIS ENTRY PRESERVES THE ORIGINAL MATRIX AND
C      INHOMOGENEOUS TERM.
C     'ASV' IS THE MATRIX
C     'BSV' IS THE INHOMOGENEOUS TERM.                                   2600   
C     'A' WILL CONTAIN THE TRIANGULAR MATRIX.  THE UPPER HALF IS THE
C      TI  TRIANGLE, THE LOWER HALF IS JUNK.
C     'B' WILL HAVE THE SOLUTIN FOR X.
C     'NMAX' IS THE ORDER OF THE MATRIX--NOTE 'NMAX' = 1 IS ALLOWED.
C     'NDIM' IS THE FIRST DIMENSION OF 'ASV', 'BSV', 'A', AND 'B'.
C     'MMAX' IS THE NUMBER OF SETS OF EQUATIONS (NUMBER OF COLUMNS IN
C            'B' AND 'BSV').
C     'DET' WILL CONTAIN
C       (SIGN OF DET) * @DET@**POWER
C     THIS IS TO HELP PREVENT UNDER/OVER FLOW IN THE DETERMINANT.        2610   
C
C     NOTE THAT THE FORTLIB ROUTINE FIHZ IS USED.  AT RETURN FIHZ WILL
C      BE ON SO THAT 10 UNDER OR OVER FLOWS WILL CAUSE  TERMINATION OF
C
C     IMPLICIT COMPLEX*16 (A-H, O-Z)
      implicit real*8 ( a-h, o-z )                                      implicit
C
C     NOTE THAT THE SECOND DEMENSION OF B, BSV IS UNNECESSARY
      DIMENSION A(NDIM, NDIM), ASV(NDIM, NDIM), B(NDIM, MMAX),
     1  BSV(NDIM, MMAX)                                                  2620   
C
C     MOVE THE STUFF OVER TO SAVE IT
      IF ((NMAX .LE. 0) .OR. (MMAX .LE. 0) .OR. (NDIM .LT. NMAX))GOTO900
      DO 10 I = 1, NMAX
      DO 10 J = 1, NMAX
   10 A(J,I) = ASV(J,I)
C
      ENTRY LINBSV (A, NMAX, NDIM, BSV, MMAX, B, POWER, DET)
C
C     'LINBSV'... THIS ENTRY PRESERVES ONLY THE INHOMOGENEOUS TERM.      2630   
C     'A' AT ENTRY THIS IS THE MATRIX, AT EXIT IT IS AS ABOVE
C     'BSV' AND 'B' ARE AS ABOVE
C     THE DIMENSIONS AND OTHER ARGS ARE AS ABOVE
      IF ((NMAX .LE. 0) .OR. (MMAX .LE. 0) .OR. (NDIM .LT. NMAX))GOTO900
  3   DO 5 M = 1, MMAX
      DO 5 I = 1, NMAX
  5   B(I,M) = BSV(I,M)
C
      ENTRY LIN (A, NMAX, NDIM, B, MMAX, POWER, DET)
C                                                                        2640   
C     'LIN'... THIS ENRRY SAVES NOTHING
C     'A' IS AS UNDER 'LINBSV'
C     'B' AT ENTRY THIS IS THE INHOMOGENEOUS TERM, AT EXIT IT IS AS
C         UNDER 'LINSAV'.
C     THE DIMENSIONS AND OTHER ARGUEMENTS ARE ALL AS ABOVE.
C
      IF ((NMAX .LE. 0) .OR. (MMAX .LE. 0) .OR. (NDIM .LT. NMAX))GOTO900
      IF (NMAX .GT. 1) GO TO 100
      DO 50 M = 1, MMAX
   50 B(1,M) = B(1,M)/A(1,1)                                             2650   
      DET = A(1,1)
      RETURN
C
C
  100 NMAXM1 = NMAX-1
C
      DO 150 K = 1, NMAXM1
      KP1 = K+1
      DO 150 J = KP1, NMAX
      PIV = -A(J,K)/A(K,K)                                               2660   
      DO 130 M = 1, MMAX
  130 B(J,M) = B(J,M) + PIV*B(K,M)
      DO 150 I = J, NMAX
  150 A(I,J) = A(I,J) + PIV*A(I,K)
C
      DO 250 M = 1, MMAX
      B(NMAX,M) = B(NMAX,M)/A(NMAX,NMAX)
      DO 250 II = 1, NMAXM1
      I = NMAX-II
      SUM = 0                                                            2670   
      IP1 = I+1
      DO 230 J = IP1, NMAX
  230 SUM = SUM + A(J,I)*B(J,M)
  250 B(I,M) = (B(I,M) - SUM)/A(I,I)
C
C******* CALLS TO FIHZ TEMPORARILY REMOVED ************
C     CALL FIHZ (12, 'OFF')
C     CALL FIHZ (13, 'OFF')
      IF ((POWER .NE. 1.D0) .AND. (POWER .NE. 0)) GO TO 400
      DET = A(1,1)                                                       2680   
      DO 350 I = 2, NMAX
  350 DET = DET*A(I,I)
      GO TO 500
  400 DET = 1.
      DO 450 I = 1, NMAX
  450 DET = DET*DSIGN(DABS(A(I,I))**POWER, A(I,I))
  500 CONTINUE
C     CALL FIHZ (12, 'ON')
C     CALL FIHZ (13, 'ON')
      RETURN                                                             2690   
C
C     HERE A(I,J) IS IN A(J,I+1); I .GE. J.  BSV AND A(J,I+1)
C       ARE SAVED.
C    *** NOTE ***  THIS WAS CHANGED FROM SAVING THE LOWER
C    TRIANGLE TO SAVING THE UPPER TRIANGLE ON 11/06/72 TO MAKE
C    THINGS FASTER ON THE /195
C
      ENTRY LINTRI (A, NMAX, NDIM, BSV, MMAX, B, POWER, DET)
      IF ((NMAX .LE. 0) .OR. (MMAX .LE. 0) .OR. (NDIM .LT. NMAX)
     1  ) GO TO 900                                                      2700   
      DO 515 I = 1, NMAX
      DO 515 J = 1, I
 515  A(I,J) = A(J,I+1)
      GO TO 3
C
  900 write (6, 903) NMAX, NDIM, MMAX
  903 FORMAT ( / '0$$$ LIN/LINSAV -- INVALID INPUT...', 3I10 )
      CALL SYSERR
C
      END                                                                2710   
C*ID* LSQPOL   C109.PTOLEMY.FORTLIB                         PER704  21:
      SUBROUTINE LSQPOL (X,Y,W,RESID,NSUB,SUM,LSUB,A,B,MSUB,NMAX,MMAX)
C
C     LEAST SQUARE POLYNOMIAL FIT
C
C     12/28/79 - REMOVE EXTRANEOUS IPIVOT - RPG
C     2/20/80 - SCALE X TO PREVENT OVERFLOW - RPG
C
C     X = INPUT ARRAY OF X VALUES, DIMENSION NSUB.
C     Y = INPUT ARRAY OF Y VALUES, DIMENSION (NMAX,LSUB).                2720   
C     W = INPUT ARRAY OF WEIGHTS, DIMENSION NSUB.
C     RESID = OUTPUT ARRAY OF RESIDUALS (FIT - INPUT Y),
C        DIMENSION (NMAX,LSUB).
C     NSUB = NUMBER OF X VALUES.
C     SUM = OUTPUT ARRAY FOR WEIGHTED SUMS OF SQUARES,
C        DIMENSION LSUB.
C     LSUB = THE NUMBER OF RIGHT-HAND-SIDE VECTORS Y TO BE FITTED.
C     A = OUTPUT ARRAY FOR THE ERROR MATRIX, DIMENSION (MMAX,MSUB).
C     B = OUTPUT ARRAY FOR THE COEFFICIENTS OF X**(I-1) IN THE FINAL
C        FIT, DIMENSION (MMAX,LSUB).                                     2730   
C     MSUB = NUMBER OF TERMS IN THE POLYNOMIAL = DEGREE - 1.
C        MSUB MUST BE <= 25 UNLESS ARRAY XPOWER (BELOW) IS EXTENDED.
C     NMAX = FIRST DIMENSION OF Y AND RESID.  MUST BE >= NSUB.
C     MMAX = FIRST DIMENSION OF A AND B.  MUST BE >= MSUB.
C
C     FOR SINGLE PRECISION- ADD C'S TO THE FOLLOWING 2 CARDS.
      implicit real*8 ( a-h, o-z )                                      implicit
      REAL*8 X,Y,W,RESID,SUM,A,B
      DIMENSION X(NSUB),Y(NMAX,LSUB),W(NSUB),RESID(NMAX,LSUB),
     1          SUM(LSUB),A(MMAX,MSUB),B(MMAX,LSUB)                      2740   
C
C     THIS COMMON LIMITS MSUB TO 25.
C
      COMMON/F402/XPOWER(100)
      EQUIVALENCE (K1,DETERM),(K2,I2),(TERM,POLY)
      N=NSUB
      L=LSUB
      M=MSUB
      M1=M+1
      M3=M+M+M                                                           2750   
      M31=M3-1
      M41=M31+M
C
C     SCALE X INTO (-1,1) TO PREVENT OVERFLOW OR UNDERFLOW.
C
      XMAX = 0.
      DO 20   K = 1, N
         XMAX = DMAX1( XMAX, DABS(X(K)) )
   20 CONTINUE
      XSCALE = 1./XMAX                                                   2760   
      DO 30  K = 1, N
         X(K) = XSCALE*X(K)
   30 CONTINUE
C
C     FORMATION AND INVERSION OF SYSTEM OF NORMAL EQUATIONS
C
      DO 100 K2=M1,M41
      XPOWER(K2)=0.0
  100 CONTINUE
      DO 200 K1=1,N                                                      2770   
      TERM=W(K1)
      DO 200 K2=M1,M31
      XPOWER(K2)=TERM+XPOWER(K2)
      TERM=X(K1)*TERM
  200 CONTINUE
      DO 300 I=1,M
      DO 300 J=1,M
      K2=I+J+M-1
      A(I,J)=XPOWER(K2)
  300 CONTINUE                                                           2780   
      DO 500 J=1,L
      DO 400 K=1,N
      TERM=W(K)*Y(K,J)
      DO 400 K2=M3,M41
      XPOWER(K2)=TERM+XPOWER(K2)
      TERM=X(K)*TERM
  400 CONTINUE
      DO 500 I=1,M
      K2=I+M31
      B(I,J)=XPOWER(K2)                                                  2790   
      IF(J.NE.L) XPOWER(K2)=0.0
  500 CONTINUE
C
C     FOR DOUBLE PRECISION- USE DOUBLE PRECISION VERSION OF MATINV
C
      CALL MATINV(A,M,B,L,DETERM,MMAX)
C
C     EVALUATION OF RESIDUALS
C
      DO 700 J=1,L                                                       2800   
      SUM(J)=0.0
      DO 700 K=1,N
      POLY=0.0
      DO 600 I2=1,M
      I=M1-I2
      POLY=X(K)*POLY+B(I,J)
  600 CONTINUE
      RESID(K,J)=POLY-Y(K,J)
      SUM(J)=SUM(J)+W(K)*RESID(K,J)**2
  700 CONTINUE                                                           2810   
C
C     NOW UN-SCALE THE COEFFICIENTS AND X.  DO IT THE SLOW WAY TO
C     PREVENT UNNECESSARY OVERFLOW OR UNDERFLOW.
C
      DO 800  J = 1, L
         DO 800  I2 = 2, M
            DO 800  I = I2, M
               B(I,J) = XSCALE*B(I,J)
  800 CONTINUE
      DO 900  K = 1, N                                                   2820   
         X(K) = XMAX*X(K)
  900 CONTINUE
      RETURN
      END
C*ID* MATINV   C109.PTOLEMY.FORTLIB                         PER704
      SUBROUTINE MATINV (A,N,B,M,DETERM,NMAX)
C
C     MATRIX INVERSION WITH ACCOMPANYING SOLUTION OF LINEAR EQUATIONS
C
C     12/28/79 - DETMAX,-MIN FOR CVA,CUN - RPG                           2830   
C
      implicit real*8 ( a-h, o-z )                                      implicit
      REAL*4 PIVOT
      DIMENSION A(NMAX,N), B(NMAX,M)
      COMMON /F402/ PIVOT(50), INDEX(50)
C     (CARDS BELOW WITH ******** IN COLUMNS 73-80 IMPLEMENT THE IMPROVED
C     COMPUTATION OF THE DETERMINANT IN THE MODIFICATION OF MARCH 1975.)
      REAL*8 DETMAX,DETMIN
      DATA DETMAX, DETMIN / 1.D300, 1.D-300 /                           ci3e
CIB      DATA DETMAX,DETMIN/Z6110000000000000,Z2110000000000000/         2840   
CVA      DATA DETMAX, DETMIN / 1.D35, 1.D-35 /
C
C
C     INITIALIZE DETERMINANT AND PIVOT ELEMENT ARRAY
C
      DETERM=1.0
      IDET=0
      DO 20 I=1,N
      PIVOT(I)=0.0
   20 INDEX(I)=0.0                                                       2850   
C
C     PERFORM SUCCESSIVE PIVOT OPERATIONS (GRAND LOOP)
C
      DO 550 I=1,N
C
C     SEARCH FOR PIVOT ELEMENT AND EXTEND DETERMINANT PARTIAL PRODUCT
C
      AMAX=0.0
      DO 105 J=1,N
      IF (PIVOT(J).NE.0.0) GO TO 105                                     2860   
      DO 100 K=1,N
      IF (PIVOT(K).NE.0.0) GO TO 100
      TEMP=DABS(A(J,K))
      IF (TEMP.LT.AMAX) GO TO 100
      IROW=J
      ICOLUM=K
      AMAX=TEMP
  100 CONTINUE
  105 CONTINUE
      INDEX(I)=4096*IROW+ICOLUM                                          2870   
      J=IROW
      AMAX=A(J,ICOLUM)
      DETERM=AMAX*DETERM
      IF (DABS(DETERM).LT.DETMAX) GO TO 130
      DETERM=DETERM*DETMIN
      IDET=IDET+1
      GO TO 140
  130 IF (DABS(DETERM).GT.DETMIN) GO TO 140
      DETERM=DETERM*DETMAX
      IDET=IDET-1                                                        2880   
  140 CONTINUE
C
C     RETURN IF MATRIX IS SINGULAR (ZERO PIVOT) AFTER COLUMN INTERCHANGE
C
      IF (DETERM.EQ.0.0) GO TO 600
C
      PIVOT(ICOLUM)=AMAX
C
C     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
C                                                                        2890   
      IF (IROW.EQ.ICOLUM) GO TO 260
      DETERM=-DETERM
      DO 200 K=1,N
      SWAP=A(J,K)
      A(J,K)=A(ICOLUM,K)
      A(ICOLUM,K)=SWAP
  200 CONTINUE
      IF (M.LE.0) GO TO 260
      DO 250 K=1,M
      SWAP=B(J,K)                                                        2900   
      B(J,K)=B(ICOLUM,K)
      B(ICOLUM,K)=SWAP
  250 CONTINUE
C
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT
C
  260 K=ICOLUM
      A(ICOLUM,K)=1.0
      DO 350 K=1,N
      A(ICOLUM,K)=A(ICOLUM,K)/AMAX                                       2910   
  350 CONTINUE
      IF (M.LE.0) GO TO 380
      DO 370 K=1,M
      B(ICOLUM,K)=B(ICOLUM,K)/AMAX
  370 CONTINUE
C
C     REDUCE NON-PIVOT ROWS
C
  380 DO 550 J=1,N
      IF (J.EQ.ICOLUM) GO TO 550                                         2920   
      T=A( J,ICOLUM)
      A( J,ICOLUM)=0.0
      DO 450 K=1,N
      A( J,K)=A( J,K)-A(ICOLUM,K)*T
  450 CONTINUE
      IF (M.LE.0) GO TO 550
      DO 500 K=1,M
      B( J,K)=B( J,K)-B(ICOLUM,K)*T
  500 CONTINUE
  550 CONTINUE                                                           2930   
C
C     INTERCHANGE COLUMNS AFTER ALL PIVOT OPERATIONS HAVE BEEN PERFORMED
C
  600 DO 710 I=1,N
      I1=N+1-I
      K=INDEX(I1)/4096
      ICOLUM=INDEX(I1)-4096*K
      IF (K.EQ.ICOLUM) GO TO 710
      DO 705 J=1,N
      SWAP=A(J,K)                                                        2940   
      A(J,K)=A(J,ICOLUM)
      A(J,ICOLUM)=SWAP
  705 CONTINUE
  710 CONTINUE
C
      PIVOT(1)=IDET
      RETURN
      END
c*id* plmsub
      SUBROUTINE PLMSUB (LMAX, MMAXin, X, PRAY)                          2950   
C
C     COMPUTES THE ASSOCIATED LEGENDRE FUNCTIONS OF THE FIRST KIND:
C            P(L,M)(X)  FOR M = 0, MMAX  ,  L = M, LMAX  AND  REAL X
C
C       THE SIGN CONVENTION OF ABRAMOWITZ AND STEIGEN (M SUPER, L SUB)
C     AND OF JACKSON ("CLASSICAL ELECTRODYNAMICS") IS USED.
C     THERE IS NOT A (-1)**M BETWEEN P(L,M) AND Y(L,M) (WHERE Y IS
C     CONDON-SHORTLEY)  SO THAT  P(L,M) AND Y(L,M)
C     HAVE THE SAME SIGN.  NOTE THAT THIS IS DIFFERENT FROM
C     MESSIAH'S CONVENTION.                                              2960   
C
C       FOR  -1 <= X <= 1,  THE "VALUE ON THE CUT" GIVEN BY ABRAMOWITZ
C     AND STEGUN  8.3.1  IS COMPUTED.  THIS IS RELATED TO P(L,M)(Z) BY
C          P(L,M)(X)  =  EXP(I .5 M PI) * P(L,M)(Z = X + I 0)
C     FOR   X  > 1,  P(L,M)(Z=X) IS COMPUTED.  THUS THE OUTPUT IS REAL
C     FOR ALL REAL X.
C
C
C     The arguments are:
C                                                                        2970   
C     LMAX - THE MAXIMUM L VALUE TO BE FOUND.
C
C     MMAX - the maximum M value.  MMAX should satisfy
C            0 <= MMAX <= LMAX  but if  MMAX > LMAX, processing
C            will continue as if MMAX = LMAX  and no error message
C            will be printed.  However,  MMAX < 0  does result in
C            A FATAL ERROR MESSAGE.  NOTE THAT THE P(L,M) ARE GENERATED
C            ONLY FOR  M > 0.
C
C     X - THE VALUE OF X AT WHICH P(L,M) IS TO BE FOUND.  SEE ABOVE      2980   
C         FOR THE ACTION TAKEN DEPENDING ON THE RANGE OF X.
C
C     PLMRAY - A ONE DIMENSIONAL ARRAY OF DIMENSION AT LEAST
C                 ( (MMAX+1) * (2*LMAX - MMAX + 2) ) / 2
C              THAT WILL RECIEVE THE P(L, M).
C              THEY ARE STORED IN PLMRAY AS (READ L AS LMAX)
C        P00, P10, P20,... PL0, P11, P21, ... PL1, P22, P32, ... ... PLL
C              THUS
C                 P(L,M) = PLMRAY(L+1+ M*(2*LMAX+1-M)/2)  for  M >= 0
C                                                                        2990   
C
C     LMAX AND MMAX ARE INTEGER INPUT ARGUMENTS.
C     X IS A REAL INPUT ARGUMENT.
C     PLMRAY IS A REAL OUTPUT ARRAY.
C
C     This is the general-purpose version.  There is a more efficent
C     version for the Cray which is less accurate for small |X-1|.
C
C     9/25/72  NEW REVISED (IMPROVED) RECURSION SCHEME - S. PIEPER
C     10/23/74 - CORRECT OUTPUT FOR  X  > 1.                             3000   
C     8/1/80 - PROGRAM IT EFFICENTLY - S.P.
C
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
C
      DIMENSION PRAY(1)
C
C
      mmax = mmaxin                                                      3010   
C
      IF (MMAX .EQ. 0) GO TO 200
C
C
      ROOT = SQRT(ABS(1-X**2))
      IF ( ABS(X) .GT. 1 ) ROOT = -ROOT
C
C
 200  PRAY(1) = 1.
C                                                                        3020   
C
      I = 1
      J = 1
      DO 499 M = 0, MMAX
      IF (M .GE. LMAX) RETURN
      J = I
C     GENERATE THE P(L=M+1, M) FOR THIS M
      I = I+1
      PRAY(I) = (2*M+1)*X*PRAY(I-1)
      MP2 = M + 2                                                        3030   
      IF (LMAX .LT. MP2) GO TO 450
C
C     NOW DO L = M+2, M+3, ..., LMAX FOR THIS M
C
      DL = MP2
      TEMP1 = (2*DL-1)*X
      TEMP2 = M-1+DL
      TEMP3 = DL-M
      IST = I+1
      IEND = I+1 + LMAX-MP2                                              3040   
C
C
      DO 399 I = IST, IEND
      PRAY(I) = ( TEMP1*PRAY(I-1) - TEMP2*PRAY(I-2) ) / TEMP3
      TEMP1 = TEMP1 + (2.*X)
      TEMP2 = TEMP2 + 1.
      TEMP3 = TEMP3 + 1.
C
 399  CONTINUE
C                                                                        3050   
C
 400  I = IEND
C
 450  IF (M .EQ. MMAX) RETURN
C
C    NOW GENERATE THE NEXT P(L=M, M)
C    NOTE THAT HERE THE DESIRED M VALUE IS REALLY M+1 SINCE
C     THE DO LOOP HAS NOT YET INCREMENTED M
C
      I = I+1                                                            3060   
C     J POINTS TO THE LAST P(M, M)
      PRAY(I) = - (2*M+1)*ROOT * PRAY(J)
 499  CONTINUE
      RETURN
C
C     THIS ENTRY COMPUTES THE LEGENDRE POLYNOMIALS
C
      ENTRY PLSUB (LMAX, X, PRAY)
      MMAX = 0
      GO TO 200                                                          3070   
C
      END
c*id* rcwfn
      SUBROUTINE RCWFN ( RHO,ETA,MINL,MAXL,FC,FCP,GC,GCP,ACCUR, IRET )
C
C     COMPUTES REGULAR AND IRREGULAR COULOMB WAVEFUNCTIONS.
C
C   *****************************************************************
C   *                                                               *
C   * !!! THE APPROPRIATE DATA STATEMENTS MUST BE ACTIVATED !!!!    *    3080   
C   *   ci3e FOR IEEE machines (RS6000, ...)                        *
C   *   CIB FOR IBM 370 TYPE MACHINES                               *
C   *   CVA FOR VAX TYPE MACHINES WITH D-FORMAT                     *
C   *   CRA FOR CRAY TYPE MACHINES - DOUBLE PRECISION MUST BE       *
C   *          DISABLED IN THE COMPLILER INVOCATION                 *
C   *   CDC FOR CDC CYBER TYPE MACHINES  - DOUBLE PREC. FUNCTIONS   *
C   *          MUST BE CHANGED TO SINGLE                            *
C   *   CUN FOR UNIVAC 1108 TYPE MACHINES                           *
C   *                                                               *
C   *****************************************************************    3090   
C
C       THIS SUBROUTINE RETURNS THE REGULAR AND IRREGULAR COULOMB
C     WAVEFUNCTIONS AND THEIR DERIVATIVES FOR SPECIFIC VALUES OF
C     ETA AND RHO AND FOR A RANGE OF L'S.  THE RANGE OF L'S (ORBITAL
C     ANGULAR MOMENTUM) NEED NOT BEGIN AT L = 0 AND, IN GENERAL,
C     IT WILL BE LESS TIME CONSUMING TO REQUEST ONLY THOSE L'S
C     THAT ARE ACTUALLY DESIRED.  THE CONVENTIONS OF ABRAMOWITZ
C     AND STEGUN (HANDBOOK OF MATHEMATICAL FUNCTIONS) ARE USED.
C
C       THIS SUBROUTINE IS AN ADAPTATION OF THE MANCHESTER SUBROUTINE    3100   
C     OF THE SAME NAME (BUT SLIGHTLY DIFFERENT ARGUMENT LIST) THAT
C     WAS PUBLISHED IN COMPUTER PHYSICS COMMUNICATIONS 8, 377 (1974).
C     BY BARNETT, FENG, STEED AND GOLDFARB.  THE CONTINUED FRACTIONS
C     INTRODUCED IN THAT ARTICLE ARE USED IN THE RANGE RHO >
C     RHO(TURN) BUT COMPLETELY DIFFERENT PROCEDURES ARE USED
C     FOR RHO < RHO(TURN).  THE NEW PROCEDURES ALLOW RESULTS
C     FOR RHO < RHO(TURN) OF ACCURACY COMPARABLE TO THOSE FOR
C     RHO > RHO(TURN) IN COMPARABLE TIME.
C
C       ONE OF THE CONTINUED FRACTIONS IS NOT CONVERGENT FOR VERY        3110   
C     SMALL RHO AND THUS THIS ROUTINE WILL FIND ONLY F AND F' FOR
C     RHO < .005.  IN ADDITION THE ROUTINE BECOMES RATHER SLOW FOR
C     RHO < .1.  FOR VERY LARGE RHO ( >> 1000 ), ANOTHER CONTINUED
C     FRACTION CONVERGES VERY SLOWLY AND THE SUBROUTINE WILL FAIL
C     TO CONVERGE IN THE MAXIMUM ALLOWED NUMBER OF ITERATIONS FOR
C     RHO > 9500 WITH THE EXCEPTION THAT LARGER RHO ARE POSSIBLE
C     IF ETA IS NEAR RHO/2.
C
C       ASIDE FROM THE ABOVE LIMITATIONS, THE SUBROUTINE HAS A VERY
C     LARGE RANGE OF ETA,  RHO, AND L FOR WHICH IT WILL RETURN           3120   
C     RELIABLE VALUES IN REASONABLE TIMES.  IT HAS BEEN TESTED FOR
C     -2000 < ETA < 5000,  1E-6 < RHO < 10000;  0 < L < 1000 .
C     THE TESTS CONSISTED OF COMPARASONS TO EXISTING ROUTINES TO
C     LOCATE CODING ERRORS AND OF COMPARASONS TO QUADRUPLE PRECISION
C     RESULTS TO DETERMINE THE NUMERICAL STABILITY OF THE ALGORITHMS.
C     HOWEVER, NO VERIFICATION OF THE CORRECTNESS OF G OR G' FOR
C     ETA < 0 HAS BEEN MADE.
C
C
C     ARGUMENTS (ALL FLOATING POINT ARGUMENTS ARE DOUBLE PRECISION) -    3130   
C
C     RHO - THE VALUE OF THE RADIAL COORDINATE AT WHICH F, F', G, AND
C           G' ARE DESIRED.  0 < RHO < 9500 IS REQUIRED.  IF
C           RHO < .005, ONLY F AND F' WILL BE FOUND.
C
C     ETA - THE VALUE OF THE CHARGE PARAMETER TO BE USED.  THE ROUTINE
C           HAS BEEN TESTED FOR -2000 < ETA < 5000.  ETA = 0 WILL
C           RESULT IN THE COMPUTATION OF THE SPHERICAL BESSEL
C           FUNCTIONS TIMES RHO.
C                                                                        3140   
C     MINL - THE MINIMUM VALUE OF L FOR WHICH THE FUNCTIONS ARE DESIRED.
C     MAXL - THE MAXIMUM VALUE OF L.
C
C     FC, FCP, GC, AND GCP - THESE ONE DIMENSIONAL ARRAYS WILL BE
C           SET TO THE COMPUTED VALUES OF F, F', G, AND G' RESPECTIVELY.
C           EACH ARRAY MUST BE OF LENGTH AT LEAST MAXL+1 AND WILL
C           BE SET AS
C             ARRAY(L+1) = FUNCTION(L)   MINL =< L =< MAXL .
C           DEPENDING ON ETA AND RHO, THE FIRST MINL ELEMENTS OF
C           EACH ARRAY MAY BE USED FOR INTERMEDIATE COMPUTATIONS AND     3150   
C           THEIR VALUES ARE NOT PREDICTABLE UPON EXIT FROM RCWFN.
C
C     ACCUR - THE DESIRED ACCURACY OF CONVERGENCE OF THE CONTINUED
C          FRACTIONS AND OTHER SERIES.  IN GENERAL THIS WILL BE
C          THE RELATIVE ACCURACY OF THE FINAL VALUES EXCEPT NEAR
C          A ZERO OF ONE OF THE FUNCTIONS.  ACCUR SHOULD BE IN
C          THE RANGE
C              1E-6 < ACCUR < 1E-15 ,
C          IF IT IS NOT, ONE OF THE ABOVE TWO VALUES WILL BE USED.
C                                                                        3160   
C     IRET - THIS INTEGER IS SET TO A RETURN CODE TO INDICATE THAT
C          RCWFN WAS SUCCESSFUL OR THE REASON FOR FAILURE:
C        0 - ALL O.K.; ALL FUNCTIONS FOUND.
C        1 - SOME O.K., THOSE FOR HIGHER L'S HAVE UNDER/OVER-FLOWED
C        2 - RHO < .005; ONLY F, F' WERE COMPUTED.
C        3 - RHO < .005; ONLY F, F' WERE COMPUTED AND SOME OF THEM
C                        UNDERFLOWED.
C        4 - ALL WILL UNDER/OVER-FLOW; NO RESULTS DEFINED.
C        5 - FAILED TO AVOID DIVIDE BY ZERO IN F'/F LOOP
C        6 - NONCONVERGENCE OF R                                         3170   
C        7 - NONCONVERGENCE OF P + IQ
C        8 - NONCONVERGENCE OF MACLAUREN SERIES
C        9 - NONCONVERGENCE OF THE TAYLOR SERIES ABOUT RHO = 2*ETA
C     IRET'S OF 5 TO 9 SHOULD BE REPORTED TO STEVE PIEPER ( X4523 )
C
C       THE FOLLOWING TABLE GIVES SAMPLE EXECUTION TIMES ON THE
C     /75 AND /195.  TIMES ARE FOR ACCUR = 1D-14 AND FOR
C     MINL = MAXL = 0.  IN GENERAL EXECUTION TIME IS WEAKLY DEPENDANT
C     ON ACCUR AND MINL.  IF MAXL > MINL, THE TIMES REQUIRED FOR
C     EACH ADDITIONAL L ARE APPROXIMATLEY 0.01 MILLISEC ON THE /195      3180   
C     AND 0.15 MILLISEC ON THE /75.  TIMES FOR ETA < 0 ARE COMPARABLE TO
C     TO THOSE FOR \ETA\.  THE FINAL COLUMN OF THE TABLE GIVES
C     THE NUMBER OF BITS LOST DUE TO ROUND OFF AND TRUNCATION
C     ERRORS AND INDICATES THE MAXIMUM PRECISION POSSIBLE ON A
C     GIVEN MACHINE.  THE /360 AND /370 HAVE 56 BITS OF PRECISION
C     (16.8 DECIMAL PLACES) AND EACH 3.3 BITS LOST REPRESENTS
C     ONE DECIMAL PLACE LOST.
C
C      ETA    RHO     TIME IN MILLISECONDS    BITS LOST
C                       /195      /75                                    3190   
C
C        0.      .01     .15       1.2            3
C        0.     1.       .13       1.2            3
C        0.   100.      1.0       11.            11
C        0.  1000.      7.6       81.
C        1.      .1     5.4       81.             9
C        1.     1.      1.0       12.5            8
C        1.    10.       .33       3.7            3
C        1.   100.      1.1       11.
C        1.  1000.      8.1       81.                                    3200   
C       10.     1.      1.5       20.             6
C       10.    10.       .8        8.3            6
C      100.    75.      1.9       21.            10
C      100.   100.      2.0       21.            10
C      100.  1000.      7.3       74.            13
C     1000.  5000.     29.       385.            16
C
C
C     MAY 23, 1976 - REVISED VERSION BY S. PIEPER.
C     NOV 19, 1976 - FIX CHOICE OF ROOT FOR RHO << ETA, ETA > 10         3210   
C     12/14/77 - CDC VERSION USING MORTRAN
C     6/27/78 - FIX CDC VRYBIG, PRINTOUTS FOR ERROR CONDITIONS - S.P.
C     1/29/80 - CUN, CVA VERSIONS, CFR%F - RPG
C     3/13/80 - BACK TO CFR%F - S.P.
C     7/9/80 - CRAY DATA STATEMENST - S.P.
C     7/3/91 - RS6000 version
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      DIMENSION FC(1),FCP(1),GC(1),GCP(1)                                3220   
C
C     VRYBIG IS THE RESULT OF AN OVERFLOW
C     BIG IS REPRESENTATIVE OF MAX NUMBER, SMALL IS ITS INVERSE
C     SMALLN IS  LOG(SMALL)  (MUST BE ACCURATE)
C     PRECIS IS ABOUT  100*(MACHINE PRECISION)
C     PRECLN IS  > \LOG(MACHINE PRECISION)\
C     PRERT3 IS LIKE THE CUBEROOT OF THE MACHINE PRECISION
C
CIB      DATA VRYBIG / Z 7FFFFFFF FFFFFFFF /
CIB      DATA BIG / Z 70800000 00000000 /                                3230   
CIB      DATA SMALL / Z 11200000 00000000 /
C     REAL*16 SMALLN / -132.39111148694955409869133519851171 Q+00 /
CIB      REAL*8 SMALLN / Z C284641FE1E589CC /
CIB      DATA PRERT3 / 1.D-7 /,  PRECIS / 1.D-15 /,  PRECLN / 44.D0 /
C
C     FOR CDC WE STOP BEFORE OVERFLOW SINCE CANNOT TRAP OVERFLOW
CDC      DATA VRYBIG / 3657 4000 0000 0000 0000 B /
CDC      DATA BIG / 3621 4000 0000 0000 0000 B /
CDC      DATA SMALL  / 0020 4000 0000 0000 0000 B /
CDC      DATA SMALLN / 6046 2632 2411 4172 5514 B /                      3240   
CDC      DATA PRERT3 / 1.E-5 /,  PRECIS /1.E-12 /,
CDC     1  PRECLN / 35.E0 /
C
CRA      DATA VRYBIG / 1.E+2450 /
CRA      DATA BIG    / 1.E+2100 /
CRA      DATA SMALL  / 1.E-2100 /
CRA      DATA SMALLN / -4835.4286 9528 749 /
CRA      DATA PRERT3 / 1.E-5 /,  PRECIS / 1.E-12 /,
CRA     1  PRECLN / 35.E0 /
C                                                                        3250   
      DATA VRYBIG / 1.79769D+308 /                                      ci3e
      DATA BIG    / 1.D+300 /                                           ci3e
      DATA SMALL  / 1.D-300 /                                           ci3e
      DATA SMALLN / -690.775527898214 /                                 ci3e
      DATA PRERT3 / 2.7E-5 /,  PRECIS / 2.E-14 /,                       ci3e
     1  PRECLN / 32.E0 /                                                ci3e
C
CUN      DATA VRYBIG / 1.D290 /, BIG / 1.D250 /, SMALL / 1.D-250 /
CUN      DATA SMALLN / -575.64627 32487 61421 0D0 /
CUN      DATA PRERT3 / 1.D-7 /, PRECIS / 1.D-16 /, PRECLN / 44.D0 /      3260   
C
CVA      DATA VRYBIG / 1.D32 /, BIG / 1.D25 /, SMALL / 1.D-25 /
CVA      DATA SMALLN / -57.564627 32487 61421 0D0 /
CVA      DATA PRERT3 / 1.D-7 /, PRECIS / 1.D-15 /, PRECLN / 44.D0 /
C
CXX      REAL*16 PI / 3.1415926535897932384626433832795028 0Q+00 /
      REAL*8  PI / 3.141592653589793238                 0D+00 /
CIB      DATA PI / Z 413243F6A8885A31 /
C
      LOGICAL FRSTSW                                                     3270   
C
C
C     COULOMB WAVEFUNCTIONS CALCULATED AT RHO = RHO BY THE
C     CONTINUED-FRACTION METHOD OF STEED   MINL,MAXL ARE ACTUAL L-VALUES
C     SEE BARNETT FENG STEED AND GOLDFARB COMPUTER PHYSICS COMMON 1974
C
C
C     HERE WE LIMIT ACCURACY TO REASONABLE VALUES FOR THE MACHINE
C
      ACC  = ACCUR                                                       3280   
      ACC = DMAX1( ACC, PRECIS )
      ACC = DMIN1( ACC, PRERT3 )
C
      LMAX = MAXL
      LMIN = MINL
      LMIN1= LMIN + 1
      XLL1 = LMIN*LMIN1
      ETA2 = ETA*ETA
C
C     DETERMIN WHICH REGION WE ARE IN                                    3290   
C
C     FOR RHO < .45, Q OF P+IQ IS POORLY DETERMINED SO WE DON'T USE IT.
C     EXCEPT THAT FOR LARGE NEGATIVE ETA THE MACLAUREN SERIES ALSO
C     HAS PROBLEMS
C
      IF ( RHO .GT. .45D0 )  GO TO 20
      IF ( ETA .GE. 0 )  GO TO 10
      IF ( -ETA*RHO .GT. 7 )  GO TO 20
C
C     FOR RHO < .005, WE ONLY RETURN F AND F' SINCE THE P+IQ RECURSION   3300   
C     IS VERY SLOWLY CONVERGENT  ( FOR VERY SMALL ETA IT IS POSSIBLE
C     TO GO TO SMALLER RHO ( RHO > 1E+4*ETA )  BUT WE IGNOR THAT HERE.
C
 10   IF ( RHO .GT. .005D0 )  GO TO 60
      IGOTO = 5
      GO TO 70
C
 20   TURN = ETA + DSQRT(ETA2 + XLL1)
      IGOTO = 1
      IF ( RHO .GE. TURN-1.D-4 )  GO TO 100                              3310   
C
C     WE ARE INSIDE THE TURNING POINT FOR MINL, CAN WE GET OUTSIDE
C     OF IT BY REDUCING MINL. (THIS IS ALWAYS POSSIBLE FOR
C     ETA < 0).
C
      IF ( RHO .LT. ETA+DABS(ETA) )  GO TO 60
C
C     YES, DO SO - THIS IS THE SAME AS  RHO > RHO(TURN)  EXCEPT
C     WE GENERATE SOME EXTRA F(L), G(L) FOR  L < MINL
C                                                                        3320   
      LMIN = .5*( DSQRT(1+4*((RHO-ETA)**2-ETA2)) - 1 )
      LMIN1 = LMIN + 1
      GO TO 80
C
C     MUST USE A DIFFERENT METHOD TO SUPPLIMENT THE BAD I Q
C     VALUE.  ALWAYS START WITH LMIN = 0 FOR SIMPLICITY.
C
C     NOTE ONLY ETA > 0 GETS TO HERE ( EXCEPT WHEN RHO < .45 )
C
 60   IGOTO = 2                                                          3330   
 70   LMIN = 0
      LMIN1 = 1
      IF ( ETA .LT. 10  .OR.  RHO .LE. ETA )  GO TO 80
      IGOTO = 3
C
 80   XLL1 = LMIN*LMIN1
C
C     HERE WE COMPUTE  F'/F  FOR L = MAXL
C     WE THEN RECURSE DOWN TO LMIN TO GENERATE THE UNNORMALIZED F'S
C     THIS SECTION IS USED FOR ALL RHO.                                  3340   
C
 100  PL   = LMAX + 1
      RHOUSE = RHO
 105  PLSAVE = PL
 110  FRSTSW = .TRUE.
C     CONTINUED FRACTION FOR  R = FP(MAXL)/F(MAXL)
      R  = ETA/PL + PL/RHOUSE
      DQ  = (ETA*RHOUSE)*2.0 + 6*PL**2
      DR = 12*PL + 6
      DEL = 0.0                                                          3350   
      D   = 0.0
      F   = 1.0
      X   = (PL*PL - PL + (ETA*RHOUSE))*(2.0*PL - 1.0)
      AI = RHOUSE*PL**2
      DI = (2*PL+1)*RHOUSE
C
C     LOOP AND CONVERGE ON R
C
      DO 139  I = 1, 100000
         H = (AI + RHOUSE*ETA2)*(RHOUSE - AI)                            3360   
         X   = X + DQ
         D = D*H + X
C
C     IF WE PASS NEAR A ZERO OF THE DIVISOR, START OVER AT
C     LARGER LMAX
C
         IF ( DABS(D) .GT. PRERT3*DABS(DR) )  GO TO 130
         PL = PL + 1
         IF ( PL .LT. PLSAVE+10 )  GO TO 110
         IRET = 5                                                        3370   
      GO TO 990
C
 130     D = 1/D
         DQ = DQ + DR
         DR = DR + 12
         AI = AI + DI
         DI = DI + 2*RHOUSE
         DEL =  DEL*(D*X - 1.0)
         IF (FRSTSW) DEL = -RHOUSE*(PL*PL + ETA2)*(PL + 1.0)*D/PL
         FRSTSW = .FALSE.                                                3380   
         R  = R + DEL
         IF(D.LT.0.0) F = -F
         IF ( DABS(DEL) .LT. DABS(R*ACC) )  GO TO 140
 139  CONTINUE
      IRET = 6
      GO TO 990
C
C     R HAS CONVERGED;  DID WE INCREASE LMAX
C
 140  IF ( PL .EQ. PLSAVE )  GO TO 160                                   3390   
C
C     RECURSE DOWN ON R TO LMAX
C     HERE THE ONLY PART OF F THAT IS OF INTEREST IS THE SIGN
C
      PL = PL-1
 150     D = ETA/PL + PL/RHOUSE
         F = (R+D)*F
         R = D - (1+ETA2/PL**2)/(R+D)
         PL = PL - 1
         IF ( PL .GT. PLSAVE )  GO TO 150                                3400   
C
C     NOW HAVE R(LMAX, RHO) OR IF IGOTO=4, R(LMIN, 2*ETA)
C
 160  IF ( IGOTO .EQ. 4 )  GO TO 210
      FC (LMAX+1) = F
      FCP(LMAX+1) = F*R
      IF( LMAX.EQ.LMIN) GO TO 200
C     DOWNWARD RECURSION TO LMIN FOR F AND FP, ARRAYS GC,GCP ARE STORAGE
      L  = LMAX
      PL = LMAX                                                          3410   
      AR = 1/RHO
      DO 189 LP  = LMIN1,LMAX
         GC (L+1) = ETA/PL + PL*AR
         GCP(L+1) = DSQRT( (ETA/PL)**2 + 1 )
         FC (L)   = (GC(L+1)*FC(L+1) + FCP(L+1))/GCP(L+1)
         FCP(L)   =  GC(L+1)*FC(L)   - GCP(L+1)*FC(L+1)
         PL = PL - 1
         L  = L - 1
C
C     IF WE ARE GETTING NEAR AN OVERFLOW, RENORMALIZE EVERYTHING DOWN    3420   
C
         IF ( DABS(FC(L+1)) .LT. BIG )  GO TO 189
         DO 179  LL = L, LMAX
            FC(LL+1) = SMALL*FC(LL+1)
            FCP(LL+1) = SMALL*FCP(LL+1)
 179     CONTINUE
 189  CONTINUE
      F  = FC (LMIN1)
      R = FCP(LMIN1)/F
C                                                                        3430   
C     HERE WE FIND
C        P + IQ  =  (G'+IF')/(G+IF)
C     THIS SECTION IS USED IN ALL CASES EXCEPT WHEN
C        15 < ETA < RHO < 2*ETA
C
 200  IF ( IGOTO .EQ. 3 )  GO TO 500
      IF ( IGOTO .EQ. 5 )  GO TO 400
C
C     NOW OBTAIN P + I.Q FOR LMIN FROM CONTINUED FRACTION (32)
C     REAL ARITHMETIC TO FACILITATE CONVERSION TO IBM USING REAL*8       3440   
 210  P  = 0.0
      Q  = RHOUSE - ETA
      PL = 0.0
      AR = -(ETA2 + XLL1)
      AI =   ETA
      BR = Q + Q
      BI = 2.0
      WI = ETA + ETA
      DR =   BR/(BR*BR + BI*BI)
      DI =  -BI/(BR*BR + BI*BI)                                          3450   
      DP = -(AR*DI + AI*DR)
      DQ =  (AR*DR - AI*DI)
C
C     LOOP AND CONVERGE ON P + IQ
C
 230     P  =  P + DP
         Q  =  Q + DQ
         PL = PL + 2.0
         AR = AR + PL
         AI = AI + WI                                                    3460   
         BI = BI + 2.0
         D  = AR*DR - AI*DI + BR
         DI = AI*DR + AR*DI + BI
         T  = 1.0/(D*D + DI*DI)
         DR =  T*D
         DI = -T*DI
         H  = BR*DR - BI*DI - 1.0
         X  = BI*DR + BR*DI
         T  = DP*H  - DQ*X
         DQ = DP*X  + DQ*H                                               3470   
         DP = T
         IF(PL.GT.46000.) GO TO 920
         IF(DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC) GO TO 230
      P  = P/RHOUSE
      Q  = Q/RHOUSE
C
C     WE NOW HAVE  R  AND  P+IQ,  IS THIS ENOUGH
C
      IF ( IGOTO .EQ. 2 )  GO TO 400
C                                                                        3480   
C     SOLVE FOR FP,G,GP AND NORMALISE F  AT L=LMIN
C
C     SINCE THIS IS FOR  RHO > RHO(TURN), F AND G ARE REASONABLE
C     NUMBERS
C
      X = (R-P)/Q
      FMAG = DSQRT( 1/(Q*(1+X**2)) )
      W = FMAG/DABS(F)
      F = W*F
      G = F*X                                                            3490   
      GP = R*G - 1/F
      IF ( IGOTO .EQ. 4 )  GO TO 600
      GO TO 800
C
C     HERE   RHO < ETA  OR  RHO < 2*ETA < 20  OR  RHO < .45
C     WE USE THE MACLAUREN SERIES TO GET  F( L=0, ETA, RHO )
C
C     FIRST COMPUTE  RHO*C(L=0, ETA)
C
 400  C = 2*PI*ETA                                                       3500   
      IF ( DABS(C) .GT. .5 )  GO TO 410
C
C     USE MACLAURIN EXPANSION OF  X / (EXP(X)-1)
C
      X = 0
      T = 1
      AR = 1
      BR = C
      AI = 1
      C = 1                                                              3510   
 405     AI = AI + 1
         AR = AR*BR/AI
         C = C + AR
         IF ( DABS(AR) .GE. ACC*C )  GO TO 405
      C = 1/C
      GO TO 430
C
C     HERE ETA IS NOT TINY.
C
 410  IF ( ETA .GT. 0 )  GO TO 420                                       3520   
      C = -C
      X = 0
      T = 1
      GO TO 425
 420  X = -SMALLN - PI*ETA
      T = SMALL
 425  IF ( C .LT. PRECLN )  C = C / (1-DEXP(-C))
 430  C = RHO*DSQRT(C)
      B1 = 1
      B2 = ETA*RHO                                                       3530   
      SUM = B1 + B2
      AI = 6
      DI = 6
      DO 449  I = 1, 10000
         B3 = ( (2*ETA*RHO)*B2 - (RHO**2)*B1 ) / AI
         AI = AI + DI
         DI = DI + 2
         SUM = SUM + B3
         STOP = DABS(B1) + DABS(B2) + DABS(B3)
         B1 = B2                                                         3540   
         B2 = B3
         IF ( DABS(SUM) .LT. BIG )  GO TO 445
         X = X - SMALLN
         SUM = SUM*SMALL
         B1 = B1*SMALL
         B2 = B2*SMALL
 445     IF ( STOP .LT. ACC*DABS(SUM) )  GO TO 450
 449  CONTINUE
      IRET = 8
      GO TO 990                                                          3550   
C
 450  SUM = ( C*DEXP(X)*SUM ) * T
C
C     DID IT UNDERFLOW
C
      IF ( SUM .EQ. 0 )  GO TO 900
C
C     WE NOW HAVE  F (=SUM),  R,  AND P  ( P ONLY IF RHO > .005 )
C     USE THE WRONSKIAN AS THE 4TH CONDITION
C                                                                        3560   
      W = SUM/F
      F = SUM
      IF ( IGOTO .EQ. 5 )  GO TO 850
      X = (R-P)*F
      IF ( DABS(X) .GT. PRERT3 )  GO TO 480
C
C     HERE F**3 AND F**4 TERMS ARE LESS THAN MACHINE PRECISION
C
      G = 1/X
      GP = P*G                                                           3570   
      GO TO 800
C
C     HERE WE MUST INCLUDE F**3, F**4; WE MUST ALSO WORRY ABOUT
C     WHICH SIGN OF THE ROOT TO USE.
C     THE POSITIVE ROOT APPLIES FOR G > F; ELSE THE NEGATIVE ROOT
C     WE USE Q IN DETERMINING WHICH IS CORRECT
C
 480  B1 = .5/X
      B2 = B1 * DSQRT(1-4*(X*F)**2)
      G = B1 + B2                                                        3580   
C     G > F IN ALL OF REGION 2 FOR ETA > 0
      IF ( ETA .GE. 0 )  GO TO 490
      SUM = 1/Q - F**2
      GP = B1 - B2
      IF ( DABS(G**2-SUM) .GT. DABS(GP**2-SUM) )  G = GP
 490  GP = P*G - X*F/G
      GO TO 800
C
C     ETA > 15  AND  ETA < RHO < 2*ETA
C                                                                        3590   
C     WE FIND G AND G' FOR LMIN, RHO=2*ETA USING THE ABOVE METHOD
C     CONSISTING OF R, P+IQ, AND W.
C
 500  RHOUSE = ETA+ETA
      PL = LMIN+1
      IGOTO = 4
      GO TO 105
C
C     NOW WE HAVE  G, G'  AT THE TURNING POINT, GO IN USING TAYLOR
C                                                                        3600   
C
 600  DEL = RHOUSE - RHO
      B1 = G
      B2  = -DEL*GP
      B3 = 0
      G = B1+B2
      ACCR = ACC/2
      DELINV = -1/DEL
      DFACTR = 3*DELINV
      X = DEL/RHOUSE                                                     3610   
      AI = X+X
      DI = AI+AI
      AR = 6
      DR = 6
      DO 639  I = 1, 10000
         S = ( AI*B3 + (X*DEL**2)*B1 ) / AR
         AR = AR + DR
         DR = DR + 2
         AI = AI + DI
         DI = DI + 2*X                                                   3620   
         G = G + S
         GP = GP + DFACTR*S
         IF ( G .GE. VRYBIG )  GO TO 900
         DFACTR = DFACTR + DELINV
         B1 = B2
         B2 = B3
         B3 = S
         IF ( S .LT. ACCR*G )  GO TO 650
 639  CONTINUE
      IRET = 9                                                           3630   
      GO TO 990
C
C     HERE WE HAVE  R = F'/F,  G,  G'
C     USE WRONSKIAN AS THE 4TH CONDITION
C
 650  F = FC(LMIN1)
      R = FCP(LMIN1)/F
      SUM = 1/(R*G-GP)
      W = SUM/F
      F = SUM                                                            3640   
C
C     WE NOW HAVE  F, R = F'/F, G, G'  AT LMIN
C
C     UPWARD RECURSION FROM GC(LMIN) AND GCP(LMIN),STORED VALUES ARE RHO
C     RENORMALISE FC,FCP FOR EACH L-VALUE
C
 800  GC (LMIN1) = G
      GCP(LMIN1) = GP
      FC(LMIN1) = F
      FCP(LMIN1) = R*F                                                   3650   
      IRET = 0
      IF(LMAX.EQ.LMIN)  RETURN
      DO  829  L = LMIN1,LMAX
         T        = GC(L+1)
         GC (L+1) = (GC(L)*GC (L+1) - GCP(L))/GCP(L+1)
         GCP(L+1) =  GC(L)*GCP(L+1) - GC(L+1)*T
         FC (L+1) = W*FC (L+1)
 829     FCP(L+1) = W*FCP(L+1)
 840  IF ( DABS(FC(LMAX+1))+DABS(FCP(LMAX+1)) .EQ. 0 )  IRET = IRET+1
      RETURN                                                             3660   
C
C     RHO < .005;  WE CANNOT FIND P OR Q AND SO RETURN ONLY F, F'.
C
 850  FC(LMIN1) = F
      FCP(LMIN1) = R*F
      IRET = 2
      IF ( LMAX .EQ. LMIN )  RETURN
      DO 859  L = LMIN1, LMAX
         FC(L+1) = W*FC(L+1)
         FCP(L+1) = W*FCP(L+1)                                           3670   
 859  CONTINUE
      GO TO 840
C
C     F AND G ARE OUT OF THE MACHINE EXPONENT RANGE FOR LMIN.
C     IT WILL BE EVEN WORSE FOR  L > LMIN  SO GIVE UP AND RETURN
C
 900  IRET = 4
      GO TO 990
C
C     P + IQ FAILED TO CONVERGE                                          3680   
 920  IRET = 7
C
C     ERROR DETECTED, PRINT THE INPUT AND SOME STUFF
C
 990  write (6, 993) IRET, RHO, ETA, MINL, MAXL, ACCUR,
     1  RHOUSE, P, Q, R, T, X, I, SUM
 993  FORMAT ( 17H0***RCWFIN IRET =, I5, 6X, 7HINPUT =,
     1    2G20.10, 2I10, G15.5 /
     2  4H ***, 6G18.8 /  I8, G18.8 )
      RETURN                                                             3690   
C
      END
C*ID* ROOTIS   C109.PTOLEMY.FORTLIB                         PER704  21:
      BLOCK DATA
C
C     BLOCK DATA FOR THE  ROOTIS  COMMON BLOCK USED BY  YLMSUB.
C
C     T___  BLOCK DATA  _______ _________ _ ______ _____
C     CONTAINING THE SQUARE ROOTS OF THE INTEGERS 0 - 601.  THE
C     _________ __ ___ ______ _____ __:                                  3700   
C
C       COMMON  /ROOTIS/  MAX, ISPACE, ROOTS(602)
C
C     _____
C
C     MAX  __ __ _______ ____ __ ___ __ ___ _______ _______ _____
C          SQUARE ROOT IS IN ROOTS (AT PRESENT MAX = 601).
C
C     ISPACE __ __ ______ _______.
C                                                                        3710   
C     ROOTS __ _ REAL*8 _____ ___ __
C           ROOTS(_+1) = ____(_),  0 <= _ <= MAX.
C
C       T___ BLOCK DATA _______ ____ ________ _ _____ ______ _____
C     _______  /ROOTIX/  _____ ___ __ ____ __ _____ ___ _______ __
C     ___  /ROOTIS/  _____ __ _________ _ _________ __ ___ ____
C        EXTERNAL  ROOTIX
C     __ ___ _______ ____ _______ __ ____ ___ __ ___  /ROOTIS/
C     ______ _____.
C                                                                        3720   
C       T__  /ROOTIS/  ______ _____ __ ____ __  YLMSUB  ___  YLM2D.
C
C     4/21/74 - S. P_____.
C     7/31/75 - EXPANDED TO 601 INTEGERS FOR PTOLEMY.
C
      implicit real*8 ( a-h, o-z )                                      implicit
C
      COMMON  /ROOTIS/  IMAX, ISPACE, ROOTI(602)
      COMMON /ROOTIX/ X
C                                                                        3730   
      DATA IMAX /601/
C
C     THE DECIMAL FORMS OF THE FOLLOWING DATA STATEMENTS MAY BE USED
C     ON NON /360 COMPUTERS.
C
      DATA  ROOTI(  1) /       0.0                      /               cni
CIB      DATA  ROOTI(  1) / Z0000000000000000 /
      DATA  ROOTI(  2) /       0.100000000000000000 0D+01 /             cni
CIB      DATA  ROOTI(  2) / Z4110000000000000 /
      DATA  ROOTI(  3) /       0.141421356237309515 0D+01 /              3740   
CIB      DATA  ROOTI(  3) / Z4116A09E667F3BCD /
      DATA  ROOTI(  4) /       0.173205080756887719 0D+01 /             cni
CIB      DATA  ROOTI(  4) / Z411BB67AE8584CAA /
      DATA  ROOTI(  5) /       0.200000000000000000 0D+01 /             cni
CIB      DATA  ROOTI(  5) / Z4120000000000000 /
      DATA  ROOTI(  6) /       0.223606797749978981 0D+01 /             cni
CIB      DATA  ROOTI(  6) / Z4123C6EF372FE950 /
      DATA  ROOTI(  7) /       0.244948974278317810 0D+01 /             cni
CIB      DATA  ROOTI(  7) / Z4127311C2812425D /
      DATA  ROOTI(  8) /       0.264575131106459049 0D+01 /              3750   
CIB      DATA  ROOTI(  8) / Z412A54FF53A5F1D3 /
      DATA  ROOTI(  9) /       0.282842712474619007 0D+01 /             cni
CIB      DATA  ROOTI(  9) / Z412D413CCCFE7799 /
      DATA  ROOTI( 10) /       0.300000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 10) / Z4130000000000000 /
      DATA  ROOTI( 11) /       0.316227766016837930 0D+01 /             cni
CIB      DATA  ROOTI( 11) / Z413298B075B4B6A5 /
      DATA  ROOTI( 12) /       0.331662479035539981 0D+01 /             cni
CIB      DATA  ROOTI( 12) / Z413510E527FADE68 /
      DATA  ROOTI( 13) /       0.346410161513775461 0D+01 /              3760   
CIB      DATA  ROOTI( 13) / Z41376CF5D0B09955 /
      DATA  ROOTI( 14) /       0.360555127546398935 0D+01 /             cni
CIB      DATA  ROOTI( 14) / Z4139B05688C2B3E7 /
      DATA  ROOTI( 15) /       0.374165738677394133 0D+01 /             cni
CIB      DATA  ROOTI( 15) / Z413BDDD422D07E92 /
      DATA  ROOTI( 16) /       0.387298334620741680 0D+01 /             cni
CIB      DATA  ROOTI( 16) / Z413DF7BD629E9DB3 /
      DATA  ROOTI( 17) /       0.400000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 17) / Z4140000000000000 /
      DATA  ROOTI( 18) /       0.412310562561766059 0D+01 /              3770   
CIB      DATA  ROOTI( 18) / Z4141F83D9ABFB41C /
      DATA  ROOTI( 19) /       0.424264068711928521 0D+01 /             cni
CIB      DATA  ROOTI( 19) / Z4143E1DB337DB366 /
      DATA  ROOTI( 20) /       0.435889894354067353 0D+01 /             cni
CIB      DATA  ROOTI( 20) / Z4145BE0CD19137E2 /
      DATA  ROOTI( 21) /       0.447213595499957939 0D+01 /             cni
CIB      DATA  ROOTI( 21) / Z41478DDE6E5FD29F /
      DATA  ROOTI( 22) /       0.458257569495584005 0D+01 /             cni
CIB      DATA  ROOTI( 22) / Z4149523AE4547A15 /
      DATA  ROOTI( 23) /       0.469041575982342951 0D+01 /              3780   
CIB      DATA  ROOTI( 23) / Z414B0BF165515A9B /
      DATA  ROOTI( 24) /       0.479583152331271956 0D+01 /             cni
CIB      DATA  ROOTI( 24) / Z414CBBB9D5DC105A /
      DATA  ROOTI( 25) /       0.489897948556635621 0D+01 /             cni
CIB      DATA  ROOTI( 25) / Z414E6238502484BA /
      DATA  ROOTI( 26) /       0.500000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 26) / Z4150000000000000 /
      DATA  ROOTI( 27) /       0.509901951359278494 0D+01 /             cni
CIB      DATA  ROOTI( 27) / Z415195957C48BFDA /
      DATA  ROOTI( 28) /       0.519615242270663180 0D+01 /              3790   
CIB      DATA  ROOTI( 28) / Z41532370B908E5FF /
      DATA  ROOTI( 29) /       0.529150262212918121 0D+01 /             cni
CIB      DATA  ROOTI( 29) / Z4154A9FEA74BE3A7 /
      DATA  ROOTI( 30) /       0.538516480713450396 0D+01 /             cni
CIB      DATA  ROOTI( 30) / Z415629A292A367CD /
      DATA  ROOTI( 31) /       0.547722557505166119 0D+01 /             cni
CIB      DATA  ROOTI( 31) / Z4157A2B748DA963C /
      DATA  ROOTI( 32) /       0.556776436283002196 0D+01 /             cni
CIB      DATA  ROOTI( 32) / Z4159159015A3070E /
      DATA  ROOTI( 33) /       0.565685424949238014 0D+01 /              3800   
CIB      DATA  ROOTI( 33) / Z415A827999FCEF32 /
      DATA  ROOTI( 34) /       0.574456264653802862 0D+01 /             cni
CIB      DATA  ROOTI( 34) / Z415BE9BA858B43C0 /
      DATA  ROOTI( 35) /       0.583095189484530052 0D+01 /             cni
CIB      DATA  ROOTI( 35) / Z415D4B9436CE8E87 /
      DATA  ROOTI( 36) /       0.591607978309961613 0D+01 /             cni
CIB      DATA  ROOTI( 36) / Z415EA843464F08B4 /
      DATA  ROOTI( 37) /       0.600000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 37) / Z4160000000000000 /
      DATA  ROOTI( 38) /       0.608276253029821978 0D+01 /              3810   
CIB      DATA  ROOTI( 38) / Z416152FECD8F70E6 /
      DATA  ROOTI( 39) /       0.616441400296897646 0D+01 /             cni
CIB      DATA  ROOTI( 39) / Z4162A17093DC1966 /
      DATA  ROOTI( 40) /       0.624499799839839831 0D+01 /             cni
CIB      DATA  ROOTI( 40) / Z4163EB83056B4E28 /
      DATA  ROOTI( 41) /       0.632455532033675860 0D+01 /             cni
CIB      DATA  ROOTI( 41) / Z41653160EB696D4A /
      DATA  ROOTI( 42) /       0.640312423743284875 0D+01 /             cni
CIB      DATA  ROOTI( 42) / Z41667332667FFC01 /
      DATA  ROOTI( 43) /       0.648074069840786016 0D+01 /              3820   
CIB      DATA  ROOTI( 43) / Z4167B11D2898498F /
      DATA  ROOTI( 44) /       0.655743852430200058 0D+01 /             cni
CIB      DATA  ROOTI( 44) / Z4168EB44A8768581 /
      DATA  ROOTI( 45) /       0.663324958071079962 0D+01 /             cni
CIB      DATA  ROOTI( 45) / Z416A21CA4FF5BCD0 /
      DATA  ROOTI( 46) /       0.670820393249936919 0D+01 /             cni
CIB      DATA  ROOTI( 46) / Z416B54CDA58FBBEF /
      DATA  ROOTI( 47) /       0.678232998312526814 0D+01 /             cni
CIB      DATA  ROOTI( 47) / Z416C846C71C3408C /
      DATA  ROOTI( 48) /       0.685565460040104413 0D+01 /              3830   
CIB      DATA  ROOTI( 48) / Z416DB0C2E0D64F99 /
      DATA  ROOTI( 49) /       0.692820323027550922 0D+01 /             cni
CIB      DATA  ROOTI( 49) / Z416ED9EBA16132AA /
      DATA  ROOTI( 50) /       0.700000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 50) / Z4170000000000000 /
      DATA  ROOTI( 51) /       0.707106781186547528 0D+01 /             cni
CIB      DATA  ROOTI( 51) / Z41712318007C2AFF /
      DATA  ROOTI( 52) /       0.714142842854284998 0D+01 /             cni
CIB      DATA  ROOTI( 52) / Z4172434A74B50F36 /
      DATA  ROOTI( 53) /       0.721110255092797869 0D+01 /              3840   
CIB      DATA  ROOTI( 53) / Z417360AD118567CE /
      DATA  ROOTI( 54) /       0.728010988928051828 0D+01 /             cni
CIB      DATA  ROOTI( 54) / Z41747B5481DBEFA5 /
      DATA  ROOTI( 55) /       0.734846922834953431 0D+01 /             cni
CIB      DATA  ROOTI( 55) / Z417593547836C717 /
      DATA  ROOTI( 56) /       0.741619848709566298 0D+01 /             cni
CIB      DATA  ROOTI( 56) / Z4176A8BFBEAB8760 /
      DATA  ROOTI( 57) /       0.748331477354788288 0D+01 /             cni
CIB      DATA  ROOTI( 57) / Z4177BBA845A0FD25 /
      DATA  ROOTI( 58) /       0.754983443527074960 0D+01 /              3850   
CIB      DATA  ROOTI( 58) / Z4178CC1F315B3D6F /
      DATA  ROOTI( 59) /       0.761577310586390821 0D+01 /             cni
CIB      DATA  ROOTI( 59) / Z4179DA34E67711BE /
      DATA  ROOTI( 60) /       0.768114574786860826 0D+01 /             cni
CIB      DATA  ROOTI( 60) / Z417AE5F9156E7B6E /
      DATA  ROOTI( 61) /       0.774596669241483382 0D+01 /             cni
CIB      DATA  ROOTI( 61) / Z417BEF7AC53D3B67 /
      DATA  ROOTI( 62) /       0.781024967590665442 0D+01 /             cni
CIB      DATA  ROOTI( 62) / Z417CF6C85D39D1A2 /
      DATA  ROOTI( 63) /       0.787400787401181113 0D+01 /              3860   
CIB      DATA  ROOTI( 63) / Z417DFBEFAE353C48 /
      DATA  ROOTI( 64) /       0.793725393319377170 0D+01 /             cni
CIB      DATA  ROOTI( 64) / Z417EFEFDFAF1D57A /
      DATA  ROOTI( 65) /       0.800000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 65) / Z4180000000000000 /
      DATA  ROOTI( 66) /       0.806225774829854958 0D+01 /             cni
CIB      DATA  ROOTI( 66) / Z4180FF01FB0DD682 /
      DATA  ROOTI( 67) /       0.812403840463596039 0D+01 /             cni
CIB      DATA  ROOTI( 67) / Z4181FC0FB1B5C05E /
      DATA  ROOTI( 68) /       0.818535277187244992 0D+01 /              3870   
CIB      DATA  ROOTI( 68) / Z4182F73477D6A456 /
      DATA  ROOTI( 69) /       0.824621125123532117 0D+01 /             cni
CIB      DATA  ROOTI( 69) / Z4183F07B357F6838 /
      DATA  ROOTI( 70) /       0.830662386291807486 0D+01 /             cni
CIB      DATA  ROOTI( 70) / Z4184E7EE6C768048 /
      DATA  ROOTI( 71) /       0.836660026534075540 0D+01 /             cni
CIB      DATA  ROOTI( 71) / Z4185DD983D657ED7 /
      DATA  ROOTI( 72) /       0.842614977317635860 0D+01 /             cni
CIB      DATA  ROOTI( 72) / Z4186D1826CAFD82E /
      DATA  ROOTI( 73) /       0.848528137423857021 0D+01 /              3880   
CIB      DATA  ROOTI( 73) / Z4187C3B666FB66CB /
      DATA  ROOTI( 74) /       0.854400374531753126 0D+01 /             cni
CIB      DATA  ROOTI( 74) / Z4188B43D4570A51C /
      DATA  ROOTI( 75) /       0.860232526704262668 0D+01 /             cni
CIB      DATA  ROOTI( 75) / Z4189A31FD1B80A70 /
      DATA  ROOTI( 76) /       0.866025403784438641 0D+01 /             cni
CIB      DATA  ROOTI( 76) / Z418A906689B97F54 /
      DATA  ROOTI( 77) /       0.871779788708134706 0D+01 /             cni
CIB      DATA  ROOTI( 77) / Z418B7C19A3226FC4 /
      DATA  ROOTI( 78) /       0.877496438739212214 0D+01 /              3890   
CIB      DATA  ROOTI( 78) / Z418C66410EB69F16 /
      DATA  ROOTI( 79) /       0.883176086632784685 0D+01 /             cni
CIB      DATA  ROOTI( 79) / Z418D4EE47B6F881C /
      DATA  ROOTI( 80) /       0.888819441731558890 0D+01 /             cni
CIB      DATA  ROOTI( 80) / Z418E360B596DC381 /
      DATA  ROOTI( 81) /       0.894427190999915878 0D+01 /             cni
CIB      DATA  ROOTI( 81) / Z418F1BBcniBFA53E /
      DATA  ROOTI( 82) /       0.900000000000000000 0D+01 /             cni
CIB      DATA  ROOTI( 82) / Z4190000000000000 /
      DATA  ROOTI( 83) /       0.905538513813741663 0D+01 /              3900   
CIB      DATA  ROOTI( 83) / Z4190E2DB86CFC11D /
      DATA  ROOTI( 84) /       0.911043357914429897 0D+01 /             cni
CIB      DATA  ROOTI( 84) / Z4191C456002CE13F /
      DATA  ROOTI( 85) /       0.916515138991168010 0D+01 /             cni
CIB      DATA  ROOTI( 85) / Z4192A475C8A8F42A /
      DATA  ROOTI( 86) /       0.921954445729288730 0D+01 /             cni
CIB      DATA  ROOTI( 86) / Z419383410C8174D1 /
      DATA  ROOTI( 87) /       0.927361849549570372 0D+01 /             cni
CIB      DATA  ROOTI( 87) / Z419460BDC99BC19F /
      DATA  ROOTI( 88) /       0.932737905308881499 0D+01 /              3910   
CIB      DATA  ROOTI( 88) / Z41953CF1D166972D /
      DATA  ROOTI( 89) /       0.938083151964685902 0D+01 /             cni
CIB      DATA  ROOTI( 89) / Z419617E2CAA2B536 /
      DATA  ROOTI( 90) /       0.943398113205660382 0D+01 /             cni
CIB      DATA  ROOTI( 90) / Z4196F19633143A0B /
      DATA  ROOTI( 91) /       0.948683298050513790 0D+01 /             cni
CIB      DATA  ROOTI( 91) / Z4197CA11611E23EF /
      DATA  ROOTI( 92) /       0.953939201416945659 0D+01 /             cni
CIB      DATA  ROOTI( 92) / Z4198A15985494D5A /
      DATA  ROOTI( 93) /       0.959166304662543912 0D+01 /              3920   
CIB      DATA  ROOTI( 93) / Z41997773ABB820B4 /
      DATA  ROOTI( 94) /       0.964365076099295493 0D+01 /             cni
CIB      DATA  ROOTI( 94) / Z419A4C64BD882A00 /
      DATA  ROOTI( 95) /       0.969535971483265802 0D+01 /             cni
CIB      DATA  ROOTI( 95) / Z419B20318222982D /
      DATA  ROOTI( 96) /       0.974679434480896401 0D+01 /             cni
CIB      DATA  ROOTI( 96) / Z419BF2DEA07CAD0C /
      DATA  ROOTI( 97) /       0.979795897113271241 0D+01 /             cni
CIB      DATA  ROOTI( 97) / Z419CC470A0490974 /
      DATA  ROOTI( 98) /       0.984885780179610482 0D+01 /              3930   
CIB      DATA  ROOTI( 98) / Z419D94EBEB1AB314 /
      DATA  ROOTI( 99) /       0.989949493661166535 0D+01 /             cni
CIB      DATA  ROOTI( 99) / Z419E6454CD7AA298 /
      DATA  ROOTI(100) /       0.994987437106619965 0D+01 /             cni
CIB      DATA  ROOTI(100) / Z419F32AF77F09B39 /
      DATA  ROOTI(101) /       0.100000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(101) / Z41A0000000000000 /
      DATA  ROOTI(102) /       0.100498756211208902 0D+02 /             cni
CIB      DATA  ROOTI(102) / Z41A0CC4A61194F81 /
      DATA  ROOTI(103) /       0.100995049383620779 0D+02 /              3940   
CIB      DATA  ROOTI(103) / Z41A197927D80E3D2 /
      DATA  ROOTI(104) /       0.101488915650922196 0D+02 /             cni
CIB      DATA  ROOTI(104) / Z41A261DC1F2B8A9A /
      DATA  ROOTI(105) /       0.101980390271855696 0D+02 /             cni
CIB      DATA  ROOTI(105) / Z41A32B2AF8917FB3 /
      DATA  ROOTI(106) /       0.102469507659595984 0D+02 /             cni
CIB      DATA  ROOTI(106) / Z41A3F382A5784C4A /
      DATA  ROOTI(107) /       0.102956301409870004 0D+02 /             cni
CIB      DATA  ROOTI(107) / Z41A4BAE6ABB4043E /
      DATA  ROOTI(108) /       0.103440804327886005 0D+02 /              3950   
CIB      DATA  ROOTI(108) / Z41A5815A7BE0543C /
      DATA  ROOTI(109) /       0.103923048454132638 0D+02 /             cni
CIB      DATA  ROOTI(109) / Z41A646E17211CBFF /
      DATA  ROOTI(110) /       0.104403065089105502 0D+02 /             cni
CIB      DATA  ROOTI(110) / Z41A70B7ED67FC9B6 /
      DATA  ROOTI(111) /       0.104880884817015154 0D+02 /             cni
CIB      DATA  ROOTI(111) / Z41A7CF35DE276598 /
      DATA  ROOTI(112) /       0.105356537528527388 0D+02 /             cni
CIB      DATA  ROOTI(112) / Z41A89209AB67B702 /
      DATA  ROOTI(113) /       0.105830052442583624 0D+02 /              3960   
CIB      DATA  ROOTI(113) / Z41A953FD4E97C74E /
      DATA  ROOTI(114) /       0.106301458127346493 0D+02 /             cni
CIB      DATA  ROOTI(114) / Z41AA1513C69681AD /
      DATA  ROOTI(115) /       0.106770782520313112 0D+02 /             cni
CIB      DATA  ROOTI(115) / Z41AAD5500154EAD0 /
      DATA  ROOTI(116) /       0.107238052947636082 0D+02 /             cni
CIB      DATA  ROOTI(116) / Z41AB94B4DC5AE6C2 /
      DATA  ROOTI(117) /       0.107703296142690081 0D+02 /             cni
CIB      DATA  ROOTI(117) / Z41AC53452546CF9B /
      DATA  ROOTI(118) /       0.108166538263919678 0D+02 /              3970   
CIB      DATA  ROOTI(118) / Z41AD11039A481BB4 /
      DATA  ROOTI(119) /       0.108627804912002157 0D+02 /             cni
CIB      DATA  ROOTI(119) / Z41ADCDF2EA954ED1 /
      DATA  ROOTI(120) /       0.109087121146357144 0D+02 /             cni
CIB      DATA  ROOTI(120) / Z41AE8A15B6DD6E2B /
      DATA  ROOTI(121) /       0.109544511501033224 0D+02 /             cni
CIB      DATA  ROOTI(121) / Z41AF456E91B52C78 /
      DATA  ROOTI(122) /       0.110000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(122) / Z41B0000000000000 /
      DATA  ROOTI(123) /       0.110453610171872607 0D+02 /              3980   
CIB      DATA  ROOTI(123) / Z41B0B9CC7955523E /
      DATA  ROOTI(124) /       0.110905365064094172 0D+02 /             cni
CIB      DATA  ROOTI(124) / Z41B172D66861F5EE /
      DATA  ROOTI(125) /       0.111355287256600439 0D+02 /             cni
CIB      DATA  ROOTI(125) / Z41B22B202B460E1C /
      DATA  ROOTI(126) /       0.111803398874989486 0D+02 /             cni
CIB      DATA  ROOTI(126) / Z41B2E2AC13EF8E8E /
      DATA  ROOTI(127) /       0.112249721603218242 0D+02 /             cni
CIB      DATA  ROOTI(127) / Z41B3997C68717BB7 /
      DATA  ROOTI(128) /       0.112694276695846449 0D+02 /              3990   
CIB      DATA  ROOTI(128) / Z41B44F9363580E84 /
      DATA  ROOTI(129) /       0.113137084989847605 0D+02 /             cni
CIB      DATA  ROOTI(129) / Z41B504F333F9DE65 /
      DATA  ROOTI(130) /       0.113578166916005472 0D+02 /             cni
CIB      DATA  ROOTI(130) / Z41B5B99DFEC63240 /
      DATA  ROOTI(131) /       0.114017542509913798 0D+02 /             cni
CIB      DATA  ROOTI(131) / Z41B66D95DD90975B /
      DATA  ROOTI(132) /       0.114455231422595971 0D+02 /             cni
CIB      DATA  ROOTI(132) / Z41B720DCDFD9DBA6 /
      DATA  ROOTI(133) /       0.114891252930760572 0D+02 /              4000   
CIB      DATA  ROOTI(133) / Z41B7D3750B168780 /
      DATA  ROOTI(134) /       0.115325625946707959 0D+02 /             cni
CIB      DATA  ROOTI(134) / Z41B885605AF2F18D /
      DATA  ROOTI(135) /       0.115758369027902255 0D+02 /             cni
CIB      DATA  ROOTI(135) / Z41B936A0C19505F0 /
      DATA  ROOTI(136) /       0.116189500386222506 0D+02 /             cni
CIB      DATA  ROOTI(136) / Z41B9E73827DBD91A /
      DATA  ROOTI(137) /       0.116619037896906010 0D+02 /             cni
CIB      DATA  ROOTI(137) / Z41BA97286D9D1D0E /
      DATA  ROOTI(138) /       0.117046999107196250 0D+02 /              4010   
CIB      DATA  ROOTI(138) / Z41BB467369E08EFD /
      DATA  ROOTI(139) /       0.117473401244707305 0D+02 /             cni
CIB      DATA  ROOTI(139) / Z41BBF51AEB19721A /
      DATA  ROOTI(140) /       0.117898261225515959 0D+02 /             cni
CIB      DATA  ROOTI(140) / Z41BCA320B75E2B63 /
      DATA  ROOTI(141) /       0.118321595661992320 0D+02 /             cni
CIB      DATA  ROOTI(141) / Z41BD50868C9E1167 /
      DATA  ROOTI(142) /       0.118743420870379173 0D+02 /             cni
CIB      DATA  ROOTI(142) / Z41BDFD4E20D58202 /
      DATA  ROOTI(143) /       0.119163752878129849 0D+02 /              4020   
CIB      DATA  ROOTI(143) / Z41BEA97922404F4A /
      DATA  ROOTI(144) /       0.119582607431013981 0D+02 /             cni
CIB      DATA  ROOTI(144) / Z41BF5509378A941F /
      DATA  ROOTI(145) /       0.120000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(145) / Z41C0000000000000 /
      DATA  ROOTI(146) /       0.120415945787922956 0D+02 /             cni
CIB      DATA  ROOTI(146) / Z41C0AA5F13B9A92F /
      DATA  ROOTI(147) /       0.120830459735945721 0D+02 /             cni
CIB      DATA  ROOTI(147) / Z41C1542803CA735F /
      DATA  ROOTI(148) /       0.121243556529821410 0D+02 /              4030   
CIB      DATA  ROOTI(148) / Z41C1FD5C5A6A18A9 /
      DATA  ROOTI(149) /       0.121655250605964393 0D+02 /             cni
CIB      DATA  ROOTI(149) / Z41C2A5FD9B1EE1CB /
      DATA  ROOTI(150) /       0.122065556157337030 0D+02 /             cni
CIB      DATA  ROOTI(150) / Z41C34E0D42E61A34 /
      DATA  ROOTI(151) /       0.122474487139158905 0D+02 /             cni
CIB      DATA  ROOTI(151) / Z41C3F58CC85B4BD1 /
      DATA  ROOTI(152) /       0.122882057274445076 0D+02 /             cni
CIB      DATA  ROOTI(152) / Z41C49C7D9BDE4E07 /
      DATA  ROOTI(153) /       0.123288280059379529 0D+02 /              4040   
CIB      DATA  ROOTI(153) / Z41C542E127B832CC /
      DATA  ROOTI(154) /       0.123693168768529818 0D+02 /             cni
CIB      DATA  ROOTI(154) / Z41C5E8B8D03F1C54 /
      DATA  ROOTI(155) /       0.124096736459908565 0D+02 /             cni
CIB      DATA  ROOTI(155) / Z41C68E05F3F9055E /
      DATA  ROOTI(156) /       0.124498995979887324 0D+02 /             cni
CIB      DATA  ROOTI(156) / Z41C732C9EBBD85BF /
      DATA  ROOTI(157) /       0.124899959967967964 0D+02 /             cni
CIB      DATA  ROOTI(157) / Z41C7D7060AD69C4F /
      DATA  ROOTI(158) /       0.125299640861416677 0D+02 /              4050   
CIB      DATA  ROOTI(158) / Z41C87ABB9F208720 /
      DATA  ROOTI(159) /       0.125698050899765348 0D+02 /             cni
CIB      DATA  ROOTI(159) / Z41C91DEBF128B266 /
      DATA  ROOTI(160) /       0.126095202129184916 0D+02 /             cni
CIB      DATA  ROOTI(160) / Z41C9C098444BC628 /
      DATA  ROOTI(161) /       0.126491106406735174 0D+02 /             cni
CIB      DATA  ROOTI(161) / Z41CA62C1D6D2DA95 /
      DATA  ROOTI(162) /       0.126885775404495205 0D+02 /             cni
CIB      DATA  ROOTI(162) / Z41CB0469E20FDA59 /
      DATA  ROOTI(163) /       0.127279220613578554 0D+02 /              4060   
CIB      DATA  ROOTI(163) / Z41CBA5919A791A31 /
      DATA  ROOTI(164) /       0.127671453348037047 0D+02 /             cni
CIB      DATA  ROOTI(164) / Z41CC463A2FC42C93 /
      DATA  ROOTI(165) /       0.128062484748656973 0D+02 /             cni
CIB      DATA  ROOTI(165) / Z41CCE664CCFFF801 /
      DATA  ROOTI(166) /       0.128452325786651291 0D+02 /             cni
CIB      DATA  ROOTI(166) / Z41CD861298AE166F /
      DATA  ROOTI(167) /       0.128840987267251255 0D+02 /             cni
CIB      DATA  ROOTI(167) / Z41CE2544B4DB83B5 /
      DATA  ROOTI(168) /       0.129228479833200856 0D+02 /              4070   
CIB      DATA  ROOTI(168) / Z41CEC3FC3F38A10F /
      DATA  ROOTI(169) /       0.129614813968157205 0D+02 /             cni
CIB      DATA  ROOTI(169) / Z41CF623A5130931F /
      DATA  ROOTI(170) /       0.130000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(170) / Z41D0000000000000 /
      DATA  ROOTI(171) /       0.130384048104052974 0D+02 /             cni
CIB      DATA  ROOTI(171) / Z41D09D4E5CCB3284 /
      DATA  ROOTI(172) /       0.130766968306220206 0D+02 /             cni
CIB      DATA  ROOTI(172) / Z41D13A2674B3A7A6 /
      DATA  ROOTI(173) /       0.131148770486040014 0D+02 /              4080   
CIB      DATA  ROOTI(173) / Z41D1D68950ED0B03 /
      DATA  ROOTI(174) /       0.131529464379659053 0D+02 /             cni
CIB      DATA  ROOTI(174) / Z41D27277F6D1A6F0 /
      DATA  ROOTI(175) /       0.131909059582729191 0D+02 /             cni
CIB      DATA  ROOTI(175) / Z41D30DF367F64CB6 /
      DATA  ROOTI(176) /       0.132287565553229529 0D+02 /             cni
CIB      DATA  ROOTI(176) / Z41D3A8FCA23DB921 /
      DATA  ROOTI(177) /       0.132664991614215995 0D+02 /             cni
CIB      DATA  ROOTI(177) / Z41D443949FEB79A1 /
      DATA  ROOTI(178) /       0.133041346956500708 0D+02 /              4090   
CIB      DATA  ROOTI(178) / Z41D4DDBC57B655E2 /
      DATA  ROOTI(179) /       0.133416640641263338 0D+02 /             cni
CIB      DATA  ROOTI(179) / Z41D57774BCDA41BE /
      DATA  ROOTI(180) /       0.133790881602596521 0D+02 /             cni
CIB      DATA  ROOTI(180) / Z41D610BEBF29DB30 /
      DATA  ROOTI(181) /       0.134164078649987382 0D+02 /             cni
CIB      DATA  ROOTI(181) / Z41D6A99B4B1F77DD /
      DATA  ROOTI(182) /       0.134536240470737103 0D+02 /             cni
CIB      DATA  ROOTI(182) / Z41D7420B49EDC5A2 /
      DATA  ROOTI(183) /       0.134907375632320414 0D+02 /              4100   
CIB      DATA  ROOTI(183) / Z41D7DA0FA190016F /
      DATA  ROOTI(184) /       0.135277492584686829 0D+02 /             cni
CIB      DATA  ROOTI(184) / Z41D871A934D9C7A8 /
      DATA  ROOTI(185) /       0.135646599662505363 0D+02 /             cni
CIB      DATA  ROOTI(185) / Z41D908D8E3868118 /
      DATA  ROOTI(186) /       0.136014705087354433 0D+02 /             cni
CIB      DATA  ROOTI(186) / Z41D99F9F8A486F75 /
      DATA  ROOTI(187) /       0.136381816969858558 0D+02 /             cni
CIB      DATA  ROOTI(187) / Z41DA35FE02D75C4B /
      DATA  ROOTI(188) /       0.136747943311773432 0D+02 /              4110   
CIB      DATA  ROOTI(188) / Z41DACBF523FEED16 /
      DATA  ROOTI(189) /       0.137113092008020883 0D+02 /             cni
CIB      DATA  ROOTI(189) / Z41DB6185C1AC9F32 /
      DATA  ROOTI(190) /       0.137477270848675199 0D+02 /             cni
CIB      DATA  ROOTI(190) / Z41DBF6B0ACFD6E3E /
      DATA  ROOTI(191) /       0.137840487520902217 0D+02 /             cni
CIB      DATA  ROOTI(191) / Z41DC8B76B44B2761 /
      DATA  ROOTI(192) /       0.138202749610852533 0D+02 /             cni
CIB      DATA  ROOTI(192) / Z41DD1FD8A3396BDF /
      DATA  ROOTI(193) /       0.138564064605510184 0D+02 /              4120   
CIB      DATA  ROOTI(193) / Z41DDB3D742C26554 /
      DATA  ROOTI(194) /       0.138924439894498044 0D+02 /             cni
CIB      DATA  ROOTI(194) / Z41DE477359432DCA /
      DATA  ROOTI(195) /       0.139283882771841194 0D+02 /             cni
CIB      DATA  ROOTI(195) / Z41DEDAADAA87EDE1 /
      DATA  ROOTI(196) /       0.139642400437689411 0D+02 /             cni
CIB      DATA  ROOTI(196) / Z41DF6D86F7D7B30A /
      DATA  ROOTI(197) /       0.140000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(197) / Z41E0000000000000 /
      DATA  ROOTI(198) /       0.140356688476181997 0D+02 /              4130   
CIB      DATA  ROOTI(198) / Z41E092197F60194B /
      DATA  ROOTI(199) /       0.140712472794702887 0D+02 /             cni
CIB      DATA  ROOTI(199) / Z41E123D42FF40FD2 /
      DATA  ROOTI(200) /       0.141067359796658844 0D+02 /             cni
CIB      DATA  ROOTI(200) / Z41E1B530C95F8B3E /
      DATA  ROOTI(201) /       0.141421356237309506 0D+02 /             cni
CIB      DATA  ROOTI(201) / Z41E2463000F855FE /
      DATA  ROOTI(202) /       0.141774468787578252 0D+02 /             cni
CIB      DATA  ROOTI(202) / Z41E2D6D289D0AC97 /
      DATA  ROOTI(203) /       0.142126704035518956 0D+02 /              4140   
CIB      DATA  ROOTI(203) / Z41E3671914C151FA /
      DATA  ROOTI(204) /       0.142478068487750071 0D+02 /             cni
CIB      DATA  ROOTI(204) / Z41E3F70450736A63 /
      DATA  ROOTI(205) /       0.142828568570857000 0D+02 /             cni
CIB      DATA  ROOTI(205) / Z41E48694E96A1E6C /
      DATA  ROOTI(206) /       0.143178210632763532 0D+02 /             cni
CIB      DATA  ROOTI(206) / Z41E515CB8A0C07B7 /
      DATA  ROOTI(207) /       0.143527000944073237 0D+02 /             cni
CIB      DATA  ROOTI(207) / Z41E5A4A8DAAC68BA /
      DATA  ROOTI(208) /       0.143874945699381587 0D+02 /              4150   
CIB      DATA  ROOTI(208) / Z41E6332D8194310E /
      DATA  ROOTI(209) /       0.144222051018559572 0D+02 /             cni
CIB      DATA  ROOTI(209) / Z41E6C15A230ACF9B /
      DATA  ROOTI(210) /       0.144568322948009602 0D+02 /             cni
CIB      DATA  ROOTI(210) / Z41E74F2F615ED3FD /
      DATA  ROOTI(211) /       0.144913767461894385 0D+02 /             cni
CIB      DATA  ROOTI(211) / Z41E7DCADDCEE6062 /
      DATA  ROOTI(212) /       0.145258390463339502 0D+02 /             cni
CIB      DATA  ROOTI(212) / Z41E869D6342F6D23 /
      DATA  ROOTI(213) /       0.145602197785610366 0D+02 /              4160   
CIB      DATA  ROOTI(213) / Z41E8F6A903B7DF4A /
      DATA  ROOTI(214) /       0.145945195193264241 0D+02 /             cni
CIB      DATA  ROOTI(214) / Z41E98326E645733A /
      DATA  ROOTI(215) /       0.146287388383277934 0D+02 /             cni
CIB      DATA  ROOTI(215) / Z41EA0F5074C57C89 /
      DATA  ROOTI(216) /       0.146628782986151802 0D+02 /             cni
CIB      DATA  ROOTI(216) / Z41EA9B26465C7C32 /
      DATA  ROOTI(217) /       0.146969384566990686 0D+02 /             cni
CIB      DATA  ROOTI(217) / Z41EB26A8F06D8E2E /
      DATA  ROOTI(218) /       0.147309198626562354 0D+02 /              4170   
CIB      DATA  ROOTI(218) / Z41EBB1D906A1AF7C /
      DATA  ROOTI(219) /       0.147648230602334005 0D+02 /             cni
CIB      DATA  ROOTI(219) / Z41EC3CB71AEEDD91 /
      DATA  ROOTI(220) /       0.147986485869487421 0D+02 /             cni
CIB      DATA  ROOTI(220) / Z41ECC743BD9F1038 /
      DATA  ROOTI(221) /       0.148323969741913260 0D+02 /             cni
CIB      DATA  ROOTI(221) / Z41ED517F7D570EC0 /
      DATA  ROOTI(222) /       0.148660687473185056 0D+02 /             cni
CIB      DATA  ROOTI(222) / Z41EDDB6AE71D2176 /
      DATA  ROOTI(223) /       0.148996644257513398 0D+02 /              4180   
CIB      DATA  ROOTI(223) / Z41EE6506865FA041 /
      DATA  ROOTI(224) /       0.149331845230680786 0D+02 /             cni
CIB      DATA  ROOTI(224) / Z41EEEE52E4FB5F41 /
      DATA  ROOTI(225) /       0.149666295470957655 0D+02 /             cni
CIB      DATA  ROOTI(225) / Z41EF77508B41FA49 /
      DATA  ROOTI(226) /       0.150000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(226) / Z41F0000000000000 /
      DATA  ROOTI(227) /       0.150332963783729083 0D+02 /             cni
CIB      DATA  ROOTI(227) / Z41F08861C882FD79 /
      DATA  ROOTI(228) /       0.150665191733193635 0D+02 /              4190   
CIB      DATA  ROOTI(228) / Z41F11076689F6AFF /
      DATA  ROOTI(229) /       0.150996688705414994 0D+02 /             cni
CIB      DATA  ROOTI(229) / Z41F1983E62B67ADF /
      DATA  ROOTI(230) /       0.151327459504215560 0D+02 /             cni
CIB      DATA  ROOTI(230) / Z41F21FBA37BBCAD6 /
      DATA  ROOTI(231) /       0.151657508881031011 0D+02 /             cni
CIB      DATA  ROOTI(231) / Z41F2A6EA673AF8EF /
      DATA  ROOTI(232) /       0.151986841535706636 0D+02 /             cni
CIB      DATA  ROOTI(232) / Z41F32DCF6F5D1C6F /
      DATA  ROOTI(233) /       0.152315462117278166 0D+02 /              4200   
CIB      DATA  ROOTI(233) / Z41F3B469CCEE237D /
      DATA  ROOTI(234) /       0.152643375224737481 0D+02 /             cni
CIB      DATA  ROOTI(234) / Z41F43AB9FB62162C /
      DATA  ROOTI(235) /       0.152970585407783546 0D+02 /             cni
CIB      DATA  ROOTI(235) / Z41F4C0C074DA3F8D /
      DATA  ROOTI(236) /       0.153297097167558916 0D+02 /             cni
CIB      DATA  ROOTI(236) / Z41F5467DB22A3D59 /
      DATA  ROOTI(237) /       0.153622914957372163 0D+02 /             cni
CIB      DATA  ROOTI(237) / Z41F5CBF22ADCF6DB /
      DATA  ROOTI(238) /       0.153948043183406524 0D+02 /              4210   
CIB      DATA  ROOTI(238) / Z41F6511E55397B99 /
      DATA  ROOTI(239) /       0.154272486205415125 0D+02 /             cni
CIB      DATA  ROOTI(239) / Z41F6D602A647CA62 /
      DATA  ROOTI(240) /       0.154596248337403066 0D+02 /             cni
CIB      DATA  ROOTI(240) / Z41F75A9F91D5813F /
      DATA  ROOTI(241) /       0.154919333848296676 0D+02 /             cni
CIB      DATA  ROOTI(241) / Z41F7DEF58A7A76CE /
      DATA  ROOTI(242) /       0.155241746962600238 0D+02 /             cni
CIB      DATA  ROOTI(242) / Z41F86305019D3D96 /
      DATA  ROOTI(243) /       0.155563491861040455 0D+02 /              4220   
CIB      DATA  ROOTI(243) / Z41F8E6CE677791CA /
      DATA  ROOTI(244) /       0.155884572681198956 0D+02 /             cni
CIB      DATA  ROOTI(244) / Z41F96A522B1AB1FE /
      DATA  ROOTI(245) /       0.156204993518133088 0D+02 /             cni
CIB      DATA  ROOTI(245) / Z41F9ED90BA73A344 /
      DATA  ROOTI(246) /       0.156524758424985280 0D+02 /             cni
CIB      DATA  ROOTI(246) / Z41FA708A824F612D /
      DATA  ROOTI(247) /       0.156843871413581220 0D+02 /             cni
CIB      DATA  ROOTI(247) / Z41FAF33FEE5EFA1D /
      DATA  ROOTI(248) /       0.157162336455017109 0D+02 /              4230   
CIB      DATA  ROOTI(248) / Z41FB75B1693B9865 /
      DATA  ROOTI(249) /       0.157480157480236220 0D+02 /             cni
CIB      DATA  ROOTI(249) / Z41FBF7DF5C6A788F /
      DATA  ROOTI(250) /       0.157797338380595000 0D+02 /             cni
CIB      DATA  ROOTI(250) / Z41FC79CA3060CD43 /
      DATA  ROOTI(251) /       0.158113883008418967 0D+02 /             cni
CIB      DATA  ROOTI(251) / Z41FCFB724C87913A /
      DATA  ROOTI(252) /       0.158429795177548596 0D+02 /             cni
CIB      DATA  ROOTI(252) / Z41FD7CD8173F4792 /
      DATA  ROOTI(253) /       0.158745078663875436 0D+02 /              4240   
CIB      DATA  ROOTI(253) / Z41FDFDFBF5E3AAF5 /
      DATA  ROOTI(254) /       0.159059737205868663 0D+02 /             cni
CIB      DATA  ROOTI(254) / Z41FE7EDE4CCF4BE9 /
      DATA  ROOTI(255) /       0.159373774505092274 0D+02 /             cni
CIB      DATA  ROOTI(255) / Z41FEFF7F7F5F1EAE /
      DATA  ROOTI(256) /       0.159687194226713121 0D+02 /             cni
CIB      DATA  ROOTI(256) / Z41FF7FDFEFF5F8FB /
      DATA  ROOTI(257) /       0.160000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(257) / Z4210000000000000 /
      DATA  ROOTI(258) /       0.160312195418813985 0D+02 /              4250   
CIB      DATA  ROOTI(258) / Z421007FE00FF6070 /
      DATA  ROOTI(259) /       0.160623784042090101 0D+02 /             cni
CIB      DATA  ROOTI(259) / Z42100FF807F60DEB /
      DATA  ROOTI(260) /       0.160934769394310813 0D+02 /             cni
CIB      DATA  ROOTI(260) / Z421017EE1Acni963 /
      DATA  ROOTI(261) /       0.161245154965970983 0D+02 /             cni
CIB      DATA  ROOTI(261) / Z42101FE03F61BAD0 /
      DATA  ROOTI(262) /       0.161554944214035103 0D+02 /             cni
CIB      DATA  ROOTI(262) / Z421027CE7B7EA376 /
      DATA  ROOTI(263) /       0.161864140562386467 0D+02 /              4260   
CIB      DATA  ROOTI(263) / Z42102FB8D4E30F48 /
      DATA  ROOTI(264) /       0.162172747402268556 0D+02 /             cni
CIB      DATA  ROOTI(264) / Z4210379F513F8570 /
      DATA  ROOTI(265) /       0.162480768092719217 0D+02 /             cni
CIB      DATA  ROOTI(265) / Z42103F81F636B80C /
      DATA  ROOTI(266) /       0.162788205960997061 0D+02 /             cni
CIB      DATA  ROOTI(266) / Z42104760C95DB310 /
      DATA  ROOTI(267) /       0.163095064303000896 0D+02 /             cni
CIB      DATA  ROOTI(267) / Z42104F3BD03C0A64 /
      DATA  ROOTI(268) /       0.163401346383681911 0D+02 /              4270   
CIB      DATA  ROOTI(268) / Z42105713104C0736 /
      DATA  ROOTI(269) /       0.163707055437449007 0D+02 /             cni
CIB      DATA  ROOTI(269) / Z42105EE68EFAD48B /
      DATA  ROOTI(270) /       0.164012194668567268 0D+02 /             cni
CIB      DATA  ROOTI(270) / Z421066B651A8AB0F /
      DATA  ROOTI(271) /       0.164316767251549827 0D+02 /             cni
CIB      DATA  ROOTI(271) / Z42106E825DA8FC2B /
      DATA  ROOTI(272) /       0.164620776331543297 0D+02 /             cni
CIB      DATA  ROOTI(272) / Z4210764AB8429C66 /
      DATA  ROOTI(273) /       0.164924225024706423 0D+02 /              4280   
CIB      DATA  ROOTI(273) / Z42107E0F66AFED07 /
      DATA  ROOTI(274) /       0.165227116418583044 0D+02 /             cni
CIB      DATA  ROOTI(274) / Z421085D06E1F0517 /
      DATA  ROOTI(275) /       0.165529453572468483 0D+02 /             cni
CIB      DATA  ROOTI(275) / Z42108D8DD3B1D9AA /
      DATA  ROOTI(276) /       0.165831239517770008 0D+02 /             cni
CIB      DATA  ROOTI(276) / Z421095479C7E6581 /
      DATA  ROOTI(277) /       0.166132477258361497 0D+02 /             cni
CIB      DATA  ROOTI(277) / Z42109CFDCD8ED009 /
      DATA  ROOTI(278) /       0.166433169770932388 0D+02 /              4290   
CIB      DATA  ROOTI(278) / Z4210A4B06BE193B9 /
      DATA  ROOTI(279) /       0.166733320005330654 0D+02 /             cni
CIB      DATA  ROOTI(279) / Z4210AC5F7C69A3C8 /
      DATA  ROOTI(280) /       0.167032930884900672 0D+02 /             cni
CIB      DATA  ROOTI(280) / Z4210B40B040E9153 /
      DATA  ROOTI(281) /       0.167332005306815113 0D+02 /             cni
CIB      DATA  ROOTI(281) / Z4210BBB307ACAFDB /
      DATA  ROOTI(282) /       0.167630546142402110 0D+02 /             cni
CIB      DATA  ROOTI(282) / Z4210C3578C15393E /
      DATA  ROOTI(283) /       0.167928556237466644 0D+02 /              4300   
CIB      DATA  ROOTI(283) / Z4210CAF8960E710D /
      DATA  ROOTI(284) /       0.168226038412607224 0D+02 /             cni
CIB      DATA  ROOTI(284) / Z4210D2962A53C75E /
      DATA  ROOTI(285) /       0.168522995463527181 0D+02 /             cni
CIB      DATA  ROOTI(285) / Z4210DA304D95FB06 /
      DATA  ROOTI(286) /       0.168819430161341337 0D+02 /             cni
CIB      DATA  ROOTI(286) / Z4210E1C7047B3B51 /
      DATA  ROOTI(287) /       0.169115345252877631 0D+02 /             cni
CIB      DATA  ROOTI(287) / Z4210E95A539F492C /
      DATA  ROOTI(288) /       0.169410743460974160 0D+02 /              4310   
CIB      DATA  ROOTI(288) / Z4210F0EA3F9397CE /
      DATA  ROOTI(289) /       0.169705627484771391 0D+02 /             cni
CIB      DATA  ROOTI(289) / Z4210F876CCDF6CD9 /
      DATA  ROOTI(290) /       0.170000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(290) / Z4211000000000000 /
      DATA  ROOTI(291) /       0.170293863659264026 0D+02 /             cni
CIB      DATA  ROOTI(291) / Z42110785DD689A29 /
      DATA  ROOTI(292) /       0.170587221092319794 0D+02 /             cni
CIB      DATA  ROOTI(292) / Z42110F086982B418 /
      DATA  ROOTI(293) /       0.170880074906350607 0D+02 /              4320   
CIB      DATA  ROOTI(293) / Z42111687A8AE14A3 /
      DATA  ROOTI(294) /       0.171172427686236901 0D+02 /             cni
CIB      DATA  ROOTI(294) / Z42111E039F40EE66 /
      DATA  ROOTI(295) /       0.171464281994822478 0D+02 /             cni
CIB      DATA  ROOTI(295) / Z4211257C5187FD09 /
      DATA  ROOTI(296) /       0.171755640373176668 0D+02 /             cni
CIB      DATA  ROOTI(296) / Z42112CF1C3C6A213 /
      DATA  ROOTI(297) /       0.172046505340852534 0D+02 /             cni
CIB      DATA  ROOTI(297) / Z42113463FA37014E /
      DATA  ROOTI(298) /       0.172336879396140858 0D+02 /              4330   
CIB      DATA  ROOTI(298) / Z42113BD2F90A1CB4 /
      DATA  ROOTI(299) /       0.172626765016320682 0D+02 /             cni
CIB      DATA  ROOTI(299) / Z4211433EC467EFFB /
      DATA  ROOTI(300) /       0.172916164657905824 0D+02 /             cni
CIB      DATA  ROOTI(300) / Z42114AA7606F8BB0 /
      DATA  ROOTI(301) /       0.173205080756887746 0D+02 /             cni
CIB      DATA  ROOTI(301) / Z4211520CD1372FEB /
      DATA  ROOTI(302) /       0.173493515728974721 0D+02 /             cni
CIB      DATA  ROOTI(302) / Z4211596F1ACC669B /
      DATA  ROOTI(303) /       0.173781471969827663 0D+02 /              4340   
CIB      DATA  ROOTI(303) / Z421160CE41341D74 /
      DATA  ROOTI(304) /       0.174068951855292120 0D+02 /             cni
CIB      DATA  ROOTI(304) / Z4211682A486ABF71 /
      DATA  ROOTI(305) /       0.174355957741626959 0D+02 /             cni
CIB      DATA  ROOTI(305) / Z42116F8334644DF9 /
      DATA  ROOTI(306) /       0.174642491965729789 0D+02 /             cni
CIB      DATA  ROOTI(306) / Z421176D9090C79A8 /
      DATA  ROOTI(307) /       0.174928556845359005 0D+02 /             cni
CIB      DATA  ROOTI(307) / Z42117E2BCA46BAB9 /
      DATA  ROOTI(308) /       0.175214154679352312 0D+02 /              4350   
CIB      DATA  ROOTI(308) / Z4211857B7BEE690D /
      DATA  ROOTI(309) /       0.175499287747842452 0D+02 /             cni
CIB      DATA  ROOTI(309) / Z42118CC821D6D3E3 /
      DATA  ROOTI(310) /       0.175783958312469473 0D+02 /             cni
CIB      DATA  ROOTI(310) / Z42119411BFCB592F /
      DATA  ROOTI(311) /       0.176068168616590093 0D+02 /             cni
CIB      DATA  ROOTI(311) / Z42119B58598F7C9F /
      DATA  ROOTI(312) /       0.176351920885483970 0D+02 /             cni
CIB      DATA  ROOTI(312) / Z4211A29BF2DEFE49 /
      DATA  ROOTI(313) /       0.176635217326556955 0D+02 /              4360   
CIB      DATA  ROOTI(313) / Z4211A9DC8F6DF104 /
      DATA  ROOTI(314) /       0.176918060129541317 0D+02 /             cni
CIB      DATA  ROOTI(314) / Z4211B11A32E8D06C /
      DATA  ROOTI(315) /       0.177200451466693494 0D+02 /             cni
CIB      DATA  ROOTI(315) / Z4211B854E0F496A0 /
      DATA  ROOTI(316) /       0.177482393492988493 0D+02 /             cni
CIB      DATA  ROOTI(316) / Z4211BF8C9D2ED1A2 /
      DATA  ROOTI(317) /       0.177763888346311774 0D+02 /             cni
CIB      DATA  ROOTI(317) / Z4211C6C16B2DB870 /
      DATA  ROOTI(318) /       0.178044938147648573 0D+02 /              4370   
CIB      DATA  ROOTI(318) / Z4211CDF34E803FD5 /
      DATA  ROOTI(319) /       0.178325545001270065 0D+02 /             cni
CIB      DATA  ROOTI(319) / Z4211D5224AAE2EE1 /
      DATA  ROOTI(320) /       0.178605710994917501 0D+02 /             cni
CIB      DATA  ROOTI(320) / Z4211DC4E63383328 /
      DATA  ROOTI(321) /       0.178885438199983184 0D+02 /             cni
CIB      DATA  ROOTI(321) / Z4211E3779B97F4A8 /
      DATA  ROOTI(322) /       0.179164728671689168 0D+02 /             cni
CIB      DATA  ROOTI(322) / Z4211EA9DF740296F /
      DATA  ROOTI(323) /       0.179443584449263618 0D+02 /              4380   
CIB      DATA  ROOTI(323) / Z4211F1C1799CA8FF /
      DATA  ROOTI(324) /       0.179722007556114285 0D+02 /             cni
CIB      DATA  ROOTI(324) / Z4211F8E226127F61 /
      DATA  ROOTI(325) /       0.180000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(325) / Z4212000000000000 /
      DATA  ROOTI(326) /       0.180277563773199461 0D+02 /             cni
CIB      DATA  ROOTI(326) / Z4212071B0ABCD838 /
      DATA  ROOTI(327) /       0.180554700852677890 0D+02 /             cni
CIB      DATA  ROOTI(327) / Z42120E33499A21A9 /
      DATA  ROOTI(328) /       0.180831413200251241 0D+02 /              4390   
CIB      DATA  ROOTI(328) / Z42121548BFE27445 /
      DATA  ROOTI(329) /       0.181107702762748346 0D+02 /             cni
CIB      DATA  ROOTI(329) / Z42121C5B70D9F824 /
      DATA  ROOTI(330) /       0.181383571472170537 0D+02 /             cni
CIB      DATA  ROOTI(330) / Z4212236B5FBE7711 /
      DATA  ROOTI(331) /       0.181659021245849495 0D+02 /             cni
CIB      DATA  ROOTI(331) / Z42122A788FC76DE5 /
      DATA  ROOTI(332) /       0.181934053986602535 0D+02 /             cni
CIB      DATA  ROOTI(332) / Z4212318304261D9A /
      DATA  ROOTI(333) /       0.182208671582885984 0D+02 /              4400   
CIB      DATA  ROOTI(333) / Z4212388AC0059C28 /
      DATA  ROOTI(334) /       0.182482875908946589 0D+02 /             cni
CIB      DATA  ROOTI(334) / Z42123F8FC68AE52B /
      DATA  ROOTI(335) /       0.182756668824970667 0D+02 /             cni
CIB      DATA  ROOTI(335) / Z421246921AD4EA49 /
      DATA  ROOTI(336) /       0.183030052177231255 0D+02 /             cni
CIB      DATA  ROOTI(336) / Z42124D91BFFCA360 /
      DATA  ROOTI(337) /       0.183303027798233593 0D+02 /             cni
CIB      DATA  ROOTI(337) / Z4212548EB9151E85 /
      DATA  ROOTI(338) /       0.183575597506858195 0D+02 /              4410   
CIB      DATA  ROOTI(338) / Z42125B89092B8FBF /
      DATA  ROOTI(339) /       0.183847763108502349 0D+02 /             cni
CIB      DATA  ROOTI(339) / Z42126280B3476096 /
      DATA  ROOTI(340) /       0.184119526395219673 0D+02 /             cni
CIB      DATA  ROOTI(340) / Z42126975BA6A3F6B /
      DATA  ROOTI(341) /       0.184390889145857741 0D+02 /             cni
CIB      DATA  ROOTI(341) / Z4212706821902E9A /
      DATA  ROOTI(342) /       0.184661853126193876 0D+02 /             cni
CIB      DATA  ROOTI(342) / Z42127757EBAF9368 /
      DATA  ROOTI(343) /       0.184932420089069289 0D+02 /              4420   
CIB      DATA  ROOTI(343) / Z42127E451BB944C3 /
      DATA  ROOTI(344) /       0.185202591774521359 0D+02 /             cni
CIB      DATA  ROOTI(344) / Z4212852FB49899CD /
      DATA  ROOTI(345) /       0.185472369909914079 0D+02 /             cni
CIB      DATA  ROOTI(345) / Z42128C17B9337834 /
      DATA  ROOTI(346) /       0.185741756210067095 0D+02 /             cni
CIB      DATA  ROOTI(346) / Z421292FD2C6A6262 /
      DATA  ROOTI(347) /       0.186010752377382751 0D+02 /             cni
CIB      DATA  ROOTI(347) / Z421299E011188575 /
      DATA  ROOTI(348) /       0.186279360101971569 0D+02 /              4430   
CIB      DATA  ROOTI(348) / Z4212A0C06A13C70B /
      DATA  ROOTI(349) /       0.186547581061776313 0D+02 /             cni
CIB      DATA  ROOTI(349) / Z4212A79E3A2CD2E6 /
      DATA  ROOTI(350) /       0.186815416922694055 0D+02 /             cni
CIB      DATA  ROOTI(350) / Z4212AE79842F2858 /
      DATA  ROOTI(351) /       0.187082869338697080 0D+02 /             cni
CIB      DATA  ROOTI(351) / Z4212B5524AE1278E /
      DATA  ROOTI(352) /       0.187349939951951932 0D+02 /             cni
CIB      DATA  ROOTI(352) / Z4212BC2891041EA7 /
      DATA  ROOTI(353) /       0.187616630392937189 0D+02 /              4440   
CIB      DATA  ROOTI(353) / Z4212C2FC595456A7 /
      DATA  ROOTI(354) /       0.187882942280559355 0D+02 /             cni
CIB      DATA  ROOTI(354) / Z4212C9CDA6892035 /
      DATA  ROOTI(355) /       0.188148877222267785 0D+02 /             cni
CIB      DATA  ROOTI(355) / Z4212D09C7B54E03E /
      DATA  ROOTI(356) /       0.188414436814167736 0D+02 /             cni
CIB      DATA  ROOTI(356) / Z4212D768DA651C63 /
      DATA  ROOTI(357) /       0.188679622641132063 0D+02 /             cni
CIB      DATA  ROOTI(357) / Z4212DE32C6628741 /
      DATA  ROOTI(358) /       0.188944436276911851 0D+02 /              4450   
CIB      DATA  ROOTI(358) / Z4212E4FA41F10C9B /
      DATA  ROOTI(359) /       0.189208879284245022 0D+02 /             cni
CIB      DATA  ROOTI(359) / Z4212EBBF4FAFDD4B /
      DATA  ROOTI(360) /       0.189472953214964157 0D+02 /             cni
CIB      DATA  ROOTI(360) / Z4212F281F2397B1D /
      DATA  ROOTI(361) /       0.189736659610102762 0D+02 /             cni
CIB      DATA  ROOTI(361) / Z4212F9422C23C47E /
      DATA  ROOTI(362) /       0.190000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(362) / Z4213000000000000 /
      DATA  ROOTI(363) /       0.190262975904404463 0D+02 /              4460   
CIB      DATA  ROOTI(363) / Z421306BB705AE7C3 /
      DATA  ROOTI(364) /       0.190525588832576496 0D+02 /             cni
CIB      DATA  ROOTI(364) / Z42130D747FBCB4B5 /
      DATA  ROOTI(365) /       0.190787840283389123 0D+02 /             cni
CIB      DATA  ROOTI(365) / Z4213142B30A929AB /
      DATA  ROOTI(366) /       0.191049731745427991 0D+02 /             cni
CIB      DATA  ROOTI(366) / Z42131ADF859F9E5E /
      DATA  ROOTI(367) /       0.191311264697089918 0D+02 /             cni
CIB      DATA  ROOTI(367) / Z42132191811B0A41 /
      DATA  ROOTI(368) /       0.191572440606680168 0D+02 /              4470   
CIB      DATA  ROOTI(368) / Z4213284125920F33 /
      DATA  ROOTI(369) /       0.191833260932508765 0D+02 /             cni
CIB      DATA  ROOTI(369) / Z42132EEE75770416 /
      DATA  ROOTI(370) /       0.192093727122985456 0D+02 /             cni
CIB      DATA  ROOTI(370) / Z421335997337FF40 /
      DATA  ROOTI(371) /       0.192353840616713434 0D+02 /             cni
CIB      DATA  ROOTI(371) / Z42133C42213EE0C9 /
      DATA  ROOTI(372) /       0.192613602842582239 0D+02 /             cni
CIB      DATA  ROOTI(372) / Z421342E881F15CC2 /
      DATA  ROOTI(373) /       0.192873015219859099 0D+02 /              4480   
CIB      DATA  ROOTI(373) / Z4213498C97B10540 /
      DATA  ROOTI(374) /       0.193132079158279666 0D+02 /             cni
CIB      DATA  ROOTI(374) / Z4213502E64DB5456 /
      DATA  ROOTI(375) /       0.193390796058137155 0D+02 /             cni
CIB      DATA  ROOTI(375) / Z421356CDEBC9B5E2 /
      DATA  ROOTI(376) /       0.193649167310370842 0D+02 /             cni
CIB      DATA  ROOTI(376) / Z42135D6B2ED19148 /
      DATA  ROOTI(377) /       0.193907194296653174 0D+02 /             cni
CIB      DATA  ROOTI(377) / Z4213640630445306 /
      DATA  ROOTI(378) /       0.194164878389475994 0D+02 /              4490   
CIB      DATA  ROOTI(378) / Z42136A9EF26F762F /
      DATA  ROOTI(379) /       0.194422220952235811 0D+02 /             cni
CIB      DATA  ROOTI(379) / Z42137135779C8DCB /
      DATA  ROOTI(380) /       0.194679223339317851 0D+02 /             cni
CIB      DATA  ROOTI(380) / Z421377C9C2114E15 /
      DATA  ROOTI(381) /       0.194935886896179262 0D+02 /             cni
CIB      DATA  ROOTI(381) / Z42137E5BD40F95A1 /
      DATA  ROOTI(382) /       0.195192212959431366 0D+02 /             cni
CIB      DATA  ROOTI(382) / Z421384EBAFD57667 /
      DATA  ROOTI(383) /       0.195448202856920652 0D+02 /              4500   
CIB      DATA  ROOTI(383) / Z42138B79579D3EAB /
      DATA  ROOTI(384) /       0.195703857907809251 0D+02 /             cni
CIB      DATA  ROOTI(384) / Z42139204CD9D81D6 /
      DATA  ROOTI(385) /       0.195959179422654231 0D+02 /             cni
CIB      DATA  ROOTI(385) / Z4213988E1409212E /
      DATA  ROOTI(386) /       0.196214168703485825 0D+02 /             cni
CIB      DATA  ROOTI(386) / Z42139F152D0F5470 /
      DATA  ROOTI(387) /       0.196468827043884993 0D+02 /             cni
CIB      DATA  ROOTI(387) / Z4213A59A1ADBB257 /
      DATA  ROOTI(388) /       0.196723155729060011 0D+02 /              4510   
CIB      DATA  ROOTI(388) / Z4213AC1CDF963908 /
      DATA  ROOTI(389) /       0.196977156035922079 0D+02 /             cni
CIB      DATA  ROOTI(389) / Z4213B29D7D635662 /
      DATA  ROOTI(390) /       0.197230829233160208 0D+02 /             cni
CIB      DATA  ROOTI(390) / Z4213B91BF663F03A /
      DATA  ROOTI(391) /       0.197484176581314976 0D+02 /             cni
CIB      DATA  ROOTI(391) / Z4213BF984CB56C77 /
      DATA  ROOTI(392) /       0.197737199332851894 0D+02 /             cni
CIB      DATA  ROOTI(392) / Z4213C6128271B923 /
      DATA  ROOTI(393) /       0.197989898732233307 0D+02 /              4520   
CIB      DATA  ROOTI(393) / Z4213CC8A99AF5453 /
      DATA  ROOTI(394) /       0.198242276015990093 0D+02 /             cni
CIB      DATA  ROOTI(394) / Z4213D30094815409 /
      DATA  ROOTI(395) /       0.198494332412792076 0D+02 /             cni
CIB      DATA  ROOTI(395) / Z4213D97474F76DF2 /
      DATA  ROOTI(396) /       0.198746069143517907 0D+02 /             cni
CIB      DATA  ROOTI(396) / Z4213DFE63D1DFF15 /
      DATA  ROOTI(397) /       0.198997487421323989 0D+02 /             cni
CIB      DATA  ROOTI(397) / Z4213E655EEFE1367 /
      DATA  ROOTI(398) /       0.199248588451712756 0D+02 /              4530   
CIB      DATA  ROOTI(398) / Z4213ECC38C9D6D4D /
      DATA  ROOTI(399) /       0.199499373432600038 0D+02 /             cni
CIB      DATA  ROOTI(399) / Z4213F32F17FE8D04 /
      DATA  ROOTI(400) /       0.199749843554381776 0D+02 /             cni
CIB      DATA  ROOTI(400) / Z4213F9989320B7F7 /
      DATA  ROOTI(401) /       0.200000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(401) / Z4214000000000000 /
      DATA  ROOTI(402) /       0.200249843945007875 0D+02 /             cni
CIB      DATA  ROOTI(402) / Z4214066560954A8F /
      DATA  ROOTI(403) /       0.200499376557634221 0D+02 /              4540   
CIB      DATA  ROOTI(403) / Z42140CC8B6D657C2 /
      DATA  ROOTI(404) /       0.200748598998847321 0D+02 /             cni
CIB      DATA  ROOTI(404) / Z4214132A04B5C969 /
      DATA  ROOTI(405) /       0.200997512422417799 0D+02 /             cni
CIB      DATA  ROOTI(405) / Z421419894C2329F0 /
      DATA  ROOTI(406) /       0.201246117974981082 0D+02 /             cni
CIB      DATA  ROOTI(406) / Z42141FE68F0AF33D /
      DATA  ROOTI(407) /       0.201494416796098861 0D+02 /             cni
CIB      DATA  ROOTI(407) / Z42142641CF569572 /
      DATA  ROOTI(408) /       0.201742410018320157 0D+02 /              4550   
CIB      DATA  ROOTI(408) / Z42142C9B0EEC7DA4 /
      DATA  ROOTI(409) /       0.201990098767241548 0D+02 /             cni
CIB      DATA  ROOTI(409) / Z421432F24FB01C7A /
      DATA  ROOTI(410) /       0.202237484161566847 0D+02 /             cni
CIB      DATA  ROOTI(410) / Z421439479381ECBD /
      DATA  ROOTI(411) /       0.202484567313165869 0D+02 /             cni
CIB      DATA  ROOTI(411) / Z42143F9ADC3F79CE /
      DATA  ROOTI(412) /       0.202731349327132939 0D+02 /             cni
CIB      DATA  ROOTI(412) / Z421445EC2BC36615 /
      DATA  ROOTI(413) /       0.202977831301844382 0D+02 /              4560   
CIB      DATA  ROOTI(413) / Z42144C3B83E57153 /
      DATA  ROOTI(414) /       0.203224014329015752 0D+02 /             cni
CIB      DATA  ROOTI(414) / Z42145288E67A7EED /
      DATA  ROOTI(415) /       0.203469899493758035 0D+02 /             cni
CIB      DATA  ROOTI(415) / Z421458D455549C1A /
      DATA  ROOTI(416) /       0.203715487874633610 0D+02 /             cni
CIB      DATA  ROOTI(416) / Z42145F1DD243060A /
      DATA  ROOTI(417) /       0.203960780543711380 0D+02 /             cni
CIB      DATA  ROOTI(417) / Z421465655F122FF6 /
      DATA  ROOTI(418) /       0.204205778566621383 0D+02 /              4570   
CIB      DATA  ROOTI(418) / Z42146BAAFD8BC921 /
      DATA  ROOTI(419) /       0.204450483002608721 0D+02 /             cni
CIB      DATA  ROOTI(419) / Z421471EEAF76C2C6 /
      DATA  ROOTI(420) /       0.204694894904587201 0D+02 /             cni
CIB      DATA  ROOTI(420) / Z42147830769755FE /
      DATA  ROOTI(421) /       0.204939015319191959 0D+02 /             cni
CIB      DATA  ROOTI(421) / Z42147E7054AF0989 /
      DATA  ROOTI(422) /       0.205182845286831927 0D+02 /             cni
CIB      DATA  ROOTI(422) / Z421484AE4B7CB793 /
      DATA  ROOTI(423) /       0.205426385841741386 0D+02 /              4580   
CIB      DATA  ROOTI(423) / Z42148AEA5CBC935F /
      DATA  ROOTI(424) /       0.205669638012031335 0D+02 /             cni
CIB      DATA  ROOTI(424) / Z421491248A282EED /
      DATA  ROOTI(425) /       0.205912602819740016 0D+02 /             cni
CIB      DATA  ROOTI(425) / Z4214975CD5768088 /
      DATA  ROOTI(426) /       0.206155281280883038 0D+02 /             cni
CIB      DATA  ROOTI(426) / Z42149D93405BE849 /
      DATA  ROOTI(427) /       0.206397674405502940 0D+02 /             cni
CIB      DATA  ROOTI(427) / Z4214A3C7CC8A358A /
      DATA  ROOTI(428) /       0.206639783197718252 0D+02 /              4590   
CIB      DATA  ROOTI(428) / Z4214A9FA7BB0AC4B /
      DATA  ROOTI(429) /       0.206881608655772027 0D+02 /             cni
CIB      DATA  ROOTI(429) / Z4214B02B4F7C0A88 /
      DATA  ROOTI(430) /       0.207123151772079801 0D+02 /             cni
CIB      DATA  ROOTI(430) / Z4214B65A49968D7F /
      DATA  ROOTI(431) /       0.207364413533277201 0D+02 /             cni
CIB      DATA  ROOTI(431) / Z4214BC876BA7F6EC /
      DATA  ROOTI(432) /       0.207605394920266946 0D+02 /             cni
CIB      DATA  ROOTI(432) / Z4214C2B2B7559234 /
      DATA  ROOTI(433) /       0.207846096908265281 0D+02 /              4600   
CIB      DATA  ROOTI(433) / Z4214C8DC2E423980 /
      DATA  ROOTI(434) /       0.208086520466848128 0D+02 /             cni
CIB      DATA  ROOTI(434) / Z4214CF03D20E5AD0 /
      DATA  ROOTI(435) /       0.208326666559996596 0D+02 /             cni
CIB      DATA  ROOTI(435) / Z4214D529A457FCFC /
      DATA  ROOTI(436) /       0.208566536146142099 0D+02 /             cni
CIB      DATA  ROOTI(436) / Z4214DB4DA6BAC4AA /
      DATA  ROOTI(437) /       0.208806130178211014 0D+02 /             cni
CIB      DATA  ROOTI(437) / Z4214E16FDACFF937 /
      DATA  ROOTI(438) /       0.209045449603668736 0D+02 /              4610   
CIB      DATA  ROOTI(438) / Z4214E790422E898F /
      DATA  ROOTI(439) /       0.209284495364563483 0D+02 /             cni
CIB      DATA  ROOTI(439) / Z4214EDAEDE6B10FE /
      DATA  ROOTI(440) /       0.209523268397569638 0D+02 /             cni
CIB      DATA  ROOTI(440) / Z4214F3CBB117DBF4 /
      DATA  ROOTI(441) /       0.209761769634030308 0D+02 /             cni
CIB      DATA  ROOTI(441) / Z4214F9E6BBC4ECB3 /
      DATA  ROOTI(442) /       0.210000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(442) / Z4215000000000000 /
      DATA  ROOTI(443) /       0.210237960416286391 0D+02 /              4620   
CIB      DATA  ROOTI(443) / Z421506177F5491BB /
      DATA  ROOTI(444) /       0.210475651798491867 0D+02 /             cni
CIB      DATA  ROOTI(444) / Z42150C2D3B4BE170 /
      DATA  ROOTI(445) /       0.210713075057054766 0D+02 /             cni
CIB      DATA  ROOTI(445) / Z42151241356CF6E0 /
      DATA  ROOTI(446) /       0.210950231097289880 0D+02 /             cni
CIB      DATA  ROOTI(446) / Z421518536F3CA675 /
      DATA  ROOTI(447) /       0.211187120819428742 0D+02 /             cni
CIB      DATA  ROOTI(447) / Z42151E63EA3D95B0 /
      DATA  ROOTI(448) /       0.211423745118659738 0D+02 /              4630   
CIB      DATA  ROOTI(448) / Z42152472A7F03F92 /
      DATA  ROOTI(449) /       0.211660104885167257 0D+02 /             cni
CIB      DATA  ROOTI(449) / Z42152A7FA9D2F8EA /
      DATA  ROOTI(450) /       0.211896201004170912 0D+02 /             cni
CIB      DATA  ROOTI(450) / Z4215308AF161F4A5 /
      DATA  ROOTI(451) /       0.212132034355964265 0D+02 /             cni
CIB      DATA  ROOTI(451) / Z4215369480174810 /
      DATA  ROOTI(452) /       0.212367605815953020 0D+02 /             cni
CIB      DATA  ROOTI(452) / Z42153C9C576AEF0B /
      DATA  ROOTI(453) /       0.212602916254693000 0D+02 /              4640   
CIB      DATA  ROOTI(453) / Z421542A278D2D036 /
      DATA  ROOTI(454) /       0.212837966537927628 0D+02 /             cni
CIB      DATA  ROOTI(454) / Z421548A6E5C2C110 /
      DATA  ROOTI(455) /       0.213072757526625161 0D+02 /             cni
CIB      DATA  ROOTI(455) / Z42154EA99FAC8A0F /
      DATA  ROOTI(456) /       0.213307290077015423 0D+02 /             cni
CIB      DATA  ROOTI(456) / Z421554AAA7FFEAAA /
      DATA  ROOTI(457) /       0.213541565040626224 0D+02 /             cni
CIB      DATA  ROOTI(457) / Z42155AAA002A9D5A /
      DATA  ROOTI(458) /       0.213775583264319486 0D+02 /              4650   
CIB      DATA  ROOTI(458) / Z421560A7A9985B93 /
      DATA  ROOTI(459) /       0.214009345590326951 0D+02 /             cni
CIB      DATA  ROOTI(459) / Z421566A3A5B2E1B1 /
      DATA  ROOTI(460) /       0.214242852856285495 0D+02 /             cni
CIB      DATA  ROOTI(460) / Z42156C9DF5E1F2DA /
      DATA  ROOTI(461) /       0.214476105895272156 0D+02 /             cni
CIB      DATA  ROOTI(461) / Z421572969B8B5CD8 /
      DATA  ROOTI(462) /       0.214709105535838880 0D+02 /             cni
CIB      DATA  ROOTI(462) / Z4215788D9812FBEB /
      DATA  ROOTI(463) /       0.214941852602046772 0D+02 /              4660   
CIB      DATA  ROOTI(463) / Z42157E82ECDABE8D /
      DATA  ROOTI(464) /       0.215174347913500128 0D+02 /             cni
CIB      DATA  ROOTI(464) / Z421584769B42A930 /
      DATA  ROOTI(465) /       0.215406592285380150 0D+02 /             cni
CIB      DATA  ROOTI(465) / Z42158A68A4A8D9F3 /
      DATA  ROOTI(466) /       0.215638586528478235 0D+02 /             cni
CIB      DATA  ROOTI(466) / Z421590590A698C4B /
      DATA  ROOTI(467) /       0.215870331449229020 0D+02 /             cni
CIB      DATA  ROOTI(467) / Z42159647CDDF1CA5 /
      DATA  ROOTI(468) /       0.216101827849743096 0D+02 /              4670   
CIB      DATA  ROOTI(468) / Z42159C34F0620BFF /
      DATA  ROOTI(469) /       0.216333076527839374 0D+02 /             cni
CIB      DATA  ROOTI(469) / Z4215A22073490377 /
      DATA  ROOTI(470) /       0.216564078277077137 0D+02 /             cni
CIB      DATA  ROOTI(470) / Z4215A80A57E8D7D1 /
      DATA  ROOTI(471) /       0.216794833886788005 0D+02 /             cni
CIB      DATA  ROOTI(471) / Z4215ADF29F948CFB /
      DATA  ROOTI(472) /       0.217025344142107066 0D+02 /             cni
CIB      DATA  ROOTI(472) / Z4215B3D94B9D5979 /
      DATA  ROOTI(473) /       0.217255609824004310 0D+02 /              4680   
CIB      DATA  ROOTI(473) / Z4215B9BE5D52A9DA /
      DATA  ROOTI(474) /       0.217485631709315470 0D+02 /             cni
CIB      DATA  ROOTI(474) / Z4215BFA1D602241C /
      DATA  ROOTI(475) /       0.217715410570772399 0D+02 /             cni
CIB      DATA  ROOTI(475) / Z4215C583B6F7AB03 /
      DATA  ROOTI(476) /       0.217944947177033690 0D+02 /             cni
CIB      DATA  ROOTI(476) / Z4215CB64017D6177 /
      DATA  ROOTI(477) /       0.218174242292714275 0D+02 /             cni
CIB      DATA  ROOTI(477) / Z4215D142B6DBADC5 /
      DATA  ROOTI(478) /       0.218403296678415551 0D+02 /              4690   
CIB      DATA  ROOTI(478) / Z4215D71FD8593CEF /
      DATA  ROOTI(479) /       0.218632111090754471 0D+02 /             cni
CIB      DATA  ROOTI(479) / Z4215DCFB673B05DF /
      DATA  ROOTI(480) /       0.218860686282392898 0D+02 /             cni
CIB      DATA  ROOTI(480) / Z4215E2D564C44CA1 /
      DATA  ROOTI(481) /       0.219089023002066448 0D+02 /             cni
CIB      DATA  ROOTI(481) / Z4215E8ADD236A58F /
      DATA  ROOTI(482) /       0.219317121994613089 0D+02 /             cni
CIB      DATA  ROOTI(482) / Z4215EE84B0D1F876 /
      DATA  ROOTI(483) /       0.219544984001001495 0D+02 /              4700   
CIB      DATA  ROOTI(483) / Z4215F45A01D483B4 /
      DATA  ROOTI(484) /       0.219772609758359110 0D+02 /             cni
CIB      DATA  ROOTI(484) / Z4215FA2DC67ADF4E /
      DATA  ROOTI(485) /       0.220000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(485) / Z4216000000000000 /
      DATA  ROOTI(486) /       0.220227155455452390 0D+02 /             cni
CIB      DATA  ROOTI(486) / Z421605D0AF9D3A44 /
      DATA  ROOTI(487) /       0.220454076850486018 0D+02 /             cni
CIB      DATA  ROOTI(487) / Z42160B9FD68A4554 /
      DATA  ROOTI(488) /       0.220680764907139100 0D+02 /              4710   
CIB      DATA  ROOTI(488) / Z4216116D75FD3E21 /
      DATA  ROOTI(489) /       0.220907220343745223 0D+02 /             cni
CIB      DATA  ROOTI(489) / Z421617398F2AAA48 /
      DATA  ROOTI(490) /       0.221133443874959816 0D+02 /             cni
CIB      DATA  ROOTI(490) / Z42161D0423457AFB /
      DATA  ROOTI(491) /       0.221359436211786544 0D+02 /             cni
CIB      DATA  ROOTI(491) / Z421622CD337F0FE8 /
      DATA  ROOTI(492) /       0.221585198061603386 0D+02 /             cni
CIB      DATA  ROOTI(492) / Z42162894C1073A17 /
      DATA  ROOTI(493) /       0.221810730128188354 0D+02 /              4720   
CIB      DATA  ROOTI(493) / Z42162E5ACD0C3EBE /
      DATA  ROOTI(494) /       0.222036033111745184 0D+02 /             cni
CIB      DATA  ROOTI(494) / Z4216341F58BADA14 /
      DATA  ROOTI(495) /       0.222261107708928698 0D+02 /             cni
CIB      DATA  ROOTI(495) / Z421639E2653E421B /
      DATA  ROOTI(496) /       0.222485954612869890 0D+02 /             cni
CIB      DATA  ROOTI(496) / Z42163FA3F3C02962 /
      DATA  ROOTI(497) /       0.222710574513200861 0D+02 /             cni
CIB      DATA  ROOTI(497) / Z421645640568C1C3 /
      DATA  ROOTI(498) /       0.222934968096079551 0D+02 /              4730   
CIB      DATA  ROOTI(498) / Z42164B229B5EBF1B /
      DATA  ROOTI(499) /       0.223159136044213966 0D+02 /             cni
CIB      DATA  ROOTI(499) / Z421650DFB6C759F4 /
      DATA  ROOTI(500) /       0.223383079036886762 0D+02 /             cni
CIB      DATA  ROOTI(500) / Z4216569B58C65239 /
      DATA  ROOTI(501) /       0.223606797749978981 0D+02 /             cni
CIB      DATA  ROOTI(501) / Z42165C55827DF1D2 /
      DATA  ROOTI(502) /       0.223830292855993918 0D+02 /             cni
CIB      DATA  ROOTI(502) / Z4216620E350F0F44 /
      DATA  ROOTI(503) /       0.224053565024080790 0D+02 /              4740   
CIB      DATA  ROOTI(503) / Z421667C57199104B /
      DATA  ROOTI(504) /       0.224276614920058037 0D+02 /             cni
CIB      DATA  ROOTI(504) / Z42166D7B3939EC6A /
      DATA  ROOTI(505) /       0.224499443206436489 0D+02 /             cni
CIB      DATA  ROOTI(505) / Z4216732F8D0E2F77 /
      DATA  ROOTI(506) /       0.224722050542442311 0D+02 /             cni
CIB      DATA  ROOTI(506) / Z421678E26E30FC21 /
      DATA  ROOTI(507) /       0.224944437584039854 0D+02 /             cni
CIB      DATA  ROOTI(507) / Z42167E93DDBC0E73 /
      DATA  ROOTI(508) /       0.225166604983954031 0D+02 /              4750   
CIB      DATA  ROOTI(508) / Z42168443DCC7BE4A /
      DATA  ROOTI(509) /       0.225388553391692881 0D+02 /             cni
CIB      DATA  ROOTI(509) / Z421689F26C6B01D0 /
      DATA  ROOTI(510) /       0.225610283453569558 0D+02 /             cni
CIB      DATA  ROOTI(510) / Z42168F9F8DBB6FE7 /
      DATA  ROOTI(511) /       0.225831795812724287 0D+02 /             cni
CIB      DATA  ROOTI(511) / Z4216954B41CD4293 /
      DATA  ROOTI(512) /       0.226053091109146287 0D+02 /             cni
CIB      DATA  ROOTI(512) / Z42169AF589B35963 /
      DATA  ROOTI(513) /       0.226274169979695223 0D+02 /              4760   
CIB      DATA  ROOTI(513) / Z4216A09E667F3BCD /
      DATA  ROOTI(514) /       0.226495033058122495 0D+02 /             cni
CIB      DATA  ROOTI(514) / Z4216A645D9411B85 /
      DATA  ROOTI(515) /       0.226715680975092688 0D+02 /             cni
CIB      DATA  ROOTI(515) / Z4216ABEBE307D6D9 /
      DATA  ROOTI(516) /       0.226936114358204328 0D+02 /             cni
CIB      DATA  ROOTI(516) / Z4216B19084E0FAF9 /
      DATA  ROOTI(517) /       0.227156333832010944 0D+02 /             cni
CIB      DATA  ROOTI(517) / Z4216B733BFD8C648 /
      DATA  ROOTI(518) /       0.227376340018041461 0D+02 /              4770   
CIB      DATA  ROOTI(518) / Z4216BCD594FA2A9A /
      DATA  ROOTI(519) /       0.227596133534820844 0D+02 /             cni
CIB      DATA  ROOTI(519) / Z4216C276054ECF79 /
      DATA  ROOTI(520) /       0.227815714997890346 0D+02 /             cni
CIB      DATA  ROOTI(520) / Z4216C81511DF145F /
      DATA  ROOTI(521) /       0.228035085019827584 0D+02 /             cni
CIB      DATA  ROOTI(521) / Z4216CDB2BBB212EB /
      DATA  ROOTI(522) /       0.228254244210266535 0D+02 /             cni
CIB      DATA  ROOTI(522) / Z4216D34F03CDA114 /
      DATA  ROOTI(523) /       0.228473193175917260 0D+02 /              4780   
CIB      DATA  ROOTI(523) / Z4216D8E9EB365354 /
      DATA  ROOTI(524) /       0.228691932520585439 0D+02 /             cni
CIB      DATA  ROOTI(524) / Z4216DE8372EF7ECE /
      DATA  ROOTI(525) /       0.228910462845191951 0D+02 /             cni
CIB      DATA  ROOTI(525) / Z4216E41B9BFB3B75 /
      DATA  ROOTI(526) /       0.229128784747791983 0D+02 /             cni
CIB      DATA  ROOTI(526) / Z4216E9B2675A6626 /
      DATA  ROOTI(527) /       0.229346898823594287 0D+02 /             cni
CIB      DATA  ROOTI(527) / Z4216EF47D60CA2C6 /
      DATA  ROOTI(528) /       0.229564805664979943 0D+02 /              4790   
CIB      DATA  ROOTI(528) / Z4216F4DBE9105E52 /
      DATA  ROOTI(529) /       0.229782505861521145 0D+02 /             cni
CIB      DATA  ROOTI(529) / Z4216FA6EA162D0F0 /
      DATA  ROOTI(530) /       0.230000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(530) / Z4217000000000000 /
      DATA  ROOTI(531) /       0.230217288664426754 0D+02 /             cni
CIB      DATA  ROOTI(531) / Z4217059005E2C01D /
      DATA  ROOTI(532) /       0.230434372436058261 0D+02 /             cni
CIB      DATA  ROOTI(532) / Z42170B1EB404B725 /
      DATA  ROOTI(533) /       0.230651251893415932 0D+02 /              4800   
CIB      DATA  ROOTI(533) / Z421710AC0B5E5E32 /
      DATA  ROOTI(534) /       0.230867927612303916 0D+02 /             cni
CIB      DATA  ROOTI(534) / Z421716380CE7039A /
      DATA  ROOTI(535) /       0.231084400165826871 0D+02 /             cni
CIB      DATA  ROOTI(535) / Z42171BC2B994CCE3 /
      DATA  ROOTI(536) /       0.231300670124407546 0D+02 /             cni
CIB      DATA  ROOTI(536) / Z4217214C125CB8B2 /
      DATA  ROOTI(537) /       0.231516738055804510 0D+02 /             cni
CIB      DATA  ROOTI(537) / Z421726D41832A0BE /
      DATA  ROOTI(538) /       0.231732604525129346 0D+02 /              4810   
CIB      DATA  ROOTI(538) / Z42172C5ACC093BB4 /
      DATA  ROOTI(539) /       0.231948270094864029 0D+02 /             cni
CIB      DATA  ROOTI(539) / Z421731E02ED21F20 /
      DATA  ROOTI(540) /       0.232163735324878004 0D+02 /             cni
CIB      DATA  ROOTI(540) / Z42173764417DC14E /
      DATA  ROOTI(541) /       0.232379000772445004 0D+02 /             cni
CIB      DATA  ROOTI(541) / Z42173CE704FB7B23 /
      DATA  ROOTI(542) /       0.232594066992260160 0D+02 /             cni
CIB      DATA  ROOTI(542) / Z421742687A3989FF /
      DATA  ROOTI(543) /       0.232808934536456320 0D+02 /              4820   
CIB      DATA  ROOTI(543) / Z421747E8A2251188 /
      DATA  ROOTI(544) /       0.233023603954620881 0D+02 /             cni
CIB      DATA  ROOTI(544) / Z42174D677DAA1D84 /
      DATA  ROOTI(545) /       0.233238075793812030 0D+02 /             cni
CIB      DATA  ROOTI(545) / Z421752E50DB3A3A2 /
      DATA  ROOTI(546) /       0.233452350598575045 0D+02 /             cni
CIB      DATA  ROOTI(546) / Z42175861532B8545 /
      DATA  ROOTI(547) /       0.233666428910958466 0D+02 /             cni
CIB      DATA  ROOTI(547) / Z42175DDC4EFA914B /
      DATA  ROOTI(548) /       0.233880311270530008 0D+02 /              4830   
CIB      DATA  ROOTI(548) / Z42176356020885CD /
      DATA  ROOTI(549) /       0.234093998214392514 0D+02 /             cni
CIB      DATA  ROOTI(549) / Z421768CE6D3C11E0 /
      DATA  ROOTI(550) /       0.234307490277199619 0D+02 /             cni
CIB      DATA  ROOTI(550) / Z42176E45917AD74E /
      DATA  ROOTI(551) /       0.234520787991171495 0D+02 /             cni
CIB      DATA  ROOTI(551) / Z421773BB6FA96C51 /
      DATA  ROOTI(552) /       0.234733891886110051 0D+02 /             cni
CIB      DATA  ROOTI(552) / Z4217793008AB5D3F /
      DATA  ROOTI(553) /       0.234946802489414601 0D+02 /              4840   
CIB      DATA  ROOTI(553) / Z42177EA35D632E43 /
      DATA  ROOTI(554) /       0.235159520326096931 0D+02 /             cni
CIB      DATA  ROOTI(554) / Z421784156EB25D05 /
      DATA  ROOTI(555) /       0.235372045918796395 0D+02 /             cni
CIB      DATA  ROOTI(555) / Z421789863D796253 /
      DATA  ROOTI(556) /       0.235584379787794944 0D+02 /             cni
CIB      DATA  ROOTI(556) / Z42178EF5CA97B3C8 /
      DATA  ROOTI(557) /       0.235796522451031905 0D+02 /             cni
CIB      DATA  ROOTI(557) / Z4217946416EBC56C /
      DATA  ROOTI(558) /       0.236008474424118937 0D+02 /              4850   
CIB      DATA  ROOTI(558) / Z421799D123530B59 /
      DATA  ROOTI(559) /       0.236220236220354316 0D+02 /             cni
CIB      DATA  ROOTI(559) / Z42179F3CF0A9FB4D /
      DATA  ROOTI(560) /       0.236431808350737782 0D+02 /             cni
CIB      DATA  ROOTI(560) / Z4217A4A77FCC0E4C /
      DATA  ROOTI(561) /       0.236643191323984645 0D+02 /             cni
CIB      DATA  ROOTI(561) / Z4217AA10D193C22D /
      DATA  ROOTI(562) /       0.236854385646540209 0D+02 /             cni
CIB      DATA  ROOTI(562) / Z4217AF78E6DA9B30 /
      DATA  ROOTI(563) /       0.237065391822593945 0D+02 /              4860   
CIB      DATA  ROOTI(563) / Z4217B4DFC079258D /
      DATA  ROOTI(564) /       0.237276210354093458 0D+02 /             cni
CIB      DATA  ROOTI(564) / Z4217BA455F46F6FD /
      DATA  ROOTI(565) /       0.237486841740758337 0D+02 /             cni
CIB      DATA  ROOTI(565) / Z4217BFA9C41AB040 /
      DATA  ROOTI(566) /       0.237697286480094263 0D+02 /             cni
CIB      DATA  ROOTI(566) / Z4217C50CEFC9FEAA /
      DATA  ROOTI(567) /       0.237907545067406367 0D+02 /             cni
CIB      DATA  ROOTI(567) / Z4217CA6EE3299D9B /
      DATA  ROOTI(568) /       0.238117617995813156 0D+02 /              4870   
CIB      DATA  ROOTI(568) / Z4217CFCF9F0D5807 /
      DATA  ROOTI(569) /       0.238327505756259690 0D+02 /             cni
CIB      DATA  ROOTI(569) / Z4217D52F244809E9 /
      DATA  ROOTI(570) /       0.238537208837531267 0D+02 /             cni
CIB      DATA  ROOTI(570) / Z4217DA8D73ABA1C4 /
      DATA  ROOTI(571) /       0.238746727726266457 0D+02 /             cni
CIB      DATA  ROOTI(571) / Z4217DFEA8E092212 /
      DATA  ROOTI(572) /       0.238956062906970423 0D+02 /             cni
CIB      DATA  ROOTI(572) / Z4217E5467430A2BB /
      DATA  ROOTI(573) /       0.239165214862027966 0D+02 /              4880   
CIB      DATA  ROOTI(573) / Z4217EAA126F15284 /
      DATA  ROOTI(574) /       0.239374184071716485 0D+02 /             cni
CIB      DATA  ROOTI(574) / Z4217EFFAA719787C /
      DATA  ROOTI(575) /       0.239582971014218771 0D+02 /             cni
CIB      DATA  ROOTI(575) / Z4217F552F5767564 /
      DATA  ROOTI(576) /       0.239791576165635973 0D+02 /             cni
CIB      DATA  ROOTI(576) / Z4217FAAA12D4C51C /
      DATA  ROOTI(577) /       0.240000000000000000 0D+02 /             cni
CIB      DATA  ROOTI(577) / Z4218000000000000 /
      DATA  ROOTI(578) /       0.240208242989286269 0D+02 /              4890   
CIB      DATA  ROOTI(578) / Z42180554BDC2DC4F /
      DATA  ROOTI(579) /       0.240416305603426146 0D+02 /             cni
CIB      DATA  ROOTI(579) / Z42180AA84CE72F89 /
      DATA  ROOTI(580) /       0.240624188310319305 0D+02 /             cni
CIB      DATA  ROOTI(580) / Z42180FFAAE35EFCB /
      DATA  ROOTI(581) /       0.240831891575845916 0D+02 /             cni
CIB      DATA  ROOTI(581) / Z4218154BE2773526 /
      DATA  ROOTI(582) /       0.241039415863879007 0D+02 /             cni
CIB      DATA  ROOTI(582) / Z42181A9BEA723AFB /
      DATA  ROOTI(583) /       0.241246761636296370 0D+02 /              4900   
CIB      DATA  ROOTI(583) / Z42181FEAC6ED614A /
      DATA  ROOTI(584) /       0.241453929352992738 0D+02 /             cni
CIB      DATA  ROOTI(584) / Z4218253878AE2E09 /
      DATA  ROOTI(585) /       0.241660919471891447 0D+02 /             cni
CIB      DATA  ROOTI(585) / Z42182A8500794E6C /
      DATA  ROOTI(586) /       0.241867732448956474 0D+02 /             cni
CIB      DATA  ROOTI(586) / Z42182FD05F129838 /
      DATA  ROOTI(587) /       0.242074368738204093 0D+02 /             cni
CIB      DATA  ROOTI(587) / Z4218351A953D0B0B /
      DATA  ROOTI(588) /       0.242280828791714349 0D+02 /              4910   
CIB      DATA  ROOTI(588) / Z42183A63A3BAD19F /
      DATA  ROOTI(589) /       0.242487113059642816 0D+02 /             cni
CIB      DATA  ROOTI(589) / Z42183FAB8B4D4315 /
      DATA  ROOTI(590) /       0.242693221990231933 0D+02 /             cni
CIB      DATA  ROOTI(590) / Z421844F24CB4E434 /
      DATA  ROOTI(591) /       0.242899156029822372 0D+02 /             cni
CIB      DATA  ROOTI(591) / Z42184A37E8B168A9 /
      DATA  ROOTI(592) /       0.243104915622864368 0D+02 /             cni
CIB      DATA  ROOTI(592) / Z42184F7C6001B446 /
      DATA  ROOTI(593) /       0.243310501211928774 0D+02 /              4920   
CIB      DATA  ROOTI(593) / Z421854BFB363DC39 /
      DATA  ROOTI(594) /       0.243515913237718422 0D+02 /             cni
CIB      DATA  ROOTI(594) / Z42185A01E395284C /
      DATA  ROOTI(595) /       0.243721152139078825 0D+02 /             cni
CIB      DATA  ROOTI(595) / Z42185F42F1521412 /
      DATA  ROOTI(596) /       0.243926218353009361 0D+02 /             cni
CIB      DATA  ROOTI(596) / Z42186482DD565022 /
      DATA  ROOTI(597) /       0.244131112314674041 0D+02 /             cni
CIB      DATA  ROOTI(597) / Z421869C1A85CC346 /
      DATA  ROOTI(598) /       0.244335834457412311 0D+02 /              4930   
CIB      DATA  ROOTI(598) / Z42186EFF531F8BAB /
      DATA  ROOTI(599) /       0.244540385212749669 0D+02 /             cni
CIB      DATA  ROOTI(599) / Z4218743BDE58000C /
      DATA  ROOTI(600) /       0.244744765010408329 0D+02 /             cni
CIB      DATA  ROOTI(600) / Z421879774ABEB0DE /
      DATA  ROOTI(601) /       0.244948974278317806 0D+02 /             cni
CIB      DATA  ROOTI(601) / Z42187EB1990B697A /
      DATA  ROOTI(602) /       0.245153013442625252 0D+02 /             cni
CIB      DATA  ROOTI(602) / Z421883EAC9F53140 /
C                                                                        4940   
C
      DATA X /0.D0/
C
      END
C*ID* SPLNCB   C109.PTOLEMY.FORTLIB                         PER704  21:
      SUBROUTINE SPLNCB (N, X, Y, B, C, D)
C
C     INTERPOLATION BY CUBIC SPLINES.
C
C       THIS ROUTINE IS A DRASTIC SIMPLIFICATION OF THE AMD ROUTINE      4950   
C     "SMOOTH".  IT WILL FIND THE COEFICENTS OF INTERPOLATING CUBICS
C     FOR A FUNCTION TABULATED ON A MONOTONIC GRID OF X'S.  THE
C     ROUTINE "INTRPC" MAY THEN BE USED TO FIND INTERPOLATED VALUES
C     OF THE FUNCTION USING THESE COEFFICENTS.  THE INTERPOLATION
C     SCHEME IS THE METHOD OF CUBIC SPLINES.  THE COEFICENTS FOUND BY
C     THIS ROUTINE ARE THE SAME AS THOSE FOUND BY THE AMD "SMOOTH"
C     ROUTINE FOR THE CASE OF   S = 0   IF THE SPLNCB ENTRY IS USED.
C
C      THE INPUT IS
C                                                                        4960   
C     N - THE NUMBER OF INPUT POINTS.  N MUST BE >= 2.
C         IF N <= 1, NOTHING IS DONE (IMMEADIATE RETURN).
C     X - THE ARRAY OF INPUT INDEPENDENT VARIABLES.  X
C         MUST BE MONOTONIC (INCREASING OR DECREASING) WITH
C         NO EQUAL POINTS.  IT MUST HAVE "N" ELEMENTS.
C     Y - THE ARRAY OF FUNCTION VALUES CORRESPONDING TO X.
C     B, C, D - THE OUTPUT ARRAYS WITH THE FITTING CUBICS.
C         THESE MUST EACH BE "N" ELEMENT ARRAYS.  THE LAST
C         ELEMENT OF EACH ARRAY WILL BE SET TO THE COEFICENTS
C         CORRESPONDING THE THE CUBIC FOR THE N-1 INTERVAL.              4970   
C
C     THE INTERPOLATING CUBICS ARE:
C
C      Y(X) = Y(I) + DEL*B(I) + DEL**2*C(I) +
C             DEL**3*D(I)
C     WHERE
C       DEL = X - X(I)
C     AND
C       X(I) <= X <= X(I+1),  1 <= I <= N-1
C                                                                        4980   
C
C     AN ALTERNATIVE ENTRY IS:
C
C     CALL  SPLNCX (N, X, Y, B, C, D, IOPT, Y1P, YNP)
C
C     WHERE THE FIRST 6 ARGUMENTS ARE AS ABOVE AND THE LAST 3 ARE
C
C     IOPT - CONTROLS OPTIONS CONCERNING THE BOUNDRY CONDITIONS:
C          = 1 MEANS THAT THE 2ND DERIVATIVES ARE ZERO AT X(1) AND
C              X(N).  THIS IS THE "NATURAL" CHOICE AND IS ALSO THE       4990   
C              RESULT OF CALLING SPLNCB.
C          = 2 MEANS THAT THE 3RD DERIVATIES ARE ZERO AT X(1) AND
C              X(N).  THIS ALLOWS AN EXACT FIT TO PURE PARABOLAS
C              WHICH IS NOT POSSIBLE WITH OPTION 1.
C          = 3 MEANS THAT THE 1ST DERIVATIVES ARE BEING SPECIFIED AT
C              X(1) AND X(N).
C     Y1P, YNP - FOR IOPT=3 THESE ARE THE REQUIRED 1ST DERIVATIVS AT
C              X(1) AND X(N).  THEY ARE NOT USED FOR IOPT = 1 OR 2.
C
C                                                                        5000   
C     7/4/73 - DRASTIC CUTS TO AMD "SMOOTH" ROUTINE. - S. PIEPER
C     9/5/73 - ADD SPLNCX ENTRY.
C     4/4/78 - PROTECT AGAINST UNDERFLOWS AT BOTH ENDS;
C              REMOVE OPTIONAL ENTRY FOR PTOLEMY VERSION - S.P.
C     5/19/78 - FIX COUNTS IN ABOVE FOR ALL ZERO ARRAYS - S.P.
C     6/28/78 - DON=T USE UNDEFINED D(NM1) AFTER 300 LOOP
C
C
      implicit real*8 ( a-h, o-z )                                      implicit
      DIMENSION X(1),Y(1),B(1),C(1),D(1)                                 5010   
C
C     NATURAL CUBIC BY DEFAULT
C
CXX 40   ASSIGN 300 TO IGO
CXX      IOPT = 1
C
  50  IF (N .LE. 1)  RETURN
C
      NM1=N-1
C                                                                        5020   
C     LONG STRINGS OF CONSTANT VALUES CAN CAUSE UNDERFLOW PROBLEMS;
C     WE REMOVE SUCH STRINGS FROM THE START.
C
      DO 74  I = 1, NM1
         B(I) = 0
         C(I) = 0
         D(I) = 0
         IF ( 1.D15*DABS(Y(I+1)-Y(I)) .GT. DABS(Y(I)) )  GO TO 80
 74   CONTINUE
      I = NM1                                                            5030   
C
 80   IBASE = I-1
      DO 84  I1 = I, NM1
         J = N+I - I1
         B(J) = 0
         C(J) = 0
         D(J) = 0
         IF ( 1.D15*DABS(Y(J-1)-Y(J)) .GT. DABS(Y(J)) )  GO TO 90
 84   CONTINUE
  90  NUSE = J - IBASE                                                   5040   
      IF ( NUSE .LE. 1 )  RETURN
      NM1 = NUSE-1
C
      H=X(IBASE+2)-X(IBASE+1)
      F=(Y(IBASE+2)-Y(IBASE+1))/H
      IF ( NUSE .EQ. 2)  GO TO 800
      DO 399 I=2,NM1
      G=H
      H=X(IBASE+I+1)-X(IBASE+I)
      E=F                                                                5050   
      F=(Y(IBASE+I+1)-Y(IBASE+I))/H
      GBY3 = G/3
      D(IBASE+I-1) = GBY3*B(IBASE+I-1)
      EPSIM1 = G+H
      RIM1B3 = F-E
C
CXX      GO TO  IGO, (200, 220, 300)
CXXC
CXXC     SPECIAL OPTIONS CHOOSEN:  END POINTS MUST BE PATCHED UP
CXXC                                                                     5060   
CXXC
CXXC     2)  Y'''(X1) = Y'''(XN) = 0
CXXC
CXX 200  IF (I .EQ. 2)  EPSIM1 = EPSIM1 + .5*G
CXX      IF (I .EQ. NM1)  EPSIM1 = EPSIM1 + .5*H
CXX      GO TO 300
CXXC
CXXC     3)  Y'(X1)  AND  Y'(XN)  FED IN
CXXC
CXX 220  IF (I .NE. 2)  GO TO 230                                        5070   
CXX      EPSIM1 = .75*G + H
CXX      RIM1B3 = F - 1.5*E - .5*Y1P
CXX      GO TO 300
CXX 230  IF (I .NE. NM1)  GO TO 300
CXX      EPSIM1 = G  +  .75*H
CXX      RIM1B3 = 1.5*F - E - .5*YNP
C
 300  B(IBASE+I) = 1/((2.D0/3.D0)*EPSIM1 - GBY3*D(IBASE+I-1))
      C(IBASE+I) = RIM1B3 - D(IBASE+I-1)*C(IBASE+I-1)
  399 CONTINUE                                                           5080   
      D(IBASE+NM1) = 0
C
      DO 500 I1=2,NM1
      I=NM1+2-I1
      C(IBASE+I)=B(IBASE+I)*C(IBASE+I)-D(IBASE+I)*C(IBASE+I+1)
  500 CONTINUE
C
CXX      GO TO (600, 520, 530),  IOPT
CXXC
CXXC     FIX UP C(IBASE+NUSE) FOR THE SPECIAL OPTIONS                    5090   
CXXC
CXX 520  C(IBASE+NUSE) = C(IBASE+NM1)
CXX      C(IBASE+1) = C(IBASE+2)
CXX      GO TO 600
CXX 530  C(IBASE+NUSE) = 1.5*(YNP-F)/H - .5*C(IBASE+NM1)
CXX      C(IBASE+1) = -1.5*(Y1P - (Y(IBASE+2)-Y(IBASE+1))/(X(IBASE+2)
CXX     1  -X(IBASE+1)))/(X(IBASE+2)-X(IBASE+1))
CXX     1   -.5*C(IBASE+2)
CXX      GO TO 600
C                                                                        5100   
 600  DO 699 I=1,NM1
      H=X(IBASE+I+1)-X(IBASE+I)
      D(IBASE+I)=(C(IBASE+I+1)-C(IBASE+I))/(3.0*H)
      B(IBASE+I)=(Y(IBASE+I+1)-Y(IBASE+I))/H-(H*D(IBASE+I)+C(IBASE+I))*H
 699  CONTINUE
C
C     DEFINE THE N'TH CUBIC AS THE N-1'TH CUBIC
C
      D(IBASE+NUSE) = D(IBASE+NM1)
      B(IBASE+NUSE) = B(IBASE+NM1) + (2*C(IBASE+NM1) + 3*D(IBASE+NM1)*H) 5110   
     1   * H
C
      RETURN
C
C     SPECIAL CASE FOR N = 2.
C
 800  D(IBASE+1) = 0
      B(IBASE+1) = F
      RETURN
C                                                                        5120   
C     ENTRY FOR EXTRA OPTIONS
C
CXX      ENTRY  SPLNCX (N, X, Y, B, C, D, IIOPT, Y1P, YNP)
CXXC
CXX      IOPT = IIOPT
CXX      IF (IOPT .LT. 1  .OR.  IOPT .GT. 3)  RETURN
CXX      GO TO (40, 920, 930), IOPT
CXX 920  ASSIGN 200 TO IGO
CXX      GO TO 50
CXX 930  ASSIGN 220 TO IGO                                               5130   
CXX      GO TO 50
C
      END
C*ID* SYSERR   C109.PTOLEMY.FORTLIB                         PER704  09:
       SUBROUTINE SYSERR
C
C      7/14/72  ARGONNE VERSION
C     7/9/80 - CRAY VERSION ADDED - S.P.
C
C      GIVES A TRACE BACK AND STOPS EXECUTION                            5140   
C
      write (6, 10)
  10  FORMAT ('0' / '0' / '0$*$*$*$*$*$  SYSERR CALLED',
     1  5X, 10('$*$*$*$*$*') )
      STOP 9876
      END
C*ID* YLMSUB   C109.PTOLEMY.FORTLIB                         PER704  21:
      SUBROUTINE YLMSUB (LMAX, MMAX, Z, YLMRAY)
C
C                                                                        5150   
C     COMPUTES THE SHERICAL HARMONICS Y(L,M)(COS THETA = Z, PHI = 0)
C      FOR M = 0, MMAX AND L = M, LMAX
C     THE CONDON-SHORTLY PHASE CONVENTION IS USED (SEE JACKSON,
C     "CLASSICAL ELECTRODYNAMICS", P65 FF).
C
C     NOTE...  THESE Y(LM) ARE FOR THE AZIMUTHAL ANGLE (PHI)
C              EQUAL TO ZERO AND HENCE ARE REAL.
C              THUS THEY ARE REALLY JUST THE P(LM) WITH DIFFERENT
C              NORMALIZATIONS
C                                                                        5160   
C     THE ARGUMENTS ARE:
C
C     LMAX - THE MAXIMUM L VALUE.  LMAX MUST SATISFY
C            0 <= LMAX    OR ELSE AN ERROR MESSAGE WILL
C            BE PRINTED AND THE STEP WILL BE TERMINATED.
C
C     MMAX - THE MAXIMUM M VALUE.  MMAX SHOULD SATISFY
C            0 <= MMAX <= LMAX  BUT IF  MMAX > LMAX, PROCESSING
C            WILL CONTINUE AS IF MMAX = LMAX  AND NO ERROR MESSAGE
C            WILL BE PRINTED.  HOWEVER,  MMAX < 0  DOES RESULT IN        5170   
C            A FATAL ERROR MESSAGE.  NOTE THAT THE Y(L,M) ARE GENERATED
C            ONLY FOR  M > 0.  VALUES OF Y(L,M) FOR M < 0 MAY BE FOUND
C           FROM
C            Y(L,-M) = (-1)**M Y(L,M)(Z,PHI=0)
C           FOR LARGE M AND |Z| THIS ROUTINE MAY UNDERFLOW SINCE
C            Y(L,M)(Z) ---> (1-Z**2)**(M/2)
C
C     Z - THE COSINE OF THE ANGLE THETA ATWHICH  Y(L,M)(THETA, PHI=0)
C         ARE TO BE GENERATED.  Z MUST SATISFY  -1 <= Z <= 1
C         (TO WITHIN 1.E-14) OR ELSE A FATAL ERROR MESSAGE IS PRINTED.   5180   
C
C     YLMRAY - A ONE DIMENSION ARRAY OF DIMENSION ATLEAST
C             ( (MMAX+1) * (2*LMAX - MMAX + 2) ) / 2
C              THAT WILL RECIEVE THE Y(L, M).
C              THEY ARE STORED IN YLMRAY AS (READ L AS LMAX)
C        Y00, Y10, Y20,... YL0, Y11, Y21, ... YL1, Y22, Y32, ... ... YLL
C     THUS
C              Y(L,M) = YLMRAY(L+1 + M*(2*LMAX+1-M)/2)
C
C     FOR A Y(L,M) WITH MORE CONVIENT (AND HENCE MORE WASTEFUL OF        5190   
C     CORE) OUTPUT SEE  YLM2D  WHICH GENERATE A 2-DIMENSIONAL
C     ARRAY OF Y(L,M), INCLUDING NEGATIVE M'S IF DESIRED.
C
C
C     NOTE THAT BOTH Z AND YLMRAY ARE  REAL*8 VARIABLES.
C
C
C     THE TIME REQUIRED FOR THIS ROUTINE IS APPROXIMATELY
C       16. + 4.9 NUM  MICROSECONDS ON THE /195 WHERE NUM IS
C     THE TOTAL NUMBER OF Y(L,M)'S TO BE FOUND (SEE MINIMUM DIMENSION    5200   
C     OF YLMRAY ABOVE).
C
C       THE RELATIVE ERRORS IN THE Y(L,M) ARE ONLY WEAKLY DEPENDANT ON
C     THE VARIABLES M, AND Z.  TYPICAL RMS AND MAXIMUM RELATIVE
C     ERRORS AS A FUNCTION OF L (AVERAGED OVER ALL ALLOWED M) ARE:
C      L    RMS     MAX
C       5   7E-16   1E-15
C      10   1E-15   1.5E-15
C      20   2E-15   5E-15
C      50   3E-14   2E-13                                                5210   
C     100  1.4E-14  8E-14
C
C     HOWEVER THERE APPEAR TO BE ISOLATED CASES OF LARGER RELATIVE ERROR
C     A FEW EXAMPLES ARE:
C       L   M   RELATIVE ERROR
C        9  3     2E-14
C       27  5     2E-12
C       98 21     2E-12
C     THE ABOVE ARE THE LARGEST ERRORS ENCOUNTERED FOR ALL (L,M)
C     LESS THAN 20, 75, AND 100, RESPECTIVELY.                           5220   
C     ALL OF THE ABOVE ERRORS ARE FOR  Z = .5  BUT SIMILAR RESULTS
C     OBTAIN FOR OTHER VALUES OF Z.
C
C     IF  LMAX <= 100  THIS ROUTINE USES A BLOCK DATA AREA THAT
C     CONTAINS THE SQUARE ROOTS OF THE FIRST 201 INTEGERS.  FOR
C     LMAX > 100 THE SQUARE ROOTS ARE COMPUTED EVERY TIME AND THUS THE
C     ROUTINE IS CONSIDERABLY SLOWER.  THIS BOUNDRY MAY BE CHANGED FROM
C     100 BY RECOMPILING THE  /ROOTIS/  COMMON BLOCK.
C
C     EXTERNAL REFERENCES:                                               5230   
C       DSQRT
C       IBCOM#  FOR ERROR MESSAGE PRINTING
C       ROOTIS  A BLOCK DATA COMMON BLOCK
C
C     BY S. PIEPER
C
C     9/25/72 - NEW, REVISED (IMPROVED) RECURSION SCHEME.
C     4/20/74 - DOUBLE PRECISION, BLOCK DATA USED.
C     10/28/74 - ALLOW ARBITRARY L.
C     2/1/77 - GREATLY IMPROVED VERSION FOR L > 100.                     5240   
C     12/14/77 - CDC STUFF IN MORTRAN
C
C
      IMPLICIT REAL*8  (A-H, O-Z)
C
      DIMENSION  YLMRAY(1)
C
C
C
C     FOLLOWING TABLE OF SQUARE ROOTS OF INTEGERS IS SETUP BY THE        5250   
C     BLOCK COMMON DRAGGED IN WITH THE DUMMY /ROOTIX/
C         ROOTI(I) = SQRT(I-1)
C
      COMMON /ROOTIS/  IMXDIM, ISPACE, ROOTI(1)
C
CIB      EXTERNAL  ROOTIX
C
      DATA  ORT4PI / .28209 47917 73878 140 0D+0 /                      cni
CIB      DATA  ORT4PI / Z4048375D 410A6DB4 /
C                                                                        5260   
      DATA IZERO /0/
C
C     CHECK THE ARGUMENTS FOR VALIDITY
C
      ZSQ = Z**2
      IF (LMAX .LT. 0
     1  .OR.  MMAX .LT. 0  .OR.  ZSQ .GT. 1.000000000002 0D+0)
     2      GO TO 120
C
      GO TO 200                                                          5270   
C
 120  write (6, 123) LMAX, MMAX, I, IMXDIM, Z
 123  FORMAT ('0$*$*$* INVALID INPUT TO YLMSUB:', 4I15, G20.10)
      CALL SYSERR
C
C
 200  IF (MMAX .NE. 0)  ROOT = DSQRT(DABS(1-ZSQ))
C
      YLMRAY(1) = ORT4PI
      IF (LMAX .EQ. 0)  RETURN                                           5280   
C
      I = 1
C
C     THE MAJOR LOOP IS ON M
C
      LMMAX = MIN0 (MMAX, LMAX-1)
C
      MMAX1 = IMXDIM/2 - 3
      LMAX1 = MIN0( LMAX, MMAX1 )
C                                                                        5290   
      DM = 0
C
      DO 899  M = IZERO, LMMAX
C
         J = I
         I = I+1
         MP2 = M + 2
         IF ( M .GT. MMAX1 )  GO TO 500
C
C     GENERATE THE Y(L=M+1, M) FOR THIS M                                5300   
C
         YLMRAY(I) = ROOTI(4+2*M)*Z*YLMRAY(I-1)
         IF (LMAX .LT. MP2)  GO TO 850
C
C     NOW DO L = M+2, M+3, ..., LMAX FOR THIS M
C
         ALAST = 1/ROOTI(2*MP2)
C
         DO 399  L = MP2, LMAX1
            I = I+1                                                      5310   
            A = ROOTI(L-M+1)*ROOTI(L+M+1) / (ROOTI(2*L)*ROOTI(2*L+2))
            YLMRAY(I) = ( Z*YLMRAY(I-1) - ALAST*YLMRAY(I-2) ) / A
            ALAST = A
 399     CONTINUE
C
         IF ( LMAX1 .GE. LMAX )  GO TO 850
         L1 = LMAX1+1
         GO TO 600
C
C                                                                        5320   
C
C     HERE WE HAVE A LARGE LMAX AND WE COMPUTE THE REQUIRED SQUARE
C     ROOTS EVERY TIME.
C
C
C     GENERATE THE Y(L=M+1, M) FOR THIS M
C
 500     YLMRAY(I) = DSQRT(3+DM+DM)*Z*YLMRAY(I-1)
         IF (LMAX .LT. MP2)  GO TO 850
         L1 = MP2                                                        5330   
C
C     NOW DO L = M+2, M+3, ..., LMAX FOR THIS M
 600     DL = L1
         ALAST = DSQRT( (DL-1-DM)*(DL-1+DM) / ((DL+DL-1)*(DL+DL-3)) )
C
         DO 799  L = L1, LMAX
            I = I+1
            A = DSQRT( (DL-DM)*(DL+DM) / ((DL+DL+1)*(DL+DL-1)) )
            YLMRAY(I) = ( Z*YLMRAY(I-1) - ALAST*YLMRAY(I-2) ) / A
            ALAST = A                                                    5340   
            DL = DL + 1
 799     CONTINUE
C
 850     IF (M .EQ. MMAX)  RETURN
C
C    NOW GENERATE THE NEXT Y(L=M, M)
C    NOTE THAT HERE THE DESIRED M VALUE IS REALLY M+1 SINCE
C     THE DO LOOP HAS NOT YET INCREMENTED M
C
         I = I+1                                                         5350   
         IF ( M .GT. MMAX1 )  GO TO 880
C     J POINTS TO THE LAST Y(M, M)
         YLMRAY(I) = -(ROOTI(4+2*M)/ROOTI(3+2*M)) * ROOT*YLMRAY(J)
         GO TO 890
 880     YLMRAY(I) = -(DSQRT(3+DM+DM)/DSQRT(2+DM+DM)) * ROOT*YLMRAY(J)
C
 890     DM = DM + 1
 899  CONTINUE
C
      RETURN                                                             5360   
C
      END
