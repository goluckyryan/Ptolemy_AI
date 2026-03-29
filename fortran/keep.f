C*ID* KEEPIT   C109.PTOLEMY.SOURCE                          PER701  22:
      SUBROUTINE KEEPIT
C
C     PERFORMS SPEAKEASY KEEP OPERATION FOR PTOLEMY
C
C     PTOLEMY INPUT IS
C
C     KEEP  OBJECTNAME  MEMBERNAME
C
C     OBJECTNAME MAY BE EITHER THE NAME OF AN OBJECT IN THE ALLOCATOR      10   
C     OR ONE OF THE FOLLOWING SPECIAL NAMES
C       ANGLEGRID - KEEPS ANGLESTEP, ANGLEMIN, ANGLEMAX
C       HEADER - KEEPS THE HEADER WITHOUT LEADING BLANKS
C       REACTION - KEEPS THE REACTION WITHOUT TRAILING BLANKS
C       ELAB
C
C     THE PDS IS REFERED TO BY THE DDNAME MYKEEP.
C
C     11/7/76 - FIRST VERSION
C     4/9/78 - CHIGNVEC ADDED TO LIST - S.P.                               20   
C     4/25/78 - LI AND LO INTERCHANGED - S.P.
C     4/27/78 - I AND ISUM ADDED TO SPECIAL NAMES - S.P.
C     3/6/79 - S MATRICES - S.P.
C     12/28/79 - STANDARD MORTRAN - RPG
C     1/16/80 - CHARACTER*8 OBNAME, MOVE CIB%F DOWN - RPG
C     12/14/18 - CROSSSECTIONS, F, ETC BY ANGLE GRID - S.P.
C     3/19/04 - speakez IV keep call - s.p.
C
      implicit real*8 ( a-h, o-z )                                      implicit
C                                                                          30   
      parameter ( namdim = 500 )                                        /allocs/
c
      INTEGER*4 Z
      COMMON /LOCptrs/ Z(namdim)
      INTEGER  FACFR4, FACFR2, FACFR1
      COMMON /ALLOCS/ FACFR4, FACFR2, FACFR1
      DIMENSION  ALLOC(1)
      character*8 dalloc(1)
      character*1  chlloc(1)
      REAL*4  ALLOC4(1)
      INTEGER    ILLOC(1)
      INTEGER*2 ILLOC2(1)
      EQUIVALENCE ( facfr4, ALLOC(1), ILLOC(1), ILLOC2(1),
     &   ALLOC4(1), chlloc(1), dalloc(1) )
 
      COMMON /LENGTH/  LENG(1)
C
      parameter ( numfloat = 153 )                                      /float/
      COMMON /FLOAT/  A, ACCURA, AI, ANGSTP, ANGMIN,
     &  ANGMAX, ASO, ASOI, ASYMPT, DELTVK,
     &  DWCUT, E, ECM, ELAB, EXS(5), GAMDIF,
     &  GAMSUM, TAU, TAUI, AM, AMA, AMB, AMBIGA, AMBIGB, AMX,
     &  AMXCS(5), AMXCGS(5), AMP,
     &  AMT, Q, R, R0, RC,
     &  RC0, RI, RI0, RSO, RSO0,
     &  RSOI, RSOI0, SPAM, SPAMP, SPAMT,
     &  STEPSZ, STEP1R, STEP1I, SUMMAX, SUMMID,
     &  SUMMIN, V, VI, VSO, VSOI,
     &  FITMUL, WRITES, ASI, RSI, RSI0,
     &  VSI, AMDMLT, STPSPR, R0E, RI0E,
     &  AE, AIE, VE, VIE, FITACC,
     &  ALMNMT, SPFACP, SPFACT, R0ESQ, RI0ESQ,
     &  AESQ, AIESQ, VESQ, VIESQ, DERIVS,
     &  FITRAT, RCP, RCT, RC0P, RC0T,
     &  SUMPTS, ACCINE, COULML, ALMXMT, EXSPT(2),
     &  POWRL, POWIM, BNDASY, SCTASY,
     &  VTR,VTRI,VTL,VTLI,VTP,VTPI,
     &  PARAM1, PARAM2, PARAM3, PARAM4, PARAM5,
     &  PAR620(15),
     &  RTR, RTRI, RTL, RTLI, RTP, RTPI,
     &  RTR0, RTRI0, RTL0, RTLI0, RTP0, RTPI0,
     &  ATR, ATRI, ATL, ATLI, ATP, ATPI,
     &  AMXGPT(2), PHIMID, GAMPHI
c
      DIMENSION PARAMS(20)
      EQUIVALENCE  ( PARAMS(1), PARAM1 )
      DIMENSION AMS(5)
      EQUIVALENCE ( AMS(1), AMA )
      DIMENSION VTEN(6), RTEN(6), RTEN0(6), ATEN(6)
      EQUIVALENCE ( VTEN(1), VTR ), (RTEN(1), RTR ),
     &  (RTEN0(1), RTR0 ), (ATEN(1), ATR )
 
      CHARACTER*1  REACT, HEADER                                        /hedcom/
      COMMON /HEDCOM/  REACT(45), HEADER(65)
 
      LOGICAL ONELSW, IDENSW                                            /inelcm/
      COMMON /INELCM/  JBPF, JBTF, LXMIN, LXMAX, ILILOR, ILILOI, MSTOP,
     &  IA12M, NMLOLX, ONELSW, NUMLX, IRDINT, IELCUP, LILOSZ, NUMLIS,
     &  ILIS, IBETAS(3), ICL2FF, LXSTEP, IDENSW, LSKIP,
     &  JPMIN, JPMAX, JPBASE, NJP, JTMIN, JTMAX, JTBASE, NJT,
     &  IINDXS, ISMATR, ISMATI, NSMAT, NLX, NSPLI, ITOCS,
     &  NASPLI, IUNITR, LIFIT, ILIFIT, NOLFIT, IDELSR, IDELSI
      INTEGER INDXJS(4,2)
      EQUIVALENCE ( INDXJS(1,1), JPMIN )
 
C
C
      CHARACTER*1  BLANK /' '/
      CHARACTER*8  OBNAME, MEMBER, MYKEEP / 'MYKEEP' /,  AS /'AS'/         40   
C
      DATA NSPECL / 16 /
      CHARACTER*8  SPECIL(16) / 'SMAG', 'SPHASE', 'LXCROSSS',
     1  'MXCROSSS', 'SIN', 'SOUT', 'F', 'HEADER', 'REACTION', 'ELAB',
     2  'ANGLEGRI', 'CHIGNVEC', 'S', 'ISUM',  'SEL', 'CROSSSEC' /
      INTEGER KINDS(16)  / 2, 2, 2, 2, 4, 4, 4, 9, 9, 2, 2, 2, 4, 4,
     1  4, 2 /
      INTEGER  ifamS(16) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2, 1, 1,
     1  1, 1 /
      INTEGER  ndimS(16) / 2, 2, 2, 2, 1, 1, 2, 1, 1, 0, 1, 2, 2, 2,       50   
     1  1, 2 /
      INTEGER ISIZES(16) / 8, 8, 8, 8,16,16,16, 1, 1, 8, 8, 8,16,16,
     1 16, 8 /
      INTEGER IGOTOS(16) / 2, 2, 3, 3, 1, 1, 3, 4, 5, 6, 7, 8, 2, 2,
     1  1, 3 /
 
      integer dims(15)
C
      CALL NXWORD ( OBNAME, *900, *900, *900 )
 50   CALL NXWORD ( MEMBER, *900, *900, *900 )                             60   
      IF ( MEMBER .EQ. AS )  GO TO 50
      do i = 1, 8
         ich = ichar(member(i:i))
         if ( ich >= ichar('A') .and. ich <= ichar('Z') )
     &     member(i:i) = char( ich - (ichar('A')-ichar('a')) )
      enddo
C
C     IS IT IN THE ALLOCATOR
C
      IOB = NAMLOC(OBNAME)                                                 70   
      LOCOB = Z(IOB)
      LENOB = LENG(IOB)
C
C     SET UP KIND, ETC
C
      KIND = 2
      KLASS = 5
      ISIZE = 8
      NCOLS = 1
      NROWS = 0                                                            80   
C
C     SEE IF WE KNOW SOMETHING SPECIAL ABOUT IT
C
      DO 189  I = 1, NSPECL
         IF ( SPECIL(I) .EQ. OBNAME )  GO TO 200
 189  CONTINUE
      GO TO 500
C
C     PROCESS SPECIAL STUFF
C                                                                          90   
 200  KIND = KINDS(I)
      ifam  =  ifamS(I)
      ndim = ndims(i)
      ISIZE = ISIZES(I)
      IGOTO = IGOTOS(I)
      GO TO ( 500, 220, 230, 240, 250, 260, 270, 280 ), IGOTO
C
C
 220  NCOLS = NMLOLX
      GO TO 500                                                           100   
C
C     SOMETHING ON THE ANGLE GRID
C
 230  NCOLS = (ANGMAX-ANGMIN)/ANGSTP + 1.5
      IF ( 8*LENOB/ISIZE .EQ. NROWS )  ndim = 1
      GO TO 500
C
C     HEADER
C
 240  DO 244  I = 1, 65                                                   110   
         IF ( HEADER(I) .NE. BLANK )  GO TO 245
 244  CONTINUE
      I = 65
 245  CALL KEEP ( HEADER(I), 9, 1, 1, 66-I, ' ', MEMBER, IER )
      NROWS = 66-I
      GO TO 265
C
C     REACTION
C
 250  DO 254  I = 1, 45                                                   120   
         IF ( REACT(46-I) .NE. BLANK )  GO TO 255
 254  CONTINUE
      I = 45
 255  CALL KEEP ( REACT, 9, 1 1, 46-I, ' ', MEMBER, IER )
      NROWS = 46-I
      GO TO 265
C
C     ELAB
C
 260  CALL KEEP ( ELAB, 2, 0, 0, 1, ' ', MEMBER, IER )                    130   
      NROWS = 1
 265  NCOLS = 1
      GO TO 800
C
C     ANGLEGRID
C
 270  CALL KEEP ( ANGSTP, 2, 1, 1, 3, ' ', MEMBER, IER )
      NROWS = 3
      GO TO 265
C                                                                         140   
C     OBJECT IS A SQUARE MATRIX OR ARRAY
C
 280  NCOLS = DSQRT( (8.*LENOB)/ISIZE + .010D0 )
      GO TO 500
C
C     OUTPUT AN OBJECT IN THE ALLOCATOR
C
 500  IF ( IOB .EQ. 0 )  GO TO 910
      LENOB = (8*LENOB)/ISIZE
      IF ( NROWS .NE. 0 )  NCOLS = LENOB/NROWS                            150   
      IF ( NCOLS .NE. 0 )  NROWS = LENOB/NCOLS
      LENOB2 = NROWS*NCOLS
      IF ( LENOB .EQ. LENOB2 )  GO TO 600
      write (6, 513) OBNAME, LENOB, NROWS,
     1  NCOLS, LENOB2, KIND, ifam, ndim, dims(1:2)
 513  FORMAT ( '0**** WARNING: ', A8, ' HAS', I6,
     1  ' ELEMENTS BUT NROWS, NCOLS & OUTPUT LENGHT =', 3I6,
     2  ' KIND, ifam, ndim =', 2I3, 2i10 )
      NROWS = LENOB
      NCOLS = 1                                                           160   
      ifam = 1
      ndim = 1
 
C
 600  if ( ndim == 2 .and. nrows == 1 ) then
         ndim = 1
         nrows = ncols
         ncols = 1
      endif
      dims(1) = nrows                                                     170   
      dims(2) = ncols
      CALL KEEP ( ALLOC(LOCOB), KIND, ifam, ndim, dims, ' ',
     1  MEMBER, IER )
C
 800  IF ( IER .NE. 0 )  GO TO 920
      write (6, 803) OBNAME, MEMBER, KIND, ifam, ndim, dims(1:2)
 803  FORMAT ( 1X, A8, ' HAS BEEN KEPT AS MEMBER ', A8,
     1  5X, 'KIND =', I2, 3X, 'family =', I2,
     2  3X, 'num. dims =', I2, 3X, 'dims =', 2I6 )
      RETURN                                                              180   
C
C     ERROR MESSAGES
C
 900  write (6, 903)
 903  FORMAT ( '0**** KEEP MUST BE FOLLOWED BY TWO NAMES.' )
      GO TO 980
C
 910  write (6, 913) OBNAME
 913  FORMAT ( '0**** ', A8, ' IS NOT A DEFINED OBJECT.' )
      RETURN                                                              190   
 920  write (6, 923) IER, OBNAME, MEMBER
 923  FORMAT ( '0**** ERROR NUMBER ', Z8,
     1  ' ENCOUNTERED WHILE KEEPING ', A8, ' AS ', A8 )
      RETURN
C
 980  write (6, 983)
 983  FORMAT ( ' **** CONTINUING WITH NEXT INPUT CARD.' )
      CALL NEWCD
      RETURN
C                                                                         200   
      END
