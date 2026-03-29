C*ID* GETLNK   C109.PTOLEMY.SOURCE                          PER701  13:
      SUBROUTINE  GETLNK ( AKEY, IRET )
C
C     PROCESSES LINKULE KEYWORDS AND LOADS LINKULE
C
C     THIS ROUTINE DETERMINES WHAT ( REAL POTENTIAL, IMAG POTENTIAL,
C     WAVE FUNCTIONS) IS TO BE SET BY THE LINKULE AND ALSO THE
C     LINKULE NAME.
C     IT THEN LOADS THE LINKULE AND SAVES ITS ADDRESS.
C                                                                          10   
C     IRET - 0 = AN ERROR; 1 = O.K.
C
C     11/24/77 - FIRST VERSION - S.P.
C     7/6/78 - CDC VERSION ADDED
C     12/7/78 - TENSOR POTENTIALS ADDED, AND USE NUMLNK - RPG
C     12/27/79 - CHANGE CDC PREFIXES TO CNI, ADD MORE NAMES - RPG
C     1/16/80 - MORE CHARACTER*8'S, REMOVE EXTRA / - RPG
C      1/21/80 - %M, CIB%F AROUND LOCF CALLS - RPG
C     2/25/80 - REMOVE 1ST INDEX FROM LNKNAM - RPG
C     3/14/80 - REMOVE BOTH ABOVE CHANGES - S.P.                           20   
C     7/28/81 - PRINT LOCATION OF LOAD - S.P.
C     6/15/83 - ADD SOME MORE TO THE CDC VERSION - S.P.
C     6/16/83 - GIVE UP ON NOS AND DO NOT USE CALLER - S.P.
C     12/19/84 - VSI MAY BE DONE BY A LINKULE - S.P.
C     4/12/85 - FIXUP CNICND RAWITSCH CARD TO TWO LINES - S.P.
C     5/17/85 - CLEANUP, ALL NON-IBM HAVE BUILT IN LINKULES - S.P.
C     3/15/03 - make CNI lines default for ease of reading;
C               If caller ever reinstated, this must be changed back - s
C               Add AV18
C     2/15/07 - add PHIFFER                                                30   
 
C
      implicit real*8 ( a-h, o-z )                                      implicit
c                                                                       /lnkblk/
c  STRUCTURE OF   /LNKBLK/
c
c  THERE IS A BLOCK OF 3 DOUBLE WORDS FOR EACH POSSIBLE USE
c  OF A LINKULE.  THE NUMLNK = NUMLINKULES USES ARE IN THE ORDER
c    REALPOTEN  IMAGPOTEN  REALVSO  IMAGVSO  COULOMB
c    WAVEFUNCTION  AND 6 TENSOR POTENTIALS
c  THE THREE-DOUBLE-WORD BLOCKS EACH HAVE THE STRUCTURE
c
c    LNKNAM(1,I)   LINKULE NAME
c
c      LNKADR(3,I)   LINKULE ADDRESS, O IF NO LINKULE IN USE.
c      LNKADR(4,I)   0 - NO L OR J DEPENDENCE IN POTENTIAL
c                    1 - L DEPENDENCE BUT NO J DEPENDENCE
c                    2 - BOTH L  AND J DEPENDENCE
c
c      LNKADR(5&6,I)  MYINTS - A TWO-ELEMENT ARRAY THAT THE
c                     LINKULE MAY USE FOR COMMUNICATION BETWEEN
c                     CALLS TO ITSELF.
c
      parameter ( NUMLINKULES = 13 )
      COMMON /LNKBLK/  NUMLNK, IUNIQU, LNKADR(6,NUMLINKULES)
      character*8  LNKNAM(3,NUMLINKULES)
      EQUIVALENCE ( LNKNAM(1,1), LNKADR(1,1) )
 
C
      CHARACTER*8  AKEY, WORD
      CHARACTER*8 AKEYS(NUMLINKULES) /
     1  'REALPOTE', 'IMAGPOTE', 'REALSOPO',
     1  'IMAGSOPO', 'COULOMBP', 'WAVEFUNC', 'REALTRPO', 'IMAGTRPO',
     2  'REALTLPO', 'IMAGTLPO', 'REALTPRP', 'IMAGTPRP', 'SIPOTENT' /       40   
C
C     FOR NON-IBM WE NEED A LIST OF SUBROUTINE (=LINKUL) NAMES
C
      parameter ( NUMNAM = 17 )
      CHARACTER*8  NAMES(numnam) / 'BKGPTELP', 'FIXEDWOO', 'GAUSSIAN',
     1   'LAGRANGE', 'LTSTELP ', 'RAWITSCH',
     2   'REID', 'SHAPE', 'SPLINE', 'TWOSHAPE', 'DEFORMED',
     3   'JDEPEN', 'JDEPENWS', 'OHTA', 'PARITWOO', 'AV18', 'PHIFFER' /
C
CIB      CHARACTER*8  LINKDD /'LINKULES'/                                  50   
      CHARACTER*8  INTERN /'INTERNAL'/
CIB      CHARACTER*8  FT10 /'FT10F001'/
C
c
C
C     FIND THE KEY
C
      DO 19  IKEY = 1, NUMLNK
         IF ( AKEY .EQ. AKEYS(IKEY) )  GO TO 50
 19   CONTINUE                                                             60   
C
C     GET THE LINKULE NAME
C
 50   CALL NXWORD ( WORD, *900, *900, *900 )
      LNKNAM(1,IKEY) = WORD
      ILOC = 0
      IF ( WORD .EQ. INTERN )  GO TO 200
C
C
C     LOAD IT                                                              70   
C
CIB      IF ( LOOKDD(LINKDD) .NE. 0 )  GO TO 125
CIB      write (6, 103) LINKDD
CIB 103  FORMAT ( '0**** ', A8, ' DD CARD MISSING.')
CIB      GO TO 950
CIB 125  IF ( LOOKDD(FT10) .NE. 0 ) GO TO 150
CIB      write (6, 103) FT10
CIB      GO TO 950
CIBC
CIB 150  ILOC = LOAD( 0, 0, LINKDD, WORD )                                 80   
CIB      IF ( ILOC .NE. 0 )  GO TO 200
CIB      write (6, 153) WORD
CIB 153  FORMAT ( '0**** COULD NOT FIND THE LINKULE ', A8 )
CIB      GO TO 950
C
C
C     IN THE NON-IBM VERSIONS THE LINKULES ARE IN THE OVERLAY
C
      DO 109  I = 1, NUMNAM
         IF ( WORD .EQ. NAMES(I) )  GO TO 140                              90   
 109  CONTINUE
      write (6, 113) WORD
 113  FORMAT ( '0*** ', A8, ' IS NOT IN THIS VERSION OF PTOLEMY.' )
      GO TO 950
C
 140  ILOC = I
C
C
 200  LNKADR(3,IKEY) = ILOC
      write (6, 203) WORD, ILOC                                           100   
 203  FORMAT ( ' LINKULE ', A8, ' IS LINKULE NUMBER/location', I8 )
      IRET = 1
      RETURN
C
 900  write (6, 903) AKEY
 903  FORMAT ( '0**** A LINKULE NAME MUST FOLLOW THE ', A8,
     1  ' KEYWORD.' )
 950  IRET = 0
      RETURN
      END                                                                 110   
C*ID* LINKUL   C109.PTOLEMY.SOURCE                          PER701  13:
      SUBROUTINE LINKUL ( LOC, ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2, ID )
C
C     LINKULE CALLER FOR PTOLEMY
C
C     THIS ROUTINE ADDS THE FIXED PART OF A LINKULE ARGUMENT LIST
C     TO THE ABOVE ARGUMENTS AND THEN CALLS THE LINKULE.
C     THE LINKULE MUST ALREADY BE LOADED IN CORE AND LOC MUST BE ITS
C     ENTRY POINT ADDRESS.                                                120   
C     UPON RETURN, THE STUFF WRITTEN IN IUNIT WILL BE COPIED TO THE
C     PRINTER FOR THE IBM VERSION.  ON OTHER VERSIONS IT IS DIRECTLY
C     PRINTED THERE.
C
C     ON NON-IBM, LINKULES ARE PART OF THE PROCESSOR AND WE JUST
C     CALL THEM DIRECTLY.
C
C     11/25/77 - FIRST VERSION - S.P.
C     7/6/78 - CDC VERSION
C     5/6/79 - PASS IPARAM ARRAY SEPARATELY - S.P.                        130   
C     12/28/79 - STANDARD MORTRAN - RPG
C     7/30/80 - IUNIT IS 6 FOR CDC - S.P.
C     6/16/83 - CALLER FAILS SOMEHOW ON NOS - S.P
C     6/24/83 - PERHAPS IT IS REALLY EXTERNAL THAT FAILS - S.P.
C     12/19/84 - PUT IN DIAGNOSITC PRINTOUT - S.P.
C     5/17/85 - SIMPLIFY NON-IBM VERSION; GIVEUP ON CALLER - S.P.
C     3/15/03 - make CNI lines default for ease of reading;
C               If caller ever reinstated, this must be changed back - s
C               Add AV18
C     2/15/07 - add PHIFFER                                               140   
C
      implicit real*8 ( a-h, o-z )                                      implicit
C
      DIMENSION  MYINTS(2), ARRAY1(1), ARRAY2(1)
      character*1 id(4)
C
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
 
      COMMON / LENGTH /  LENG(1)
C
      COMMON / FLOAT / F(116)                                             150   
      COMMON / INTGER / I(46)
      EQUIVALENCE ( IPRINT, I(17) )
      COMMON / JBLOCK / JB(1)
      COMMON / SWITCH / IS(1)
      COMMON / CNSTNT / C(1)
      COMMON / TEMPVS / T(1)
      COMMON / WAVCOM / WAV(1)
      COMMON / KANDM / KM(1)
      COMMON / INTRNL / IN(1)
C                                                                         160   
      LOGICAL DBUGSW
      EXTERNAL  NALLOC, NAMLOC
C
CIB      DATA IUNIT /10/
      DATA IUNIT /6/
C
      CHARACTER*1  LINE(133)
C
c
C                                                                         170   
      DBUGSW = MOD(IPRINT, 10) .EQ. 9
C
      IF (DBUGSW)  write (6, 23) LOC, ALIAS, MYINTS, IPOTYP, IREQUE, ID
 23   FORMAT (' CALLING LINKULE AT', i8, 1X, A8, 2I8, 2I4,1X, 4a )
C     NOTE:
C
C     F(112) IS BEGINING OF PARAM'S
C     I(46) IS BEGINNING OF IPARAM'S
C
CIB      CALL CALLER ( LOC, ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,        180   
CIB     1  IUNIT, NUMOUT,
CIB     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
CIB     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
CIB     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
CIB     5  I(46) )
C
      GO TO ( 101, 102, 103, 104, 105, 106,  107, 108, 109,
     1    110, 111, 112, 113, 114, 115, 116, 117 ), LOC
C
 101  CALL bkgptelp ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,              190   
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 102  CALL fixedwoo ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,                     200   
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 103  CALL gaussian ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,                210   
     5  I(46) )
      GO TO 190
C
 104  CALL lagrange ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190                                                           220   
C
 105  CALL ltstelp  ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 106  CALL RAWITSch ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,              230   
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 107  CALL REID ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,                     240   
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 108  CALL shape ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,                250   
     5  I(46) )
      GO TO 190
C
 109  CALL spline ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190                                                           260   
C
 110  CALL twoshape ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 111   CALL DEFORMEd ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,             270   
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 112  CALL jdepen   ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,                     280   
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 113  CALL jdepenws ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,                290   
     5  I(46) )
      GO TO 190
C
 114  CALL ohta     ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190                                                           300   
C
 115  CALL paritwoo ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 116  CALL av18 ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,                  310   
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 117  CALL phiffer ( ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
     1  IUNIT, NUMOUT,
     2  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,                     320   
     3  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
     4  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
     5  I(46) )
      GO TO 190
C
 190  CONTINUE
      IF (DBUGSW)  write (6, 193) IRET, MYINTS, IRETUR, NUMOUT
 193  FORMAT ( ' LINKULE RETURNS', 5I10 )
CIB      IF ( NUMOUT .EQ. 0 )  RETURN
CIBC                                                                      330   
CIB 300  READ ( IUNIT, 303, END=400 )  LINE
CIB 303  FORMAT ( 133A1 )
CIB      write (6, 303) LINE
CIB      GO TO 300
CIB 400  REWIND IUNIT
      RETURN
 900  write (6, 903) LOC
 903  FORMAT ( '0**** LINKULE', I10, ' NOT AVAILABLE.' )
      STOP 1122
      END                                                                 340   
