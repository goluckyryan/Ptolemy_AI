C*ID* IALLOC   C109.PTOLEMY.SOURCE                          PER701  09:
      FUNCTION IALLOC(IWRDS)
C
C     THIS IS THE ENTIRE NUMBERED STORAGE PACKAGE
C
C-----------------------------------------------------------------------
C  THIS ROUTINE ASSIGNS SPACE IN THE ALLOCATED REGION WHICH BEGINS
C  AT ALLOC(IBASE) AND EXTENDS TO ALLOC(IBASE+NWORDS)
C
C  THE ALLOCATED ARRAY IS LABELED BY A NUMBER- RETURNED AS THE             10   
C  FUNCTIONAL VALUE FROM IALLOC
C
C  SPACE IS ASSIGNED TO UP TO NAMDIM ARRAYS.  EACH ALLOCATION ICREMENTS
C  THE SCALAR IAM BY 1.  THE LOCATION OF THE FIRST BYTE OF THE ARRAY
C  IS KEPT IN L(IAM) AND ITS LENGTH (IWRDS) IS SAVED IN LENG(IAM).
C  AN ARRAY IS FREED FOR FUTURE RE-USE BY THE OPERATION
C                L(NUMBER)=-L(NUMBER)
C  WHERE NUMBER IS THE IDENTIFYING ARRAY NUMBER (FUNCTION VALUE OF
C  IALLOC) AND L IS THE INTEGER*4 ARRAY IN THE LABELED COMMON LOC
C                                                                          20   
C
C  IF ALLOCATION IS IMPOSSIBLE IALLOC RETURNS A FUNCTION VALUE 0
C     UNLESS IF   STOPSW = .TRUE.  (SEE BELOW)  INWHICH CASE THE JOB IS
C     TERMINATED WITH AN ERROR TRACE.
C
C  TO CHANGE THE MAXIMUM NUMBER OF NAMES RE-COMPILE WITH
C  COMMON AREAS NAMES,LOC AND LENGTH WITH NEW DIMENSIONS
C  AND WITH THE NEW SIZE IN THE NAMDIM DATA STATEMENT BELOW.
C
C-----------------------------------------------------------------------   30   
C
C
C  THIS VERSION 7/14/73    D.GLOECKNER
C     2/3/75 - STOPSW ADDED, NEW STATISTIC IN NSCAT; NSSTAT ENTRY - S.P.
C              ALSO ADD SOME STUFF FOR NAMED-NUMBERED-STORAGE -- ALL
C              SUCH ADDITIONS CONTAIN A REFERENCE TO "NAMES" OR "BLANK"
C     2/16/75 - MANY CHANGES - S.P.
C     2/24/75 - FIX BUG IN COMPRESSING (IGOT UNDEFINED) S.P.
C     8/15/75 - ADD GIVEAL AND CHOPIT ENTRIES D.G.
C     11/7/75 - MAKE REFRESHABLE;  STOP FOR =< 0 SIZE.                     40   
C     6/6/76 - WORD SIZE FACTORS IN ALLOCS
C     9/10/76 - ADD IREDEF
C     12/28/76 - FIX BUG IN ABOVE
C     4/17/77 - ALLOW OVERLAYS BETWEEN GIVEAL AND CHOPIT
C     5/10/77 - FIX GIVEAL BUG WHEN NAMESTACK IS FULL
C     6/11/77 - FIX COMPRESS BUG WHEN IT IS ALREADY COMPRESSED.
C     1/25/78 - FIX BUG FOR EXACT FIT COMPRESS.
C     4/25/78 - NO ENTRIES FOR CRAY VERSION - S.P.
C     4/26/78 - PRINT DIAGNOSTICS BASED ON NSPTSW - S.P.
C     7/3/78 - AVOID CDC OVERFLOW ON INTEGRETIY HCHECK                     50   
C     7/24/79 - NEW NSDUMP MADE A SEPARATE SUBROUTINE - S.P.
C     11/4/79 - INCREASE TO 100 NAMES - S.P.
C     12/31/79 - STANDARD MORTRAN; CVA, CUN VERSIONS - RPG
C     1/24/80 - CHARACTER*8 THENAM; LIB%STOP FOR CVA; CORRECT URWALK
C        FOR CUN; 9950 TO 9550; NAMLCA (NOT NAMLOC) ARG TO IREDEF - RPG
C     6/15/80 - CORRECT BOB'S IREDEF FIX FOR CDC - S.P.
C     7/14/80 - CRAY #2 VERSION -- MUCH CLOSER TO IBM FORTRAN - S.P.
C     7/23/80 - CRAY VERSION CAN ALOS HAVE EXPANDING ALLOCATOR - S.P.
C     7/28/80 - ALOW COMPRESS LOOP TO VECTORIZE - S.P.
C     8/20/85 - ALLOW WHOLE THING TO MOVE WHEN EXPANDING - S.P.            60   
C
      implicit real*8 ( a-h, o-z )                                      implicit
      parameter ( NUMSWITCHES = 34 )                                    /switch/
      INTEGER  PROBLM, R0TYPE, CONVRG
      COMMON /SWITCH/  IASYMP, ICHECK, IECHO, IELAST, IFIT,
     &  IVRTEX, IRLWAV, ITSO, LINEAR, NUCONL,
     &  NEXT, IFITER, IEPOW, ISAVMX, ISAVB,
     &  LABANG, PROBLM, MASTYP, R0TYPE, IBSPAS,
     &  ISW1, ISW2, ISW3, ISCRFL, NOBFAC,
     &  CONVRG, IEXTYP, IBRNSB, MEORD, KFRMTP,
     &  IWRTWV, LORNTZ, IRELPT, ISAVSM
 
C
      INTEGER GIVEAL, CHOPIT                                            cnd
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
 
C
      COMMON /LENGTH/  LENG(namdim)
      CHARACTER*8  NAMES                                                   70   
      COMMON /NAMCOM/  NAMES(namdim)
C
C
C     FOLLOWING COMMON CONTROLS NUMBERED STORAGE
C
C     IBASE - OFFSET (IN 8-BYTE WORDS) TO START OF ALLOCATOR
C     NWORDS - NUMBER OF 8-BYTE WORDS ASSIGNED TO ALLOCATOR.
C     IMAX - SIZE OF NAME STACK
C     INAME - HIGHEST LOCATION REACHED IN NAME STACK
C     IGOT - CURRENT NUMBER OF "ACTIVE" SLOTS IN NAME STACK                80   
C     NLOC - FIRST FREE WORD AT END OF ALLOCATOR
C     IPEAK - PEAK ACTUALLY USED SPACE IN UNITS OF 8-BYTES.
C     NAMPEK - PEAK NUMBER OF NAMES ACTUALLY IN USE.
C     NUMCOM - NUMBER OF COMPRESSES THAT HAVE OCCURED.
C     ICHEC - ADDRESS (IN 8-BYTE WORDS) OF END OF ALLOCATOR.
C     THENAM - NAME OF OBJECT BEING ALLOCATED (SET BY NALLOC).
C     IWASAT - GIVEAL OBJECT LOCATION
C     NAMLC2 - GIVEAL NAME STACK LOCATION
C
      LOGICAL NSPTSW                                                    /usage/
      COMMON/USAGE/IBASE, NWORDS, IMAX, INAME, IGOT, NLOC,
     & IPEAK, NAMPEK, NUMCOM, ICHEC, THENAM, IWASAT, NAMLC2,
     & NSPTSW
 
 
C
C
C     STRUCTURE OF EACH ITEM IN THE ALLOCATOR IS
C
C     WORD         IBM ETC.           CDC ETC.
C       1  NNNN 0000 FFFF 0000    00000 NNNNN
C     2,3..        DATA               DATA
C
C             WORD = 8 BYTES       WORD = 60 BITS
C     IN BOTH CASES NNNN = NUMBER OF OBJECT (USED FOR VALIDITY            100   
C     CHECKS).  IN IBM THE HEXADECIMAL STRUCTURE IS SHOWN;
C     IN THE CDC THE NUMBER IS JUST A FULL-WORD INTEGER.
C
CCCcni      character*8  ANAME / Z 0000 0000 FFFF 0000 /
      character*8  ANAME                                                cni
      integer*8 ianame                                                  cni
      equivalence ( aname, ianame )                                     cni
      data  iANAME / Z'00000000FFFF0000' /                              cni
CIB      character*8  ANAME / Z 0000 0000 FFFF 0000 /
CVA      REAL*8  ANAME / '00000000FFFF0000'X /                            110   
      INTEGER*2 NUMBER(4)
      EQUIVALENCE(NUMBER(1),ANAME)
C
      LOGICAL  STOPSW /.TRUE./, REDFSW
      LOGICAL ERRSW,  PRTSW, FONDSW, ALLSW
      CHARACTER*8  BLANK /' '/, THENAM
C
c
C
C                                                                         120   
      REDFSW = .FALSE.
      IF ( NSPTSW )  write (6, 3) THENAM, IWRDS
  3   FORMAT ( ' IALLOC/NALLOC: ', A8, I15 )
C
  5   IF(IWRDS.GT.0)GO TO 100
      write (6, 13) THENAM, IWRDS
 13   FORMAT('0**** ', A8, '  - REQUEST FOR ', I15, ' WORDS. ****')
      GO TO 9500
C
C                                                                         130   
C     PREPARE TO ALLOCATE AN OBJECT.
C
C     FIRST GO THROUGH THE LIST AND FIND USED SPACE AND NUMBERS.
C     DURING THIS PASS WE ACCUMULATE INFO FOR STATISTICAL PURPOSES
C     AND CHECK THE INTEGRETY OF THE ALLOCATOR.  ALSO WE LOOK FOR AN
C     EXACT FIT TO A FREED SPOT WHICH WILL BE USED IF FOUND.
C
C     WE ALSO LOOK FOR AN AVAILABLE NAME SPOT.
C
 100  FONDSW = .FALSE.                                                    140   
      NAMUSE = 0
      LUSE = 0
      NAMLOC = 0
      IF ( INAME .EQ. 0 )  GO TO 150
      DO 149  N = 1, INAME
C
C     IS THIS AN EMPTY SLOT
C
         IF (z(N) .NE. 0)  GO TO 120
         IF ( NAMLOC .EQ. 0 )  NAMLOC = N                                 150   
         GO TO 149
C
C     IT IS IN USE OR IS A FREED SPOT.  CHECK INTEGRETY
C
 120     NUMBER(1) = N
      IF ( dALLOC(IABS(z(N))-1) .EQ. ANAME )  GO TO 130
         write (6, 123) N, z(N), LENG(N), NAMES(N),
     1      z(N-1), LENG(N-1), NAMES(N-1)
 123  FORMAT ('0*** HEADER OF OBJECT', I4, ' HAS BEEN DESTROYED.' /
     1   ' THIS AND PREVIOUS OBJECT LOCATION, LENGTH AND NAME ARE:' /     160   
     2   ( 2I15, 2X, A8 ) )
C
C     IS IT FREE SPACE OR A STILL IN USE OBJECT
C
 130     IF ( z(N) .LT. 0 )  GO TO 140
         NAMUSE = NAMUSE + 1
         LUSE = LUSE + LENG(N) + 1
         GO TO 149
C
 140     IF (FONDSW)  GO TO 149                                           170   
         IF ( IWRDS .NE. LENG(N) )  GO TO 149
         NAMLOC = N
         ILOC = -z(NAMLOC) - 1
         FONDSW = .TRUE.
C
 149  CONTINUE
C
C     SOME STATISTICS
C
 150  IPEAK = MAX0( IPEAK, LUSE+IWRDS+1 )                                 180   
      NAMPEK = MAX0( NAMPEK, NAMUSE+1 )
C
C     DID WE FIND AN EXACT FIT
C
      IF ( FONDSW )  GO TO 250
C
C     NO, DID WE ATLEAST FIND A NAME SLOT
C
      IF ( NAMLOC .NE. 0 )  GO TO 200
C                                                                         190   
C     NO, WE MUST GET A NEW NAME SLOT IF THERE ARE MORE AVAILABLE
C
      IF(INAME.LT.IMAX)GO TO 170
C
C     WE HAVE GONE THROUGH THE WHOLE STACK OF NAMES, ARE THERE HOLES
C
      IF ( NAMUSE .LT. IMAX )  GO TO 200
 160  write (6, 163) THENAM
 163  FORMAT ('0*** ', A8, 2X, '-- NO MORE NAME SLOTS AVAILABLE ***')
      GO TO 9500                                                          200   
C
C     ADD TO END OF NAME LIST
C
 170  INAME = INAME + 1
      NAMLOC = INAME
C
C     IS THERE SPACE AT THE END OF THE ALLOCATOR
C
 200  IF ( NLOC+IWRDS .LE. ICHEC )  GO TO 220
C                                                                         210   
C     NO, WILL COMPRESSION SOLVE IT
C
      IF ( LUSE+IWRDS .LT. NWORDS )  GO TO 500
crayC
crayC     ON CDC & CRAY ATTEMPT TO EXPAND THE ALLOCATOR
crayC
cray      NNEED = LUSE + IWRDS + 1
cray      NASK  = MAX0( NNEED, LUSE+10000 )
cray      IBOLD = IBASE
cray      CALL GRAB ( ALLOC, FACFR1, IBASE, NASK )                        220   
crayC
crayC     INCREASING THE SIZE MAY HAVE CAUSED THE ENTIRE MESS TO MOVE
crayC
cray      IF ( IBASE .EQ. IBOLD )  GO TO 209
cray      IDEL = IBASE - IBOLD
cray      NLOC = NLOC + IBOLD
cray      DO 207  N = 1, INAME
cray          IF ( z(N) .GT. 0 )  z(N) = z(N) + IDEL
cray          IF ( z(N) .LT. 0 )  z(N) = z(N) - IDEL
cray 207  CONTINUE                                                        230   
cray 209  IF ( NASK .LT. NNEED ) GO TO 210
cray      NWORDS = NASK
cray      ICHEC = IBASE + NWORDS - 1
cray      GO TO 100
C
 210  write (6, 213) THENAM, IWRDS, IPEAK
213   FORMAT ('0*** ', A8, 2X, '-- NOT ENOUGH ALLOCATOR.  ATTEMPTING',
     1  ' TO DEFINE A', I12, ' WORD OBJECT.'  /
     2  5X, 'TOTAL ALLOCATOR SIZE MUST BE AT LEASET', I13, ' WORDS.' )
      GO TO 9500                                                          240   
C
C     THERE IS ROOM TO ALLOCATE AND ALSO A NAME SLOT IS AVAILABLE.
C     IF NAMLOC IS NOT YET ASSIGNED THEN WE MUST COMPRESS TO GET THE
C     SLOT.
C
 220  IF ( NAMLOC .EQ. 0 )  GO TO 500
      ILOC = NLOC
      NLOC = NLOC + IWRDS + 1
      IGOT = IGOT + 1
C                                                                         250   
C     DEFINE THE OBJECT
C
 250  IF ( NSPTSW )  write (6, 253) NAMLOC, ILOC, NLOC
 253  FORMAT ( ' ASSIGNED', I4, I12, 5X, 'NLOC =', I12 )
C
      z(NAMLOC)=ILOC+1
C
C  PUT NAME (NUMBER) IN FIRST INTEGER*2 PART OF ARRAY
C
      NUMBER(1)=NAMLOC                                                    260   
      dALLOC(ILOC)=ANAME
      LENG(NAMLOC)=IWRDS
      NAMES(NAMLOC) = BLANK
      IALLOC=NAMLOC
      NAMLC2 = NAMLOC
      IF ( REDFSW )  GO TO 620
      RETURN
C
      ENTRY GIVEAL(IWRDS)
C                                                                         270   
C
C     COMPRESS THE ALLOCATOR AND
C  GIVE LARGEST CONTIGUOUS PIECE FROM TOP OF ALLOCATOR TO END
C  RETURNS NUMBER OF WORDS YOU GOT
C  A CALL TO THIS ENTRY SHOULD BE IMMEDIATELY FOLLOWED WITH A CALL TO
C  CHOPIT
C
      IF ( NSPTSW )  write (6, 303)
 303  FORMAT ( ' GIVEAL CALLED' )
      REDFSW = .FALSE.                                                    280   
      ALLSW = .TRUE.
      IF ( NLOC .GT. IBASE )  GO TO 505
 300  IWRDS=ICHEC-NLOC-2
      IWASAT=NLOC
      IF ( IWRDS .GE. 2 )  GO TO 310
      WRITE(6,301)
301   FORMAT(' ***ERROR/NOTHING TO GIVE IN GIVEAL  CALL')
      GO TO 9500
C
C     FIND A NAME STACK LOCATION FOR GIVEAL                               290   
C
 310  IF ( INAME .EQ. 0 )  GO TO 170
      DO 304  N = 1, INAME
         IF ( z(N) .NE. 0 )  GO TO 304
         NAMLOC = N
         GO TO 200
 304  CONTINUE
      IF ( INAME .LT. IMAX )  GO TO 170
      GO TO 160
C                                                                         300   
      ENTRY CHOPIT(IWRDS)
C
C
C  THIS CALL MUST FOLLOW GIVEAL--IT TRUNCATES N.STORAGE OBJECT
C
C  PASS NUMBER OF WORDS YOU REALLY WANT
C  CHECK IF THIS CALL IMMEDIATELY FOLLOWS GIVEAL
C
      IF ( NSPTSW )  write (6, 323) IWRDS
 323  FORMAT ( ' CHOPIT CALLED:', I15 )                                   310   
      ANAME=dALLOC(IWASAT)
      IHAD=NUMBER(1)
      IPEAK=IPEAK+IWRDS+1
      IF(NAMLC2.NE.IHAD)GO TO 320
      NAMPEK=NAMPEK+1
      NLOC=IWASAT+IWRDS+1
      LENG(NAMLC2)=IWRDS
      RETURN
320   WRITE (6,321)
321   FORMAT(' ***ERROR CALL TO CHOPIT MUST IMMEDIATELY FOLLOW GIVEAL')   320   
      GO TO 9500
C
C
C
C     COMPRESSION IS NECESSARY AND IT WILL SOLVE OUR PROBLEMS
C     OR
C     GIVEAL HAS BEEN CALLED AND WE WANT TO GET THE MOST POSSIBLE
C
C
C  THIS SECTION FOR COMPRESSION                                           330   
C
C
 500  ALLSW = .FALSE.
C
      IF ( NAMLOC .NE. 0 )  WRITE(6,7432)  THENAM
7432  FORMAT( '0 ', A8, '  --  COMPRESSION IS OCCURING NOW.' )
      IF ( NAMLOC .EQ. 0 )  write (6, 7433) THENAM
7433  FORMAT( '0 ', A8, '  --  NAME STACK COMPRESSION IS OCCURING NOW.'
     1  )
C                                                                         340   
 505  LL=IBASE
      LCNT=0
      LSUM = 0
      NUMCOM = NUMCOM + 1
C
C     FIND FIRST FREE OBJECT
C
 510  LLOW = LL
      IF ( LCNT .GE. IGOT )  GO TO 600
C                                                                         350   
      ANAME=dALLOC(LL)
C
C  NEGATIVE LOCATION INDICATES A FREE LOCATION
C
C
C  STRUCTURE CHECK
      NAMLOC = NUMBER(1)
      IF(IABS(z(NAMLOC)).NE.(LL+1))GO TO 8171
      IF(z(NAMLOC).LT.0)GO TO 520
      LL=LL+LENG(NAMLOC)+1                                                360   
      LCNT=LCNT+1
      GO TO 510
8171  WRITE(6,7170)LL, ANAME, NAMLOC, LENG(NAMLOC), z(NAMLOC),
     1  NAMES(NAMLOC)
7170  FORMAT(' **** STRUCTURE CHECK: LOC, HEADER, NAMLOC, LENG, L,',
     1    ' NAME:',
     2  /   I13, Z18, 3I13, 2X, A8 )
      GO TO 9500
C
C     FIRST FREE LOCATION FOUND                                           370   
C
C
C  LOWEST FREE LOCATION IS AT LLOW-- CORRESPONDING TO THE ARRAY NUMBER
C  NAMLOC... THE LCNT(TH) ONE UP FROM THE BOTTOM OF THE ALLOCATOR
C
 520  IHCNT=LCNT
C
C
C     HAVE FOUND A FREE LOCATION.  MARK IT NOT IN USE AND GO ON.
C                                                                         380   
 530  z(NAMLOC)=0
      LSUM = LSUM + LENG(NAMLOC) + 1
C
C  HAVE FOUND FIRST FREE LOCATION--NOW FIND NEXT FILLED ONE
C
 540  IF(IHCNT.GE.(IGOT-1))GO TO 600
      IHCNT=IHCNT+1
      LL=LL+LENG(NAMLOC)+1
      IF(LL.GE.NLOC)GO TO 600
      ANAME=dALLOC(LL)                                                    390   
      NAMLOC = NUMBER(1)
C
C  INTERNAL STRUCTURE CHECK
C
      IF(IABS(z(NAMLOC)).NE.(LL+1))GO TO 8171
      IF(z(NAMLOC).LT.0)GO TO 530
C
C
C     FOUND A STILL-IN-USE OBJECT, COMPRESS
C                                                                         400   
      N=LENG(NAMLOC)+1
CRACDIR$ IVDEP
      DO 579 I=1,N
579   dALLOC(LLOW+I-1) = dALLOC(LL+I-1)
C
C  NOTE WE HAVE ALSO MOVED THE IDENTIFYING NUMBER
C
      LCNT=LCNT+1
      z(NAMLOC)=LLOW+1
      LLOW = LLOW+N                                                       410   
      GO TO 540
C
C     COMPRESS ALL DONE.  UPDATE NUMBER OF ACTIVE SLOTS AND TOP
C     OF ALLOCATOR "IN USE".
C
 600  IGOT=LCNT
      NLOC = LLOW
      IF ( ALLSW )  GO TO 300
C
C                                                                         420   
      WRITE(6,5101)IGOT, LSUM
5101  FORMAT(' ',I5,' ARRAYS COMPRESSED', I8, ' WORDS FREED' / )
      GO TO 100
C
C
C     CALL IREDEF ( NEWSIZE, NUMBER )   REDEFINES THE LENGTH
C     OF AN EXISTING OBJECT AND DUPLICATES THE EXISTING INFORMATION
C     INTO IT.  THE NUMBER AND NAME ARE NOT CHANGED BUT THE
C     LOCATION IN THE ALLOCATOR IS CHANGED.  SPACE MUST BE
C     AVAILABLE TO MAKE A TEMPORARY DUPLICATE COPY (EVEN FRO              430   
C     REDUCTIONS IN THE SIZE OF AN OBJECT).
C
C
      ENTRY IREDEF_F ( IWRDS, NAMLCA )
C
C
      NAMLOC = NAMLCA
      NAMOLD = NAMLOC
      THENAM = NAMES(NAMOLD)
      IF ( NSPTSW )  write (6, 613) THENAM, IWRDS                         440   
 613  FORMAT ( ' IREDEF: ', A8, I15 )
      IWDOLD = LENG(NAMOLD)
      REDFSW = .TRUE.
      GO TO 5
C
C     COPY OLD TO NEW INCLUDING HEADER
C
 620  N = MIN0( IWRDS, IWDOLD ) + 1
      LL = z(NAMOLD) - 1
      DO 629  I = 1, N                                                    450   
         dALLOC(ILOC-1+I) = dALLOC(LL-1+I)
 629  CONTINUE
C
C     CHANGE OLD TO POINT TO NEW
C
      z(NAMOLD) = ILOC+1
      LENG(NAMOLD) = IWRDS
      NAMES(NAMOLD) = THENAM
C
C     CHANGE NEW TO POINT TO OLD AND FREE IT                              460   
C
      dALLOC(LL) = ANAME
      LENG(NAMLOC) = IWDOLD
      z(NAMLOC) = -(LL+1)
      IALLOC = NAMOLD
      NAMLOC = NAMOLD
C
      RETURN
C
C                                                                         470   
C
C
      ENTRY ISTART_F (NBYTES)
CDC      NBYTES = IWRDS
C
C
C  INITIALIZE ALLOCATION
C
C     THE INPUT IS THE NUMBER OF BYTES TO ALLOCATE, OR IF NEGATIVE, THE
C     NUMBER OF BYTES TO LEAVE FREE.                                      480   
C
C
      NUMCOM = 0
      FACFR4 = 2
      FACFR2 = 4
      FACFR1 = 8
CRA      FACFR4 = 1
CRA      FACFR2 = 1
CRA      FACFR1 = 1
C                                                                         490   
      NWORDS = NBYTES
      IBASE = 0
      CALL GRAB(ALLOC,FACFR1,IBASE,NWORDS)
c
      if ( ibase .lt. 0 ) then
         write (0,*) ' ibase < 0:', ibase, nwords
!!!         write (6,*) ' ibase < 0 is not allowed!!!', ibase
!!!         stop 8765
      endif
c                                                                         500   
      ICHEC=IBASE+NWORDS-1
C
C  ALLOC(IBASE) IS FIRST BYTE OF ALLOCATED AREA--- ALLOC(IBASE+NWORDS)
C  IS LAST BYTE
C
      IMAX= NAMDIM
      IPEAK = 0
      NAMPEK = 0
      NUMCOM = 0
C                                                                         510   
 710  NLOC = IBASE
      INAME = 0
      IGOT = 0
      THENAM = BLANK
      DO 719 I=1, IMAX
      NAMES(I) = BLANK
719   z(I)=0
C
      IF ( NSPTSW )  write (6, 723) IBASE, ICHEC, IMAX
 723  FORMAT ( ' ALLOCATOR RANGE:', 2I12, 5X, '# NAMES =', I5 )           520   
      IALLOC=IBASE
      RETURN
C
C
C
C
      ENTRY ICLEAR_F (IDUM)
C
C
C  THIS ENTRY CLEARS ALL DATA                                             530   
C
      WRITE(6,8429)
 8429 FORMAT(/'   -------  ALL NUMBERED OBJECTS CLEARED  ------'/)
      GO TO 710
C
C     NSSTAT  PRINTS NUMBERED STORAGE STATISTICS
C
       ENTRY NSSTAT_F (IDUM)
      PRTSW = .FALSE.
      ALLSW = .FALSE.                                                     540   
      GO TO 900
C
C
C
C  NSCAT GIVES A LISTING OF NUMBERED OBJECT,THEIR LENGTHS AND LOCATIONS
C
C     IF  IDUM = -1  THEN FREE AREAS ARE ALSO LISTED
C
C
      ENTRY NSCAT_F (IDUM)                                                550   
C
CDC      IDUM = IWRDS
C
C
      PRTSW = .TRUE.
      ALLSW = IDUM .EQ. -1
 900  ERRSW = .FALSE.
C
 910  LUSE = 0
      NAMUSE = 0                                                          560   
C
      IF (PRTSW)  WRITE(6,923 )
923   FORMAT ( '             ---- NUMBERED STORAGE CATALOG ----' )
      IF ( NSPTSW .OR.  PRTSW .AND. ALLSW )  write (6, 924) IBASE, NLOC,
     1  ICHEC, IGOT
 924  FORMAT ( '0    FROM', I8, 3X, 'TO', I8, 6X,
     1  'TOP OF ALLOC =', I8, I8, ' ACTIVE NAME SLOTS' )
      IF ( PRTSW )  write (6, 927)
 927  FORMAT (  '0  NUMBER   NAME     LENGTH      START        END')
      DO  949 III=1,INAME                                                 570   
         LLOW = IABS(z(III))
         LL = LLOW + LENG(III)+1
         IF(z(III).LE.0)GO TO  940
         LUSE = LUSE + LENG(III) + 1
         NAMUSE = NAMUSE + 1
         IF (PRTSW)  WRITE(6,  933)III, NAMES(III), LENG(III), LLOW, LL
 933     FORMAT(' ', I8, 2X,  A8,  I8, I12, I11 )
         GO TO 949
 940     IF ( ALLSW  .AND.  LLOW .GT. 0 )  write (6, 943) III,
     1      NAMES(III), LENG(III), LLOW, LL                               580   
 943     FORMAT(' ', I8, 2X,  A8,  I8, I12, I11, 5X, 'FREED' )
 949  CONTINUE
C
      LSIZE = NLOC - IBASE
      BSIZE = LSIZE/128.
      BAVAIL = NWORDS/128.
      BUSE = LUSE/128.
      BPEAK = IPEAK/128.
      write (6, 5427) NWORDS, BAVAIL, LSIZE, BSIZE, LUSE, BUSE, IPEAK,
     1  BPEAK, IMAX, INAME, NAMUSE, NAMPEK, NUMCOM                        590   
 5427 FORMAT ( '0', T17, 'ALLOCATOR STATISTICS' /
     1    T41, 'WORDS   KILOBYTES' /
     2  ' ALLOCATOR SIZE:',  T38, I8, F10.1  /
     3  ' CURRENT UNCOMPRESSED SIZE:',  T38, I8, F10.1  /
     4  ' CURRENT COMPRESSED SIZE:', T38, I8, F10.1  /
     4  ' PEAK COMPRESSED SIZE:', T38, I8, F10.1 /
     5  '0NAME STACK SIZE:', T38, I8  /
     6  ' CURRENT HIGHEST NAME:',  T38, I8  /
     7  ' CURRENT NAMES IN USE:',  T38, I8   /
     9  ' PEAK NAMES IN USE:', T38, I8 /                                  600   
     A  ' NUMBER OF ALLOCATOR COMPRESSES:', T38, I8 //)
C
      IF (ERRSW)  GO TO 9550
C
      RETURN
C
C     NO ROOM FOR ALLOCATION, STOP OR RETURN A 0
C
 9500 IALLOC = 0
      IF ( .NOT. STOPSW )  RETURN                                         610   
C
      write (6, 9503)
 9503 FORMAT ('-****** EXECUTION IS BEING TERMINATED BY NUMBERED',
     1  ' STORAGE ******' / )
      ERRSW = .TRUE.
      PRTSW = .TRUE.
      ALLSW = .TRUE.
      GO TO 910
C
cccCIB 9550 CALL ERRTRA                                                   620   
 9550 CALL SYSERR
CRA 9550 CALL TRBK( 101 )
CRA      CALL TRBK
CUN 9550 CALL URWALK ( 'IALLOC' )
CVA 9550 CALL LIB%STOP ( %VAL(4) )
      STOP 1234
C
C
      END
C                                                                         630   
C
C*ID* IREDEF
      SUBROUTINE IREDEF ( N1, N2 )
C
C   These are things inside IALLOC which are called instead
C   of being invoked as functions -- violates the standard
C   to do it directly
C
      IDUM = IREDEF_F ( N1, N2 )
      RETURN                                                              640   
C
      ENTRY ISTART ( N1 )
      IDUM = ISTART_F ( N1 )
      RETURN
C
      ENTRY ICLEAR ( N1 )
      IDUM = ICLEAR_F ( N1 )
      RETURN
C
      ENTRY NSSTAT ( N1 )                                                 650   
      IDUM = NSSTAT_F ( N1 )
      RETURN
C
      ENTRY NSCAT ( N1 )
      IDUM = NSCAT_F ( N1 )
      RETURN
      END
C
C
C*ID* NALLOC   C109.PTOLEMY.SOURCE                          PER701  15:   660   
      FUNCTION  NALLOC( IWRDS, NAME )
C
C     NAMED-NUMBERED-STORAGE
C
C     2/3/75 - S. PIEPER
C     11/7/75 - STORE NAME IN COMMON FOR IALLOC ERROR MESSAGES.
C     4/25/78 - CRAY VERSION WITH NO ENTRIES - S.P.
C     12/28/79 - STANDARD MORTRAN - RPG
C     1/16/80 - CHARACTER*8 THENAM - RPG
C     2/25/80 - USE MORTRAN MACRO  USAGE  - RPG                           670   
C     6/13/80 - CHAR*8 NAME FOR CDC - S.P.
C     7/14/80 - CRAY VERSION #2 -- MUCH CLOSER TO IBM - S.P.
C
C     USE AS
C
C     1)   ILOC = NALLOC( IWRDS, 'NAME    ' )
C       WHERE ILOC WILL BE AS IS RETURNED FROM IALLOC.  'NAME'
C       (WHICH NEED NOT BE ON AN 8-BYTE BOUNDRY) IS A NAME THAT WILL
C       BE ASSOCIATED WITH THE ILOC'TH LOCATION IN THE QUERY STACK.
C       IF 'NAME' IS ALREADY IN THE QUERY STACK THEN IT WILL BE FREED     680   
C       AND AN ATTEMPT MADE TO REUSE THE SPACE.  IN THIS CASE IALLOC
C       IS NOT CALLED.
C     2)   ILOC = NAMLOC( 'NAME    ' )
C       RETURNS THE QUERY STACK LOCATION OF AN ALREADY DEFINED NAME.
C       IF 'NAME' IS NOT IN THE QUERY STACK AS A NON-FREED ELEMENT,
C       THEN 0 IS RETURNED.
C
c   Note:  if 'NAME    ' is coded as a literal, it must be 8 chars!
C
      IMPLICIT  REAL*8 (A-H, O-Z)                                         690   
C
      CHARACTER*(*)  NAME
C
      CHARACTER*8  NAMES
      COMMON /NAMCOM/  NAMES(1)
      COMMON /LOCptrs/  L(1)
      COMMON /LENGTH/  LENG(1)
C
      LOGICAL NSPTSW                                                    /usage/
      COMMON/USAGE/IBASE, NWORDS, IMAX, INAME, IGOT, NLOC,
     & IPEAK, NAMPEK, NUMCOM, ICHEC, THENAM, IWASAT, NAMLC2,
     & NSPTSW
 
 
C                                                                         700   
      LOGICAL FNDSW
      CHARACTER*8  BLANK/' '/, THENAM
C
      CHARACTER*1  NAME1(8),  BLNK1 /' '/
      CHARACTER*8  NAME8
      EQUIVALENCE (NAME8, NAME1(1))
C
c
C
C                                                                         710   
C     IS IT DEFINED
C
      ASSIGN 100 TO IGO
      GO TO 800
 100  IF ( .NOT. FNDSW )  GO TO 150
C
C     FOUND - FREE IT AND TRY TO REUSE
C
      L(IL) = IABS(L(IL))
      IF ( L(IL) .EQ. IWRDS )  GOTO 180                                   720   
      L(IL) = -L(IL)
      NAMES(IL) = BLANK
C
 150  THENAM = NAME8
      IL = IALLOC(IWRDS)
      THENAM = BLANK
      IF ( IL .EQ. 0 )  GO TO 200
 180  NAMES(IL) = NAME8
 200  NALLOC = IL
      RETURN                                                              730   
C
C
C     LOCATE A NAME
C
      ENTRY  NAMLOC ( NAME )
C
      ASSIGN  300 TO IGO
      GO TO 800
 300  NALLOC = 0
      IF ( .NOT. FNDSW )  RETURN                                          740   
      IF ( L(IL) .LT. 0 )  RETURN
      NALLOC = IL
      RETURN
C
C     SEARCH FOR A NAME, GET THE NAME FROM THE ARGUMENT LIST
C
 800  FNDSW = .FALSE.
      IL = 0
C
C     PAD IT OUT WITH BLANKS FOR IBM ETC.                                 750   
C
      NAME8 = NAME
C
 810  DO 839  I = 1, IMAX
         IF ( NAMES(I) .NE. NAME8 )  GO TO 839
         IF ( L(I) .GT. 0 )  GO TO 850
C     SEMI FOUND - PERHAPS THERE IS A BETTER ENTRY
      NAMES(I) = BLANK
      IL = I
         FNDSW = .TRUE.                                                   760   
 839  CONTINUE
      GO TO 870
C
C     REALLY FOUND
 850  IL = I
      FNDSW = .TRUE.
 870  GO TO IGO, (100, 300)
C
C
      END                                                                 770   
C*ID* NSDUMP   C109.PTOLEMY.SOURCE                          PER701  15:
      SUBROUTINE  NSDUMP ( IDUMP, NSKEY )
C
C  NSDUMP PRINTS OBJECTS IN NUMBERED STORAGE
C
C     IDUMP - NUMBER OF OBJECT TO BE PRINTED;
C             IF ZERO, ALL OBJECTS ARE PRINTED.
C
C     NSKEY - TYPE OF PRINTOUT DESIRED -
C       1 - AUTOMATIC (DEFAULT REAL ON OTHER MACHINES)                    780   
C       2 - REAL*4
c       3 - integer*4
C          ( FOR 2&3 REAL VERSUS INT IS DETERMINED NUMBER BY NUMBER)
C       4 - INTEGER*2
C       5 - HEXADECIMAL (OCTAL AS APPROPRIATE)
C       6 - REAL*8
C
C     7/24/79 - NSDUMP SEPARATED FROM IALLOC - S.P.
C     2/7/80 - DETERMIN TYPE OF DATA BEING PRINTED - S.P.
C     3/18/80 - TYPE 6 FORCES REAL*8 - S.P.                               790   
C     3/24/80 - BASE 4-BYTE DETECTION ON NUMBER OF OCCURANCES - S.P.
C     6/3/80 - FIXES FOR CDC VERSION - S.P.
C     7/14/80 - CRAY VERSION - S.P.
C     7/16/80 - CDM, CDO PREFIXES - S.P.
C     7/22/80 - CAN'T MIX I-FMT WITH REAL DATA ON CRAY - S.P.
C
C
      implicit real*8 ( a-h, o-z )                                      implicit
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
 
C
      CHARACTER*8  NAMES
      COMMON /NAMCOM/  NAMES(1)
      COMMON /LENGTH/  LENG(1)
C
      COMMON /INTRNL/  UNDEF, NOTDEF
C
      CHARACTER*8  FMTS(22) / '(', 20*' ', '1x )' /
      CHARACTER*8  IFMT / 'I13,' /                                      cnr
CRA      CHARACTER*8  IFMT / 'F13.0,'/                                    810   
      CHARACTER*8  NOTFMT(2) /'5X,4HNot', ' 4HDef.,'/,
     1  GFMT /'G13.5,'/
C
      REAL*4  BUFFER(10)
      INTEGER  IBUF(10)
      EQUIVALENCE ( BUFFER(1), IBUF(1) )
C
CIB      LOGICAL R8SW
cccCIB      LOGICAL*1  EXPMSK /Z7F/,  ZEREXP /Z40/
cccCIB      REAL*8  DOUBLE                                                820   
      INTEGER*4  INTGER(2)
      REAL*4  SINGLE(2)
cccCIB      LOGICAL*1  ONEBYT(8)
      EQUIVALENCE ( INTGER(1), SINGLE(1) )
cccCIB      EQUIVALENCE ( DOUBLE, INTGER(1), ONEBYT(1) )
C
C
      LOGICAL NSPTSW                                                    /usage/
      COMMON/USAGE/IBASE, NWORDS, IMAX, INAME, IGOT, NLOC,
     & IPEAK, NAMPEK, NUMCOM, ICHEC, THENAM, IWASAT, NAMLC2,
     & NSPTSW
 
 
C
C     SETUP SIZE OF EACH WORD TO PRINT                                    830   
C
      GO TO ( 110, 120, 120, 130, 120, 110 ), NSKEY
C
 110  IFAC = 1
      GO TO 200
 120  IFAC = FACFR4
      GO TO 200
 130  IFAC = FACFR2
      GO TO 200
C                                                                         840   
C
 200  NSTRT = 1
      NEND = INAME
      IF ( IDUMP .EQ. 0 )  GO TO 230
      NSTRT = IDUMP
      NEND = IDUMP
C
 230  DO 879  NUMLOC = NSTRT, NEND
         LX = Z(NUMLOC)
         IF ( LX .LE. 0 )  GO TO 879                                      850   
         LX = IFAC*LX - IFAC
         LLX=LENG(NUMLOC)
         LLX = IFAC*LLX
c automatic determination does not work
ccc         IF ( NSKEY .LE. 3  .OR.  NSKEY .EQ. 6 )  GO TO 400
         WRITE(6, 257)NUMLOC, NAMES(NUMLOC), LLX
 257     FORMAT(/'  ----NUMBERED OBJECT',I6,'  NAMED ', A8,
     1     '  WITH', I6, ' ELEMENTS----' /)
C
         GO TO ( 310, 320, 330, 340, 350, 310 ), NSKEY                    860   
C
C
C ccc not true    NOTE...  FOR NSKEY=1, 2 3 OR 6 WE NEVER REALLY GET HER
C ccc not true    AND THUS 310 THROUGH 330 ARE REALLY NOT USED ANYMORE.
C
 310     write (6, 313) ( ALLOC(LX+I), I = 1, LLX )
 313     FORMAT ( 1X, 10G13.5 )
         GO TO 879
 320     write (6, 313) ( ALLOC4(LX+I), I = 1, LLX )
         GO TO 879                                                        870   
 330     write (6, 333) ( ILLOC(LX+I), I = 1, LLX )
 333     FORMAT ( 1X, I10, 9I13 )
         GO TO 879
 340     write (6, 333) ( ILLOC2(LX+I), I = 1, LLX )
         GO TO 879
 350     write (6, 353) ( ILLOC(LX+I), I = 1, LLX )
 353     FORMAT ( 5( Z10, Z9 ) )
         GO TO 879
C
C     DUMP, DUMP8, DUMPR4, DUMPI4.  HERE WE ATTEMPT TO DETERMIN THE       880   
C     NATURE OF THE NUMBER AND PRINT IT ACCORDINGALLY.
C     IF IT HAS BEEN DECLAIRED AS 4- OR 8-BYTE DATA, WE ACCEPT
C     THIS.  ELSE WE TRY TO SEE IF IT SHOULD BE TREATED AS
C     FOUR-BYTE OR EIGHT-BYTE.
C
C     ON CDC & CRAY THERE IS, OF COURSE NO DIFFERENCE
C     THIS HAS NOT YET BEEN THOUGHT ABOUT FOR THE VAX.
C
 400     ASSIGN 510 TO IGOTO
         IF ( NSKEY .EQ. 6 )  ASSIGN 600 TO IGOTO                         890   
cccCIB         I4 = 0
cccCIB         IZERO = 0
cccCIB         R8SW = .FALSE.
cccCIB         IF ( NSKEY .NE. 1 )  GO TO 500
cccCIB         DO 449  I = 1, LLX
cccCIB            DOUBLE = ALLOC(LX+I)
cccCIB            IF ( INTGER(1) .NE. NOTDEF )  GO TO 410
cccCIBC
cccCIBC     FIRST WORD IS NOTDEF; INORDER TO BE DOUBLE PRECISION,
cccCIBC     SECOND MUST BE ZERO.                                          900   
cccCIBC
cccCIB            IF ( INTGER(2) .NE. 0 )  GO TO 470
cccCIB            GO TO 449
cccCIBC
cccCIBC     IF IT IS A DOUBLE WORD OF ZEROS THEN WE CAN'T TELL ANYTHING
cccCIBC     WE USE FLOATING POINT ZEROS CHECKS TO CATCH FLOATING -0
cccCIBC
cccCIB 410        IF ( SINGLE(1) .NE. 0 )  GO TO 420
cccCIB            IF ( SINGLE(2) .NE. 0 )  GO TO 470
cccCIB            IZERO = IZERO + 1                                       910   
cccCIB            GO TO 449
cccCIBC
cccCIBC     IS IT AN INTEGER  (I.E. FLOATING POINT EXPONENT OF 00 OR FF)
cccCIBC
cccCIB 420        IF ( IABS(INTGER(1)) .LT. 16000000 )  GO TO 470
cccCIBC
cccCIBC     IT IS REAL; IS SECOND WORD A REASONABLE INTEGER.  HERE WE
cccCIBC     ASSUME THAT THE UNLIKELY SEQUENCE OF 16 0'S OR 16 1'S
cccCIBC     DOES NOT APPEAR IN THE MIDDLE OF THE DOUBLE PRECISION NUMBER
cccCIBC     IF IT IS ALL ZEROS, WE MUST ASSUME IT IS THE TAIL OF AN       920   
cccCIBC     EXACT FLOATING POINT NUMBER.
cccCIBC
cccCIB            IF ( INTGER(2) .EQ. 0 )  GO TO 449
cccCIB            IF ( IABS(INTGER(2)) .LT. 32000 )  GO TO 430
cccCIB            IF ( INTGER(2) .EQ. NOTDEF )  GO TO 470
cccCIBC
cccCIBC     NOW CHECK FOR A REASONALBLE FLOATING POINT NUMBER.
cccCIBC     IS IT NORMALIZED?  IS EXPONENT MODERATE?
cccCIBC
cccCIB            IF ( ONEBYT(6) .LT. 16 )  GO TO 440                     930   
cccCIB            IF ( IABS( iand(EXPMSK, ONEBYT(5))-ZEREXP ) .GT. 10 )
cccCIB     1         GO TO 449
cccCIB 430        I4 = I4+1
cccCIB            GO TO 449
cccCIB 440        R8SW = .TRUE.
cccCIB 449     CONTINUE
cccCIBC
cccCIBC     WE NEVER FOUND A DEFINET INDICATOR; IS IT LIKELY TO BE
cccCIBC     REAL*4?
cccCIBC     IF 50% OF THE NON-ZERO WORDS ARE REASONABLE 4-BYTE            940   
cccCIBC     INFORMATION, THEN WE USE THEM.
cccCIBC
cccCIB         ASSIGN 600 TO IGOTO
cccCIB         IF ( R8SW .OR.  2*I4 .LT. LLX-IZERO )
cccCIB     1      GO TO 500
cccCIBC
cccCIBC     IT APPEARS TO BE 4-BYTE INFORMATION
cccCIBC
cccCIB 470     LLX = 2*LLX
cccCIB         LX = 2*LX                                                  950   
cccCIB         ASSIGN 510 TO IGOTO
C
C     NOW PRINT ARRAY CHECKING EACH ELEMENT FOR INTEGER VS REAL
C
 500     write (6, 257)  NUMLOC, NAMES(NUMLOC), LLX
         I1 = 2
         I2 = 1
        write (6,503) (alloc(lx+i), i=1,llx)
503     format ( 5g15.5 )
         DO 699  I = 1, LLX                                               960   
            GO TO IGOTO, ( 510, 600 )
 510        INTGER(1) = ILLOC(LX+I)
C
C     DETERMIN HOW TO PRINT THIS NUMBER
C     GET RID OF FLOATING POINT -0
C
            IF ( INTGER(1) .EQ. NOTDEF )  GO TO 550
            IF ( SINGLE(1) .EQ. 0 )  SINGLE(1) = 0
            IF ( IABS(INTGER(1)) .LT. 16000000 )  GO TO 530
C                                                                         970   
C     IT IS FLOATING POINT
C
            FMTS(I1) = GFMT
            GO TO 535
C     IT IS INTEGER
 530        FMTS(I1) = IFMT
CCRA
CCRA     ON CRAY MUST HAVE AGREEMENT BETWEEN FMT AND ARRAY
CCRA
CRA      SINGLE(1) = INTGER(1)                                            980   
 535        I1 = I1 + 1
            IBUF(I2) = INTGER(1)
            I2 = I2+1
            GO TO 650
C
C     IT IS NOT DEFINED
C
 550        FMTS(I1) = NOTFMT(1)
            FMTS(I1+1) = NOTFMT(2)
            I1 = I1+2                                                     990   
            GO TO 650
C
C     REAL*8
C
 600        IF ( ALLOC(LX+I) .EQ. UNDEF )  GO TO 550
            BUFFER(I2) = ALLOC(LX+I)
            FMTS(I1) = GFMT
            I2 = I2+1
            I1 = I1+1
C                                                                        1000   
C     HAVE WE FILLED UP AN OUTPUT LINE
C
 650        IF ( MOD(I,10) .NE. 0  .AND.  I .NE. LLX )  GO TO 699
            FMTS(I1) = FMTS(22)
            I2 = I2-1
            IF ( I2 .EQ. 0 )  GO TO 670
            WRITE ( 6, FMTS )  ( BUFFER(II), II = 1, I2 )
            GO TO 680
 670        WRITE ( 6, FMTS )
 680        I1 = 2                                                       1010   
            I2 = 1
 699     CONTINUE
C
 879  CONTINUE
      RETURN
      END
