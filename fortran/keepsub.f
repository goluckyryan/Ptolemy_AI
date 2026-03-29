      SUBROUTINE KEEP ( OBJ,KIND,IFAM,NDIM,IDIM,LIB,NAME,IER )
C
C   Subroutine to create a Speakeasy KEEP object
!
!   This is the  Speakeasy  keep.mor, processed by Mortran for
!   unix machines with the resulting .f file named keepsub.mor for
!   Ptolemy compilations (confusing)!
C
C    Use as
C                                                                          10   
C      CALL KEEP ( OBJ, KIND, IFAM, NDIM, IDIM, LIB, NAME, IER )
C
C     With
C
C     OBJ - the Fortran array to be output
C     KIND - the Speakeasy kind describing the data:
C       1 - integer*4
C       2 - real*8
C       4 - complex*16
C       5 - real*4                                                         20   
C       6 - name literal  (8 characters in character*8)
C       7 - integer*2
C       8 - logical*1 or integer*1
C       9 - character literal (character*1)
C     Note:  If KIND = 1, 5, or 7 objects are recovered in
C     Speakeasy, they should be immeadiately converted to
C     standard Real*8 objects with a statement of the form
C         X = REAL8( X )
C     because almost all Speakeasy statements cannot process
C     such kinds.                                                          30   
C     IFAM - the Speakeasy family describing the structure of
C            the object:
C       0 - scalar
C       1 - arrays (1 to 15 dimensional; only 1 or 2 may be
C                   processed by Speakeasy for now)
C       2 - vectors or matrices (depening on 1 or 2 dimensions)
C       3 - set
C     other values of KIND or IFAM may not be used.
C     NDIM - number of dimensions.  This must be 0 for scalars (IFAM=0),
C            1 for sets (IFAM=3)  and  (for all practical purposes)        40   
C            1 or 2 for arrays or vectors/matrices
C     IDIM - an array of dimension NDIM containing the dimensions.
C            This is not refered to for scalars.
C     LIB - the "library name" on which to write.
C           This is used as follows:
C           VAX - a character*20 object containing a logical name that
C                 must already be assigned to the directory
C                 in which to write.
C           IBM - a character*8 scalar containing the ddname or
C                 filename that refers to the output dataset               50   
C                 or file.  The allocation (or filedef) must
C                 already be set up.
C           UNIX- a character array (up to 247 characters) containing th
C                 full path of the object, i.e. "/usr/fred" for the obje
C                 /usr/fred/object.
C           Win32- a character array (up to 247 characters) containing t
C                 full path of the object, i.e. "\usr\fred" for the obje
C                 \usr\fred\object.
C     NAME - a character*8 containing the member name to be written.
C            This is used as follows:                                      60   
C           VAX - the file name written is   membername.SPK  in
C                 the directory specified by LIB.
C           IBM -
C             1)  if NAME is blank or 0.d0, then a sequential
C                 dataset (or file) is written.
C             2)  if NAME is not blank or 0.d0, then it is the
C                 member name to be put in the PDS (or CMS
C                 simulated PDS) pointed to by LIB.
C           UNIX- the file name written is  object in the directory
C                 specified by LIB.                                        70   
C     IER - return code.  0 means everything done o.k.; nonzero
C           means an error occured as follows:
C       -1 - KIND, NDIM or IDIM are invalid
C      VAX -
C       +1 - error during attempted open of the file.
C      IBM -
C       2 - the ddname or filename refered to by LIB is not defined
C       3 - (CMS only)  cannot read first record of an existing
C           file to check the type.
C       4 - bad or undefined RECFM, LRECL or BLKSIZE.                      80   
C       5 - I/O error (also means full disk in CMS)
C       6 - (CMS) attempt to write on a read-only disk,
C           (OS) No room in directory
C       7 - I/O error updating the directory
C       8 - Not enough memory for I/O buffers
C       9 - OBTAIN failed for KEEPCK call.
C      10 - LIB does not point to a PDS (KEEPCK call only).
C      11 - Could not open pds to read directory (KEEPCK)
C      12 - Error reading directory (KEEPCK)
C      13 - PDS last record written pointer has been corrupted --          90   
C           do not attempt to write into this pds anymore!
C      UNIX -
C       +1 - error during attempted open of the file.
C      Win32 -
C       +1 - error during attempted open of the file.
C
C     In addition in OS, IER may have the hexadecimal format
C     of  00ccc0rr  where ccc is a system completion code and rr is
C     the associated return code.  The most common cases are
C     ccc = B37, D37 or E37 (all with rr = 04) which mean the disk        100   
C     is full.  Because of this possible hexadecimal format,
C     it is best to always print IER with Z8 format.
C
C
C      CALL KEEPCK ( OBJ, KIND, IFAM, NDIM, IDIM, LIB, NAME, IER )
C
C     This entry may be used to keep an object with checking of
C     the validity of the PDS pointers.  There have been a few
C     reported problems in which the PDS last record pointer has
C     been messed up in a way that causes new information to              110   
C     overwrite existing data in the PDS.  The extra code invoked
C     by this call attempts to detect this situation before the PDS
C     is damaged.  If the error is detected, diagnostic messages
C     are printed and nothing is written.  There is no guarantee
C     that this call will detect the damage; it is possible that
C     the PDS will be destroyed even by this call.
C
C     All of the arguments have the same meanings as for   KEEP
C
C     The special processing of this entry is done only in the            120   
C     OS version.  For all other versions, this entry is equivalent
C     to the KEEP entry.
C
C
C   5/1/83 - WRITTEN FOR VAX SPEAKEASY IV  - JEG
C   8/24/84 - UPGRADED TO BE COMPATIBLE WITH IBM VERSIONS.   W.T.
C   1/10/85 - FINISH FIRST JOINT IBM&VAX VERSIONS - S.P.
C   1/30/86 - ADD TEST FOR CORRUPTED PDS POINTERS - S.P.
C   9/10/90 - convert to Fortran-77 - cleanup |IBM| errors - s.p.
C   3/9/89 - ADDED SUPPORT FOR UNIX VERSIONS - T.H.                       130   
C   4/17/91 - allow AMODE=31 on IBM - s.p.
C   4/29/91 - MAKE OBJ CHARACTER*8 FOR ALL MACHINES - T.H.
C   11/10/92 - Change CVA -> VMS for new macros.   W.T.
C   8/23/99 - Add code for Windows.  W.T.
C   3/20/04 - make specific unix fortran version for Ptolemy - s.p.
C
      CHARACTER*8 OBJ(1)
      INTEGER IDIM(1)
C
      CHARACTER*1 LIB(1)                                                  140   
      CHARACTER*1 FILE(257)
      CHARACTER*257 FFILE
      EQUIVALENCE(FILE(1),FFILE)
C
      CHARACTER*8 NAME, NAME2
      CHARACTER*1 NAME1(8)
      EQUIVALENCE(NAME2,NAME1(1))
C
C   FOR THE VARIOUS DATA TYPES ON THE CARD
C                                                                         150   
      CHARACTER*8 CARD(10)
      REAL*8 RCARD(10)
      INTEGER*4 CARD4(20)
      INTEGER*2 CARD2(40)
      CHARACTER*1 CARD1(80)
      EQUIVALENCE (CARD1(1),CARD(1))
      EQUIVALENCE (CARD2(1),CARD(1))
      EQUIVALENCE (CARD4(1),CARD(1))
      EQUIVALENCE (RCARD(1),CARD(1))
C                                                                         160   
C   THE SIZES (IN BYTES) AND BOUNDARIES AS A FUNCTION OF SPEAKEASY KIND
C
      INTEGER  ISIZE(9)
      INTEGER  IBOUND(9)
C
      CHARACTER*8  TIME
      CHARACTER*8  DATE
      CHARACTER*8  ENDDAT
      CHARACTER*8  ZERNAM
      integer*8 izernam                                                   170   
      equivalence ( zernam, izernam )
      INTEGER*4  ZFFA3
      INTEGER*4  VER
C
      DATA ISIZE /4,8,0,16,4,8,2,1,1/
      DATA IBOUND /4,8,0,8,4,8,2,1,1/
      DATA TIME /' '/
      DATA DATE /' '/
      DATA ENDDAT /'ENDOFDAT'/
      DATA iZERNAM /Z'0000000000000000'/                                  180   
      DATA ZFFA3 /Z'FFA30000'/
      DATA VER /4061000/
C
C
      GO TO 20
C
      ENTRY KEEPCK ( OBJ,KIND,IFAM,NDIM,IDIM,LIB,NAME,IER )
C
C
C                                                                         190   
C     READ THE DIRECTORY AND GT MAX TTR
C
C
C     THE MAX. DIRECTORY TTR MUST BE LESS THAN THE CURRENT TTR
C
C
C
C
C
 20   IER = 0                                                             200   
      NAME2 = NAME
C
C     VALIDATE SOME OF THE INPUT AND GET THE NUMBER OF ELEMENTS
C     TO BE PROCESSED; COMPLEX COUNTS AS 2 REALS.
C
      IF ( KIND .LE. 0 .OR. KIND .GE. 10 ) GO TO 900
      IF ( KIND .EQ. 3 ) GO TO 900
      IF ( NDIM .LT. 0 .OR. NDIM .GT. 15 ) GO TO 900
      NELS = 1
      IF ( NDIM .EQ. 0 ) GO TO 30                                         210   
      DO 29 I = 1, NDIM
         IF ( IDIM(I) .LT. 0 ) GO TO 900
         NELS = NELS*IDIM(I)
  29  CONTINUE
  30  IF ( KIND .EQ. 4 ) NELS = 2*NELS
C
C   OPEN THE FILE FOR WRITING
C
C
C   Construct the full path name of the file.                             220   
      DO 310 K = 1,247
         L = K
         IF (LIB(K).EQ.' ') GO TO 48
         FILE(K) = LIB(K)
310   CONTINUE
48    IF (FILE(L-1) .NE. '/' .and. l /= 1 ) THEN
         FILE(L) = '/'
      ELSE
         L = L-1
      END IF                                                              230   
      DO J = 1,8
         K = J
         IF (NAME1(J).EQ.' ') GO TO 66
         FILE(L+J) = NAME1(J)
      enddo
      k = 9
66    FILE(L+k) = '.'
      FILE(L+k+1) = 'k'
      FILE(L+k+2) = 'e'
      FILE(L+k+3) = 'e'                                                   240   
      FILE(L+k+4) = 'p'
      FILE(L+k+5) = CHAR(0)
      CALL SWOPEN(IDCB,ICODE,FFILE,1,0)
      IF (ICODE .NE. 0) GO TO 920
C
C
C   SET UP THE FIRST CARD
C
 200  CARD4(1) = ZFFA3
      CARD4(2) = 1                                                        250   
      CARD4(3) = 0
      CARD4(4) = 1
      CARD4(5) = VER
      CARD4(6) = 0
      CARD4(6) = 0
      CARD4(7) = ZFFA3
      CARD4(8) = 4
      CARD4(9) = 0
      CARD4(10) = 2
      CARD(6) = TIME                                                      260   
      CARD(7) = DATE
      CARD4(15) = ZFFA3
      CARD4(16) = 1
      CARD4(17) = 0
C
      NN = 0
      IFACT = 8/IBOUND(KIND)
      IF (((NELS/IFACT)*IFACT) .NE. NELS) NN = 1
      LHEAD = 2+NDIM/2
      NTWORD = (NELS * (IBOUND(KIND)))/8 + LHEAD + NN+1                   270   
      CARD4(18) = NTWORD
      CARD(10) = NAME
C
C   WRITE OUT THE FIRST CARD
C
      CALL SWRITE(IDCB,ICODE,CARD(1),80)
      IF (ICODE .NE. 0) GO TO 920
C
C   SET UP THE SECOND CARD INFORMATION
C                                                                         280   
      CARD1(3) = CHAR(LHEAD)
      CARD1(4) = CHAR(0)
      CARD2(3) = KIND
      CARD1(7) = CHAR(IFAM)
      CARD1(8) = CHAR(NDIM)
      CARD1(9) = CHAR(ISIZE(KIND))
      CARD1(10) = CHAR(IBOUND(KIND))
C
      IF ( NDIM .EQ. 0 ) GO TO 330
      DO 329 K = 1,NDIM                                                   290   
 329  CARD4(K+3) = IDIM(K)
C
 330  IOUT = 0
      NWORD = NTWORD - LHEAD-1
      NLEFT = NWORD
C
C     LOOP TO OUTPUT THE DATA, STARTING SOMEWHERE INSIDE THE
C     SECOND CARD.
C
 350  IF ( NLEFT .EQ. 0 ) GO TO 400                                       300   
      LMAX = MIN0(NLEFT,10-LHEAD)
      DO 369 K=1,LMAX
         CARD(K+LHEAD) = OBJ(K+IOUT)
 369  CONTINUE
      IOUT = IOUT + LMAX
      LHEAD = LHEAD+LMAX
      NLEFT = NLEFT-LMAX
C
      IF ( LHEAD .LT. 10 ) GO TO 400
C                                                                         310   
C   WRITE OUT THE CARD
C
      CALL SWRITE(IDCB,ICODE,CARD(1),80)
      IF (ICODE .NE. 0) GO TO 920
C
C     SUBSEQUENT CARDS ARE FULLY DATA
C
      LHEAD = 0
      GO TO 350
C                                                                         320   
C   WRITE OUT THE LAST CARD, CLOSE AND RETURN
C
 400  CARD(LHEAD+1) = ENDDAT
      LHEAD = LHEAD+1
      IF ( LHEAD .GT. 10 ) GO TO 420
      DO 409 I = LHEAD, 10
         RCARD(I) = 0
 409  CONTINUE
 420  CONTINUE
      CALL SWRITE(IDCB,ICODE,CARD(1),80)                                  330   
      IF (ICODE .NE. 0) GO TO 920
      CALL SWCLSE(IDCB)
      RETURN
 920  IER = 1
      RETURN
C
C
C    CLEAN UP THE IBM-ERROR RETURNS TO MAKE EASIER TO PROCESS
C
C                                                                         340   
 900  IER = -1
      RETURN
C
      END
