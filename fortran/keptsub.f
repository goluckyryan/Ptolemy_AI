      SUBROUTINE KEPT(OBJ,KIND,IFAM,NDIMS,IDIM,LIB,NAME,IER1, MBYTE1)
C
C   Subroutine to read a Speakeasy KEEP object
!
!   This is the  Speakeasy  keep.mor, processed by Mortran for
!   unix machines with the resulting .f file named keptsub.mor for
!   Ptolemy compilations (confusing)!
C
C    Use as
C                                                                          10   
C      CALL KEPT ( OBJ, KIND, IFAM, NDIM, IDIM, LIB, NAME, IER, MBYTE )
C
C     With
C
C     OBJ - the Fortran array into which data is to be placed
C     KIND - set to the Speakeasy kind describing the data:
C       1 - integer*4
C       2 - real*8
C       4 - complex*16
C       5 - real*4                                                         20   
C       6 - name literal  (8 characters in real*8)
C       7 - integer*2
C       8 - logical*1
C       9 - character literal (character in logical*1)
C     IFAM - set to the Speakeasy family describing the structure of
C            the object:
C       0 - scalar
C       1 - arrays (1 to 15 dimensional; only 1 or 2 may be
C                   processed by Speakeasy for now)
C       2 - vectors or matrices (depening on 1 or 2 dimensions)            30   
C       3 - set
C     other values of KIND or IFAM may be recovered but will
C     generally not be usable in Fortran programs.
C     NDIM - set to number of dimensions.  Will be 0 for scalars (IFAM=0
C            1 for sets (IFAM=3)  and  (for all practical purposes)
C            1 or 2 for arrays or vectors/matrices
C     IDIM - an array of dimension 15 that will get the dimensions.
C            This is not refered to for scalars; only the first NDIM
C            elements are defined.
C     LIB - the "library name" from which to read.                         40   
C           This is used as follows:
C           VAX - a character*20 object containing a logical name that
C                 must already be assigned to the directory
C                 from which to read.
C           IBM - a character*8 scalar containing the ddname or
C                 filename that refers to the input dataset
C                 or file.  The allocation (or filedef) must
C                 already be set up.
C           UNIX- a character array (up to 247 characters) containing th
C                 full path of the object, i.e. "/usr/fred" for the obje   50   
C                 /usr/fred/object.
C     NAME - a character*8 containing the member name to be read.
C            This is used as follows:
C           VAX - the file name read is   membername.SPK  in
C                 the directory specified by LIB.
C           IBM -
C             1)  if NAME is blank or 0.D0, then a sequential
C                 dataset (or file) is read.
C             2)  if NAME is not blank or 0.D0, then it is the
C                 member name to be read from the PDS (or CMS              60   
C                 simulated PDS) pointed to by LIB.
C           UNIX- the file name read is  object in the directory
C                 specified by LIB.
C     IER - return code.  0 means everything done o.k.; nonzero
C           means an error occured as follows:
C     -1 - not a Speakeasy Keep type object
C     -2 - Speakeasy KEEPLIST or CHECKPOINT
C     -3 - Invalid KLASS in a Speakeasy-III KEEP
C     -4 - MBYTE/8 is smaller than the double-word size of the data;
C          only the first MBYTE/8 double words were recovered.             70   
C      VAX -
C      +1 - error during attempted open of the file.
C       5 - I/O error
C       9 - unexpected end of file or file is completely empty.
C      UNIX -
C      +1 - error during attempted open of the file.
C       5 - I/O error
C       9 - unexpected end of file or file is completely empty.
C      IBM -
C      +2 - the ddname or filename refered to by LIB is not defined        80   
C       3 - Member name not found in PDS.
C       4 - bad or undefined RECFM, LRECL or BLKSIZE.
C       5 - I/O error.
C       8 - Not enough memory for I/O buffers
C       9 - unexpected end of file or file is completely empty.
C     In addition in OS, IER may have the hexadecimal format
C     of  00ccc0rr  where ccc is a system completion code and rr is
C     the associated return code.  For this reason it may be
C     best to always print IER with Z8 format.
C                                                                          90   
C     MBYTE - size of the array OBJ in bytes.  This should be a
C           multiple of 8 since the data is always moved in groups
C           of 8 bytes;  thus a nine-element character array
C           requires MBYTE=16 to be fully recovered.
C           If MBYTE specifies too small an array, then the array
C           will be filled and IER=-4 will be returned.
C
C        A special case occurs if MBYTE is 0:
C           then only the KIND,...,IDIM information is returned.
C           Following the KEPT call, space may be allocated for           100   
C           OBJ and KEPTCN called to continue the reading, or
C           KEPTCL may be used to cancel the reading;  note that
C           after calls to KEPT with MBYTE=0, one of KEPTCN or
C           KEPTCL must be called.  These calls are described next.
C
C
C     CALL KEPTCN ( OBJ, IER, MBYTE )
C
C     This call may be used after a call to KEPT that was made
C     with MBYTE=0.  The processing of the KEPT call will be              110   
C     concluded by loading the data into OBJ.  The three
C     arguments have the same meanings as above except that
C     MBYTE=0 makes no sense.
C
C
C     CALL KEPTCL
C
C     This call may be used after a call to KEPT that was made
C     with MBYTE=0 to cancel further processing of the object.
C     There are no arguments.                                             120   
C
C
C   5/1/83 - WRITTEN FOR VAX SPEAKEASY IV - JEG
C   9/5/84 - UPGRADED TO BE COMPATIBLE WITH IBM VERSIONS.   W.T.
C   1/15/85 - FINISH UP MERGE OF IBM AND VAX VERSIONS - S.P.
C   3/25/91 - SUPPORT FOR UNIX VERSIONS ADDED - T.H.
C   4/18/91 - ALLOW AMODE=31 on IBM - s.p.
C   4/29/91 - MAKE OBJ CHARACTER*8 - T.H.
C   11/10/92 - Change CVA -> VMS for new macros.   W.T.
C   3/20/04 - make specific unix fortran version for Ptolemy - s.p.       130   
C
      CHARACTER*8 OBJ(1)
      INTEGER KIND,IDIM(1)
C
C   FOR THE LIBRARY NAME
C
      CHARACTER*1 LIB(1)
      CHARACTER*1 FILE(257)
      CHARACTER*257 FFILE
      EQUIVALENCE(FILE(1),FFILE)                                          140   
C
C   FOR USING THE MEMBER NAME
C
      CHARACTER*8 NAME
      CHARACTER*8 NAME2
      CHARACTER*1 NAME1(8)
      REAL*8 RNAME
      EQUIVALENCE (RNAME, NAME2, NAME1(1))
C
C   FOR THE VARIOUS DATA TYPES ON THE CARD                                150   
C
      REAL*8 RCARD(10)
      CHARACTER*8 CARD(10)
      INTEGER*4 CARD4(20)
      INTEGER*2 CARD2(40)
      CHARACTER*1 CARD1(80)
      EQUIVALENCE (RCARD(1),CARD(1))
      EQUIVALENCE (CARD1(1),CARD(1))
      EQUIVALENCE (CARD2(1),CARD(1))
      EQUIVALENCE (CARD4(1),CARD(1))                                      160   
      INTEGER*4  ZFFA3
C
      SAVE IOFF, NWDS                                                   cni
C
C
      ICALL = 1
      IER = 0
C
      NAME2 = NAME
      DO 310 K = 1,247                                                    170   
         L = K
         IF (LIB(K).EQ.' ') GO TO 48
         FILE(K) = LIB(K)
310   CONTINUE
48    IF (FILE(L-1) .NE. '/' .and. l /= 1 ) THEN
         FILE(L) = '/'
      ELSE
         L = L-1
      END IF
      DO J = 1,8                                                          180   
         K = J
         IF (NAME1(J).EQ.' ') GO TO 66
         FILE(L+J) = NAME1(J)
      enddo
      k = 9
66    FILE(L+k) = '.'
      FILE(L+k+1) = 'k'
      FILE(L+k+2) = 'e'
      FILE(L+k+3) = 'e'
      FILE(L+k+4) = 'p'                                                   190   
      FILE(L+k+5) = CHAR(0)
      CALL SROPEN(IDCB,ICODE,FFILE)
      IF (ICODE .NE. 0) GO TO 915
      CALL SRREAD(IDCB,ICODE,CARD(1),NUMCHR,80,1)
      IF (ICODE .EQ. 3) GO TO 910
      IF (ICODE .NE. 0) GO TO 920
C
C   OPEN THE FILE AND READ THE FIRST CARD
C
C                                                                         200   
C     CHECK FOR SPEAKEASY-III KEEP
C
C
C     PROCESS 1ST CARD OF THE KEEP FILE
C
      IF (CARD4(1).NE.ZFFA3.OR.CARD4(15).NE.ZFFA3) GO TO 900
      IF (CARD4(2).NE.1) GO TO 900
      IF (CARD4(3).NE.0) GO TO 900
      IF (CARD4(6).NE.0) GO TO 900
C                                                                         210   
      IF (CARD4(16).NE.1) GO TO 905
      NWDS = CARD4(18)
C
      CALL SRREAD(IDCB,ICODE,CARD(1),NUMCHR,80,1)
      IF (ICODE .EQ. 3) GO TO 910
      IF (ICODE .NE. 0) GO TO 920
      LHEAD = ICHAR(CARD1(3))
      NWDS = NWDS-LHEAD-1
C
      KIND=CARD2(3)                                                       220   
      IFAM = ICHAR(CARD1(7))
      NDIMS = ICHAR(CARD1(8))
C
      ISIZE = ICHAR(CARD1(9))
      IBOUND = ICHAR(CARD1(10))
C
      IF ( NDIMS .EQ. 0 ) GO TO 140
      DO 129 K = 1,NDIMS
 129  IDIM(K) = CARD4(K+3)
C                                                                         230   
C     SKIP TO END OF THE HEADER IF NECESSARY
C
 140  ISKIP = LHEAD/10
      IF ( ISKIP .EQ. 0 ) GO TO 170
      DO 159 MM = 1,ISKIP
      CALL SRREAD(IDCB,ICODE,CARD(1),NUMCHR,80,1)
      IF (ICODE .EQ. 3) GO TO 910
      IF (ICODE .NE. 0) GO TO 920
      LHEAD = LHEAD-10
 159  CONTINUE                                                            240   
 170  IOFF = LHEAD
      GO TO 300
C
C     PROCESS SPEAKEASY-III OBJECT
C
C
C     READY TO PROCESS THE DATA
C
C     NWDS IS THE NUMBER OF DOUBLE WORDS OF DATA IN THE FILE
C     IOFF IS THE NUMBER OF DOUBLE WORDS IN THE PRESENT CARD              250   
C          IMAGE PRECEEDING THE DATA
C
 300  IER1 = 0
      IF ( MBYTE1 .EQ. 0 ) RETURN
      MBYTE = MBYTE1
      GO TO 350
C
C     THIS ENTRY CONTINUES A READ
C
      ENTRY KEPTCN ( OBJ, IER2, MBYTE2 )                                  260   
      ICALL = 2
      MBYTE = MBYTE2
C
 350  IER = 0
      IF ( NWDS .EQ. 0 ) GO TO 990
      MAXWRD = MBYTE/8
      IF ( MAXWRD .LT. NWDS ) IER = -4
      IF ( MAXWRD .EQ. 0 ) GO TO 990
      NWDS = MIN( MAXWRD, NWDS )
      L = 0                                                               270   
C
C     PROCESS THE DATA ON THIS CARD
C
 330  NDO = MIN( NWDS, 10-IOFF )
      DO 339 I = 1, NDO
      OBJ(L+I) = CARD(IOFF+I)
 339  CONTINUE
      L = L+NDO
      NWDS = NWDS-NDO
      IOFF = 0                                                            280   
      IF ( NWDS .EQ. 0 ) GO TO 990
      CALL SRREAD(IDCB,ICODE,CARD(1),NUMCHR,80,1)
      IF (ICODE .EQ. 3) GO TO 910
      IF (ICODE .NE. 0) GO TO 920
      GO TO 330
C
C     THIS ENTRY CLOSES WITHOUT FURTHER READING
C
      ENTRY KEPTCL
      ICALL = 3                                                           290   
      GO TO 990
C
C**************************************ERROR MESSAGES
C
 900  IER = -1
      GO TO 990
 905  IER = -2
      GO TO 990
 910  IER = +9
      GO TO 990                                                           300   
 915  IER = 1
      GO TO 990
 920  IER = 5
      GO TO 990
 930  IER = -3
      GO TO 990
 980  IF ( IER .EQ. 3 ) IER = 9
      IF ( IER .EQ. 4 ) IER = 5
C
 990  CALL SRCLSE(IDCB)                                                   310   
      IF ( ICALL .EQ. 1 ) IER1 = IER
      IF ( ICALL .EQ. 2 ) IER2 = IER
      RETURN
      DATA ZFFA3 /Z'FFA30000'/
      END
