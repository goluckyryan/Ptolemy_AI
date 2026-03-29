C*ID* PTOLEMY  C109.PTOLEMY.SOURCE                          PER701  16:
CRA      PROGRAM PTOLEMY
C
C     MAIN PROGRAM FOR THE DWBA CODE
C
C     THIS ACTS AS A SWITCHER FOR CONTRL.  CONTRL IS CALLED AND
C     READS INPUT.  EACH TIME A SUBSTANTIAL CALCULATION IS TO BE DONE
C     CONTRL SETS  IGOTO  AND  JGOTO  AND RETURNS TO HERE.  IGOTO
C     IS THEN USED TO SWITCH TO THE APPROPRIATE ROUTINE.  UPON
C     COMPLETION OF THE CALCULATION, CONTRL IS AGAIN CALLED.  THIS MODE    10   
C     OF OPERATION IS USED TO ALLOW  CONTRL  TO BE OVERLAYED AGAINST THE
C     LARGE CALCULATIONAL BLOCKS.
C
C     1/18/75 - MADE FROM CONTRL - S. PIEPER
C     1/18/75 - REARRANGED TO BE OVERLAYED - S. P.
C     1/24/75 - FOR SPLIT UP INELDC
C     2/5/75 - FOR SPLIT UP WAVEF
C     2/7/75 - FIRST ADDITION OF GRDSET
C     2/10/75 - SEPARATE GRIDSET, INELDC, DSIGMA; TEMP1,2,3
C     2/25/75 - RETURN CODES FROM WAVEF, WAVSET                            20   
C     4/17/75 - NEW COMMONS
C     4/23/75 - CALL WAVPOT AFTER GRDSET
C     5/5/75 - COMPLETE DWBA PATH WITHOUT SEPARATE CALLS
C     5/22/75 - STAND ALONE SWITH OFF ON WAVSET CALLS
C     5/28/75 - PRBPRT AND ANGSET SEPARATED FROM GRDSET.
C     6/9/75 - CALL LINTRP WHEN NEEDED
C     7/8/75 - CALL GETSCT PRIOR TO GRDSET.
C     8/7/75 -  CALL FITELE ADDED
C     10/13/75 - EXTERNALS FOR OVERLAY PROCESSING
C     11/6/75 - FITEL INTO TWO PIECES                                      30   
C     9/1/76 - ADD VCSQ12 TO EXTERNALS
C     12/26/76 - INELASTIC SCATTERING ADDED
C     4/18/77 - CHANGE PARM TO PRM TO AVOID COMMON CONFLICTS
C     4/26/77 - FITSW IN WAVSET
C     10/14/77 - PLMSUB ADDED TO EXTERNALS FOR OVERLAY; DATAN2 OUT
C     11/20/77 - NEW WAVPOT CALLING LOCATION.
C     12/1/77 - INTRPC ADDED TO EXTERNALS - S.P.
C     12/20/77 - FIX USEHS ERROR MADE ON 11/20 - S.P.
C     5/8/78 - CALL COUPL2; SPLIT OFF LINTRP CALL - S.P.
C     6/27/78 - CALL TO TABLES FOR CDC                                     40   
C     5/2/79 - ADD COUPL2 STUFF AND WAVETC TO EXTERNAL LIST - S.P.
C     6/10/79 - CMS VERSIN - S.P.
C     7/5/79 - ADD SDERIV TO EXTERNALS - S.P.
C     12/28/79 - ADD CLEBSH, 3J, 6J, 9J TO EXTERNALS - S.P.
C     12/27/79 - CFR$F, MORE COMMENTS - RPG
C     1/15/80 - COUPLED CHANNELS CALLS - S.P.
C     2/8/80 - CALL COUPLN - S.P.
C     3/21/80 - CALL FFACST - S.P.
C     6/3/80 - CALL TABLES ONLY IN THE BROOKHAVEN VERSION - S.P.
C     6/25/80 - COUPL2 ONLY IN IBM, DOESN'T WORK IN CDC - S.P.             50   
C     7/10/80 - CHANGE CDC BUFFER SIZES, CIB ON EXTERNALS,
C               REMOVE LOTS IN CRAY VERSION WITH CNR - S.P.
C     7/15/80 - COMMONS FOR THE CRAY VERSION - S.P.
C     7/23/80 - VCSBLK COMMON - S.P.
C     7/28/80 - GRABCOM WITH DATA STATEMENT FOR CRAY - S.P.
C     7/30/80 - READ TO END OF FILE ON CRAY TO CEAR UP I.O. - S.P.
C     8/13/80 - REMOVE CDC TAPE10 DEFINITION - S.P.
C     12/17/80 - REMOVE CRA FROM STMNT 94; CUPROT IN EXTERNALS - S.P.
C     4/3/81 - REMOVE COUPL2 - S.P.
C     9/21/81 - CUPROT EXTERNAL REMOVED, ADD CUPSPN, CUPAB - S.P.          60   
C     12/14/81 - GETSCT NOW CALLED ALWAYS - S.P.
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
 
      LOGICAL*1 WASSET                                                  /intrnl/
      COMMON /INTRNL/  UNDEF, NOTDEF, ICHANB, ICHANW, ICNTT, EBNDS(2),
     &  IDONE, ISTRIP, IHSAVE, NEXTSV, IPOT, IWFN, LIMOST, IEXCIT,
     &  R0MASS, WASSET(40), RATMAS,
     &  LSPECS(4), NODESP(4), LSPCPT(2), NODEPT(2),ICD3I
      real*4 undef4
      equivalence ( undef4, notdef )
 
C
CCRA      THE CRAY LOADER IS PRIMATIVE AND WE NAME ALL COMMONS IN
CCRA      MAIN PROGRAM JUST TO BE SAFE.
CCRA
CRA      ALLOCS; USAGE; FLOAT; INTGER; JBLOCK; LNKBLK; HEDCOM; KANDM       70   
CRA      CNSTNT; FORMF; INELCM; CCBLK; FTIMEI; WAVCOM; TEMPVS
CRA      GRIDCMFORINEL; INPBUF; LOCFITFORDWBA
CRA      COMMON /LENGTH/ LENG(100);   COMMON /NAMCOM/ NAMES(100)
CRA      COMMON /FACTRL/ IFAC1, IFAC2, FTB(97)
CRA      COMMON /LOGFAC/ ILOG1, ILOG2, FLTB(101)
CRA      COMMON /VCSBLK/ VCSQ12(3,11)
CCRA
CCRA      FOLLOWING COMMON IS FOR RETURNING RETURN CODES FROM CRAY'S
CCRA      CRUMMY OVERLAY PROCESSOR.
CCRA                                                                       80   
CRA      COMMON /IRETURN/  IPRM1, IGOTO, JGOTO
CCRA
CCRA      FOLLOWING ALLOWS GRAB TO KNOW WHEN HLM IN JC TABLE HAS CHANGED
CCRA
CRA      COMMON /GRABCOM/  LASTHLM
C
      LOGICAL ONESW
CRA      DATA LASTHLM /0/
C
c                                                                          90   
C     INITIALIZE WITH ALL THE DEFAULTS
C
      JGOTO = 3
C
C     START OR CONTINUE THE INPUT PROCESS
C
 90   CALL CONTRL ( IGOTO, JGOTO, IPRM1, IPRM2, IPRM3, PRM1 )
C
      IF ( IGOTO .NE. 0 )  GO TO 95
CCRA                                                                      100   
CCRA     CRAY SEEMS TO NOT LIKE TO HAVE THE INPUT FILE LEFT WITH
CCRA     UNREAD RECORDS OR END OF FILE.  WE FINISH IT OFF NOW.
CCRA
CRA  92  READ ( 5, 93, END=94 ) AWORD
CRA  93  FORMAT ( A8 )
CRA      GO TO 92
  94  STOP
C
C
C     WHAT IS TO BE DONE NEXT                                             110   
C
C     IGOTO = 4 AVAILABLE FOR RECYCLING
C
  95  ONESW = .TRUE.
      GO TO ( 100, 200, 300, 90, 500, 600, 620, 650, 670, 700, 720, 740,
     1  800, 850, 870 ), IGOTO
C                                                                       cnr
C     IGOTO = 1:  ; AFTER PROJECTILE OR TARGET.                         cnr
C                                                                       cnr
 100  CALL BOUND ( IPRM1, IPRM2, IPRM3 )                                  120   
      GO TO 90                                                          cnr
C                                                                       cnr
C     IGOTO = 2:  ; AFTER FIT OR FITELASTIC                             cnr
C                                                                       cnr
200   CALL FITEL ( IPRM1 )                                              cnr
      IF ( IPRM1 .EQ. 0 )  GO TO 90                                     cnr
      CALL FITEL2 ( IPRM1 )                                             cnr
      GO TO 90                                                          cnr
C                                                                       cnr
C     IGOTO = 3:  ; AFTER CALLTSO.                                        130   
C                                                                       cnr
COS 300  CALL TSO ( PRM1 )
 300  write (6, 303)                                                    cni
 303  FORMAT ( '0**** TSO KEYWORD NOT AVAILABLE IN this VERSION')       cni
      GO TO 90                                                          cnr
C
C     IGOTO = 5:  ; AFTER INCOMING OR OUTGOING.
C
 500  CALL WAVSET ( IPRM1, .FALSE., .FALSE. )
      GO TO 90                                                            140   
C
C     GRIDSETUP
C
C     IGOTO = 6:  ; AFTER GRIDSETUP.
C
 600  CALL PRBPRT ( IPRM1 )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      IF ( PROBLM .NE. 24 )  GO TO 605
C
C     COUPLED CHANNELS CALLS                                              150   
C
      CALL BASLBL ( IPRM1  )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      CALL BASCPL ( IPRM1 )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
C
 605  CALL GETSCT ( IPRM1 )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      IF ( ISTRIP .NE. 0 )  CALL GRDSET ( IPRM1 )                       cnr
C                                                                         160   
C     FOLLOWING DONE FOR ALL CALCULATIONS
C
      IF ( ISTRIP .EQ. 0 )  CALL INGRST ( IPRM1 )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      CALL ANGSET ( IPRM1 )
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      IF ( ISTRIP .EQ. 0 )  GO TO 615
      IF ( ONESW )  GO TO 90
      GO TO 620
C                                                                         170   
C     SETUP COULOMB INTEGRALS FOR INELASTIC AND CC
C
 615  CALL COULST ( IPRM1 )
C
      IF ( IPRM1 .EQ. 0  .OR.  ONESW )  GO TO 90
C
C     IGOTO = 7:  ; AFTER RADIALINTEGRALS
C
 620  CONTINUE
      IF ( PROBLM .EQ. 20 )  CALL INELDC ( IPRM1 )                        180   
      IF ( PROBLM .EQ. 22 )  CALL INRDIN ( IPRM1 )                      cnr
      IF ( PROBLM .EQ. 24 )  CALL COUPLN ( IPRM1 )
      IF ( IPRM1 .EQ. 0  .OR.  ONESW )  GO TO 90
C
C     IGOTO = 8:  ; AFTER LINTERPOLATION
C
 650  IPRM1 = 1
      IF ( PROBLM .EQ. 24 )  CALL FFACST (IPRM1)
      IF ( IPRM1 .EQ. 0 )  GO TO 90
      CALL LINTRP                                                         190   
      IF (ONESW)  GO TO 90
C
C     IGOTO = 9:  ; AFTER CROSSSECTION.
C
 670  CALL  XSECTN ( IPRM1 )
      GO TO 90
C                                                                       cnr
C     IGOTO = 10, 11, 12:  ; AFTER DOTEMP1, DOTEMP2, DOTEMP3.           cnr
C                                                                       cnr
 700  CALL TEMP1                                                          200   
      GO TO 90                                                          cnr
 720  CALL TEMP2                                                        cnr
      GO TO 90                                                          cnr
 740  CALL TEMP3                                                        cnr
      GO TO 90                                                          cnr
C                                                                       cnr
C     IGOTO = 13:  ; AFTER SCATTERING.                                  cnr
C                                                                       cnr
 800  CALL WAVEF ( IPRM1 )                                              cnr
      GO TO 90                                                            210   
C
C     IGOTO = 14:  ; AFTER ALL CHANNELS ARE SET UP:
C     FULL TRANSFER OR INELASTIC CALCULATION
C
 850  ONESW = .FALSE.
      GO TO 600
C
C     IGOTO = 15:  ; AFTER CROSSSECTION IF LINTRPOLATION IS NOT DONE:
C     BOTH LINTRP AND XSECTN
C                                                                         220   
 870  ONESW = .FALSE.
      GO TO 650
C
CRA 100  CONTINUE
CRA 200  CONTINUE
CRA 300  CONTINUE
CRA 700  CONTINUE
CRA 720  CONTINUE
CRA 740  CONTINUE
CRA 800  write (6, 903)                                                   230   
CRA 903  FORMAT ('0*** REQUESTED PROGRAM NOT AVAILABLE.' )
CRA      IPRM1 = 0
CRA      GO TO 90
C
      END
