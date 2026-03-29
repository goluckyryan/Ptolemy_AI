      PROGRAM DUMPMT
      IMPLICIT REAL*8 (A-H,O-Z)
      
      DO 10 IZ = 0, 119
        DO 10 IA = 1, 300
          CALL MASEXX(IZ, IA, EX, ER, IERR)
          IF (IERR .EQ. 0) THEN
            WRITE(6,'(I4,I5,F15.6,F12.6)') IZ, IA, EX, ER
          ENDIF
   10 CONTINUE
      END
