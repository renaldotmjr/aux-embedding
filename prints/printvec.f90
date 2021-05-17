     SUBROUTINE PRINTVEC(FILE_I, VEC, VSIZE, N_col)
!
!     Purpose: Print the vector VEC(VSIZE) in FILE_NAME with N_col
!               number of colunms
!
!     History: - Creation (11.06.12, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!
!     List of local variables:
!
!     ******************************************************************

      IMPLICIT NONE

      ! List of local variables
      INTEGER        :: i, VSIZE, N_col, FILE_I

      CHARACTER (80) :: INPUT_NAME

      REAL (8)        :: VEC(VSIZE)


      DO i=1, VSIZE

          IF ( i .EQ. VSIZE) THEN

              WRITE (FILE_I,'(F18.6)') VEC(i)

          ELSE

              IF ( MOD(i,N_col) .NE. 0 ) WRITE (FILE_I,'(F18.6,$)') VEC(i)
              IF ( MOD(i,N_col) .EQ. 0 ) WRITE (FILE_I,'(F18.6)') VEC(i)

          ENDIF

      ENDDO

!
!     *** End of SUBROUTINE PRINTVEC ***
!
      END SUBROUTINE PRINTVEC
