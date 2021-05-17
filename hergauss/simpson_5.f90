     SUBROUTINE SIMPSON_5(XY, NP, INT_SIMP)
!
!     Purpose: Integrates a function using the Simpson's method
!              *** Without the line element ***
!
!     History: - Creation (16.04.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:

!
!     ******************************************************************
!

      IMPLICIT NONE

      INTEGER       :: I, NP
      REAL (8)      :: XY(NP), INT_SIMP


      INT_SIMP = 0.0D0

      DO I=1, (NP - 4) , 4

          INT_SIMP = INT_SIMP + 14.0D0*XY(I) + 64.0D0*XY(I+1) + 24.0D0*XY(I+2) + &
                   & 64.0D0*XY(I+3) + 14.0D0*XY(I+4)

      ENDDO

      INT_SIMP = INT_SIMP/45.0D0

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE SIMPSON_5
