     SUBROUTINE SIMPSON(XY, NP, INT_SIMP)
!
!     Purpose: Integrates a function using the Simpson's method with 5 points
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

      DO I=1, (NP - 2) , 2

          INT_SIMP = INT_SIMP + XY(I) + 4.0D0*XY(I+1) + XY(I+2)

      ENDDO

      INT_SIMP = INT_SIMP/3.0D0

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE SIMPSON
