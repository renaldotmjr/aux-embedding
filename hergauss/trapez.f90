     SUBROUTINE TRAPEZ(XY, NP, INT_TRAP)
!
!     Purpose: Integrates a function using the Trapezoid's method
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
      REAL (8)      :: XY(NP), INT_TRAP


      INT_TRAP = 0.0D0

      DO I=1, (NP - 1)

          INT_TRAP = INT_TRAP + XY(I) + XY(I+1)

      ENDDO

      INT_TRAP = 0.5D0*INT_TRAP

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE TRAPEZ
