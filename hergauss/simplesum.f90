     SUBROUTINE SIMPLESUM(XY, NP, INT_SS)
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
      REAL (8)      :: XY(NP), INT_SS


      INT_SS = 0.0D0

      DO I=1, NP

          INT_SS = INT_SS + XY(I)

      ENDDO

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE SIMPLESUM
