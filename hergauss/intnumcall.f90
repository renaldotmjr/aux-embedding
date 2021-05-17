     SUBROUTINE INTNUMCALL(XY, NP, INT_NUM)
!
!     Purpose: Call the numerical integration subroutines
!
!     History: - Creation (16.04.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:

!
!     ******************************************************************
!

      IMPLICIT NONE

      INTEGER       :: NP
      REAL (8)      :: XY(NP), INT_NUM


      INT_NUM = 0.0D0

      !CALL SIMPSON(XY, NP, INT_NUM)
      !CALL SIMPSON_5(XY, NP, INT_NUM)
      !CALL TRAPEZ(XY, NP, INT_NUM)
      CALL SIMPLESUM(XY, NP, INT_NUM)

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE INTNUMCALL
