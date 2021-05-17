      SUBROUTINE OVERLAP_SS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_SS)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [s_A|s_B]
!
!     History: - Creation (28.06.11, CO)
!
!     ******************************************************************
!
!     List of input variables:
!     ZETA     : Gaussian exponent of auxiliary function A.
!     ZETB     : Gaussian exponent of auxiliary function B.
!     XA,YA,ZA : Atomic coordinates of center A.
!     XB,YB,ZB : Atomic coordinates of center B.
!
!     List of output variable:
!     OVER_SS(1,1) = [s_A|s_B]
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!
!     ------------------------------------------------------------------
!
!     List of input variables:
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB
!
!     List of output variable:
      REAL(8) :: OVER_SS(1,1)
!
!     List of local variables:
      REAL(8) :: ZETAB, FACTOR, rAB2
!
!     ------------------------------------------------------------------
!
      ZETAB = (ZETA*ZETB)/(ZETA+ZETB)
!
      FACTOR = (PI/(ZETA+ZETB))**(1.5)
!
      rAB2 = (XA-XB)**2 + (YA-YB)**2 + (ZA-ZB)**2
!
!     ------------------------------------------------------------------
!
      OVER_SS(1,1) = FACTOR*EXP(-ZETAB*rAB2)   ! [s_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_SS ***
!
      END SUBROUTINE OVERLAP_SS

