      SUBROUTINE OVERLAP_SP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_SP)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [s_A|pi_B]
!
!      ---------------------------------------------
!        i    AX, AY, AZ | BX, BY, BZ |     pi
!      ---------------------------------------------
!        1    0   0   0  | 1   0   0  | [s_A|px_B]
!        2    0   0   0  | 0   1   0  | [s_A|py_B]
!        3    0   0   0  | 0   0   1  | [s_A|pz_B]
!      ---------------------------------------------
!
!     AX, AY, AZ: Angular momentum index of auxiliary function centered
!                 in atom A.
!
!     BX, BY, BZ: Angular momentum index of auxiliary function centered
!                 in atom B.
!
!     History: - Creation (02.09.11, MJR)
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
!     OVER_SP(1,1) = [s_A|px_B]
!     OVER_SP(1,2) = [s_A|py_B]
!     OVER_SP(1,3) = [s_A|pz_B]
!
!     ******************************************************************
!
      IMPLICIT NONE

!
!     ------------------------------------------------------------------
!
!     List of input variables:
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB
!
!     List of output variable:
      REAL(8) :: OVER_SP(1,3)
!
!     List of local variables:
      REAL(8) :: OVER_PS(3,1)
!
!     ------------------------------------------------------------------
!     This integrals uses the [pi_A|s_B] integrals
!
      CALL OVERLAP_PS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_PS)
!
!     ------------------------------------------------------------------
!
      OVER_SP(1,1) = -OVER_PS(1,1)     ! [s_A|px_B] = -[px_A|s_B]
      OVER_SP(1,2) = -OVER_PS(2,1)     ! [s_A|py_B] = -[py_A|s_B]
      OVER_SP(1,3) = -OVER_PS(3,1)     ! [s_A|pz_B] = -[pz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_SP ***
!
      END SUBROUTINE OVERLAP_SP


