      SUBROUTINE OVERLAP_PP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_PP)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [pi_A|pj_B]
!
!      ----------------------------------------------------------
!       i/j   AX, AY, AZ | BX, BY, BZ  | BX, BY, BZ  | BX, BY, BZ |
!      ----------------------------------------------------------
!                        | 1   0   0   | 0   1   0   | 0   0   1  |
!      ----------------------------------------------------------
!        1    1   0   0  | [px_A|px_B] | [px_A|py_B] | [px_A|pz_B]
!        2    0   1   0  | [py_A|px_B] | [py_A|py_B] | [py_A|pz_B]
!        3    0   0   1  | [pz_A|px_B] | [pz_A|py_B] | [pz_A|pz_B]
!      ----------------------------------------------------------
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
!     OVER_PP(1,1) = [px_A|px_B]
!     OVER_PP(1,2) = [px_A|py_B]
!     OVER_PP(1,3) = [px_A|pz_B]
!     OVER_PP(2,1) = [py_A|px_B]
!     OVER_PP(2,2) = [py_A|py_B]
!     OVER_PP(2,3) = [py_A|pz_B]
!     OVER_PP(3,1) = [pz_A|px_B]
!     OVER_PP(3,2) = [pz_A|py_B]
!     OVER_PP(3,3) = [pz_A|pz_B]
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
      REAL(8) :: OVER_PP(3,3)
!
!     List of local variables:
      REAL(8) :: OVER_DS(6,1)
!
!     ------------------------------------------------------------------
!     This integrals uses the [di_A|s_B] integrals
!
      CALL OVERLAP_DS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_DS)
!
!     ------------------------------------------------------------------
!
      OVER_PP(1,1) = -OVER_DS(1,1)     ! [px_A|px_B] = -[dxx_A|s_B]
      OVER_PP(1,2) = -OVER_DS(2,1)     ! [px_A|py_B] = -[dxy_A|s_B]
      OVER_PP(1,3) = -OVER_DS(3,1)     ! [px_A|pz_B] = -[dxz_A|s_B]
      OVER_PP(2,1) = -OVER_DS(2,1)     ! [py_A|px_B] = -[dxy_A|s_B]
      OVER_PP(2,2) = -OVER_DS(4,1)     ! [py_A|py_B] = -[dyy_A|s_B]
      OVER_PP(2,3) = -OVER_DS(5,1)     ! [py_A|pz_B] = -[dyz_A|s_B]
      OVER_PP(3,1) = -OVER_DS(3,1)     ! [pz_A|px_B] = -[dxz_A|s_B]
      OVER_PP(3,2) = -OVER_DS(5,1)     ! [pz_A|py_B] = -[dyz_A|s_B]
      OVER_PP(3,3) = -OVER_DS(6,1)     ! [pz_A|pz_B] = -[dzz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_PP ***
!
      END SUBROUTINE OVERLAP_PP


