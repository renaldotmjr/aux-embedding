      SUBROUTINE OVERLAP_SD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_SD)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [s_A|di_B]
!
!      ---------------------------------------------
!        i    AX, AY, AZ | BX, BY, BZ |     pi
!      ---------------------------------------------
!        1    0   0   0  | 2   0   0  | [s_A|dxx_B]
!        2    0   0   0  | 1   1   0  | [s_A|dxy_B]
!        3    0   0   0  | 1   0   1  | [s_A|dxz_B]
!        4    0   0   0  | 0   2   0  | [s_A|dyy_B]
!        5    0   0   0  | 0   1   1  | [s_A|dyz_B]
!        6    0   0   0  | 0   0   2  | [s_A|dzz_B]
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
!     OVER_SD(1,1) = [s_A|dxx_B]
!     OVER_SD(1,2) = [s_A|dxy_B]
!     OVER_SD(1,3) = [s_A|dxz_B]
!     OVER_SD(1,4) = [s_A|dyy_B]
!     OVER_SD(1,5) = [s_A|dyz_B]
!     OVER_SD(1,6) = [s_A|dzz_B]
!
!     ******************************************************************
!
      IMPLICIT NONE

!     ------------------------------------------------------------------

!     List of input variables:
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB

!     List of output variable:
      REAL(8) :: OVER_SD(1,6)

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
      OVER_SD(1,1) = OVER_DS(1,1)     ! [s_A|dxx_B] = [dxx_A|s_B]
      OVER_SD(1,2) = OVER_DS(2,1)     ! [s_A|dxy_B] = [dxy_A|s_B]
      OVER_SD(1,3) = OVER_DS(3,1)     ! [s_A|dxz_B] = [dxz_A|s_B]
      OVER_SD(1,4) = OVER_DS(4,1)     ! [s_A|dyy_B] = [dxy_A|s_B]
      OVER_SD(1,5) = OVER_DS(5,1)     ! [s_A|dyz_B] = [dyy_A|s_B]
      OVER_SD(1,6) = OVER_DS(6,1)     ! [s_A|dzz_B] = [dyz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_SD ***
!
      END SUBROUTINE OVERLAP_SD


