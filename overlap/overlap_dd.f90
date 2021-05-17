      SUBROUTINE OVERLAP_DD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_DD)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [di_A|dj_B]
!
!--------------------------------------------------------------------------------------------------------|
!  i/j  AX, AY, AZ |  BX, BY, BZ  |  BX, BY, BZ |  BX, BY, BZ |  BX, BY, BZ |  BX, BY, BZ |  BX, BY, BZ  |
!                  |  2   0   0   |  1   1   0  |  1   0   1  |  0   2   0  |  0   1   1  |   0   0   2  |
!--------------------------------------------------------------------------------------------------------|
!  1    2   0   0  | [dxx_A|dxx_B] [dxx_A|dxy_B] [dxx_A|dxz_B] [dxx_A|dyy_B] [dxx_A|dyz_B] [dxx_A|dzz_B] |
!  2    1   1   0  | [dxy_A|dxx_B] [dxy_A|dxy_B] [dxy_A|dxz_B] [dxy_A|dyy_B] [dxy_A|dyz_B] [dxy_A|dzz_B] |
!  3    1   0   1  | [dxz_A|dxx_B] [dxz_A|dxy_B] [dxz_A|dxz_B] [dxz_A|dyy_B] [dxz_A|dyz_B] [dxz_A|dzz_B] |
!  4    0   2   0  | [dyy_A|dxx_B] [dyy_A|dxy_B] [dyy_A|dxz_B] [dyy_A|dyy_B] [dyy_A|dyz_B] [dyy_A|dzz_B] |
!  5    0   1   1  | [dyz_A|dxx_B] [dyz_A|dxy_B] [dyz_A|dxz_B] [dyz_A|dyy_B] [dyz_A|dyz_B] [dyz_A|dzz_B] |
!  6    0   0   2  | [dzz_A|dxx_B] [dzz_A|dxy_B] [dzz_A|dxz_B] [dzz_A|dyy_B] [dzz_A|dyz_B] [dzz_A|dzz_B] |
!--------------------------------------------------------------------------------------------------------
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
!     OVER_DD(1,1) = [dxx_A|dxx_B]
!     OVER_DD(1,2) = [dxx_A|dxy_B]
!     OVER_DD(1,3) = [dxx_A|dxz_B]
!     OVER_DD(1,4) = [dxx_A|dyy_B]
!     OVER_DD(1,5) = [dxx_A|dyz_B]
!     OVER_DD(1,6) = [dxx_A|dzz_B]
!     OVER_DD(2,1) = [dxy_A|dxx_B]
!     OVER_DD(2,2) = [dxy_A|dxy_B]
!     OVER_DD(2,3) = [dxy_A|dxz_B]
!     OVER_DD(2,4) = [dxy_A|dyy_B]
!     OVER_DD(2,5) = [dxy_A|dyz_B]
!     OVER_DD(2,6) = [dxy_A|dzz_B]
!     OVER_DD(3,1) = [dxz_A|dxx_B]
!     OVER_DD(3,2) = [dxz_A|dxy_B]
!     OVER_DD(3,3) = [dxz_A|dxz_B]
!     OVER_DD(3,4) = [dxz_A|dyy_B]
!     OVER_DD(3,5) = [dxz_A|dyz_B]
!     OVER_DD(3,6) = [dxz_A|dzz_B]
!     OVER_DD(4,1) = [dyy_A|dxx_B]
!     OVER_DD(4,2) = [dyy_A|dxy_B]
!     OVER_DD(4,3) = [dyy_A|dxz_B]
!     OVER_DD(4,4) = [dyy_A|dyy_B]
!     OVER_DD(4,5) = [dyy_A|dyz_B]
!     OVER_DD(4,6) = [dyy_A|dzz_B]
!     OVER_DD(5,1) = [dyz_A|dxx_B]
!     OVER_DD(5,2) = [dyz_A|dxy_B]
!     OVER_DD(5,3) = [dyz_A|dxz_B]
!     OVER_DD(5,4) = [dyz_A|dyy_B]
!     OVER_DD(5,5) = [dyz_A|dyz_B]
!     OVER_DD(5,6) = [dyz_A|dzz_B]
!     OVER_DD(6,1) = [dzz_A|dxx_B]
!     OVER_DD(6,2) = [dzz_A|dxy_B]
!     OVER_DD(6,3) = [dzz_A|dxz_B]
!     OVER_DD(6,4) = [dzz_A|dyy_B]
!     OVER_DD(6,5) = [dzz_A|dyz_B]
!     OVER_DD(6,6) = [dzz_A|dzz_B]
!
!     ******************************************************************
!
      IMPLICIT NONE

!     ------------------------------------------------------------------

!     List of input variables:
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB
!
!     List of output variable:
      REAL(8) :: OVER_DD(6,6)

!     List of local variables:
      REAL(8) :: OVER_GS(15,1)

!     ------------------------------------------------------------------
!     This integrals uses the [gi_A|s_B] integrals

      CALL OVERLAP_GS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_GS)

!     ------------------------------------------------------------------

      OVER_DD(1,1) = OVER_GS(1,1)   ! [dxx_A|dxx_B] = [gxxxx_A|s_B]
      OVER_DD(1,2) = OVER_GS(2,1)   ! [dxx_A|dxy_B] = [gxxxy_A|s_B]
      OVER_DD(1,3) = OVER_GS(3,1)   ! [dxx_A|dxz_B] = [gxxxz_A|s_B]
      OVER_DD(1,4) = OVER_GS(4,1)   ! [dxx_A|dyy_B] = [gxxyy_A|s_B]
      OVER_DD(1,5) = OVER_GS(5,1)   ! [dxx_A|dyz_B] = [gxxyz_A|s_B]
      OVER_DD(1,6) = OVER_GS(6,1)   ! [dxx_A|dzz_B] = [gxxzz_A|s_B]
      OVER_DD(2,1) = OVER_GS(2,1)   ! [dxy_A|dxx_B] = [gxxxy_A|s_B]
      OVER_DD(2,2) = OVER_GS(4,1)   ! [dxy_A|dxy_B] = [gxxyy_A|s_B]
      OVER_DD(2,3) = OVER_GS(5,1)   ! [dxy_A|dxz_B] = [gxxyz_A|s_B]
      OVER_DD(2,4) = OVER_GS(7,1)   ! [dxy_A|dyy_B] = [gxyyy_A|s_B]
      OVER_DD(2,5) = OVER_GS(8,1)   ! [dxy_A|dyz_B] = [gxyyz_A|s_B]
      OVER_DD(2,6) = OVER_GS(9,1)   ! [dxy_A|dzz_B] = [gxyzz_A|s_B]
      OVER_DD(3,1) = OVER_GS(3,1)   ! [dxz_A|dxx_B] = [gxxxz_A|s_B]
      OVER_DD(3,2) = OVER_GS(5,1)   ! [dxz_A|dxy_B] = [gxxyz_A|s_B]
      OVER_DD(3,3) = OVER_GS(6,1)   ! [dxz_A|dxz_B] = [gxxzz_A|s_B]
      OVER_DD(3,4) = OVER_GS(8,1)   ! [dxz_A|dyy_B] = [gxyyz_A|s_B]
      OVER_DD(3,5) = OVER_GS(9,1)   ! [dxz_A|dyz_B] = [gxyzz_A|s_B]
      OVER_DD(3,6) = OVER_GS(10,1)  ! [dxz_A|dzz_B] = [gxzzz_A|s_B]
      OVER_DD(4,1) = OVER_GS(4,1)   ! [dyy_A|dxx_B] = [gxxyy_A|s_B]
      OVER_DD(4,2) = OVER_GS(7,1)   ! [dyy_A|dxy_B] = [gxyyy_A|s_B]
      OVER_DD(4,3) = OVER_GS(8,1)   ! [dyy_A|dxz_B] = [gxyyz_A|s_B]
      OVER_DD(4,4) = OVER_GS(11,1)  ! [dyy_A|dyy_B] = [gyyyy_A|s_B]
      OVER_DD(4,5) = OVER_GS(12,1)  ! [dyy_A|dyz_B] = [gyyyz_A|s_B]
      OVER_DD(4,6) = OVER_GS(13,1)  ! [dyy_A|dzz_B] = [gyyzz_A|s_B]
      OVER_DD(5,1) = OVER_GS(5,1)   ! [dyz_A|dxx_B] = [gxxyz_A|s_B]
      OVER_DD(5,2) = OVER_GS(8,1)   ! [dyz_A|dxy_B] = [gxyyz_A|s_B]
      OVER_DD(5,3) = OVER_GS(9,1)   ! [dyz_A|dxz_B] = [gxyzz_A|s_B]
      OVER_DD(5,4) = OVER_GS(12,1)  ! [dyz_A|dyy_B] = [gyyyz_A|s_B]
      OVER_DD(5,5) = OVER_GS(13,1)  ! [dyz_A|dyz_B] = [gyyzz_A|s_B]
      OVER_DD(5,6) = OVER_GS(14,1)  ! [dyz_A|dzz_B] = [gyzzz_A|s_B]
      OVER_DD(6,1) = OVER_GS(6,1)   ! [dzz_A|dxx_B] = [gxxzz_A|s_B]
      OVER_DD(6,2) = OVER_GS(9,1)   ! [dzz_A|dxy_B] = [gxyzz_A|s_B]
      OVER_DD(6,3) = OVER_GS(10,1)  ! [dzz_A|dxz_B] = [gxzzz_A|s_B]
      OVER_DD(6,4) = OVER_GS(13,1)  ! [dzz_A|dyy_B] = [gyyzz_A|s_B]
      OVER_DD(6,5) = OVER_GS(14,1)  ! [dzz_A|dyz_B] = [gyzzz_A|s_B]
      OVER_DD(6,6) = OVER_GS(15,1)  ! [dzz_A|dzz_B] = [gzzzz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_DP***
!
      END SUBROUTINE OVERLAP_DD
