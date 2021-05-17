      SUBROUTINE OVERLAP_DP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_DP)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [di_A|pj_B]
!
!      -----------------------------------------------------------
!       i/j   AX, AY, AZ |  BX, BY, BZ | BX, BY, BZ |  BX, BY, BZ |
!      -----------------------------------------------------------|
!                        |  1   0   0  | 0   1   0  |  0   0   1  |
!      -----------------------------------------------------------|
!        1    2   0   0  | [dxx_A|px_B] [dxx_A|py_B] [dxx_A|pz_B] |
!        2    1   1   0  | [dxy_A|px_B] [dxy_A|py_B] [dxy_A|pz_B] |
!        3    1   0   1  | [dxz_A|px_B] [dxz_A|py_B] [dxz_A|pz_B] |
!        4    0   2   0  | [dyy_A|px_B] [dyy_A|py_B] [dyy_A|pz_B] |
!        5    0   1   1  | [dyz_A|px_B] [dyz_A|py_B] [dyz_A|pz_B] |
!        6    0   0   2  | [dzz_A|px_B] [dzz_A|py_B] [dzz_A|pz_B] |
!      -----------------------------------------------------------|
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
!     OVER_DP(1,1) = [dxx_A|px_B]
!     OVER_DP(1,2) = [dxx_A|py_B]
!     OVER_DP(1,3) = [dxx_A|pz_B]
!     OVER_DP(2,1) = [dxy_A|px_B]
!     OVER_DP(2,2) = [dxy_A|py_B]
!     OVER_DP(2,3) = [dxy_A|pz_B]
!     OVER_DP(3,1) = [dxz_A|px_B]
!     OVER_DP(3,2) = [dxz_A|py_B]
!     OVER_DP(3,3) = [dxz_A|pz_B]
!     OVER_DP(4,1) = [dyy_A|px_B]
!     OVER_DP(4,2) = [dyy_A|py_B]
!     OVER_DP(4,3) = [dyy_A|pz_B]
!     OVER_DP(5,1) = [dyz_A|px_B]
!     OVER_DP(5,2) = [dyz_A|py_B]
!     OVER_DP(5,3) = [dyz_A|pz_B]
!     OVER_DP(6,1) = [dzz_A|px_B]
!     OVER_DP(6,2) = [dzz_A|py_B]
!     OVER_DP(6,3) = [dzz_A|pz_B]
!
!     ******************************************************************
!
      IMPLICIT NONE
!
!     ------------------------------------------------------------------
!
!     List of input variables:
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB

!     List of output variable:
      REAL(8) :: OVER_DP(6,3)

!     List of local variables:
      REAL(8) :: OVER_FS(10,1)
!
!     ------------------------------------------------------------------
!     This integrals uses the [fi_A|s_B] integrals
!
      CALL OVERLAP_FS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_FS)
!
!     ------------------------------------------------------------------
!
      OVER_DP(1,1) = -OVER_FS(1,1)    ! [dxx_A|px_B] = -[fxxx_A|s_B]
      OVER_DP(1,2) = -OVER_FS(2,1)    ! [dxx_A|py_B] = -[fxxy_A|s_B]
      OVER_DP(1,3) = -OVER_FS(3,1)    ! [dxx_A|pz_B] = -[fxxz_A|s_B]
      OVER_DP(2,1) = -OVER_FS(2,1)    ! [dxy_A|px_B] = -[fxxy_A|s_B]
      OVER_DP(2,2) = -OVER_FS(4,1)    ! [dxy_A|py_B] = -[fxyy_A|s_B]
      OVER_DP(2,3) = -OVER_FS(5,1)    ! [dxy_A|pz_B] = -[fxyz_A|s_B]
      OVER_DP(3,1) = -OVER_FS(3,1)    ! [dxz_A|px_B] = -[fxxz_A|s_B]
      OVER_DP(3,2) = -OVER_FS(5,1)    ! [dxz_A|py_B] = -[fxyz_A|s_B]
      OVER_DP(3,3) = -OVER_FS(6,1)    ! [dxz_A|pz_B] = -[fxzz_A|s_B]
      OVER_DP(4,1) = -OVER_FS(4,1)    ! [dyy_A|px_B] = -[fxyy_A|s_B]
      OVER_DP(4,2) = -OVER_FS(7,1)    ! [dyy_A|py_B] = -[fyyy_A|s_B]
      OVER_DP(4,3) = -OVER_FS(8,1)    ! [dyy_A|pz_B] = -[fyyz_A|s_B]
      OVER_DP(5,1) = -OVER_FS(5,1)    ! [dyz_A|px_B] = -[fxyz_A|s_B]
      OVER_DP(5,2) = -OVER_FS(8,1)    ! [dyz_A|py_B] = -[fyyz_A|s_B]
      OVER_DP(5,3) = -OVER_FS(9,1)    ! [dyz_A|pz_B] = -[fyzz_A|s_B]
      OVER_DP(6,1) = -OVER_FS(6,1)    ! [dzz_A|px_B] = -[fxzz_A|s_B]
      OVER_DP(6,2) = -OVER_FS(9,1)    ! [dzz_A|py_B] = -[fyzz_A|s_B]
      OVER_DP(6,3) = -OVER_FS(10,1)   ! [dzz_A|pz_B] = -[fzzz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_DP***
!
      END SUBROUTINE OVERLAP_DP


