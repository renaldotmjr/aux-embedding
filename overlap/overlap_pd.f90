      SUBROUTINE OVERLAP_PD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_PD)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [pi_A|dj_B]
!
!---------------------------------------------
!  i/j  AX, AY, AZ |  BX, BY, BZ | BX, BY, BZ | BX, BY, BZ | BX, BY, BZ | BX, BY, BZ | BX, BY, BZ  |
!                  |  2   0   0  | 1   1   0  | 1   0   1  | 0   2   0  | 0   1   1  |  0   0   2  |
!--------------------------------------------------------------------------------------------------|
!  1    1   0   0  | [px_A|dxx_B] [px_A|dxy_B] [px_A|dxz_B] [px_A|dyy_B] [px_A|dyz_B] [px_A|dzz_B] |
!  2    0   1   0  | [py_A|dxx_B] [py_A|dxy_B] [py_A|dxz_B] [py_A|dyy_B] [py_A|dyz_B] [py_A|dzz_B] |
!  3    0   0   1  | [pz_A|dxx_B] [pz_A|dxy_B] [pz_A|dxz_B] [pz_A|dyy_B] [pz_A|dyz_B] [pz_A|dzz_B] |
!--------------------------------------------------------------------------------------------------|
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
!     OVER_PD(1,1) = [px_A|dxx_B]
!     OVER_PD(1,2) = [px_A|dxy_B]
!     OVER_PD(1,3) = [px_A|dxz_B]
!     OVER_PD(1,4) = [px_A|dyy_B]
!     OVER_PD(1,5) = [px_A|dyz_B]
!     OVER_PD(1,6) = [px_A|dzz_B]
!     OVER_PD(2,1) = [py_A|dxx_B]
!     OVER_PD(2,2) = [py_A|dxy_B]
!     OVER_PD(2,3) = [py_A|dxz_B]
!     OVER_PD(2,4) = [py_A|dyy_B]
!     OVER_PD(2,5) = [py_A|dyz_B]
!     OVER_PD(2,6) = [py_A|dzz_B]
!     OVER_PD(3,1) = [pz_A|dxx_B]
!     OVER_PD(3,2) = [pz_A|dxy_B]
!     OVER_PD(3,3) = [pz_A|dxz_B]
!     OVER_PD(3,4) = [pz_A|dyy_B]
!     OVER_PD(3,5) = [pz_A|dyz_B]
!     OVER_PD(3,6) = [pz_A|dzz_B]
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
      REAL(8) :: OVER_PD(3,6)

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
      OVER_PD(1,1) = OVER_FS(1,1)    ! [px_A|dxx_B] = [fxxx_A|s_B]
      OVER_PD(1,2) = OVER_FS(2,1)    ! [px_A|dxy_B] = [fxxy_A|s_B]
      OVER_PD(1,3) = OVER_FS(3,1)    ! [px_A|dxz_B] = [fxxz_A|s_B]
      OVER_PD(1,4) = OVER_FS(4,1)    ! [px_A|dyy_B] = [fxyy_A|s_B]
      OVER_PD(1,5) = OVER_FS(5,1)    ! [px_A|dyz_B] = [fxyz_A|s_B]
      OVER_PD(1,6) = OVER_FS(6,1)    ! [px_A|dzz_B] = [fxzz_A|s_B]
      OVER_PD(2,1) = OVER_FS(2,1)    ! [py_A|dxx_B] = [fxxy_A|s_B]
      OVER_PD(2,2) = OVER_FS(4,1)    ! [py_A|dxy_B] = [fxyy_A|s_B]
      OVER_PD(2,3) = OVER_FS(5,1)    ! [py_A|dxz_B] = [fxyz_A|s_B]
      OVER_PD(2,4) = OVER_FS(7,1)    ! [py_A|dyy_B] = [fyyy_A|s_B]
      OVER_PD(2,5) = OVER_FS(8,1)    ! [py_A|dyz_B] = [fyyz_A|s_B]
      OVER_PD(2,6) = OVER_FS(9,1)    ! [py_A|dzz_B] = [fyzz_A|s_B]
      OVER_PD(3,1) = OVER_FS(3,1)    ! [pz_A|dxx_B] = [fxxz_A|s_B]
      OVER_PD(3,2) = OVER_FS(5,1)    ! [pz_A|dxy_B] = [fxyz_A|s_B]
      OVER_PD(3,3) = OVER_FS(6,1)    ! [pz_A|dxz_B] = [fxzz_A|s_B]
      OVER_PD(3,4) = OVER_FS(8,1)    ! [pz_A|dyy_B] = [fyyz_A|s_B]
      OVER_PD(3,5) = OVER_FS(9,1)    ! [pz_A|dyz_B] = [fyzz_A|s_B]
      OVER_PD(3,6) = OVER_FS(10,1)   ! [pz_A|dzz_B] = [fzzz_A|s_B]
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_PD ***
!
      END SUBROUTINE OVERLAP_PD


