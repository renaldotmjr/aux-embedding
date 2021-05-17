      SUBROUTINE OVERLAP_FS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_FS)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [fi_A|s_B]
!
!         --------------------------
!           i    AX, AY, AZ    fi
!         --------------------------
!           1    3   0   0    f_xxx
!           2    2   1   0    f_xxy
!           3    2   0   1    f_xxz
!           4    1   2   0    f_xyy
!           5    1   1   1    f_xyz
!           6    1   0   2    f_xzz
!           7    0   3   0    f_yyy
!           8    0   2   1    f_yyz
!           9    0   1   2    f_yzz
!          10    0   0   3    f_zzz
!         --------------------------
!
!     AX, AY, AZ: Angular momentum index of auxiliary function.
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
!     OVER_FS( 1,1) = [f_xxx_A|s_B]
!     OVER_FS( 2,1) = [f_xxy_A|s_B]
!     OVER_FS( 3,1) = [f_xxz_A|s_B]
!     OVER_FS( 4,1) = [f_xyy_A|s_B]
!     OVER_FS( 5,1) = [f_xyz_A|s_B]
!     OVER_FS( 6,1) = [f_xzz_A|s_B]
!     OVER_FS( 7,1) = [f_yyy_A|s_B]
!     OVER_FS( 8,1) = [f_yyz_A|s_B]
!     OVER_FS( 9,1) = [f_yzz_A|s_B]
!     OVER_FS(10,1) = [f_zzz_A|s_B]
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
      REAL(8) :: OVER_FS(10,1)

!     List of local variables:
      REAL(8) :: ZETAB, FACTOR1, FACTOR2, OVER_SS(1,1)
      REAL(8) :: DIFX, DIFY, DIFZ
      REAL(8) :: FUNC_F1, FUNC_F2, FUNC_F3
      INTEGER :: I
!
!     ------------------------------------------------------------------
!
      ZETAB = (ZETA*ZETB)/(ZETA+ZETB)

      DIFX = XB - XA
      DIFY = YB - YA
      DIFZ = ZB - ZA

      CALL OVERLAP_SS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_SS)
!
!     ------------------------------------------------------------------
!
      OVER_FS( 1,1) = FUNC_F1(ZETAB, DIFX)               ! [f_xxx_A|s_B]

      OVER_FS( 2,1) = FUNC_F2(ZETAB, DIFX, DIFY)         ! [f_xxy_A|s_B]

      OVER_FS( 3,1) = FUNC_F2(ZETAB, DIFX, DIFZ)         ! [f_xxz_A|s_B]

      OVER_FS( 4,1) = FUNC_F2(ZETAB, DIFY, DIFX)         ! [f_xyy_A|s_B]

      OVER_FS( 5,1) = FUNC_F3(ZETAB, DIFX, DIFY, DIFZ)   ! [f_xyz_A|s_B]

      OVER_FS( 6,1) = FUNC_F2(ZETAB, DIFZ, DIFX)         ! [f_xzz_A|s_B]

      OVER_FS( 7,1) = FUNC_F1(ZETAB, DIFY)               ! [f_yyy_A|s_B]

      OVER_FS( 8,1) = FUNC_F2(ZETAB, DIFY, DIFZ)         ! [f_yyz_A|s_B]

      OVER_FS( 9,1) = FUNC_F2(ZETAB, DIFZ, DIFY)         ! [f_yzz_A|s_B]

      OVER_FS(10,1) = FUNC_F1(ZETAB, DIFZ)               ! [f_zzz_A|s_B]

!     Multiplying everything by OVER_SS
      OVER_FS = OVER_FS*OVER_SS(1,1)
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_FS ***
!
      END SUBROUTINE OVERLAP_FS



      REAL(8) FUNCTION FUNC_F1(ZET, DIF1)
      REAL(8) :: ZET, DIF1
!     FUNC_F1 => f_iii , i=x,y,z

         FUNC_F1 = 8.0D0*(ZET**3)*(DIF1**3) - 12.0D0*(ZET**2)*DIF1
!
      END FUNCTION FUNC_F1



      REAL(8) FUNCTION FUNC_F2(ZET, DIF1, DIF2)
      REAL(8) :: ZET, DIF1, DIF2
!     FUNC_F2 => f_iij , i=x,y,z , j=x,y,z , i.NE.j

         FUNC_F2 = 8.0D0*(ZET**3)*(DIF1**2)*DIF2 - 4.0D0*(ZET**2)*DIF2
!
      END FUNCTION FUNC_F2



      REAL(8) FUNCTION FUNC_F3(ZET, DIF1, DIF2, DIF3)
      REAL(8) :: ZET, DIF1, DIF2, DIF3
!     FUNC_F3 => f_ijk , i=x,y,z , j=x,y,z , k=x,y,z , i.NE.j.NE.k

         FUNC_F3 = 8.0D0*(ZET**3)*DIF1*DIF2*DIF3

      END FUNCTION FUNC_F3

