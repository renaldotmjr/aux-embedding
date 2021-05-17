      SUBROUTINE OVERLAP_DS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_DS)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [di_A|s_B]
!
!         --------------------------
!           i    AX, AY, AZ    di
!         --------------------------
!           1    2   0   0    d_xx
!           2    1   1   0    d_xy
!           3    1   0   1    d_xz
!           4    0   2   0    d_yy
!           5    0   1   1    d_yz
!           6    0   0   2    d_zz
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
!     OVER_DS(1,1) = [d_xx_A|s_B]
!     OVER_DS(2,1) = [d_xy_A|s_B]
!     OVER_DS(3,1) = [d_xz_A|s_B]
!     OVER_DS(4,1) = [d_yy_A|s_B]
!     OVER_DS(5,1) = [d_yz_A|s_B]
!     OVER_DS(6,1) = [d_zz_A|s_B]
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
      REAL(8) :: OVER_DS(6,1)

!     List of local variables:
      REAL(8) :: ZETAB, OVER_SS(1,1)
      REAL(8) :: DIFX, DIFY, DIFZ
      REAL(8) :: FUNC_D1, FUNC_D2
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
      OVER_DS(1,1) = FUNC_D1(ZETAB, DIFX)        ! [d_xx_A|s_B]

      OVER_DS(2,1) = FUNC_D2(ZETAB, DIFX, DIFY)  ! [d_xy_A|s_B]

      OVER_DS(3,1) = FUNC_D2(ZETAB, DIFX, DIFZ)  ! [d_xz_A|s_B]

      OVER_DS(4,1) = FUNC_D1(ZETAB, DIFY)        ! [d_yy_A|s_B]

      OVER_DS(5,1) = FUNC_D2(ZETAB, DIFY, DIFZ)  ! [d_yz_A|s_B]

      OVER_DS(6,1) = FUNC_D1(ZETAB, DIFZ)        ! [d_zz_A|s_B]

!     Multiplying everything by OVER_SS
      OVER_DS = OVER_DS*OVER_SS(1,1)

!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_DS ***
!
      END SUBROUTINE OVERLAP_DS



      REAL(8) FUNCTION FUNC_D1(ZET, DIF1)
      REAL(8) :: ZET, DIF1
!     FUNC_D1 => d_ii , i=x,y,z

         FUNC_D1 = 4.0D0*(ZET**2)*(DIF1**2) - 2.0D0*ZET

      END FUNCTION FUNC_D1



      REAL(8) FUNCTION FUNC_D2(ZET, DIF1, DIF2)
      REAL(8) :: ZET, DIF1, DIF2
!     FUNC_D2 => d_ij , i=x,y,z , j=x,y,z , i.NE.j

         FUNC_D2 = 4.0D0*(ZET**2)*(DIF1*DIF2)

      END FUNCTION FUNC_D2

