      SUBROUTINE OVERLAP_PS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_PS)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [pi_A|s_B]
!
!         ------------------------
!           i    AX, AY, AZ   pi
!         ------------------------
!           1    1   0   0    px
!           2    0   1   0    py
!           3    0   0   1    pz
!         ------------------------
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
!     OVER_PS(1,1) = [px_A|s_B]
!     OVER_PS(2,1) = [py_A|s_B]
!     OVER_PS(3,1) = [pz_A|s_B]
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
      REAL(8) :: OVER_PS(3,1)

!     List of local variables:
      REAL(8) :: ZETAB, OVER_SS(1,1)
      REAL(8) :: DIFX, DIFY, DIFZ
      REAL(8) :: FUNC_P1
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
      OVER_PS(1,1) = FUNC_P1(ZETAB, DIFX) ! [px_A|s_B]

      OVER_PS(2,1) = FUNC_P1(ZETAB, DIFY) ! [py_A|s_B]

      OVER_PS(3,1) = FUNC_P1(ZETAB, DIFZ) ! [pz_A|s_B]

!     Multiplying everything by OVER_SS
      OVER_PS = OVER_PS*OVER_SS(1,1)
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_PS ***

      END SUBROUTINE OVERLAP_PS


      REAL(8) FUNCTION FUNC_P1(ZET, DIF1)
      REAL(8) :: ZET, DIF1

      FUNC_P1 = 2.0D0*ZET*DIF1

      END FUNCTION FUNC_P1

