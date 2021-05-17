      SUBROUTINE OVERLAP_GS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVER_GS)
!
!     Purpose: Calculation of unnormalized primitive overlap integrals
!              [gi_A|s_B]
!
!         ----------------------------
!           i    AX, AY, AZ      gi
!         ----------------------------
!           1    4   0   0    g_xxxx
!           2    3   1   0    g_xxxy
!           3    3   0   1    g_xxxz
!           4    2   2   0    g_xxyy
!           5    2   1   1    g_xxyz
!           6    2   0   2    g_xxzz
!           7    1   3   0    g_xyyy
!           8    1   2   1    g_xyyz
!           9    1   1   2    g_xyzz
!          10    1   0   3    g_xzzz
!          11    0   4   0    g_yyyy
!          12    0   3   1    g_yyyz
!          13    0   2   2    g_yyzz
!          14    0   1   3    g_yzzz
!          15    0   0   4    g_zzzz
!         ----------------------------
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
!     OVER_GS( 1,1) = [g_xxxx_A|s_B]
!     OVER_GS( 2,1) = [g_xxxy_A|s_B]
!     OVER_GS( 3,1) = [g_xxxz_A|s_B]
!     OVER_GS( 4,1) = [g_xxyy_A|s_B]
!     OVER_GS( 5,1) = [g_xxyz_A|s_B]
!     OVER_GS( 6,1) = [g_xxzz_A|s_B]
!     OVER_GS( 7,1) = [g_xyyy_A|s_B]
!     OVER_GS( 8,1) = [g_xyyz_A|s_B]
!     OVER_GS( 9,1) = [g_xyzz_A|s_B]
!     OVER_GS(10,1) = [g_xzzz_A|s_B]
!     OVER_GS(11,1) = [g_yyyy_A|s_B]
!     OVER_GS(12,1) = [g_yyyz_A|s_B]
!     OVER_GS(13,1) = [g_yyzz_A|s_B]
!     OVER_GS(14,1) = [g_yzzz_A|s_B]
!     OVER_GS(15,1) = [g_zzzz_A|s_B]
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
      REAL(8) :: OVER_GS(15,1)
!
!     List of local variables:
      REAL(8) :: ZETAB, FACTOR1, FACTOR2, OVER_SS(1,1)
      REAL(8) :: DIFX, DIFY, DIFZ
      REAL(8) :: FUNC_G1, FUNC_G2, FUNC_G3, FUNC_G4
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
      OVER_GS( 1,1) = FUNC_G1(ZETAB, DIFX)               ! [g_xxxx_A|s_B]

      OVER_GS( 2,1) = FUNC_G2(ZETAB, DIFX, DIFY)         ! [g_xxxy_A|s_B]

      OVER_GS( 3,1) = FUNC_G2(ZETAB, DIFX, DIFZ)         ! [g_xxxz_A|s_B]

      OVER_GS( 4,1) = FUNC_G3(ZETAB, DIFX, DIFY)         ! [g_xxyy_A|s_B]

      OVER_GS( 5,1) = FUNC_G4(ZETAB, DIFX, DIFY, DIFZ)   ! [g_xxyz_A|s_B]

      OVER_GS( 6,1) = FUNC_G3(ZETAB, DIFX, DIFZ)         ! [g_xxzz_A|s_B]

      OVER_GS( 7,1) = FUNC_G2(ZETAB, DIFY, DIFX)         ! [g_xyyy_A|s_B]

      OVER_GS( 8,1) = FUNC_G4(ZETAB, DIFY, DIFX, DIFZ)   ! [g_xyyz_A|s_B]

      OVER_GS( 9,1) = FUNC_G4(ZETAB, DIFZ, DIFX, DIFY)   ! [g_xyzz_A|s_B]

      OVER_GS(10,1) = FUNC_G2(ZETAB, DIFZ, DIFX)         ! [g_xzzz_A|s_B]

      OVER_GS(11,1) = FUNC_G1(ZETAB, DIFY)               ! [g_yyyy_A|s_B]

      OVER_GS(12,1) = FUNC_G2(ZETAB, DIFY, DIFZ)         ! [g_yyyz_A|s_B]

      OVER_GS(13,1) = FUNC_G3(ZETAB, DIFY, DIFZ)         ! [g_yyzz_A|s_B]

      OVER_GS(14,1) = FUNC_G2(ZETAB, DIFZ, DIFY)         ! [g_yzzz_A|s_B]

      OVER_GS(15,1) = FUNC_G1(ZETAB, DIFZ)               ! [g_zzzz_A|s_B]

      OVER_GS = OVER_GS*OVER_SS(1,1)
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_GS ***
!
      END SUBROUTINE OVERLAP_GS


      REAL(8) FUNCTION FUNC_G1(ZET, DIF1)
      REAL(8) :: ZET, DIF1
!     FUNC_G1 => g_iiii , i=x,y,z

         FUNC_G1 = 16.0D0*(ZET**4)*(DIF1**4) - 48.0D0*(ZET**3)*(DIF1**2) + &
         & 12.0D0*(ZET**2)

      END FUNCTION FUNC_G1



      REAL(8) FUNCTION FUNC_G2(ZET, DIF1, DIF2)
      REAL(8) :: ZET, DIF1, DIF2
!     FUNC_G2 => g_iiij , i=x,y,z , j=x,y,z , i.NE.j

         FUNC_G2 = 16.0D0*(ZET**4)*(DIF1**3)*DIF2 - 24.0D0*(ZET**3)*DIF1*DIF2

      END FUNCTION FUNC_G2



      REAL(8) FUNCTION FUNC_G3(ZET, DIF1, DIF2)
      REAL(8) :: ZET, DIF1, DIF2
!     FUNC_G3 => g_iijj , i=x,y,z , j=x,y,z , i.NE.j

         FUNC_G3 = 16.0D0*(ZET**4)*(DIF1**2)*(DIF2**2) - &
         & 8.0D0*(ZET**3)*( (DIF1**2) + (DIF2**2) ) + 4.0D0*(ZET**2)

      END FUNCTION FUNC_G3



      REAL(8) FUNCTION FUNC_G4(ZET, DIF1, DIF2, DIF3)
      REAL(8) :: ZET, DIF1, DIF2, DIF3
!     FUNC_G4 => g_iijk , i=x,y,z , j=x,y,z , k=x,y,z , i.NE.j.NE.k

         FUNC_G4 = 16.0D0*(ZET**4)*(DIF1**2)*DIF2*DIF3 - &
         & 8.0D0*(ZET**3)*DIF2*DIF3

      END FUNCTION FUNC_G4

