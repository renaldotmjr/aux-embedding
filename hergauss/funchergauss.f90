
!     Purpose: Hermite Gaussian functions <s_A|  <p_A|  <d_A|
!
!     History: - Creation (14.02.12, MJR)
!
!     ******************************************************************
!     List of FUNCTIONS:
!
!     HERGA_S(NA, ZETA, r_eA2)                : s-type Hermite-Gaussian
!     HERGA_P(NA, ZETA, DIF, r_eA2)           : p-type Hermite-Gaussian
!     HERGA_D1(NA, ZETA, DIF, r_eA2)          : d-type Hermite-Gaussian
!     HERGA_D2(NA, ZETA, DIF1, DIF2, r_eA2)   : d-type Hermite-Gaussian
!
!     ******************************************************************
!

!
!     ********* Definition of <s_A| function *********
      REAL(8) FUNCTION HERGA_S(NA, ZETA, r_eA2)
      REAL(8) :: NA, ZETA, r_eA2

!     <s_A| = Na*exp( -ZETA*(r-A)^2 )
      HERGA_S = NA * exp( -ZETA*r_eA2 )

      END FUNCTION HERGA_S
!     ************************************************

!
!     ********* Definition of <p_A| function *********
      REAL(8) FUNCTION HERGA_P(NA, ZETA, DIF, r_eA2)
      REAL(8) :: NA, ZETA, DIF, r_eA2

!     <p_A| = -2*Na_pi*ZETA*(re-rA)*exp( -ZETA*(r-A)^2 )
      HERGA_P = -2.0D0*NA*ZETA * DIF * exp( -ZETA*r_eA2 )

      END FUNCTION HERGA_P
!     ************************************************

!
!     ******************** Definition of <d_A| function 1 ********************
      REAL(8) FUNCTION HERGA_D1(NA, ZETA, DIF, r_eA2)
      REAL(8) :: NA, ZETA, DIF, r_eA2

!   1 <d_A| = [4*Na_di*ZETA^2*(re-rA)^2 - 2*Na_di*ZETA]*exp( -ZETA*(r-A)^2 )
      HERGA_D1 = (4.0D0*NA*(ZETA**2) * (DIF**2) - 2*NA*ZETA ) * exp( -ZETA*r_eA2 )

      END FUNCTION HERGA_D1
!     ************************************************************************

!
!     ******************** Definition of <d_A| function 2 ********************
      REAL(8) FUNCTION HERGA_D2(NA, ZETA, DIF1, DIF2, r_eA2)
      REAL(8) :: NA, ZETA, DIF1, DIF2, r_eA2

!   2 <d_A| = 4*Na_di*ZETA^2*(rei-rAi)*(rej-rAj)*exp( -ZETA*(r-A)^2 )
      HERGA_D2 = 4.0D0*NA*(ZETA**2) * DIF1*DIF2 * exp( -ZETA*r_eA2 )

      END FUNCTION HERGA_D2
!     ************************************************************************
