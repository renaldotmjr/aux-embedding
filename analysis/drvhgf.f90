     SUBROUTINE DRVHGF(LTOTSHL, COEF, ZETA, NA, xe, ye, ze, r_eA2, fit_rho)
!
!     Purpose: Calculate RHO in 'ri' position with Hermite-Gaussians
!               centered on atoms
!
!     History: - Creation (04.06.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:

!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     List of Local variables
      INTEGER       :: LTOTSHL
      REAL (8)      :: fit_rho, ZETA, HG_f, xe, ye, ze, r_eA2
      REAL (8)      :: HERGA_S, HERGA_P, HERGA_D1, HERGA_D2

!     List of local ALLOCATABLE variables
      ! In this case, NA will not be allocated for optimize this calculation
      REAL (8), DIMENSION(14) :: NA, COEF  ! Normalization coef.


      ! ##################  S-TYPE block ##################
      IF ( LTOTSHL == 0 ) THEN
          ! Calculating the Hermite-Gaussian function
          HG_f = COEF(1) * HERGA_S(NA(1), ZETA, r_eA2)
      ! ############### END of S-TYPE block ###############

      ! ##################  P-TYPE block ##################
      ELSEIF ( LTOTSHL == 1 ) THEN
          ! Calculating the Hermite-Gaussian function
          HG_f = COEF(1) * HERGA_P(NA(1), ZETA, xe, r_eA2) + &
               & COEF(2) * HERGA_P(NA(2), ZETA, ye, r_eA2) + &
               & COEF(3) * HERGA_P(NA(3), ZETA, ze, r_eA2)
          ! ############### END of P-TYPE block ###############

          ! ##################  D-TYPE block ##################
      ELSEIF ( LTOTSHL == 2 ) THEN
          ! Calculating the Hermite-Gaussian function
          HG_f = COEF(1) * HERGA_D1(NA(1), ZETA, xe, r_eA2) + &
               & COEF(2) * HERGA_D2(NA(2), ZETA, xe, ye, r_eA2) + &
               & COEF(3) * HERGA_D2(NA(3), ZETA, xe, ze, r_eA2) + &
               & COEF(4) * HERGA_D1(NA(4), ZETA, ye, r_eA2) + &
               & COEF(5) * HERGA_D2(NA(5), ZETA, ye, ze, r_eA2) + &
               & COEF(6) * HERGA_D1(NA(6), ZETA, ze, r_eA2)

      ENDIF
      ! ############### END of D-TYPE integrals block ###############

      fit_rho = HG_f

!
!     *** End of SUBROUTINE DRVHGF ***
!
      END SUBROUTINE DRVHGF
