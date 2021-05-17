     SUBROUTINE Nelec_FITATM(iatom, Nel_FIT)
!
!     Purpose: Calculate the number of electrons using the fitted
!                coefficients in 'iatom'
!
!     History: - Creation (25.06.12, MJR)
!
!     ******************************************************************
!     List of OUTPUT variables:
!     Nel_FIT                 : Normalization coefficient

!
!     List of local variables:

!
!     List of USED variables from EMBPARAM:

!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER                 :: ISET, iatom

      REAL (8)                :: NORMA, ZETA, Nel_FIT, NORM_S

!     Starting the variables
      Nel_FIT = 0.0D0

!     Loop over each ZET()
      DO ISET=LAUXPTR(iatom), UAUXPTR(iatom)

!        Getting the NORMA and ZETA values
         NORMA = EMBNORMA(LAUXSHL(ISET))
         ZETA = ZET(ISET)

!        If a /= s, <a> = 0
         IF (LTOTAL(ISET) == 0) THEN

            ! Calculating the integral
            CALL NORMALIZATION_S(NORMA, ZETA, NORM_S)

            ! Calculating the number of electrons for fitted coefficients
            Nel_FIT = Nel_FIT + EMB_COEF( LAUXSHL(ISET) ,1) * NORM_S

         ENDIF
      ENDDO




!
!     *** End of SUBROUTINE FITERROR ***
!
      END SUBROUTINE Nelec_FITATM
