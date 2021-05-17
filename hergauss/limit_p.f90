     SUBROUTINE LIMIT_P(Na, ZETA, TOL, D)
!
!     Purpose: Find the index in grid nearest the atom
!
!     History: - Creation (19.10.11, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     Na           : Normalization coefficient
!     ZETA         : basis exponents
!
!     List of OUTPUT variables
!     D            : The Limit of numerical integration
!
!     ******************************************************************
!
      USE EMBTOL

      IMPLICIT NONE

      INTEGER       :: I, IMAX

      REAL (8)      :: Na, ZETA, D, Di, TOL


!     Initial values of I, Di and D
      Di = 10

!     Testing if TOL is in the existence condition
      IF ( TOL < 2*Na*ZETA*Di ) THEN

!        Calculating D iterativelly
         do I = 0, IMAX

            D = SQRT( (-1.0D0/ZETA)*DLOG( TOL/(2*Na*ZETA*Di) ) )

            ! Exit if XTOL tolerance was achieved
            IF ( DABS(D - Di) <= XTOL ) EXIT

            ! Passing D to Di
            Di = D
         enddo

      ELSE

         write(6, *) 'Low numerical integration tolerance!'
         write(6, *) 'Tolerance must be .LT. ', 2*Na*ZETA*Di
         STOP

      ENDIF

!
!     *** End of SUBROUTINE LIMIT_P ***
!
      END SUBROUTINE LIMIT_P
