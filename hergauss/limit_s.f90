     SUBROUTINE LIMIT_S(Na, ZETA, TOL, D)
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
      IMPLICIT NONE

      REAL (8)      :: Na, ZETA, D, TOL


!     Testing if TOL is in the existence condition
      IF ( TOL < Na ) THEN

!        Calculating D
         D = SQRT( (-1.0D0/ZETA)*DLOG(TOL/Na) )

      ELSE

         write(6, *) 'Low numerical integration tolerance!'
         write(6, *) 'Tolerance must be .LT. ', Na
         STOP

      ENDIF

!
!     *** End of SUBROUTINE LIMIT_S ***
!
      END SUBROUTINE LIMIT_S

