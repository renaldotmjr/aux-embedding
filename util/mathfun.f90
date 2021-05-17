     SUBROUTINE MATHFUN
!
!     Purpose: Calculation of mathematical constants and functions
!
!     History: - Creation (27.09.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     I
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     PI                      : Pi number
!     DEGRAD                  : DEG/RAD conversion factors
!     RADDEG                  : DEG/RAD conversion factors
!     FAC(A)                  : Vector with factorial of A
!     DFAC(A)                 : Vector with double factorial of A
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER   :: I

!     *** Definition of Pi and DEG/RAD conversion factors ***

      PI = 2.0*ACOS(0.0)

      DEGRAD = PI/180.0
      RADDEG = 180.0/PI
!
!     *** Definition of factorial function ***
!
      FAC(0) = 1.0

      DO I=1, MAXFAC
         FAC(I) = DFLOAT(I)*FAC(I-1)
      ENDDO
!
!     *** Definition of double factorial function ***
!
      DFAC(-1) = 1.0
      DFAC(0) = 1.0
      DFAC(1) = 1.0
      DFAC(2) = 2.0

      DO I=3, 2*MAXFAC+1
         DFAC(I) = DFLOAT(I)*DFAC(I-2)
      ENDDO

!
!     *** End of SUBROUTINE MOUNT_ZET ***
!
      END SUBROUTINE MATHFUN
