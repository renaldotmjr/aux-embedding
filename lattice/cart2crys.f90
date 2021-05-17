      SUBROUTINE CART2CRYS(Xa, Ya, Za)
!
!     Purpose: Convert coords from cartesian to crystalline
!
!     History: - Creation (07.10.11, MJR)
!
!     ******************************************************************
!
!     List of input/output variables:
!     X, Y, Z           : input -> cartesian coordinates
!                         output -> crystalline coordinates
!
!     List of USED variables from EMBPARAM:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!
!     ******************************************************************

      USE EMBPARAM

      IMPLICIT NONE

      REAL(8)  :: CART(3,1)
      REAL(8)  :: CRYS(3,1)
      REAL(8)  :: Xa, Ya, Za

      CART(1,1) = Xa
      CART(2,1) = Ya
      CART(3,1) = Za

      CRYS = MATMUL(MLAT_inv, CART)

      Xa = CRYS(1,1)
      Ya = CRYS(2,1)
      Za = CRYS(3,1)

!
!     *** End of SUBROUTINE CART2CRYS ***
!
      END SUBROUTINE CART2CRYS
