      SUBROUTINE CRYS2CART(Xa, Ya, Za)
!
!     Purpose: Convert coords from crystalline to cartesian
!
!     History: - Creation (07.10.11, MJR)
!
!     ******************************************************************
!
!     List of input/output variables:
!     X, Y, Z           : input -> crystalline coordinates
!                         output -> cartesian coordinates
!
!     List of USED variables from EMBPARAM:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!
!     ******************************************************************

      USE EMBPARAM

      IMPLICIT NONE

      REAL(8)  :: CART(3,1), CRYS(3,1)
      REAL(8)  :: Xa, Ya, Za

      CRYS(1,1) = Xa
      CRYS(2,1) = Ya
      CRYS(3,1) = Za

      CART = MATMUL(MLAT,CRYS)

      Xa = CART(1,1)
      Ya = CART(2,1)
      Za = CART(3,1)

!
!     *** End of SUBROUTINE CART2CRYS ***
!
      END SUBROUTINE CRYS2CART
