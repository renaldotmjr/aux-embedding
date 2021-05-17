     SUBROUTINE INDEX2CRYS(Ix,Iy,Iz, Xa,Ya,Za)
!
!     Purpose: Calculate the index associated with the crystalline coords.
!
!     History: - Creation (17.10.11, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     Ix,Iy,Iz                : Index of the point in rho
!
!     List of USED variables from EMBPARAM:
!     del_x, del_y, del_z     : Steps of cryst. coords.
!
!     List of OUTPUT variables
!     X,Y,Z                   : point crystalline coordinates
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER       :: Ix,Iy,Iz

      REAL (8)      :: Xa, Ya, Za

!     Calculating the crystalline coordinates
      Xa = (DFLOAT(Ix-1))*del_x
      Ya = (DFLOAT(Iy-1))*del_y
      Za = (DFLOAT(Iz-1))*del_z

!
!     *** End of SUBROUTINE INDEX2CRYS ***
!
      END SUBROUTINE INDEX2CRYS
