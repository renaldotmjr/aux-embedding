     SUBROUTINE CRYS2INDEX(Xa,Ya,Za, Ix,Iy,Iz)
!
!     Purpose: Calculate the index associated with the crystalline coords.
!
!     History: - Creation (17.10.11, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     X,Y,Z                   : point crystalline coordinates
!
!     List of USED variables from EMBPARAM:
!     del_x, del_y, del_z     : Steps of cryst. coords.
!
!     List of OUTPUT variables
!     Ix,Iy,Iz                : Index of the point in rho
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER       :: Ix,Iy,Iz

      REAL (8)      :: Xa, Ya, Za

!     Calculating the index in grid
      Ix = IDNINT( (Xa/del_x) + 1.0 )
      Iy = IDNINT( (Ya/del_y) + 1.0 )
      Iz = IDNINT( (Za/del_z) + 1.0 )

!
!     *** End of SUBROUTINE CRYS2INDEX ***
!
      END SUBROUTINE CRYS2INDEX
