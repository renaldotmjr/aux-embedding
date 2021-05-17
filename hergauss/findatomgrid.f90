     SUBROUTINE FINDATOMGRID(XC,YC,ZC)
!
!     Purpose: Find the index in grid nearest the atom
!
!     History: - Creation (17.10.11, MJR)
!
!     ******************************************************************
!     List of INPUT/OUTPUT variables:
!     XC,YC,ZC          : INPUT  -> Atomic CARTESIAN coordinates
!                         OUTPUT -> Atomic CRYSTAL coordinates of the
!                                   nearest point in grid
!
!     List of local variables
!     Ix,Iy,Iz          : Index of the point in rho
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER       :: Ix,Iy,Iz

      REAL (8)      :: XC, YC, ZC


!     Changing coordinates from cartesian to crystalline
      CALL CART2CRYS(XC, YC, ZC)

!     Calculating the index in grid nearest the atom position
      CALL CRYS2INDEX(XC,YC,ZC, Ix,Iy,Iz)

!     Calculating the crystal coordinates of the nearest point in grid
      CALL INDEX2CRYS(Ix,Iy,Iz, XC,YC,ZC)
!
!     *** End of SUBROUTINE FINDATOMGRID ***
!
      END SUBROUTINE FINDATOMGRID

