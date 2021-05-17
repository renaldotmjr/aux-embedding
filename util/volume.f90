      SUBROUTINE VOLUME
!
!     Purpose: Calculate the unit cell volume and the volume element
!
!     History: - Creation (27.09.11, MJR)
!
!     ******************************************************************
!
!     List of USED variables in EMBPARAM:
!     NX, NY, NZ              : Grid in a, b and c
!     A_LATx, A_LATy, A_LATz  : Lattice parameter a
!     B_LATx, B_LATy, B_LATz  : Lattice parameter b
!     C_LATx, C_LATy, C_LATz  : Lattice parameter c
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     VOLUME                  : Unit Cell volume (CONVERTED TO A.U.)
!     VOLUME_EL               : Volume element for numerical integration
!
!     List of local variables:
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      REAL(8)    ::     aXb_x, aXb_y, aXb_z

!
!     V = | (a x b) * c |   -> Where 'a', 'b' and 'c' are vectors
!
!     (a x b) = (aXb_x, aXb_y, aXb_z)
      aXb_x = A_LATy*B_LATz - A_LATz*B_LATy
      aXb_y = A_LATz*B_LATx - A_LATx*B_LATz
      aXb_z = A_LATx*B_LATy - A_LATy*B_LATx

!     V = | (a x b) * c |
      VOLUME_UC = DABS( aXb_x*C_LATx + aXb_y*C_LATy + aXb_z*C_LATz )

!     The volume element = V/grid
      VOLUME_EL = VOLUME_UC/DFLOAT(NX*NY*NZ)
!
!     *** End of SUBROUTINE VOLUME ***
!
      END SUBROUTINE VOLUME
