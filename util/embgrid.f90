      SUBROUTINE EMB_GRID
!
!     Purpose: Calculate the crystaline coordinates of the mapped grid
!
!     History: - Creation (20.09.11, MJR)
!              - Modified (17.10.11, MJR) - Grid in crystal coords
!
!     ******************************************************************
!
!     List of input variables:
!     FILE_NAME    : PWScf output file name
!MODULE
!     List of local variables:
!
!     List of USED variables from EMBPARAM:
!     NX, NY, NZ                     : Grid in a, b and c
!     A_LATx, A_LATy, A_LATz         : Lattice parameter a
!     B_LATx, B_LATy, B_LATz         : Lattice parameter b
!     C_LATx, C_LATy, C_LATz         : Lattice parameter c
!     
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     X(), Y(), Z()                  : Grid mapped coordinate vectors
!     del_x, del_y, del_z            : Steps of cryst. coords.
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER           :: IX, IY, IZ

!     Allocating the vectors with grid coords
      ALLOCATE ( X(NX), Y(NY), Z(NZ) )

!     Steps for the crystalline coordinates of mapped grid
      del_x = 1.0D0/DFLOAT(NX-1)
      del_y = 1.0D0/DFLOAT(NY-1)
      del_z = 1.0D0/DFLOAT(NZ-1)

!
!     Calculating the X, Y and Z vectors (Maybe this will not be used)
      DO IX=1,NX
         X(IX) = DFLOAT(IX-1)*del_x
      ENDDO
      DO IY=1,NY
         Y(IY) = DFLOAT(IY-1)*del_y
      ENDDO
      DO IZ=1,NZ
         Z(IZ) = DFLOAT(IZ-1)*del_z
      ENDDO

!
!     *** End of SUBROUTINE EMB_GRID ***
!
      END SUBROUTINE EMB_GRID
