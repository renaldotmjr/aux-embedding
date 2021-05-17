      SUBROUTINE UNITCONV
!
!     Purpose: Convert the units to atomic units
!
!     History: - Creation (17.10.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!
!     List of USED variables from EMBPARAM:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!
!     ******************************************************************
!
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

!     Conveting the lattice parameters
      A_LATx = A_LATx*ANG2AU
      B_LATx = B_LATx*ANG2AU
      C_LATx = C_LATx*ANG2AU
      A_LATy = A_LATy*ANG2AU
      B_LATy = B_LATy*ANG2AU
      C_LATy = C_LATy*ANG2AU
      A_LATz = A_LATz*ANG2AU
      B_LATz = B_LATz*ANG2AU
      C_LATz = C_LATz*ANG2AU

!     Converting the atomic coordinates
      X0 = X0*ANG2AU
      Y0 = Y0*ANG2AU
      Z0 = Z0*ANG2AU



!
!     *** End of SUBROUTINE UNITCONV ***
!
      END SUBROUTINE UNITCONV

