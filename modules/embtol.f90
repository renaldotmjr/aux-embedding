      MODULE EMBTOL
!     TOLERANCE DEFINITIONS FOR EMBEDDING CODE
!
!     Purpose: EMBedding TOLerance definitions
!
!     History: - Creation (19.10.11, MJR)
!
!     ******************************************************************
!
!     List of Variables
!     INTTOL           : Tolerance for numerical integrations
!     XTOL             : Tolerance for iterativelly solving the
!                          numerical integrations limits
!     XTMAX            : Maximum of iterations in limits calculation
!
!     ******************************************************************
!
      REAL(8) :: INTTOL, XTOL

      INTEGER :: XTMAX

!
!     *** End of MODULE EMBTOL ***
!
      END MODULE EMBTOL
