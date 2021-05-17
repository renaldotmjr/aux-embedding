     SUBROUTINE DEALLOC_EMB
!
!     Purpose: Deallocate everything
!
!     History: - Creation (28.09.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     Deallocating everything
      DEALLOCATE( SYMB, ZATOM, X0, Y0, Z0, RHO, X , Y , Z, ZET, LTOTAL, &
                & ORBTOTAL, LAUXPTR, UAUXPTR, LAUXSHL, UAUXSHL, AUXL, &
                & EMBNORMA, VEC_b )


!
!     *** End of SUBROUTINE DEALLOC_EMB ***
!
      END SUBROUTINE DEALLOC_EMB
