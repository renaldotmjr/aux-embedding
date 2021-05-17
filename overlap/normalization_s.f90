      SUBROUTINE NORMALIZATION_S(NORMA, ZETA, NORMA_S)
!
!     Purpose: Calculation of normalization integral <a>
!
!     History: - Creation (03.10.11, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     NORMA    : Normalization coefficient of a
!     ZETA     : Gaussian exponent of auxiliary function A.
!
!     List of output variable:
!     NORMA_S = <a>
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE
!
!     ------------------------------------------------------------------
!
!     List of input variables:
      REAL(8) :: NORMA, ZETA
!
!     List of output variable:
      REAL(8) :: NORMA_S
!
!     ------------------------------------------------------------------
!
      NORMA_S =  NORMA*((PI/ZETA)**1.5)  ! <a> = Na*(PI/ZETA)^1.5
!
!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_SS ***
!
      END SUBROUTINE NORMALIZATION_S
