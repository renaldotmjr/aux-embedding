      SUBROUTINE EMBNORM
!
!     Purpose: Calculate the AUX primitive norms
!
!     History: - Creation (27.09.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     ZETA
!     FACTOR
!     ISET
!     ISHL
!     LC, LX, LY, LZ
!
!     List of USED variables in EMBPARAM:
!     NSET                    : Number of aux basis functions
!     LAUXSHL               : Lower pointer to aux function shells
!     UAUXSHL               : Upper pointer to aux function shells
!     AUXL(,3)              : Matrix with all primitives hermite
!                             gaussian function angular momentum index:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     EMBNORMA()            : Normalization coefficients of each shell
!
!     ******************************************************************
      USE EMBPARAM

      IMPLICIT NONE

      REAL (8)             :: ZETA, FACTOR, LC
      INTEGER              :: ISET, ISHL, LX, LY, LZ

!     Allocating tha EMBNORMA vector
      ALLOCATE( EMBNORMA(NSHL) )

!     Loop over each ZET()
      DO ISET=1, NSET

!        Getting the ZETA and LC values
         ZETA = ZET(ISET)
         LC = DFLOAT( LTOTAL(ISET) )
         FACTOR = DSQRT(2.0D0) * ((PI/ZETA)**(2.5D0)) * ( (ZETA**LC)/(2.0D0*LC+1.0D0) )

!        Loop over each AUX shell
         do ISHL=LAUXSHL(ISET), UAUXSHL(ISET)

            LX = AUXL(ISHL,1)
            LY = AUXL(ISHL,2)
            LZ = AUXL(ISHL,3)
            EMBNORMA(ISHL) = DFAC(2*LX-1) * DFAC(2*LY-1) * DFAC(2*LZ-1) * FACTOR
            EMBNORMA(ISHL) = 1.0D0/DSQRT(EMBNORMA(ISHL))

         enddo

      ENDDO

!
!     *** End of SUBROUTINE EMBNORM ***
!
      END SUBROUTINE EMBNORM
