     SUBROUTINE LATTICEVEC
!
!     Purpose: Organize ALATs in a field and calculate the norms
!
!     History: - Creation (31.05.12, MJR)
!
!     ******************************************************************
!     List of USED variables in EMBPARAM:
!     A_LATx, A_LATy, A_LATz      : Lattice parameter a        (REAL)
!     B_LATx, B_LATy, B_LATz      : Lattice parameter b        (REAL)
!     C_LATx, C_LATy, C_LATz      : Lattice parameter c        (REAL)
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     LATVEC(3,3)                 : The lattice vectors in a filed (REAL)
!     LATVECNORM(3)               : The lattice vec. normas    (REAL)
!
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER        :: I, J, K, CC

!     Constructing the LATVEC field
      LATVEC(1,1) = A_LATx
      LATVEC(1,2) = A_LATy
      LATVEC(1,3) = A_LATz

      LATVEC(2,1) = B_LATx
      LATVEC(2,2) = B_LATy
      LATVEC(2,3) = B_LATz

      LATVEC(3,1) = C_LATx
      LATVEC(3,2) = C_LATy
      LATVEC(3,3) = C_LATz

!     Calculating the norm (squared)
      LATVECNORM(1) = DOT_PRODUCT( LATVEC(1,:), LATVEC(1,:))
      LATVECNORM(2) = DOT_PRODUCT( LATVEC(2,:), LATVEC(2,:))
      LATVECNORM(3) = DOT_PRODUCT( LATVEC(3,:), LATVEC(3,:))

      MAXVECNORM = MAXVAL(LATVECNORM)
      MINVECNORM = MINVAL(LATVECNORM)

      NTRANS3 = ( (2*NTRANS) +1)**3

      ALLOCATE( TRANSLATVEC(NTRANS3, 3) )

      CC = 1

      DO I=-NTRANS, NTRANS
          DO J=-NTRANS, NTRANS
              DO K=-NTRANS, NTRANS

                  TRANSLATVEC(CC,1) = DFLOAT(I)
                  TRANSLATVEC(CC,2) = DFLOAT(J)
                  TRANSLATVEC(CC,3) = DFLOAT(K)

                  CC = CC + 1

              ENDDO
          ENDDO
      ENDDO


!
!     *** End of SUBROUTINE LATTICEVEC ***
!
      END SUBROUTINE LATTICEVEC
