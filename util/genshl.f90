      SUBROUTINE GENSHL
!
!     Purpose: Mount the vectors: ORBTOTAL, LAUXSHL, UAUXSHL, AUXL
!
!     History: - Creation (26.09.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     ISET                  :
!     ISHL                  :
!     AL                    :
!     COUNTL                :
!
!     List of USED variables in EMBPARAM:
!     NSET                  : Number of aux basis functions
!     LTOTAL()              : Aux. basis total angular momentum of
!                             each atomic orbital
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     NSHL                  : Total number of primitives hermite
!                              gaussian functions
!     ORBTOTAL()            : Number of primitive hermite gaussian
!                              functions of each aux func
!     LAUXSHL               : Lower pointer to aux function shells
!     UAUXSHL               : Upper pointer to aux function shells
!     AUXL(,3)              : Matrix with all primitives hermite
!                              gaussian function angular momentum index:
!                               (1) -> Lax
!                               (2) -> Lay
!                               (3) -> Laz
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      INTEGER              :: ISET, ISHL, AL, COUNTL

!     Initializing NSHL
      NSHL = 0

!     Mounting ORBTOTAL
      ALLOCATE( ORBTOTAL(NSET) )
!     Calculating the number of shells of each aux function
      DO ISET=1, NSET

         AL = LTOTAL(ISET) + 1
         ORBTOTAL(ISET) = (AL*(AL+1))/2.0D0

         ! Accounting to NSHL
         NSHL = NSHL + ORBTOTAL(ISET)
      ENDDO

!     Generating LAUXSHL and UAUXSHL
      ALLOCATE( LAUXSHL(NSET), UAUXSHL(NSET), AUXL(NSHL, 3) )

      COUNTL = 1

      DO ISET=1, NSET

!        If the orbital is an S-TYPE
         IF (LTOTAL(ISET) == 0) THEN

            ! Passing lower and upper index
            LAUXSHL(ISET) = COUNTL
            UAUXSHL(ISET) = COUNTL + (ORBTOTAL(ISET) - 1)

            ! Passing Lax, Lay and Laz of each AO aux function
            AUXL(COUNTL  ,1) = 0 ; AUXL(COUNTL  ,2) = 0 ; AUXL(COUNTL  ,3) = 0

            COUNTL = COUNTL + ORBTOTAL(ISET)

         ENDIF

!        If the orbital is an P-TYPE
         IF (LTOTAL(ISET) == 1) THEN

            ! Passing lower and upper index
            LAUXSHL(ISET) = COUNTL
            UAUXSHL(ISET) = COUNTL + (ORBTOTAL(ISET) - 1)

            ! Passing Lax, Lay and Laz of each AO aux function
            AUXL(COUNTL  ,1) = 1 ; AUXL(COUNTL  ,2) = 0 ; AUXL(COUNTL  ,3) = 0
            AUXL(COUNTL+1,1) = 0 ; AUXL(COUNTL+1,2) = 1 ; AUXL(COUNTL+1,3) = 0
            AUXL(COUNTL+2,1) = 0 ; AUXL(COUNTL+2,2) = 0 ; AUXL(COUNTL+2,3) = 1

            COUNTL = COUNTL + ORBTOTAL(ISET)

         ENDIF

!        If the orbital is an D-TYPE
         IF (LTOTAL(ISET) == 2) THEN

            ! Passing lower and upper index
            LAUXSHL(ISET) = COUNTL
            UAUXSHL(ISET) = COUNTL + (ORBTOTAL(ISET) - 1)

            ! Passing Lax, Lay and Laz of each AO aux function
            AUXL(COUNTL  ,1) = 2 ; AUXL(COUNTL  ,2) = 0 ; AUXL(COUNTL  ,3) = 0
            AUXL(COUNTL+1,1) = 1 ; AUXL(COUNTL+1,2) = 1 ; AUXL(COUNTL+1,3) = 0
            AUXL(COUNTL+2,1) = 1 ; AUXL(COUNTL+2,2) = 0 ; AUXL(COUNTL+2,3) = 1
            AUXL(COUNTL+3,1) = 0 ; AUXL(COUNTL+3,2) = 2 ; AUXL(COUNTL+3,3) = 0
            AUXL(COUNTL+4,1) = 0 ; AUXL(COUNTL+4,2) = 1 ; AUXL(COUNTL+4,3) = 1
            AUXL(COUNTL+5,1) = 0 ; AUXL(COUNTL+5,2) = 0 ; AUXL(COUNTL+5,3) = 2

            COUNTL = COUNTL + ORBTOTAL(ISET)

         ENDIF

      ENDDO

!
!     *** End of SUBROUTINE GETBASIS ***
!
      END SUBROUTINE GENSHL
