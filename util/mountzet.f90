      SUBROUTINE MOUNT_ZET(ZETA, LTOT)
!
!     Purpose: Mount the vector with auxiliary basis set zeta values
!                -> This subroutine reallocates some vectors
!
!     History: - Creation (23.09.11, MJR)
!              - Inclusion of deMon standard
!                         (06.08.12, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     ZETA        : Auxiliary basis zeta value
!     LTOT        : Total angular momentum of the aux basis
!
!     List of local variables:
!     LOC_ZET     : The local zeta vector
!     LOC_LTOT    : The local total L vector
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     NSET                  : Number of aux basis functions
!     ZET()                 : Aux. basis zeta of each atomic orbital
!     LTOTAL()              : Aux. basis total angular momentum of
!                               each atomic orbital
!
!     ******************************************************************
!
      USE EMBPARAM
      USE INPDEFVARS

      IMPLICIT NONE

!     Local variables
      REAL (8), DIMENSION(:),      ALLOCATABLE :: LOC_ZET
      INTEGER , DIMENSION(:),      ALLOCATABLE :: LOC_LTOT

      REAL (8)                                 :: ZETA
      INTEGER           :: I, J, LTOT, INDEX_NSET

!     Enter in the first read of the basis in file!
      IF (NSET == 0) THEN

         IF (AUXIS_STD .EQ. 'OPEN') THEN

            NSET = NSET + 1
!           Here initialize the matrixes ZET and LTOTAL
            ALLOCATE(  ZET(NSET), LTOTAL(NSET) )
            ZET(NSET) = ZETA
            LTOTAL(NSET) = LTOT

         ELSEIF (AUXIS_STD .EQ. 'deMon') THEN

            NSET = NSET + (1+LTOT)
!           Here initialize the matrixes ZET and LTOTAL
            ALLOCATE(  ZET(NSET), LTOTAL(NSET) )

            DO J=0, LTOT
               ZET(J+1) = ZETA
               LTOTAL(J+1) = J
            ENDDO

         ENDIF

      ELSE

         IF (AUXIS_STD .EQ. 'OPEN') THEN

            ! Allocating the local temporary matrix to preserv the old data
            ALLOCATE( LOC_ZET(NSET), LOC_LTOT(NSET) )

            ! Copying the old data to local variables
            do I=1, NSET
               LOC_ZET(I) = ZET(I)
               LOC_LTOT(I) = LTOTAL(I)
            enddo

            ! deallocating the old vector
            DEALLOCATE( ZET, LTOTAL)

            NSET = NSET + 1

            ! Allocating the new vector
            ALLOCATE(  ZET(NSET), LTOTAL(NSET) )
            ! Passing the old data to the new vectors
            do I=1, (NSET-1)
               ZET(I) = LOC_ZET(I)
               LTOTAL(I) = LOC_LTOT(I)
            enddo
            ! Passing the new values to the end of vectors
            ZET(NSET) = ZETA
            LTOTAL(NSET) = LTOT

            ! Deallocating the local variables
            DEALLOCATE( LOC_ZET, LOC_LTOT )

         ELSEIF (AUXIS_STD .EQ. 'deMon') THEN


            ! Allocating the local temporary matrix to preserv the old data
            ALLOCATE( LOC_ZET(NSET), LOC_LTOT(NSET) )
            ! Copying the old data to local variables
            do I=1, NSET
               LOC_ZET(I) = ZET(I)
               LOC_LTOT(I) = LTOTAL(I)
            enddo

            ! deallocating the old vector
            DEALLOCATE( ZET, LTOTAL)

            NSET = NSET + (1+LTOT)

            ! Allocating the new vector
            ALLOCATE(  ZET(NSET), LTOTAL(NSET) )
            ! Passing the old data to the new vectors
            do I=1, ( NSET - (1+LTOT) )
               ZET(I) = LOC_ZET(I)
               LTOTAL(I) = LOC_LTOT(I)
            enddo

            ! Passing the new values to the end of vectors
            DO J=0, LTOT
               INDEX_NSET = NSET - (1+LTOT) + (J+1)
               ZET(INDEX_NSET) = ZETA
               LTOTAL(INDEX_NSET) = J
            ENDDO

            ! Deallocating the local variables
            DEALLOCATE( LOC_ZET, LOC_LTOT )

         ENDIF


      ENDIF

!
!     *** End of SUBROUTINE MOUNT_ZET ***
!
      END SUBROUTINE MOUNT_ZET
