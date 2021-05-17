     SUBROUTINE MOUNTTVECS
!
!     Purpose: Mount a set of translational vectors
!
!     History: - Creation (13.06.12, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!
!     List of local variables:
!
!     ******************************************************************

      USE EMBPARAM
      USE INPDEFVARS

      IMPLICIT NONE

      ! List of local variables
      INTEGER        :: I, J, K, L, CC

      LOGICAL        :: LETUC

      REAL (8)       :: LOC_LATVEC(3,3), base_in(3), MOD_base

!     Constructing the LOC_LATVEC field
      LOC_LATVEC(1,1) = A_LATx
      LOC_LATVEC(1,2) = A_LATy
      LOC_LATVEC(1,3) = A_LATz

      LOC_LATVEC(2,1) = B_LATx
      LOC_LATVEC(2,2) = B_LATy
      LOC_LATVEC(2,3) = B_LATz

      LOC_LATVEC(3,1) = C_LATx
      LOC_LATVEC(3,2) = C_LATy
      LOC_LATVEC(3,3) = C_LATz


      ! Calculating the number of translation vectors, not including the (0,0,0)
      SIZE_TV = (REPNUM(1,2) - REPNUM(1,1) + 1)* &
             &  (REPNUM(2,2) - REPNUM(2,1) + 1)* &
             &  (REPNUM(3,2) - REPNUM(3,1) + 1)  &
             & - NTAKEUCOUT

      ! Allocating the field TRANSVECS
      ALLOCATE( TRANSVECS(SIZE_TV, 3) )

      WRITE(6,*)
      WRITE(6,*) 'Take out the vectors in deMon.cub file'
      do i = 1, NTAKEUCOUT
          WRITE(6,"(3I5)") TAKEUCOUTvecs(i,1), TAKEUCOUTvecs(i,2), TAKEUCOUTvecs(i,3)
      enddo
      WRITE(6,*)


      ! Calculating the tranlation vectors
      CC = 1
      DO I=REPNUM(1,1), REPNUM(1,2)
          DO J=REPNUM(2,1), REPNUM(2,2)
              DO K=REPNUM(3,1), REPNUM(3,2)

                  base_in(1) = DFLOAT(I)
                  base_in(2) = DFLOAT(J)
                  base_in(3) = DFLOAT(K)
!                  MOD_base = DOT_PRODUCT( base_in ,base_in )

                  LETUC = .TRUE.
                  DO L=1, NTAKEUCOUT

                      IF ( ((I .EQ. TAKEUCOUTvecs(L,1)) .AND. &
                          & (J .EQ. TAKEUCOUTvecs(L,2)) .AND. &
                          & (K .EQ. TAKEUCOUTvecs(L,3))) ) THEN

                          LETUC = .FALSE.
                          EXIT

                      ENDIF

                  ENDDO

                  ! IF to take out the original unit cell vector (0,0,0)
                  !IF ( DABS(MOD_base) > 1.0D-6 ) THEN
                  IF ( LETUC ) THEN

                      DO L=1,3
                          TRANSVECS(CC,L) = base_in(1)*LOC_LATVEC(1,L) + &
                                        &   base_in(2)*LOC_LATVEC(2,L) + &
                                        &   base_in(3)*LOC_LATVEC(3,L)
                      ENDDO

                      CC = CC + 1

                  ENDIF

              ENDDO
          ENDDO
      ENDDO

!
!     *** End of SUBROUTINE MOUNTTVECS ***
!
      END SUBROUTINE MOUNTTVECS
