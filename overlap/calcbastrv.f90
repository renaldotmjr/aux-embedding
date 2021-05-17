     SUBROUTINE CALCBASTRV(RNIJ)
!
!     Purpose: CALCulate the BASe TRanslation Vector
!
!     History: - Creation (31.05.12, MJR)
!
!     ******************************************************************
!     List of input/output variables:
!     RNIJ               : INPUT  -> RIJ coordinates
!                          OUTPUT -> Ni vector (base translation vector)
!
!     List of USED variables in EMBPARAM:
!     LATVEC(3,3)        : The lattice vectors in a filed (REAL)
!     LATVECNORM(3)      : The lattice vec. normas    (REAL)
!
!
!     ******************************************************************
!
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

!     List of local variables:
      INTEGER       :: I
      REAL (8)      :: RNIJ(3), PROJRLAT(3), BASEVEC(3)


      IF (DEBUG) THEN

          WRITE(6,*) 'XYZ_ATM = ', RNIJ(1), RNIJ(2), RNIJ(3)

      ENDIF

      ! Projecting the RNIJ at all lattice vectors
      DO I=1, 3

          PROJRLAT(I) = (DOT_PRODUCT( RNIJ, LATVEC(I,:) ) ) / LATVECNORM(I)

      ENDDO

      IF (DEBUG) THEN

          WRITE(6,*) 'PROJRLAT = ', PROJRLAT(1), PROJRLAT(2), PROJRLAT(3)

      ENDIF

      ! Analysing PROJRLAT(3) to construct the base translation vector
      DO I=1, 3

          IF ( DABS(PROJRLAT(I)) < (0.5D0 - DIFTHRESHOLD) ) THEN    ! -0.5 < Proj(i) < 0.5

              BASEVEC(I) = 0.0D0

          ELSE IF ( (PROJRLAT(I) < (-0.5D0 - DIFTHRESHOLD) ) .OR.  &  ! Proj(i) <= -0.5
              &    ( (PROJRLAT(I) + 0.5D0) < DIFTHRESHOLD )   ) THEN

              BASEVEC(I) = 1.0D0

          ELSE IF ( (PROJRLAT(I) > (0.5D0 + DIFTHRESHOLD) ) .OR.  &  ! Proj(i) >= 0.5
              &    ( (PROJRLAT(I) - 0.5D0) < DIFTHRESHOLD )   ) THEN

              BASEVEC(I) = -1.0D0

          ELSE  ! Should never enter here!

              WRITE(6,*)
              WRITE(6,*) 'Problems with the translation vectors... Aborting!'

              STOP

          ENDIF


      ENDDO

      RNIJ = BASEVEC

      IF (DEBUG) THEN

          WRITE(6,*) 'BASEVEC = ', BASEVEC(1), BASEVEC(2), BASEVEC(3)

      ENDIF


!
!     *** End of SUBROUTINE PRJRIJLAT ***
!
      END SUBROUTINE CALCBASTRV
