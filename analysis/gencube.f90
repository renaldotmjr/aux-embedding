     SUBROUTINE GENCUBE(FILNAM)
!
!     Purpose: Generate file.cube of densities
!
!     History: - Creation (27.03.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     NA                  : Normalization coefficient

!
!     List of local variables:

!
!     List of USED variables from EMBPARAM:
!     RHO(,,,)                : Mapped density RHO(X,Y,Z)

!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      CHARACTER (80)          :: FILNAM

      INTEGER                 :: I, IX, IY, IZ

      REAL (8)                :: ERROR, RHO_F, xe, ye, ze
      REAL (8), DIMENSION(3)  :: re

!     Starting the variables
      RHO_F = 0.0D0
      ERROR = 0.0D0



      OPEN(2, FILE=FILNAM, STATUS='UNKNOWN')

      WRITE(2, *) 'TITLE 1'
      WRITE(2, *) 'TITLE 2'

      WRITE(2, '(I5, 3F14.8)') NATOM, 0.0, 0.0, 0.0

      WRITE(2, '(I5,3F14.8)') NX, A_LATx/DFLOAT(NX-1),   0.0,   0.0
      WRITE(2, '(I5,3F14.8)') NY,   0.0, B_LATy/DFLOAT(NY-1),   0.0
      WRITE(2, '(I5,3F14.8)') NZ,   0.0,   0.0, C_LATz/DFLOAT(NZ-1)

      DO I=1, NATOM
          WRITE(2, '(I5,4F14.8)') ZATOM(I), 0.0, X0(I), Y0(I), Z0(I)
      ENDDO

      DO IX=1, NX
          DO IY=1, NY
              DO IZ=1, NZ

!!!!!!!!!!!!!!!!!!!
                  IF ( (IZ .EQ. 1) .AND. (IY .EQ. 1) ) THEN
!!!!!!!!!!!!!!!!!!!

                  ! Calculating the coords of electron in crys coords
                  ! This coordinates are always in the mapped grid.
                  xe = del_x * DFLOAT(IX - 1)
                  ye = del_y * DFLOAT(IY - 1)
                  ze = del_z * DFLOAT(IZ - 1)

                  !Switching to cartesian coordinates
                  CALL CRYS2CART(xe, ye, ze)

                  ! Calculating the fitted RHO in the point (xe, ye, ze)
                  re(1) = xe
                  re(2) = ye
                  re(3) = ze
                  CALL RHOFIT(re, RHO_F)

                  ERROR =  RHO_F - RHO(IX,IY,IZ)


!!!!!!!!!!!!!!!!!!
!                  IF ( (IZ .EQ. 1) .AND. (IY .EQ. 1) ) THEN

                      WRITE(6, '(4F20.6)') xe, RHO(IX,IY,IZ), RHO_F, ERROR

!                  ENDIF
!!!!!!!!!!!!!!!!!!


                  IF ( IZ .EQ. NZ) THEN

                      WRITE(2, '(F15.6)') ERROR

                  ELSE

                      IF ( MOD(IZ,6) .NE. 0 ) WRITE(2, '(F15.6,$)') ERROR

                      IF ( MOD(IZ,6) .EQ. 0 ) WRITE(2, '(F15.6)') ERROR

                  ENDIF


!!!!!!!!!!!!!!!!!!
                  ENDIF
!!!!!!!!!!!!!!!!!!


              ENDDO
          ENDDO
      ENDDO

       CLOSE(2)

!
!     *** End of SUBROUTINE GENCUBE ***
!
      END SUBROUTINE GENCUBE
