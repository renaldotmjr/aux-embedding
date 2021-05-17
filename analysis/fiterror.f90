     SUBROUTINE FITERROR(p_ERROR, ME, SD)
!
!     Purpose: Calculate the error function between
!                the input and fitted rho
!
!     History: - Creation (13.02.12, MJR)
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

      INTEGER                 :: IX, IY, IZ

      REAL (8)                :: p_ERROR, ME, Erro, SD, RHO_F, xe, ye, ze
      REAL (8), DIMENSION(3)  :: re

!     Starting the variables
      RHO_F = 0.0D0
      ME = 0.0D0
      p_ERROR = 0.0D0

      DO IX=1, NX
          DO IY=1, NY
              DO IZ=1, NZ

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

                  p_ERROR = p_ERROR + ( RHO(IX,IY,IZ) - RHO_F )**2

                  ! Calculating the Mean Error
                  ME = ME + ( RHO_F - RHO(IX,IY,IZ) )

              ENDDO
          ENDDO
      ENDDO

      ! The error with unit: number of electrons
      p_ERROR = DSQRT(p_ERROR) * VOLUME_EL

      ! Calculating the Mean Error
      ME = (1.0D0/(NX*NY*NZ)) * ME

      !
      ! CALCULATING THE SAMPLE STANDARD DEVIATION
      !
      Erro = 0.0D0
      RHO_F = 0.0D0
      SD = 0.0D0

      DO IX=1, NX
          DO IY=1, NY
              DO IZ=1, NZ

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

                  Erro = RHO_F - RHO(IX,IY,IZ)

                  ! Calculating the Standard Deviation
                  SD = SD + (Erro - ME)**2

              ENDDO
          ENDDO
      ENDDO

      ! Calculating the Standard Deviation
      SD = DSQRT( (1.0D0/(NX*NY*NZ - 1)) * SD )

!
!     *** End of SUBROUTINE FITERROR ***
!
      END SUBROUTINE FITERROR
