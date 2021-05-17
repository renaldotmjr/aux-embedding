     SUBROUTINE HERGAUSS_P(NA, ZETA, XA, YA, ZA, RESULT, IGRID)
!
!     Purpose: Calculation of numerical integrals <p_A|rho>
!              THE NUCLEAR COORDINATES MUST COME IN ATOMIC UNITS (Bohr)
!
!     History: - Creation (28.09.11, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     NA                  : Normalization coefficients ***
!                            -> NA(1) = NA_px
!                            -> NA(2) = NA_py
!                            -> NA(3) = NA_pz
!                         *** In this case, all NA are equal, but the
!                             implementation will be done separatelly
!                             to preserve the generality in the CALL
!
!     ZETA                : Basis exponent
!     XA, YA, ZA          : Nuclear cartesian coordinates
!
!     List of local variables:
!     r_eA2, DIFX, DIFY, DIFZ
!
!     List of USED variables from EMBPARAM:
!     RHO(,,,)                : Mapped density RHO(X,Y,Z)
!     X() , Y() , Z()         : Mapped grid coordinates
!     VOLUME_EL               : Volume element
!
!     List of OUTPUT variables
!     HG_P(1) = <px_A|rho>
!     HG_P(2) = <py_A|rho>
!     HG_P(3) = <pz_A|rho>
!
!     ******************************************************************
!
      USE EMBPARAM
      USE EMBTOL

      IMPLICIT NONE

!     List of local variables:
      INTEGER       :: IGRID, IX,IY,IZ, NI(3)
      INTEGER       :: N_int                    ! Number of points in the integration algorithm
      INTEGER       :: IRx, IRy, IRz            ! Index to get RHO

      REAL (8)      :: NA(3), ZETA, XA,YA,ZA, HG_P(3)
      REAL (8)      :: XC,YC,ZC, r_eA2, rhoI, LIMIT
      REAL (8)      :: RESULT(3)
      REAL (8)      :: xe,ye,ze, xec,yec,zec    ! The electron position

      REAL (8)      :: HERGA_P

      REAL (8), DIMENSION(:), ALLOCATABLE  :: INT_X1, INT_Y1, INT_Z1
      REAL (8), DIMENSION(:), ALLOCATABLE  :: INT_X2, INT_Y2, INT_Z2
      REAL (8), DIMENSION(:), ALLOCATABLE  :: INT_X3, INT_Y3, INT_Z3

      REAL (8)      :: INT_TEMP

!     Initializing HG_P and IGRID values
      RESULT = 0.0D0
      IGRID = 0

!     Calculating the limit distance (LIMIT) for a s-type hermite gaussian:
      CALL LIMIT_P(NA, ZETA, INTTOL, LIMIT)

!     Calculating the number of steps in each direction
      CALL CALCNSTEP(LIMIT, del_x, del_y, del_z, NI)


      ! *******************************************************************
      ! NEED TO DEFINE THE SIZE OF THE NUMERICAL INTEGRATION ALGORITHM HERE
      ! *******************************************************************
      !N_int = 3      ! For Simpson
      !N_int = 5      ! For Simpson_5
      N_int = 2      ! For Trapezoid or simple sum

      ! Making NI() a multiple of N_int
      NI(1) = NI(1) + ( N_int - MOD(NI(1),N_int) ) + N_int
      NI(2) = NI(2) + ( N_int - MOD(NI(2),N_int) ) + N_int
      NI(3) = NI(3) + ( N_int - MOD(NI(3),N_int) ) + N_int

      ! Allocating the vectors with the line integration results
      ! Allocate twice the NI lenght. If NI is a multiple of N_int, 2*NI is also
      ALLOCATE( INT_X1( 2*NI(1) ), INT_X2( 2*NI(1) ), INT_X3( 2*NI(1) ) )
      ALLOCATE( INT_Y1( 2*NI(2) ), INT_Y2( 2*NI(2) ), INT_Y3( 2*NI(2) ) )
      ALLOCATE( INT_Z1( 2*NI(3) ), INT_Z2( 2*NI(3) ), INT_Z3( 2*NI(3) ) )

      ! *******************************************************************


!     Calculating the crystal coordinates of the point in grid nearest the atom.
      ! Passing the atomic cartesian coordinates to (XC,YC,ZC)
      XC = XA
      YC = YA
      ZC = ZA
      CALL FINDATOMGRID(XC,YC,ZC)   ! Now (XC,YC,ZC) is a point in grid nearest
                                    ! the atom and is in crystalline coordinates

!     Integrating over the grid!
      DO IZ=1, 2*NI(3)

         DO IY=1, 2*NI(2)

            DO IX=1, 2*NI(1)

              ! Calculating the coords of electron in crys coords
              ! This coordinates are always in the mapped grid.
              xe = XC + del_x*DFLOAT(IX - NI(1))
              ye = YC + del_y*DFLOAT(IY - NI(2))
              ze = ZC + del_z*DFLOAT(IZ - NI(3))

              ! Converting the electron coordinates from crystal to cartesian
              xec = xe
              yec = ye
              zec = ze
              CALL CRYS2CART(xec, yec, zec)

              ! Centering in atom
              xec = xec - XA
              yec = yec - YA
              zec = zec - ZA

              ! Calculating the hermite gaussian values
              r_eA2 = (xec)**2 + (yec)**2 + (zec)**2

!             <px_A| = -2*Na_px*ZETA*(xe-XA)*exp( -ZETA*(r-A)^2 )
              HG_P(1) = HERGA_P(NA(1), ZETA, xec, r_eA2)

!             <py_A| = -2*Na_py*ZETA*(ye-YA)*exp( -ZETA*(r-A)^2 )
              HG_P(2) = HERGA_P(NA(2), ZETA, yec, r_eA2)

!             <pz_A| = -2*Na_pz*ZETA*(ze-ZA)*exp( -ZETA*(r-A)^2 )
              HG_P(3) = HERGA_P(NA(3), ZETA, zec, r_eA2)


              ! Putting the elec. coords. in unitary cube
              ! ***** This subroutine has IF structures inside!!
              CALL PUTINGRID(xe, ye, ze)

              ! Calculating the index of elec coords. in cube to get the RHO value
              CALL CRYS2INDEX(xe,ye,ze, IRx,IRy,IRz)

              ! Getting the RHO value
              rhoI = RHO(IRx, IRy, IRz)

              ! Passing the current calculation to the total result
              INT_X1(IX) = rhoI*HG_P(1)
              INT_X2(IX) = rhoI*HG_P(2)
              INT_X3(IX) = rhoI*HG_P(3)

              ! Accounting the number of inspections of the integration
              IGRID = IGRID +1
            ENDDO

            ! INTEGRATING X_Px
            CALL INTNUMCALL(INT_X1, 2*NI(1), INT_TEMP)
            INT_Y1(IY) = INT_TEMP
            ! INTEGRATING X_Py
            CALL INTNUMCALL(INT_X2, 2*NI(1), INT_TEMP)
            INT_Y2(IY) = INT_TEMP
            ! INTEGRATING X_Pz
            CALL INTNUMCALL(INT_X3, 2*NI(1), INT_TEMP)
            INT_Y3(IY) = INT_TEMP


         ENDDO

         ! INTEGRATING Y_Px
         CALL INTNUMCALL(INT_Y1, 2*NI(2), INT_TEMP)
         INT_Z1(IZ) = INT_TEMP
         ! INTEGRATING Y_Py
         CALL INTNUMCALL(INT_Y2, 2*NI(2), INT_TEMP)
         INT_Z2(IZ) = INT_TEMP
         ! INTEGRATING Y_Pz
         CALL INTNUMCALL(INT_Y3, 2*NI(2), INT_TEMP)
         INT_Z3(IZ) = INT_TEMP

      ENDDO

      ! INTEGRATING Z_Px
      CALL INTNUMCALL(INT_Z1, 2*NI(3), INT_TEMP)
      RESULT(1) = INT_TEMP
      ! INTEGRATING Z_Px
      CALL INTNUMCALL(INT_Z2, 2*NI(3), INT_TEMP)
      RESULT(2) = INT_TEMP
      ! INTEGRATING Z_Px
      CALL INTNUMCALL(INT_Z3, 2*NI(3), INT_TEMP)
      RESULT(3) = INT_TEMP


!     Multiplying by the volume element
      RESULT(1) = VOLUME_EL*RESULT(1)
      RESULT(2) = VOLUME_EL*RESULT(2)
      RESULT(3) = VOLUME_EL*RESULT(3)

      ! Deallocating ...
      DEALLOCATE(INT_X1, INT_Y1, INT_Z1)
      DEALLOCATE(INT_X2, INT_Y2, INT_Z2)
      DEALLOCATE(INT_X3, INT_Y3, INT_Z3)

!
!     *** End of SUBROUTINE HERGAUSS_P ***
!
      END SUBROUTINE HERGAUSS_P



