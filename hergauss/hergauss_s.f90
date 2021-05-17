     SUBROUTINE HERGAUSS_S(NA, ZETA, XA, YA, ZA, RESULT, IGRID)
!
!     Purpose: Calculation of numerical integral <s_A|rho>
!              THE NUCLEAR COORDINATES MUST COME IN ATOMIC UNITS (Bohr)
!
!     History: - Creation (27.09.11, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     NA                  : Normalization coefficient
!     ZETA                : Basis exponent
!     XA, YA, ZA          : Nuclear cartesian coordinates
!
!     List of local variables:
!     r_eA2
!
!     List of USED variables from EMBPARAM:VOLUME_EL
!     RHO(,,,)                : Mapped density RHO(X,Y,Z)
!     X() , Y() , Z()         : Mapped grid coordinates
!     VOLUME_EL               : Volume element
!
!     List of OUTPUT variables
!     RESULT = <s_A|rho>
!     IGRID                   : Number of inspectioned points.
!
!     ******************************************************************
!
      USE EMBPARAM
      USE EMBTOL

      IMPLICIT NONE

      INTEGER       :: IGRID, IX,IY,IZ, NI(3)
      INTEGER       :: N_int                    ! Number of points in the integration algorithm
      INTEGER       :: IRx, IRy, IRz            ! Index to get RHO

      REAL (8)      :: NA, ZETA, XA,YA,ZA
      REAL (8)      :: XC,YC,ZC, r_eA2, rhoI, LIMIT
      REAL (8)      :: HERGA_S, HG_S, RESULT(1)
      REAL (8)      :: xe,ye,ze, xec,yec,zec    ! The electron position

      REAL (8)      :: XCR,YCR,ZCR

      REAL (8), DIMENSION(:), ALLOCATABLE  :: INT_X, INT_Y, INT_Z
      REAL (8)      :: INT_TEMP


!     Initializing HG_S and IGRID values
      LIMIT = 0.0D0
      RESULT = 0.0D0
      IGRID = 0
      INT_TEMP = 0.0D0

!     Calculating the limit distance (LIMIT) for a s-type hermite gaussian:
      CALL LIMIT_S(NA, ZETA, INTTOL, LIMIT)

!     Calculating the number of steps in each direction
      CALL CALCNSTEP(LIMIT, del_x, del_y, del_z, NI)


      ! *******************************************************************
      ! NEED TO DEFINE THE SIZE OF THE NUMERICAL INTEGRATION ALGORITHM HERE
      ! *******************************************************************
      !N_int = 3      ! For Simpson
      !N_int = 5      ! For Simpson_5
      N_int = 2      ! For Trapezoid

      ! Making NI() a multiple of N_int
      NI(1) = NI(1) + ( N_int - MOD(NI(1),N_int) ) + N_int
      NI(2) = NI(2) + ( N_int - MOD(NI(2),N_int) ) + N_int
      NI(3) = NI(3) + ( N_int - MOD(NI(3),N_int) ) + N_int

      ! Allocating the vectors with the line integration results
      ! Allocate twice the NI lenght. If NI is a multiple of N_int, 2*NI is also
      ALLOCATE( INT_X( 2*NI(1) ) )
      ALLOCATE( INT_Y( 2*NI(2) ) )
      ALLOCATE( INT_Z( 2*NI(3) ) )

      ! *******************************************************************


!     Calculating the crystal coordinates of the point in grid nearest the atom.
      ! Passing the atomic cartesian coordinates to (XC,YC,ZC)
      XC = XA
      YC = YA
      ZC = ZA
      CALL FINDATOMGRID(XC,YC,ZC)   ! Now (XC,YC,ZC) is a point in grid nearest
                                    ! the atom and is in crystalline coordinates

!      WRITE(6,'(A,3F10.6)') 'Atomic Coords:  ', XA, YA, ZA
!      WRITE(6,*) '-----------------------------'
!      WRITE(6,'(A,F10.6,A,F10.6,A,F10.6)') 'S_limit:  ', LIMIT, '  Na:  ', NA, '  ZETA:   ', ZETA
!      WRITE(6,'(A,3I5)') 'Nx  Ny  NZ  (steps):    ', NI(1), NI(2), NI(3)
      XCR = XC
      YCR = YC
      ZCR = ZC
      CALL CRYS2CART(XCR,YCR,ZCR)
!      WRITE(6,'(A,3F10.6)') 'XCR  YCR  ZCR  (cart) :    ', XCR, YCR, ZCR
!      WRITE(6,*)
!      WRITE(6,*) 'Points inspected in grid:'

!     Integrating over the grid!  #### POOR METHOD ####
      DO IZ=1, 2*NI(3)

         DO IY=1, 2*NI(2)

            DO IX=1, 2*NI(1)

              ! Calculating the coords of electron in crys coords
              ! This coordinates are always in the mapped grid.
              xe = XC + del_x*DFLOAT(IX - NI(1))
              ye = YC + del_y*DFLOAT(IY - NI(2))
              ze = ZC + del_z*DFLOAT(IZ - NI(3))

!              WRITE(6,'(A,3F10.6)') 'elec-cryst:   ', xe, ye, ze

              ! Converting the electron coordinates from crystal to cartesian
              xec = xe
              yec = ye
              zec = ze
              CALL CRYS2CART(xec, yec, zec)

!              WRITE(6,'(A,3F10.6)') 'elec-Carte:   ', xec, yec, zec

              ! Centering in atom
              xec = xec - XA
              yec = yec - YA
              zec = zec - ZA

              ! Calculating the hermite gaussian value
              r_eA2 = (xec)**2 + (yec)**2 + (zec)**2
              HG_S = HERGA_S(NA, ZETA, r_eA2)

              ! Putting the elec. coords. in unitary cube
              ! ***** This subroutine has IF structures inside!!
              CALL PUTINGRID(xe, ye, ze)

              ! Calculating the index of elec coords. in cube to get the RHO value
              CALL CRYS2INDEX(xe,ye,ze, IRx,IRy,IRz)

!              WRITE(6,'(A,3F10.6,A,3I5)') 'xe ye ze:  ', xe,ye,ze, 'IRx IRy IRz:  ', IRx,IRy,IRz

              ! Getting the RHO value
              rhoI = RHO(IRx, IRy, IRz)

!              CALL CRYS2CART(xe,ye,ze)
!              WRITE(6,'(A,3I10)') 'Grid-Idx      ', IRx, IRy, IRz
!              WRITE(6,'(A,6F20.10)') 'Cart-Grid:    ', xec, yec, zec, HG_S*rhoI, HG_S, rhoI
!              WRITE(6,*)

              ! Passing the current calculation to the total result
              INT_X(IX) = HG_S*rhoI

              ! Accounting the number of inspections of the integration
              IGRID = IGRID +1
            ENDDO

            ! INTEGRATING X
            CALL INTNUMCALL(INT_X, 2*NI(1), INT_TEMP)
            INT_Y(IY) = INT_TEMP

         ENDDO

         ! INTEGRATING Y
         CALL INTNUMCALL(INT_Y, 2*NI(2), INT_TEMP)
         INT_Z(IZ) = INT_TEMP

      ENDDO

      ! INTEGRATING Z
      CALL INTNUMCALL(INT_Z, 2*NI(3), INT_TEMP)
      RESULT = INT_TEMP


!     Multiplying by the volume element
      RESULT = RESULT*VOLUME_EL

      ! Deallocating ...
      DEALLOCATE(INT_X, INT_Y, INT_Z)

!
!     *** End of SUBROUTINE HERGAUSS_S ***
!
      END SUBROUTINE HERGAUSS_S



