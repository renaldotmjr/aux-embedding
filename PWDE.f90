      PROGRAM PWDE
!
!  Program: Plane-Waves Density Embedding
!  Date:
!  Authors:
!           Renaldo Tenorio de Moura Jr.
!           Claudio de Oliveira
!           Helio Anderson Duarte
!           Andreas M Koster
!
!  Objective:
!   This program uses the density calculated from PWSCF and evaluated
!   in an uniform grid to be projected on the auxiliary basis set
!   (Hermite Gaussian functions).
!
!  Modules used:

      USE EMBPARAM
      USE EMBTOL
      USE global_variables
      USE INPDEFVARS

      IMPLICIT NONE

!  Local variables:
!
!     Temporary variables to FILE_NAME reading

      REAL (8) :: DELTA_X, DELTA_Y, DELTA_Z, DELTA_VOLUME, SUM_RHO, p_ERROR, ME, SD, Nel_FIT, Nel_FIT_atm
      REAL (8), DIMENSION(:,:),    ALLOCATABLE :: INT_NUM_ANALY
      INTEGER  :: I, J, K, IX, IY, IZ, countprt, LEN_STR

      ! Getting the input file name
      CALL getarg ( 1, INPUT_NAME )

      WRITE(6,*) 'Input file:', INPUT_NAME

      DEBUG = .FALSE.
      DEBUG_INV = .FALSE.


!     Defining the default options
      CALL DEFAULTS

!     Reading the input file and setting the options
      CALL READINPUT(INPUT_NAME)

!     Calling the physical constants
      CALL PHYSCON

!     Calling tha MATHFUN
      CALL MATHFUN

!     Reading RHO, the grid, the lattice parameters and atomic coordinates
!         from PWscf output.
!     This subroutine allocates dinamically RHO(,,,), X(), Y() and Z().
      WRITE(6,*) 'PW_OUT_NAME file:', PW_OUT_NAME

      CALL RHO_PWSCF(PW_OUT_NAME)

!     Checking the unit system!
!     -----------------------------------------------------------------------------------
      LEN_STR = LEN_TRIM(UNITSY)
      IF (UNITSY(1:LEN_STR) .EQ. 'ANGSTROM') THEN
          WRITE(6,*) 'Input unit system: ', UNITSY(1:LEN_STR)
          WRITE(6,*) 'Converting to: BOHR'
          WRITE(6,*)
          ! Converting units to au
          CALL UNITCONV
      ELSEIF (UNITSY(1:LEN_STR) .EQ. 'BOHR') THEN
          WRITE(6,*) 'Input unit system: ', UNITSY(1:LEN_STR)
      ELSE
          WRITE(6,*) 'Not supported unit system: ', UNITSY(1:LEN_STR)
      ENDIF
!     -----------------------------------------------------------------------------------


!     Calculating the matrices that will be used in crystal-cartesian conversions
      CALL LATCONV

!     Calculating the grid map coordinates (coordinates in a.u.)
      CALL EMB_GRID

!     Finding the atomi numbers
      CALL FINDZATOM

!     This subroutine reads the basis in BASIS_NAME file and
!          mounts LAUXPTR(), UAUXPTR(), ZET(), LTOTAL()
      IF ( (AUXIS_STD .EQ. 'OPEN') .OR. (AUXIS_STD .EQ. 'deMon') ) THEN

          LEN_STR = LEN_TRIM(AUXIS_NAME)
          WRITE(6,*)
          WRITE(6,*) 'Reading the AUXIS functions from: ', AUXIS_NAME(1:LEN_STR)
          CALL GETBASIS(AUXIS_NAME)

      ELSE
          LEN_STR = LEN_TRIM(AUXIS_STD)
          WRITE(6,*)
          WRITE(6,*) 'Not supported AUXIS type file: ', AUXIS_STD(1:LEN_STR)

      ENDIF



      CALL GENSHL

!     Normalize the AO aux functions
      CALL EMBNORM

!     Calculate the unit cell vomule and the volume element (in a.u.)
      CALL VOLUME

!     Calculating the numerical integrals
      WRITE(6,*)
      WRITE(6,*) 'Integrating over mapped density to obtain '
      WRITE(6,*) ' the total mapped number of electrons... '
      WRITE(6,*)

      !#########################################################################################
!     Integracao de RHO - somente soma (pobre)
      SUM_RHO = 0.0D0
      DO IX=1,NX
       DO IY=1,NY
        DO IZ=1,NZ
         SUM_RHO = SUM_RHO + RHO(IX,IY,IZ)
        ENDDO
       ENDDO
      ENDDO
      SUM_RHO = SUM_RHO*VOLUME_EL

!     Integracao de RHO - somente soma (pobre)
!      DO IX=1,NX
!         DO IY=1,NY
!            DO IZ=1,NZ
!               RHO(IX,IY,IZ) = ( DFLOAT(TNelec)/SUM_RHO )*RHO(IX,IY,IZ)
!            ENDDO
!         ENDDO
!      ENDDO

!     Integracao de RHO - somente soma (pobre)
!      SUM_RHO = 0.0D0
!      DO IX=1,NX
!       DO IY=1,NY
!        DO IZ=1,NZ
!         SUM_RHO = SUM_RHO + RHO(IX,IY,IZ)
!        ENDDO
!       ENDDO
!      ENDDO
!      SUM_RHO = SUM_RHO*VOLUME_EL

      !#########################################################################################
      WRITE(6,*) 'Done!'

      WRITE(6,*)
      WRITE(6,*) 'Calculating the numerical overlaps between the  '
      WRITE(6,*) ' mapped density and the HG functions... '
      WRITE(6,*)

      ! Calculating the numerical integrals
      CALL INTAUXRHO

      WRITE(6,*) 'Done!'


!     Calculating the NORMALIZED overlap matrix
      CALL INTOVERLAP

!     Solving the Linear Equation System
      CALL RESLIJAC


!     Doing the FIT TESTS
!     -----------------------------------------------------------------------------------
      LEN_STR = LEN_TRIM(FITTEST)
      IF (FITTEST(1:LEN_STR) .EQ. 'YES') THEN

          WRITE(6,*)
          WRITE(6,*) 'Doing the tests of fitting coefficients... '

          ! The number of expansions to the fit analysis
          NTRANS = 4
          ! Doing the fit error analysis
          CALL FITERROR(p_ERROR, ME, SD)
          ! Calculating the analytical number of electrons
          CALL Nelec_FIT(Nel_FIT)

          WRITE(6,*) 'Done!'

          WRITE(6,*)
          WRITE(6,*) '------------------------ FIT TEST ANALYSIS ------------------------'
          WRITE(6,*)
          WRITE(6,*) 'EMB_COEF used'
          DO I=1, (NSHL+1)
              WRITE(6,'(I5, F25.20)') I, EMB_COEF(I,1)
          ENDDO
          WRITE(6,*)
          WRITE(6,*) 'Fit test analysis:'
          WRITE(6,"(A,F15.9)") 'p_ERROR = ', p_ERROR
          WRITE(6,"(A,F15.9)") 'Mean Error = ', ME
          WRITE(6,"(A,F15.9)") 'SD Error = ', SD
          WRITE(6,*)
          WRITE(6,*) 'Analytical integration of number of electrons:'
          WRITE(6,"(A,F27.20)") 'Fitted N_elec = ', Nel_FIT
          WRITE(6,*)
          WRITE(6,*) '-------------------------------------------------------------------'

      ELSEIF (FITTEST(1:LEN_STR) .EQ. 'NO') THEN
          WRITE(6,*) 'Any fit test will be done!'
      ELSE
          WRITE(6,*) 'Invalid option for fit tests: ', FITTEST(1:LEN_STR)
      ENDIF
!     -----------------------------------------------------------------------------------



!##################################################################################################
      Nel_FIT = 0

      WRITE(6,*)
      WRITE(6,*) 'ATOMIC COORDINATES AND FITTED CHARGES'
      ! Calculating the Nelec for each atom!
      DO I=1, NATOM

          CALL Nelec_FITATM(I, Nel_FIT_atm)
          Nel_FIT = Nel_FIT + Nel_FIT_atm

          WRITE(6,"(I3,3X,A2,I5,3F15.9,F10.5)") I, SYMB(I),ZATOM(I),X0(I),Y0(I),Z0(I), Nel_FIT_atm

      ENDDO

      WRITE(6,*)
      WRITE(6,*) 'Total fitted charge', Nel_FIT

!##################################################################################################




      CALL PRINTDEMONCUBE(F_shr)










      ! The number of expansions to the fit analysis
      NTRANS = 4

      WRITE(6,*)
      WRITE(6,*) 'Printing .cube file'
      CUBE_NAME = 'CUBE.cube'
      CALL GENCUBE(CUBE_NAME)

      ALLOCATE( INT_NUM_ANALY(NSHL+1, 1) )
      INT_NUM_ANALY = MATMUL(OVERLAP_AB, EMB_COEF)
      WRITE(6,*)
      WRITE(6,*) 'Numerical integrals calculated and from LE system'
      DO I=1, (NSHL+1)
         WRITE(6,'(I5, 3F25.20, F7.2)') I, VEC_b(I), INT_NUM_ANALY(I,1), &
            & ( DABS(VEC_b(I)) - DABS(INT_NUM_ANALY(I,1))), &
            & 100.0D0*( DABS(VEC_b(I)) - DABS(INT_NUM_ANALY(I,1))) / DABS(INT_NUM_ANALY(I,1))
      ENDDO

!
      WRITE(6,*)
      WRITE(6,"(A,F15.6,A)") 'Volume  = ', VOLUME_UC, &
     &   ' [a.u.]^3'
      WRITE(6,*)
      WRITE(6,"(A,3F15.9,A,I5)") 'A_LAT = {', A_LATx, A_LATy, A_LATz, '}   NX = ', NX
      WRITE(6,"(A,3F15.9,A,I5)") 'B_LAT = {', B_LATx, B_LATy, B_LATz, '}   NY = ', NY
      WRITE(6,"(A,3F15.9,A,I5)") 'C_LAT = {', C_LATx, C_LATy, C_LATz, '}   NZ = ', NZ
      WRITE(6,*)

      WRITE(6,"(A,F15.9,A)") 'SUM_RHO = ', SUM_RHO, ' electrons'

!
      WRITE(6,"(A,I5)") 'NATOM = ', NATOM
      DO I=1,NATOM
         WRITE(6,"(I3,3X,A2,I5,3F15.9,2I5)") I, SYMB(I),ZATOM(I),X0(I),Y0(I),Z0(I), LAUXPTR(I), UAUXPTR(I)
      ENDDO
!
      WRITE(6,*)
      WRITE(6,"(A,I10)") 'NSHL = ', NSHL
      WRITE(6,*)
      DO I=1,NSET
         WRITE(6,"(I5,3X,F15.9,4I7)") I, ZET(I), LTOTAL(I), ORBTOTAL(I), LAUXSHL(I), UAUXSHL(I)
      ENDDO

      DO I=1,NSHL
         WRITE(6,"(I5,3X,3I10,3F20.14, I20)") I, AUXL(I,1), AUXL(I,2), AUXL(I,3), &
                                          & EMBNORMA(I), VEC_b(I), VEC_b_GRID(I)
      ENDDO

      WRITE(6,"(A,F15.9)") 'VEC_b(NSHL+1) = ', VEC_b(NSHL+1)

      WRITE(6,*) 'Overlap Matrix:'
      WRITE(6,*)

      DO I=1, NSHL
         DO J=1, NSHL
            WRITE(6,"(I5,I5,F20.15)") I, J, OVERLAP_AB(I,J)
         ENDDO
      ENDDO





      CALL EMBCUBREADAUXIS




      ! Deallocating everything
      CALL DEALLOC_EMB

      END PROGRAM PWDE
