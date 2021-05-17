      SUBROUTINE RHO_PWSCF(FILE_NAME)
!
!     Purpose: Read the rho values in the PWscf output
!
!     History: - Creation (28.06.11, CO)
!              - Modification (19.09.11, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     FILE_NAME    : PWscf output file name
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     NATOM                   : Number of atoms
!     NX, NY, NZ              : Grid in a, b and c
!     A_LATx, A_LATy, A_LATz  : Lattice parameter a
!     B_LATx, B_LATy, B_LATz  : Lattice parameter b
!     C_LATx, C_LATy, C_LATz  : Lattice parameter c
!     ATOM                    : Atomic symbols
!     X0(), Y0(), Z0()        : Atomic coordinates
!     RHO(,,,)                : Mapped density RHO(X,Y,Z)
!     X() , Y() , Z()         : Mapped grid coordinates
!
!     List of local variables:
!     TMP                     : The line in the PWscf output reading
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     Temporary variable to FILE_NAME reading
      CHARACTER (80)                          :: FILE_NAME
      CHARACTER (80)                          :: TMP

      INTEGER  :: I, IX, IY, IZ

      ! Opening the PWscf output
      OPEN(1,FILE=FILE_NAME, STATUS='OLD')  ! pwscf output - post process

   10 READ(1,"(A)",END=20) TMP

!     Reading the lattice parameters
      IF(TMP(2:8).EQ.'PRIMVEC') THEN
         READ(1,*) A_LATx, A_LATy, A_LATz
         READ(1,*) B_LATx, B_LATy, B_LATz
         READ(1,*) C_LATx, C_LATy, C_LATz
      ENDIF
!
!     Reading the symbol and atomic coordinates of the atoms in the unit cell
      IF(TMP(2:10).EQ.'PRIMCOORD') THEN
         READ(1,*) NATOM
         ALLOCATE (X0(NATOM), Y0(NATOM), Z0(NATOM))
         ALLOCATE (SYMB(NATOM),ZATOM(NATOM))
         DO I=1,NATOM
            READ(1,"(A4,2X,3F15.9)") &
            &  SYMB(I), X0(I), Y0(I), Z0(I)
        ENDDO
      ENDIF
!
!     Reading RHO and the grid
      IF(TMP(1:19).EQ.'DATAGRID_3D_UNKNOWN') THEN
         READ(1,*) NX, NY, NZ
         READ(1,*)
         READ(1,*)
         READ(1,*)
         READ(1,*)
         
         ALLOCATE (RHO(NX,NY,NZ))
         
         READ(1,*) ( ( (RHO(IX,IY,IZ), IX=1,NX ) , IY=1,NY), IZ=1,NZ)

      ENDIF
!
      GOTO 10
!
   20 CONTINUE
      
      ! Closing the PWscf file
      CLOSE(1)
!
!     *** End of SUBROUTINE RHO_PWSCF ***
!
      END SUBROUTINE RHO_PWSCF


