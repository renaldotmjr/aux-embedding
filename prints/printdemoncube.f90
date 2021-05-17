     SUBROUTINE PRINTDEMONCUBE(F_shrink)
!
!     Purpose: Print the <deMon.cub> (with the used-defined name)
!
!     History: - Creation (11.06.12, MJR)
!              - Inclusion of F_shrink (30/04/2013, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!
!     List of local variables:
!     DIF_ZATOMS(NATOM)      : Atomic number of diferent atoms
!     NDIFATOM               : Number of diferent atoms
!
!     ******************************************************************

      USE EMBPARAM
      USE INPDEFVARS

      IMPLICIT NONE

      ! List of local variables
      INTEGER     :: DIF_ZATOMS(NATOM), DIF_IATOM(NATOM)
      INTEGER     :: NDIFATOM

      INTEGER     :: I, LEN_STR

      REAL (8)    :: F_shrink             ! The F shrinkage factor

      ! Initializing the variables


      ! Opening the deMon_CUB_NAME file for writing
      ! If the file already exists, it will be deleted
      OPEN(24,FILE=deMon_CUB_NAME, STATUS='REPLACE')

      WRITE (24,"(A)") '#'
      WRITE (24,"(A)") 'AUXILIARY FUNCTIONS'

      ! Closing deMon_CUB_NAME
      CLOSE(24)

      ! Generating a vector with the non-equal atomic numbers
      CALL CHECKSAMEATM(DIF_ZATOMS, DIF_IATOM, NDIFATOM)

      ! Printing the AUXIS functions in .cub file
      CALL PRINTAUXIS(DIF_ZATOMS, DIF_IATOM, NDIFATOM, F_shrink)

      ! Opening the deMon_CUB_NAME file for writing
      ! If the file already exists, it will be deleted
      OPEN(24,FILE=deMon_CUB_NAME, STATUS='OLD', POSITION='APPEND')


      LEN_STR = LEN_TRIM(UNITSY)
      WRITE (24,"(A)") '#'
      WRITE (24,"(2A)") 'EMBEDDING DENSITIES ', UNITSY(1:LEN_STR)

      ! Printing the riplicated coodinates and coefficients
      ! This function do not print the coordinates of oniginal unit cell
      CALL PRINTEXTEMBCH(F_shrink)

      ! Closing deMon_CUB_NAME
      CLOSE(24)

!
!     *** End of SUBROUTINE PRINTDEMONCUBE ***
!
      END SUBROUTINE PRINTDEMONCUBE
