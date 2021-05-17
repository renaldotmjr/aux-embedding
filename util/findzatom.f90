      SUBROUTINE FINDZATOM
!
!     Purpose: Get the atomic number from physcon.f90
!
!     History: - Creation (02.02.2012, MJR)
!
!     ******************************************************************
!
!
!
!     ******************************************************************

      use global_variables
      USE EMBPARAM
      use String_Manipulation

      IMPLICIT NONE

!     Local variables
      INTEGER       :: I, J, C_test
      CHARACTER(4)  :: Str_atom      ! The local name of input atom
      CHARACTER(4)  :: ElemSymb      ! The local name of parameter atom
      ! Initializing ZATOM to zero
      ZATOM = 0

      ! Do over the atoms
      DO I = 1, NATOM

          ! Starting the C_test for testing if the atom was found
          C_test = -1

          Str_atom = SYMB(I)

          ! ### Pre-processing the Str_atom string ###
          ! - Adjusting to the left
          Str_atom = ADJUSTL(Str_atom)
          ! - Converting from upper to lower case
          CALL LowCase(Str_atom)
          ! - Converting eventual numbers to space character
          CALL NumSpace(Str_atom)

          ! Do over the atoms in the ELSYM() array
          do J = 0, 103

              ! Copying ELSYM(J) to ElemSymb
              ElemSymb = ELSYM(J)
              ! - Converting from upper to lower case
              CALL LowCase(ElemSymb)

              ! Testing if the strings are equal
              IF( Str_atom .EQ. ElemSymb ) THEN
                  ZATOM(I) = J
                  C_test = 1
                  EXIT
              END IF

          end do

          ! Testing if the atom was found
          IF( C_test == -1 ) THEN

              WRITE(*,*) 'ATOM ', I, SYMB(I), '(', Str_atom ,')', 'NOT FOUND!'
              STOP 'ATOM NOT FOUND'

          END IF

      END DO

      END SUBROUTINE FINDZATOM
