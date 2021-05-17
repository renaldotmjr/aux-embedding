     SUBROUTINE CHECKSAMEATM(DIF_ZATOM, DIF_IATOM, NDIFATOMS)
!
!     Purpose: Check the same atom types inside the unit cell
!
!     History: - Creation (11.06.12, MJR)
!              - Inclusion of DIF_IATOM (30/04/2013, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!     NATOM                   : Number of atoms            (INTEGER)
!     ZATOM                   : Atomic numbers             (INTEGER)
!
!     List of local variables:
!     DIF_ZATOM(NATOM)        : Atomic number of diferent atoms
!     NDIFATOMS               : Number of diferent atoms
!
!     ******************************************************************

      USE EMBPARAM

      IMPLICIT NONE

      ! List of local variables
      INTEGER     :: DIF_ZATOM(NATOM), DIF_IATOM(NATOM)
      INTEGER     :: NDIFATOMS

      INTEGER     :: I_ZATOM, J_DIFATOM, EQUIV_count

      !initializing the variables
      NDIFATOMS = 0

      ! passing the first Zatom to DIF_ZATOM(1)
      DIF_ZATOM(1) = ZATOM(1)
      DIF_IATOM(1) = 1
      NDIFATOMS = 1
      J_DIFATOM = 1

      ! Loop starting in the seccond atom in ZATOM
      DO I_ZATOM = 2, NATOM

          ! initializing the equivalence counter
          EQUIV_count = 0

          DO J_DIFATOM = 1, NDIFATOMS

              IF ( ZATOM(I_ZATOM) .EQ. DIF_ZATOM(J_DIFATOM) ) THEN
                  EQUIV_count = 1
                  EXIT
              ENDIF

          ENDDO

          IF ( EQUIV_count .EQ. 0 ) THEN
              NDIFATOMS = NDIFATOMS + 1
              DIF_ZATOM(NDIFATOMS) = ZATOM(I_ZATOM)
              DIF_IATOM(NDIFATOMS) = I_ZATOM

          ENDIF

      ENDDO


!
!     *** End of SUBROUTINE CHECKSAMEATM ***
!
      END SUBROUTINE CHECKSAMEATM
