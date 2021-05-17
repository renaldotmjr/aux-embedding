     SUBROUTINE PRINTAUXIS(DIF_ZATOM, DIF_IATOM, NDIFATOMS, F_shrink)
!
!     Purpose: Check the same atom types inside the unit cell
!
!     History: - Creation (11.06.12, MJR)
!              - Inclusion of F_shrink (30/04/2013, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!
!     List of local variables:
!     DIF_ZATOM(NATOM)        : Atomic number of diferent atoms
!     DIF_IATOM(NATOM)        : Index of diferent atoms, the first in order
!     NDIFATOMS               : Number of diferent atoms
!
!     ******************************************************************

      USE EMBPARAM
      USE INPDEFVARS
      USE global_variables

      IMPLICIT NONE

      ! List of local variables
      INTEGER            :: DIF_ZATOM(NDIFATOMS), DIF_IATOM(NDIFATOMS)
      INTEGER            :: NDIFATOMS

      CHARACTER (80)     :: TMP

      REAL (8)           :: ZETA, Nel_FIT_atm
      REAL (8)           :: Qd_fit               ! Atomic charge fitted originally
      REAL (8)           :: F_shrink             ! The F shrinkage factor
      REAL (8)           :: Zd_F                 ! The new Zd

      INTEGER     :: I, J, K, Zcurr, NEatom, NFUNC, LTOT
      INTEGER     :: IN                 ! Number of AUXIS function blocks


      ! Opening the deMon_CUB_NAME file for writing
      OPEN(24,FILE=deMon_CUB_NAME, STATUS='OLD', POSITION='APPEND')

      ! Loop over the non-equal atoms
      DO I=1, NDIFATOMS

         ! Opening the basis file
         OPEN(1,FILE=AUXIS_NAME, STATUS='OLD')
!         if( FILE NOT EXIST?????)

      10 READ(1,"(A)",END=20) TMP

         IF(TMP(1:1).NE.'#') THEN
            ! Converting string to integer
            read(TMP(1:3), "(I5)") Zcurr

!           Searching the atomic number in basis file
            IF( Zcurr == DIF_ZATOM(I) ) THEN

               ! Reading the number of electrons treated by the basis
               read(TMP(9:11), "(I5)") NEatom

!
!              CALCULATING THE NEW N ELEC TO PRINT IN .cub
!
               ! Calculating the number of electrons treated by the auxis
               CALL Nelec_FITATM(I, Nel_FIT_atm)
               ! Calculating the original fitted charge Qd_fit
               Qd_fit = DFLOAT(NEatom) - Nel_FIT_atm
               ! Calculating the new Zd (NEatom)
               Zd_F = Qd_fit + F_shrink*Nel_FIT_atm
!
!              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

               ! Writing at the deMon_CUB_NAME file
               ! Writing the atom label and the number of treated electrons
               WRITE(24, "(A,F8.4)") ELSYM(DIF_ZATOM(I)), Zd_F

               ! Reading the number of L-set
               READ(1,*) IN

               ! Writing at the deMon_CUB_NAME file
               ! Writing the number of AUXIS blocks
               WRITE(24, "(I5)") IN

               do J=1, IN
                  ! Reading the number of L-functions and the L value
                  READ(1,*) NFUNC, LTOT

                  ! Writing at the deMon_CUB_NAME file
                  ! Writing the total angular momentum and the number of functions in this block
                  WRITE(24, "(2I5)") LTOT, NFUNC

                  do K=1, NFUNC
                     ! Reading Zeta
                     READ(1,*) ZETA

                     ! Writing at the deMon_CUB_NAME file
                     ! Writing the the ZETA values
                     WRITE(24, "(F20.9)") ZETA

                  enddo

               enddo

            ENDIF ! End of IF that finds the atom in file

         ENDIF ! End of IF that check the '#' comment chacartere

         GOTO 10
      20 CONTINUE

         ! Closing the basis file
         CLOSE(1)

      ENDDO  ! End of non-equal atoms DO

      ! Closing the deMon_CUB_NAME file
      CLOSE(24)

!
!     *** End of SUBROUTINE PRINTAUXIS ***
!
      END SUBROUTINE PRINTAUXIS
