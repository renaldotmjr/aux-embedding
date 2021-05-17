     SUBROUTINE PRINTEXTEMBCH(F_shrink)
!
!     Purpose: Print the Extended embedding charges
!               The file with index 24 must be open in this subroutone call
!
!     History: - Creation (11.06.12, MJR)
!              - Inclusion of F_shrink (30/04/2013, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!
!     List of local variables:
!
!     ******************************************************************

      USE EMBPARAM
      USE INPDEFVARS
      USE global_variables

      IMPLICIT NONE

      ! List of local variables
      INTEGER        :: i, ITRANS, IATOM, ISET, ishl, iprt, LEN_STR

      REAL (8)       :: XYZ_ATOMA(3), VEC_PRT(NSHL), UNITFAC
      REAL (8)       :: F_shrink             ! The F shrinkage factor

      CHARACTER (80) :: FNAME


      ! Mounting the translation vectors
      CALL MOUNTTVECS

      ! Checking the unit system!
!     -----------------------------------------------------------------------------------
      LEN_STR = LEN_TRIM(UNITSY)
      IF (UNITSY(1:LEN_STR) .EQ. 'ANGSTROM') THEN

          UNITFAC = AU2ANG

      ELSEIF (UNITSY(1:LEN_STR) .EQ. 'BOHR') THEN

          UNITFAC = 1.0D0

      ELSE
          WRITE(6,*) 'Not supported unit system: ', UNITSY(1:LEN_STR)
          STOP
      ENDIF
!     -----------------------------------------------------------------------------------


      ! Opening the deMon_CUB_NAME file for writing
      OPEN(24,FILE=deMon_CUB_NAME, STATUS='OLD', POSITION='APPEND')

!!!!!!!!!!!!
      FNAME = 'CHARGE.cub'
      OPEN(64,FILE=FNAME, STATUS='REPLACE')
!!!!!!!!!!!

      DO ITRANS=1, SIZE_TV

          DO IATOM=1, NATOM

              ! Translating each atom in the unit cell
              XYZ_ATOMA(1) = X0(IATOM)
              XYZ_ATOMA(2) = Y0(IATOM)
              XYZ_ATOMA(3) = Z0(IATOM)
              do i=1,3
                  XYZ_ATOMA(i) = XYZ_ATOMA(i) + TRANSVECS(ITRANS, i)
              enddo

              ! Printing the Label and coordinates
              WRITE (24,"(A5,3F20.9)") ELSYM( ZATOM(IATOM) ), UNITFAC*XYZ_ATOMA(1),&
                                      & UNITFAC*XYZ_ATOMA(2), UNITFAC*XYZ_ATOMA(3)

!!!!!!!!!!!!!!!!!!!!!!!
              ! Printing the Label and coordinates
              WRITE (64,"(3F20.9, A5)") UNITFAC*XYZ_ATOMA(1), UNITFAC*XYZ_ATOMA(2), UNITFAC*XYZ_ATOMA(3), &
                                       & ELSYM( ZATOM(IATOM) )
!!!!!!!!!!!!!!!!!!!!!!!

              iprt = 0

              ! Run over the AUXPTR of IATOMA
              DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                  ! Run over the AUXSHL of UXPTR in IATOMA
                  do ishl=LAUXSHL(ISET), UAUXSHL(ISET)

                      iprt = iprt +1
                      ! Including the F_shrink...
                      VEC_PRT(iprt) = F_shrink * EMBNORMA(ishl)*EMB_COEF(ishl,1)

                  enddo
              ENDDO

              CALL PRINTVEC(24, VEC_PRT, iprt, 1)


          ENDDO

      ENDDO

      ! Closing the deMon_CUB_NAME file
      CLOSE(24)


!!!!!!!!!!!!!!!!!!!
      CLOSE(64)
!!!!!!!!!!!!!!!!!!!

!
!     *** End of SUBROUTINE PRINTEXTEMBCH ***
!
      END SUBROUTINE PRINTEXTEMBCH
