      SUBROUTINE GETBASIS(BASIS_NAME)
!
!     Purpose: Mount the vector with auxiliary basis set zeta values
!                -> This subroutine reallocates some vectors
!
!     History: - Creation (23.09.11, MJR)
!              - Inclusion of deMon standard in MOUNTZET.F90
!                         (06.08.12, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     BASIS_NAME  : Auxiliary basis set name
!
!     List of local variables:
!     BASIS_NAME
!     TMP
!     ZETA
!     Zcurr, IN, NFUNC, LTOT
!     I J, K,
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     LAUXPTR               : Lower pointer to aux function
!     UAUXPTR               : Upper pointer to aux function
!
!     ******************************************************************
!     Called subroutines:
!     MOUNT_ZET(ZETA, LTOT  : in mountzet.f90
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     Temporary variable to BASIS_NAME reading
      CHARACTER (80)               :: BASIS_NAME
      CHARACTER (80)               :: TMP
 
      REAL (8)                     :: ZETA

      INTEGER                      :: I, J, K, Zcurr, IN, NFUNC, LTOT

!     Initializing NSET and Zcurr
      NSET = 0
      Zcurr = -1

!     Allocating LAUXPTR UAUXPTR
      ALLOCATE( LAUXPTR(NATOM), UAUXPTR(NATOM), Nelec(NATOM) )

      DO I=1, NATOM

         ! Passing the LAUXPTR vallue
         LAUXPTR(I) = NSET + 1

         ! Opening the basis file
         OPEN(1,FILE=BASIS_NAME, STATUS='OLD')
!         if( FILE NOT EXIST?????)

      10 READ(1,"(A)",END=20) TMP

         IF(TMP(1:1).NE.'#') THEN
            ! Converting string to integer
            read(TMP(1:3), "(I5)") Zcurr

!           Searching the atomic number in basis file
            IF( Zcurr == ZATOM(I) ) THEN

               ! Reading the number of electrons treated by the basis
               read(TMP(9:11), "(I5)") Nelec(I)

!               WRITE(6, "(A,2I3)") 'Found: I, Zcurr =', I, Zcurr
               ! Reading the number of L-set
               READ(1,*) IN
               
               do J=1, IN
                  ! Reading the number of L-functions and the L value
                  READ(1,*) NFUNC, LTOT

                  do K=1, NFUNC
                     ! Reading Zeta
                     READ(1,*) ZETA
                     ! Passing Zeta and Ltot to MOUNT_ZET
                     ! This subroutine increments the NSET value! The NSET value
                     !    will be passed below to UAUXPTR
                     CALL MOUNT_ZET(ZETA, LTOT)
                  enddo

               enddo

            ENDIF ! End of IF that finds the atom in file

         ENDIF ! End of IF that check the '#' comment chacartere

         GOTO 10
      20 CONTINUE

         ! Closing the basis file
         CLOSE(1)

         ! Testing if the current atom has a basis in file! If no, STOP!
         IF(Zcurr .EQ. -1) THEN

             WRITE(6,"(A,A5)") 'Missing AUXIS for: ', ZATOM(I)
             STOP

         ENDIF

         ! Passing the UAUXPTR vallue
         UAUXPTR(I) = NSET

      ENDDO  ! End of NATOM DO

      ! Calculating the total number of electrons
      TNelec = SUM(Nelec)
!
!     *** End of SUBROUTINE GETBASIS ***
!
      END SUBROUTINE GETBASIS

