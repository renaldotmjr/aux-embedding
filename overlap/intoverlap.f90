     SUBROUTINE INTOVERLAP
!
!     Purpose: INTegral OVERLAP matrix calculation
!
!     History: - Creation (30.09.11, MJR)
!              - Including CCM build of OVERLAP (31.05.12, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     ZETA, ZETB       :
!     XA,YA,ZA         : Atomic coordinates of center A.
!     XB,YB,ZB         : Atomic coordinates of center B.
!     IATOM, JATOM     :
!     ISET, JSET       :
!     ishl, jshl, i, j :
!     NORM_S           : The normalization integral value
!     RIJ(3)
!     RIJNORM
!
!     List of USED variables in EMBPARAM:
!     NATOM            : Number of atoms
!     ZET()            : Aux. basis zeta of each atomic orbital
!     LTOTAL()         : Aux. basis total angular momentum of
!                        each atomic orbital
!     ORBTOTAL()       : Number of primitive hermite gaussian
!                        functions of each aux func
!     LAUXPTR          : Lower pointer to aux function
!     UAUXPTR          : Upper pointer to aux function
!     LAUXSHL          : Lower pointer to aux function shells
!     UAUXSHL          : Upper pointer to aux function shells
!     X0, Y0, Z0       : Atomic coordinates
!     EMBNORMA()       : Normalization coefficients of each shell
!
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     OVERLAP_AB(,)    : Matrix with normalized overlap integrals
!
!     ******************************************************************
!
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

!     List of Local variables
      INTEGER       :: IATOM, JATOM, ISET, JSET, ishl, jshl, i, j
      REAL (8)      :: ZETA, ZETB, XA,YA,ZA, XB,YB,ZB, NORMA, NORM_S

      ! Local for the diagonalization
      INTEGER                               :: N, NDIM, IEGEN, NR
      INTEGER                               :: II, NVAL
      INTEGER , DIMENSION(:),   ALLOCATABLE :: IQ
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: H, EIGVEC
      REAL (8), DIMENSION(:),   ALLOCATABLE :: EIGVAL, XXX


      REAL (8)      :: RIJ(3), RIJNORM, XYZ_TATOM(3)

!     The matrix that will receive the overlap result in each block [i_A|j_B]
      REAL (8), DIMENSION(:,:), ALLOCATABLE ::  OVERLAP


!     Allocating the overlap integral matrix
      ALLOCATE( OVERLAP_AB(NSHL+1, NSHL+1) )

!     Creating the lattice vectors and the norms
      CALL LATTICEVEC

!     Initializing the matrix with zero
      OVERLAP_AB = 0.0D0

      ! Run over Atom I
      DO IATOM=1, NATOM

         ! Passing I atomic coordinates
         XA = X0(IATOM)
         YA = Y0(IATOM)
         ZA = Z0(IATOM)

         ! Run over Atom J
         DO JATOM=1, NATOM

             ! Passing J atomic coordinates
             XB = X0(JATOM)
             YB = Y0(JATOM)
             ZB = Z0(JATOM)

             !*********************************************************************
             !**** DOING THE CYCLIC CLUSTER MODEL BUILND OF THE OVERLAP MATRIX ****

             ! Generating the distance vector and calculating the norm
             RIJ(1) = XB - XA
             RIJ(2) = YB - YA
             RIJ(3) = ZB - ZA
             RIJNORM = RIJ(1)**2 + RIJ(2)**2 + RIJ(3)**2   ! Squared distance


             ! The first if: Same atoms or if Rij < |a|/2
             IF ( (IATOM .EQ. JATOM) .OR. (RIJNORM .LT. ((MINVECNORM/4.0D0)-DIFTHRESHOLD) ) ) THEN

                 ! Run over the AUXPTR of Atom I
                 DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                     ! Getting the ZET value for ISET of atom I
                     ZETA = ZET(ISET)

                     ! Run over the AUXPTR of Atom J
                     DO JSET=LAUXPTR(JATOM), UAUXPTR(JATOM)

                         ! Getting the ZET value for JSET of atom J
                         ZETB = ZET(JSET)

                         ! Allocating the result matrix OVERLAP
                         ALLOCATE( OVERLAP( ORBTOTAL(ISET), ORBTOTAL(JSET) ) )

                         ! Calling the subroutine that calls the integration
                         CALL DRVOVERLAP(LTOTAL(ISET), LTOTAL(JSET), ZETA,ZETB, &
                                       & XA,YA,ZA, XB,YB,ZB, OVERLAP, &
                                       & ORBTOTAL(ISET), ORBTOTAL(JSET) )

                         !------ Passing the result to OVERLAP_AB ------
                         !  *** Passing the result and normalizing ***
                         ! Run over ISET and JSET in AUXSHL
                         i = 1
                         do ishl=LAUXSHL(ISET), UAUXSHL(ISET)
                         j = 1
                         do jshl=LAUXSHL(JSET), UAUXSHL(JSET)
                            OVERLAP_AB(ishl, jshl) = EMBNORMA(ishl)*EMBNORMA(jshl)*OVERLAP(i,j)
                            j = j +1
                         enddo
                         i = i + 1
                         enddo
                         !----------------------------------------------

!                        ******************* DEBUG print! *******************
                         IF (DEBUG) THEN
                             write (6, '(A, A,I5, A, A,I5)') '  IATOM  ', SYMB(IATOM), &
                                              & IATOM, '  JATOM  ', SYMB(JATOM), JATOM
                             write (6, *) 'ISET  JSET  LT(ISET) LT(JSET)  ORBT(ISET) ORBT(JSET)'
                             write (6, '(I4,I6,I7,I9,I11,I10)') ISET, JSET, LTOTAL(ISET), &
                                          & LTOTAL(JSET), ORBTOTAL(ISET), ORBTOTAL(JSET)
                             do i=1, ORBTOTAL(ISET)
                                 do j=1, ORBTOTAL(JSET)
                                     write (6, "(A, I5, I5, F25.20)") 'i   j   OVERLAP  ', &
                                              & i, j, OVERLAP(i,j)
                                 enddo
                             enddo
                             write (6,*)
                             write (6,*)
                         ENDIF
!                        **************** END-OF DEBUG print! ***************

                         ! Deallocating OVERLAP for this block
                         DEALLOCATE( OVERLAP )

                     ENDDO ! End-do of AUXPTR in JATOM
                 ENDDO  ! End-do of AUXPTR in IATOM

             ! The seccond if: Rij = |a|/2 where 'a' is anyone of LATVEC
             ELSEIF ( ( (RIJNORM - (LATVECNORM(1)/4.0D0)) < DIFTHRESHOLD ) .OR. &
             &         ( (RIJNORM - (LATVECNORM(2)/4.0D0)) < DIFTHRESHOLD ) .OR. &
             &         ( (RIJNORM - (LATVECNORM(3)/4.0D0)) < DIFTHRESHOLD )   ) THEN

                 ! Run over the AUXPTR of Atom I
                 DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                     ! Getting the ZET value for ISET of atom I
                     ZETA = ZET(ISET)

                     ! Run over the AUXPTR of Atom J
                     DO JSET=LAUXPTR(JATOM), UAUXPTR(JATOM)

                         ! Getting the ZET value for JSET of atom J
                         ZETB = ZET(JSET)

                         ! Allocating the result matrix OVERLAP
                         ALLOCATE( OVERLAP( ORBTOTAL(ISET), ORBTOTAL(JSET) ) )

                         ! Calculating the atom at real unit cell
                         ! Calling the subroutine that calls the integration
                         CALL DRVOVERLAP(LTOTAL(ISET), LTOTAL(JSET), ZETA,ZETB, &
                                       & XA,YA,ZA, XB,YB,ZB, OVERLAP, &
                                       & ORBTOTAL(ISET), ORBTOTAL(JSET) )

                         !------ Passing the result to OVERLAP_AB ------
                         !  *** Passing the result and normalizing ***
                         ! Run over ISET and JSET in AUXSHL
                         i = 1
                         do ishl=LAUXSHL(ISET), UAUXSHL(ISET)
                             j = 1
                             do jshl=LAUXSHL(JSET), UAUXSHL(JSET)
                                 OVERLAP_AB(ishl, jshl) = EMBNORMA(ishl)*EMBNORMA(jshl)*(0.5D0*OVERLAP(i,j))
                                 j = j +1
                             enddo
                             i = i + 1
                         enddo
                         !----------------------------------------------
                         ! Deallocating OVERLAP for this block
                         DEALLOCATE( OVERLAP )
                     ENDDO ! End-do of AUXPTR in JATOM
                 ENDDO  ! End-do of AUXPTR in IATOM

                 ! Calculating the atom at the virtual cell

                 ! Calculating the translation base vector!
                 CALL CALCBASTRV(RIJ)

                 XYZ_TATOM(1) = XB
                 XYZ_TATOM(2) = YB
                 XYZ_TATOM(3) = ZB

                 ! Translating the atom coords. using the base translation vector RIJ
                 CALL TRANSATOM(XYZ_TATOM, RIJ)

                 XB = XYZ_TATOM(1)
                 YB = XYZ_TATOM(2)
                 ZB = XYZ_TATOM(3)

                 ! Run over the AUXPTR of Atom I
                 DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                     ! Getting the ZET value for ISET of atom I
                     ZETA = ZET(ISET)

                     ! Run over the AUXPTR of Atom J
                     DO JSET=LAUXPTR(JATOM), UAUXPTR(JATOM)

                         ! Getting the ZET value for JSET of atom J
                         ZETB = ZET(JSET)

                         ! Allocating the result matrix OVERLAP
                         ALLOCATE( OVERLAP( ORBTOTAL(ISET), ORBTOTAL(JSET) ) )

                         ! Calculating the overlap between the IATOM and translated JATOM
                         CALL DRVOVERLAP(LTOTAL(ISET), LTOTAL(JSET), ZETA,ZETB, &
                                       & XA,YA,ZA, XB,YB,ZB, OVERLAP, &
                                       & ORBTOTAL(ISET), ORBTOTAL(JSET) )

                         !------ Passing the result to OVERLAP_AB ------
                         !  *** Passing the result and normalizing ***
                         ! Run over ISET and JSET in AUXSHL
                         i = 1
                         do ishl=LAUXSHL(ISET), UAUXSHL(ISET)
                             j = 1
                             do jshl=LAUXSHL(JSET), UAUXSHL(JSET)
                                OVERLAP_AB(ishl, jshl) = OVERLAP_AB(ishl, jshl) + &
                                                       & EMBNORMA(ishl)*EMBNORMA(jshl)*(0.5D0*OVERLAP(i,j))
                                j = j +1
                             enddo
                             i = i + 1
                         enddo
                         !----------------------------------------------

!                        ******************* DEBUG print! *******************
                         IF (DEBUG) THEN
                             write (6, '(A, A,I5, A, A,I5)') '  IATOM  ', SYMB(IATOM), &
                                              & IATOM, '  JATOM  ', SYMB(JATOM), JATOM
                             write (6, *) 'ISET  JSET  LT(ISET) LT(JSET)  ORBT(ISET) ORBT(JSET)'
                             write (6, '(I4,I6,I7,I9,I11,I10)') ISET, JSET, LTOTAL(ISET), &
                                              & LTOTAL(JSET), ORBTOTAL(ISET), ORBTOTAL(JSET)
                             do i=1, ORBTOTAL(ISET)
                                 do j=1, ORBTOTAL(JSET)
                                     write (6, "(A, I5, I5, F25.20)") 'i   j   OVERLAP  ', &
                                              & i, j, OVERLAP(i,j)
                                 enddo
                             enddo
                             write (6,*)
                             write (6,*)
                         ENDIF
!                        **************** END-OF DEBUG print! ***************

                         ! Deallocating OVERLAP for this block
                         DEALLOCATE( OVERLAP )
                     ENDDO ! End-do of AUXPTR in JATOM
                 ENDDO  ! End-do of AUXPTR in IATOM

             ELSE    ! If Rij > |a|/2

                 ! Calculating the atom at the virtual cell

                 ! Calculating the translation base vector!
                 CALL CALCBASTRV(RIJ)

                 XYZ_TATOM(1) = XB
                 XYZ_TATOM(2) = YB
                 XYZ_TATOM(3) = ZB

                 ! Translating the atom coords. using the base translation vector RIJ
                 CALL TRANSATOM(XYZ_TATOM, RIJ)

                 XB = XYZ_TATOM(1)
                 YB = XYZ_TATOM(2)
                 ZB = XYZ_TATOM(3)

                 ! Run over the AUXPTR of Atom I
                 DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                     ! Getting the ZET value for ISET of atom I
                     ZETA = ZET(ISET)

                     ! Run over the AUXPTR of Atom J
                     DO JSET=LAUXPTR(JATOM), UAUXPTR(JATOM)

                         ! Getting the ZET value for JSET of atom J
                         ZETB = ZET(JSET)

                         ! Allocating the result matrix OVERLAP
                         ALLOCATE( OVERLAP( ORBTOTAL(ISET), ORBTOTAL(JSET) ) )

                         ! Calculating the overlap between the IATOM and translated JATOM
                         CALL DRVOVERLAP(LTOTAL(ISET), LTOTAL(JSET), ZETA,ZETB, &
                                       & XA,YA,ZA, XB,YB,ZB, OVERLAP, &
                                       & ORBTOTAL(ISET), ORBTOTAL(JSET) )

                         !------ Passing the result to OVERLAP_AB ------
                         !  *** Passing the result and normalizing ***
                         ! Run over ISET and JSET in AUXSHL
                         i = 1
                         do ishl=LAUXSHL(ISET), UAUXSHL(ISET)
                             j = 1
                             do jshl=LAUXSHL(JSET), UAUXSHL(JSET)
                                 OVERLAP_AB(ishl, jshl) = EMBNORMA(ishl)*EMBNORMA(jshl)*OVERLAP(i,j)
                                 j = j +1
                             enddo
                             i = i + 1
                         enddo
                         !----------------------------------------------

!                        ******************* DEBUG print! *******************
                         IF (DEBUG) THEN
                             write (6, '(A, A,I5, A, A,I5)') '  IATOM  ', SYMB(IATOM), &
                                              & IATOM, '  JATOM  ', SYMB(JATOM), JATOM
                             write (6, *) 'ISET  JSET  LT(ISET) LT(JSET)  ORBT(ISET) ORBT(JSET)'
                             write (6, '(I4,I6,I7,I9,I11,I10)') ISET, JSET, LTOTAL(ISET), &
                                              & LTOTAL(JSET), ORBTOTAL(ISET), ORBTOTAL(JSET)
                             do i=1, ORBTOTAL(ISET)
                                 do j=1, ORBTOTAL(JSET)
                                     write (6, "(A, I5, I5, F25.20)") 'i   j   OVERLAP  ', &
                                              & i, j, OVERLAP(i,j)
                                 enddo
                             enddo
                             write (6,*)
                             write (6,*)
                         ENDIF
!                        **************** END-OF DEBUG print! ***************

                         ! Deallocating OVERLAP for this block
                         DEALLOCATE( OVERLAP )
                     ENDDO ! End-do of AUXPTR in JATOM
                 ENDDO  ! End-do of AUXPTR in IATOM

             ENDIF


         ENDDO ! End-do of JATOM
     ENDDO ! End-do of IATOM

!     Calculating the normalization integrals <a>
!      These integrals appear in the last row and col of OVERLAP_AB
!     ** This loop could be done in the nested loop above, but I prefered to put separated

!     Loop over each ZET()
      DO ISET=1, NSET

!        Getting the NORMA, ZETA and LC values
         NORMA = EMBNORMA(LAUXSHL(ISET))
         ZETA = ZET(ISET)

!        If a /= s, <a> = 0
         IF (LTOTAL(ISET) == 0) THEN

            ! Calculating the integral
            CALL NORMALIZATION_S(NORMA, ZETA, NORM_S)

            ! Passing the value to the line
            OVERLAP_AB(NSHL+1, LAUXSHL(ISET)) = NORM_S

            ! Passing the value to the column
            OVERLAP_AB(LAUXSHL(ISET), NSHL+1) = OVERLAP_AB(NSHL+1, LAUXSHL(ISET))

         ENDIF
      ENDDO


      ! Initializing variables
      NR = 0
      EIGVAL = 0.0D0
      EIGVEC = 0.0D0
      IQ = 0.0D0

      ! Passing the parameters for use in JACOBI
      N = NSHL
      NDIM = N
      IEGEN = 0

      ALLOCATE( H(NSHL, NSHL) )
      DO i=1, NSHL
          DO j=1, NSHL

              H(i,j) = OVERLAP_AB(i,j)

          ENDDO
      ENDDO

      ! Allocating the eigenvalues and eigenvectors vectors
      ALLOCATE( EIGVAL(NSHL), EIGVEC(NSHL, NSHL) )
      ALLOCATE( XXX(NSHL), IQ(NSHL) )

!     Diagonalizing OVERLAP_AB and getting eigenvalues and eigenvectors
      CALL JACOBI(H, N, NDIM, IEGEN, EIGVAL, NR, IQ, EIGVEC, XXX, .FALSE.)

      ! Testing if we some negative eigenvalue
      NVAL=0
      DO II=1, NDIM

          IF (EIGVAL(II) < 0.0D0) THEN
              NVAL = NVAL + 1
          ENDIF

      ENDDO

      IF (NVAL .NE. 0) THEN
          WRITE(6,*) 'Overlap matrix not positive defined!'
          WRITE(6,*) 'Means that it has ', NVAL, 'negative eigenvalue(s)'
          WRITE(6,*) '*** Consider to expand the unit cell! ***'
      ELSE
          WRITE(6,*) 'Overlap matrix is positive defined!'
      ENDIF

!     ******************* DEBUG print! *******************
      IF (DEBUG_INV) THEN

          WRITE(6,*) '*************************************************'
          WRITE(6,*) '**** Overlap diagonalization (JACOBI Method) ****'
          WRITE(6,*)
          WRITE(6,*) 'DIAGONALIZATION:'
          WRITE(6,*) ' Eigenvalues:'
          WRITE(6,*)
          DO II=1, NDIM
              WRITE(6,"(F18.10)") EIGVAL(II)
          ENDDO
          WRITE(6,*)


      ENDIF
!     **************** END-OF DEBUG print! ***************

      ! Deallocating local variables
      DEALLOCATE(H, EIGVAL, EIGVEC, XXX, IQ)


!
!     *** End of SUBROUTINE INTOVERLAP ***
!
      END SUBROUTINE INTOVERLAP
