     SUBROUTINE INTAUXRHO
!
!     Purpose: Calculate numerical integrals in mapped grid electron
!              density
!
!     History: - Creation (28.09.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     IATOM            : Run over the atoms
!     ISET             : Run over the hermite gaussian sets
!     I, inorm         : run over the normalization coefficient of each set
!     iint             : Run over the results to indexation
!     ZETA             : The exponent
!     NORMA()          : Receives the normalization coef.
!     INTEGRAL()       : Receives the integral results
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
!
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     VEC_b()  = <i_A|rho>    : Vector with the numerical integrals over
!                                unit cell GRID
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     List of Local variables
      INTEGER       :: IATOM, ISET, inorm, I, iint, IGRID
      REAL (8)      :: ZETA

!     List of local ALLOCATABLE variables
      REAL (8), DIMENSION(:), ALLOCATABLE  :: NORMA, &     ! Normalization coef.
                                              &  INTEGRAL


!     Allocating VEC_b
      ALLOCATE( VEC_b(NSHL+1) )
      ALLOCATE( VEC_b_GRID(NSHL) )

      ! Starting VEC_b with zero
      VEC_b = 0.0D0
      ! The NSHL+1 vallue in VEC_b is 'TNelec', the total electron number
      VEC_b(NSHL+1) = DFLOAT(TNelec)

!     Loop over the atoms
      DO IATOM=1, NATOM

         DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

            ZETA = ZET(ISET)

!           Allocating NORMA and INTEGRAL with the size in ORBTOTAL(ISET)
            ALLOCATE( NORMA(ORBTOTAL(ISET)), INTEGRAL(ORBTOTAL(ISET)) )


            ! Starting the with zero the vectors
            INTEGRAL = 0.0D0

            ! Passing the normalization coefficients to NORMA
            I=1
            do inorm=LAUXSHL(ISET), UAUXSHL(ISET)
               NORMA(I) = EMBNORMA(inorm)
               I = I+1
            enddo
            !------------------------------------------------

            ! ##################  S-TYPE integrals block ##################
            IF ( LTOTAL(ISET) == 0 ) THEN
               ! Calculating the integrall
               CALL HERGAUSS_S(NORMA, ZETA, &
               & X0(IATOM), Y0(IATOM), Z0(IATOM), INTEGRAL, IGRID)
            ! ############### END of S-TYPE integrals block ###############

            ! ##################  P-TYPE integrals block ##################
            ELSEIF ( LTOTAL(ISET) == 1 ) THEN
               ! Calculating the integrall
               CALL HERGAUSS_P(NORMA, ZETA, &
               & X0(IATOM), Y0(IATOM), Z0(IATOM), INTEGRAL, IGRID)
            ! ############### END of P-TYPE integrals block ###############

            ! ##################  D-TYPE integrals block ##################
            ELSEIF ( LTOTAL(ISET) == 2 ) THEN
               ! Calculating the integrall
               CALL HERGAUSS_D(NORMA, ZETA, &
               & X0(IATOM), Y0(IATOM), Z0(IATOM), INTEGRAL, IGRID)
            ENDIF
            ! ############### END of D-TYPE integrals block ###############


            ! --- Passing the result to VEC_b ---
            I=1
            do iint=LAUXSHL(ISET), UAUXSHL(ISET)
               ! Passing the result
               VEC_b(iint) = INTEGRAL(I)

               ! Passing the number of
               VEC_b_GRID(iint) = IGRID

               I = I+1
            enddo
            !------------------------------------

!           Deallocating NORMA and INTEGRAL
            DEALLOCATE(NORMA, INTEGRAL)

         ENDDO

      ENDDO

!
!     *** End of SUBROUTINE INTAUXRHO ***
!
      END SUBROUTINE INTAUXRHO
