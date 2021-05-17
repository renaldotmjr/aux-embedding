     SUBROUTINE RHOFIT(ri, fitrho)
!
!     Purpose: Calculate RHO in 'ri' position with Hermite-Gaussians
!               centered on atoms
!
!     History: - Creation (14.02.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     ri(xe,ye,ze)        : The electron position to calculate the rho
!     IATOM            : Run over the atoms
!     ISET             : Run over the hermite gaussian sets
!     I, inorm         : run over the normalization coefficient of each set
!     ZETA             : The exponent
!     NA
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

!     List of OUTPUT variables:
!     fitrho         : The fitted RHO calculated in 'ri'

!     List of local variables:
!
!     List of USED variables from EMBPARAM:
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     List of Local variables
      INTEGER       :: IATOM, ISET, inorm, I, ii
      REAL (8)      :: fitrho, rho_f, ZETA, HG_f, xe, ye, ze, r_eA2
      REAL (8)      :: HERGA_S, HERGA_P, HERGA_D1, HERGA_D2
      REAL (8)      :: COORD_elec_UC(NTRANS3,4), XYZ_ATOM(3)

!     List of local ALLOCATABLE variables
      REAL (8), DIMENSION(3)  :: ri
      ! In this case, NA will not be allocated for optimize this calculation
      REAL (8), DIMENSION(14) :: NA, COEF  ! Normalization coef.

      rho_f = 0.0D0

!     Loop over the atoms
      DO IATOM=1, NATOM

         DO ii=1, NTRANS3

             XYZ_ATOM(1) = X0(IATOM)
             XYZ_ATOM(2) = Y0(IATOM)
             XYZ_ATOM(3) = Z0(IATOM)

             CALL TRANSATOM(XYZ_ATOM, TRANSLATVEC(ii,:))

             COORD_elec_UC(ii,1) = ri(1) - XYZ_ATOM(1)
             COORD_elec_UC(ii,2) = ri(2) - XYZ_ATOM(2)
             COORD_elec_UC(ii,3) = ri(3) - XYZ_ATOM(3)

             COORD_elec_UC(ii,4) = COORD_elec_UC(ii,1)**2 + &
                                 & COORD_elec_UC(ii,2)**2 + &
                                 & COORD_elec_UC(ii,3)**2
         ENDDO


         DO ii=1, NTRANS3

             xe = COORD_elec_UC(ii,1)
             ye = COORD_elec_UC(ii,2)
             ze = COORD_elec_UC(ii,3)
             r_eA2 = COORD_elec_UC(ii,4)

             DO ISET=LAUXPTR(IATOM), UAUXPTR(IATOM)

                ZETA = ZET(ISET)

                ! Passing the normalization coefficients to NA
                I=1
                do inorm=LAUXSHL(ISET), UAUXSHL(ISET)
                    NA(I) = EMBNORMA(inorm)
                    COEF(I) = EMB_COEF(inorm,1)
                    I = I+1
                enddo
                !------------------------------------------------

                HG_f = 0.0D0

                CALL DRVHGF(LTOTAL(ISET), COEF, ZETA, NA, xe, ye, ze, r_eA2, HG_f)

                rho_f = rho_f + HG_f

             ENDDO

         ENDDO

      ENDDO

      fitrho = rho_f

!
!     *** End of SUBROUTINE RHOFIT ***
!
      END SUBROUTINE RHOFIT
