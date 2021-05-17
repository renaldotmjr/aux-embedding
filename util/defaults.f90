     SUBROUTINE DEFAULTS
!
!     Purpose: Set the default options.
!
!     History: - Creation (06.06.12, MJR)
!
!     ******************************************************************
!
!     List of USED variables from EMBPARAM:
!     DIFTHRESHOLD     : Defining a zero for the calculation of differences

!
!
!     List of USED variables from EMBTOL:
!     INTTOL           : Tolerance for numerical integrations!
!     XTOL             : Tolerance for iterativelly solving the
!                          numerical integrations limits
!     XTMAX            : Maximum of iterations in limits calculation
!
!     List of USED variables from global_variables:
!
!
!     List of USED variables from INPDEFVARS:
!     UNITSY               : Input unit system option
!     FITTEST              : Fitting tests option
!     EQUIVALENCY          : Equivalency between atoms input option
!
!
!
!     ******************************************************************

      USE EMBPARAM
      USE EMBTOL
      USE global_variables
      USE INPDEFVARS

      IMPLICIT NONE

      DIFTHRESHOLD = 1.0D-6         ! Defining a zero for the calculation of differences

      INTTOL = 1.0D-14              ! Tolerance for numerical integrations!
      XTOL = 1.0D-6                 ! Tolerance for iterativelly solving the
                                    !  numerical integrations limits
      XTMAX = 10000                 ! Maximum of iterations in limits calculation
      NTRANS = 4                    ! The number of expansions to the fit analysis
      F_shr = 1.0                   ! The shrinkage factor!

      UNITSY = 'ANGSTROM'           ! OPTIONS < BOHR ; ANGSTROM >
      FITTEST = 'NO'                ! OPTIONS < YES ; NO >
      EQUIVALENCE = 'SAMEEQ'        ! OPTIONS < ALLNONEQ ; SAMEEQ ; GEOMLABEL >

      PRT_deMon_CUBE = 'NO'         ! OPTIONS < YES ; NO>
      REPLICATION = 'NO'            ! OPTIONS < YES ; NO>

      AUXIS_STD = 'OPEN'            ! OPTIONS < deMon ; OPEN >

      REPNUM = 0                    ! Initializing the field with the number of
                                    !  replications in each lattice vec.

      BTAKEUCOUT = 'NO'             ! OPTIONS < YES ; NO >
      NTAKEUCOUT = 0                ! Initializing the number of UC that will be token out

!
!     *** End of SUBROUTINE DEFAULTS ***
!
      END SUBROUTINE DEFAULTS
