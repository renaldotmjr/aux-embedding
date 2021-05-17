      MODULE INPDEFVARS
!
!     Purpose: INPut DEFinitions of some VARiables
!
!     History: - Creation (29.05.12, MJR)
!
!     ******************************************************************
!
!     List of Variables
!     INPUT_NAME           : The PWDE input file name
!     PW_OUT_NAME          : The PWScf output file name
!     BASIS_NAME           : The basis file name    
!     CUBE_NAME            : The .cube file name    
!     deMon_CUB_NAME       : The deMon.cub file name
!
!     UNITSY               : Input unit system option
!     FITTEST              : Fitting tests option
!     EQUIVALENCY          : Equivalency between atoms input option
!
!     ******************************************************************
!

      CHARACTER (80) :: INPUT_NAME
      CHARACTER (80) :: PW_OUT_NAME
      CHARACTER (80) :: AUXIS_NAME
      CHARACTER (80) :: CUBE_NAME
      CHARACTER (80) :: deMon_CUB_NAME

      CHARACTER (80) :: UNITSY
      CHARACTER (80) :: FITTEST
      CHARACTER (80) :: EQUIVALENCE
      CHARACTER (80) :: PRT_deMon_CUBE
      CHARACTER (80) :: REPLICATION

      CHARACTER (80) :: AUXIS_STD

      CHARACTER (80) :: BTAKEUCOUT

      REAL(8)        :: F_shr

      INTEGER        :: REPNUM(3,2), SIZE_TV, NTAKEUCOUT
      REAL (8), DIMENSION(:,:),    ALLOCATABLE ::  TRANSVECS
      INTEGER , DIMENSION(:,:),    ALLOCATABLE ::  TAKEUCOUTvecs

!
!     *** End of MODULE DEFVARS ***
!
      END MODULE INPDEFVARS
