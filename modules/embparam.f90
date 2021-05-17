      MODULE EMBPARAM
!     PARAMETER DEFINITION FOR EMBEDDING CODE
!
!     Purpose: EMBedding PARAMeters definitions
!
!     History: - Creation (19.09.11, MJR)
!
!     ******************************************************************
!
!     List of PARAMETER
!     AU2ANG           : Conversion factor from a.u. to angs
!     ANG2AU           : Conversion factor from angs to a.u.
!     MAXFAC           : Maximum factor argument
!     DIFTHRESHOLD     : Defining a zero for the calculation of differences
!
!     List of static variables
!     PI                      : Pi number
!     DEGRAD                  : DEG/RAD conversion factors
!     RADDEG                  : DEG/RAD conversion factors
!     FAC(A)                  : Vector with factorial of A
!     DFAC(A)                 : Vector with double factorial of A
!     NATOM                   : Number of atoms            (INTEGER)
!     NX, NY, NZ              : Grid in a, b and c         (INTEGER)
!     A_LATx, A_LATy, A_LATz  : Lattice parameter a        (REAL)
!     B_LATx, B_LATy, B_LATz  : Lattice parameter b        (REAL)
!     C_LATx, C_LATy, C_LATz  : Lattice parameter c        (REAL)
!     LATVEC(3,3)             : The lattice vectors in a filed (REAL)
!     LATVECNORM(3)           : The lattice vec. normas    (REAL)
!     MAXVECNORM              : The greater value          (REAL)
!     MINVECNORM              : The smaler value          (REAL)
!     del_x, del_y, del_z     : Steps of cryst. coords.    (REAL)
!     MLAT(3,3),MLAT_inv(3,3) : Lattice matrices for conversions
!     VOLUME_UC               : Unit Cell volume           (REAL)
!     VOLUME_EL               : Volume element for numerical integration (REAL)
!     TNelec                  : Total number of electrons accounted
!                                in current basis set
!
!     NSET                    : Number of aux basis functions (INTEGER)
!     NSHL                    : Total number of primitives hermite
!                                gaussian functions           (INTEGER)
!
!     List of ALLOCATABLE variables
!     SYMB                    : Atomic symbols             (CHARACTER)
!     ZATOM                   : Atomic numbers             (INTEGER)
!     Nelec                   : Number of electrons accounted by
!                                the basis set of each atom
!     X0, Y0, Z0              : Atomic coordinates         (REAL)
!     RHO(,,,)                : Mapped density RHO(X,Y,Z)  (REAL)
!     X() , Y() , Z()         : Mapped grid coordinates    (REAL)
!
!       -> Auxiliary basis block variables
!       ZET()                 : Aux. basis zeta of each atomic orbital
!       LTOTAL()              : Aux. basis total angular momentum of
!                               each atomic orbital
!       ORBTOTAL()            : Number of primitive hermite gaussian
!                               functions of each aux func
!       LAUXPTR               : Lower pointer to aux function
!       UAUXPTR               : Upper pointer to aux function
!       LAUXSHL               : Lower pointer to aux function shells
!       UAUXSHL               : Upper pointer to aux function shells
!       AUXL(,3)              : Matrix with all primitives hermite
!                               gaussian function angular momentum index:
!                                (1) -> Lax
!                                (2) -> Lay
!                                (3) -> Laz
!       EMBNORMA()            : Normalization coefficients of each shell
!
!     List of the EMBEDDING MATRIX and VECTORS
!     VEC_b()  = <i_A|rho>    : Vector with the numerical integrals over
!                                unit cell GRID
!     VEC_b_GRID()            : Vector with the number of inspectioned
!                                points in grid for each integral
!     EMB_COEF                : The Embedding coefficients
!     OVERLAP_AB(,)           : Matrix with normalized overlap integrals
!
!
!     ******************************************************************
!
      REAL(8)      :: AU2ANG, ANG2AU

      PARAMETER(AU2ANG = 0.529177249D0)   ! Conversion factor from a.u. to angs
      PARAMETER(ANG2AU = 1.0D0/AU2ANG)    ! Conversion factor from angs to a.u.

      REAL(8)      :: DIFTHRESHOLD        ! defining a zero for the calculation of differences

      REAL(8)      :: PI, DEGRAD, RADDEG

      ! Maximum double factor argument
      INTEGER :: MAXFAC
      PARAMETER(MAXFAC = 20)

      REAL(8)      :: FAC(0:MAXFAC), DFAC(-1:2*MAXFAC+1)

!     Number of atoms and electrons in the Embedding cluster
      INTEGER      :: NATOM, TNelec

!     Number of aux basis functions
      INTEGER      :: NSET, &       ! Number of basis functions
                    & NSHL          ! Number of primitives hermite gaussian functions

!     Grid parameters and boundary condition parameter
      INTEGER      :: NX, NY, NZ

!     Lattice parameters
      REAL (8)     :: A_LATx, A_LATy, A_LATz
      REAL (8)     :: B_LATx, B_LATy, B_LATz
      REAL (8)     :: C_LATx, C_LATy, C_LATz
      REAL (8)     :: LATVEC(3,3), LATVECNORM(3)
      REAL (8)     :: MAXVECNORM, MINVECNORM
      REAL (8)     :: del_x, del_y, del_z
      REAL (8)     :: MLAT(3,3), MLAT_inv(3,3)
      REAL (8)     :: VOLUME_UC, VOLUME_EL

      INTEGER      :: NTRANS, NTRANS3
      REAL (8), DIMENSION(:,:),    ALLOCATABLE ::  TRANSLATVEC

!     Atomic symbols and atomic coordinates in unit cell 
      CHARACTER (4), DIMENSION(:), ALLOCATABLE :: SYMB              ! Atomic symbol
      INTEGER      , DIMENSION(:), ALLOCATABLE :: ZATOM, Nelec      ! Atomic number
      REAL (8)     , DIMENSION(:), ALLOCATABLE :: X0, Y0, Z0        ! Atomic coordinates

!     Mapped density on grid (x,y,z)
      REAL (8), DIMENSION(:,:,:),  ALLOCATABLE :: RHO               ! RHO(X,Y,Z)
      REAL (8), DIMENSION(:),      ALLOCATABLE :: X , Y , Z         ! Mapped grid

!     Zeta vector and total angular momentum vector
      REAL (8), DIMENSION(:),      ALLOCATABLE :: ZET               ! Zeta of each aux function
      INTEGER , DIMENSION(:),      ALLOCATABLE :: LTOTAL, ORBTOTAL  ! Total L of each aux function

!     Pointers to the lower and upper ref of aux function of each atom
      INTEGER , DIMENSION(:),      ALLOCATABLE :: LAUXPTR, UAUXPTR  ! Total L of each aux function

!     Pointers to the lower and upper ref of aux function of each atom
      INTEGER , DIMENSION(:),      ALLOCATABLE :: LAUXSHL, UAUXSHL  ! Total L of each aux function

!     Matrix with all primitives hermite gaussian function angular momentum index:
      INTEGER , DIMENSION(:,:),    ALLOCATABLE :: AUXL              ! Lax, Lay and Laz of each aux
!                                                                     function

!     Vector with normalization coefficients of each shell
      REAL (8), DIMENSION(:),      ALLOCATABLE :: EMBNORMA          ! Zeta of each aux function

!     List of the EMBEDDING MATRIX and VECTORS
      INTEGER , DIMENSION(:),      ALLOCATABLE :: VEC_b_GRID        ! Number of inspectioned points
      REAL (8), DIMENSION(:),      ALLOCATABLE :: VEC_b             ! Vector with the numerical integrals
      REAL (8), DIMENSION(:,:),    ALLOCATABLE :: OVERLAP_AB        ! Vector with the numerical integrals
      REAL (8), DIMENSION(:,:),    ALLOCATABLE :: EMB_COEF
!
!     *** End of MODULE EMBPARAM ***
!
      END MODULE EMBPARAM


