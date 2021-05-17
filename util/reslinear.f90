    SUBROUTINE RESLINEAR
!
!     Purpose: Resolution of linear system equations
!
!     History: - Creation (04.10.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!
!     List of USED variables in EMBPARAM:
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!     VEC_b()  = <i_A|rho>    : Vector with the numerical integrals over
!                                unit cell GRID
!
!     ******************************************************************
!     *** SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ):
!     Arguments
!     =========
!
!     N       (input) INTEGER
!             The number of linear equations, i.e., the order of the
!             matrix A.  N >= 0.
!
!     NRHS    (input) INTEGER
!             The number of right hand sides, i.e., the number of columns
!             of the matrix B.  NRHS >= 0.
!
!     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!             On entry, the N-by-N coefficient matrix A.
!             On exit, the factors L and U from the factorization
!             A = P*L*U; the unit diagonal elements of L are not stored.
!
!     LDA     (input) INTEGER
!             The leading dimension of the array A.  LDA >= max(1,N).
!
!     IPIV    (output) INTEGER array, dimension (N)
!             The pivot indices that define the permutation matrix P;
!             row i of the matrix was interchanged with row IPIV(i).
!
!     B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!             On entry, the N-by-NRHS matrix of right hand side matrix B.
!             On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!     LDB     (input) INTEGER
!             The leading dimension of the array B.  LDB >= max(1,N).
!
!     INFO    (output) INTEGER
!             = 0:  successful exit
!             < 0:  if INFO = -i, the i-th argument had an illegal value
!             > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                   has been completed, but the factor U is exactly
!                   singular, so the solution could not be computed.
!
!     ******************************************************************
!
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

      ! Local variables:
      INTEGER                               :: I, N, NRHS, LDA, LDB, INFO
      INTEGER , DIMENSION(:),   ALLOCATABLE :: IPIV
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: A

      ! Allocating EMB_COEF
      ALLOCATE( EMB_COEF(NSHL+1, 1), A(NSHL+1, NSHL+1), IPIV(NSHL+1) )

      ! Passing the VEC_b to EMB_COEF
      do I=1, NSHL+1
         EMB_COEF(I,1) = VEC_b(I)
      enddo
      ! Passing the parameters for use in DGESV
      N = NSHL+1
      NRHS = 1
      LDA = NSHL+1
      LDB = NSHL+1

      A = OVERLAP_AB

      ! Calling the subroutine
      CALL DGESV( N, NRHS, A, LDA, IPIV, EMB_COEF, LDB, INFO )

!     ******************* DEBUG print! *******************
      IF (DEBUG_INV) THEN

         WRITE(6,*) '***********************************************'
         WRITE(6,*) '**** The Linear Equation System Resolution ****'
         WRITE(6,*)
         WRITE(6,*) 'INFO =  ', INFO
         WRITE(6,*) 'N =     ', NSHL+1
         WRITE(6,*) 'NRHS =  ', 1
         WRITE(6,*) 'LDA =   ', NSHL+1
         WRITE(6,*) 'LDB =   ', NSHL+1
         WRITE(6,*)
         WRITE(6,*) 'EMB_COEF    Pivoting'
         DO I=1, N
            WRITE(6,'(I5, F25.20, I5)') I, EMB_COEF(I,1), IPIV(I)
         ENDDO
         WRITE(6,*)
         WRITE(6,*) '**** END-OF Linear Equation System Resolution ****'
         WRITE(6,*) '**************************************************'

      ENDIF
!     **************** END-OF DEBUG print! ***************

!
!     *** End of SUBROUTINE RESLINEAR ***
!
      END SUBROUTINE RESLINEAR
