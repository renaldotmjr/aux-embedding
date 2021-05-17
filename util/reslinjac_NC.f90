    SUBROUTINE RESLIJAC_NC
!
!     Purpose: Resolution of linear system equations using jacobi routine
!              (WITHOUT CONSTRAIN CONDITION: S_rho = N_electrons)
!
!     History: - Creation (13.02.12, MJR, CO)
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
!    ######### THE SUBROUTINE JACOBI WAS TAKEN IN deMok-2K code #########
!
!    SUBROUTINE JACOBI(H,N,NDIM,IEGEN,EIGVAL,NR,IQ,U,X)
!C
!C     Purpose: JACOBI computes eigenvalues and eigenvectors of a real
!C              symmetric matrix using the Jacobi method.
!C
!C     History: - Creation (June 1994, B. Ahlswede; the original program
!C                          was written by F. Corbato and M. Merwin)
!C
!C              - Implemented in deMon (19.10.00, AMK)
!C
!C     ******************************************************************
!C
!C     On input:
!C
!C             H: Contains the matrix to be diagonalized.
!C                Only the upper triangle of the matrix need to
!C                be supplied. This part is changed on exit,
!C                but the diagonal elements are unchanged.
!C             N: Order of matrix H.
!C          NDIM: Dimension of H as declared in the calling program.
!C         IEGEN: Must be 0 if eigenvectors are to be computed.
!C
!C     On output:
!C
!C        EIGVAL: Contains the eigenvalues in ascending order.
!C             U: Contains the corresponding eigenvectors.
!C            NR: Number of Jacobi-rotations.
!
!C     *** X(I) contains the largest element in ith column           ***
!C     *** IQ(I) holds second subscript defining position of element ***
!C
!C     ------------------------------------------------------------------
!C
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

      !
      ! Local variables:
      !
      INTEGER                               :: N, NDIM, IEGEN, NR
      INTEGER                               :: I
      INTEGER , DIMENSION(:),   ALLOCATABLE :: IQ
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: H, EIGVEC
      REAL (8), DIMENSION(:),   ALLOCATABLE :: EIGVAL, XXX

      ! The transformed matrix
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: EIGVEC_T
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: EMB_COEF_NC
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: B_trf     ! Transformed VEC_b
      REAL (8), DIMENSION(:,:), ALLOCATABLE :: X_trf     ! Transformed EMB_COEF

      ! Linear resolution


      ! Allocating EMB_COEF and the matriz to be diagonalized
!      ALLOCATE( EMB_COEF(NSHL+1, 1))
      ALLOCATE( H(NSHL, NSHL) )
      ! Allocating the eigenvalues and eigenvectors vectors
      ALLOCATE( EIGVAL(NSHL), EIGVEC(NSHL, NSHL) )
      ! Allocating transformed vectors
      ALLOCATE( B_trf(NSHL, 1), EIGVEC_T(NSHL, NSHL) )

      ! Allocating EMB_COEF
      ALLOCATE( EMB_COEF_NC(NSHL, 1) )
      ALLOCATE( X_trf(NSHL, 1) )


      ALLOCATE( XXX(NSHL), IQ(NSHL) )

      ! Initializing variables
      NR = 0
      EIGVAL = 0.0D0
      EIGVEC = 0.0D0
      X = 0.0D0
      IQ = 0.0D0

      ! Passing the parameters for use in JACOBI
      N = NSHL
      NDIM = N
      IEGEN = 0

      H = OVERLAP_AB

!     Diagonalizing OVERLAP_AB and getting eigenvalues and eigenvectors
      CALL JACOBI(H, N, NDIM, IEGEN, EIGVAL, NR, IQ, EIGVEC, XXX, .FALSE.)

!     ******************* DEBUG print! *******************
      IF (DEBUG_INV) THEN

          WRITE(6,*) '***************************************************************'
          WRITE(6,*) '**** The Linear Equation System Resolution (JACOBI Method) ****'
          WRITE(6,*) '  **** WITHOUT CONSTRAIN CONDITION-> S_rho = N_electrons ****'
          WRITE(6,*)
          WRITE(6,*) 'DIAGONALIZATION:'
          WRITE(6,*) '   Eigenvectors:'
          WRITE(6,*)
          WRITE(6,"(3F18.10)") EIGVEC(1,1), EIGVEC(1,2), EIGVEC(1,3)
          WRITE(6,"(3F18.10)") EIGVEC(2,1), EIGVEC(2,2), EIGVEC(2,3)
          WRITE(6,"(3F18.10)") EIGVEC(3,1), EIGVEC(3,2), EIGVEC(3,3)
          WRITE(6,*)
          WRITE(6,*) '   Eigenvvalues:'
          WRITE(6,*)
          WRITE(6,"(3F18.10)") EIGVAL(1), EIGVAL(2), EIGVAL(3)
          WRITE(6,*)

      ENDIF
!     **************** END-OF DEBUG print! ***************

      !
      ! Solving Linear equations
      !

      ! Transposing EIGVEC to EIGVEC_T
      EIGVEC_T = TRANSPOSE(EIGVEC)

      ! Transforming VEC_b
      B_trf(:,1) = VEC_b(:)
      B_trf = MATMUL(EIGVEC_T, B_trf)           ! b_transf = U_transf * b

      ! Solving diagonal equations
      do I = 1, N                               !
          X_trf(I,1) = B_trf(I,1)/EIGVAL(I)     ! X_transf = D_inv * B_transf
      end do                                    !

      ! Back transformation: Obtaining the ENB_COEF values
      EMB_COEF_NC = 0.0D0
      EMB_COEF_NC = MATMUL(EIGVEC, X_trf)          ! X = U * X_transf


!     ******************* DEBUG print! *******************
      IF (DEBUG_INV) THEN

         WRITE(6,*)
         WRITE(6,*) '         EMB_COEF            X                IQ'
         WRITE(6,*)
         DO I=1, N
            WRITE(6,'(I5, F16.10, F15.10, I10)') I, EMB_COEF_NC(I,1), XXX(I), IQ(I)
         ENDDO
         WRITE(6,*)
         WRITE(6,*) '**** END-OF Linear Equation System Resolution (JACOBI Method) ****'
         WRITE(6,*) '******************************************************************'
      ENDIF
!     **************** END-OF DEBUG print! ***************


      ! Deallocating local variables
      DEALLOCATE(H, EIGVAL, EIGVEC, B_trf, EIGVEC_T, X_trf, XXX, IQ)

      EMB_COEF = EMB_COEF_NC

!
!     *** End of SUBROUTINE RESLIJAC ***
!
      END SUBROUTINE RESLIJAC_NC
