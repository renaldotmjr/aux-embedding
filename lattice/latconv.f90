      SUBROUTINE LATCONV
!
!     Purpose: Creates the matrices to crystaline and cartesian conversions
!
!     History: - Creation (07.10.11, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     invPIV(3),
!     TRF_INF, TRI_INF
!     WORK(9)
!
!     List of USED variables from EMBPARAM:
!     A_LATx, A_LATy, A_LATz         : Lattice parameter a
!     B_LATx, B_LATy, B_LATz         : Lattice parameter b
!     C_LATx, C_LATy, C_LATz         : Lattice parameter c
!
!     List of USED-and-MODIFIED variables in EMBPARAM:
!
!     ******************************************************************
!
      USE EMBPARAM
      USE global_variables

      IMPLICIT NONE

      INTEGER     :: I, invPIV(3), TRF_INF, TRI_INF
      REAL(8)     :: WORK(9)


!     Generating the matrix MLAT, that converts cystaline to cartesian coods.
      MLAT(1,1) = A_LATx
      MLAT(1,2) = B_LATx
      MLAT(1,3) = C_LATx
      MLAT(2,1) = A_LATy
      MLAT(2,2) = B_LATy
      MLAT(2,3) = C_LATy
      MLAT(3,1) = A_LATz
      MLAT(3,2) = B_LATz
      MLAT(3,3) = C_LATz

!     The inverse matrix converts from cartesian to crystaline coods.
!     Inverting MLAT
      MLAT_inv = MLAT
!     DGETRF computes an LU factorization of a general M-by-N matrix A
!        using partial pivoting with row interchanges.
      CALL DGETRF( 3, 3, MLAT_inv, 3, invPIV, TRF_INF )

      if (TRF_INF == 0) then

!        DGETRI computes the inverse of a matrix using the LU factorization
!           computed by DGETRF.
!        Calculating the inverse
         CALL DGETRI( 3, MLAT_inv, 3, invPIV, WORK, 9, TRI_INF )

      else
         write(6,*) 'Problems with MLAT inversion'
         stop
      endif

! Checking if the inversion was ok
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
      if (TRI_INF /= 0) then

         write(6,*) 'Promlems with MLAT inversion'
         stop

      endif

      IF (DEBUG) THEN

         write(6,*) 'MLAT:'
         do I=1, 3
            write(6,'(3F20.13)') MLAT(i,1), MLAT(i,2), MLAT(i,3)
         enddo
         write(6,*)
         write(6,*) 'MLAT:'
         do I=1, 3
            write(6,'(3F20.13)') MLAT_inv(i,1), MLAT_inv(i,2), MLAT_inv(i,3)
         enddo

      ENDIF

!
!     *** End of SUBROUTINE LATCONV ***
!
      END SUBROUTINE LATCONV

