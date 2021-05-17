      REAL FUNCTION PYTHAG(A,B)
!C
!C     Purpose: PYTHAG finds SQRT(A**2 + B**2) without overflow or
!C              destructive underflow. This function is used by the
!C              jacobi diagonalization routine.
!C
!C     History: - Implemented (24.02.95, MK)
!C
!C     ******************************************************************
!C
      IMPLICIT NONE
!C
      REAL A,B,P,R,S,T,U
!C
!C     ------------------------------------------------------------------
!C
!C     *** Calculate (A**2 + B**2)**1/2 by binomial expansion ***
!C
      P = MAX(ABS(A),ABS(B))
      IF (P.EQ.0.0) GO TO 20
      R = (MIN(ABS(A),ABS(B))/P)**2
!C
   10 CONTINUE
!C
      T = 4.0 + R
      IF (T.EQ.4.0) GO TO 20
      S = R/T
      U = 1.0 + 2.0*S
      P = U*P
      R = (S/U)**2*R
      GO TO 10
!C
   20 CONTINUE
!C
      PYTHAG = P
!C
!C     ------------------------------------------------------------------
!C
!C     *** End of REAL FUNCTION PYTHAG ***
!C
      END
