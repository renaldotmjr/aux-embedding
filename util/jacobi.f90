      SUBROUTINE JACOBI(H,N,NDIM,IEGEN,EIGVAL,NR,IQ,U,X, SORT)
!C
!C     Purpose: JACOBI computes eigenvalues and eigenvectors of a real
!C              symmetric matrix using the Jacobi method.
!C
!C     History: - Creation (June 1994, B. Ahlswede; the original program
!C                          was written by F. Corbato and M. Merwin)
!C
!C              - Implemented in deMon (19.10.00, AMK)
!
!               - Inclusion of the SORTING option (02.02.2012, MJR, CO)
!                  -> SORT = .TRUE. : Sort the eigenvalues
!                  -> SORT = .FALSE. : Do not sort the eigenvalues
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
!C
!C     ------------------------------------------------------------------
!C
      IMPLICIT NONE

      LOGICAL :: SORT
!C
      REAL EPS
      PARAMETER (EPS = 1.E-10)
!C
      INTEGER I,IEGEN,IMI1,IPIV,J,JPIV,K,N,NDIM,NR
      REAL COSINE,HDIMIN,HDTEST,HII,HTEMP,P,PYTHAG,SINE,TANG,XDIF,XMAX
!C
      REAL EIGVAL(NDIM),H(NDIM,NDIM)
!C
      INTEGER IQ(NDIM)
      REAL U(NDIM,NDIM),X(NDIM)
!C
!C     ------------------------------------------------------------------
!C
!C     *** Save the diagonal elements ***
!C
      DO 10 I=1,N
        EIGVAL(I) = H(I,I)
   10 CONTINUE
!C
!C     *** Initialize eigenvectors ***
!C
      IF (IEGEN.EQ.0) THEN
        DO 610 I=1,N
          DO 600 J=1,N
            IF (I.EQ.J) THEN
              U(J,I) = 1.0
            ELSE
              U(J,I) = 0.0
            END IF
  600     CONTINUE
  610   CONTINUE
      END IF
!C
      NR = 0
      IF (N.LE.1) GO TO 1100
!C
!C     *** Scan for largest off diagonal element in each column      ***
!C     *** X(I) contains the largest element in ith column           ***
!C     *** IQ(I) holds second subscript defining position of element ***
!C
      DO 35 I=2,N
        X(I) = 0.0
        IMI1 = I - 1
        DO 30 J=1,IMI1
          IF (X(I).LE.ABS(H(J,I))) THEN
            X(I) = ABS(H(J,I))
            IQ(I) = J
          END IF
   30   CONTINUE
   35 CONTINUE
!C
!C     *** Set indicator for shut-off ***
!C
      HDTEST = 1.0E38
!C
!C     *** Find maximum of X(I)'s for pivot element ***
!C     *** and test for end of problem              ***
!C
   40 CONTINUE
!C
      XMAX = 0.0
      DO 70 I=2,N
        IF (XMAX.LT.X(I)) THEN
          XMAX = X(I)
          IPIV = IQ(I)
          JPIV = I
        END IF
   70 CONTINUE
!C
      IF (XMAX.LE.0.0) GO TO 1000
      IF ((HDTEST.GT.0.0).AND.(XMAX.GT.HDTEST)) GO TO 148
      HDIMIN = ABS(EIGVAL(1))
!C
!C     *** Search for smallest diagonalelement ***
!C
      DO 110 I=2,N
        IF (HDIMIN.GT.ABS(EIGVAL(I))) THEN
          HDIMIN = ABS(EIGVAL(I))
        END IF
  110 CONTINUE
!C
!C     *** End rotations if MAX(H(I,J)).LT.(EPS)*ABS(H(K,K)-MIN) ***
!C
      HDTEST = HDIMIN*EPS
      IF (HDTEST.GE.XMAX) GO TO 1000
!C
  148 CONTINUE
!C
      NR = NR + 1
!C
!C     *** Compute TANG, SINE AND COSINE, H(I,I), H(J,J) ***
!C     *** Formel: Quantenchemie Bd. 5,  S. 473          ***
!C
      XDIF = EIGVAL(IPIV) - EIGVAL(JPIV)
      P = XDIF/(2.0*H(IPIV,JPIV))
!C
!C     *** PYTHAG berechnet SQRT(A**2 + B**2) - aus Eispack ***
!C
      TANG = 1.0/(P+SIGN(1.0,P)*PYTHAG(P,1.0))
      COSINE = 1.0/PYTHAG(1.0,TANG)
      SINE = TANG*COSINE
      HII = EIGVAL(IPIV)
      EIGVAL(IPIV) = COSINE*COSINE*(HII + TANG*(2.0*H(IPIV,JPIV) + &
     &               TANG*EIGVAL(JPIV)))
      EIGVAL(JPIV) = COSINE*COSINE*(EIGVAL(JPIV) - &
     &               TANG*(2.0*H(IPIV,JPIV) - TANG*HII))
      H(IPIV,JPIV) = 0.0
!C
!C     *** Inspect the IQ's between I+1 and N-1 to determine   ***
!C     *** wether a new maximum value should be computed since ***
!C     *** the present maximum is in the I's or J's column     ***
!C
      DO 350 I=2,N
        IF ((I.EQ.IPIV).OR.(I.EQ.JPIV)) GO TO 350
        IF ((IQ(I).EQ.IPIV).OR.(IQ(I).EQ.JPIV)) THEN
          K = IQ(I)
          HTEMP = H(K,I)
          H(K,I) = 0.0
          IMI1 = I - 1
          X(I) = 0.0
!C
!C     *** Search in depleted column new maximum ***
!C
          DO 320 J=1,IMI1
            IF (X(I).GT.ABS(H(J,I))) THEN
              X(I) = ABS(H(J,I))
              IQ(I) = J
            END IF
  320     CONTINUE
          H(K,I) = HTEMP
        END IF
  350 CONTINUE
!C
!C     *** Change the order of elements of H ***
!C
!C     *** I < IPIV, I < JPIV ***
!C
      X(IPIV) = 0.0
      X(JPIV) = 0.0
      DO 400 I=1,IPIV-1
        HTEMP = H(I,IPIV)
        H(I,IPIV) = COSINE*HTEMP + SINE*H(I,JPIV)
        IF (X(IPIV).LT.ABS(H(I,IPIV))) THEN
          X(IPIV) = ABS(H(I,IPIV))
          IQ(IPIV) = I
        END IF
        H(I,JPIV) = - SINE*HTEMP + COSINE*H(I,JPIV)
        IF (X(JPIV).LT.ABS(H(I,JPIV))) THEN
          X(JPIV) = ABS(H(I,JPIV))
          IQ(JPIV) = I
        END IF
  400 CONTINUE
!C
!C     *** I > IPIV, I < JPIV ***
!C
      DO 410 I=IPIV+1,JPIV-1
        HTEMP = H(IPIV,I)
        H(IPIV,I) = COSINE*HTEMP + SINE*H(I,JPIV)
        IF (X(I).LT.ABS(H(IPIV,I))) THEN
          X(I) = ABS(H(IPIV,I))
          IQ(I) = IPIV
        END IF
        H(I,JPIV) = - SINE*HTEMP + COSINE*H(I,JPIV)
        IF (X(JPIV).LT.ABS(H(I,JPIV))) THEN
          X(JPIV) = ABS(H(I,JPIV))
          IQ(JPIV) = I
        END IF
  410 CONTINUE
!C
!C     *** I > IPIV, I > JPIV ***
!C
      DO 420 I=JPIV+1,N
        HTEMP = H(IPIV,I)
        H(IPIV,I) = COSINE*HTEMP + SINE*H(JPIV,I)
        IF (X(I).LT.ABS(H(IPIV,I))) THEN
          X(I) = ABS(H(IPIV,I))
          IQ(I) = IPIV
        END IF
        H(JPIV,I) = - SINE*HTEMP + COSINE*H(JPIV,I)
        IF (X(I).LT.ABS(H(JPIV,I))) THEN
          X(I) = ABS(H(JPIV,I))
          IQ(I) = JPIV
        END IF
  420 CONTINUE
!C
!C     *** Compute new eigenvectors ***
!C
      IF (IEGEN.EQ.0) THEN
        DO 550 I=1,N
          HTEMP = U(I,IPIV)
          U(I,IPIV) = COSINE*HTEMP + SINE*U(I,JPIV)
          U(I,JPIV) = - SINE*HTEMP + COSINE*U(I,JPIV)
  550  CONTINUE
      END IF
!C
!C     *** GO TO NEXT CYCLE ***
!C
      GO TO 40
!C
 1000 CONTINUE

!C
!C     *** Sort eigenvalues ***
!C
! Inclusion of the SORTING option (02.02.2012, MJR, CO)
      IF(SORT .EQV. .TRUE.) THEN

          DO 920 I=1,N-1
            K = I
            P = EIGVAL(I)
            DO 910 J=I+1,N
              IF (EIGVAL(J).LT.P) THEN
                K = J
                P = EIGVAL(J)
              END IF
  910       CONTINUE
            IF (K.NE.I) THEN
              EIGVAL(K) = EIGVAL(I)
              EIGVAL(I) = P
              DO 900 J=1,N
                P = U(J,I)
                U(J,I) = U(J,K)
                U(J,K) = P
  900         CONTINUE
            END IF
  920     CONTINUE

      END IF

!C
 1100 CONTINUE
!C
!C     ------------------------------------------------------------------
!C
!C     *** End of SUBROUTINE JACOBI ***
!C
      END
