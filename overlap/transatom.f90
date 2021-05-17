     SUBROUTINE TRANSATOM(XYZ_ATOM, BASEVEC)
!
!     Purpose: Do a translation in an atom
!
!     History: - Creation (01.06.12, MJR)
!
!     ******************************************************************
!     List of input/output variables:
!     BASEVEC            : INPUT  -> Base translation vector
!     XYZ_ATOM           : INPUT  -> Original coordinates
!                          OUTPUT -> Translated coordinates
!
!     List of USED variables in EMBPARAM:
!     LATVEC(3,3)        : The lattice vectors in a filed (REAL)
!
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     List of local variables:
      INTEGER       :: I
      REAL (8)      :: XYZ_ATOM(3), BASEVEC(3), TEMP_VEC(3)


      DO I=1,3

          TEMP_VEC(I) = BASEVEC(1)*LATVEC(1,I) + &
          &             BASEVEC(2)*LATVEC(2,I) + &
          &             BASEVEC(3)*LATVEC(3,I)

      ENDDO

      DO I=1,3

          XYZ_ATOM(I) = XYZ_ATOM(I) + TEMP_VEC(I)

      ENDDO

!
!     *** End of SUBROUTINE TRANSATOM ***
!
      END SUBROUTINE TRANSATOM
