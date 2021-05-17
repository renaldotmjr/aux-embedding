     SUBROUTINE PUTINGRID(Xa,Ya,Za)
!
!     Purpose: Put X,Y,Z crystal coordinates inside the unitary cube.
!
!     History: - Creation (20.10.11, MJR)
!
!     ******************************************************************
!     List of INPUT/OUTPUT variables:
!     X, Y, Z     : INPUT  -> Crystalline coords in or outside the
!                             unitary cube
!                   OUTPUT -> Crystalline coords INSIDE the unitary cube
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

      REAL (8)      :: Xa,Ya,Za

!     Here, the GOTO inside the IF structures is to avoid the posterior
!      unnecessary ELSEIF inspections

!     Testing if X is in outside of grid
      IF ( Xa < 0.0D0 ) THEN
         Xa = Xa - DFLOAT( INT(Xa - 1) )
         GOTO 10
      ELSEIF ( (DABS(Xa - 1.0D0) <= 1.0D-6) .OR. (Xa > 1.0D0) ) THEN
         Xa = Xa - DFLOAT( INT(Xa) )
      ENDIF

   10 CONTINUE

!     Testing if Y is in outside of grid
      IF ( Ya < 0.0D0 ) THEN
         Ya = Ya - DFLOAT( INT(Ya - 1) )
         GOTO 20
      ELSEIF ( (DABS(Ya - 1.0D0) <= 1.0D-6) .OR. (Ya > 1.0D0) ) THEN
         Ya = Ya - DFLOAT( INT(Ya) )
      ENDIF

   20 CONTINUE

!     Testing if Z is in outside of grid
      IF ( Za < 0.0D0 ) THEN
         Za = Za - DFLOAT( INT(Za - 1) )
         GOTO 30
      ELSEIF ( (DABS(Za - 1.0D0) <= 1.0D-6) .OR. (Za > 1.0D0) ) THEN
         Za = Za - DFLOAT( INT(Za) )
      ENDIF

   30 CONTINUE

!     If X >= 0.0 .OR. X < 1.0, don't need to do anything!
!     If Y >= 0.0 .OR. Y < 1.0, don't need to do anything!
!     If Z >= 0.0 .OR. Z < 1.0, don't need to do anything!
!     The point is already on the grid

!
!     *** End of SUBROUTINE PUTXGRID ***
!
      END SUBROUTINE PUTINGRID
