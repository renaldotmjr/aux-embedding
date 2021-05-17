      SUBROUTINE CALCNSTEP(D, dx, dy, dz, NI)
!
!     Purpose: Calculate the number of steps in distance D
!
!     History: - Creation (19.10.11, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     D                 : The distance
!     dx, dy, dz        : Vector of steps (only in one coord)
!                         -> dx,  0,  0
!                         ->  0, dy,  0
!                         ->  0,  0, dz
!
!     List of local variables:
!     X, Y, Z           : Will receive dx, dy and dz
!     MOD_R             : Modulus of (Xa, Ya, Za) in cartesian coord.
!
!     List of output variables:
!     NI                : The number of steps in D in directions x,y,z
!
!     List of USED variables from EMBPARAM:
!
!
!     ******************************************************************

      IMPLICIT NONE

      INTEGER  :: NI(3)

      REAL(8)  :: D, dx, dy, dz, Xa, Ya, Za, MOD_R


!     *** CALCULATING NI(1) = NIx
      ! Converting to cartesian
      Xa = dx
      Ya = 0.0
      Za = 0.0
      CALL CRYS2CART(Xa, Ya, Za)
      ! Calculating the modulus
      MOD_R = SQRT(Xa**2 + Ya**2 + Za**2)
      ! Calculating the number of steps
      NI(1) = ( IDNINT(D/MOD_R) + 1 )

!     *** CALCULATING NI(2) = NIy
      ! Converting to cartesian
      Xa = 0.0
      Ya = dy
      Za = 0.0
      CALL CRYS2CART(Xa, Ya, Za)
      ! Calculating the modulus
      MOD_R = SQRT(Xa**2 + Ya**2 + Za**2)
      ! Calculating the number of steps
      NI(2) = ( IDNINT(D/MOD_R) + 1 )

!     *** CALCULATING NI(3) = NIz
      ! Converting to cartesian
      Xa = 0.0
      Ya = 0.0
      Za = dz
      CALL CRYS2CART(Xa, Ya, Za)
      ! Calculating the modulus
      MOD_R = SQRT(Xa**2 + Ya**2 + Za**2)
      ! Calculating the number of steps
      NI(3) = ( IDNINT(D/MOD_R) + 1 )


!
!     *** End of SUBROUTINE CALCNSTEPS ***
!
      END SUBROUTINE CALCNSTEP
