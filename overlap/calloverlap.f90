      SUBROUTINE CALLOVERLAP(LTA, LTB, ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP)
!
!     Purpose: Call the OVERLAP_ij subroutines
!
!     History: - Creation (30.09.11, MJR)
!
!     ******************************************************************
!
!     List of input variables:
!     LCA      : Total angular momentum of function A
!     LCB      : Total angular momentum of function B
!     ZETA     : Gaussian exponent of auxiliary function A.
!     ZETB     : Gaussian exponent of auxiliary function B.
!     XA,YA,ZA : Atomic coordinates of center A.
!     XB,YB,ZB : Atomic coordinates of center B.
!
!     ******************************************************************
!
      IMPLICIT NONE

!     ------------------------------------------------------------------

!     Local variables:
      INTEGER :: I, J

!     List of input variables:
      INTEGER :: LTA, LTB
      REAL(8) :: ZETA,ZETB, XA,YA,ZA, XB,YB,ZB

!     List of output variable:
      REAL (8), DIMENSION(:,:), ALLOCATABLE ::  OVERLAP_tmp


      ! Allocating the result matrix OVERLAP
      ALLOCATE( OVERLAP_tmp(LTA, LTB) )
      
      !******************** Call [s_A|s_B] block ********************
      IF ( LTA == 1) .AND. (LTB == 1) ) THEN
         CALL OVERLAP_SS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [s_A|p_B] block ********************
      ELSEIF ( (LTA == 1) .AND. (LTB == 3) ) THEN
         CALL OVERLAP_SP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [s_A|d_B] block ********************
      ELSEIF ( (LTA == 1) .AND. (LTB == 6) ) THEN
         CALL OVERLAP_SD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [p_A|s_B] block ********************
      ELSEIF ( (LTA == 3) .AND. (LTB == 1) ) THEN
         CALL OVERLAP_PS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [p_A|p_B] block ********************
      ELSEIF ( (LTA == 3) .AND. (LTB == 3) ) THEN
         CALL OVERLAP_PP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [p_A|d_B] block ********************
      ELSEIF ( (LTA == 3) .AND. (LTB == 6) ) THEN
         CALL OVERLAP_PD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [d_A|s_B] block ********************
      ELSEIF ( (LTA == 6) .AND. (LTB == 1) ) THEN
         CALL OVERLAP_DS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [d_A|p_B] block ********************
      ELSEIF ( (LTA == 6) .AND. (LTB == 3) ) THEN
         CALL OVERLAP_DP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************

      !******************** Call [d_A|d_B] block ********************
      ELSEIF ( (LTA == 6) .AND. (LTB == 6) ) THEN
         CALL OVERLAP_DD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, OVERLAP_tmp)
      !**************************************************************
      ENDIF

!     ------------------------------------------------------------------

      DO I=1, LTA
         DO J=1, LTB




         ENDDO
      ENDDO

!     ------------------------------------------------------------------
!
!     *** End of SUBROUTINE OVERLAP_SD ***
!
      END SUBROUTINE CALLOVERLAP
