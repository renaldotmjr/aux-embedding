     SUBROUTINE DRVOVERLAP(IL,JL, ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP, IO,JO)
!
!     Purpose: Driver routine to calculate OVERLAP INTegral
!
!     History: - Creation (31.05.12, MJR)
!
!     ******************************************************************
!
!     List of local variables:
!     IL, JL           : Angular momentum of the shells
!     OL, OL           : Number of functions per shell in atoms I and J
!     ZETA, ZETB       : Exponents of each shell
!     XA,YA,ZA         : Atomic coordinates of center A.
!     XB,YB,ZB         : Atomic coordinates of center B.
!
!
!     ******************************************************************
!
      USE EMBPARAM

      IMPLICIT NONE

!     List of Local variables
      INTEGER       :: IL, JL, IO, JO
      REAL (8)      :: ZETA, ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP(IO,JO)


      ! Calling the subroutine for integrations

      !******************** Call [s_A|s_B] block ********************
      IF ( (IL == 0) .AND. (JL == 0) ) THEN
         CALL OVERLAP_SS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [s_A|p_B] block ********************
      ELSEIF ( (IL == 0) .AND. (JL == 1) ) THEN
         CALL OVERLAP_SP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [s_A|d_B] block ********************
      ELSEIF ( (IL == 0) .AND. (JL == 2) ) THEN
         CALL OVERLAP_SD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [p_A|s_B] block ********************
      ELSEIF ( (IL == 1) .AND. (JL == 0) ) THEN
         CALL OVERLAP_PS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [p_A|p_B] block ********************
      ELSEIF ( (IL == 1) .AND. (JL == 1) ) THEN
         CALL OVERLAP_PP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [p_A|d_B] block ********************
      ELSEIF ( (IL == 1) .AND. (JL == 2) ) THEN
         CALL OVERLAP_PD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [d_A|s_B] block ********************
      ELSEIF ( (IL == 2) .AND. (JL == 0) ) THEN
         CALL OVERLAP_DS(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [d_A|p_B] block ********************
      ELSEIF ( (IL == 2) .AND. (JL == 1) ) THEN
         CALL OVERLAP_DP(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************

      !******************** Call [d_A|d_B] block ********************
      ELSEIF ( (IL == 2) .AND. (JL == 2) ) THEN
         CALL OVERLAP_DD(ZETA,ZETB, XA,YA,ZA, XB,YB,ZB, LOC_OVERLAP)
      !**************************************************************
      ENDIF

!
!     *** End of SUBROUTINE DRVOVERLAP ***
!
      END SUBROUTINE DRVOVERLAP
