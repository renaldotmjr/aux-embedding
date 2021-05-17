     SUBROUTINE READINPUT(F_NAME)
!
!     Purpose: Read the input file and set the variables to the program
!              runing
!
!     History: - Creation (06.06.12, MJR)
!
!     ******************************************************************
!     List of INPUT variables:
!     F_NAME              : The name of the input file
!
!     List of local variables:
!
!
!     List of USED variables from EMBPARAM:
!
!     ******************************************************************

      USE INPDEFVARS
      USE EMBTOL

      IMPLICIT NONE

      INTEGER               :: i

      CHARACTER (80)        :: F_NAME
      CHARACTER (200)       :: TMP_STR, SUB_TMP_STR_1
      CHARACTER (200)       :: SUB_A1, SUB_A2, SUB_A3, SUB_A4, SUB_A5
      CHARACTER             :: SUB_AA(200,3), AAA

      INTEGER               :: L_F_NAME, LT_TMP_STR, LT_TMP
      INTEGER               :: SUB_D1, SUB_D2, SUB_D3

      REAL(8)               :: SUB_F1


      ! Opening the PWDE input
      F_NAME = ADJUSTL(F_NAME)
      F_NAME = TRIM(F_NAME)
      L_F_NAME = LEN_TRIM(F_NAME)
      WRITE(6,*) 'F_NAME file:   ', '|',F_NAME(1:L_F_NAME),'|'
      WRITE(6,*)
      OPEN(99,FILE=F_NAME(1:L_F_NAME), STATUS='OLD')

   10 READ(99,"(A)",END=20) TMP_STR

      ! Adjusting to the left and calculating the lenght of TMP_STR
      TMP_STR = ADJUSTL(TMP_STR)
      LT_TMP_STR = LEN_TRIM(TMP_STR)

      ! Initializing SUB_A1
      SUB_A1 = ' '

      ! Ignore blank lines!
      IF (LT_TMP_STR .GT. 0) THEN
          read(TMP_STR, *) SUB_A1
          WRITE(6,*) TMP_STR(1:LT_TMP_STR)
      ENDIF


!     ##################################################################
!     Starting the search for the keywords


      !Searching PW_OUT name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'PW_OUT'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          PW_OUT_NAME = SUB_A2(1:LT_TMP)

      ENDIF
      !---------------------------------------------------------

      !Searching AUXIS name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'AUXIS'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2, SUB_A3
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          AUXIS_NAME = SUB_A2(1:LT_TMP)

          ! getting the type of AUXIS file!
          LT_TMP = LEN_TRIM(SUB_A3)
          IF (LT_TMP .GT. 0) THEN
              AUXIS_STD = SUB_A3(1:LT_TMP)
          ENDIF

      ENDIF
      !---------------------------------------------------------

      !Searching INTEGRATION name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'INTEGRATION'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          read(SUB_A2, "(F15.9)") INTTOL

      ENDIF
      !---------------------------------------------------------

      !Searching UNIT name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'UNIT'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          UNITSY = SUB_A2(1:LT_TMP)

      ENDIF
      !---------------------------------------------------------

      !Searching FITTESTS name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'FITTESTS'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          FITTEST = SUB_A2(1:LT_TMP)

      ENDIF
      !---------------------------------------------------------

      !Searching EQUIVALENCE name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'EQUIVALENCE'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          EQUIVALENCE = SUB_A2(1:LT_TMP)

      ENDIF
      !---------------------------------------------------------

      !Searching GENDEMONCUBE name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'GENDEMONCUBE'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2, SUB_F1
          SUB_A2 = TRIM(SUB_A2)
          LT_TMP = LEN_TRIM(SUB_A2)
          deMon_CUB_NAME = SUB_A2(1:LT_TMP)

          F_shr = SUB_F1

          PRT_deMon_CUBE = 'YES'

      ENDIF
      !---------------------------------------------------------

      !Searching REPLICATION name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'REPLICATION'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          do i = 1, 3
              READ(99,*) SUB_D1, SUB_D2
              REPNUM(i,1) = SUB_D1
              REPNUM(i,2) = SUB_D2
          enddo

          REPLICATION = 'YES'

      ENDIF
      !---------------------------------------------------------

      !Searching TAKEUCOUT name at input file!
      !---------------------------------------------------------
      SUB_TMP_STR_1 = 'TAKEUCOUT'
      SUB_A1 = TRIM(SUB_A1)

      IF (SUB_A1 .EQ. SUB_TMP_STR_1) THEN

          read(TMP_STR, *) SUB_A1, SUB_A2
          SUB_A2 = TRIM(SUB_A2)
          read(SUB_A2, *) NTAKEUCOUT

          ALLOCATE( TAKEUCOUTvecs(NTAKEUCOUT, 3) )

          do i = 1, NTAKEUCOUT
              READ(99,*) SUB_D1, SUB_D2, SUB_D3
              TAKEUCOUTvecs(i,1) = SUB_D1
              TAKEUCOUTvecs(i,2) = SUB_D2
              TAKEUCOUTvecs(i,3) = SUB_D3
          enddo

          BTAKEUCOUT = 'YES'

      ENDIF
      !---------------------------------------------------------


      GOTO 10

   20 CONTINUE


      ! Closing the PWDE file
      CLOSE(99)

      WRITE(6,*)
      WRITE(6,*) 'End of ', F_NAME(1:L_F_NAME), ' reading!'
      WRITE(6,*)

!
!     *** End of SUBROUTINE REDINPUT ***
!
      END SUBROUTINE READINPUT
