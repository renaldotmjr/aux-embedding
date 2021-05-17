MODULE String_Manipulation
!
!     Purpose: Some routines for string manipulation.
!
!     History: - Creation (02.02.2012, MJR)
!
!     Functions based on the book: Upgrading to Fortran 90. By Cooper
!                                  Redwine (1995). pp.80
!
!     ******************************************************************
!
!
!
!     ******************************************************************

   IMPLICIT NONE

   ! The last character is the blank (space) that is not changed by the routins
   CHARACTER(*), PRIVATE, PARAMETER   :: LOWER = 'abcdefghijklmnopqrstuvwxyz '
   CHARACTER(*), PRIVATE, PARAMETER   :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
   CHARACTER(*), PRIVATE, PARAMETER   :: NUMB  = '1234567890 '

!  The content of the module
   CONTAINS

!
!  SUBROUTINE TO CONVERT LOWER TO UPPER CASE
!
   SUBROUTINE UpCase(Input)

!    Input argument and the output result
     CHARACTER(*)                :: Input
     CHARACTER( LEN(Input) )     :: Output

!    Local variables
     INTEGER    :: i, n

     ! Initializating the output arument with the same content of input
     Output = Input

     ! Loop over output string elements
     DO i = 1, LEN(Output)

       ! Find the location of an output letter in LOWER parameter string
       n = INDEX( LOWER, Output(i:i) )

       ! If current substring is a lower case letter, make it upper case
       IF (n /= 0) Output(i:i) = UPPER(n:n)

     END DO

     Input = Output

   END SUBROUTINE UpCase

!
!  SUBROUTINE TO CONVERT UPPER TO LOWER CASE
!
   SUBROUTINE LowCase(Input)

!    Input argument and the output result
     CHARACTER(*)                :: Input
     CHARACTER( LEN(Input) )     :: Output

!    Local variables
     INTEGER    :: i, n

     ! Initializating the output arument with the same content of input
     Output = Input

     ! Loop over OUTPUT string elements
     DO i = 1, LEN(Output)

       ! Find the location of an output letter in UPPER parameter string
       n = INDEX( UPPER, Output(i:i) )

       ! If current substring is an upper case letter, make it lower case
       IF ( n /= 0 ) Output(i:i) = LOWER(n:n)

     END DO

     Input = Output

   END SUBROUTINE LowCase

!
!  SUBROUTINE TO CONVERT NUMBER CARACTERS IN BLANK SPACE
!
   SUBROUTINE NumSpace(Input)

!    Input argument and the output result
     CHARACTER(*)                :: Input
     CHARACTER( LEN(Input) )     :: Output

!    Local variables
     INTEGER    :: i, n

     ! Initializating the output arument with the same content of input
     Output = Input

     ! Loop over OUTPUT string elements
     DO i = 1, LEN(Output)

       ! Find the location of an output number letter in NUMB parameter string
       n = INDEX( NUMB, Output(i:i) )

       ! If current substring is a number letter, make it a space
       IF ( n /= 0 ) Output(i:i) = ' '

     END DO

     Input = Output

   END SUBROUTINE NumSpace

  END MODULE String_Manipulation
