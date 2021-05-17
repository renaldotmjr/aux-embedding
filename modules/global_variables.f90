module global_variables

  implicit none

  public COVRAD, STDMATOM, ELSYM, DEBUG 
  REAL :: COVRAD(0:103),STDMATOM(0:103),PBCBOX(3,3),PBCINVBOX(3,3)
  CHARACTER(4) :: ELSYM(0:103)
  LOGICAL :: DEBUG, DEBUG_INV

end module global_variables
