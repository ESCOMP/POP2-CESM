module ecosys_constants

  use kinds_mod, only : int_kind

  implicit none

  public

!-----------------------------------------------------------------------------
! number of ecosystem tracers (also in ecosys_constants.F90)
!-----------------------------------------------------------------------------
  integer(int_kind), parameter :: ecosys_tracer_cnt = ECOSYS_NT

!-----------------------------------------------------------------------------
! number of ecosystem constituents and grazing interactions
!-----------------------------------------------------------------------------
  integer (KIND=int_kind), parameter :: &
       zooplankton_cnt = ZOOPLANKTON_CNT, &
       autotroph_cnt   = AUTOTROPH_CNT,   &
       grazer_prey_cnt = GRAZER_PREY_CNT

end module ecosys_constants
