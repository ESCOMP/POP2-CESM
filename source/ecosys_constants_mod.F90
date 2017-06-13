module ecosys_constants_mod

  ! !DESCRIPTION:
  !  This module provides a place to set parameters (in the Fortran sense) that
  !  are used across the various ecosystem subroutines
  !
  !  Written by: Michael Levy, NCAR, June 2017

  use kinds_mod, only : int_kind

  Implicit None
  public

  ! # of tracers expected from MARBL
  integer(kind=int_kind), parameter :: marbl_tracer_cnt = MARBL_NT

  ! Parameters handling size of namelists
  integer(kind=int_kind), parameter :: ecosys_nl_in_size     = 262144
  integer(kind=int_kind), parameter :: ecosys_nl_cnt         = 256
  integer(kind=int_kind), parameter :: ecosys_nl_buffer_size = 32768

  ! Surface forcing output variables
  integer(kind=int_kind), parameter :: sfo_cnt = 2
  integer(kind=int_kind)            :: flux_co2_id
  integer(kind=int_kind)            :: totalChl_id

end module ecosys_constants_mod
