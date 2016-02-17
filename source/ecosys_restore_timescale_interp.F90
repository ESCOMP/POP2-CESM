module ecosys_restore_timescale_interp
  !
  ! Control spatial variability of restoring timescale with
  ! interpolation based on depth
  !

  use marbl_kinds_mod, only : r8, log_kind, int_kind
  use marbl_constants_mod, only : c0, c2, c1000
  use domain_size, only : km ! FIXME

  implicit none

  private

  type, public :: ecosys_restore_timescale_interp_type
     real (r8), dimension(km), public :: &
          inv_restoring_time_scale ! inverse restoring time scale for nutrients (1/secs)

     real (r8), private :: rest_time_inv_surf ! inverse restoring timescale at surface
     real (r8), private :: rest_time_inv_deep ! inverse restoring timescale at depth
     real (r8), private :: rest_z0 ! shallow end of transition regime
     real (r8), private :: rest_z1 ! deep end of transition regime

   contains
     procedure, public :: init
     procedure, private :: read_namelist
     procedure, private :: interpolate_restoring_timescale

  end type ecosys_restore_timescale_interp_type

  real (r8), parameter, private :: default_rest_time_inv_surf = c0
  real (r8), parameter, private :: default_rest_time_inv_deep = c0
  real (r8), parameter, private :: default_rest_z0 = c1000
  real (r8), parameter, private :: default_rest_z1 = c2 * c1000

contains


!*****************************************************************************

subroutine init(this, nl_buffer, zt, marbl_status_log)

  use marbl_kinds_mod, only : i4, r8
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : error_msg

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_interp_type) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(marbl_log_type), intent(inout) :: marbl_status_log
  real (r8), dimension(km), intent(in) :: zt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(*), parameter :: subname = 'ecosys_restore_timescale_interp:init'

  !-----------------------------------------------------------------------
  call this%read_namelist(nl_buffer, marbl_status_log)
  if (marbl_status_log%labort_marbl) then
    error_msg = "error code returned from this%read_namelist"
    call marbl_status_log%log_error(error_msg, subname)
    return
  end if
  call this%interpolate_restoring_timescale(zt)
  
end subroutine init

!*****************************************************************************

subroutine read_namelist(this, nl_buffer, marbl_status_log)

  use marbl_kinds_mod, only : r8, i4
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use marbl_namelist_mod    , only : marbl_namelist
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : error_msg

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_interp_type) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: nml_error
  real(r8) :: rest_time_inv_surf, rest_time_inv_deep, rest_z0, rest_z1
  character(*), parameter :: subname = 'ecosys_restore_timescale_interp:read_namelist'
  character(len=marbl_nl_buffer_size) :: tmp_nl_buffer

  !-----------------------------------------------------------------------
  namelist /ecosys_restore_timescale_interp_nml/ &
       rest_time_inv_surf, &
       rest_time_inv_deep, &
       rest_z0, rest_z1

  rest_time_inv_surf = default_rest_time_inv_surf
  rest_time_inv_deep = default_rest_time_inv_deep
  rest_z0 = default_rest_z0
  rest_z1 = default_rest_z1

  tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_restore_timescale_interp_nml')
  read(tmp_nl_buffer, nml=ecosys_restore_timescale_interp_nml, iostat=nml_error)
  if (nml_error /= 0) then
     error_msg = "Error reading ecosys_restore_timescale_interp_nml"
     call marbl_status_log%log_error(error_msg, subname)
     return
  else
    call marbl_status_log%log_namelist('ecosys_restore_timescale_interp_nml', &
                                       tmp_nl_buffer, subname)
  end if

  this%rest_time_inv_surf = rest_time_inv_surf
  this%rest_time_inv_deep = rest_time_inv_deep
  this%rest_z0 = rest_z0
  this%rest_z1 = rest_z1

end subroutine read_namelist

!*****************************************************************************

subroutine interpolate_restoring_timescale(this, zt)
  !
  ! Initialize the spatial variability of the restoring time scale
  ! with an interpolation.
  !
  use marbl_kinds_mod, only : int_kind, r8
  use marbl_constants_mod, only : p5
  use domain_size, only : km ! FIXME

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_interp_type) :: this
  real (r8), dimension(km) :: zt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: k

  do k = 1, km
     if (zt(k) < this%rest_z0) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf
     else if (zt(k) > this%rest_z1) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_deep
     else if (this%rest_z1 == this%rest_z0) then
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf + p5 * &
             (this%rest_time_inv_deep - this%rest_time_inv_surf)
     else
        this%inv_restoring_time_scale(k) = this%rest_time_inv_surf + &
             (zt(k) - this%rest_z0) / (this%rest_z1 - this%rest_z0) * &
             (this%rest_time_inv_deep - this%rest_time_inv_surf)
     endif
  end do
end subroutine interpolate_restoring_timescale

!*****************************************************************************

end module ecosys_restore_timescale_interp
