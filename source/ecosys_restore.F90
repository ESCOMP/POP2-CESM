module ecosys_restore_mod
  !
  ! Module to generalize restoring any non-autotroph tracer
  !

  use marbl_kinds_mod      , only : r8, log_kind, int_kind
  use marbl_constants_mod  , only : c0, c2, c1000
  use marbl_interface_types, only : marbl_single_restoring_field_type
  use marbl_interface_types, only : marbl_domain_type
  use marbl_sizes          , only : ecosys_tracer_cnt

  implicit none

  private

  type, public :: marbl_restore_type
    logical :: lrestore_any
    type(marbl_single_restoring_field_type), dimension(ecosys_tracer_cnt) :: tracer_restore
  contains
     procedure, public :: init
!     procedure, public :: restore_tracers
!     procedure, private :: set_tracer_read_metadata
  end type marbl_restore_type

  real (r8), parameter, private :: default_rest_time_inv_surf = c0
  real (r8), parameter, private :: default_rest_time_inv_deep = c0
  real (r8), parameter, private :: default_rest_z0 = c1000
  real (r8), parameter, private :: default_rest_z1 = c2 * c1000

contains

!*****************************************************************************

subroutine init(this, nl_buffer, domain, marbl_status_log)

  ! initialize marbl_restore instance to default values, then read
  ! namelist and setup tracers that need to be restored

  use marbl_kinds_mod   , only : char_len, int_kind, i4, log_kind
  use marbl_logging     , only : marbl_log_type
  use marbl_logging     , only : error_msg
  use marbl_logging     , only : status_msg
  use marbl_namelist_mod, only : marbl_nl_cnt
  use marbl_namelist_mod, only : marbl_nl_buffer_size
  use marbl_namelist_mod, only : marbl_namelist

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(marbl_restore_type), intent(inout) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(marbl_domain_type), intent(in) :: domain

  type(marbl_log_type), intent(inout) :: marbl_status_log
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_file_varnames
  real(r8) :: rest_time_inv_surf, rest_time_inv_deep, rest_z0, rest_z1
  integer(int_kind) :: nml_error, t
  character(len=marbl_nl_buffer_size) :: tmp_nl_buffer
  character(*), parameter :: subname = 'ecosys_restore:Init'

  !-----------------------------------------------------------------------

  namelist /ecosys_restore_nml/       &
       restore_short_names,           &
       restore_filenames,             &
       restore_file_varnames,         &
       rest_time_inv_surf,            &
       rest_time_inv_deep,            &
       rest_z0, rest_z1

  this%lrestore_any = .false.

!  do t = 1, ecosys_tracer_cnt
!     call this%set_tracer_read_metadata(t)
!  end do

  ! initialize namelist variables to default values
  restore_short_names = ''
  restore_filenames = ''
  restore_file_varnames = ''

  rest_time_inv_surf = default_rest_time_inv_surf
  rest_time_inv_deep = default_rest_time_inv_deep
  rest_z0 = default_rest_z0
  rest_z1 = default_rest_z1

  tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_restore_nml')
  call marbl_status_log%log_noerror(tmp_nl_buffer, "MNL MNL MNL")
  read(tmp_nl_buffer, nml=ecosys_restore_nml, iostat=nml_error)
  if (nml_error /= 0) then
     error_msg = "Error reading ecosys_restore_nml"
     call marbl_status_log%log_error(error_msg, subname)
     return
  else
    call marbl_status_log%log_namelist('ecosys_restore_nml', tmp_nl_buffer, subname)
  end if

  ! FIXME(bja, 2014-10) assert(len(restore_short_names) == len(restore_filenames))
  write(status_msg, "(A)") "Found restore variables : "
  call marbl_status_log%log_noerror(status_msg, subname)
  do t = 1, size(restore_short_names)
     if (len(trim(restore_short_names(t))) > 0) then
        write(status_msg, "(6A)") trim(restore_short_names(t)), " --> ", &
             trim(restore_filenames(t)), " [ ", trim(restore_file_varnames(t)), " ]"
        call marbl_status_log%log_noerror(status_msg, subname)
     end if
     ! MNL TO DO: allocate appropriate memory in this%tracer_restore
  end do

end subroutine Init

#if 0
!*****************************************************************************

subroutine restore_tracers(this)
  !
  !  restore a variable if required
  !
  use marbl_kinds_mod       , only : r8, int_kind, log_kind
  use marbl_constants_mod   , only : c0

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(marbl_restore_type), intent(inout) :: this

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: n
  !-----------------------------------------------------------------------

!    do n=1,tracer_cnt
!    enddo
  end associate

end subroutine restore_tracers

!*****************************************************************************

subroutine set_tracer_read_metadata(this, index, &
     mod_varname, filename, file_varname, file_fmt, &
     scale_factor, default_val)

  use marbl_kinds_mod    , only : char_len, r8
  use marbl_constants_mod, only : c0, c1

  implicit none

  !  initialize a tracer_read type to common default values.

  class(marbl_restore_type), intent(inout) :: this
  integer(int_kind)              , intent(in) :: index
  character(char_len) , optional , intent(in) :: mod_varname
  character(char_len) , optional , intent(in) :: filename
  character(char_len) , optional , intent(in) :: file_varname
  character(char_len) , optional , intent(in) :: file_fmt
  real(r8)            , optional , intent(in) :: scale_factor
  real(r8)            , optional , intent(in) :: default_val


  this%tracers(index)%restore_file_info%mod_varname  = 'unknown'
  this%tracers(index)%restore_file_info%filename     = 'unknown'
  this%tracers(index)%restore_file_info%file_varname = 'unknown'
  this%tracers(index)%restore_file_info%file_fmt     = 'bin'
  this%tracers(index)%restore_file_info%scale_factor = c1
  this%tracers(index)%restore_file_info%default_val  = c0

  if (present(mod_varname  )) this%tracers(index)%restore_file_info%mod_varname  = mod_varname
  if (present(filename     )) this%tracers(index)%restore_file_info%filename     = filename
  if (present(file_varname )) this%tracers(index)%restore_file_info%file_varname = file_varname
  if (present(file_fmt     )) this%tracers(index)%restore_file_info%file_fmt     = file_fmt
  if (present(scale_factor )) this%tracers(index)%restore_file_info%scale_factor = scale_factor
  if (present(default_val  )) this%tracers(index)%restore_file_info%default_val  = default_val

end subroutine set_tracer_read_metadata
#endif

!*****************************************************************************

end module ecosys_restore_mod
