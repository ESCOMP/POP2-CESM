module ecosys_restore_timescale_file
  !
  ! Control spatial variability of restoring timescale with an input file
  !
  use marbl_kinds_mod, only : r8, int_kind

  implicit none

  private

  type, public :: ecosys_restore_timescale_file_type
     ! inverse restoring timescale for variable interior restoring
     real (r8), dimension(:,:,:), allocatable, public :: restore_rtau
     ! maximum level for applying variable interior restoring
     integer (int_kind), dimension(:,:,:), allocatable, public :: restore_max_level

   contains
     procedure, public :: init
     procedure, private :: read_namelist
     procedure, private :: read_restoring_timescale_from_file

  end type ecosys_restore_timescale_file_type

contains

!*****************************************************************************

subroutine init(this, nl_buffer, marbl_status_log)

  use marbl_kinds_mod, only : char_len, i4
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : error_msg

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(char_len) :: file_format
  character(char_len) :: file_name
  character(*), parameter :: subname = 'ecosys_restore_timescale_file:init'

  !-----------------------------------------------------------------------
  call this%read_namelist(nl_buffer, file_name, file_format, marbl_status_log)
  if (marbl_status_log%labort_marbl) then
    error_msg = "error code returned from this%read_namelist"
    call marbl_status_log%log_error(error_msg, subname)
    return
  end if
  call this%read_restoring_timescale_from_file(file_name, file_format)
  
end subroutine init

!*****************************************************************************

subroutine read_namelist(this, nl_buffer, file_name, file_format, marbl_status_log)

  use marbl_kinds_mod, only : int_kind, i4, char_len
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use marbl_namelist_mod    , only : marbl_namelist
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : error_msg

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(char_len), intent(out) :: file_format
  character(char_len), intent(out) :: file_name
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: nml_error
  character(char_len) :: restore_timescale_file_format
  character(char_len) :: restore_timescale_file_name
  character(*), parameter :: subname = 'ecosys_restore_timescale_file:read_namelist'
  character(len=marbl_nl_buffer_size) :: tmp_nl_buffer

  !-----------------------------------------------------------------------
  namelist /ecosys_restore_timescale_file_nml/ &
       restore_timescale_file_format, &
       restore_timescale_file_name

  restore_timescale_file_format = ''
  restore_timescale_file_name = ''

  tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_restore_timescale_file_nml')
  read(tmp_nl_buffer, nml=ecosys_restore_timescale_file_nml, iostat=nml_error)
  if (nml_error /= 0) then
     error_msg = "Error reading ecosys_restore_timescale_file_nml"
     call marbl_status_log%log_error(error_msg, subname)
     return
  else
    call marbl_status_log%log_namelist('ecosys_restore_timescale_file_nml', &
                                       tmp_nl_buffer, subname)
  end if

  file_format = restore_timescale_file_format
  file_name = restore_timescale_file_name

end subroutine read_namelist

!*****************************************************************************

subroutine read_restoring_timescale_from_file(this, file_format, file_name)
  !
  ! Initialize the spatially variable restoring timescale from the the
  ! user specified file
  !
  use marbl_kinds_mod, only : char_len
  use blocks, only : nx_block, ny_block ! FIXME
  use prognostic, only : max_blocks_clinic ! FIXME
  use time_management, only : seconds_in_day ! FIXME
  use passive_tracer_tools, only : read_field ! FIXME

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_timescale_file_type) :: this
  character(len=char_len), intent(in) :: file_format
  character(len=char_len), intent(in) :: file_name

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

    allocate(this%restore_rtau(nx_block, ny_block, max_blocks_clinic))
    allocate(this%restore_max_level(nx_block, ny_block, max_blocks_clinic))

    call read_field(file_format, file_name, &
         'NUTR_RESTORE_MAX_LEVEL', this%restore_rtau)

    this%restore_max_level = nint(this%restore_rtau)

    call read_field(file_format, file_name, &
         'NUTR_RESTORE_RTAU', this%restore_rtau)

    this%restore_rtau = this%restore_rtau / seconds_in_day ! convert days to secs

  end subroutine read_restoring_timescale_from_file

!*****************************************************************************

end module
