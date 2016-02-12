module ecosys_restore_mod
  !
  ! Module to generalize restoring any non-autotroph tracer
  !

  use marbl_kinds_mod                 , only : r8, log_kind, int_kind
  use marbl_sizes                     , only : ecosys_tracer_cnt
  use marbl_interface_types           , only : tracer_read => marbl_tracer_read_type
  use ecosys_restore_timescale_file   , only : ecosys_restore_timescale_file_type
  use ecosys_restore_timescale_interp , only : ecosys_restore_timescale_interp_type

  implicit none

  private

  type, private :: restore_info
     logical(log_kind) :: restore            ! flag indicating if this tracer should be restored
     type(tracer_read) :: restore_file_info  ! info about file containing restoring field
     integer(int_kind) :: tavg_restore_index ! index for tavg output
     real (r8), dimension(:, :, :, :), allocatable :: data ! restoring field
  end type restore_info

  type, public :: ecosys_restore_type
     logical(log_kind), public :: restore_any_tracer ! true if we are restoring any tracer
     type(restore_info), allocatable, private :: tracers(:)

     ! true if geographically varying nutrient restoring is read from
     ! a file (formally lnutr_variable_restore)
     logical(log_kind), private :: spatial_variability_from_file
     type(ecosys_restore_timescale_file_type), private :: timescale_file
     type(ecosys_restore_timescale_interp_type), private :: timescale_interp

   contains

     procedure, public :: init
     procedure, public :: read_restoring_fields
     procedure, public :: restore_tracers
     procedure, public :: initialize_restoring_timescale
     procedure, private :: read_namelist
     procedure, private :: initialize_restore_read_vars
     procedure, private :: set_tracer_read_metadata

  end type ecosys_restore_type

contains

!*****************************************************************************

subroutine init(this, nl_buffer, ind_name_table, marbl_status_log)

  ! initialize ecosys_restore instance to default values, then read
  ! namelist and setup tracers that need to be restored

  use marbl_kinds_mod       , only : char_len, int_kind, i4, log_kind
  use marbl_kinds_mod       , only : c0, c2
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : error_msg
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use passive_tracer_tools  , only : ind_name_pair, name_to_ind ! FIXME
  use marbl_interface_types , only : tracer_read  => marbl_tracer_read_type

  implicit none

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  type(ind_name_pair) , dimension(:), intent(in) :: ind_name_table

  type(marbl_log_type), intent(inout) :: marbl_status_log
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer (int_kind) :: t
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt) :: restore_file_varnames
  character(*), parameter :: subname = 'ecosys_restore:Init'

  !-----------------------------------------------------------------------

  ! initialize class data to default values
  allocate(this%tracers(ecosys_tracer_cnt))

  this%restore_any_tracer = .false.
  this%tracers(:)%restore = .false.
  this%tracers(:)%tavg_restore_index = -1

  do t = 1, ecosys_tracer_cnt
     call this%set_tracer_read_metadata(t)
  end do

  call this%read_namelist(nl_buffer, restore_short_names, restore_filenames, &
                          restore_file_varnames, marbl_status_log)
  if (marbl_status_log%labort_marbl) then
    error_msg = "error code returned from this%read_namelist"
    call marbl_status_log%log_error(error_msg, subname)
    return
  end if

  call this%initialize_restore_read_vars(restore_short_names, restore_filenames, &
       restore_file_varnames, ind_name_table, marbl_status_log)
  if (marbl_status_log%labort_marbl) then
    error_msg = "error code returned from this%initialize_restore_read_vars"
    call marbl_status_log%log_error(error_msg, subname)
    return
  end if

end subroutine Init


!*****************************************************************************
subroutine read_namelist(this, nl_buffer, restore_short_names,                &
                         restore_filenames, restore_file_varnames,            &
                         marbl_status_log)

  ! Read the ecosys_restore namelist and broadcast to all
  ! processes. Store results in the ecosys_restore_vars

  use marbl_kinds_mod, only : char_len, char_len, int_kind, i4
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size
  use marbl_namelist_mod    , only : marbl_namelist
  use marbl_logging         , only : marbl_log_type
  use marbl_logging         , only : status_msg

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(out) :: restore_file_varnames
  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(*), parameter :: subname = "ecosys_restore:read_namelist"
  integer(int_kind) :: nml_error, t
  logical(log_kind) :: spatial_variability_from_file
  character(len=marbl_nl_buffer_size) :: tmp_nl_buffer

  !-----------------------------------------------------------------------

  namelist /ecosys_restore_nml/ &
       restore_short_names, &
       restore_filenames, &
       restore_file_varnames, &
       spatial_variability_from_file

  ! initialize namelist variables to default values
  restore_short_names = ''
  restore_filenames = ''
  restore_file_varnames = ''
  spatial_variability_from_file = .false.

  tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_restore_nml')
  read(tmp_nl_buffer, nml=ecosys_restore_nml, iostat=nml_error)
  if (nml_error /= 0) then
     call marbl_status_log%log_error("Error reading ecosys_restore_nml", subname)
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
  end do

  ! assign namelist variables to corresponding instance variables
  this%spatial_variability_from_file = spatial_variability_from_file

end subroutine read_namelist

!*****************************************************************************

subroutine initialize_restore_read_vars(this, restore_short_names, restore_filenames, &
     restore_file_varnames, ind_name_table, marbl_status_log)
  !
  ! Read the ecosys_restore namelist and broadcast to all
  ! processes. Store results in the ecosys_restore_vars%tracers(i)%restore_file_info
  !
  ! FIXME(bja, 2014-10) this%restore_any_tracer is set as a
  ! side-effect of this function! This isn't clear and is a Bad Thing (tm)
  !
  ! NOTE(bja, 2014-10) assumes that restore file is ALWAYS netcdf!
  use marbl_kinds_mod, only : char_len, int_kind, log_kind
  use passive_tracer_tools, only : ind_name_pair, name_to_ind ! FIXME
  use marbl_logging, only : marbl_log_type
  use marbl_logging, only : status_msg

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_short_names
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_filenames
  character(len=char_len), dimension(ecosys_tracer_cnt), intent(in) :: restore_file_varnames
  type(ind_name_pair), dimension(:), intent(in) :: ind_name_table

  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind)       :: t, tracer_index
  character(len=char_len) :: file_format
  character(len=char_len) :: message
  character(*), parameter :: subname='ecosys_restore:initialize_restore_read_vars'

  !-----------------------------------------------------------------------

  file_format = 'nc'

  this%restore_any_tracer = .false.

  ! reinitialize any restore vars requested by the user
  do t = 1, size(restore_short_names)
     if (len(trim(restore_short_names(t))) > 0) then

        ! attempt to map the user specified tracer name to tracer index
        tracer_index = name_to_ind(trim(restore_short_names(t)), ind_name_table)
        if (tracer_index > 0) then

           this%restore_any_tracer = .true.
           this%tracers(tracer_index)%restore = .true.
           call this%set_tracer_read_metadata(tracer_index, &
                mod_varname=restore_short_names(t),         &
                filename=restore_filenames(t),              &
                file_varname=restore_file_varnames(t),      &
                file_fmt=file_format)
           write(status_msg, "(3A)") "Setting up restoring for '", trim(restore_short_names(t)), "'"
           call marbl_status_log%log_noerror(status_msg, subname)

        else

           write(message, *) "ERROR: Could not find user requested restore variable '", &
                              trim(restore_short_names(t)), "'"
           call marbl_status_log%log_error(message, subname)
           return

        end if
     end if
  end do

end subroutine initialize_restore_read_vars

!*****************************************************************************

subroutine initialize_restoring_timescale(this, nl_buffer, zt,     &
                                          marbl_status_log)
  !
  ! Initialize the spatially varying restoring timescale by 
  !
  use marbl_kinds_mod, only : i4, r8
  use marbl_logging, only : marbl_log_type
  use marbl_logging, only : error_msg
  use marbl_namelist_mod    , only : marbl_nl_cnt
  use marbl_namelist_mod    , only : marbl_nl_buffer_size

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
  real (r8), dimension(:), intent(in) :: zt
  type(marbl_log_type), intent(inout) :: marbl_status_log

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  character(*), parameter :: subname = 'ecosys_restore:initialize_restoring_timescale'

  if (this%restore_any_tracer) then
     if (this%spatial_variability_from_file) then
        call this%timescale_file%init(nl_buffer, marbl_status_log)
     else
        call this%timescale_interp%init(nl_buffer, zt, marbl_status_log)
     end if
     if (marbl_status_log%labort_marbl) then
       error_msg = "error code returned from this%timescale_[file|interp]%init"
       call marbl_status_log%log_error(error_msg, subname)
       return
     end if
     call marbl_status_log%log_noerror("Restoring timescale set.", subname)
  end if

end subroutine initialize_restoring_timescale

!*****************************************************************************

subroutine read_restoring_fields(this, LAND_MASK)
  !
  !  load restoring fields if required
  !
  use marbl_kinds_mod, only : int_kind, c0
  use blocks, only : nx_block, ny_block ! FIXME
  use domain_size, only : km ! FIXME
  use domain, only : nblocks_clinic ! FIXME
  use prognostic, only : max_blocks_clinic ! FIXME
  use grid, only : KMT ! FIXME
  use passive_tracer_tools, only : read_field ! FIXME

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  logical (log_kind), dimension(:,:,:), intent(in) :: LAND_MASK

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer (int_kind) :: i,j    ! indoces for looping over horizontal grid
  integer (int_kind) :: k      ! index for looping over levels
  integer (int_kind) :: iblock ! index for looping over blocks
  integer (int_kind) :: tracer_index

  !-----------------------------------------------------------------------

  do tracer_index = 1, size(this%tracers)
     associate(tracer => this%tracers(tracer_index))

       if (tracer%restore) then
!!$          write(stdout, *) "Reading restore data for ", trim(tracer%restore_file_info%mod_varname)
          allocate(tracer%data(nx_block, ny_block, km, max_blocks_clinic))

          call read_field(tracer%restore_file_info%file_fmt, &
               tracer%restore_file_info%filename, &
               tracer%restore_file_info%file_varname, &
               tracer%data)

          do iblock = 1, nblocks_clinic
             do k = 1, km
                do j=1,ny_block
                   do i=1,nx_block
                      if (LAND_MASK(i,j,iblock) .and. (k.le.KMT(i,j,iblock))) then
                         tracer%data(i,j,k,iblock) = tracer%data(i,j,k,iblock) * tracer%restore_file_info%scale_factor
                      else
                         tracer%data(i,j,k,iblock) = c0
                      end if
                   end do
                end do
             end do
          end do
       endif
     end associate
  end do

end subroutine read_restoring_fields

!*****************************************************************************

subroutine restore_tracers(this, tracer_cnt, vert_level, x_index, y_index, &
     block_id, local_data, restore_data)
  !
  !  restore a variable if required
  !
  use marbl_kinds_mod       , only : c0, r8, int_kind, log_kind

  implicit none
  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------
  class(ecosys_restore_type), intent(inout) :: this
  integer(int_kind), intent(in) :: tracer_cnt
  integer(int_kind), intent(in) :: vert_level
  integer(int_kind), intent(in) :: x_index
  integer(int_kind), intent(in) :: y_index
  integer(int_kind), intent(in) :: block_id
  real(r8), intent(in) :: local_data(tracer_cnt)
  real(r8), intent(out) :: restore_data(tracer_cnt)
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind) :: n
  !-----------------------------------------------------------------------

  associate( &
       tracer => this%tracers(:), &
       rtau => this%timescale_file%restore_rtau, &
       restore_max_level => this%timescale_file%restore_max_level, &
       inv_restoring_time_scale => this%timescale_interp%inv_restoring_time_scale &
       )

    do n=1,tracer_cnt
       if (tracer(n)%restore) then
          if (this%spatial_variability_from_file) then
             restore_data(n) = rtau(x_index, y_index, block_id) * &
                  merge((tracer(n)%data(x_index, y_index, vert_level, block_id) - local_data(n)), &
                  c0, vert_level <= restore_max_level(x_index, y_index, block_id))
          else
             restore_data(n) = (tracer(n)%data(x_index, y_index, vert_level, block_id) - local_data(n)) * &
                  inv_restoring_time_scale(vert_level)
          endif
       else
          restore_data(n) = c0
       endif
    enddo
  end associate

end subroutine restore_tracers

!*****************************************************************************

subroutine set_tracer_read_metadata(this, index, &
     mod_varname, filename, file_varname, file_fmt, &
     scale_factor, default_val)

  use marbl_kinds_mod , only : char_len, r8, c0, c1

  implicit none

  !  initialize a tracer_read type to common default values.

  class(ecosys_restore_type), intent(inout) :: this
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

!*****************************************************************************

end module ecosys_restore_mod
