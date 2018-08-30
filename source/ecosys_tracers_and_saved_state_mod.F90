module ecosys_tracers_and_saved_state_mod

  ! !DESCRIPTION:
  !  This module provides routines to read data to initialize the tracers
  !  introduced to POP via the MARBL library. It also sets up data-types to
  !  populate MARBL's saved_state type (used to pass data between the surface
  !  and interior computations).

  use constants, only : c0, c1

  use kinds_mod, only : r8, int_kind, log_kind, char_len, char_len_long

  use io_types, only : stdout

  use io_tools, only : document

  use exit_mod, only : sigAbort, exit_POP

  use passive_tracer_tools, only : tracer_read

  use communicate, only : my_task, master_task

  implicit none
  private

  ! this struct is necessary because there is some global state
  ! that needs to be preserved for from one time step to the next
  type, public :: ecosys_saved_state_type
     integer :: rank
     character(len=char_len) :: short_name
     character(len=char_len) :: file_varname
     character(len=char_len) :: units

     ! (nx_block, ny_block, max_blocks_clinic)
     real (r8), allocatable, dimension(:,:,:)   :: field_2d
     ! (km, nx_block, ny_block, max_blocks_clinic) [column-first for MARBL!]
     real (r8), allocatable, dimension(:,:,:,:) :: field_3d
  end type ecosys_saved_state_type

  type(ecosys_saved_state_type), allocatable, target, dimension(:), public :: surface_flux_saved_state
  type(ecosys_saved_state_type), allocatable, target, dimension(:), public :: interior_tendency_saved_state
  ! Transpose field_3d into this array before writing restart
  ! (nx_block, ny_block, km, max_blocks_clinic)
  real(r8), allocatable, target, dimension(:,:,:,:), public :: saved_state_field_3d

  !-----------------------------------------------------------------------
  ! public variables
  !-----------------------------------------------------------------------

  ! # of tracers expected from MARBL
  integer(int_kind), parameter, public :: marbl_tracer_cnt = MARBL_NT

  ! Indices of tracers needed for virtual flux or river fluxes
  integer(int_kind), public :: dic_ind
  integer(int_kind), public :: alk_ind
  integer(int_kind), public :: dic_alt_co2_ind
  integer(int_kind), public :: alk_alt_co2_ind
  integer(int_kind), public :: di13c_ind
  integer(int_kind), public :: di14c_ind
  integer(int_kind), public :: o2_ind
  integer(int_kind), public :: no3_ind
  integer(int_kind), public :: po4_ind
  integer(int_kind), public :: don_ind
  integer(int_kind), public :: donr_ind
  integer(int_kind), public :: dop_ind
  integer(int_kind), public :: dopr_ind
  integer(int_kind), public :: sio3_ind
  integer(int_kind), public :: fe_ind
  integer(int_kind), public :: doc_ind
  integer(int_kind), public :: docr_ind
  integer(int_kind), public :: do13ctot_ind
  integer(int_kind), public :: do14ctot_ind

  !---------------------------------------------------------------------
  !  Private variables read in via &ecosys_tracer_init_nml
  !---------------------------------------------------------------------

  character(char_len), target :: init_ecosys_option           ! namelist option for initialization of bgc
  character(char_len), target :: init_ecosys_init_file        ! filename for option 'file'
  character(char_len), target :: init_ecosys_init_file_fmt    ! file format for option 'file'
  type(tracer_read),   target :: tracer_init_ext(MARBL_NT) ! namelist variable for initializing tracers

  character(char_len), target :: ciso_init_ecosys_option        ! option for initialization of bgc
  character(char_len), target :: ciso_init_ecosys_init_file     ! filename for option 'file'
  character(char_len), target :: ciso_init_ecosys_init_file_fmt ! file format for option 'file'
  type(tracer_read),   target :: ciso_tracer_init_ext(MARBL_NT) ! namelist variable for initializing tracers

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ecosys_tracers_and_saved_state_init
  public :: ecosys_saved_state_setup
  public :: ecosys_saved_state_write_restart
  public :: set_defaults_tracer_read

  !-----------------------------------------------------------------------

Contains

  !-----------------------------------------------------------------------

  subroutine ecosys_tracers_and_saved_state_init(marbl_instance,              &
                                                 ecosys_driver_ind_begin,     &
                                                 ciso_on,                     &
                                                 init_ts_file_fmt,            &
                                                 read_restart_filename,       &
                                                 tracer_d_module,             &
                                                 module_name,                 &
                                                 tracer_nml,                  &
                                                 land_mask,                   &
                                                 TRACER_MODULE,               &
                                                 errorCode)

    use POP_ErrorMod, only : POP_Success, POP_ErrorSet

    use io_read_fallback_mod, only : io_read_fallback_register_tracer
    use io_read_fallback_mod, only : io_read_fallback_register_field

    use passive_tracer_tools, only : rest_read_tracer_block
    use passive_tracer_tools, only : file_read_single_tracer
    use passive_tracer_tools, only : read_field
    use passive_tracer_tools, only : ind_name_pair

    use prognostic, only : curtime
    use prognostic, only : oldtime
    use prognostic, only : newtime
    use prognostic, only : tracer_field_type => tracer_field

    use time_management, only : check_time_flag
    use time_management, only : eval_time_flag

    use grid, only : fill_points
    use grid, only : n_topo_smooth
    use grid, only : KMT

    use domain, only : nblocks_clinic
    use domain_size, only : km, nt

    use constants, only : delim_fmt, char_blank, ndelim_fmt

    use ecosys_forcing_saved_state_mod, only : ecosys_forcing_saved_state_init

    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_state_init

    use marbl_interface, only : marbl_interface_class

    type(marbl_interface_class),          intent(in)    :: marbl_instance
    integer (int_kind),                   intent(in)    :: ecosys_driver_ind_begin ! starting index of ecosys tracers in global tracer
    logical,                              intent(in)    :: ciso_on
    character(len=*),                     intent(in)    :: init_ts_file_fmt        ! format (bin or nc) for input file
    character(len=*),                     intent(in)    :: read_restart_filename   ! file name for restart file
    type(tracer_field_type),              intent(in)    :: tracer_d_module(:)      ! descriptors for each tracer
    character(len=*),   dimension(:),     intent(in)    :: module_name
    character(len=*),                     intent(in)    :: tracer_nml
    logical(log_kind) , dimension(:,:,:), intent(in)    :: land_mask
    real(r8),                             intent(inout) :: tracer_module(:,:,:,:,:,:)
    integer(int_kind),                    intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_tracers_and_saved_state_mod:ecosys_tracers_and_saved_state_init'

    integer                      :: n, k, iblock
    character(len=char_len)      :: init_option, init_file_fmt
    character(len=char_len)      :: ecosys_restart_filename            ! modified file name for restart file
    integer(int_kind)            :: nml_error                          ! error flag for nml read
    character(len=char_len_long) :: ioerror_msg
    type(tracer_read), dimension(marbl_tracer_cnt) :: tracer_inputs   ! metadata about file to read

    !-----------------------------------------------------------------------

    namelist /ecosys_tracer_init_nml/                                         &
         tracer_init_ext, init_ecosys_option,                                 &
         init_ecosys_init_file, init_ecosys_init_file_fmt,                    &
         ciso_init_ecosys_option, ciso_init_ecosys_init_file,                 &
         ciso_init_ecosys_init_file_fmt, ciso_tracer_init_ext

    !-----------------------------------------------------------------------
    !  Set defaults for &ecosys_tracer_init_nml
    !-----------------------------------------------------------------------

    do n = 1, size(tracer_init_ext)
       call set_defaults_tracer_read(tracer_init_ext(n))
    end do
    if (ciso_on) then
      do n = 1, size(ciso_tracer_init_ext)
         call set_defaults_tracer_read(ciso_tracer_init_ext(n))
      end do
    end if

    init_ecosys_option = 'unknown'
    init_ecosys_init_file = 'unknown'
    init_ecosys_init_file_fmt = 'bin'

    ciso_init_ecosys_option                 = 'unknown'
    ciso_init_ecosys_init_file              = 'unknown'
    ciso_init_ecosys_init_file_fmt          = 'bin'

    !-----------------------------------------------------------------------
    !  Read &ecosys_tracer_init_nml
    !-----------------------------------------------------------------------

    read(tracer_nml, nml=ecosys_tracer_init_nml, iostat=nml_error, iomsg=ioerror_msg)
    if (nml_error /= 0) then
       write(stdout, *) subname, ": process ", my_task, ": namelist read error: ", nml_error, " : ", ioerror_msg
       call exit_POP(sigAbort, 'ERROR reading ecosys_tracer_init_nml from buffer.')
    end if

    if (my_task == master_task) then
       write(stdout,*)
       write(stdout,ndelim_fmt)
       write(stdout,*)
       write(stdout,*) ' ecosys_tracer_init_nml namelist settings:'
       write(stdout,*)
       write(stdout,ecosys_tracer_init_nml)
       write(stdout,*)
       write(stdout,delim_fmt)
    endif

    !-----------------------------------------------------------------------
    ! Set default values for tracer_inputs
    !-----------------------------------------------------------------------

    do n=1,marbl_tracer_cnt
      tracer_inputs(n)%mod_varname  = tracer_d_module(n)%short_name
      tracer_inputs(n)%file_varname = tracer_d_module(n)%short_name
      tracer_inputs(n)%scale_factor = c1
      tracer_inputs(n)%default_val  = c0
      select case (trim(module_name(n)))
        case('ecosys')
          tracer_inputs(n)%filename = init_ecosys_init_file
          tracer_inputs(n)%file_fmt = init_ecosys_init_file_fmt
        case('ciso')
          tracer_inputs(n)%filename = ciso_init_ecosys_init_file
          tracer_inputs(n)%file_fmt = ciso_init_ecosys_init_file_fmt
        case DEFAULT
          call document(subname, 'unknown module_name', trim(module_name(n)))
          call exit_POP(sigAbort, 'Stopping in ' // subname)
      end select
    end do

    !-----------------------------------------------------------------------
    ! Update tracer_inputs from tracer_init_ext and ciso_tracer_init_ext
    !-----------------------------------------------------------------------

    call set_tracer_read(tracer_init_ext, tracer_inputs)
    if (ciso_on) then
      call set_tracer_read(ciso_tracer_init_ext, tracer_inputs)
    end if

    !-----------------------------------------------------------------------
    !  Enable tracers and saved state to be read from older IC/restart files
    !-----------------------------------------------------------------------

    call io_read_fallback_register_tracer(tracername='DIC_ALT_CO2', &
       fallback_opt='alt_field', alt_tracername='DIC')

    call io_read_fallback_register_tracer(tracername='ALK_ALT_CO2', &
       fallback_opt='alt_field', alt_tracername='ALK')

    call io_read_fallback_register_tracer(tracername='DOCr', &
       fallback_opt='const', const_val=38.0_r8)

    call io_read_fallback_register_tracer(tracername='Lig', &
       fallback_opt='const', const_val=2.0e-3_r8)

    call io_read_fallback_register_tracer(tracername='spP', &
       fallback_opt='alt_field', alt_tracername='spC', scalefactor=c1/117.0_r8)

    call io_read_fallback_register_tracer(tracername='diatP', &
       fallback_opt='alt_field', alt_tracername='diatC', scalefactor=c1/117.0_r8)

    call io_read_fallback_register_tracer(tracername='diazP', &
       fallback_opt='alt_field', alt_tracername='diazC', scalefactor=0.32_r8/117.0_r8)

    call io_read_fallback_register_field(fieldname='MARBL_PH_SURF', &
       fallback_opt='const', const_val=c0)

    call io_read_fallback_register_field(fieldname='MARBL_PH_SURF_ALT_CO2', &
       fallback_opt='const', const_val=c0)

    call io_read_fallback_register_field(fieldname='MARBL_PH_3D', &
       fallback_opt='const', const_val=c0)

    call io_read_fallback_register_field(fieldname='MARBL_PH_3D_ALT_CO2', &
       fallback_opt='const', const_val=c0)

    !-----------------------------------------------------------------------
    !  initialize saved state
    !-----------------------------------------------------------------------

    ecosys_restart_filename = char_blank

    select case (trim(init_ecosys_option))

    case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

       init_file_fmt = init_ecosys_init_file_fmt
       if (init_ecosys_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call document(subname, 'no restart file to read ecosys from')
             call exit_POP(sigAbort, 'Stopping in ' // subname)
          endif
          ecosys_restart_filename = read_restart_filename
          init_file_fmt = init_ts_file_fmt
       else  ! do not read from TS restart file
          ecosys_restart_filename = trim(init_ecosys_init_file)
       endif

       call read_restart_saved_state(init_file_fmt, ecosys_restart_filename,  &
            surface_flux_saved_state)
       call read_restart_saved_state(init_file_fmt, ecosys_restart_filename,  &
            interior_tendency_saved_state)

    case ('file', 'ccsm_startup')
       call set_saved_state_zero(surface_flux_saved_state)
       call set_saved_state_zero(interior_tendency_saved_state)

    case default
       call document(subname, 'unknown init_ecosys_option', init_ecosys_option)
       call exit_POP(sigAbort, 'Stopping in ' // subname)

    end select

    call ecosys_forcing_saved_state_init(ecosys_restart_filename)
    call ecosys_running_mean_state_init(marbl_instance, ecosys_restart_filename)

    !-----------------------------------------------------------------------
    !  initialize tracers
    !-----------------------------------------------------------------------

    do n=1,marbl_tracer_cnt

      ! Is tracer read from restart file or initial condition?
      ! What is the file name and format?
      select case (trim(module_name(n)))
        case('ecosys')
          init_option = init_ecosys_option
          ecosys_restart_filename = trim(init_ecosys_init_file)
          init_file_fmt = init_ecosys_init_file_fmt

        case('ciso')
          init_option = ciso_init_ecosys_option
          ecosys_restart_filename = trim(ciso_init_ecosys_init_file)
          init_file_fmt = ciso_init_ecosys_init_file_fmt

        case DEFAULT
          call document(subname, 'unknown module_name', trim(module_name(n)))
          call exit_POP(sigAbort, 'Stopping in ' // subname)

      end select

      select case (trim(init_option))

        ! For restart run, either read from specified file or TS restart file
        case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')
           if (ecosys_restart_filename == 'same_as_TS') then
              if (read_restart_filename == 'undefined') then
                 call document(subname, 'no restart file to read ',           &
                      trim(module_name(n)))
                 call exit_POP(sigAbort, 'Stopping in ' // subname)
              endif
              ecosys_restart_filename = read_restart_filename
              init_file_fmt = init_ts_file_fmt
           endif

           call rest_read_tracer_block(ecosys_driver_ind_begin+n-1, &
                                       init_file_fmt,               &
                                       ecosys_restart_filename,     &
                                       tracer_d_module(n:n),        &
                                       TRACER_MODULE(:,:,:,n:n,:,:))

        case ('file', 'ccsm_startup')
           call file_read_single_tracer(tracer_inputs, TRACER_MODULE, n)

           if (n_topo_smooth > 0) then
             do k=1, km
                call fill_points(k, TRACER_MODULE(:, :, k, n, oldtime, :), errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(oldtime)')
                   return
                endif

                call fill_points(k, TRACER_MODULE(:, :, k, n, curtime, :), errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(newtime)')
                   return
                endif

             enddo
           endif

        case default
         call document(subname, 'unknown init_option', init_option)
         call exit_POP(sigAbort, 'Stopping in ' // subname)

      end select

    end do

    do iblock = 1, nblocks_clinic
       do n = 1, marbl_tracer_cnt
          do k = 1, km
             where (.not. land_mask(:, :, iblock) .or. k > KMT(:, :, iblock))
                TRACER_MODULE(:, :, k, n, curtime, iblock) = c0
                TRACER_MODULE(:, :, k, n, oldtime, iblock) = c0
             end where
          end do
       end do
    end do

  end subroutine ecosys_tracers_and_saved_state_init

  !-----------------------------------------------------------------------

  subroutine set_defaults_tracer_read(tracer_read_var, mod_varname, filename, &
             file_varname, scale_factor, default_val, file_fmt)

    type(tracer_read),          intent(inout) :: tracer_read_var
    character(len=*), optional, intent(in)    :: mod_varname
    character(len=*), optional, intent(in)    :: filename
    character(len=*), optional, intent(in)    :: file_varname
    real(r8),         optional, intent(in)    :: scale_factor
    real(r8),         optional, intent(in)    :: default_val
    character(len=*), optional, intent(in)    :: file_fmt

    if (present(mod_varname)) then
      tracer_read_var%mod_varname = trim(mod_varname)
    else
      tracer_read_var%mod_varname = 'unknown'
    end if

    if (present(filename)) then
      tracer_read_var%filename = trim(filename)
    else
      tracer_read_var%filename = 'unknown'
    end if

    if (present(file_varname)) then
      tracer_read_var%file_varname = trim(file_varname)
    else
      tracer_read_var%file_varname = 'unknown'
    end if

    if (present(scale_factor)) then
      tracer_read_var%scale_factor = scale_factor
    else
      tracer_read_var%scale_factor = c1
    end if

    if (present(default_val)) then
      tracer_read_var%default_val = default_val
    else
      tracer_read_var%default_val = c0
    end if

    if (present(file_fmt)) then
      tracer_read_var%file_fmt = trim(file_fmt)
    else
      tracer_read_var%file_fmt = 'nc'
    end if

  end subroutine set_defaults_tracer_read

  !-----------------------------------------------------------------------

  subroutine read_restart_saved_state(file_fmt, filename, state)

    use passive_tracer_tools, only : read_field
    use blocks              , only : nx_block, ny_block
    use domain_size         , only : km, max_blocks_clinic

    character(len=*), intent(in) :: file_fmt
    character(len=*), intent(in) :: filename
    real(r8) :: tmp_field_3d(nx_block, ny_block, km, max_blocks_clinic)
    type(ecosys_saved_state_type), dimension(:), intent(inout) :: state

    integer :: k, n

    do n=1,size(state)
      select case (state(n)%rank)
        case (2)
          call read_field(file_fmt, filename, state(n)%file_varname,        &
               state(n)%field_2d(:,:,:))
        case (3)
          call read_field(file_fmt, filename, state(n)%file_varname,        &
               tmp_field_3d)
          ! Transpose data after reading
          do k=1, km
            state(n)%field_3d(k,:,:,:) = tmp_field_3d(:,:,k,:)
          end do
      end select
    end do

  end subroutine read_restart_saved_state

  !-----------------------------------------------------------------------

  subroutine set_saved_state_zero(state)

    type(ecosys_saved_state_type), dimension(:), intent(inout) :: state
    integer :: n

    do n=1,size(state)
      select case (state(n)%rank)
        case (2)
          state(n)%field_2d = c0
        case (3)
          state(n)%field_3d = c0
      end select
    end do

  end subroutine set_saved_state_zero

  !-----------------------------------------------------------------------

  subroutine ecosys_saved_state_setup(state, marbl_state)

    use blocks, only : nx_block, ny_block
    use domain_size, only : km, max_blocks_clinic
    use marbl_interface_public_types, only : marbl_saved_state_type
    use io_tools, only : document

    type(ecosys_saved_state_type), allocatable, intent(out) :: state(:)
    type(marbl_saved_state_type),               intent(in)  :: marbl_state

    character(len=*), parameter :: subname =                                  &
                  'ecosys_tracers_and_saved_state_mod:ecosys_saved_state_setup'
    integer :: num_fields
    integer :: n

    num_fields = marbl_state%saved_state_cnt
    allocate(state(num_fields))
    do n=1, num_fields
      state(n)%rank = marbl_state%state(n)%rank
      select case (state(n)%rank)
        case (2)
          allocate(state(n)%field_2d(nx_block, ny_block, max_blocks_clinic))
        case (3)
          select case (trim(marbl_state%state(n)%vertical_grid))
            case ('layer_avg')
              allocate(state(n)%field_3d(km, nx_block, ny_block, max_blocks_clinic))
            case DEFAULT
          call document(subname, 'n', n)
          call document(subname, 'marbl_state(n)%vertical_grid',              &
                        trim(marbl_state%state(n)%vertical_grid))
          call document(subname, 'Invalid vert grid from MARBL saved state')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
          end select
        case DEFAULT
          call document(subname, 'n', n)
          call document(subname, 'state(n)%rank', state(n)%rank)
          call document(subname, 'Invalid rank from MARBL saved state')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
      end select
      state(n)%short_name = trim(marbl_state%state(n)%short_name)
      write(state(n)%file_varname, "(2A)") 'MARBL_', trim(state(n)%short_name)
      state(n)%units = trim(marbl_state%state(n)%units)
    end do

  end subroutine ecosys_saved_state_setup

  !-----------------------------------------------------------------------

  subroutine ecosys_saved_state_write_restart(restart_file, action, state, iodesc)

    use domain     , only : nblocks_clinic
    use domain_size, only : nx_global
    use domain_size, only : ny_global
    use domain_size, only : km
    use constants  , only : field_loc_center
    use constants  , only : field_type_scalar
    use io         , only : data_set
    use io         , only : datafile
    use io_types   , only : io_dim
    use io_types   , only : io_field_desc
    use io_types   , only : add_attrib_file
    use io_types   , only : construct_io_dim
    use io_types   , only : construct_io_field
    use io_tools   , only : document

    implicit none

    type (datafile),                                     intent(inout) :: restart_file
    character(len=*),                                    intent(in)    :: action
    type(ecosys_saved_state_type), target, dimension(:), intent(in)    :: state
    type(io_field_desc),                   dimension(:), intent(inout) :: iodesc
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname='ecosys_tracers_and_saved_state_mod:ecosys_saved_state_write_restart'
    character(len=char_len)     :: log_message

    integer(int_kind) :: k, n
    type(io_dim)      :: i_dim, j_dim ! dimension descriptors
    type(io_dim)      :: k_dim        ! dimension descriptor for vertical levels

    if (trim(action) .eq. 'define') then
      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)
      k_dim = construct_io_dim('k', km)

      do n=1,size(state)
        select case (state(n)%rank)
          case (2)
            iodesc(n) = construct_io_field(trim(state(n)%file_varname), i_dim, &
                        j_dim, units = trim(state(n)%units), grid_loc = '2110',&
                        field_loc  = field_loc_center,                         &
                        field_type = field_type_scalar,                        &
                        d2d_array  = state(n)%field_2d(:,:,1:nblocks_clinic))
          case (3)
            iodesc(n) = construct_io_field(trim(state(n)%file_varname), i_dim, &
                        j_dim, k_dim,                                          &
                        units = trim(state(n)%units), grid_loc = '3111',       &
                        field_loc  = field_loc_center,                         &
                        field_type = field_type_scalar,                        &
                        d3d_array  = saved_state_field_3d(:,:,:,1:nblocks_clinic))
        end select
        write(log_message, "(3A)") "Setting up IO field for ",                  &
                                   trim(state(n)%file_varname), " (in restart)"
        call document(subname, log_message)
        call data_set (restart_file, 'define', iodesc(n))
      end do
    else if (trim(action) .eq. 'write') then
      do n=1,size(state)
        if (state(n)%rank .eq. 3) then
          do k=1,km
            saved_state_field_3d(:,:,k,:) = state(n)%field_3d(k, :,:,:)
          end do
        end if
        call data_set (restart_file, 'write', iodesc(n))
      end do
    end if

  end subroutine ecosys_saved_state_write_restart

  !-----------------------------------------------------------------------

  subroutine set_tracer_read(tracer_ext, tracer_inputs)

    type(tracer_read), dimension(:), intent(in)    :: tracer_ext
    type(tracer_read), dimension(:), intent(inout) :: tracer_inputs

    character(len=*), parameter :: subname = 'ecosys_tracers_and_saved_state_mod:set_tracer_read'
    integer :: n, ind, tracer_ind
    character(len=char_len) :: err_msg

    do n=1,size(tracer_ext)
      if (trim(tracer_ext(n)%mod_varname).ne.'unknown') then
        tracer_ind = 0
        do ind = 1, size(tracer_inputs)
          if (trim(tracer_ext(n)%mod_varname).eq.                             &
              trim(tracer_inputs(ind)%mod_varname)) then
            tracer_ind = ind
            exit
          end if
        end do

        if (tracer_ind.eq.0) then
          write(err_msg, "(A,1X,A)") 'No tracer defined with name',           &
                                     trim(tracer_ext(n)%mod_varname)
          call document(subname, err_msg)
          call exit_POP(sigAbort, 'Stopping in ' // subname)
        end if

        if (trim(tracer_ext(n)%filename).ne.'unknown')                        &
          tracer_inputs(ind)%filename = tracer_ext(n)%filename
        if (trim(tracer_ext(n)%file_varname).ne.'unknown')                    &
          tracer_inputs(ind)%file_varname = tracer_ext(n)%file_varname
        tracer_inputs(ind)%scale_factor = tracer_ext(n)%scale_factor
        tracer_inputs(ind)%default_val  = tracer_ext(n)%default_val
      end if
    end do

  end subroutine set_tracer_read

  !-----------------------------------------------------------------------

end module ecosys_tracers_and_saved_state_mod
