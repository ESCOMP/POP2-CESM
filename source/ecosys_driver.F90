! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_driver

  ! !DESCRIPTION:
  !  This module provides support for the ecosystem module and dependend tracer modules
  !  The base model calls subroutines in passive_tracers, which then call
  !  this module if ecosys_on is true. Ecosys_driver then calls subroutines
  !  in individual ecosystem modules (so far marbl_mod and marbl_ciso_mod)
  !
  !  Written by: Alexandra Jahn, NCAR, Nov/Dec 2012

  use POP_KindsMod
  use POP_ErrorMod
  use POP_IOUnitsMod

  use kinds_mod                 , only : r8, int_kind, log_kind, char_len, char_len_long
  use blocks                    , only : block, nx_block, ny_block
  use domain                    , only : nblocks_clinic
  use domain                    , only : distrb_clinic
  use domain_size               , only : max_blocks_clinic, km, nt
  use io_types                  , only : stdout, nml_in, nml_filename
  use io_tools                  , only : document
  use exit_mod                  , only : sigAbort, exit_POP
  use constants                 , only : c0, c1, p5, delim_fmt, char_blank, ndelim_fmt
  use communicate               , only : my_task, master_task

  use marbl_sizes               , only : num_surface_forcing_fields 
  use marbl_sizes               , only : num_interior_forcing_fields 

  use marbl_config_mod          , only : lflux_gas_co2

  use marbl_logging             , only : marbl_log_type

  use marbl_interface           , only : marbl_interface_class

  use marbl_namelist_mod        , only : marbl_nl_in_size
  use marbl_namelist_mod        , only : marbl_nl_cnt
  use marbl_namelist_mod        , only : marbl_nl_buffer_size
  use marbl_namelist_mod        , only : marbl_nl_split_string
  use marbl_namelist_mod        , only : marbl_namelist

  use ecosys_tavg               , only : ecosys_tavg_init
  use ecosys_tavg               , only : ecosys_tavg_accumulate
  use ecosys_tavg               , only : ecosys_tavg_accumulate_flux

  use passive_tracer_tools      , only : tracer_read

  use timers                    , only : timer_start
  use timers                    , only : timer_stop
  use timers                    , only : get_timer

  use ecosys_tracers_and_saved_state_mod, only : ecosys_tracers_and_saved_state_init
  use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_setup
  use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_construct_io_fields
  use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_type
  use ecosys_tracers_and_saved_state_mod, only : saved_state_surf
  use ecosys_tracers_and_saved_state_mod, only : saved_state_interior
  use ecosys_tracers_and_saved_state_mod, only : dic_ind, alk_ind, dic_alt_co2_ind
  use ecosys_tracers_and_saved_state_mod, only : di13c_ind, di14c_ind

  ! Provide marbl_tracer_cnt to passive_tracers.F90 (as ecosys_tracer_cnt)
  use ecosys_tracers_and_saved_state_mod, only : ecosys_tracer_cnt => marbl_tracer_cnt
  ! Provide ecosys_forcing_tracer_ref_val to passive_tracers.F90
  ! (as ecosys_driver_tracer_ref_val)
  use ecosys_forcing_mod, only : ecosys_driver_tracer_ref_val => ecosys_forcing_tracer_ref_val

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_driver_init
  public :: ecosys_driver_set_interior
  public :: ecosys_driver_set_global_scalars
  public :: ecosys_driver_comp_global_averages
  public :: ecosys_driver_set_sflux
  public :: ecosys_driver_tavg_forcing
  public :: ecosys_driver_write_restart
  public :: ecosys_driver_unpack_source_sink_terms

  ! These are public so passive_tracers.F90 only uses ecosys_driver
  public :: ecosys_driver_tracer_ref_val
  public :: ecosys_tracer_cnt

  private :: ecosys_driver_init_rmean_var
  private :: ecosys_driver_update_scalar_rmeans

  !*****************************************************************************

  !-----------------------------------------------------------------------
  ! timers
  !-----------------------------------------------------------------------

  integer (int_kind) :: ecosys_interior_timer
  integer (int_kind) :: ecosys_set_sflux_timer
  integer (int_kind) :: ecosys_comp_CO3terms_timer
  
  !-----------------------------------------------------------------------
  ! module variables 
  !-----------------------------------------------------------------------

  type(marbl_interface_class) :: marbl_instances(max_blocks_clinic)

  integer (int_kind)  :: totChl_surf_nf_ind = 0                ! total chlorophyll in surface layer 
  integer (int_kind)  :: sflux_co2_nf_ind   = 0                ! air-sea co2 gas flux 
  integer (int_kind)  :: num_elements  = nx_block*ny_block     ! number of surface elements passed to marbl

  character (char_len)                       :: ecosys_tadvect_ctype                  ! advection method for ecosys tracers
  logical   (log_kind) , public              :: ecosys_qsw_distrb_const
  logical   (log_kind)                       :: ciso_on 
  logical   (log_kind) , allocatable         :: land_mask(:, :, :)
  real      (r8)       , allocatable         :: surface_forcing_diags(:, :, :, :)
  real      (r8) :: surface_forcing_outputs(nx_block, ny_block, max_blocks_clinic, 2)
  integer :: flux_co2_id
  integer :: totalChl_id

  ! Variables related to global averages

  real (r8), allocatable, target, dimension(:,:,:,:) :: glo_avg_fields_interior
  real (r8), allocatable, target, dimension(:,:,:,:) :: glo_avg_fields_surface
  real (r8), allocatable, dimension(:,:,:)           :: glo_avg_area_masked
  real (r8)                                          :: glo_avg_norm_fact

  ! FIXME : move the following variables into MARBL

  integer (int_kind), allocatable, target, dimension(:) :: glo_avg_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:) :: glo_avg_rmean_ind_surface
  integer (int_kind), allocatable, target, dimension(:) :: glo_scalar_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:) :: glo_scalar_rmean_ind_surface

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_driver_init( &
       ciso_active_flag,         &
       init_ts_file_fmt,         &
       read_restart_filename,    &
       tracer_d_module,          &
       tracer_module,            &
       tadvect_ctype,            &
       errorCode)

    ! !DESCRIPTION:
    !  Initialize ecosys_driver passive tracers. This involves:
    !  1) setting ecosys and ecosys_ciso module index bounds
    !  2) calling ecosys and ecosys_ciso module init subroutine

    use grid                  , only : dz
    use grid                  , only : zt
    use grid                  , only : zw
    use grid                  , only : region_mask
    use grid                  , only : TAREA
    use global_reductions     , only : global_sum
    use constants             , only : field_loc_center
    use broadcast             , only : broadcast_scalar
    use time_management       , only : init_time_flag
    use passive_tracer_tools  , only : set_tracer_indices
    use named_field_mod       , only : named_field_register
    use named_field_mod       , only : named_field_set
    use prognostic            , only : curtime
    use prognostic            , only : oldtime
    use prognostic            , only : tracer_field_type => tracer_field
    use mcog                  , only : mcog_nbins
    use registry              , only : registry_match
    use named_field_mod       , only : named_field_get_index
    use named_field_mod       , only : named_field_get
    use named_field_mod       , only : named_field_register
    use named_field_mod       , only : named_field_set
    use running_mean_mod      , only : running_mean_get_var
    use ecosys_forcing_mod    , only : ecosys_forcing_init
    use marbl_logging         , only : marbl_log_type

    implicit none

    logical                  , intent(in)    :: ciso_active_flag      ! set ciso_on
    character (*)            , intent(in)    :: init_ts_file_fmt      ! format (bin or nc) for input file
    character (*)            , intent(in)    :: read_restart_filename ! file name for restart file
    type (tracer_field_type) , intent(inout) :: tracer_d_module(:)    ! descriptors for each tracer
    real (r8)                , intent(inout) :: tracer_module(:,:,:,:,:,:)
    character (char_len)     , intent(out)   :: tadvect_ctype(:)      ! advection method for ecosys tracers
    integer (POP_i4)         , intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    ! FIXME (MNL): we should treat namelists in the driver in a POP-ish manner
    !              rather than a MARBL-ish manner, but the MARBL process is
    !              less error-prone (reading namelist on each task rather
    !              than relying on broadcasts to ensure every task has every
    !              namelist variable set correctly).
    type(marbl_log_type)                :: ecosys_status_log
    character(*), parameter             :: subname = 'ecosys_driver:ecosys_driver_init'
    character(char_len)                 :: log_message
    integer (int_kind)                  :: cumulative_nt, n, bid, k, i, j
    integer (int_kind)                  :: num_fields
    integer (int_kind)                  :: nml_error                          ! error flag for nml read
    integer (int_kind)                  :: iostat                             ! io status flag
    character (char_len)                :: sname, lname, units, coordinates
    character (4)                       :: grid_loc
    character(len=marbl_nl_buffer_size) :: nl_buffer(marbl_nl_cnt)
    character(len=marbl_nl_buffer_size) :: tmp_nl_buffer
    character(len=marbl_nl_in_size)     :: nl_str
    character(char_len_long)            :: ioerror_msg
    integer (int_kind)                  :: auto_ind                           ! autotroph functional group index
    integer (int_kind)                  :: iblock                             ! index for looping over blocks
    character (char_len)                :: ecosys_restart_filename            ! modified file name for restart file
    character (char_len)                :: init_file_fmt                      ! file format for option 'file'
    logical   (log_kind)                :: lmarginal_seas                     ! is ecosystem active in marginal seas?
    integer(int_kind)                   :: marbl_actual_tracer_cnt            ! # of tracers actually in MARBL
    integer (int_kind)                  :: glo_avg_field_cnt
    real (r8)                           :: rmean_val
    real (r8)                           :: fe_frac_dust
    real (r8)                           :: fe_frac_bc

    !-----------------------------------------------------------------------
    !  read in ecosys_driver namelist, to set namelist parameters that
    !  should be the same for all ecosystem-related modules
    !-----------------------------------------------------------------------

    namelist /ecosys_driver_nml/ &
         lmarginal_seas, ecosys_tadvect_ctype, ecosys_qsw_distrb_const

    errorCode = POP_Success
    ciso_on   = ciso_active_flag

    lmarginal_seas        = .true.
    ecosys_tadvect_ctype  = 'base_model'
    ecosys_qsw_distrb_const = .true.

    nl_buffer(:) = ''
    nl_str    = ''
    if (my_task == master_task) then
       ! read the namelist file into a buffer.
       open(unit=nml_in, file=nml_filename, action='read', access='stream', &
            form='unformatted', iostat=nml_error)
       if (nml_error == 0) then
          read(unit=nml_in, iostat=nml_error, iomsg=ioerror_msg) nl_str

          ! we should always reach the EOF to capture the entire file...
          if (.not. is_iostat_end(nml_error)) then
             write(stdout, '(a, a, i8)') subname, &
                  ": IO ERROR reading namelist file into buffer: ", nml_error
             write(stdout, '(a)') ioerror_msg
          else
             write(stdout, '(a, a, a)') "Read '", trim(nml_filename), "' until EOF."
          end if

          write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
               len_trim(nl_str), " characters."

          write(stdout, '(a)') "  If it looks like part of the namelist is missing, "
          write(stdout, '(a)') "  compare the number of characters read to the actual "
          write(stdout, '(a)') "  size of your file ($ wc -c pop2_in) and increase "
          write(stdout, '(a)') "  the buffer size if necessary."
       else
          write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
               "opening namelist file : ", trim(nml_filename)
       end if
       close(nml_in)
    end if

    call broadcast_scalar(nml_error, master_task)
    if (.not. is_iostat_end(nml_error)) then
       ! NOTE(bja, 2015-01) assuming that eof is the only proper exit
       ! code from the read.
       call document(subname, 'ERROR reading pop namelist file into buffer.')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    ! broadcast the namelist string
    call broadcast_scalar(nl_str, master_task)

    ! process namelist string to store in nl_buffer
    call marbl_nl_split_string(nl_str, nl_buffer)

    ! now every process reads the namelists from the buffer
    ioerror_msg=''
    call ecosys_status_log%construct()

    ! ecosys_driver_nml
    tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_driver_nml',ecosys_status_log)
    if (ecosys_status_log%labort_marbl) then
      call ecosys_status_log%log_error_trace('marbl_namelist', subname)
      call print_marbl_log(ecosys_status_log, 1)
    end if
    read(tmp_nl_buffer, nml=ecosys_driver_nml, iostat=nml_error, iomsg=ioerror_msg)
    if (nml_error /= 0) then
       write(stdout, *) subname, ": process ", my_task, ": namelist read error: ", nml_error, " : ", ioerror_msg
       call exit_POP(sigAbort, 'ERROR reading ecosys_driver_nml from buffer.')
    end if

    if (my_task == master_task) then
       write(stdout,*)
       write(stdout,ndelim_fmt)
       write(stdout,*)
       write(stdout,*) ' ecosys_driver:'
       write(stdout,*)
       write(stdout,*) ' ecosys_driver_nml namelist settings:'
       write(stdout,*)
       write(stdout,ecosys_driver_nml)
       write(stdout,*)
       write(stdout,delim_fmt)
    endif

    !-----------------------------------------------------------------------
    !  timer init
    !-----------------------------------------------------------------------

    call get_timer(ecosys_interior_timer     , 'ECOSYS_INTERIOR'   , nblocks_clinic , distrb_clinic%nprocs)
    call get_timer(ecosys_set_sflux_timer    , 'ECOSYS_SET_SFLUX'  , 1              , distrb_clinic%nprocs)
    call get_timer(ecosys_comp_CO3terms_timer, 'comp_CO3terms'     , nblocks_clinic , distrb_clinic%nprocs)

    !--------------------------------------------------------------------
    !  Initialize module variable land mask
    !--------------------------------------------------------------------

    allocate(land_mask(nx_block, ny_block, nblocks_clinic) )
    if (lmarginal_seas) then
       land_mask = region_mask /= 0
    else
       land_mask = region_mask > 0
    endif

    !--------------------------------------------------------------------
    !  MARBL setup is 3 steps:
    !  1) Configure (set variables that affect tracer count / other parms)
    !  2) Initialize ("lock" config vars so the aren't changed during init
    !     or in the time loop; write config vars to log; set parameters)
    !  3) Complete setup ("lock" parmameters so they aren't changed during
    !     time loop; write parameters to log)
    !
    !  POP can set up saved state, tracers, and forcing fields after (3)
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !  Configure marbl 
    !--------------------------------------------------------------------

    do iblock=1, nblocks_clinic

       call marbl_instances(iblock)%config(lgcm_has_global_ops = .true.,      &
                                           gcm_nl_buffer = nl_buffer)
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, ")%config()"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

    end do

    !--------------------------------------------------------------------
    !  Initialize marbl 
    !--------------------------------------------------------------------

    do iblock=1, nblocks_clinic

       call marbl_instances(iblock)%init(                                     &
            gcm_num_levels = km,                                              & 
            gcm_num_PAR_subcols = mcog_nbins,                                 &
            gcm_num_elements_interior_forcing = 1,                            & 
            gcm_num_elements_surface_forcing = num_elements,                  &
            gcm_dz = dz,                                                      &
            gcm_zw = zw,                                                      &
            gcm_zt = zt,                                                      &
            gcm_nl_buffer = nl_buffer,                                        &
            marbl_tracer_cnt = marbl_actual_tracer_cnt)

       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, ")%init()"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

       ! Make sure MARBL tracer count lines up with what POP expects
       if (marbl_actual_tracer_cnt.ne.ecosys_tracer_cnt) then
         write(log_message,"(A,I0,A,I0)") 'MARBL is computing tendencies for ', &
                    marbl_actual_tracer_cnt, ' tracers, but POP is expecting ', &
                    ecosys_tracer_cnt
         call document(subname, log_message)
         call exit_POP(sigAbort, 'Stopping in ' // subname)
       end if

    end do

    !--------------------------------------------------------------------
    !  Complete marbl setup
    !--------------------------------------------------------------------

    do iblock=1, nblocks_clinic

       call marbl_instances(iblock)%complete_config_and_init()
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, ")%complete_init_and_config()"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

    end do

    !--------------------------------------------------------------------
    !  Initialize ecosys tracers and saved state
    !--------------------------------------------------------------------

    call ecosys_saved_state_setup(saved_state_surf,                    &
         marbl_instances(1)%surface_saved_state)
    call ecosys_saved_state_setup(saved_state_interior,                &
         marbl_instances(1)%interior_saved_state)

    ! Initialize tracer_d_module input argument (needed before reading
    ! tracers from restart file)
    do n = 1, ecosys_tracer_cnt
       tracer_d_module(n)%short_name       = marbl_instances(1)%tracer_metadata(n)%short_name      
       tracer_d_module(n)%long_name        = marbl_instances(1)%tracer_metadata(n)%long_name       
       tracer_d_module(n)%units            = marbl_instances(1)%tracer_metadata(n)%units           
       tracer_d_module(n)%tend_units       = marbl_instances(1)%tracer_metadata(n)%tend_units      
       tracer_d_module(n)%flux_units       = marbl_instances(1)%tracer_metadata(n)%flux_units      
       tracer_d_module(n)%lfull_depth_tavg = marbl_instances(1)%tracer_metadata(n)%lfull_depth_tavg
       tracer_d_module(n)%scale_factor     = c1
    end do

    ! pass ecosys_tracer_init_nml to
    ! ecosys_tracers_and_saved_state_init()
    tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_tracer_init_nml',       &
                                   ecosys_status_log)
    if (ecosys_status_log%labort_marbl) then
      call ecosys_status_log%log_error_trace('marbl_namelist', subname)
      call print_marbl_log(ecosys_status_log, 1)
    end if

    call ecosys_tracers_and_saved_state_init(                    &
       ciso_on,                                                  &
       init_ts_file_fmt,                                         &
       read_restart_filename,                                    &
       tracer_d_module(:),                                       &
       marbl_instances(1)%tracer_metadata(:)%tracer_module_name, &
       tmp_nl_buffer,                                            &
       land_mask,                                                &
       tracer_module(:,:,:,:,:,:),                               &
       ecosys_restart_filename,                                  &
       errorCode)       

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_driver_init')
       return
    endif

    !--------------------------------------------------------------------
    !  Initialize ecosys forcing fields
    !--------------------------------------------------------------------

    ! Determine MARBL tracer indices for tracers that use virtual fluxes
    ! (must be done prior to call to ecosys_forcing_init, when vflux is set)
    dic_ind   = marbl_instances(1)%get_tracer_index('DIC')
    call document(subname, 'dic_ind', dic_ind)
    alk_ind   = marbl_instances(1)%get_tracer_index('ALK')
    call document(subname, 'alk_ind', alk_ind)
    dic_alt_co2_ind = marbl_instances(1)%get_tracer_index('DIC_ALT_CO2')
    call document(subname, 'dic_alt_co2_ind', dic_alt_co2_ind)
    di13c_ind = marbl_instances(1)%get_tracer_index('DI13C')
    call document(subname, 'di13c_ind', di13c_ind)
    di14c_ind = marbl_instances(1)%get_tracer_index('DI14C')
    call document(subname, 'di14c_ind', di14c_ind)

    ! forcing module requires two MARBL parameter values (set during init)
    call marbl_instances(1)%parameters%get('iron_frac_in_dust', fe_frac_dust, &
                                           ecosys_status_log)
    if (ecosys_status_log%labort_marbl) then
      call ecosys_status_log%log_error_trace('parameters%get(iron_frac_in_dust)', subname)
      call print_marbl_log(ecosys_status_log, 1)
    end if

    call marbl_instances(1)%parameters%get('iron_frac_in_bc', fe_frac_bc,     &
                                           ecosys_status_log)
    if (ecosys_status_log%labort_marbl) then
      call ecosys_status_log%log_error_trace('parameters%get(iron_frac_in_bc)', subname)
      call print_marbl_log(ecosys_status_log, 1)
    end if

    ! pass ecosys_forcing_data_nml
    ! to ecosys_forcing_init()
    ! Also pass marbl_instance%surface_forcing_metadata
    tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_forcing_data_nml',      &
                                   ecosys_status_log)
    if (ecosys_status_log%labort_marbl) then
      call ecosys_status_log%log_error_trace('marbl_namelist', subname)
      call print_marbl_log(ecosys_status_log, 1)
    end if

    call ecosys_forcing_init(ciso_on,                                         &
                             num_elements,                                    &
                             land_mask,                                       &
                             fe_frac_dust,                                    &
                             fe_frac_bc,                                      &
                             marbl_instances(1)%surface_input_forcings,       &
                             marbl_instances(1)%interior_input_forcings,      &
                             tmp_nl_buffer)

    !--------------------------------------------------------------------
    !  Initialize ecosys_driver module variables
    !--------------------------------------------------------------------

    associate(diag_cnt => marbl_instances(1)%surface_forcing_diags%diag_cnt)
    allocate(surface_forcing_diags(nx_block, ny_block, diag_cnt, nblocks_clinic))
    end associate

    tadvect_ctype(1:ecosys_tracer_cnt) = ecosys_tadvect_ctype

    !--------------------------------------------------------------------
    ! Initialize tavg ids (need only do this using first block)
    !--------------------------------------------------------------------

    call ecosys_tavg_init(marbl_instances(1))

    !--------------------------------------------------------------------
    ! Register and set Chl field for short-wave absorption
    ! Register and set fco2, the  air-sea co2 gas flux
    !--------------------------------------------------------------------

    do iblock=1, nblocks_clinic
       ! Register flux_co2 with MARBL surface forcing outputs
       call marbl_instances(iblock)%surface_forcing_output%add_sfo(           &
              num_elements = num_elements,                                    &
              field_name   = "flux_co2",                                      &
              sfo_id       = flux_co2_id,                                     &
              marbl_status_log = marbl_instances(iblock)%StatusLog)
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, &
                                     ")%surface_forcing_output%add_sfo(flux_co2)"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

       ! Register totalChl with MARBL surface forcing outputs
       call marbl_instances(iblock)%surface_forcing_output%add_sfo(           &
              num_elements = num_elements,                                    &
              field_name   = "totalChl",                                      &
              sfo_id       = totalChl_id,                                     &
              marbl_status_log = marbl_instances(iblock)%StatusLog)
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, &
                                     ")%surface_forcing_output%add_sfo(totalChl)"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()
    end do
    call named_field_register('SFLUX_CO2'        , sflux_co2_nf_ind)
    call named_field_register('model_chlorophyll', totChl_surf_nf_ind)

    !--------------------------------------------------------------------
    ! allocate space for fields for which global averages are to be computed
    !--------------------------------------------------------------------

    glo_avg_field_cnt = size(marbl_instances(1)%glo_avg_fields_interior, dim=1)
    allocate(glo_avg_fields_interior(nx_block, ny_block, nblocks_clinic, glo_avg_field_cnt))

    glo_avg_field_cnt = size(marbl_instances(1)%glo_avg_fields_surface, dim=2)
    allocate(glo_avg_fields_surface(nx_block, ny_block, nblocks_clinic, glo_avg_field_cnt))

    ! initialize to zero so that values not set at runtime don't cause problems in global sum function
    glo_avg_fields_interior(:,:,:,:) = c0
    glo_avg_fields_surface(:,:,:,:) = c0

    if ((size(glo_avg_fields_interior, dim=4) /= 0) .or. (size(glo_avg_fields_surface, dim=4) /= 0)) then
       allocate(glo_avg_area_masked(nx_block, ny_block, nblocks_clinic))
       where (land_mask(:,:,:))
          glo_avg_area_masked(:,:,:) = TAREA(:,:,:)
       else where
          glo_avg_area_masked(:,:,:) = c0
       end where

       glo_avg_norm_fact = c1 / global_sum(glo_avg_area_masked(:,:,:), distrb_clinic, field_loc_center)
    end if

    ! FIXME : move the setup of running means of global averages into MARBL

    call ecosys_driver_init_rmean_var(marbl_instances(1)%glo_avg_rmean_interior, &
                                      ecosys_restart_filename, glo_avg_rmean_ind_interior)

    call ecosys_driver_init_rmean_var(marbl_instances(1)%glo_avg_rmean_surface, &
                                      ecosys_restart_filename, glo_avg_rmean_ind_surface)

    call ecosys_driver_init_rmean_var(marbl_instances(1)%glo_scalar_rmean_interior, &
                                      ecosys_restart_filename, glo_scalar_rmean_ind_interior)

    call ecosys_driver_init_rmean_var(marbl_instances(1)%glo_scalar_rmean_surface, &
                                      ecosys_restart_filename, glo_scalar_rmean_ind_surface)

    ! copy values from POP's running mean to MARBL interface

    do n = 1, size(glo_avg_rmean_ind_interior(:))
       call running_mean_get_var(glo_avg_rmean_ind_interior(n), vals_0d=rmean_val)
       do iblock = 1, nblocks_clinic
          marbl_instances(iblock)%glo_avg_rmean_interior(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_avg_rmean_ind_surface(:))
       call running_mean_get_var(glo_avg_rmean_ind_surface(n), vals_0d=rmean_val)
       do iblock = 1, nblocks_clinic
          marbl_instances(iblock)%glo_avg_rmean_surface(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_scalar_rmean_ind_interior(:))
       call running_mean_get_var(glo_scalar_rmean_ind_interior(n), vals_0d=rmean_val)
       do iblock = 1, nblocks_clinic
          marbl_instances(iblock)%glo_scalar_rmean_interior(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_scalar_rmean_ind_surface(:))
       call running_mean_get_var(glo_scalar_rmean_ind_surface(n), vals_0d=rmean_val)
       do iblock = 1, nblocks_clinic
          marbl_instances(iblock)%glo_scalar_rmean_surface(n)%rmean = rmean_val
       end do
    end do

  end subroutine ecosys_driver_init

  !-----------------------------------------------------------------------

  subroutine ecosys_driver_init_rmean_var(marbl_running_mean_var, ecosys_restart_filename, rmean_ind)

    use marbl_interface_types, only : marbl_running_mean_0d_type
    use running_mean_mod     , only : running_mean_define_var
    use running_mean_mod     , only : running_mean_init_var

    type(marbl_running_mean_0d_type), intent(in)  :: marbl_running_mean_var(:)
    character(char_len)             , intent(in)  :: ecosys_restart_filename
    integer (int_kind), allocatable , intent(out) :: rmean_ind(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: rmean_var_cnt
    integer (int_kind) :: n

    !-----------------------------------------------------------------------

    rmean_var_cnt = size(marbl_running_mean_var, 1)
    allocate(rmean_ind(rmean_var_cnt))

    do n = 1, rmean_var_cnt
       call running_mean_define_var(name=trim(marbl_running_mean_var(n)%sname), rank=0, &
          timescale=marbl_running_mean_var(n)%timescale, index=rmean_ind(n))

       if (marbl_running_mean_var(n)%linit_by_val) then
          call running_mean_init_var(rmean_ind(n), vals_0d=marbl_running_mean_var(n)%init_val)
       else
          call running_mean_init_var(rmean_ind(n), filename=ecosys_restart_filename)
       end if
    end do

  end subroutine ecosys_driver_init_rmean_var

  !***********************************************************************

  subroutine ecosys_driver_set_interior(     &
       FRACR_BIN, QSW_RAW_BIN, QSW_BIN,      &
       TEMP_OLD, TEMP_CUR,                   &
       SALT_OLD, SALT_CUR,                   &
       TRACER_MODULE_OLD, TRACER_MODULE_CUR, &
       DTRACER_MODULE, this_block)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute source-sink terms
    !  accumulate commnon tavg fields related to source-sink terms

    use constants          , only : salt_to_ppt
    use grid               , only : KMT
    use grid               , only : DZT
    use grid               , only : partial_bottom_cells
    use mcog               , only : mcog_nbins
    use state_mod          , only : ref_pressure
    use ecosys_forcing_mod , only : interior_forcing_fields
    use ecosys_forcing_mod , only : ecosys_forcing_set_interior_time_varying_forcing_data
    use ecosys_forcing_mod , only : dustflux_ind
    use ecosys_forcing_mod , only : PAR_col_frac_ind
    use ecosys_forcing_mod , only : surf_shortwave_ind
    use ecosys_forcing_mod , only : temperature_ind
    use ecosys_forcing_mod , only : salinity_ind
    use ecosys_forcing_mod , only : pressure_ind

    implicit none

    real (r8), dimension(nx_block, ny_block, mcog_nbins) , intent(in)    :: FRACR_BIN         ! fraction of cell occupied by mcog bin
    real (r8), dimension(nx_block, ny_block, mcog_nbins) , intent(in)    :: QSW_RAW_BIN       ! raw (directly from cpl) shortwave into each mcog column (W/m^2)
    real (r8), dimension(nx_block, ny_block, mcog_nbins) , intent(in)    :: QSW_BIN           ! shortwave into each mcog bin, potentially modified by coszen factor (W/m^2)
    real (r8), dimension(nx_block, ny_block, km)         , intent(in)    :: TEMP_OLD          ! old potential temperature (C)
    real (r8), dimension(nx_block, ny_block, km)         , intent(in)    :: TEMP_CUR          ! current potential temperature (C)
    real (r8), dimension(nx_block, ny_block, km)         , intent(in)    :: SALT_OLD          ! old salinity (msu)
    real (r8), dimension(nx_block, ny_block, km)         , intent(in)    :: SALT_CUR          ! current salinity (msu)
    real (r8), dimension(:,:,:,:)                        , intent(in)    :: TRACER_MODULE_OLD ! old tracer values
    real (r8), dimension(:,:,:,:)                        , intent(in)    :: TRACER_MODULE_CUR ! current tracer values
    type (block)                                         , intent(in)    :: this_block        ! block information for this block
    real (r8), dimension(:,:,:,:)                        , intent(inout) :: DTRACER_MODULE    ! computed source/sink terms

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_driver:ecosys_driver_set_interior'
    character(len=char_len) :: log_message
    integer (int_kind)      :: i   ! nx_block loop index
    integer (int_kind)      :: c   ! ny_block / column loop index
    integer (int_kind)      :: k   ! vertical level index
    integer (int_kind)      :: bid ! local block address for this block
    integer (int_kind)      :: n, d, ncols
    real (r8), dimension(nx_block, ny_block, km) :: temperature, salinity
    real (r8), dimension(nx_block, ny_block, km) :: pressure
    !-----------------------------------------------------------------------

    bid = this_block%local_id

    ! NTOE: gcm dependent quantities (i.e. time stepping). need to be
    ! averaged into a single quantity for marbl ecosys

    call timer_start(ecosys_interior_timer, block_id=bid)

    associate(marbl_interior_forcings => marbl_instances(bid)%interior_input_forcings)

    !-----------------------------------------------------------------------
    ! Set input surface forcing data and surface saved state data
    !-----------------------------------------------------------------------

    temperature = p5*(temp_old + temp_cur)
    salinity = p5*(salt_old + salt_cur)*salt_to_ppt
    do k=1,km
      ! NOTE: ref_pressure is a function, not an array
       pressure(:,:,k) = ref_pressure(k)
    end do
    call ecosys_forcing_set_interior_time_varying_forcing_data(FRACR_bin,     &
                            QSW_RAW_BIN, QSW_BIN, temperature, salinity,      &
                            pressure, ecosys_qsw_distrb_const, bid)

    do c = this_block%jb,this_block%je
       do i = this_block%ib,this_block%ie

          if (land_mask(i,c,bid)) then

             !-----------------------------------------------------------
             ! Copy data from slab to column
             !-----------------------------------------------------------

             ! --- set marbl_domain kmt and if partial bottom cells then also delta_z ---

             marbl_instances(bid)%domain%kmt = KMT(i, c, bid)
             if (partial_bottom_cells) then
                marbl_instances(bid)%domain%delta_z(:) = DZT(i, c, :, bid)
             end if

             ! --- set forcing fields ---

             do n = 1, num_interior_forcing_fields
               if (allocated(interior_forcing_fields(n)%field_0d)) then
                 marbl_instances(bid)%interior_input_forcings(n)%field_0d(1) = &
                      interior_forcing_fields(n)%field_0d(i,c,bid)
               else
                 marbl_instances(bid)%interior_input_forcings(n)%field_1d(1,:) = &
                      interior_forcing_fields(n)%field_1d(i,c,:,bid)
               end if
             end do

             ! --- set column tracers ---

             do n = 1, ecosys_tracer_cnt
                marbl_instances(bid)%column_tracers(n, :) = p5*(tracer_module_old(i, c, :, n) + tracer_module_cur(i, c, :, n))
             end do 

             ! --- copy data from slab to column for marbl_saved_state ---
             do n=1,size(saved_state_interior)
               marbl_instances(bid)%interior_saved_state%state(n)%field_3d(:,1) = &
                 saved_state_interior(n)%field_3d(i,c,:,bid)
             end do

             !-----------------------------------------------------------
             !  compute time derivatives for ecosystem state variables
             !-----------------------------------------------------------

             call marbl_instances(bid)%set_interior_forcing()
             if (marbl_instances(bid)%StatusLog%labort_marbl) then
                write(log_message,"(A,I0,A)") "marbl_instances(", bid, &
                                              ")%set_surface_forcing()"
                call marbl_instances(bid)%StatusLog%log_error_trace(log_message, subname)
             end if
             call print_marbl_log(marbl_instances(bid)%StatusLog, bid)
             call marbl_instances(bid)%StatusLog%erase()

             !-----------------------------------------------------------
             ! copy marbl column data back to slab
             !-----------------------------------------------------------

             do n=1,size(saved_state_interior)
               saved_state_interior(n)%field_3d(i,c,:,bid) =               &
                 marbl_instances(bid)%interior_saved_state%state(n)%field_3d(:,1)
             end do
                
             do n = 1, ecosys_tracer_cnt
                dtracer_module(i, c, :, n) = marbl_instances(bid)%column_dtracers(n, :)
             end do

             ! copy values to be used in computing requested global averages
             ! arrays have zero extent if none are requested
             glo_avg_fields_interior(i, c, bid, :) = marbl_instances(bid)%glo_avg_fields_interior(:)

          !-----------------------------------------------------------
          ! Update pop tavg diags
          !-----------------------------------------------------------

             call ecosys_tavg_accumulate((/i/), (/c/), bid,                                   &
                  marbl_interior_forcing_diags = marbl_instances(bid)%interior_forcing_diags, &
                  marbl_interior_restore_diags = marbl_instances(bid)%interior_restore_diags)
 
          end if ! end if land_mask > 0

       end do ! do i
    end do ! do c
    
    end associate

    call timer_stop(ecosys_interior_timer, block_id=bid)

  end subroutine ecosys_driver_set_interior

  !***********************************************************************

  subroutine ecosys_driver_set_sflux( &
       u10_sqr,                       &
       ifrac,                         &
       press,                         &
       dust_flux,                     &
       black_carbon_flux,             &
       sst,                           &
       sss,                           &
       surface_vals_old,              &
       surface_vals_cur,              &
       stf_module)

    ! DESCRIPTION: 
    ! Sets surface input forcing data and calls marbl to
    ! compute surface tracer fluxes

    use named_field_mod      , only : named_field_set
    use time_management      , only : check_time_flag
    use domain               , only : nblocks_clinic
    use ecosys_forcing_mod   , only : ecosys_forcing_set_surface_time_varying_forcing_data
    use ecosys_forcing_mod   , only : surface_forcing_fields

    implicit none

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: u10_sqr           ! 10m wind speed squared (cm/s)**2
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: ifrac             ! sea ice fraction (non-dimensional)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: press             ! sea level atmospheric pressure (dyne/cm**2)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: dust_flux         ! dust flux (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: black_carbon_flux ! black carbon flux (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: sst               ! sea surface temperature (c)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: sss               ! sea surface salinity (psu)
    real (r8), dimension(:,:,:,:)                             , intent(in)    :: surface_vals_old
    real (r8), dimension(:,:,:,:)                             , intent(in)    :: surface_vals_cur  ! module tracers
    real (r8), dimension(:,:,:,:)                             , intent(inout) :: stf_module

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_driver:ecosys_driver_set_sflux'
    character(char_len) :: log_message
    integer (int_kind) :: index_marbl                                 ! marbl index
    integer (int_kind) :: i, j, iblock, n                             ! pop loop indices
    integer (int_kind) :: glo_scalar_cnt
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Set input surface forcing data and surface saved state data
    !-----------------------------------------------------------------------

    call ecosys_forcing_set_surface_time_varying_forcing_data( &
         ciso_on,                         &
         land_mask,                       &
         u10_sqr,                         &
         ifrac,                           &
         press,                           &
         dust_flux,                       &
         black_carbon_flux,               &
         sst,                             &
         sss)

    !-----------------------------------------------------------------------
    ! Set output surface tracer fluxes
    !-----------------------------------------------------------------------

    call timer_start(ecosys_set_sflux_timer)

    call ecosys_driver_set_global_scalars('surface')

    ! FIXME : add OMP directive to this loop
    do iblock = 1, nblocks_clinic

       !-----------------------------------------------------------------------
       ! Copy data from slab data structure to column input for marbl
       !-----------------------------------------------------------------------

       do j = 1, ny_block
          do i = 1,nx_block
             index_marbl = i + (j-1)*nx_block

             do n = 1,num_surface_forcing_fields
                marbl_instances(iblock)%surface_input_forcings(n)%field_0d(index_marbl) = &
                     surface_forcing_fields(n)%field_0d(i,j,iblock)
             end do

             do n = 1,ecosys_tracer_cnt
                marbl_instances(iblock)%surface_vals(index_marbl,n) = &
                     p5*(surface_vals_old(i,j,n,iblock) + surface_vals_cur(i,j,n,iblock))

             end do

             do n=1,size(saved_state_surf)
               marbl_instances(iblock)%surface_saved_state%state(n)%field_2d(index_marbl) = &
                 saved_state_surf(n)%field_2d(i,j,iblock)
             end do

          end do
       end do

       !-----------------------------------------------------------------------
       ! Determine surface forcing flux - marbl
       !-----------------------------------------------------------------------

       call marbl_instances(iblock)%set_surface_forcing()

       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
          write(log_message,"(A,I0,A)") "marbl_instances(", iblock, &
                                        ")%set_surface_forcing()"
          call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

       !-----------------------------------------------------------------------
       ! Copy data from marbl output column to pop slab data structure
       !-----------------------------------------------------------------------

       do j = 1, ny_block
          do i = 1,nx_block
             index_marbl = i + (j-1)*nx_block

             do n=1,size(saved_state_surf)
               saved_state_surf(n)%field_2d(i,j,iblock) = &
                 marbl_instances(iblock)%surface_saved_state%state(n)%field_2d(index_marbl)
             end do

             do n=1,2
               surface_forcing_outputs(i,j,iblock,n) = &
                  marbl_instances(iblock)%surface_forcing_output%sfo(n)%forcing_field(index_marbl)
             end do
             
             do n = 1,ecosys_tracer_cnt
                stf_module(i,j,n,iblock) = &
                     marbl_instances(iblock)%surface_tracer_fluxes(index_marbl,n)  
             end do

             do n=1,marbl_instances(1)%surface_forcing_diags%diag_cnt
                surface_forcing_diags(i,j,n,iblock) = &
                     marbl_instances(iblock)%surface_forcing_diags%diags(n)%field_2d(index_marbl)
             end do

             ! copy values to be used in computing requested global averages
             ! arrays have zero extent if none are requested
             glo_avg_fields_surface(i,j,iblock,:) = marbl_instances(iblock)%glo_avg_fields_surface(index_marbl,:)
          enddo
       end do

    enddo ! end loop over iblock

    call ecosys_driver_comp_global_averages('surface')

    call timer_stop(ecosys_set_sflux_timer)

    !-----------------------------------------------------------------------
    ! Update named fields
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       call named_field_set(totChl_surf_nf_ind, iblock, surface_forcing_outputs(:, :, iblock, totalChl_id))
       
       !  set air-sea co2 gas flux named field, converting units from
       !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
       if (lflux_gas_co2) then
          call named_field_set(sflux_co2_nf_ind, iblock, 44.0e-8_r8 * surface_forcing_outputs(:, :, iblock, flux_co2_id))
       end if
    end do
    
    ! Note: out of this subroutine rest of pop needs stf_module, ph_prev_surf, named_state, diagnostics

  end subroutine ecosys_driver_set_sflux

  !***********************************************************************

  subroutine ecosys_driver_comp_global_averages(field_source)

    ! DESCRIPTION: 
    ! perform global operations

    use global_reductions, only : global_sum_prod
    use constants        , only : field_loc_center
    use running_mean_mod , only : running_mean_update_var
    use running_mean_mod , only : running_mean_get_var

    character (*), intent(in) :: field_source ! 'interior' or 'surface'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real (r8), pointer          :: glo_avg_fields(:,:,:,:)
    integer (int_kind), pointer :: glo_avg_rmean_ind(:)
    real (r8), allocatable      :: glo_avg(:)
    real (r8), allocatable      :: glo_avg_rmean(:)
    integer (int_kind)          :: glo_avg_field_cnt
    integer (int_kind)          :: n, iblock
    !-----------------------------------------------------------------------

    if (field_source .eq. 'interior') then
       glo_avg_fields    => glo_avg_fields_interior(:,:,:,:)
       glo_avg_rmean_ind => glo_avg_rmean_ind_interior(:)
    else
       glo_avg_fields    => glo_avg_fields_surface(:,:,:,:)
       glo_avg_rmean_ind => glo_avg_rmean_ind_surface(:)
    end if

    glo_avg_field_cnt = size(glo_avg_fields, dim=4)

    if (glo_avg_field_cnt /= 0) then
       allocate(glo_avg(glo_avg_field_cnt))
       allocate(glo_avg_rmean(glo_avg_field_cnt))

       ! compute global means, and update their running means
       do n = 1, glo_avg_field_cnt
          glo_avg(n) = glo_avg_norm_fact * global_sum_prod(glo_avg_fields(:,:,:,n), &
             glo_avg_area_masked(:,:,:), distrb_clinic, field_loc_center)
          call running_mean_update_var(glo_avg_rmean_ind(n), vals_0d=glo_avg(n))
          call running_mean_get_var(glo_avg_rmean_ind(n), vals_0d=glo_avg_rmean(n))
       end do

       ! store global means, and their running means, into appropriate component of marbl_instances
       if (field_source .eq. 'interior') then
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_avg_averages_interior(:)    = glo_avg(:)
             marbl_instances(iblock)%glo_avg_rmean_interior(:)%rmean = glo_avg_rmean(:)
          end do
       else
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_avg_averages_surface(:)    = glo_avg(:)
             marbl_instances(iblock)%glo_avg_rmean_surface(:)%rmean = glo_avg_rmean(:)
          end do
       end if

       deallocate(glo_avg_rmean)
       deallocate(glo_avg)
    end if

  end subroutine ecosys_driver_comp_global_averages

  !***********************************************************************

  subroutine ecosys_driver_set_global_scalars(field_source)

    character (*), intent(in) :: field_source ! 'interior' or 'surface'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind)          :: iblock

    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       call marbl_instances(iblock)%set_global_scalars(field_source)
    end do

    call ecosys_driver_update_scalar_rmeans(field_source)

  end subroutine ecosys_driver_set_global_scalars

  !***********************************************************************

  subroutine ecosys_driver_update_scalar_rmeans(field_source)

    ! DESCRIPTION: 
    ! update running means of scalar variables

    use running_mean_mod , only : running_mean_update_var
    use running_mean_mod , only : running_mean_get_var
    use io_types         , only : stdout

    character (*), intent(in) :: field_source ! 'interior' or 'surface'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname   = 'ecosys_driver:ecosys_driver_update_scalar_rmeans'
    character(*), parameter :: fmt_str   = '(A,1X,A)'
    character(*), parameter :: fmt_str_i = '(A,1X,A,1X,I0)'
    character(*), parameter :: fmt_str_e = '(A,1X,A,1X,E23.16)'
    real (r8)               :: rmean_val
    integer (int_kind)      :: n, iblock
    !-----------------------------------------------------------------------

    if (field_source .eq. 'interior') then
       do n = 1, size(glo_scalar_rmean_ind_interior(:))
          ! verify that all instances have same value of glo_scalar_interior
          do iblock = 2, nblocks_clinic
             if (marbl_instances(iblock)%glo_scalar_interior(n) /= marbl_instances(1)%glo_scalar_interior(n)) then
                write(stdout, fmt_str)   subname, 'mismatch in glo_scalar_interior values across MARBL instances'
                write(stdout, fmt_str_i) subname, 'rmean index', n
                write(stdout, fmt_str_e) subname, 'iblock 1 value', marbl_instances(1)%glo_scalar_interior(n)
                write(stdout, fmt_str_i) subname, 'iblock', iblock
                write(stdout, fmt_str_e) subname, 'mismatched value', marbl_instances(iblock)%glo_scalar_interior(n)
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             end if
          end do

          call running_mean_update_var(glo_scalar_rmean_ind_interior(n), vals_0d=marbl_instances(1)%glo_scalar_interior(n))
          call running_mean_get_var(glo_scalar_rmean_ind_interior(n), vals_0d=rmean_val)
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_scalar_rmean_interior(n)%rmean = rmean_val
          end do
       end do
    else
       do n = 1, size(glo_scalar_rmean_ind_surface(:))
          ! verify that all instances have same value of glo_scalar_surface
          do iblock = 2, nblocks_clinic
             if (marbl_instances(iblock)%glo_scalar_surface(n) /= marbl_instances(1)%glo_scalar_surface(n)) then
                write(stdout, fmt_str)   subname, 'mismatch in glo_scalar_surface values across MARBL instances'
                write(stdout, fmt_str_i) subname, 'rmean index', n
                write(stdout, fmt_str_e) subname, 'iblock 1 value', marbl_instances(1)%glo_scalar_surface(n)
                write(stdout, fmt_str_i) subname, 'iblock', iblock
                write(stdout, fmt_str_e) subname, 'mismatched value', marbl_instances(iblock)%glo_scalar_surface(n)
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             end if
          end do

          call running_mean_update_var(glo_scalar_rmean_ind_surface(n), vals_0d=marbl_instances(1)%glo_scalar_surface(n))
          call running_mean_get_var(glo_scalar_rmean_ind_surface(n), vals_0d=rmean_val)
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_scalar_rmean_surface(n)%rmean = rmean_val
          end do
       end do
    end if

  end subroutine ecosys_driver_update_scalar_rmeans

  !***********************************************************************

  subroutine ecosys_driver_tavg_forcing()

    ! !DESCRIPTION:
    !  accumulate common tavg fields for tracer surface fluxes

    implicit none 

    call ecosys_tavg_accumulate_flux(surface_forcing_diags, marbl_instances(:))

  end subroutine ecosys_driver_tavg_forcing

  !***********************************************************************

  subroutine ecosys_driver_unpack_source_sink_terms(source, destination)

    real(r8), dimension(:, :, :), intent(in)  :: source
    real(r8), dimension(:, :, :), intent(out) :: destination

    destination(:, :, :) = source(:, :, :)

  end subroutine ecosys_driver_unpack_source_sink_terms

  !*****************************************************************************

  subroutine ecosys_driver_write_restart(restart_file, action)

    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files

    use io_types, only : io_field_desc
    use io      , only : data_set
    use io      , only : datafile

    implicit none

    character(*)                 , intent(in)    :: action
    type (datafile)              , intent(inout) :: restart_file

    type (io_field_desc), dimension(:), allocatable, save :: surf_iodesc
    type (io_field_desc), dimension(:), allocatable, save :: col_iodesc
    integer :: n
    !-----------------------------------------------------------------------

    if (trim(action) == 'define') then

       allocate(surf_iodesc(size(saved_state_surf)))
       surf_iodesc = ecosys_saved_state_construct_io_fields(restart_file,     &
            saved_state_surf, size(saved_state_surf))

       allocate(col_iodesc(size(saved_state_interior)))
       col_iodesc = ecosys_saved_state_construct_io_fields(restart_file,      &
           saved_state_interior, size(saved_state_interior))

    else if (trim(action) == 'write') then

       do n=1,size(saved_state_surf)
          call data_set (restart_file, 'write', surf_iodesc(n))
       end do
       do n=1,size(saved_state_interior)
          call data_set (restart_file, 'write', col_iodesc(n))
       end do
       deallocate(surf_iodesc)
       deallocate(col_iodesc)

    endif

  end subroutine ecosys_driver_write_restart

  !*****************************************************************************

  subroutine print_marbl_log(log_to_print, iblock)

    use marbl_logging             , only : marbl_status_log_entry_type

    class(marbl_log_type), intent(in) :: log_to_print
    integer,               intent(in) :: iblock

    character(len=*), parameter :: subname = 'ecosys_driver:print_marbl_log'
    type(marbl_status_log_entry_type), pointer :: tmp
    logical :: iam_master

    ! FIXME (02-2016,mnl): need better logic on which items to print
    iam_master = (my_task.eq.master_task).and.(iblock.eq.1)
    tmp => log_to_print%FullLog
    do while (associated(tmp))
      if (tmp%lall_tasks.or.iam_master) then
        if (tmp%lall_tasks.and.(.not.iam_master)) then
          write(stdout, "(A,I0,A,I0,2A)") "(Task ", my_task, ', block ', iblock, &
                '): ', trim(tmp%LogMessage)
        else
          write(stdout, "(A)") trim(tmp%LogMessage)
        end if
      end if
      tmp => tmp%next
    end do

    if (log_to_print%labort_marbl) then
      call document(subname, 'ERROR reported from MARBL library')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

  end subroutine print_marbl_log

  !*****************************************************************************

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

