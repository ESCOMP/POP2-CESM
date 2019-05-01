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

  use shr_infnan_mod            , only : shr_infnan_isnan

  use marbl_interface           , only : marbl_interface_class

  use namelist_from_str_mod     , only : namelist_split_by_line
  use namelist_from_str_mod     , only : namelist_split_by_nl
  use namelist_from_str_mod     , only : namelist_find

  use ecosys_tavg               , only : ecosys_tavg_init
  use ecosys_tavg               , only : ecosys_tavg_accumulate_interior
  use ecosys_tavg               , only : ecosys_tavg_accumulate_surface
  use ecosys_tavg               , only : ecosys_tavg_accumulate_scalar_rmeans

  use passive_tracer_tools      , only : tracer_read

  use timers                    , only : timer_start
  use timers                    , only : timer_stop
  use timers                    , only : get_timer

  use ecosys_tracers_and_saved_state_mod, only : ecosys_tracers_and_saved_state_init
  use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_setup
  use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_type
  use ecosys_tracers_and_saved_state_mod, only : surface_flux_saved_state
  use ecosys_tracers_and_saved_state_mod, only : interior_tendency_saved_state
  use ecosys_tracers_and_saved_state_mod, only : dic_ind, alk_ind, dic_alt_co2_ind, alk_alt_co2_ind
  use ecosys_tracers_and_saved_state_mod, only : di13c_ind, di14c_ind
  use ecosys_tracers_and_saved_state_mod, only : o2_ind, no3_ind, po4_ind, don_ind, donr_ind, dop_ind, dopr_ind
  use ecosys_tracers_and_saved_state_mod, only : sio3_ind, fe_ind, doc_ind, docr_ind, do13ctot_ind, do14ctot_ind

  ! Provide marbl_tracer_cnt to passive_tracers.F90 (as ecosys_tracer_cnt)
  use ecosys_tracers_and_saved_state_mod, only : ecosys_tracer_cnt => marbl_tracer_cnt
  ! Provide ecosys_forcing_tracer_ref_val to passive_tracers.F90
  ! (as ecosys_driver_tracer_ref_val)
  use ecosys_forcing_mod, only : ecosys_driver_tracer_ref_val => ecosys_forcing_tracer_ref_val

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_driver_init
  public :: ecosys_driver_set_interior_forcing
  public :: ecosys_driver_set_interior
  public :: ecosys_driver_set_global_scalars
  public :: ecosys_driver_comp_global_averages
  public :: ecosys_driver_set_sflux_forcing
  public :: ecosys_driver_set_sflux
  public :: ecosys_driver_post_set_sflux
  public :: ecosys_driver_tavg_forcing
  public :: ecosys_driver_write_restart
  public :: ecosys_driver_unpack_source_sink_terms
  public :: ecosys_driver_print_marbl_timers

  ! These are public so passive_tracers.F90 only uses ecosys_driver
  public :: ecosys_driver_tracer_ref_val
  public :: ecosys_tracer_cnt

  private :: ecosys_driver_update_scalar_rmeans

  !*****************************************************************************

  !-----------------------------------------------------------------------
  ! timers
  !-----------------------------------------------------------------------

  integer (int_kind) :: ecosys_interior_forcing_timer
  integer (int_kind) :: ecosys_interior_pop_to_marbl
  integer (int_kind) :: ecosys_interior_marbl_to_pop
  integer (int_kind) :: ecosys_interior_marbl_tavg
  integer (int_kind) :: ecosys_interior_timer
  integer (int_kind), target :: ecosys_surface_global_sum_timer
  integer (int_kind), target :: ecosys_interior_global_sum_timer
  integer (int_kind) :: ecosys_set_sflux_timer

  !-----------------------------------------------------------------------
  ! module variables
  !-----------------------------------------------------------------------

  type(marbl_interface_class) :: marbl_instances(max_blocks_clinic)

  integer (int_kind)  :: totChl_surf_nf_ind = 0 ! total chlorophyll in surface layer
  integer (int_kind)  :: sflux_co2_nf_ind   = 0 ! air-sea co2 gas flux

  character (char_len)               :: ecosys_tadvect_ctype    ! advection method for ecosys tracers
  logical   (log_kind) , public      :: ecosys_qsw_distrb_const
  logical   (log_kind)               :: ciso_on
  logical   (log_kind) , allocatable :: land_mask(:, :, :)
  real      (r8)       , allocatable :: surface_flux_diags(:, :, :, :)

  integer   (int_kind)               :: sfo_cnt
  integer   (int_kind)               :: flux_co2_id
  integer   (int_kind)               :: totalChl_id
  real      (r8)       , allocatable :: surface_flux_outputs(:, :, :, :)


  ! Variables related to global averages

  real (r8), allocatable, target, dimension(:,:,:,:) :: glo_avg_fields_interior
  real (r8), allocatable, target, dimension(:,:,:,:) :: glo_avg_fields_surface
  real (r8), allocatable, dimension(:,:,:)           :: glo_avg_area_masked
  real (r8)                                          :: glo_avg_norm_fact

  ! Variables used to track POP index of columns passed to MARBL
  ! (surface fluxes are computed for multiple columns simultaneously)
  integer, allocatable, dimension(:)   :: marbl_col_cnt                          ! nblocks_clinic
  integer, allocatable, dimension(:,:) :: marbl_col_to_pop_i, marbl_col_to_pop_j ! max(marbl_col_cnt) x nblocks_clinic

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_driver_init( &
       ecosys_driver_ind_begin,  &
       init_ts_file_fmt,         &
       read_restart_filename,    &
       tracer_d_module,          &
       tracer_module,            &
       tadvect_ctype,            &
       lhas_riv_flux,            &
       errorCode)

    ! !DESCRIPTION:
    !  Initialize ecosys_driver passive tracers. This involves:
    !  1) setting ecosys and ecosys_ciso module index bounds
    !  2) calling ecosys and ecosys_ciso module init subroutine

    use grid, only : dz
    use grid, only : zt
    use grid, only : zw
    use grid, only : region_mask
    use grid, only : TAREA

    use prognostic , only : curtime
    use prognostic , only : oldtime
    use prognostic , only : tracer_field_type => tracer_field

    use global_reductions, only : global_sum
    use constants, only : field_loc_center
    use broadcast, only : broadcast_scalar
    use time_management, only : init_time_flag
    use passive_tracer_tools, only : set_tracer_indices
    use mcog, only : mcog_nbins
    use registry, only : registry_match
    use named_field_mod, only : named_field_register
    use ecosys_forcing_mod, only : ecosys_forcing_init
    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_saved_state_get_var_vals

    integer (int_kind)       , intent(in)    :: ecosys_driver_ind_begin ! starting index of ecosys tracers in global tracer
                                                                        ! array, passed through to rest_read_tracer_block
    character (*)            , intent(in)    :: init_ts_file_fmt      ! format (bin or nc) for input file
    character (*)            , intent(in)    :: read_restart_filename ! file name for restart file
    type (tracer_field_type) , intent(inout) :: tracer_d_module(:)    ! descriptors for each tracer
    real (r8)                , intent(inout) :: tracer_module(:,:,:,:,:,:)
    character (char_len)     , intent(out)   :: tadvect_ctype(:)      ! advection method for ecosys tracers
    logical                  , intent(out)   :: lhas_riv_flux(:)      ! true if a tracer has a riverine flux
    integer (POP_i4)         , intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    ! FIXME (MNL): we should treat namelists in the driver in a POP-ish manner
    !              rather than a MARBL-ish manner, but the MARBL process is
    !              less error-prone (reading namelist on each task rather
    !              than relying on broadcasts to ensure every task has every
    !              namelist variable set correctly).
    integer(int_kind), parameter :: pop_in_tot_len    = 262144
    integer(int_kind), parameter :: pop_in_nl_max_len = 32768
    integer(int_kind), parameter :: pop_in_nl_cnt     = 256

    character(len=*), parameter      :: subname = 'ecosys_driver:ecosys_driver_init'
    character(char_len)              :: log_message
    integer (int_kind)               :: cumulative_nt, n, bid, k, i, j
    integer (int_kind)               :: num_fields
    integer (int_kind)               :: nml_error                          ! error flag for nml read
    integer (int_kind)               :: iostat                             ! io status flag
    character (char_len)             :: sname, lname, units, coordinates
    character (4)                    :: grid_loc
    integer (int_kind)               :: auto_ind                           ! autotroph functional group index
    integer (int_kind)               :: iblock                             ! index for looping over blocks
    character (char_len)             :: init_file_fmt                      ! file format for option 'file'
    logical   (log_kind)             :: lmarginal_seas                     ! is ecosystem active in marginal seas?
    integer(int_kind)                :: marbl_actual_tracer_cnt            ! # of tracers actually in MARBL
    integer (int_kind)               :: glo_avg_field_cnt
    real (r8), allocatable           :: rmean_vals(:)
    ! Variables for processing namelists
    character (char_len)             :: marbl_settings_file                  ! name of file containing MARBL settings
    character(len=pop_in_nl_max_len) :: nl_buffer(pop_in_nl_cnt)
    character(len=pop_in_nl_max_len) :: tmp_nl_buffer
    character(len=256)               :: marbl_nl_line
    character(len=pop_in_tot_len)    :: nl_str
    character(char_len_long)         :: ioerror_msg
    logical(log_kind)                :: lscalar ! true => glo_scalar_rmean_ind, false => glo_avg_rmean_ind

    !-----------------------------------------------------------------------
    !  read in ecosys_driver namelist, to set namelist parameters that
    !  should be the same for all ecosystem-related modules
    !-----------------------------------------------------------------------

    namelist /ecosys_driver_nml/ &
         lmarginal_seas, ecosys_tadvect_ctype, ecosys_qsw_distrb_const, &
         marbl_settings_file

    errorCode = POP_Success

    lmarginal_seas            = .true.
    ecosys_tadvect_ctype      = 'base_model'
    ecosys_qsw_distrb_const   = .true.
    marbl_settings_file       = 'marbl_in'

    ! -----------------------------
    ! Read pop namelist into string
    ! -----------------------------

    ! (a) initialize
    nl_buffer(:) = ''
    nl_str       = ''

    ! (b) Only master task reads
    if (my_task == master_task) then
       ! read the POP namelist file into a buffer.
       open(unit=nml_in, file=nml_filename, action='read', access='stream', &
            form='unformatted', iostat=nml_error)
       if (nml_error == 0) then
          read(unit=nml_in, iostat=nml_error, iomsg=ioerror_msg) nl_str

          ! we should always reach the EOF to capture the entire file...
          if (.not. is_iostat_end(nml_error)) then
             write(stdout, '(4a, i0)') subname, &
                  ": IO ERROR reading ", trim(nml_filename), " into buffer: ", nml_error
             write(stdout, '(a)') ioerror_msg
          else
             write(stdout, '(a, a, a)') "Read '", trim(nml_filename), "' until EOF."
          end if

          write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
               len_trim(nl_str), " characters."

          write(stdout, '(a)') "  If it looks like part of the namelist is missing, "
          write(stdout, '(a)') "  compare the number of characters read to the actual "
          write(stdout, '(a)') "  size of your file ($ wc -c pop_in) and increase "
          write(stdout, '(a)') "  the buffer size if necessary."
       else
          write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
               "opening namelist file : ", trim(nml_filename)
       end if
       close(nml_in)
    end if

    ! (c) if error reading, let all tasks know / abort
    call broadcast_scalar(nml_error, master_task)
    if (.not. is_iostat_end(nml_error)) then
       ! NOTE(bja, 2015-01) assuming that eof is the only proper exit
       ! code from the read.
       call document(subname, 'ERROR reading pop namelist file into buffer.')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    ! (d) Otherwise broadcast namelist contents to all tasks
    call broadcast_scalar(nl_str, master_task)

    ! (e) process namelist string to store in nl_buffer
    call namelist_split_by_nl(nl_str, nl_buffer)

    ! (f) read ecosys_driver_nml
    tmp_nl_buffer = namelist_find(nl_buffer, 'ecosys_driver_nml')
    read(tmp_nl_buffer, nml=ecosys_driver_nml, iostat=nml_error, iomsg=ioerror_msg)
    if (nml_error /= 0) then
       write(stdout, *) subname, ": process ", my_task, ": namelist read error: ", nml_error, " : ", ioerror_msg
       call exit_POP(sigAbort, 'ERROR reading ecosys_driver_nml from buffer.')
    end if

    ! (g) Output namelist to log
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

    call get_timer(ecosys_interior_forcing_timer   , 'ECOSYS_INTERIOR_FORCING'          , 1             , distrb_clinic%nprocs)
    call get_timer(ecosys_interior_pop_to_marbl    , 'ECOSYS_INTERIOR_POP_TO_MARBL_XFER', nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_interior_marbl_to_pop    , 'ECOSYS_INTERIOR_MARBL_TO_POP_XFER', nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_interior_marbl_tavg      , 'ECOSYS_INTERIOR_MARBL_TAVG'       , nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_interior_timer           , 'ECOSYS_INTERIOR'                  , nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_surface_global_sum_timer , 'ECOSYS_SURFACE_GLOBAL_SUM'        , 1             , distrb_clinic%nprocs)
    call get_timer(ecosys_interior_global_sum_timer, 'ECOSYS_INTERIOR_GLOBAL_SUM'       , 1             , distrb_clinic%nprocs)
    call get_timer(ecosys_set_sflux_timer          , 'ECOSYS_SET_SFLUX'                 , nblocks_clinic, distrb_clinic%nprocs)

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
    !  Initialize MARBL
    !--------------------------------------------------------------------

    ! (1) Read marbl input file and call put_setting()

    ! (1a) only master task opens file
    if (my_task == master_task) &
       ! read the marbl_in into buffer
       open(unit=nml_in, file=marbl_settings_file, iostat=nml_error)
    call broadcast_scalar(nml_error, master_task)
    if (nml_error .ne. 0) then
      write(stdout, '(a, a, i8, a, a)') subname, ": IO ERROR ", nml_error, &
           "opening namelist file : ", trim(marbl_settings_file)
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! (1b) master task reads file and broadcasts line-by-line
    marbl_nl_line = ''
    do
      ! i. Read next line on master, iostat value out
      !    (Exit loop if read is not successful; either read error or end of file)
      if (my_task .eq. master_task) read(nml_in, "(A)", iostat=nml_error) marbl_nl_line
      call broadcast_scalar(nml_error, master_task)
      if (nml_error .ne. 0) exit

      ! ii. Broadcast line just read in on master_task to all tasks
      call broadcast_scalar(marbl_nl_line, master_task)

      ! iii. Call put_setting from each block
      do iblock=1, max(1,nblocks_clinic)
         ! call marbl_instance%put line by line
         call marbl_instances(iblock)%put_setting(marbl_nl_line)
      end do
    end do

    ! (1c) we should always reach the EOF to capture the entire file...
    if (.not. is_iostat_end(nml_error)) then
       write(stdout, '(4a, i0)') subname, &
            ": IO ERROR reading ", trim(marbl_settings_file), ": ", nml_error
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    else
       if (my_task == master_task) &
         write(stdout, '(a, a, a)') "Read '", trim(marbl_settings_file), "' until EOF."
    end if
    if (my_task == master_task) close(nml_in)

    ! (2) Set up mapping between columns we pass to MARBL and the POP indices
    call gen_marbl_to_pop_index_mapping()

    ! (3) call marbl%init()
    do iblock=1, max(1,nblocks_clinic)

       call marbl_instances(iblock)%init(                                     &
            gcm_num_levels = km,                                              &
            gcm_num_PAR_subcols = mcog_nbins,                                 &
            gcm_num_elements_surface_flux = marbl_col_cnt(iblock),            &
            gcm_delta_z = dz,                                                 &
            gcm_zw = zw,                                                      &
            gcm_zt = zt,                                                      &
            lgcm_has_global_ops = .true.)

       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, ")%init()"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()

    end do

    ! Is ciso enabled in this run?
    call marbl_instances(1)%get_setting('ciso_on', ciso_on)
    marbl_actual_tracer_cnt = size(marbl_instances(1)%tracer_metadata)

    ! Make sure MARBL tracer count lines up with what POP expects
    if (marbl_actual_tracer_cnt.ne.ecosys_tracer_cnt) then
      write(log_message,"(A,I0,A,I0)") 'MARBL is computing tendencies for ', &
                 marbl_actual_tracer_cnt, ' tracers, but POP is expecting ', &
                 ecosys_tracer_cnt
      call document(subname, log_message)
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    !--------------------------------------------------------------------
    !  Initialize ecosys forcing fields
    !--------------------------------------------------------------------

    ! Determine MARBL tracer indices for tracers that use virtual fluxes
    ! (must be done prior to call to ecosys_forcing_init, when vflux is set)
    dic_ind   = marbl_instances(1)%get_tracer_index('DIC')
    alk_ind   = marbl_instances(1)%get_tracer_index('ALK')
    dic_alt_co2_ind = marbl_instances(1)%get_tracer_index('DIC_ALT_CO2')
    alk_alt_co2_ind = marbl_instances(1)%get_tracer_index('ALK_ALT_CO2')
    di13c_ind = marbl_instances(1)%get_tracer_index('DI13C')
    di14c_ind = marbl_instances(1)%get_tracer_index('DI14C')

    o2_ind       = marbl_instances(1)%get_tracer_index('O2')
    no3_ind      = marbl_instances(1)%get_tracer_index('NO3')
    po4_ind      = marbl_instances(1)%get_tracer_index('PO4')
    don_ind      = marbl_instances(1)%get_tracer_index('DON')
    donr_ind     = marbl_instances(1)%get_tracer_index('DONr')
    dop_ind      = marbl_instances(1)%get_tracer_index('DOP')
    dopr_ind     = marbl_instances(1)%get_tracer_index('DOPr')
    sio3_ind     = marbl_instances(1)%get_tracer_index('SiO3')
    fe_ind       = marbl_instances(1)%get_tracer_index('Fe')
    doc_ind      = marbl_instances(1)%get_tracer_index('DOC')
    docr_ind     = marbl_instances(1)%get_tracer_index('DOCr')
    do13ctot_ind = marbl_instances(1)%get_tracer_index('DO13Ctot')
    do14ctot_ind = marbl_instances(1)%get_tracer_index('DO14Ctot')

    ! pass ecosys_forcing_data_nml
    ! to ecosys_forcing_init()
    tmp_nl_buffer = namelist_find(nl_buffer, 'ecosys_forcing_data_nml')

    call ecosys_forcing_init(ciso_on,                                         &
                             land_mask,                                       &
                             marbl_instances(1)%surface_flux_forcings,        &
                             marbl_instances(1)%interior_tendency_forcings,   &
                             tmp_nl_buffer,                                   &
                             lhas_riv_flux)

    !--------------------------------------------------------------------
    !  Initialize ecosys tracers and saved state
    !--------------------------------------------------------------------

    call ecosys_saved_state_setup(surface_flux_saved_state, &
         marbl_instances(1)%surface_flux_saved_state)
    call ecosys_saved_state_setup(interior_tendency_saved_state, &
         marbl_instances(1)%interior_tendency_saved_state)

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
    tmp_nl_buffer = namelist_find(nl_buffer, 'ecosys_tracer_init_nml')

    call ecosys_tracers_and_saved_state_init(                    &
       marbl_instances(1),                                       &
       ecosys_driver_ind_begin,                                  &
       ciso_on,                                                  &
       init_ts_file_fmt,                                         &
       read_restart_filename,                                    &
       tracer_d_module(:),                                       &
       marbl_instances(1)%tracer_metadata(:)%tracer_module_name, &
       tmp_nl_buffer,                                            &
       land_mask,                                                &
       tracer_module(:,:,:,:,:,:),                               &
       errorCode)

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_driver_init')
       return
    endif

    ! copy values from POP's running mean to MARBL interface
    allocate(rmean_vals(size(marbl_instances(1)%glo_avg_rmean_interior_tendency)))
    lscalar = .false.
    call ecosys_running_mean_saved_state_get_var_vals('interior_tendency', lscalar, rmean_vals(:))
    do n = 1, size(rmean_vals)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_avg_rmean_interior_tendency(n)%rmean = rmean_vals(n)
       end do
    end do
    deallocate(rmean_vals)

    allocate(rmean_vals(size(marbl_instances(1)%glo_avg_rmean_surface_flux)))
    lscalar = .false.
    call ecosys_running_mean_saved_state_get_var_vals('surface_flux', lscalar, rmean_vals(:))
    do n = 1, size(rmean_vals)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_avg_rmean_surface_flux(n)%rmean = rmean_vals(n)
       end do
    end do
    deallocate(rmean_vals)

    allocate(rmean_vals(size(marbl_instances(1)%glo_scalar_rmean_interior_tendency)))
    lscalar = .true.
    call ecosys_running_mean_saved_state_get_var_vals('interior_tendency', lscalar, rmean_vals(:))
    do n = 1, size(rmean_vals)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_scalar_rmean_interior_tendency(n)%rmean = rmean_vals(n)
       end do
    end do
    deallocate(rmean_vals)

    allocate(rmean_vals(size(marbl_instances(1)%glo_scalar_rmean_surface_flux)))
    lscalar = .true.
    call ecosys_running_mean_saved_state_get_var_vals('surface_flux', lscalar, rmean_vals(:))
    do n = 1, size(rmean_vals)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_scalar_rmean_surface_flux(n)%rmean = rmean_vals(n)
       end do
    end do
    deallocate(rmean_vals)

    !--------------------------------------------------------------------
    !  Initialize ecosys_driver module variables
    !--------------------------------------------------------------------

    associate(diags => marbl_instances(1)%surface_flux_diags%diags)
      allocate(surface_flux_diags(nx_block, ny_block, size(diags), nblocks_clinic))
    end associate

    surface_flux_diags = c0

    tadvect_ctype(1:ecosys_tracer_cnt) = ecosys_tadvect_ctype

    !--------------------------------------------------------------------
    ! Initialize tavg ids (need only do this using first block)
    !--------------------------------------------------------------------

    call ecosys_tavg_init(marbl_instances(1))

    !--------------------------------------------------------------------
    ! Register and set Chl field for short-wave absorption
    ! Register and set fco2, the  air-sea co2 gas flux
    !--------------------------------------------------------------------

    sfo_cnt = 0

    ! Register flux_co2 with MARBL surface flux outputs
    sfo_cnt = sfo_cnt + 1
    do iblock=1, max(1,nblocks_clinic)
       call marbl_instances(iblock)%surface_flux_output%add_sfo(              &
              num_elements = marbl_col_cnt(iblock),                           &
              field_name   = "flux_co2",                                      &
              sfo_id       = flux_co2_id,                                     &
              marbl_status_log = marbl_instances(iblock)%StatusLog)
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, &
                                     ")%surface_flux_output%add_sfo(flux_co2)"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()
    end do

    ! Register totalChl with MARBL surface flux outputs
    sfo_cnt = sfo_cnt + 1
    do iblock=1, max(1,nblocks_clinic)
       call marbl_instances(iblock)%surface_flux_output%add_sfo(              &
              num_elements = marbl_col_cnt(iblock),                           &
              field_name   = "totalChl",                                      &
              sfo_id       = totalChl_id,                                     &
              marbl_status_log = marbl_instances(iblock)%StatusLog)
       if (marbl_instances(iblock)%StatusLog%labort_marbl) then
         write(log_message,"(A,I0,A)") "marbl(", iblock, &
                                     ")%surface_flux_output%add_sfo(totalChl)"
         call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
       end if
       call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
       call marbl_instances(iblock)%StatusLog%erase()
    end do

    allocate(surface_flux_outputs(nx_block,ny_block,nblocks_clinic,sfo_cnt))
    surface_flux_outputs = c0

    call named_field_register('SFLUX_CO2'        , sflux_co2_nf_ind)
    call named_field_register('model_chlorophyll', totChl_surf_nf_ind)

    !--------------------------------------------------------------------
    ! allocate space for fields for which global averages are to be computed
    !--------------------------------------------------------------------

    glo_avg_field_cnt = size(marbl_instances(1)%glo_avg_fields_interior_tendency, dim=1)
    allocate(glo_avg_fields_interior(nx_block, ny_block, nblocks_clinic, glo_avg_field_cnt))

    glo_avg_field_cnt = size(marbl_instances(1)%glo_avg_fields_surface_flux, dim=2)
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

  end subroutine ecosys_driver_init

  !***********************************************************************

  subroutine ecosys_driver_set_interior_forcing(FRACR_BIN, QSW_RAW_BIN, QSW_BIN)

    ! !DESCRIPTION:
    !  prepare forcing for computation of interior source-sink terms

    use mcog               , only : mcog_nbins
    use ecosys_forcing_mod , only : ecosys_forcing_set_interior_time_varying_forcing_data

    real (r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic) , intent(in)    :: FRACR_BIN         ! fraction of cell occupied by mcog bin
    real (r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic) , intent(in)    :: QSW_RAW_BIN       ! raw (directly from cpl) shortwave into each mcog column (W/m^2)
    real (r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic) , intent(in)    :: QSW_BIN           ! shortwave into each mcog bin, potentially modified by coszen factor (W/m^2)

    call timer_start(ecosys_interior_forcing_timer)

    call ecosys_forcing_set_interior_time_varying_forcing_data( &
         FRACR_BIN, QSW_RAW_BIN, QSW_BIN, ecosys_qsw_distrb_const, land_mask)

    call timer_stop(ecosys_interior_forcing_timer)

  end subroutine ecosys_driver_set_interior_forcing

  !***********************************************************************

  subroutine ecosys_driver_set_interior(     &
       TRACER_MODULE_OLD, TRACER_MODULE_CUR, &
       DTRACER_MODULE, this_block)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute source-sink terms
    !  accumulate commnon tavg fields related to source-sink terms

    use grid               , only : KMT
    use grid               , only : DZT
    use grid               , only : partial_bottom_cells
    use grid               , only : TLOND, TLATD
    use ecosys_forcing_mod , only : interior_tendency_forcings

    real (r8), dimension(:,:,:,:), intent(in)    :: TRACER_MODULE_OLD ! old tracer values
    real (r8), dimension(:,:,:,:), intent(in)    :: TRACER_MODULE_CUR ! current tracer values
    type (block)                 , intent(in)    :: this_block        ! block information for this block
    real (r8), dimension(:,:,:,:), intent(inout) :: DTRACER_MODULE    ! computed source/sink terms

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_driver:ecosys_driver_set_interior'
    character(len=char_len)     :: log_message

    integer (int_kind) :: i   ! nx_block loop index
    integer (int_kind) :: c   ! ny_block / column loop index
    integer (int_kind) :: bid ! local block address for this block
    integer (int_kind) :: n, d, ncols, k
    !-----------------------------------------------------------------------

    bid = this_block%local_id

    call timer_start(ecosys_interior_timer, block_id=bid)

    do c = this_block%jb,this_block%je
       do i = this_block%ib,this_block%ie

          if (land_mask(i,c,bid)) then

             !-----------------------------------------------------------
             ! Copy data from slab to column
             !-----------------------------------------------------------

             call timer_start(ecosys_interior_pop_to_marbl, block_id=bid)

             ! --- set marbl_domain kmt and if partial bottom cells then also delta_z ---

             marbl_instances(bid)%domain%kmt = KMT(i, c, bid)
             if (partial_bottom_cells) then
                marbl_instances(bid)%domain%delta_z(:) = DZT(i, c, :, bid)
             end if

             ! --- set forcing fields ---

             do n = 1, size(interior_tendency_forcings)
               if (interior_tendency_forcings(n)%rank == 2) then
                 marbl_instances(bid)%interior_tendency_forcings(n)%field_0d(1) = &
                      interior_tendency_forcings(n)%field_0d(i,c,bid)
               else
                 marbl_instances(bid)%interior_tendency_forcings(n)%field_1d(1,:) = &
                      interior_tendency_forcings(n)%field_1d(i,c,:,bid)
               end if
             end do

             ! --- set column tracers, averaging 2 time levels into 1 ---

             do n = 1, ecosys_tracer_cnt
                marbl_instances(bid)%tracers(n, :) = p5*(tracer_module_old(i, c, :, n) + tracer_module_cur(i, c, :, n))
             end do

             ! --- copy data from slab to column for marbl_saved_state ---
             do n=1,size(interior_tendency_saved_state)
               marbl_instances(bid)%interior_tendency_saved_state%state(n)%field_3d(:,1) = &
                 interior_tendency_saved_state(n)%field_3d(:,i,c,bid)
             end do
             call timer_stop(ecosys_interior_pop_to_marbl, block_id=bid)

             !-----------------------------------------------------------
             !  compute time derivatives for ecosystem state variables
             !-----------------------------------------------------------

             call marbl_instances(bid)%interior_tendency_compute()
             if (marbl_instances(bid)%StatusLog%labort_marbl) then
                write(log_message,"(A,I0,A)") "marbl_instances(", bid, &
                                              ")%set_interior_forcing()"
                call marbl_instances(bid)%StatusLog%log_error_trace(log_message, subname)
             end if
             call print_marbl_log(marbl_instances(bid)%StatusLog, bid, i, c)
             call marbl_instances(bid)%StatusLog%erase()

             !-----------------------------------------------------------
             ! copy marbl column data back to slab
             !-----------------------------------------------------------

             call timer_start(ecosys_interior_marbl_to_pop, block_id=bid)

             do n=1,size(interior_tendency_saved_state)
               interior_tendency_saved_state(n)%field_3d(:,i,c,bid) =               &
                 marbl_instances(bid)%interior_tendency_saved_state%state(n)%field_3d(:,1)
             end do

             !-----------------------------------------------------------
             ! before copying tendencies, check to see if any are NaNs
             !-----------------------------------------------------------

             do k = 1, KMT(i, c, bid)
                if (any(shr_infnan_isnan(marbl_instances(bid)%interior_tendencies(:, k)))) then
                   write(stdout, *) subname, ': NaN in dtracer_module, (i,j,k)=(', &
                      this_block%i_glob(i), ',', this_block%j_glob(c), ',', k, ')'
                   write(stdout, *) '(lon,lat)=(', TLOND(i,c,bid), ',', TLATD(i,c,bid), ')'
                   do n = 1, ecosys_tracer_cnt
                      write(stdout, *) trim(marbl_instances(1)%tracer_metadata(n)%short_name), ' ', &
                         marbl_instances(bid)%tracers(n, k), ' ', &
                         marbl_instances(bid)%interior_tendencies(n, k)
                   end do
                   do n = 1, size(interior_tendency_forcings)
                      associate (forcing_field => interior_tendency_forcings(n))
                         write(stdout, *) trim(forcing_field%metadata%marbl_varname)
                         if (forcing_field%rank == 2) then
                            write(stdout, *) forcing_field%field_0d(i,c,bid)
                         else
                            if (forcing_field%ldim3_is_depth) then
                               write(stdout, *) forcing_field%field_1d(i,c,k,bid)
                            else
                               write(stdout, *) forcing_field%field_1d(i,c,:,bid)
                            end if
                         end if
                      end associate
                   end do
                   call exit_POP(sigAbort, 'Stopping in ' // subname)
                end if
             end do

             do n = 1, ecosys_tracer_cnt
                dtracer_module(i, c, 1:KMT(i, c, bid), n) = marbl_instances(bid)%interior_tendencies(n, 1:KMT(i, c, bid))
             end do

             ! copy values to be used in computing requested global averages
             ! arrays have zero extent if none are requested
             glo_avg_fields_interior(i, c, bid, :) = marbl_instances(bid)%glo_avg_fields_interior_tendency(:)
             call timer_stop(ecosys_interior_marbl_to_pop, block_id=bid)

             !-----------------------------------------------------------
             ! Update pop tavg diags
             !-----------------------------------------------------------

             call timer_start(ecosys_interior_marbl_tavg, block_id=bid)

             call ecosys_tavg_accumulate_interior(i, c, marbl_instances(bid), bid)

             call timer_stop(ecosys_interior_marbl_tavg, block_id=bid)

          end if ! end if land_mask > 0

       end do ! do i
    end do ! do c

    call timer_stop(ecosys_interior_timer, block_id=bid)

  end subroutine ecosys_driver_set_interior

  !***********************************************************************

  subroutine ecosys_driver_set_sflux_forcing( &
       u10_sqr,                       &
       ifrac,                         &
       press,                         &
       atm_fine_dust_flux,            &
       atm_coarse_dust_flux,          &
       seaice_dust_flux,              &
       atm_black_carbon_flux,         &
       seaice_black_carbon_flux,      &
       sst,                           &
       sss)

    ! DESCRIPTION:
    ! Sets surface flux forcing data

    use ecosys_forcing_mod   , only : ecosys_forcing_set_surface_time_varying_forcing_data

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: u10_sqr                  ! 10m wind speed squared (cm/s)**2
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: ifrac                    ! sea ice fraction (non-dimensional)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: press                    ! sea level atmospheric pressure (dyne/cm**2)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: atm_fine_dust_flux       ! fine dust flux from atm (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: atm_coarse_dust_flux     ! coarse dust flux from atm (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: seaice_dust_flux         ! dust flux from seaice (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: atm_black_carbon_flux    ! black carbon flux from atm (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: seaice_black_carbon_flux ! black carbon flux from seaice (g/cm**2/s)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: sst                      ! sea surface temperature (c)
    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) , intent(in)    :: sss                      ! sea surface salinity (psu)

    !-----------------------------------------------------------------------
    ! Set surface flux forcing data
    !-----------------------------------------------------------------------

    call ecosys_forcing_set_surface_time_varying_forcing_data( &
         ciso_on,                         &
         land_mask,                       &
         u10_sqr,                         &
         ifrac,                           &
         press,                           &
         atm_fine_dust_flux,              &
         atm_coarse_dust_flux,            &
         seaice_dust_flux,                &
         atm_black_carbon_flux,           &
         seaice_black_carbon_flux,        &
         sst,                             &
         sss)

  end subroutine ecosys_driver_set_sflux_forcing


  !***********************************************************************

  subroutine ecosys_driver_set_sflux( &
       tracers_at_surface_old,        &
       tracers_at_surface_cur,        &
       stf_module,                    &
       stf_riv_module,                &
       iblock)

    ! DESCRIPTION:
    ! Calls marbl to compute surface tracer fluxes

    use ecosys_forcing_mod   , only : surface_flux_forcings
    use ecosys_forcing_mod   , only : ecosys_forcing_comp_stf_riv
    use blocks               , only : get_block
    use domain               , only : blocks_clinic
    use grid                 , only : TLOND, TLATD

    real (r8), dimension(:,:,:), intent(in)    :: tracers_at_surface_old
    real (r8), dimension(:,:,:), intent(in)    :: tracers_at_surface_cur  ! module tracers
    real (r8), dimension(:,:,:), intent(inout) :: stf_module
    real (r8), dimension(:,:,:), intent(inout) :: stf_riv_module
    integer (int_kind)         , intent(in)    :: iblock

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_driver:ecosys_driver_set_sflux'
    character(len=char_len)     :: log_message

    integer (int_kind) :: index_marbl  ! marbl index
    integer (int_kind) :: i, j, n      ! pop loop indices
    type(block) :: this_block
    !-----------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    if (marbl_col_cnt(iblock) .eq. 0) return

    call timer_start(ecosys_set_sflux_timer, iblock)

    ! Set up block id
    this_block = get_block(blocks_clinic(iblock), iblock)

    !-----------------------------------------------------------------------
    ! Copy data from slab data structure to column input for marbl
    !-----------------------------------------------------------------------

    do index_marbl = 1, marbl_col_cnt(iblock)
       i = marbl_col_to_pop_i(index_marbl,iblock)
       j = marbl_col_to_pop_j(index_marbl,iblock)

       do n = 1,size(surface_flux_forcings)
          marbl_instances(iblock)%surface_flux_forcings(n)%field_0d(index_marbl) = &
               surface_flux_forcings(n)%field_0d(i,j,iblock)
       end do

       do n = 1,ecosys_tracer_cnt
          marbl_instances(iblock)%tracers_at_surface(index_marbl,n) = &
               p5*(tracers_at_surface_old(i,j,n) + tracers_at_surface_cur(i,j,n))
       end do

       do n=1,size(surface_flux_saved_state)
         marbl_instances(iblock)%surface_flux_saved_state%state(n)%field_2d(index_marbl) = &
           surface_flux_saved_state(n)%field_2d(i,j,iblock)
       end do

    end do

    !-----------------------------------------------------------------------
    ! Compute surface fluxes in MARBL
    !-----------------------------------------------------------------------

    call marbl_instances(iblock)%surface_flux_compute()

    if (marbl_instances(iblock)%StatusLog%labort_marbl) then
       write(log_message,"(A,I0,A)") "marbl_instances(", iblock, &
                                     ")%surface_flux_compute()"
       call marbl_instances(iblock)%StatusLog%log_error_trace(log_message, subname)
    end if
    call print_marbl_log(marbl_instances(iblock)%StatusLog, iblock)
    call marbl_instances(iblock)%StatusLog%erase()

    !-----------------------------------------------------------------------
    ! Copy data from marbl output column to pop slab data structure
    !-----------------------------------------------------------------------

    do index_marbl = 1, marbl_col_cnt(iblock)
       i = marbl_col_to_pop_i(index_marbl,iblock)
       j = marbl_col_to_pop_j(index_marbl,iblock)

       do n=1,size(surface_flux_saved_state)
         surface_flux_saved_state(n)%field_2d(i,j,iblock) = &
           marbl_instances(iblock)%surface_flux_saved_state%state(n)%field_2d(index_marbl)
       end do

       do n=1,sfo_cnt
         surface_flux_outputs(i,j,iblock,n) = &
            marbl_instances(iblock)%surface_flux_output%sfo(n)%forcing_field(index_marbl)
       end do

       !-----------------------------------------------------------
       ! before copying surface fluxes, check to see if any are NaNs
       !-----------------------------------------------------------

       if (any(shr_infnan_isnan(marbl_instances(iblock)%surface_fluxes(index_marbl,:)))) then
          write(stdout, *) subname, ': NaN in stf_module, (i,j)=(', &
             this_block%i_glob(i), ',', this_block%j_glob(j), ')'
          write(stdout, *) '(lon,lat)=(', TLOND(i,j,iblock), ',', TLATD(i,j,iblock), ')'
          do n = 1, ecosys_tracer_cnt
             write(stdout, *) trim(marbl_instances(1)%tracer_metadata(n)%short_name), ' ', &
                marbl_instances(iblock)%tracers_at_surface(index_marbl,n), ' ', &
                marbl_instances(iblock)%surface_fluxes(index_marbl,n)
          end do
          do n = 1, size(surface_flux_forcings)
             associate (forcing_field => surface_flux_forcings(n))
                write(stdout, *) trim(forcing_field%metadata%marbl_varname)
                if (forcing_field%rank == 2) then
                   write(stdout, *) forcing_field%field_0d(i,j,iblock)
                else
                   write(stdout, *) forcing_field%field_1d(i,j,:,iblock)
                end if
             end associate
          end do
          call exit_POP(sigAbort, 'Stopping in ' // subname)
       end if

       do n = 1,ecosys_tracer_cnt
          stf_module(i,j,n) = &
               marbl_instances(iblock)%surface_fluxes(index_marbl,n)
       end do

       do n=1,size(marbl_instances(1)%surface_flux_diags%diags)
          surface_flux_diags(i,j,n,iblock) = &
               marbl_instances(iblock)%surface_flux_diags%diags(n)%field_2d(index_marbl)
       end do

       ! copy values to be used in computing requested global averages
       ! arrays have zero extent if none are requested
       glo_avg_fields_surface(i,j,iblock,:) = marbl_instances(iblock)%glo_avg_fields_surface_flux(index_marbl,:)
    end do

    !-----------------------------------------------------------------------
    ! compute river fluxes
    !-----------------------------------------------------------------------

    call ecosys_forcing_comp_stf_riv(ciso_on, stf_riv_module(:,:,:), iblock)

    call timer_stop(ecosys_set_sflux_timer, iblock)

  end subroutine ecosys_driver_set_sflux

  !***********************************************************************

  subroutine ecosys_driver_post_set_sflux()

    ! DESCRIPTION:
    ! sflux related code that occurs after the threaded region needs to be exited

    use POP_HaloMod          , only : POP_HaloUpdate
    use POP_GridHorzMod      , only : POP_gridHorzLocCenter
    use POP_FieldMod         , only : POP_fieldKindScalar
    use POP_ErrorMod         , only : POP_Success
    use domain               , only : POP_haloClinic
    use named_field_mod      , only : named_field_set
    use ecosys_forcing_saved_state_mod , only : ecosys_forcing_saved_state_update

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_driver:ecosys_driver_post_set_sflux'

    integer (int_kind) :: iblock    ! block loop index
    integer (int_kind) :: errorCode ! halo update error code
    !-----------------------------------------------------------------------

    ! Added halo update for totalChl surface flux output because sw_absorption
    ! had issue with prognostic Chl and we didn't track down why the halo cells
    ! mattered in that case
    ! klindsay thought -- KPP smooths boundary layer depth; CVMix computes KPP
    ! at all points (not just phys points) and then smooths boundary layer depth
    ! so computing vertical mixing coefficients only on phys points will require
    ! halo update for HBLT (and we may be able to remove this halo update as a
    ! result)
    call POP_HaloUpdate(surface_flux_outputs(:,:,:,totalChl_id), POP_haloclinic, &
         POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
    if (errorCode /= POP_Success) then
       call document(subname, 'error updating SFO halo for totalChl from MARBL')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    !-----------------------------------------------------------------------
    ! Update ecosys_forcing_saved_state fields
    !-----------------------------------------------------------------------

    call ecosys_forcing_saved_state_update(surface_flux_outputs(:,:,:,flux_co2_id))

    !-----------------------------------------------------------------------
    ! Update named fields
    !-----------------------------------------------------------------------

    call named_field_set(totChl_surf_nf_ind, surface_flux_outputs(:,:,:,totalChl_id))

    !  set air-sea co2 gas flux named field, converting units from
    !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
    do iblock = 1, nblocks_clinic
       call named_field_set(sflux_co2_nf_ind, iblock, 44.0e-8_r8 * surface_flux_outputs(:,:,iblock,flux_co2_id))
    end do

  end subroutine ecosys_driver_post_set_sflux

  !***********************************************************************

  subroutine ecosys_driver_comp_global_averages(field_source)

    ! DESCRIPTION:
    ! perform global operations

    use global_reductions, only : global_sum_prod
    use constants        , only : field_loc_center
    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_saved_state_update
    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_saved_state_get_var_vals
    use POP_CommMod      , only : POP_Barrier

    character (*), intent(in) :: field_source ! 'interior_tendency' or 'surface_flux'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real(r8), pointer          :: glo_avg_fields(:,:,:,:)
    real(r8), allocatable      :: glo_avg(:)
    real(r8), allocatable      :: glo_avg_rmean(:)
    integer(int_kind)          :: glo_avg_field_cnt
    integer(int_kind)          :: n, iblock
    integer(int_kind), pointer :: timer
    logical(log_kind)          :: lscalar ! true => glo_scalar_rmean_ind, false => glo_avg_rmean_ind
    !-----------------------------------------------------------------------

    lscalar = .false. ! all global means in this routine are glo_avg_rmean
    if (trim(field_source) == 'interior_tendency') then
       timer             => ecosys_interior_global_sum_timer
       glo_avg_fields    => glo_avg_fields_interior(:,:,:,:)
    else
       timer             => ecosys_surface_global_sum_timer
       glo_avg_fields    => glo_avg_fields_surface(:,:,:,:)
    end if

    glo_avg_field_cnt = size(glo_avg_fields, dim=4)

    if (glo_avg_field_cnt /= 0) then
       allocate(glo_avg(glo_avg_field_cnt))
       allocate(glo_avg_rmean(glo_avg_field_cnt))

       ! call barrier to get consistent timing info
       ! overhead is low because we're already calling a global_sum right after this
       call POP_Barrier()

       call timer_start(timer)

       ! compute global means, and update their running means
       do n = 1, glo_avg_field_cnt
          glo_avg(n) = glo_avg_norm_fact * global_sum_prod(glo_avg_fields(:,:,:,n), &
             glo_avg_area_masked(:,:,:), distrb_clinic, field_loc_center)
       end do
       call ecosys_running_mean_saved_state_update(field_source, lscalar, glo_avg)
       call ecosys_running_mean_saved_state_get_var_vals(field_source, lscalar, glo_avg_rmean)

       call timer_stop(timer)

       ! store global means, and their running means, into appropriate component of marbl_instances
       if (trim(field_source) == 'interior_tendency') then
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_avg_averages_interior_tendency(:)    = glo_avg(:)
             marbl_instances(iblock)%glo_avg_rmean_interior_tendency(:)%rmean = glo_avg_rmean(:)
          end do
       else
          do iblock = 1, nblocks_clinic
             marbl_instances(iblock)%glo_avg_averages_surface_flux(:)    = glo_avg(:)
             marbl_instances(iblock)%glo_avg_rmean_surface_flux(:)%rmean = glo_avg_rmean(:)
          end do
       end if

       deallocate(glo_avg_rmean)
       deallocate(glo_avg)
    end if

  end subroutine ecosys_driver_comp_global_averages

  !***********************************************************************

  subroutine ecosys_driver_set_global_scalars(field_source)

    character (*), intent(in) :: field_source ! 'interior_tendency' or 'surface_flux'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind)          :: iblock

    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       call marbl_instances(iblock)%set_global_scalars(field_source)
    end do

    call ecosys_driver_update_scalar_rmeans(field_source)

    call ecosys_tavg_accumulate_scalar_rmeans(marbl_instances(1), field_source)

  end subroutine ecosys_driver_set_global_scalars

  !***********************************************************************

  subroutine ecosys_driver_update_scalar_rmeans(field_source)

    ! DESCRIPTION:
    ! update running means of scalar variables

    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_saved_state_update
    use ecosys_running_mean_saved_state_mod, only : ecosys_running_mean_saved_state_get_var_vals
    use io_types         , only : stdout

    character (*), intent(in) :: field_source ! 'interior_tendency' or 'surface_flux'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname   = 'ecosys_driver:ecosys_driver_update_scalar_rmeans'
    character(len=*), parameter :: fmt_str   = '(A,1X,A)'
    character(len=*), parameter :: fmt_str_i = '(A,1X,A,1X,I0)'
    character(len=*), parameter :: fmt_str_e = '(A,1X,A,1X,E23.16)'

    real(r8), allocatable :: rmean_vals(:)
    integer(int_kind) :: n, iblock
    logical(log_kind) :: lscalar ! true => glo_scalar_rmean_ind, false => glo_avg_rmean_ind
    !-----------------------------------------------------------------------

    lscalar = .true. ! all global means in this subroutine are glo_scalar_rmean
    if (trim(field_source) == 'interior_tendency') then
      allocate(rmean_vals(size(marbl_instances(1)%glo_scalar_interior_tendency)))
       do n = 1, size(rmean_vals)
          ! verify that all instances have same value of glo_scalar_interior
          do iblock = 2, nblocks_clinic
             if (marbl_instances(iblock)%glo_scalar_interior_tendency(n) /= &
                 marbl_instances(1)%glo_scalar_interior_tendency(n)) then
                write(stdout, fmt_str)   subname, 'mismatch in glo_scalar_interior values across MARBL instances'
                write(stdout, fmt_str_i) subname, 'rmean index', n
                write(stdout, fmt_str_e) subname, 'iblock 1 value', marbl_instances(1)%glo_scalar_interior_tendency(n)
                write(stdout, fmt_str_i) subname, 'iblock', iblock
                write(stdout, fmt_str_e) subname, 'mismatched value', marbl_instances(iblock)%glo_scalar_interior_tendency(n)
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             end if
          end do
       end do

       call ecosys_running_mean_saved_state_update(field_source, lscalar, marbl_instances(1)%glo_scalar_interior_tendency)
       call ecosys_running_mean_saved_state_get_var_vals(field_source, lscalar, rmean_vals)
       do iblock = 1, nblocks_clinic
          do n=1, size(rmean_vals)
             marbl_instances(iblock)%glo_scalar_rmean_interior_tendency(n)%rmean = rmean_vals(n)
          end do
      end do
    else
       allocate(rmean_vals(size(marbl_instances(1)%glo_scalar_surface_flux)))
       do n = 1, size(rmean_vals)
          ! verify that all instances have same value of glo_scalar_surface
          do iblock = 2, nblocks_clinic
             if (marbl_instances(iblock)%glo_scalar_surface_flux(n) /= marbl_instances(1)%glo_scalar_surface_flux(n)) then
                write(stdout, fmt_str)   subname, 'mismatch in glo_scalar_surface values across MARBL instances'
                write(stdout, fmt_str_i) subname, 'rmean index', n
                write(stdout, fmt_str_e) subname, 'iblock 1 value', marbl_instances(1)%glo_scalar_surface_flux(n)
                write(stdout, fmt_str_i) subname, 'iblock', iblock
                write(stdout, fmt_str_e) subname, 'mismatched value', marbl_instances(iblock)%glo_scalar_surface_flux(n)
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             end if
          end do
       end do

       call ecosys_running_mean_saved_state_update(field_source, lscalar, marbl_instances(1)%glo_scalar_surface_flux)
       call ecosys_running_mean_saved_state_get_var_vals(field_source, lscalar, rmean_vals)
       do iblock = 1, nblocks_clinic
          do n=1, size(rmean_vals)
             marbl_instances(iblock)%glo_scalar_rmean_surface_flux(n)%rmean = rmean_vals(n)
          end do
       end do
    end if

  end subroutine ecosys_driver_update_scalar_rmeans

  !***********************************************************************

  subroutine ecosys_driver_tavg_forcing(stf_module, bid)

    ! !DESCRIPTION:
    !  accumulate common tavg fields for tracer surface fluxes
    real (r8), dimension(:,:,:), intent(in) :: stf_module
    integer,                     intent(in) :: bid

    call ecosys_tavg_accumulate_surface(marbl_col_to_pop_i(1:marbl_col_cnt(bid),bid), &
                                        marbl_col_to_pop_j(1:marbl_col_cnt(bid),bid), &
                                        stf_module(:,:,:), marbl_instances(bid), bid)

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

    use ecosys_tracers_and_saved_state_mod, only : ecosys_saved_state_write_restart
    use ecosys_tracers_and_saved_state_mod, only : saved_state_field_3d
    use io_types, only : io_field_desc
    use io      , only : datafile

    character(len=*), intent(in)    :: action
    type (datafile),  intent(inout) :: restart_file

    type (io_field_desc), dimension(:), allocatable, save :: surface_iodesc
    type (io_field_desc), dimension(:), allocatable, save :: interior_iodesc
    integer :: n
    !-----------------------------------------------------------------------

    if (.not. allocated(surface_iodesc)) then
       allocate(surface_iodesc(size(surface_flux_saved_state)))
       allocate(interior_iodesc(size(interior_tendency_saved_state)))
       allocate(saved_state_field_3d(nx_block, ny_block, km, max_blocks_clinic))
    end if

    call ecosys_saved_state_write_restart(restart_file, action,               &
                                          surface_flux_saved_state, surface_iodesc)
    call ecosys_saved_state_write_restart(restart_file, action,               &
                                          interior_tendency_saved_state, interior_iodesc)

    if (trim(action) == 'write') then
       deallocate(saved_state_field_3d)
       deallocate(surface_iodesc)
       deallocate(interior_iodesc)
    endif

  end subroutine ecosys_driver_write_restart

  !*****************************************************************************

  subroutine ecosys_driver_print_marbl_timers(stats)

    use timers, only : timer_format, stats_fmt1, stats_fmt2, stats_fmt3, stats_fmt4
    use constants, only : bignum, blank_fmt, delim_fmt
    use global_reductions, only : global_maxval, global_minval, global_sum
    use domain, only : distrb_clinic

    ! Accumulate timing data from MARBL structure
    ! * marbl_instances(:)%timer_summary%cumulative_runtimes(:)
    ! Process:
    ! (1) On each "node" (MPI task), compute total runtime across blocks
    !     -- max over blocks in threaded regions, otherwise sum over blocks
    ! (2) Master task will write the max across all nodes to ocn.log
    ! (3) If additional stats are requested:
    !     -- Pass node and block times to master task
    !     -- Master task writes min, max, avg over both nodes and blocks to log

    logical, optional, intent(in) :: stats

    real(r8), allocatable, dimension(:,:) :: run_time
    real(r8), allocatable, dimension(:) :: block_max, block_min, block_tot
    real(r8) :: node_time, node_min, node_max, node_avg
    integer :: iblock, nblocks, num_timers, timer_id
    logical :: stats_local

    if (present(stats)) then
      stats_local = stats
    else
      stats_local = .false.
    end if

    do iblock=1, max(1, nblocks_clinic)
      call marbl_instances(iblock)%extract_timing()
    end do

    num_timers = marbl_instances(1)%timer_summary%num_timers
    allocate(run_time(num_timers, nblocks_clinic))
    allocate(block_min(num_timers))
    allocate(block_max(num_timers))
    allocate(block_tot(num_timers))
    do iblock=1,nblocks_clinic
      associate(timers => marbl_instances(iblock)%timer_summary)
        run_time(:,iblock) =timers%cumulative_runtimes(:)
      end associate
    end do
    if (nblocks_clinic > 0) then
      block_min(:) = minval(run_time, dim=2)
      block_max(:) = maxval(run_time, dim=2)
      block_tot(:) = sum(run_time, dim=2)
    else
      block_min(:) = c0
      block_max(:) = c0
      block_tot(:) = c0
    end if

    nblocks   = global_sum(nblocks_clinic, distrb_clinic)

    do timer_id=1,num_timers

      if (marbl_instances(1)%timer_summary%is_threaded(timer_id)) then
        ! FIXME: to account for situation where nblocks > num_threads > 1
        !        we need to compute block stats better
        ! call threaded_stats(run_time, block_min, block_max, block_tot)
        node_time = block_max(timer_id)
      else
        node_time = block_tot(timer_id)
      end if

      node_max = global_maxval(node_time)

      if (my_task.eq.master_task) then
        write(stdout, timer_format) timer_id, node_max,                       &
                  trim(marbl_instances(1)%timer_summary%names(timer_id))
      end if

      if (stats_local) then
        node_min = global_minval(node_time)
        node_avg = global_sum(node_time, distrb_clinic)/real(distrb_clinic%nprocs, r8)

        block_min(timer_id) = global_minval(block_min(timer_id))
        block_max(timer_id) = global_maxval(block_max(timer_id))
        block_tot(timer_id) = global_sum(block_tot(timer_id), distrb_clinic)

        if (my_task.eq.master_task) then
          if (marbl_instances(1)%timer_summary%is_threaded(timer_id)) then
            write(stdout, "(A)") "Per node stats come from max per block"
          else
            write(stdout, "(A)") "Per node stats come from sum per block"
          end if ! is_threaded

          write(stdout, stats_fmt1) node_min
          write(stdout, stats_fmt2) node_max
          write(stdout, stats_fmt3) node_avg
          write(stdout, stats_fmt4) block_min(timer_id)
          write(stdout, stats_fmt2) block_max(timer_id)
          write(stdout, stats_fmt3) block_tot(timer_id)/real(nblocks, r8)

        end if ! master task
      end if ! stats_local
    end do ! timer loop

    if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
    endif
    deallocate(run_time, block_min, block_max, block_tot)

  end subroutine ecosys_driver_print_marbl_timers

  !*****************************************************************************

  subroutine gen_marbl_to_pop_index_mapping

    use blocks,        only : get_block
    use domain,        only : blocks_clinic

    character(len=*), parameter :: subname = 'ecosys_driver:gen_marbl_to_pop_index_mapping'

    type(block) :: this_block
    integer     :: index_marbl, i, j, iblock

    allocate(marbl_col_cnt(max(1,nblocks_clinic)))
    if (nblocks_clinic .eq. 0) marbl_col_cnt(1) = 0
    do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock), iblock)
      ! marbl_col_cnt is the number of ocean grid cells in the phys grid
      ! (no halo region)
      marbl_col_cnt(iblock) = count(land_mask(this_block%ib:this_block%ie, &
                                              this_block%jb:this_block%je, &
                                              iblock))
    end do

    allocate(marbl_col_to_pop_i(maxval(marbl_col_cnt), nblocks_clinic))
    allocate(marbl_col_to_pop_j(maxval(marbl_col_cnt), nblocks_clinic))

    do iblock = 1, nblocks_clinic
      index_marbl = 0
      this_block = get_block(blocks_clinic(iblock), iblock)
      do j = this_block%jb, this_block%je
        do i = this_block%ib, this_block%ie
          if (land_mask(i,j,iblock)) then
            index_marbl = index_marbl + 1
            if (index_marbl .gt. marbl_col_cnt(iblock)) then
              ! Exceeding expected number of columns
              write(stdout,"(3A)") 'INTERNAL ERROR (', subname, &
                                   '): index_marbl exceeds marbl_col_cnt'
              write(stdout,"(A,I0)") 'marbl_col_cnt = ', marbl_col_cnt(iblock)
              write(stdout,"(A,I0)") 'index_marbl = ', index_marbl
              call exit_POP(sigAbort, 'Stopping in ' // subname)
            end if
            marbl_col_to_pop_i(index_marbl, iblock) = i
            marbl_col_to_pop_j(index_marbl, iblock) = j
          end if
        end do
      end do
      if (index_marbl .lt. marbl_col_cnt(iblock)) then
        ! Exceeding expected number of columns
        write(stdout,"(3A)") 'INTERNAL ERROR (', subname, &
                             '): index_marbl is less than marbl_col_cnt'
        write(stdout,"(A,I0)") 'marbl_col_cnt = ', marbl_col_cnt(iblock)
        write(stdout,"(A,I0)") 'index_marbl = ', index_marbl
        call exit_POP(sigAbort, 'Stopping in ' // subname)
      end if
    end do

  end subroutine gen_marbl_to_pop_index_mapping

  !*****************************************************************************

  subroutine print_marbl_log(log_to_print, iblock, i, j)

    use marbl_logging, only : marbl_status_log_entry_type
    use marbl_logging, only : marbl_log_type
    use grid,          only : TLATD, TLOND
    use blocks,        only : get_block
    use domain,        only : blocks_clinic

    class(marbl_log_type), intent(in) :: log_to_print
    integer,               intent(in) :: iblock
    integer,     optional, intent(in) :: i, j

    character(len=*), parameter :: subname = 'ecosys_driver:print_marbl_log'
    character(len=char_len)     :: message_prefix, message_location
    type(marbl_status_log_entry_type), pointer :: tmp
    integer :: i_loc, j_loc, elem_old
    logical :: iam_master
    type(block) :: this_block

    ! Set up block id
    if (nblocks_clinic .eq. 0) return
    this_block = get_block(blocks_clinic(iblock), iblock)
    elem_old = -1
    write(message_prefix, "(A,I0,A,I0,A)") '(Task ', my_task, ', block ', iblock, ')'

    ! FIXME (02-2016,mnl): need better logic on which items to print
    iam_master = (my_task.eq.master_task).and.(iblock.eq.1)
    tmp => log_to_print%FullLog
    do while (associated(tmp))
      ! 1) Do I need to write this message? Yes, if all tasks should write this
      !    or if I am master_task
      if ((.not. tmp%lonly_master_writes) .or. iam_master) then
        ! 2) Print message location? (only if ElementInd changed and is positive)
        if (tmp%ElementInd .ne. elem_old) then
          if (tmp%ElementInd .gt. 0) then
            if (present(i) .and. present(j)) then
              write(message_location, "(A,F8.3,A,F7.3,A,I0,A,I0,A,I0)") &
                   'Message from (lon, lat) (', TLOND(i, j, iblock), ', ', &
                   TLATD(i, j, iblock), '), which is global (i,j) (', &
                   this_block%i_glob(i), ', ', this_block%j_glob(j), &
                   '). Level: ', tmp%ElementInd
            else
              i_loc = marbl_col_to_pop_i(tmp%ElementInd, iblock)
              j_loc = marbl_col_to_pop_j(tmp%ElementInd, iblock)
              write(message_location, "(A,F8.3,A,F7.3,A,I0,A,I0,A)") &
                   'Message from (lon, lat) (', TLOND(i_loc, j_loc, iblock), &
                   ", ", TLATD(i_loc, j_loc, iblock), '), which is global (i,j) (', &
                   this_block%i_glob(i_loc), ', ', this_block%j_glob(j_loc), ')'
            end if ! i,j present

            ! master task does not need prefix
            if (.not. iam_master) then
              write(stdout, "(A,1X,A)") trim(message_prefix), trim(message_location)
            else
              write(stdout, "(A)") trim(message_location)
            end if ! print message prefix?

          end if   ! ElementInd > 0
          elem_old = tmp%ElementInd
        end if     ! ElementInd /= elem_old

        ! master task does not need prefix
        if (.not. iam_master) then
          write(stdout, "(A,1X,A)") trim(message_prefix), trim(tmp%LogMessage)
        else
          write(stdout, "(A)") trim(tmp%LogMessage)
        end if     ! print message prefix?
      end if       ! write the message?
      tmp => tmp%next
    end do

    if (log_to_print%labort_marbl) then
      write(stdout, "(A)") 'ERROR reported from MARBL library'
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

  end subroutine print_marbl_log

  !*****************************************************************************

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
