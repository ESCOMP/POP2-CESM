! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_driver

  ! !DESCRIPTION:
  !  This module provides support for the ecosystem module and dependend tracer modules
  !  The base model calls subroutines in passive_tracers, which then call
  !  this module if ecosys_on is true. Ecosys_driver then calls subroutines
  !  in individual ecosystem modules (so far ecosys_mod and marbl_ciso_mod)
  !
  !  Written by: Alexandra Jahn, NCAR, Nov/Dec 2012

  ! !USES:

  use POP_KindsMod
  use POP_ErrorMod
  use POP_IOUnitsMod

  use kinds_mod                 , only : r8, int_kind, log_kind, char_len, char_len_long
  use blocks                    , only : block, nx_block, ny_block
  use domain                    , only : nblocks_clinic
  use domain                    , only : distrb_clinic
  use domain_size               , only : max_blocks_clinic, km, nt
  use io_types                  , only : stdout, nml_in, nml_filename
  use prognostic                , only : TRACER
  use exit_mod                  , only : sigAbort, exit_pop
  use constants                 , only : c0, c1, p5, delim_fmt, char_blank, ndelim_fmt

  use marbl_share_mod           , only : autotroph_cnt, zooplankton_cnt
  use marbl_share_mod           , only : ecosys_tracer_cnt
  use marbl_share_mod           , only : ecosys_ciso_tracer_cnt
  use marbl_share_mod           , only : ecosys_ind_begin, ecosys_ind_end
  use marbl_share_mod           , only : ecosys_ciso_ind_begin, ecosys_ciso_ind_end
  use marbl_share_mod           , only : ecosys_surface_share_type 

  use marbl_interface           , only : marbl_interface_class
  use marbl_interface           , only : marbl_sizes_type
  use marbl_interface_types     , only : marbl_diagnostics_type
  use marbl_interface_types     , only : photosynthetically_available_radiation_type
  use marbl_share_mod           , only : marbl_forcing_share_type
  use marbl_interface_types     , only : marbl_forcing_input_type
  use marbl_interface_types     , only : marbl_forcing_output_type
  use ecosys_mod                , only : marbl_ecosys_set_interior
  use ecosys_mod                , only : marbl_ecosys_compute_totalChl
  use ecosys_diagnostics_mod    , only : max_forcing_diags

  use ecosys_tavg               , only : ecosys_tavg_init
  use ecosys_tavg               , only : ecosys_tavg_accumulate
  use ecosys_tavg               , only : ecosys_tavg_accumulate_flux

  use ecosys_restore_mod        , only : ecosys_restore_type
  use passive_tracer_tools      , only : ind_name_pair

  use timers                    , only : timer_start
  use timers                    , only : timer_stop
  use timers                    , only : get_timer
  use mcog                      , only : mcog_nbins

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_driver_tracer_cnt_init
  public :: ecosys_driver_init
  public :: ecosys_driver_set_interior
  public :: ecosys_driver_set_sflux
  public :: ecosys_driver_tracer_ref_val
  public :: ecosys_driver_tavg_forcing
  public :: ecosys_driver_write_restart
  public :: ecosys_driver_unpack_source_sink_terms

  private :: ecosys_driver_init_tracers_and_saved_state
  private :: ecosys_driver_init_interior_restore 
  private :: ecosys_driver_init_sflux
  private :: ecosys_driver_read_sflux
  private :: ecosys_driver_write_ecosys_restart

  type :: ecosys_saved_state_type
     ! this struct is necessary because there is some global state
     ! that needs to be preserved for from one time step to the next
     real (r8) , dimension(:, :, :)         , allocatable :: dust_FLUX_IN       ! dust flux not stored in STF since dust is not prognostic
     real (r8) , dimension(:, :, :, :)      , allocatable :: PH_PREV_3D         ! computed pH_3D from previous time step
     real (r8) , dimension(:, :, :, :)      , allocatable :: PH_PREV_ALT_CO2_3D ! computed pH_3D from previous time step, alternative CO2
     real (r8) , dimension(:, :, :)         , allocatable :: ph_surf            ! computed ph from previous time step
     real (r8) , dimension(:, :, :)         , allocatable :: ph_surf_alt_co2    ! computed ph from previous time step, alternative CO2
  contains
     procedure :: construct => ecosys_saved_state_constructor
  end type ecosys_saved_state_type

  !-----------------------------------------------------------------------
  !  module variables required by forcing_passive_tracer
  !-----------------------------------------------------------------------

  integer (int_kind), public :: ecosys_driver_tracer_cnt

  !-----------------------------------------------------------------------
  ! timers
  !-----------------------------------------------------------------------

  integer (int_kind) :: ecosys_interior_timer
  integer (int_kind) :: ecosys_shr_strdata_advance_timer
  integer (int_kind) :: ecosys_pre_sflux_timer
  integer (int_kind) :: ecosys_set_sflux_timer
  integer (int_kind) :: ecosys_comp_CO3terms_timer
  integer (int_kind) :: ecosys_ciso_interior_timer
  integer (int_kind) :: ecosys_ciso_set_sflux_timer
  
  !--------------------------------------------------------------------
  ! removed from marbl_share because they are read from pop's
  ! namelist and passed into marbl (lmarginal_seas) or not used in
  ! marbl at all (tadvect, ecosys_qsw_distrb_const).
  ! --------------------------------------------------------------------

  logical(log_kind)           :: lmarginal_seas       ! Is ecosystem active in marginal seas ?
  character(char_len)         :: ecosys_tadvect_ctype ! advection method for ecosys tracers
  logical (log_kind) , public :: ecosys_qsw_distrb_const

  !-----------------------------------------------------------------------
  !  data storage for interaction with the marbl bgc library
  !-----------------------------------------------------------------------

  type(marbl_interface_class)     , dimension(max_blocks_clinic) :: marbl
  type(marbl_diagnostics_type)    , dimension(max_blocks_clinic) :: marbl_interior_diags
  type(marbl_diagnostics_type)    , dimension(max_blocks_clinic) :: marbl_restore_diags
  type(marbl_diagnostics_type)    , dimension(max_blocks_clinic) :: marbl_forcing_diags
  type(marbl_forcing_input_type)  , dimension(max_blocks_clinic) :: marbl_forcing_input
  type(marbl_forcing_output_type) , dimension(max_blocks_clinic) :: marbl_forcing_output
  type(marbl_forcing_share_type)  , dimension(max_blocks_clinic) :: marbl_forcing_share

  ! Computed diagnostics for surface fluxes
  real (r8) :: flux_diags(nx_block, ny_block, max_forcing_diags, max_blocks_clinic)

  !-----------------------------------------------------------------------
  ! pop data storage for interaction with marbl
  !-----------------------------------------------------------------------

  type(ecosys_saved_state_type) :: ecosys_saved_state 
  type(ecosys_restore_type)     :: ecosys_restore
  logical(log_kind), dimension(:, :, :), allocatable :: land_mask

  !-----------------------------------------------------------------------
  !  average surface tracer value related variables
  !  used as reference value for virtual flux computations
  !-----------------------------------------------------------------------

  logical (log_kind)  :: ciso_on 
  logical (log_kind)  :: vflux_flag(ecosys_tracer_cnt)     ! which tracers get virtual fluxes applied
  real (r8)           :: surf_avg(ecosys_tracer_cnt)       ! average surface tracer values
  logical (log_kind)  :: ciso_vflux_flag(ecosys_ciso_tracer_cnt)
  real (r8)           :: ciso_surf_avg(ecosys_ciso_tracer_cnt)  ! average surface tracer values

  !-----------------------------------------------------------------------
  ! needed as a module variable because of interface to ecosys_write_restart
  !-----------------------------------------------------------------------

  type(ind_name_pair) :: ind_name_table(ecosys_tracer_cnt) ! derived type & parameter for tracer index lookup
  type(ind_name_pair) :: ciso_ind_name_table(ecosys_ciso_tracer_cnt)

  !-----------------------------------------------------------------------
  !  PAR variable for each thread
  !  FIXME(ktl) move to appropriate marbl derived type and allocate and initialize it in marbl_init
  !-----------------------------------------------------------------------

  type(photosynthetically_available_radiation_type) :: PAR_instances(max_blocks_clinic)

  !-----------------------------------------------------------------------
  !  named field indices
  !-----------------------------------------------------------------------

  integer (int_kind) :: totChl_surf_nf_ind = 0 ! total chlorophyll in surface layer 
  integer (int_kind) :: sflux_co2_nf_ind   = 0 ! air-sea co2 gas flux 
  integer (int_kind) :: atm_co2_nf_ind     = 0 ! atmospheric co2

  !-----------------------------------------------------------------------
  !  module variables related to ph computations
  !-----------------------------------------------------------------------

  real (r8), allocatable, target :: iron_patch_flux(:, :, :) ! localized iron patch flux

  !-----------------------------------------------------------------------
  !  restoring climatologies for nutrients
  !-----------------------------------------------------------------------

  real (r8), allocatable, target :: fesedflux(:, :, :, :)      !  sedimentary Fe inputs

  !-----------------------------------------------------------------------
  !  define array for holding flux-related quantities that need to be time-averaged
  !-----------------------------------------------------------------------

  integer (int_kind) :: num_elements  = nx_block !TODO - make this an input namelist value

  type(ecosys_surface_share_type) :: surface_share

  !-----------------------------------------------------------------------
  !  ciso_data_ind_d13c is the index for the D13C data for the
  !  current timestep
  !  Note that ciso_data_ind_d13c is always less than ciso_atm_d13c_data_nbval.
  !  To enable OpenMP parallelism, duplicating data_ind for each block
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(max_blocks_clinic) :: ciso_data_ind_d13c = -1 ! data index for D13C data
  integer (int_kind), dimension(max_blocks_clinic) :: ciso_data_ind_d14c = -1 ! data index for D14C data

  real (r8) :: D13C(nx_block, ny_block, max_blocks_clinic)  ! atm 13co2 value
  real (r8) :: D14C(nx_block, ny_block, max_blocks_clinic)  ! atm 14co2 value
  real (r8) :: D14C_glo_avg  ! global average D14C over the ocean, computed from current D14C field

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_driver_tracer_cnt_init(ciso_active_flag)

    ! !DESCRIPTION:
    !  Zero-level initialization of ecosys_driver,
    !  which involves setting the ecosys_driver_tracer_cnt

    ! !INPUT PARAMETERS:

    logical (kind=log_kind), intent(in)  ::  &
         ciso_active_flag             ! ecosys ciso is on
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    ! Determine module variable ciso_on
    !-----------------------------------------------------------------------

    if (ciso_active_flag) then
       ciso_on = .true.
    else
       ciso_on = .false.
    end if

    !-----------------------------------------------------------------------
    !  Determine ecosys_driver_tracer_cnt, depending on whether only ecosys
    !  or also other modules are on
    !-----------------------------------------------------------------------

    ecosys_driver_tracer_cnt = ecosys_tracer_cnt

    if (ciso_on) then
       ecosys_driver_tracer_cnt = ecosys_driver_tracer_cnt + ecosys_ciso_tracer_cnt
    end if

  end subroutine ecosys_driver_tracer_cnt_init

  !***********************************************************************

  subroutine ecosys_driver_init( &
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

    use marbl_namelist_mod    , only : marbl_nl_in_size
    use marbl_namelist_mod    , only : marbl_nl_cnt
    use marbl_namelist_mod    , only : marbl_nl_buffer_size
    use marbl_namelist_mod    , only : marbl_nl_split_string
    use marbl_namelist_mod    , only : marbl_namelist
    use marbl_interface_types , only : tracer_field => marbl_tracer_metadata_type
    use marbl_interface_types , only : tracer_read  => marbl_tracer_read_type
    use marbl_share_mod       , only : autotrophs
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : comp_surf_avg_flag 
    use marbl_share_mod       , only : comp_surf_avg_freq_iopt
    use marbl_share_mod       , only : comp_surf_avg_freq
    use marbl_share_mod       , only : ciso_comp_surf_avg_flag 
    use marbl_share_mod       , only : ciso_comp_surf_avg_freq_iopt
    use marbl_share_mod       , only : ciso_comp_surf_avg_freq
    use grid                  , only : region_mask
    use broadcast             , only : broadcast_scalar
    use io_tools              , only : document
    use time_management       , only : init_time_flag
    use passive_tracer_tools  , only : set_tracer_indices
    use passive_tracer_tools  , only : read_field
    use communicate           , only : my_task, master_task
    use named_field_mod       , only : named_field_register
    use named_field_mod       , only : named_field_set
    use prognostic            , only : curtime
    use prognostic            , only : oldtime


    character (*)        , intent(in)    :: init_ts_file_fmt      ! format (bin or nc) for input file
    character (*)        , intent(in)    :: read_restart_filename ! file name for restart file
    type (tracer_field)  , intent(inout) :: tracer_d_module(:)    ! descriptors for each tracer
    real (r8)            , intent(inout) :: tracer_module(:,:,:,:,:,:)
    character (char_len) , intent(out)   :: tadvect_ctype(:)      ! advection method for ecosys tracers
    integer (POP_i4)     , intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter             :: subname = 'ecosys_driver:ecosys_driver_init'
    integer (int_kind)                  :: cumulative_nt, n, bid, k, i, j
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
    real(r8)                            :: work(nx_block, ny_block)           
    real (r8)                           :: surface_vals(ecosys_tracer_cnt)
    character (char_len)                :: ecosys_restart_filename            ! modified file name for restart file
    character (char_len)                :: init_ecosys_option                 ! namelist option for initialization of bgc
    character (char_len)                :: init_ecosys_init_file              ! filename for option 'file'
    character (char_len)                :: init_ecosys_init_file_fmt          ! file format for option 'file'
    type(tracer_read)                   :: tracer_init_ext(ecosys_tracer_cnt) ! namelist variable for initializing tracers
    type(marbl_sizes_type)              :: marbl_sizes

    !-----------------------------------------------------------------------
    !  read in ecosys_driver namelist, to set namelist parameters that
    !  should be the same for all ecosystem-related modules
    !-----------------------------------------------------------------------

    namelist /ecosys_driver_nml/ &
         lmarginal_seas, ecosys_tadvect_ctype, ecosys_qsw_distrb_const

    errorCode = POP_Success

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
       call exit_POP(sigAbort, 'ERROR reading pop namelist file into buffer.')
    endif

    ! broadcast the namelist string
    call broadcast_scalar(nl_str, master_task)
    ! process namelist string to store in nl_buffer
    call marbl_nl_split_string(nl_str, nl_buffer)

    ! now every process reads the namelist from the buffer
    ioerror_msg=''
    tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_driver_nml')
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
    call get_timer(ecosys_pre_sflux_timer    , 'ECOSYS_PRE_SFLUX'  , 1              , distrb_clinic%nprocs)
    call get_timer(ecosys_comp_CO3terms_timer, 'comp_CO3terms'     , nblocks_clinic , distrb_clinic%nprocs)
    if (ndep_data_type == 'shr_stream') then
       call get_timer(ecosys_shr_strdata_advance_timer, 'ecosys_shr_strdata_advance', 1, distrb_clinic%nprocs)
    endif
    if (ciso_on) then
       call get_timer(ecosys_ciso_interior_timer , 'ECOSYS_CISO_INTERIOR' , nblocks_clinic, distrb_clinic%nprocs)
       call get_timer(ecosys_ciso_set_sflux_timer, 'ECOSYS_CISO_SET_SFLUX', 1             , distrb_clinic%nprocs)
    end if

    !-----------------------------------------------------------------------
    !  Set up indices for ecosys modules that are on
    !-----------------------------------------------------------------------

    !  These indices are relative to ecosys_driver_ind_begin/end
    !  This means they start at 1 and end at ecosys_driver_tracer_cnt
    !  ecosys_driver then passes the ecosys module tracers back to passive
    !  tracers within TRACER(:,:,:,ecosys_driver_ind_beg,ecosys_driver_ind_end)

    cumulative_nt = 0

    call set_tracer_indices('ECOSYS', ecosys_tracer_cnt, cumulative_nt,  &
         ecosys_ind_begin, ecosys_ind_end)

    if (ciso_on) then
       call set_tracer_indices('CISO', ecosys_ciso_tracer_cnt, cumulative_nt,  &
            ecosys_ciso_ind_begin, ecosys_ciso_ind_end)
    end if

    if (cumulative_nt /= ecosys_driver_tracer_cnt) then
       call document(subname, 'ecosys_driver_tracer_cnt', ecosys_driver_tracer_cnt)
       call document(subname, 'cumulative_nt', cumulative_nt)
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: ecosys_driver_tracer_cnt does not match cumulative nt')
    end if

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
    !  Initialize marbl 
    !--------------------------------------------------------------------

    do iblock=1, nblocks_clinic
       call marbl(iblock)%init(                                &
            ciso_on,                                           &
            nl_buffer,                                         &
            num_elements_interior=1,                           & 
            num_elements_forcing=num_elements,                 & 
            num_levels=km,                                     & 
            marbl_tracer_metadata=tracer_d_module,             &
            marbl_sizes=marbl_sizes,                           &         
            marbl_interior_diags=marbl_interior_diags(iblock), &
            marbl_restore_diags=marbl_restore_diags(iblock),   &
            marbl_forcing_diags=marbl_forcing_diags(iblock),   &
            marbl_forcing_input=marbl_forcing_input(iblock),   &
            marbl_forcing_output=marbl_forcing_output(iblock), &
            marbl_forcing_share=marbl_forcing_share(iblock))

       ! Check marbl(iblock)%StatusLog, abort on error
       ! This might require a global collection of labort_marbl, abort
       ! if any are true
    end do

    !--------------------------------------------------------------------
    !  Determine advection type
    !--------------------------------------------------------------------

    tadvect_ctype(1:ecosys_driver_tracer_cnt) = ecosys_tadvect_ctype

    !--------------------------------------------------------------------
    !  Initialize surface average times
    !--------------------------------------------------------------------

    call init_time_flag('ecosys_comp_surf_avg', &
         comp_surf_avg_flag,                    &
         default=.false.,                       &
         freq_opt=comp_surf_avg_freq_iopt,      &
         freq=comp_surf_avg_freq,               &
         owner='ecosys_init')

    if (ciso_on) then
       call init_time_flag('ciso_ecosys_comp_surf_avg', &
            ciso_comp_surf_avg_flag,                    &
            default=.false.,                            &
            freq_opt=ciso_comp_surf_avg_freq_iopt,      &
            freq=ciso_comp_surf_avg_freq,               &
            owner='ciso_ecosys_init')
    end if

    !--------------------------------------------------------------------
    !  Initialize module variables ind_name_table and ciso_ind_name_table
    !--------------------------------------------------------------------

    ! NOTE (mvertens, 2015-11), this currently has to be a module variable 
    ! since ecosys_driver_write_restart is called by passive_tracers and 
    ! does not pass in tracer_d_module into the interface

    do n = 1, ecosys_tracer_cnt
       ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
    end do

    if (ciso_on) then
       do n = 1, ecosys_ciso_tracer_cnt
          ciso_ind_name_table(n) = ind_name_pair(n, tracer_d_module(ecosys_ciso_ind_begin+n-1)%short_name)
       end do
    end if
       
    !--------------------------------------------------------------------
    !  Initialize ecosys tracers
    !--------------------------------------------------------------------

    surf_avg(:) = 0._r8

    call ecosys_driver_init_tracers_and_saved_state( &
       init_ts_file_fmt,                             &
       read_restart_filename,                        &
       tracer_d_module(1:ecosys_ind_end),            &
       tracer_module(:,:,:,1:ecosys_ind_end,:,:),    &
       ecosys_restart_filename,                      &
       errorCode)       

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_driver_init')
       return
    endif

    if (ciso_on) then
       call ecosys_driver_ciso_init_tracers(                                    &
            init_ts_file_fmt,                                                   &
            read_restart_filename,                                              &
            tracer_d_module(ecosys_ciso_ind_begin:ecosys_ciso_ind_end),         &
            tracer_module(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:,:), &
            errorCode)

       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_ciso_init')
          return
       endif

       call ecosys_driver_ciso_init_atm_D13_D14()
    end if

    !--------------------------------------------------------------------
    !  Initialize surface forcing flux
    !--------------------------------------------------------------------

    call ecosys_driver_init_sflux()   

    !--------------------------------------------------------------------
    !  Initialize interior restoring
    !--------------------------------------------------------------------

    call ecosys_restore%init(nml_filename, nml_in, ind_name_table)
    call ecosys_driver_init_interior_restore(ecosys_restore)

    !--------------------------------------------------------------------
    ! Initialize tavg ids (need only do this using first block)
    !--------------------------------------------------------------------

    call ecosys_tavg_init(marbl_interior_diags(1), marbl_restore_diags(1), marbl_forcing_diags(1))

    !--------------------------------------------------------------------
    ! Register Chl field for short-wave absorption
    !--------------------------------------------------------------------

    !  apply land mask to tracers
    !  set Chl field for short-wave absorption

    call named_field_register('model_chlorophyll', totChl_surf_nf_ind)
    !$OMP PARALLEL DO PRIVATE(iblock, work, surface_vals, i, j)
    do iblock=1, nblocks_clinic
       do i = 1, nx_block
         do j = 1, ny_block
            surface_vals(1:ecosys_tracer_cnt) =p5*( &
                 tracer_module(i, j, 1, 1:ecosys_ind_end, oldtime, iblock) +   &
                 tracer_module(i, j, 1, 1:ecosys_ind_end, curtime, iblock))

           work(i,j) = marbl_ecosys_compute_totalChl( tracer_in=surface_vals(:), nb=1, ne=ecosys_ind_end )
         end do
       end do
       call named_field_set(totChl_surf_nf_ind, iblock, WORK)
    enddo
    !$OMP END PARALLEL DO

  end subroutine ecosys_driver_init

  !-----------------------------------------------------------------------

  subroutine ecosys_driver_init_tracers_and_saved_state(&
       init_ts_file_fmt, read_restart_filename, tracer_d_module, TRACER_MODULE, &
       ecosys_restart_filename, errorCode)       

    use marbl_parms           , only : dic_ind
    use marbl_parms           , only : alk_ind
    use marbl_parms           , only : dic_alt_co2_ind
    use marbl_parms           , only : di13c_ind
    use marbl_parms           , only : di14c_ind
    use marbl_share_mod       , only : init_ecosys_option
    use marbl_share_mod       , only : init_ecosys_init_file
    use marbl_share_mod       , only : init_ecosys_init_file_fmt
    use marbl_share_mod       , only : use_nml_surf_vals         
    use marbl_share_mod       , only : tracer_init_ext
    use marbl_share_mod       , only : surf_avg_dic_const
    use marbl_share_mod       , only : surf_avg_alk_const
    use marbl_share_mod       , only : comp_surf_avg_flag
    use passive_tracer_tools  , only : extract_surf_avg
    use passive_tracer_tools  , only : comp_surf_avg
    use passive_tracer_tools  , only : rest_read_tracer_block
    use passive_tracer_tools  , only : file_read_tracer_block
    use passive_tracer_tools  , only : field_exists_in_file
    use passive_tracer_tools  , only : read_field
    use passive_tracer_tools  , only : ind_name_pair
    use marbl_interface_types , only : tracer_field => marbl_tracer_metadata_type
    use prognostic            , only : curtime
    use prognostic            , only : oldtime
    use prognostic            , only : newtime
    use time_management       , only : check_time_flag
    use time_management       , only : eval_time_flag
    use ecosys_constants      , only : ecosys_tracer_cnt
    use io_tools              , only : document
    use grid                  , only : fill_points
    use grid                  , only : n_topo_smooth
    use grid                  , only : KMT
    use exit_mod              , only : exit_POP

    implicit none
    
    ! !INPUT PARAMETERS:

    character (*)                , intent(in)    :: init_ts_file_fmt                   ! format (bin or nc) for input file
    character (*)                , intent(in)    :: read_restart_filename              ! file name for restart file
    type(tracer_field)           , intent(in)    :: tracer_d_module(:)                 ! descriptors for each tracer

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8)                    , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)

    ! !OUTPUT PARAMETERS:

    character(char_len)          , intent(out)   :: ecosys_restart_filename     ! modified file name for restart file
    integer (POP_i4)             , intent(out)   :: errorCode
    
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = &
        'ecosys_driver_mod:ecosys_driver_init_tracers_and_saved_state'
    integer :: n, k, bid

    !-----------------------------------------------------------------------
    !  allocate memory for saved state
    !-----------------------------------------------------------------------
    call ecosys_saved_state%construct()

    !-----------------------------------------------------------------------
    !  initialize tracers
    !-----------------------------------------------------------------------

    ! initialize module variable - virtual flux flag array
    ! FIXME(mnl,2016-01): do these really belong here, or should they be in
    !                     MARBL? Decision needs to be made
    vflux_flag(:) = .false.
    vflux_flag(dic_ind) = .true.
    vflux_flag(alk_ind) = .true.
    vflux_flag(dic_alt_co2_ind) = .true.

    if (ecosys_ciso_tracer_cnt > 0) then
       ciso_vflux_flag(:) = .false.
       ciso_vflux_flag(di13c_ind) = .true.
       ciso_vflux_flag(di14c_ind) = .true.
    end if

    !-----------------------------------------------------------------------
    !  allocate components in PAR derived type
    !  FIXME(ktl) eventually is allocation for a single instance
    !  FIXME(ktl) use PAR_nsubcols instead of MCOG_nbins
    !  FIXME(mnl) not tracer related, should not be in
    !             ecosys_driver_init_tracers_and_saved_state
    !-----------------------------------------------------------------------

    do bid=1, nblocks_clinic
       allocate(PAR_instances(bid)%col_frac(mcog_nbins))
       allocate(PAR_instances(bid)%interface(0:km,mcog_nbins))
       allocate(PAR_instances(bid)%avg(km,mcog_nbins))
       allocate(PAR_instances(bid)%KPARdz(km))
    enddo

    select case (init_ecosys_option)
       
    case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

       ecosys_restart_filename = char_blank
       
       if (init_ecosys_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call document(subname, 'no restart file to read ecosys from')
             call exit_POP(sigAbort, 'stopping in ' // subname)
          endif
          ecosys_restart_filename = read_restart_filename
          init_ecosys_init_file_fmt = init_ts_file_fmt
          
       else  ! do not read from TS restart file
          
          ecosys_restart_filename = trim(init_ecosys_init_file)
          
       endif
       
       call rest_read_tracer_block(init_ecosys_init_file_fmt, &
                                   ecosys_restart_filename,   &
                                   tracer_d_module,           &
                                   TRACER_MODULE)
       
       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_SURF')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_SURF', ecosys_saved_state%ph_surf)
       else
          call document(subname, 'PH_SURF does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_SURF to 0')
          ecosys_saved_state%ph_surf = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_SURF_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_SURF_ALT_CO2', ecosys_saved_state%ph_surf_alt_co2)
       else
          call document(subname, 'PH_SURF_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2 to 0')
          ecosys_saved_state%ph_surf_alt_co2 = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_3D')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D', ecosys_saved_state%PH_PREV_3D)
       else
          call document(subname, 'PH_3D does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_3D to 0')
          ecosys_saved_state%PH_PREV_3D  = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_3D_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D_ALT_CO2', ecosys_saved_state%PH_PREV_ALT_CO2_3D)
       else
          call document(subname, 'PH_3D_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2_3D to 0')
          ecosys_saved_state%PH_PREV_ALT_CO2_3D = c0
       endif

       if (use_nml_surf_vals) then
          surf_avg(:)               = c0
          surf_avg(dic_ind)         = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind)         = surf_avg_alk_const
       else
          call extract_surf_avg(&
               init_ecosys_init_file_fmt,     &
               ecosys_restart_filename,       &
               ecosys_tracer_cnt, vflux_flag, &
               ind_name_table, surf_avg)
       endif

       ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do
       call eval_time_flag(comp_surf_avg_flag) 

       if (check_time_flag(comp_surf_avg_flag)) then
          call comp_surf_avg(&
               TRACER_MODULE(:, :, 1, :, oldtime, :), &
               TRACER_MODULE(:, :, 1, :, curtime, :), &
               ecosys_tracer_cnt, vflux_flag, surf_avg)
       end if

    case ('file', 'ccsm_startup')
       call document(subname, 'ecosystem vars being read from separate files')

       call file_read_tracer_block(init_ecosys_init_file_fmt, &
                                   init_ecosys_init_file,     &
                                   tracer_d_module,           &
                                   ind_name_table,            &
                                   tracer_init_ext,           &
                                   TRACER_MODULE)

       if (n_topo_smooth > 0) then
          do n = 1, ecosys_tracer_cnt
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
          enddo
       endif

       ecosys_saved_state%ph_surf = c0
       ecosys_saved_state%ph_surf_alt_co2 = c0
       ecosys_saved_state%ph_prev_3d = c0
       ecosys_saved_state%ph_prev_alt_co2_3d = c0

       if (use_nml_surf_vals) then
          surf_avg                  = c0
          surf_avg(dic_ind)         = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind)         = surf_avg_alk_const
       else
          call comp_surf_avg(&
               TRACER_MODULE(:, :, 1, :, oldtime, :), &
               TRACER_MODULE(:, :, 1, :, curtime, :), &
               ecosys_tracer_cnt, vflux_flag, surf_avg)
       endif

    case default
       call document(subname, 'init_ecosys_option', init_ecosys_option)
       call exit_POP(sigAbort, 'unknown init_ecosys_option')

    end select

    do bid=1, nblocks_clinic
       do n = 1, ecosys_tracer_cnt
          do k = 1, km
             where (.not. land_mask(:, :, bid) .or. k > KMT(:, :, bid))
                TRACER_MODULE(:, :, k, n, curtime, bid) = c0
                TRACER_MODULE(:, :, k, n, oldtime, bid) = c0
             end where
          end do
       end do
    end do

  end subroutine ecosys_driver_init_tracers_and_saved_state

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

    use constants             , only : salt_to_ppt
    use grid                  , only : KMT
    use grid                  , only : DZT
    use grid                  , only : dz
    use grid                  , only : partial_bottom_cells
    use grid                  , only : zt
    use grid                  , only : zw
    use marbl_interface_types , only : marbl_column_domain_type
    use marbl_interface_types , only : marbl_gcm_state_type

    ! !INPUT PARAMETERS:

    real (r8), dimension(nx_block, ny_block, mcog_nbins), intent(in) :: &
         FRACR_BIN,    &! fraction of cell occupied by mcog bin
         QSW_RAW_BIN,  &! raw (directly from cpl) shortwave into each mcog column (W/m^2)
         QSW_BIN        ! shortwave into each mcog bin, potentially modified by coszen factor (W/m^2)

    real (r8), dimension(nx_block, ny_block, km), intent(in) :: &
         TEMP_OLD,          &! old potential temperature (C)
         TEMP_CUR,          &! current potential temperature (C)
         SALT_OLD,          &! old salinity (msu)
         SALT_CUR            ! current salinity (msu)

    real (r8), dimension(:,:,:,:), intent(in) :: &
         TRACER_MODULE_OLD, &! old tracer values
         TRACER_MODULE_CUR   ! current tracer values

    ! !OUTPUT PARAMETERS:

    real (r8), dimension(:,:,:,:), intent(inout) :: &
         DTRACER_MODULE      ! computed source/sink terms

    type (block), intent(in) :: this_block   ! block information for this block

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: i ! nx_block loop index
    integer (int_kind) :: c ! ny_block / column loop index
    integer (int_kind) :: k ! vertical level index
    integer (int_kind) :: n, d, ncols
    integer (int_kind) :: bid ! local block address for this block

    type(marbl_column_domain_type) :: marbl_domain
    type(marbl_gcm_state_type)     :: marbl_gcm_state

    real (r8) :: column_tracer_module(ecosys_driver_tracer_cnt, km)
    real (r8) :: column_dtracer(ecosys_driver_tracer_cnt, km)
    real (r8) :: column_restore_local(ecosys_tracer_cnt, km) ! local restoring terms for nutrients (mmol ./m^3/sec)
    real (r8) :: column_dust_flux_in
    real (r8) :: column_fesedflux(km)
    real (r8) :: column_ph_prev_3d(km)
    real (r8) :: column_ph_prev_alt_co2_3d(km)

    !-----------------------------------------------------------------------

    bid = this_block%local_id

    ! FIXME(bja, 2015-07) one time copy of global marbl_domain
    ! related memory from slab to column ordering. move entire
    ! copy to ecosys_driver_init.
    marbl_domain%PAR_nsubcols = mcog_nbins

    marbl_domain%km = km
    allocate(marbl_domain%dz(marbl_domain%km))
    allocate(marbl_domain%delta_z(marbl_domain%km))
    allocate(marbl_domain%zw(marbl_domain%km))
    allocate(marbl_domain%zt(marbl_domain%km))
    allocate(marbl_gcm_state%temperature(marbl_domain%km))
    allocate(marbl_gcm_state%salinity(marbl_domain%km))
    allocate(marbl_domain%PAR_col_frac(marbl_domain%PAR_nsubcols))
    allocate(marbl_domain%surf_shortwave(marbl_domain%PAR_nsubcols))

    ! NTOE: gcm dependent quantities (i.e. time stepping). need to be
    ! averaged into a single quantity for marbl ecosys

    call timer_start(ecosys_interior_timer, block_id=bid)

    do c = this_block%jb,this_block%je
       do i = this_block%ib,this_block%ie

          marbl_domain%land_mask = land_mask(i, c, bid)
          marbl_domain%kmt = KMT(i, c, bid)

          if (marbl_domain%land_mask) then

             !-----------------------------------------------------------
             ! copy data from slab to column for marbl
             !-----------------------------------------------------------

             do ncols = 1,marbl_domain%PAR_nsubcols
                marbl_domain%PAR_col_frac(ncols) = FRACR_BIN(i, c, ncols)

                ! select short-wave forcing
                if (ecosys_qsw_distrb_const) then
                   marbl_domain%surf_shortwave(ncols) = QSW_RAW_BIN(i, c, ncols)
                else
                   marbl_domain%surf_shortwave(ncols) = QSW_BIN(i, c, ncols)
                end if
             end do

             do k = 1, marbl_domain%km
                marbl_domain%dz(k) = dz(k)
                if (partial_bottom_cells) then
                   marbl_domain%delta_z(k) = DZT(i, c, k, bid)
                else
                   marbl_domain%delta_z(k) = dz(k)
                end if
                marbl_domain%zw(k) = zw(k)
                marbl_domain%zt(k) = zt(k)

                marbl_gcm_state%temperature(k) = p5*(TEMP_OLD(i, c, k) + TEMP_CUR(i, c, k))

                marbl_gcm_state%salinity(k)    = p5*(SALT_OLD(i, c, k) + SALT_CUR(i, c, k))*salt_to_ppt

                do n = 1, ecosys_ind_end
                   column_tracer_module(n, k) = p5*(TRACER_MODULE_OLD(i, c, k, n) + TRACER_MODULE_CUR(i, c, k, n))
                end do
                if (ciso_on) then
                   do n = ecosys_ciso_ind_begin, ecosys_ciso_ind_end
                      column_tracer_module(n, k) = p5*(TRACER_MODULE_OLD(i, c, k, n) + TRACER_MODULE_CUR(i, c, k, n))
                   end do
                end if

             end do ! do k

             if (marbl_domain%kmt > 0) then 

                !  set tracer restore fields
                do k = 1, marbl_domain%km
                   call ecosys_restore%restore_tracers(ecosys_tracer_cnt,     &
                        vert_level=k, x_index=i, y_index=c, block_id=bid,     &
                        local_data=column_tracer_module(:,k),                 &
                        restore_data=column_restore_local(:, k))
                end do

                column_dust_flux_in  = ecosys_saved_state%dust_FLUX_IN(i, c, bid)
                do k = 1, marbl_domain%km
                   column_ph_prev_3d(k)         = ecosys_saved_state%ph_prev_3d(i, c, k, bid)
                   column_ph_prev_alt_co2_3d(k) = ecosys_saved_state%ph_prev_alt_co2_3d(i, c, k, bid)
                   ! FIXME(mnl,2016-01) column_fesedflux could be in interior forcing field datatype
                   column_fesedflux(k)          = fesedflux(i, c, k, bid)
                end do

                !-----------------------------------------------------------
                !  compute time derivatives for ecosystem state variables
                !-----------------------------------------------------------

                call marbl_ecosys_set_interior( &
                     ciso_on,                   &
                     marbl_domain,              &
                     marbl_gcm_state,           &
                     marbl_interior_diags(bid), &
                     marbl_restore_diags(bid),  &
                     column_restore_local,      &
                     PAR_instances(bid),        &
                     column_fesedflux,          &
                     column_dust_flux_in,       &
                     column_tracer_module,      &
                     column_ph_prev_3d,         &
                     column_ph_prev_alt_co2_3d, &
                     column_dtracer)

                !-----------------------------------------------------------
                ! copy marbl column data back to slab
                !-----------------------------------------------------------

                do k = 1, marbl_domain%km
                   ecosys_saved_state%ph_prev_3d(i, c, k, bid)         = column_ph_prev_3d(k)
                   ecosys_saved_state%ph_prev_alt_co2_3d(i, c, k, bid) = column_ph_prev_alt_co2_3d(k)

                   do n = 1, ecosys_ind_end
                      dtracer_module(i, c, k, n) = column_dtracer(n, k)
                   end do 
                   if (ciso_on) then
                      do n = ecosys_ciso_ind_begin, ecosys_ciso_ind_end
                         dtracer_module(i, c, k, n) = column_dtracer(n, k)
                      end do
                   end if
                end do ! do k

             end if ! end if marbl_domain%kmt > 0

          end if ! end if marbl_domain%land_mask > 0

          call ecosys_tavg_accumulate((/i/), (/c/), bid,                     &
               marbl_interior_diags=marbl_interior_diags(bid),               &
               marbl_restore_diags=marbl_restore_diags(bid))

       end do ! do i
    end do ! do c
    
    call timer_stop(ecosys_interior_timer, block_id=bid)

    deallocate(marbl_domain%delta_z)
    deallocate(marbl_domain%dz)
    deallocate(marbl_domain%zw)
    deallocate(marbl_domain%zt)
    deallocate(marbl_gcm_state%temperature)
    deallocate(marbl_gcm_state%salinity)
    deallocate(marbl_domain%PAR_col_frac)
    deallocate(marbl_domain%surf_shortwave)

  end subroutine ecosys_driver_set_interior

  !***********************************************************************

  subroutine ecosys_driver_set_sflux( &
       U10_SQR,IFRAC,PRESS,SST,SSS,   &
       SURFACE_VALS_OLD,              &
       SURFACE_VALS_CUR,              &
       STF_MODULE)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute surface fluxes

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use constants            , only : xkw_coeff
    use named_field_mod      , only : named_field_set
    use named_field_mod      , only : named_field_set
    use marbl_parms          , only : parm_Fe_bioavail
    use marbl_share_mod      , only : gas_flux_forcing_iopt
    use marbl_share_mod      , only : gas_flux_forcing_iopt_drv
    use marbl_share_mod      , only : gas_flux_forcing_iopt_file
    use marbl_share_mod      , only : lflux_gas_co2
    use marbl_share_mod      , only : lflux_gas_o2
    use marbl_share_mod      , only : autotrophs
    use marbl_share_mod      , only : comp_surf_avg_flag
    use marbl_share_mod      , only : ciso_comp_surf_avg_flag 
    use ecosys_mod           , only : marbl_set_sflux
    use marbl_ciso_mod       , only : marbl_ciso_set_sflux
    use passive_tracer_tools , only : comp_surf_avg
    use time_management      , only : check_time_flag
    use domain               , only : nblocks_clinic

    ! !INPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
         U10_SQR,      &! 10m wind speed squared (cm/s)**2
         IFRAC,        &! sea ice fraction (non-dimensional)
         PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
         SST,          &! sea surface temperature (C)
         SSS            ! sea surface salinity (psu)

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8), dimension(:,:,:,:), intent(inout) :: &
         STF_MODULE

    ! OUTPUT PARAMETERS:

    real (r8), dimension(:,:,:,:), intent(in) :: &
         surface_vals_old, &
         surface_vals_cur ! module tracers

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: &
         i, j, iblock, n, & ! loop indices
         auto_ind,        & ! autotroph functional group index
         errorCode          ! errorCode from HaloUpdate call

    real (r8)             , dimension(nx_block, ny_block, max_blocks_clinic) :: &
         ifrac_file_input , & ! file input ice fraction (non-dimensional)
         xkw_file_input   , & ! file input portion of piston velocity (cm/s)
         ap_file_input    , & ! used atm pressure (converted from dyne/cm**2 to atm)
         iron_flux_in     , & ! iron flux
         ifrac_used       , & ! used ice fraction (non-dimensional)
         xco2             , & ! atmospheric co2 conc. (dry-air, 1 atm)
         xco2_alt_co2     , & ! atmospheric alternative CO2 (dry-air, 1 atm)
         flux_co2

    real (r8), dimension(nx_block, ny_block) :: &
         work1 ! temporaries for averages

    real (r8), dimension(nx_block, ny_block, 13, max_blocks_clinic) :: &
         input_forcing_data !FIXME - don't hardwire 13

    real (r8), dimension(ecosys_driver_tracer_cnt) :: &
         surface_vals
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Set surf_avg (module variable)
    !-----------------------------------------------------------------------

    if (check_time_flag(comp_surf_avg_flag))  then
       call comp_surf_avg( &
            surface_vals_old(:, :, 1:ecosys_ind_end, :), &
            surface_vals_cur(:, :, 1:ecosys_ind_end, :), &
            ecosys_tracer_cnt, vflux_flag, surf_avg)
    end if

    if (ciso_on) then
       if (check_time_flag(ciso_comp_surf_avg_flag)) then
          call comp_surf_avg( &
               surface_vals_old(:, :, ecosys_ciso_ind_begin:ecosys_ciso_ind_end, :), &
               surface_vals_cur(:, :, ecosys_ciso_ind_begin:ecosys_ciso_ind_end, :), &
               ecosys_ciso_tracer_cnt, ciso_vflux_flag, ciso_surf_avg)
       end if
    end if

    !-----------------------------------------------------------------------
    ! Set input surface forcing data
    !-----------------------------------------------------------------------

    call ecosys_driver_read_sflux(                                &
         input_forcing_data,                                      &
         XCO2, XCO2_ALT_CO2,                                      &
         ifrac_file_input, xkw_file_input, ap_file_input, iron_flux_in)

    !-----------------------------------------------------------------------
    ! Apply OCMIP ice fraction mask when input is from a file.
    !-----------------------------------------------------------------------

    ! FIXME(mnl,2016-01): ecosys_mod::ecosys_update_surface_forcing_fields()
    !                     that applies this fix and also multiplies iron_flux_in
    !                     by parm_Fe_bioavail ; similar multiplication of
    !                     dust_flux by 0.98 in read_sflux should also be moved

    if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
       ifrac_used = ifrac_file_input
       where (ifrac_used < 0.2000_r8) ifrac_used = 0.2000_r8
       where (ifrac_used > 0.9999_r8) ifrac_used = 0.9999_r8
    else if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv) then
       ifrac_used = ifrac
       where (ifrac_used < c0) ifrac_used = c0
       where (ifrac_used > c1) ifrac_used = c1
    endif

    !-----------------------------------------------------------------------
    ! Update carbon isotope atmosphere deltas
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ecosys_driver_ciso_update_atm_D13C_D14C()
    end if

    !-----------------------------------------------------------------------
    ! Set surface flux 
    !-----------------------------------------------------------------------

    call timer_start(ecosys_set_sflux_timer)

    do iblock = 1, nblocks_clinic
       do j = 1, ny_block

          !-----------------------------------------------------------------------
          ! Copy data from slab data structure to column input for marbl
          !-----------------------------------------------------------------------

          marbl_forcing_input(iblock)%u10_sqr(:)         = u10_sqr(:,j,iblock)
          marbl_forcing_input(iblock)%sst(:)             = sst(:,j,iblock)
          marbl_forcing_input(iblock)%sss(:)             = sss(:,j,iblock)
          marbl_forcing_input(iblock)%xco2(:)            = xco2(:,j,iblock)
          marbl_forcing_input(iblock)%xco2_alt_co2(:)    = xco2_alt_co2(:,j,iblock)
          marbl_forcing_input(iblock)%ifrac(:)           = ifrac_used(:,j,iblock)
          marbl_forcing_input(iblock)%ph_prev(:)         = ecosys_saved_state%ph_surf(:,j,iblock)         
          marbl_forcing_input(iblock)%ph_prev_alt_co2(:) = ecosys_saved_state%ph_surf_alt_co2(:,j,iblock) 
          marbl_forcing_input(iblock)%iron_flux(:)       = iron_flux_in(:,j,iblock)
          marbl_forcing_input(iblock)%dust_flux(:)       = ecosys_saved_state%dust_flux_in(:,j,iblock)
          marbl_forcing_input(iblock)%land_mask(:)       = land_mask(:,j,iblock)

          if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv) then
             marbl_forcing_input(iblock)%xkw(:) = xkw_coeff * u10_sqr(:,j,iblock)
          else
             marbl_forcing_input(iblock)%xkw(:) = xkw_file_input(:,j,iblock)
          end if

          if (lflux_gas_o2 .or. lflux_gas_co2) then
             !  assume PRESS is in cgs units (dyne/cm**2) since that is what is
             !    required for pressure forcing in barotropic
             !  want units to be atmospheres
             !  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
             !  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8

             marbl_forcing_input(iblock)%atm_press(:) = press(:,j,iblock) / 101.325e+4_r8
          end if

          marbl_forcing_input(iblock)%input_forcings(:,:) = input_forcing_data(:,j,:,iblock)

          if (ciso_on) then
             marbl_forcing_input(iblock)%d13c(:) = D13C(:,j,iblock)
             marbl_forcing_input(iblock)%d14c(:) = D14C(:,j,iblock)
          end if

          do n = 1,ecosys_driver_tracer_cnt
             marbl_forcing_input(iblock)%surface_vals(:,n) = &
                  p5*(surface_vals_old(:, j, n, iblock) &
                    + surface_vals_cur(:, j, n, iblock))
          end do


          !-----------------------------------------------------------------------
          ! Determine surface flux - marbl
          !-----------------------------------------------------------------------

          ! FIXME - remove marbl_forcing_share from below and call marbl_ciso_set_sflux
          ! from marbl_set_sflux

          call marbl_set_sflux(              &
               ciso_on,                      &
               num_elements,                 &
               marbl_forcing_input(iblock),  &
               marbl_forcing_output(iblock), &
               marbl_forcing_share(iblock),  &
               marbl_forcing_diags(iblock))

          if (ciso_on) then
             call marbl_ciso_set_sflux(                                                                  &
                  num_elements,                                                                          &
                  ecosys_ciso_tracer_cnt,                                                                &
                  marbl_forcing_input(iblock)%land_mask,                                                 &
                  marbl_forcing_input(iblock)%sst,                                                       &
                  marbl_forcing_input(iblock)%d13c,                                                      &
                  marbl_forcing_input(iblock)%d14c,                                                      &
                  d14c_glo_avg,                                                                          &
                  marbl_forcing_input(iblock)%surface_vals(:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end), &
                  marbl_forcing_share(iblock),                                                           &
                  marbl_forcing_output(iblock)%stf_module(:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),  &
                  marbl_forcing_diags(iblock))
          end if

          !-----------------------------------------------------------------------
          ! Copy data from marbl output column to pop slab data structure
          !-----------------------------------------------------------------------

          ! FIXME(mnl,2016-01): should these be saved state?
          ecosys_saved_state%ph_surf(:,j,iblock)         = marbl_forcing_output(iblock)%ph_prev(:)
          ecosys_saved_state%ph_surf_alt_co2(:,j,iblock) = marbl_forcing_output(iblock)%ph_prev_alt_co2(:)

          ! FIXME(mnl,2016-01): do we need a data structure to handle this? future might include additional fluxes...
          flux_co2(:,j,iblock) = marbl_forcing_output(iblock)%flux_co2(:)

          do n=1,marbl_forcing_diags(iblock)%diag_cnt
            flux_diags(:,j,n,iblock) = marbl_forcing_diags(iblock)%diags(n)%field_2d(:)
          end do

          do n = 1,ecosys_driver_tracer_cnt
             stf_module(:,j,n,iblock) = marbl_forcing_output(iblock)%stf_module(:,n)  
          end do

       enddo
    enddo

    call timer_stop(ecosys_set_sflux_timer)

    !-----------------------------------------------------------------------
    ! Update named fields
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       do i = 1, nx_block
          do j = 1, ny_block
             surface_vals(1:ecosys_tracer_cnt) =p5*( &
                  surface_vals_old(i,j,1:ecosys_tracer_cnt,iblock) + &
                  surface_vals_cur(i,j,1:ecosys_tracer_cnt,iblock))
             work1(i,j) = marbl_ecosys_compute_totalChl(   &
                  tracer_in=surface_vals(:), nb=1, ne=ecosys_tracer_cnt)
          end do
       end do
       call named_field_set(totChl_surf_nf_ind, iblock, work1)
       
       !  set air-sea co2 gas flux named field, converting units from
       !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
       if (lflux_gas_co2) then
          call named_field_set(sflux_co2_nf_ind, iblock, 44.0e-8_r8 * flux_co2(:, :, iblock))
       end if
    end do
    
    ! Note: out of this subroutine rest of pop needs stf_module, ph_surf, named_state, diagnostics

  end subroutine ecosys_driver_set_sflux

  !***********************************************************************

  subroutine ecosys_driver_write_restart(restart_file, action)

    ! !DESCRIPTION:
    !  call restart routines for each tracer module that
    !  write fields besides the tracers themselves

    use io_types, only : datafile

    implicit none 

    character(*)            , intent(in)     :: action
    type (datafile)         , intent (inout) :: restart_file

    call ecosys_driver_write_ecosys_restart(restart_file, action)
    if (ciso_on) then
       call ecosys_driver_ciso_write_restart(restart_file, action)
    end if

  end subroutine ecosys_driver_write_restart

  !***********************************************************************

  subroutine ecosys_driver_tavg_forcing()

    ! !DESCRIPTION:
    !  accumulate common tavg fields for tracer surface fluxes

    implicit none 

    call ecosys_tavg_accumulate_flux(flux_diags, marbl_forcing_diags)

  end subroutine ecosys_driver_tavg_forcing

  !***********************************************************************

  function ecosys_driver_tracer_ref_val(ind)
    !
    ! !DESCRIPTION:
    !  return reference value for tracer with global tracer index ind
    !  this is used in virtual flux computations

    implicit none

    integer(int_kind) , intent(in) :: ind
    real(r8) :: ecosys_driver_tracer_ref_val

    !  default value for reference value is 0

    ecosys_driver_tracer_ref_val = c0

    if (ind >= 1 .and. ind <= ecosys_ind_end) then
       if (vflux_flag(ind)) then
          ecosys_driver_tracer_ref_val = surf_avg(ind)
       else
          ecosys_driver_tracer_ref_val = c0
       endif
    endif

    if (ciso_on) then
       if (ind >= ecosys_ciso_ind_begin .and. ind <= ecosys_ciso_ind_end) then
          if (ciso_vflux_flag(ind-ecosys_ciso_ind_begin+1)) then
             ecosys_driver_tracer_ref_val = ciso_surf_avg(ind-ecosys_ciso_ind_begin+1)
          else
             ecosys_driver_tracer_ref_val = c0
          endif
       end if
    end if
       
  end function ecosys_driver_tracer_ref_val

  !***********************************************************************

  subroutine ecosys_driver_unpack_source_sink_terms(source, destination)

    real(r8), dimension(:, :, :), intent(in)  :: source
    real(r8), dimension(:, :, :), intent(out) :: destination

    destination(:, :, :) = source(:, :, :)

  end subroutine ecosys_driver_unpack_source_sink_terms

  !*****************************************************************************

  subroutine ecosys_driver_write_ecosys_restart(restart_file, action)

    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files

    use domain_size           , only : nx_global
    use domain_size           , only : ny_global
    use constants             , only : field_loc_center
    use constants             , only : field_type_scalar
    use io                    , only : data_set
    use io                    , only : datafile
    use io_types              , only : io_dim
    use io_types              , only : io_field_desc
    use io_types              , only : add_attrib_file
    use io_types              , only : construct_io_dim
    use io_types              , only : construct_io_field
    use passive_tracer_tools  , only : ind_name_pair

    implicit none

    character(*)                 , intent(in)    :: action
    type (datafile)              , intent(inout) :: restart_file

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (char_len)       :: short_name   ! tracer name temporaries
    integer (int_kind)         :: n
    type (io_dim)              :: i_dim, j_dim ! dimension descriptors
    type (io_dim)              :: k_dim        ! dimension descriptor for vertical levels
    type (io_field_desc), save :: ph_surf_iodesc
    type (io_field_desc), save :: ph_surf_alt_co2_iodesc
    type (io_field_desc), save :: ph_3d_alt_co2_iodesc
    type (io_field_desc), save :: ph_3d_iodesc
    !-----------------------------------------------------------------------

    if (trim(action) == 'add_attrib_file') then

       short_name = char_blank
       do n=1, ecosys_tracer_cnt
          if (vflux_flag(n)) then   
             short_name = 'surf_avg_' // ind_name_table(n)%name
             call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
          endif
       end do

    else if (trim(action) == 'define') then

       i_dim = construct_io_dim('i', nx_global)
       j_dim = construct_io_dim('j', ny_global)
       k_dim = construct_io_dim('k', km)

       ph_surf_iodesc = construct_io_field('PH_SURF', i_dim, j_dim,                    &
            long_name  ='surface pH at current time',                                  &
            units      ='pH', grid_loc='2110',                                         &
            field_loc  = field_loc_center,                                             &
            field_type = field_type_scalar,                                            &
            d2d_array  = ecosys_saved_state%ph_surf(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', ph_surf_iodesc)

       ph_surf_alt_co2_iodesc = construct_io_field('PH_SURF_ALT_CO2', i_dim, j_dim,    &
            long_name  ='surface pH, alternate CO2, at current time',                  &
            units      ='pH', grid_loc='2110',                                         &
            field_loc  = field_loc_center,                                             &
            field_type = field_type_scalar,                                            &
            d2d_array  = ecosys_saved_state%ph_surf_alt_co2(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', ph_surf_alt_co2_iodesc)

       ph_3d_alt_co2_iodesc = construct_io_field('PH_3D_ALT_CO2', i_dim, j_dim, k_dim, &
            long_name  ='3D pH, alternate CO2, at current time',                       &
            units      ='pH', grid_loc='3111',                                         &
            field_loc  = field_loc_center,                                             &
            field_type = field_type_scalar,                                            &
            d3d_array  = ecosys_saved_state%ph_prev_alt_co2_3d(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', ph_3d_alt_co2_iodesc)

       ph_3d_iodesc = construct_io_field('PH_3D', i_dim, j_dim, k_dim,                 &
            long_name  ='3D pH at current time',                                       &
            units      ='pH', grid_loc='3111',                                         &
            field_loc  = field_loc_center,                                             &
            field_type = field_type_scalar,                                            &
            d3d_array  = ecosys_saved_state%ph_prev_3d(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', ph_3d_iodesc)

    else if (trim(action) == 'write') then

       call data_set (restart_file, 'write', ph_surf_iodesc)
       call data_set (restart_file, 'write', ph_surf_alt_co2_iodesc)
       call data_set (restart_file, 'write', ph_3d_iodesc)
       call data_set (restart_file, 'write', ph_3d_alt_co2_iodesc)

    endif

  end subroutine ecosys_driver_write_ecosys_restart

  !*****************************************************************************

  subroutine ecosys_driver_ciso_write_restart(restart_file, action)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files
    !-----------------------------------------------------------------------

    use constants , only : char_blank
    use constants , only : field_loc_center
    use constants , only : field_type_scalar
    use io        , only : datafile
    use io_types  , only : add_attrib_file

    implicit none

    character(*)   , intent(in)     :: action
    type (datafile), intent (inout) :: restart_file

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (char_len) :: short_name   ! tracer name temporaries
    integer (int_kind)   :: n
    !-----------------------------------------------------------------------
    
    if (trim(action) == 'add_attrib_file') then
       short_name = char_blank
       do n = 1,ecosys_ciso_tracer_cnt
          if (ciso_vflux_flag(n)) then
             short_name = 'surf_avg_' // ciso_ind_name_table(n)%name
             call add_attrib_file(restart_file, trim(short_name), ciso_surf_avg(n))
          endif
       end do
    endif
    
  end subroutine ecosys_driver_ciso_write_restart

  !*****************************************************************************

  subroutine ecosys_driver_init_interior_restore(ecosys_restore)

    ! !DESCRIPTION:
    !  Initialize interior restoring computations for ecosys tracer module.

    use ecosys_restore_mod    , only : ecosys_restore_type
    use grid                  , only : KMT
    use grid                  , only : zt
    use io_types              , only : nml_in, nml_filename
    use blocks                , only : nx_block, ny_block
    use domain_size           , only : max_blocks_clinic, km
    use marbl_share_mod       , only : fesedflux_input
    use passive_tracer_tools  , only : read_field

    implicit none

    type(ecosys_restore_type)   , intent(inout) :: ecosys_restore

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: &
         k,               & ! index for looping over levels
         i, j,            & ! index for looping over horiz. dims.
         iblock             ! index for looping over blocks

    real (r8) :: &
         subsurf_fesed      ! sum of subsurface fesed values
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  initialize restoring timescale (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%initialize_restoring_timescale(nml_filename, nml_in, zt)

    !-----------------------------------------------------------------------
    !  load restoring fields (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%read_restoring_fields(land_mask)

    !-----------------------------------------------------------------------
    !  load fesedflux
    !  add subsurface positives to 1 level shallower, to accomodate overflow pop-ups
    !-----------------------------------------------------------------------

    allocate(fesedflux(nx_block, ny_block, km, max_blocks_clinic))

    call read_field(fesedflux_input%file_fmt, &
         fesedflux_input%filename, &
         fesedflux_input%file_varname, &
         fesedflux)

    do iblock=1, nblocks_clinic
       do j=1, ny_block
          do i=1, nx_block
             if (KMT(i, j, iblock) > 0 .and. KMT(i, j, iblock) < km) then
                subsurf_fesed = c0
                do k=KMT(i, j, iblock)+1, km
                   subsurf_fesed = subsurf_fesed + fesedflux(i, j, k, iblock)
                enddo
                fesedflux(i, j, KMT(i, j, iblock), iblock) = fesedflux(i, j, KMT(i, j, iblock), iblock) + subsurf_fesed
             endif
          enddo
       enddo

       do k = 1, km
          where (.not. land_mask(:, :, iblock) .or. k > KMT(:, :, iblock)) &
               fesedflux(:, :, k, iblock) = c0
          fesedflux(:, :, k, iblock) = fesedflux(:, :, k, iblock) * fesedflux_input%scale_factor
       enddo
    end do

  end subroutine ecosys_driver_init_interior_restore

  !***********************************************************************

  subroutine ecosys_driver_init_sflux()

    ! !DESCRIPTION:
    !  Initialize surface flux computations for ecosys tracer module.

    use registry              , only : registry_match
    use passive_tracer_tools  , only : read_field
    use forcing_tools         , only : find_forcing_times
    use io_tools              , only : document
    use named_field_mod       , only : named_field_get_index
    use named_field_mod       , only : named_field_get
    use named_field_mod       , only : named_field_register
    use named_field_mod       , only : named_field_set
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : atm_co2_iopt
    use marbl_share_mod       , only : atm_co2_iopt_drv_prog
    use marbl_share_mod       , only : atm_co2_iopt_drv_diag
    use marbl_share_mod       , only : gas_flux_forcing_iopt
    use marbl_share_mod       , only : gas_flux_forcing_iopt_file
    use marbl_share_mod       , only : gas_flux_forcing_iopt_drv
    use marbl_share_mod       , only : gas_flux_forcing_file
    use marbl_share_mod       , only : fice_file
    use marbl_share_mod       , only : xkw_file
    use marbl_share_mod       , only : ap_file
    use marbl_share_mod       , only : dust_flux
    use marbl_share_mod       , only : iron_flux
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : nox_flux_monthly
    use marbl_share_mod       , only : nhy_flux_monthly
    use marbl_share_mod       , only : din_riv_flux
    use marbl_share_mod       , only : dip_riv_flux
    use marbl_share_mod       , only : don_riv_flux
    use marbl_share_mod       , only : dop_riv_flux
    use marbl_share_mod       , only : dsi_riv_flux
    use marbl_share_mod       , only : dfe_riv_flux
    use marbl_share_mod       , only : dic_riv_flux
    use marbl_share_mod       , only : doc_riv_flux
    use marbl_share_mod       , only : alk_riv_flux
    use marbl_share_mod       , only : liron_patch  
    use marbl_share_mod       , only : iron_patch_flux_filename  
    use marbl_share_mod       , only : iron_patch_month  

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_init_sflux'

    logical (log_kind) :: &
         luse_INTERP_WORK     ! does INTERP_WORK need to be allocated

    integer (int_kind) :: &
         n,                 & ! index for looping over tracers
         iblock               ! index for looping over blocks

    real (r8), dimension (nx_block, ny_block) :: &
         WORK

    real (r8), dimension (nx_block, ny_block, 12, max_blocks_clinic), target :: &
         WORK_READ            ! temporary space to read in fields

    real (r8) :: INTERP_WORK(nx_block, ny_block, max_blocks_clinic, 1) ! temp array for interpolate_forcing output

    !-----------------------------------------------------------------------

    luse_INTERP_WORK = .false.

    !-----------------------------------------------------------------------
    !  read gas flux forcing (if required)
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

       luse_INTERP_WORK = .true.

       !-----------------------------------------------------------------------
       !  first, read ice file
       !-----------------------------------------------------------------------

       allocate(fice_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(fice_file%input%filename) == 'unknown') &
            fice_file%input%filename = gas_flux_forcing_file

       call read_field(fice_file%input%file_fmt, &
            fice_file%input%filename, &
            fice_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             fice_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  fice_file%DATA(:, :, iblock, 1, n) = c0
             fice_file%DATA(:, :, iblock, 1, n) = &
                  fice_file%DATA(:, :, iblock, 1, n) * fice_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(fice_file%data_time, &
            fice_file%data_inc, fice_file%interp_type, &
            fice_file%data_next, fice_file%data_time_min_loc, &
            fice_file%data_update, fice_file%data_type)

       !-----------------------------------------------------------------------
       !  next, read piston velocity file
       !-----------------------------------------------------------------------

       allocate(xkw_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(xkw_file%input%filename) == 'unknown') &
            xkw_file%input%filename = gas_flux_forcing_file

       call read_field(xkw_file%input%file_fmt, &
            xkw_file%input%filename, &
            xkw_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             xkw_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  xkw_file%DATA(:, :, iblock, 1, n) = c0
             xkw_file%DATA(:, :, iblock, 1, n) = &
                  xkw_file%DATA(:, :, iblock, 1, n) * xkw_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(xkw_file%data_time, &
            xkw_file%data_inc, xkw_file%interp_type, &
            xkw_file%data_next, xkw_file%data_time_min_loc, &
            xkw_file%data_update, xkw_file%data_type)

       !-----------------------------------------------------------------------
       !  last, read atmospheric pressure file
       !-----------------------------------------------------------------------

       allocate(ap_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(ap_file%input%filename) == 'unknown') &
            ap_file%input%filename = gas_flux_forcing_file

       call read_field(ap_file%input%file_fmt, &
            ap_file%input%filename, &
            ap_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             ap_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  ap_file%DATA(:, :, iblock, 1, n) = c0
             ap_file%DATA(:, :, iblock, 1, n) = &
                  ap_file%DATA(:, :, iblock, 1, n) * ap_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(ap_file%data_time, &
            ap_file%data_inc, ap_file%interp_type, &
            ap_file%data_next, ap_file%data_time_min_loc, &
            ap_file%data_update, ap_file%data_type)

    endif

    !-----------------------------------------------------------------------
    !  load dust flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(dust_flux%input%filename) /= 'none' .and. &
         trim(dust_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dust_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(dust_flux%input%filename) == 'unknown') &
            dust_flux%input%filename = gas_flux_forcing_file

       call read_field(dust_flux%input%file_fmt, &
            dust_flux%input%filename, &
            dust_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dust_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dust_flux%DATA(:, :, iblock, 1, n) = c0
             dust_flux%DATA(:, :, iblock, 1, n) = &
                  dust_flux%DATA(:, :, iblock, 1, n) * dust_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dust_flux%data_time, &
            dust_flux%data_inc, dust_flux%interp_type, &
            dust_flux%data_next, dust_flux%data_time_min_loc, &
            dust_flux%data_update, dust_flux%data_type)

       dust_flux%has_data = .true.
    else
       dust_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  load iron flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(iron_flux%input%filename) /= 'none' .and. &
         trim(iron_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(iron_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(iron_flux%input%filename) == 'unknown') &
            iron_flux%input%filename = gas_flux_forcing_file

       call read_field(iron_flux%input%file_fmt, &
            iron_flux%input%filename, &
            iron_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             iron_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  iron_flux%DATA(:, :, iblock, 1, n) = c0
             iron_flux%DATA(:, :, iblock, 1, n) = &
                  iron_flux%DATA(:, :, iblock, 1, n) * iron_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(iron_flux%data_time, &
            iron_flux%data_inc, iron_flux%interp_type, &
            iron_flux%data_next, iron_flux%data_time_min_loc, &
            iron_flux%data_update, iron_flux%data_type)

       iron_flux%has_data = .true.
    else
       iron_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  load nox & noy flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(ndep_data_type) /= 'none' .and. &
         trim(ndep_data_type) /= 'monthly-calendar' .and. &
         trim(ndep_data_type) /= 'shr_stream') then
       call document(subname, 'ndep_data_type', ndep_data_type)
       call exit_POP(sigAbort, 'unknown ndep_data_type')
    endif

    if (trim(ndep_data_type) == 'monthly-calendar' .and. &
         trim(nox_flux_monthly%input%filename) /= 'none' .and. &
         trim(nox_flux_monthly%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(nox_flux_monthly%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(nox_flux_monthly%input%filename) == 'unknown') &
            nox_flux_monthly%input%filename = gas_flux_forcing_file

       call read_field(nox_flux_monthly%input%file_fmt, &
            nox_flux_monthly%input%filename, &
            nox_flux_monthly%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             nox_flux_monthly%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  nox_flux_monthly%DATA(:, :, iblock, 1, n) = c0
             nox_flux_monthly%DATA(:, :, iblock, 1, n) = &
                  nox_flux_monthly%DATA(:, :, iblock, 1, n) * nox_flux_monthly%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(nox_flux_monthly%data_time, &
            nox_flux_monthly%data_inc, nox_flux_monthly%interp_type, &
            nox_flux_monthly%data_next, nox_flux_monthly%data_time_min_loc, &
            nox_flux_monthly%data_update, nox_flux_monthly%data_type)

       nox_flux_monthly%has_data = .true.
    else
       nox_flux_monthly%has_data = .false.
    endif

    if (trim(ndep_data_type) == 'monthly-calendar' .and. &
         trim(nhy_flux_monthly%input%filename) /= 'none' .and. &
         trim(nhy_flux_monthly%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(nhy_flux_monthly%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(nhy_flux_monthly%input%filename) == 'unknown') &
            nhy_flux_monthly%input%filename = gas_flux_forcing_file

       call read_field(nhy_flux_monthly%input%file_fmt, &
            nhy_flux_monthly%input%filename, &
            nhy_flux_monthly%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             nhy_flux_monthly%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  nhy_flux_monthly%DATA(:, :, iblock, 1, n) = c0
             nhy_flux_monthly%DATA(:, :, iblock, 1, n) = &
                  nhy_flux_monthly%DATA(:, :, iblock, 1, n) * nhy_flux_monthly%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(nhy_flux_monthly%data_time, &
            nhy_flux_monthly%data_inc, nhy_flux_monthly%interp_type, &
            nhy_flux_monthly%data_next, nhy_flux_monthly%data_time_min_loc, &
            nhy_flux_monthly%data_update, nhy_flux_monthly%data_type)

       nhy_flux_monthly%has_data = .true.
    else
       nhy_flux_monthly%has_data = .false.
    endif

    !-----------------------------------------------------------------------

    if (trim(ndep_data_type) == 'shr_stream') then
       call document(subname, 'ndep_data_type', ndep_data_type)
       call exit_POP(sigAbort, &
            'shr_stream option only supported when CCSMCOUPLED is defined')
    endif

    !-----------------------------------------------------------------------
    !  load river nutrient flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(din_riv_flux%input%filename) /= 'none' .and. &
        trim(din_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(din_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(din_riv_flux%input%file_fmt, &
            din_riv_flux%input%filename, &
            din_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             din_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  din_riv_flux%DATA(:, :, iblock, 1, n) = c0
             din_riv_flux%DATA(:, :, iblock, 1, n) = &
                  din_riv_flux%DATA(:, :, iblock, 1, n) * din_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(din_riv_flux%data_time, &
            din_riv_flux%data_inc, din_riv_flux%interp_type, &
            din_riv_flux%data_next, din_riv_flux%data_time_min_loc, &
            din_riv_flux%data_update, din_riv_flux%data_type)

       din_riv_flux%has_data = .true.
    else
       din_riv_flux%has_data = .false.
    endif


    if (trim(dip_riv_flux%input%filename) /= 'none' .and. &
        trim(dip_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dip_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dip_riv_flux%input%file_fmt, &
            dip_riv_flux%input%filename, &
            dip_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dip_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dip_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dip_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dip_riv_flux%DATA(:, :, iblock, 1, n) * dip_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dip_riv_flux%data_time, &
            dip_riv_flux%data_inc, dip_riv_flux%interp_type, &
            dip_riv_flux%data_next, dip_riv_flux%data_time_min_loc, &
            dip_riv_flux%data_update, dip_riv_flux%data_type)

       dip_riv_flux%has_data = .true.
    else
       dip_riv_flux%has_data = .false.
    endif


    if (trim(don_riv_flux%input%filename) /= 'none' .and. &
        trim(don_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(don_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(don_riv_flux%input%file_fmt, &
            don_riv_flux%input%filename, &
            don_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             don_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  don_riv_flux%DATA(:, :, iblock, 1, n) = c0
             don_riv_flux%DATA(:, :, iblock, 1, n) = &
                  don_riv_flux%DATA(:, :, iblock, 1, n) * don_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(don_riv_flux%data_time, &
            don_riv_flux%data_inc, don_riv_flux%interp_type, &
            don_riv_flux%data_next, don_riv_flux%data_time_min_loc, &
            don_riv_flux%data_update, don_riv_flux%data_type)

       don_riv_flux%has_data = .true.
    else
       don_riv_flux%has_data = .false.
    endif


    if (trim(dop_riv_flux%input%filename) /= 'none' .and. &
        trim(dop_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dop_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dop_riv_flux%input%file_fmt, &
            dop_riv_flux%input%filename, &
            dop_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dop_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dop_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dop_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dop_riv_flux%DATA(:, :, iblock, 1, n) * dop_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dop_riv_flux%data_time, &
            dop_riv_flux%data_inc, dop_riv_flux%interp_type, &
            dop_riv_flux%data_next, dop_riv_flux%data_time_min_loc, &
            dop_riv_flux%data_update, dop_riv_flux%data_type)

       dop_riv_flux%has_data = .true.
    else
       dop_riv_flux%has_data = .false.
    endif

    if (trim(dsi_riv_flux%input%filename) /= 'none' .and. &
        trim(dsi_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dsi_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dsi_riv_flux%input%file_fmt, &
            dsi_riv_flux%input%filename, &
            dsi_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dsi_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dsi_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dsi_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dsi_riv_flux%DATA(:, :, iblock, 1, n) * dsi_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dsi_riv_flux%data_time, &
            dsi_riv_flux%data_inc, dsi_riv_flux%interp_type, &
            dsi_riv_flux%data_next, dsi_riv_flux%data_time_min_loc, &
            dsi_riv_flux%data_update, dsi_riv_flux%data_type)

       dsi_riv_flux%has_data = .true.
    else
       dsi_riv_flux%has_data = .false.
    endif


    if (trim(dfe_riv_flux%input%filename) /= 'none' .and. &
        trim(dfe_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dfe_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dfe_riv_flux%input%file_fmt, &
            dfe_riv_flux%input%filename, &
            dfe_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dfe_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dfe_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dfe_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dfe_riv_flux%DATA(:, :, iblock, 1, n) * dfe_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dfe_riv_flux%data_time, &
            dfe_riv_flux%data_inc, dfe_riv_flux%interp_type, &
            dfe_riv_flux%data_next, dfe_riv_flux%data_time_min_loc, &
            dfe_riv_flux%data_update, dfe_riv_flux%data_type)

       dfe_riv_flux%has_data = .true.
    else
       dfe_riv_flux%has_data = .false.
    endif


    if (trim(dic_riv_flux%input%filename) /= 'none' .and. &
        trim(dic_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dic_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dic_riv_flux%input%file_fmt, &
            dic_riv_flux%input%filename, &
            dic_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dic_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  dic_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dic_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dic_riv_flux%DATA(:, :, iblock, 1, n) * dic_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dic_riv_flux%data_time, &
            dic_riv_flux%data_inc, dic_riv_flux%interp_type, &
            dic_riv_flux%data_next, dic_riv_flux%data_time_min_loc, &
            dic_riv_flux%data_update, dic_riv_flux%data_type)

       dic_riv_flux%has_data = .true.
    else
       dic_riv_flux%has_data = .false.
    endif


    if (trim(alk_riv_flux%input%filename) /= 'none' .and. &
        trim(alk_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(alk_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(alk_riv_flux%input%file_fmt, &
            alk_riv_flux%input%filename, &
            alk_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             alk_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  alk_riv_flux%DATA(:, :, iblock, 1, n) = c0
             alk_riv_flux%DATA(:, :, iblock, 1, n) = &
                  alk_riv_flux%DATA(:, :, iblock, 1, n) * alk_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(alk_riv_flux%data_time, &
            alk_riv_flux%data_inc, alk_riv_flux%interp_type, &
            alk_riv_flux%data_next, alk_riv_flux%data_time_min_loc, &
            alk_riv_flux%data_update, alk_riv_flux%data_type)

       alk_riv_flux%has_data = .true.
    else
       alk_riv_flux%has_data = .false.
    endif

    if (trim(doc_riv_flux%input%filename) /= 'none' .and. &
        trim(doc_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(doc_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(doc_riv_flux%input%file_fmt, &
            doc_riv_flux%input%filename, &
            doc_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             doc_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. land_mask(:, :, iblock)) &
                  doc_riv_flux%DATA(:, :, iblock, 1, n) = c0
             doc_riv_flux%DATA(:, :, iblock, 1, n) = &
                  doc_riv_flux%DATA(:, :, iblock, 1, n) * doc_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(doc_riv_flux%data_time, &
            doc_riv_flux%data_inc, doc_riv_flux%interp_type, &
            doc_riv_flux%data_next, doc_riv_flux%data_time_min_loc, &
            doc_riv_flux%data_update, doc_riv_flux%data_type)

       doc_riv_flux%has_data = .true.
    else
       doc_riv_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  load iron PATCH flux fields (if required)
    !-----------------------------------------------------------------------

    if (liron_patch) then

       !maltrud iron patch
       !  assume patch file has same normalization and format as deposition file

       allocate(IRON_PATCH_FLUX(nx_block, ny_block, max_blocks_clinic))

       if (trim(iron_flux%input%filename) == 'unknown') &
            iron_flux%input%filename = gas_flux_forcing_file

       call read_field(iron_flux%input%file_fmt, &
            iron_flux%input%filename, &
            iron_patch_flux_filename, &
            IRON_PATCH_FLUX)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             where (.not. land_mask(:, :, iblock)) &
                  IRON_PATCH_FLUX(:, :, iblock) = c0
             iron_flux%DATA(:, :, iblock, 1, n) = &
                  IRON_PATCH_FLUX(:, :, iblock) * iron_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    !-----------------------------------------------------------------------
    !  register and set
    !     fco2, the air-sea co2 gas flux
    !-----------------------------------------------------------------------

    call named_field_register('SFLUX_CO2', sflux_co2_nf_ind)
    !$OMP PARALLEL DO PRIVATE(iblock, WORK)
    do iblock=1, nblocks_clinic
       WORK = c0
       call named_field_set(sflux_co2_nf_ind, iblock, WORK)
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !  verify running coupled if gas fluxes use coupler forcing
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
         (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv .or. &
         atm_co2_iopt == atm_co2_iopt_drv_prog .or. &
         atm_co2_iopt == atm_co2_iopt_drv_diag) .and. &
         .not. registry_match('lcoupled')) then
       call exit_POP(sigAbort, 'ecosys_init: ecosys module requires the ' /&
            &/ 'flux coupler when gas_flux_forcing_opt=drv')
    endif

    !-----------------------------------------------------------------------
    !  get named field index for atmospheric CO2, if required
    !-----------------------------------------------------------------------

    if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_prog) then
       call named_field_get_index('ATM_CO2_PROG', atm_co2_nf_ind, &
            exit_on_err=.false.)
       if (atm_co2_nf_ind == 0) then
          call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
               &/ 'atmopsheric CO2 from the flux coupler ' /&
               &/ 'and it is not present')
       endif
    endif

    if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_diag) then
       call named_field_get_index('ATM_CO2_DIAG', atm_co2_nf_ind, &
            exit_on_err=.false.)
       if (atm_co2_nf_ind == 0) then
          call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
               &/ 'atmopsheric CO2 from the flux coupler ' /&
               &/ 'and it is not present')
       endif
    endif

  end subroutine ecosys_driver_init_sflux

  !*****************************************************************************

  subroutine ecosys_driver_read_sflux( &
       input_forcing_data,             &
       XCO2, XCO2_ALT_CO2,             &
       IFRAC_FILE_INPUT, XKW_FILE_INPUT, AP_FILE_INPUT, IRON_FLUX_IN)

    use shr_pio_mod           , only : shr_pio_getiotype, shr_pio_getiosys
    use POP_IOUnitsMod        , only : inst_name
    use POP_HaloMod           , only : POP_HaloUpdate 
    use POP_GridHorzMod       , only : POP_gridHorzLocCenter 
    use POP_CommMod           , only : POP_communicator 
    use POP_FieldMod          , only : POP_fieldKindScalar
    use domain                , only : blocks_clinic
    use domain                , only : POP_haloClinic
    use blocks                , only : get_block
    use constants             , only : field_loc_center
    use constants             , only : field_type_scalar
    use constants             , only : hflux_factor
    use passive_tracer_tools  , only : forcing_monthly_every_ts
    use forcing_tools         , only : interpolate_forcing
    use forcing_tools         , only : update_forcing_data
    use named_field_mod       , only : named_field_get
    use time_management       , only : isecond
    use time_management       , only : iminute
    use time_management       , only : ihour
    use time_management       , only : iday
    use time_management       , only : imonth
    use time_management       , only : iyear
    use time_management       , only : thour00
    use shr_strdata_mod       , only : shr_strdata_type
    use shr_strdata_mod       , only : shr_strdata_advance 
    use strdata_interface_mod , only : strdata_input_type
    use strdata_interface_mod , only : POP_strdata_create
    use strdata_interface_mod , only : POP_strdata_advance

    use marbl_oxygen          , only : o2sat
    use marbl_parms           , only : f_qsw_par
    use marbl_parms           , only : ind_nox_flux
    use marbl_parms           , only : ind_nhy_flux
    use marbl_parms           , only : ind_nh4_flux     
    use marbl_parms           , only : ind_no3_flux     
    use marbl_parms           , only : ind_din_riv_flux 
    use marbl_parms           , only : ind_dip_riv_flux 
    use marbl_parms           , only : ind_don_riv_flux 
    use marbl_parms           , only : ind_dop_riv_flux 
    use marbl_parms           , only : ind_dsi_riv_flux 
    use marbl_parms           , only : ind_dfe_riv_flux
    use marbl_parms           , only : ind_dic_riv_flux
    use marbl_parms           , only : ind_alk_riv_flux
    use marbl_parms           , only : ind_doc_riv_flux

    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : gas_flux_forcing_iopt
    use marbl_share_mod       , only : gas_flux_forcing_iopt_file
    use marbl_share_mod       , only : atm_co2_iopt
    use marbl_share_mod       , only : atm_co2_iopt_const
    use marbl_share_mod       , only : atm_co2_iopt_drv_prog
    use marbl_share_mod       , only : atm_co2_iopt_drv_diag
    use marbl_share_mod       , only : atm_co2_const
    use marbl_share_mod       , only : atm_alt_co2_const
    use marbl_share_mod       , only : atm_alt_co2_iopt
    use marbl_share_mod       , only : ndep_data_type 
    use marbl_share_mod       , only : ndep_shr_stream_var_cnt
    use marbl_share_mod       , only : ndep_shr_stream_no_ind
    use marbl_share_mod       , only : ndep_shr_stream_nh_ind
    use marbl_share_mod       , only : ndep_shr_stream_file
    use marbl_share_mod       , only : ndep_shr_stream_year_first
    use marbl_share_mod       , only : ndep_shr_stream_year_last
    use marbl_share_mod       , only : ndep_shr_stream_year_align
    use marbl_share_mod       , only : dust_flux        
    use marbl_share_mod       , only : iron_flux        
    use marbl_share_mod       , only : fice_file        
    use marbl_share_mod       , only : xkw_file         
    use marbl_share_mod       , only : ap_file          
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
    use marbl_share_mod       , only : din_riv_flux     
    use marbl_share_mod       , only : dip_riv_flux     
    use marbl_share_mod       , only : don_riv_flux     
    use marbl_share_mod       , only : dop_riv_flux     
    use marbl_share_mod       , only : dsi_riv_flux     
    use marbl_share_mod       , only : dfe_riv_flux     
    use marbl_share_mod       , only : dic_riv_flux     
    use marbl_share_mod       , only : alk_riv_flux     
    use marbl_share_mod       , only : doc_riv_flux     
    use marbl_share_mod       , only : liron_patch  
    use marbl_share_mod       , only : iron_patch_month  
    use passive_tracer_tools  , only : comp_surf_avg

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.

    real (r8), dimension(:, :, :)    , intent(out) :: XCO2             ! atmospheric co2 conc. (dry-air, 1 atm)
    real (r8), dimension(:, :, :)    , intent(out) :: XCO2_ALT_CO2     ! atmospheric alternative CO2 (dry-air, 1 atm)
    real (r8), dimension(:, :, :)    , intent(out) :: IFRAC_FILE_INPUT ! used ice fraction (non-dimensional)
    real (r8), dimension(:, :, :)    , intent(out) :: XKW_FILE_INPUT   ! portion of piston velocity (cm/s)
    real (r8), dimension(:, :, :)    , intent(out) :: AP_FILE_INPUT    ! used atm pressure (converted from dyne/cm**2 to atm)
    real (r8), dimension(:, :, :)    , intent(out) :: IRON_FLUX_IN     ! iron flux
    real (r8), dimension(:, :, :, :) , intent(out) :: INPUT_FORCING_DATA

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_driver:ecosys_driver_read_sflux'

    logical (log_kind) :: first_call = .true.

    type (block) :: this_block      ! block info for the current block

    integer (int_kind) :: i, j, iblock, n ! loop indices
    integer (int_kind) :: auto_ind        ! autotroph functional group index
    integer (int_kind) :: mcdate, sec     ! date vals for shr_strdata_advance
    integer (int_kind) :: errorCode       ! errorCode from HaloUpdate call

    real (r8) :: scalar_temp
    real (r8) :: INTERP_WORK(nx_block, ny_block, max_blocks_clinic, 1) ! temp array for interpolate_forcing output
    real (r8) :: SHR_STREAM_WORK(nx_block, ny_block, max_blocks_clinic)

    character (char_len)     :: tracer_data_label       ! label for what is being updated
    character (char_len)     :: tracer_data_names(1)    ! short names for input data fields
    integer (int_kind)       :: tracer_bndy_loc(1)      ! location   for ghost tracer_bndy_type cell updates
    integer (int_kind)       :: tracer_bndy_type(1)     ! field type for ghost tracer_bndy_type cell updates
    character (char_len)     :: ndep_shr_stream_fldList ! label for what is being updated
    type(shr_strdata_type)   :: ndep_sdat               ! input data stream for ndep
    type(strdata_input_type) :: ndep_inputlist
    !-----------------------------------------------------------------------

    call timer_start(ecosys_pre_sflux_timer)

    !-----------------------------------------------------------------------
    !  fluxes initially set to 0
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       input_forcing_data(:, :, :, iblock) = c0
    enddo

    !-----------------------------------------------------------------------
    !  Interpolate gas flux forcing data if necessary
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
         gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(fice_file%data_time,      &
               fice_file%data_time_min_loc,  fice_file%interp_type,    &
               fice_file%data_next,          fice_file%data_update,    &
               fice_file%data_type,          fice_file%data_inc,       &
               fice_file%DATA(:, :, :, :, 1:12), fice_file%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:, :, :, :, 1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_FILE_INPUT = INTERP_WORK(:, :, :, 1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(xkw_file%data_time,      &
               xkw_file%data_time_min_loc,  xkw_file%interp_type,    &
               xkw_file%data_next,          xkw_file%data_update,    &
               xkw_file%data_type,          xkw_file%data_inc,       &
               xkw_file%DATA(:, :, :, :, 1:12), xkw_file%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            xkw_file%DATA(:, :, :, :, 1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_FILE_INPUT = INTERP_WORK(:, :, :, 1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(ap_file%data_time,    &
               ap_file%data_time_min_loc,  ap_file%interp_type,    &
               ap_file%data_next,          ap_file%data_update,    &
               ap_file%data_type,          ap_file%data_inc,       &
               ap_file%DATA(:, :, :, :, 1:12), ap_file%data_renorm,    &
               tracer_data_label,          tracer_data_names,      &
               tracer_bndy_loc,            tracer_bndy_type,       &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:, :, :, :, 1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_FILE_INPUT = INTERP_WORK(:, :, :, 1)

    endif

    !-----------------------------------------------------------------------
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  compute CO2 flux, computing disequilibrium one row at a time
    !-----------------------------------------------------------------------

    if (lflux_gas_co2) then

       !-----------------------------------------------------------------------
       !  Set XCO2
       !-----------------------------------------------------------------------

       select case (atm_co2_iopt)
       case (atm_co2_iopt_const)
          XCO2 = atm_co2_const
       case (atm_co2_iopt_drv_prog, atm_co2_iopt_drv_diag)
          do iblock = 1, nblocks_clinic
             call named_field_get(atm_co2_nf_ind, iblock, XCO2(:, :, iblock))
          end do
       end select

       select case (atm_alt_co2_iopt)
       case (atm_co2_iopt_const)
          XCO2_ALT_CO2 = atm_alt_co2_const
       end select

    endif  !  lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    if (iron_flux%has_data) then
       if (thour00 >= iron_flux%data_update) then
          tracer_data_names = iron_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Iron Flux'
          call update_forcing_data(iron_flux%data_time,    &
               iron_flux%data_time_min_loc,  iron_flux%interp_type,    &
               iron_flux%data_next,          iron_flux%data_update,    &
               iron_flux%data_type,          iron_flux%data_inc,       &
               iron_flux%DATA(:, :, :, :, 1:12), iron_flux%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               iron_flux%filename,           iron_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            iron_flux%DATA(:, :, :, :, 1:12), &
            iron_flux%data_time,         iron_flux%interp_type, &
            iron_flux%data_time_min_loc, iron_flux%interp_freq, &
            iron_flux%interp_inc,        iron_flux%interp_next, &
            iron_flux%interp_last,       0)
       if (liron_patch .and. imonth == iron_patch_month) then
          IRON_FLUX_IN = INTERP_WORK(:, :, :, 1) + IRON_PATCH_FLUX
       else
          IRON_FLUX_IN = INTERP_WORK(:, :, :, 1)
       endif
    else
       IRON_FLUX_IN = c0
    endif

    if (dust_flux%has_data) then
       if (thour00 >= dust_flux%data_update) then
          tracer_data_names = dust_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Dust Flux'
          call update_forcing_data(dust_flux%data_time,    &
               dust_flux%data_time_min_loc,  dust_flux%interp_type,    &
               dust_flux%data_next,          dust_flux%data_update,    &
               dust_flux%data_type,          dust_flux%data_inc,       &
               dust_flux%DATA(:, :, :, :, 1:12), dust_flux%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               dust_flux%filename,           dust_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            dust_flux%DATA(:, :, :, :, 1:12),    &
            dust_flux%data_time,         dust_flux%interp_type, &
            dust_flux%data_time_min_loc, dust_flux%interp_freq, &
            dust_flux%interp_inc,        dust_flux%interp_next, &
            dust_flux%interp_last,       0)

       ecosys_saved_state%dust_FLUX_IN = INTERP_WORK(:, :, :, 1)

       !-----------------------------------------------------------------------
       !  Reduce surface dust flux due to assumed instant surface dissolution
       !  Can't use parm_fe_bioavail when using solFe input files
       !-----------------------------------------------------------------------

       !     dust_FLUX_IN = dust_FLUX_IN * (c1 - parm_Fe_bioavail)
       ecosys_saved_state%dust_FLUX_IN = ecosys_saved_state%dust_FLUX_IN * 0.98_r8
    else
       ecosys_saved_state%dust_FLUX_IN = c0
    endif

    !JW TODO: dust_FLUX_IN in through structure?

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------

    if (nox_flux_monthly%has_data) then
       if (thour00 >= nox_flux_monthly%data_update) then
          tracer_data_names = nox_flux_monthly%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'NOx Flux'
          call update_forcing_data(nox_flux_monthly%data_time,    &
               nox_flux_monthly%data_time_min_loc, nox_flux_monthly%interp_type, &
               nox_flux_monthly%data_next, nox_flux_monthly%data_update, &
               nox_flux_monthly%data_type,          nox_flux_monthly%data_inc, &
               nox_flux_monthly%DATA(:, :, :, :, 1:12), nox_flux_monthly%data_renorm, &
               tracer_data_label,                   tracer_data_names, &
               tracer_bndy_loc,                     tracer_bndy_type, &
               nox_flux_monthly%filename, nox_flux_monthly%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            nox_flux_monthly%DATA(:, :, :, :, 1:12), &
            nox_flux_monthly%data_time,         nox_flux_monthly%interp_type, &
            nox_flux_monthly%data_time_min_loc, nox_flux_monthly%interp_freq, &
            nox_flux_monthly%interp_inc,        nox_flux_monthly%interp_next, &
            nox_flux_monthly%interp_last,       0)
        input_forcing_data(:,:, ind_nox_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (nhy_flux_monthly%has_data) then
       if (thour00 >= nhy_flux_monthly%data_update) then
          tracer_data_names = nhy_flux_monthly%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'NHy Flux'
          call update_forcing_data(nhy_flux_monthly%data_time,    &
               nhy_flux_monthly%data_time_min_loc, nhy_flux_monthly%interp_type, &
               nhy_flux_monthly%data_next, nhy_flux_monthly%data_update, &
               nhy_flux_monthly%data_type,          nhy_flux_monthly%data_inc, &
               nhy_flux_monthly%DATA(:, :, :, :, 1:12), nhy_flux_monthly%data_renorm, &
               tracer_data_label,                   tracer_data_names, &
               tracer_bndy_loc,                     tracer_bndy_type, &
               nhy_flux_monthly%filename, nhy_flux_monthly%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            nhy_flux_monthly%DATA(:, :, :, :, 1:12), &
            nhy_flux_monthly%data_time,         nhy_flux_monthly%interp_type, &
            nhy_flux_monthly%data_time_min_loc, nhy_flux_monthly%interp_freq, &
            nhy_flux_monthly%interp_inc,        nhy_flux_monthly%interp_next, &
            nhy_flux_monthly%interp_last,       0)
        input_forcing_data(:,:, ind_nhy_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (trim(ndep_data_type) == 'shr_stream') then
       if (first_call) then

          ndep_inputlist%field_name = 'ndep data'
          ndep_inputlist%short_name = 'ndep'
          ndep_inputlist%year_first = ndep_shr_stream_year_first
          ndep_inputlist%year_last  = ndep_shr_stream_year_last
          ndep_inputlist%year_align = ndep_shr_stream_year_align
          ndep_inputlist%file_name  = ndep_shr_stream_file
          ndep_inputlist%field_list = ' '
          do n = 1, ndep_shr_stream_var_cnt
             if (n == ndep_shr_stream_no_ind) then
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) // 'NOy_deposition'
               end if
             if (n == ndep_shr_stream_nh_ind) then
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) // 'NHx_deposition'
               end if
             if (n < ndep_shr_stream_var_cnt) then
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) // ':'
               end if
          end do

          call POP_strdata_create(ndep_inputlist)

          first_call = .false.
       endif

       mcdate = iyear*10000 + imonth*100 + iday
       sec = isecond + 60 * (iminute + 60 * ihour)

       call timer_start(ecosys_shr_strdata_advance_timer)
       call shr_strdata_advance(ndep_sdat, mcdate, sec, POP_communicator, 'ndep')
       call timer_stop(ecosys_shr_strdata_advance_timer)

       ! process NO3 flux, store results in SHR_STREAM_WORK array
       ! instead of directly into STF_MODULE
       ! to avoid argument copies in HaloUpdate calls

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock), iblock)
          do j=this_block%jb, this_block%je
             do i=this_block%ib, this_block%ie
                n = n + 1
                SHR_STREAM_WORK(i, j, iblock) = ndep_sdat%avs(1)%rAttr(ndep_shr_stream_no_ind, n)
             enddo
          enddo
       enddo

       call POP_HaloUpdate(SHR_STREAM_WORK, POP_haloClinic, &
            POP_gridHorzLocCenter,          &
            POP_fieldKindScalar, errorCode, &
            fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call exit_POP(sigAbort, subname /&
               &/ ': error updating halo for Ndep fields')
       endif

       do iblock = 1, nblocks_clinic
          where (land_mask(:, :, iblock))
             input_forcing_data(:,:, ind_no3_flux,iblock) = SHR_STREAM_WORK(:, :, iblock)
          endwhere
       enddo

       ! process NH4 flux, store results in SHR_STREAM_WORK array
       ! instead of directly into STF_MODULE
       ! to avoid argument copies in HaloUpdate calls

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock), iblock)
          do j=this_block%jb, this_block%je
             do i=this_block%ib, this_block%ie
                n = n + 1
                SHR_STREAM_WORK(i, j, iblock) = &
                     ndep_sdat%avs(1)%rAttr(ndep_shr_stream_nh_ind, n)
             enddo
          enddo
       enddo

       call POP_HaloUpdate(SHR_STREAM_WORK, POP_haloClinic, &
            POP_gridHorzLocCenter,          &
            POP_fieldKindScalar, errorCode, &
            fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call exit_POP(sigAbort, subname /&
               &/ ': error updating halo for Ndep fields')
       endif

       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1, nblocks_clinic
          where (land_mask(:, :, iblock))
             input_forcing_data(:,:, ind_nh4_flux,iblock) = SHR_STREAM_WORK(:, :, iblock)
          endwhere
       enddo
       !$OMP END PARALLEL DO

    endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------

    if (din_riv_flux%has_data) then
       if (thour00 >= din_riv_flux%data_update) then
          tracer_data_names = din_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DIN River Flux'
          call update_forcing_data(din_riv_flux%data_time,    &
               din_riv_flux%data_time_min_loc,  din_riv_flux%interp_type,    &
               din_riv_flux%data_next,          din_riv_flux%data_update,    &
               din_riv_flux%data_type,          din_riv_flux%data_inc,       &
               din_riv_flux%DATA(:, :, :, :, 1:12), din_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               din_riv_flux%filename,           din_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            din_riv_flux%DATA(:, :, :, :, 1:12), &
            din_riv_flux%data_time,         din_riv_flux%interp_type, &
            din_riv_flux%data_time_min_loc, din_riv_flux%interp_freq, &
            din_riv_flux%interp_inc,        din_riv_flux%interp_next, &
            din_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_din_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (dip_riv_flux%has_data) then
       if (thour00 >= dip_riv_flux%data_update) then
          tracer_data_names = dip_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DIP River Flux'
          call update_forcing_data(dip_riv_flux%data_time,    &
               dip_riv_flux%data_time_min_loc,  dip_riv_flux%interp_type,    &
               dip_riv_flux%data_next,          dip_riv_flux%data_update,    &
               dip_riv_flux%data_type,          dip_riv_flux%data_inc,       &
               dip_riv_flux%DATA(:, :, :, :, 1:12), dip_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dip_riv_flux%filename,           dip_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dip_riv_flux%DATA(:, :, :, :, 1:12), &
            dip_riv_flux%data_time,         dip_riv_flux%interp_type, &
            dip_riv_flux%data_time_min_loc, dip_riv_flux%interp_freq, &
            dip_riv_flux%interp_inc,        dip_riv_flux%interp_next, &
            dip_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_dip_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (don_riv_flux%has_data) then
       if (thour00 >= don_riv_flux%data_update) then
          tracer_data_names = don_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DON River Flux'
          call update_forcing_data(don_riv_flux%data_time,    &
               don_riv_flux%data_time_min_loc,  don_riv_flux%interp_type,    &
               don_riv_flux%data_next,          don_riv_flux%data_update,    &
               don_riv_flux%data_type,          don_riv_flux%data_inc,       &
               don_riv_flux%DATA(:, :, :, :, 1:12), don_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               don_riv_flux%filename,           don_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            don_riv_flux%DATA(:, :, :, :, 1:12), &
            don_riv_flux%data_time,         don_riv_flux%interp_type, &
            don_riv_flux%data_time_min_loc, don_riv_flux%interp_freq, &
            don_riv_flux%interp_inc,        don_riv_flux%interp_next, &
            don_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_don_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (dop_riv_flux%has_data) then
       if (thour00 >= dop_riv_flux%data_update) then
          tracer_data_names = dop_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(dop_riv_flux%data_time,    &
               dop_riv_flux%data_time_min_loc,  dop_riv_flux%interp_type,    &
               dop_riv_flux%data_next,          dop_riv_flux%data_update,    &
               dop_riv_flux%data_type,          dop_riv_flux%data_inc,       &
               dop_riv_flux%DATA(:, :, :, :, 1:12), dop_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dop_riv_flux%filename,           dop_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dop_riv_flux%DATA(:, :, :, :, 1:12), &
            dop_riv_flux%data_time,         dop_riv_flux%interp_type, &
            dop_riv_flux%data_time_min_loc, dop_riv_flux%interp_freq, &
            dop_riv_flux%interp_inc,        dop_riv_flux%interp_next, &
            dop_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_dop_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (dsi_riv_flux%has_data) then
       if (thour00 >= dsi_riv_flux%data_update) then
          tracer_data_names = dsi_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DSI River Flux'
          call update_forcing_data(dsi_riv_flux%data_time,    &
               dsi_riv_flux%data_time_min_loc,  dsi_riv_flux%interp_type,    &
               dsi_riv_flux%data_next,          dsi_riv_flux%data_update,    &
               dsi_riv_flux%data_type,          dsi_riv_flux%data_inc,       &
               dsi_riv_flux%DATA(:, :, :, :, 1:12), dsi_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dsi_riv_flux%filename,           dsi_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dsi_riv_flux%DATA(:, :, :, :, 1:12), &
            dsi_riv_flux%data_time,         dsi_riv_flux%interp_type, &
            dsi_riv_flux%data_time_min_loc, dsi_riv_flux%interp_freq, &
            dsi_riv_flux%interp_inc,        dsi_riv_flux%interp_next, &
            dsi_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_dsi_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (dfe_riv_flux%has_data) then
       if (thour00 >= dfe_riv_flux%data_update) then
          tracer_data_names = dfe_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DFE River Flux'
          call update_forcing_data(dfe_riv_flux%data_time,    &
               dfe_riv_flux%data_time_min_loc,  dfe_riv_flux%interp_type,    &
               dfe_riv_flux%data_next,          dfe_riv_flux%data_update,    &
               dfe_riv_flux%data_type,          dfe_riv_flux%data_inc,       &
               dfe_riv_flux%DATA(:, :, :, :, 1:12), dfe_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dfe_riv_flux%filename,           dfe_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dfe_riv_flux%DATA(:, :, :, :, 1:12), &
            dfe_riv_flux%data_time,         dfe_riv_flux%interp_type, &
            dfe_riv_flux%data_time_min_loc, dfe_riv_flux%interp_freq, &
            dfe_riv_flux%interp_inc,        dfe_riv_flux%interp_next, &
            dfe_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_dfe_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (dic_riv_flux%has_data) then
       if (thour00 >= dic_riv_flux%data_update) then
          tracer_data_names = dic_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DIC River Flux'
          call update_forcing_data(dic_riv_flux%data_time,    &
               dic_riv_flux%data_time_min_loc,  dic_riv_flux%interp_type,    &
               dic_riv_flux%data_next,          dic_riv_flux%data_update,    &
               dic_riv_flux%data_type,          dic_riv_flux%data_inc,       &
               dic_riv_flux%DATA(:, :, :, :, 1:12), dic_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dic_riv_flux%filename,           dic_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dic_riv_flux%DATA(:, :, :, :, 1:12), &
            dic_riv_flux%data_time,         dic_riv_flux%interp_type, &
            dic_riv_flux%data_time_min_loc, dic_riv_flux%interp_freq, &
            dic_riv_flux%interp_inc,        dic_riv_flux%interp_next, &
            dic_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_dic_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (alk_riv_flux%has_data) then
       if (thour00 >= alk_riv_flux%data_update) then
          tracer_data_names = alk_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'ALK River Flux'
          call update_forcing_data(alk_riv_flux%data_time,    &
               alk_riv_flux%data_time_min_loc,  alk_riv_flux%interp_type,    &
               alk_riv_flux%data_next,          alk_riv_flux%data_update,    &
               alk_riv_flux%data_type,          alk_riv_flux%data_inc,       &
               alk_riv_flux%DATA(:, :, :, :, 1:12), alk_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               alk_riv_flux%filename,           alk_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            alk_riv_flux%DATA(:, :, :, :, 1:12), &
            alk_riv_flux%data_time,         alk_riv_flux%interp_type, &
            alk_riv_flux%data_time_min_loc, alk_riv_flux%interp_freq, &
            alk_riv_flux%interp_inc,        alk_riv_flux%interp_next, &
            alk_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_alk_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    if (doc_riv_flux%has_data) then
       if (thour00 >= doc_riv_flux%data_update) then
          tracer_data_names = doc_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'PP River Flux'
          call update_forcing_data(doc_riv_flux%data_time,    &
               doc_riv_flux%data_time_min_loc,  doc_riv_flux%interp_type,    &
               doc_riv_flux%data_next,          doc_riv_flux%data_update,    &
               doc_riv_flux%data_type,          doc_riv_flux%data_inc,       &
               doc_riv_flux%DATA(:, :, :, :, 1:12), doc_riv_flux%data_renorm, &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               doc_riv_flux%filename,           doc_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            doc_riv_flux%DATA(:, :, :, :, 1:12), &
            doc_riv_flux%data_time,         doc_riv_flux%interp_type, &
            doc_riv_flux%data_time_min_loc, doc_riv_flux%interp_freq, &
            doc_riv_flux%interp_inc,        doc_riv_flux%interp_next, &
            doc_riv_flux%interp_last,       0)
       input_forcing_data(:,:, ind_doc_riv_flux,:) = INTERP_WORK(:, :, :, 1)
    endif

    call timer_stop(ecosys_pre_sflux_timer)

  end subroutine ecosys_driver_read_sflux

  !*****************************************************************************

  subroutine ecosys_driver_ciso_init_tracers(&
       init_ts_file_fmt, &
       read_restart_filename, &
       tracer_d_module, &
       tracer_module, &
       errorCode)

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initialize ecosys_ciso tracer module. This involves setting metadata, reading
    ! the module namelist, setting initial conditions, setting up forcing,
    ! and defining additional tavg variables.
    !---------------------------------------------------------------------

    use marbl_interface_types , only : tracer_field => marbl_tracer_metadata_type
    use marbl_share_mod       , only : ciso_comp_surf_avg_freq_iopt
    use marbl_share_mod       , only : ciso_comp_surf_avg_freq
    use marbl_share_mod       , only : ciso_use_nml_surf_vals
    use marbl_share_mod       , only : ciso_surf_avg_di13c_const
    use marbl_share_mod       , only : ciso_surf_avg_di14c_const
    use marbl_share_mod       , only : ciso_tracer_init_ext
    use marbl_share_mod       , only : ciso_lecovars_full_depth_tavg
    use marbl_share_mod       , only : ciso_init_ecosys_option
    use marbl_share_mod       , only : ciso_init_ecosys_init_file
    use marbl_share_mod       , only : ciso_init_ecosys_init_file_fmt
    use marbl_share_mod       , only : ciso_comp_surf_avg_flag
    use marbl_parms           , only : di13c_ind
    use marbl_parms           , only : di14c_ind
    use passive_tracer_tools  , only : comp_surf_avg
    use passive_tracer_tools  , only : extract_surf_avg
    use passive_tracer_tools  , only : rest_read_tracer_block
    use passive_tracer_tools  , only : file_read_tracer_block
    use passive_tracer_tools  , only : read_field
    use time_management       , only : check_time_flag
    use time_management       , only : eval_time_flag
    use grid                  , only : fill_points
    use grid                  , only : n_topo_smooth
    use grid                  , only : n_topo_smooth
    use grid                  , only : KMT
    use prognostic            , only : curtime
    use prognostic            , only : oldtime
    use io_tools              , only : document

    implicit none

    character (*)       , intent(in)    :: init_ts_file_fmt      ! format (bin or nc) for input file
    character (*)       , intent(in)    :: read_restart_filename ! file name for restart file
    type (tracer_field) , intent(inout) :: tracer_d_module(:)    ! descriptors for each tracer
    real (r8)           , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)
    integer (POP_i4)    , intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ciso_mod:ecosys_ciso_init'

    integer (int_kind) :: &
         auto_ind,                        & ! autotroph functional group index
         n,                               & ! index for looping over tracers
         k,                               & ! index for looping over depth levels
         l,                               & ! index for looping over time levels
         ind,                             & ! tracer index for tracer name from namelist
         iblock                             ! index for looping over blocks

    character (char_len) :: &
         ecosys_ciso_restart_filename  ! modified file name for restart file
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  initialize tracers
    !-----------------------------------------------------------------------

    select case (ciso_init_ecosys_option)

    case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

       ecosys_ciso_restart_filename = char_blank

       if (ciso_init_ecosys_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call document(subname, 'no restart file to read ciso vars from')
             call exit_POP(sigAbort, 'stopping in ' /&
                  &/ subname)
          endif
          ecosys_ciso_restart_filename = read_restart_filename
          ciso_init_ecosys_init_file_fmt = init_ts_file_fmt

       else  ! do not read from TS restart file

          ecosys_ciso_restart_filename = trim(ciso_init_ecosys_init_file)

       endif

       call rest_read_tracer_block(         &
            ciso_init_ecosys_init_file_fmt, &
            ecosys_ciso_restart_filename,   &
            tracer_d_module,                &
            tracer_module)

       if (ciso_use_nml_surf_vals) then
          surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
          surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
       else
          call extract_surf_avg(               &
               ciso_init_ecosys_init_file_fmt, &
               ecosys_ciso_restart_filename,   &
               ecosys_ciso_tracer_cnt,         &
               ciso_vflux_flag,                &
               ciso_ind_name_table,            &
               ciso_surf_avg)
       endif

       ! evaluate time_flag(ciso_comp_surf_avg_flag)%value via time_to_do
       call eval_time_flag(ciso_comp_surf_avg_flag) 

       if (check_time_flag(ciso_comp_surf_avg_flag)) then
          call comp_surf_avg(                    &
               tracer_module(:,:,1,:,oldtime,:), &
               tracer_module(:,:,1,:,curtime,:), &
               ecosys_ciso_tracer_cnt,           &
               ciso_vflux_flag,                  &
               ciso_surf_avg)
       end if

    case ('file', 'ccsm_startup')
       call document(subname, 'ciso vars being read from separate files')

       call file_read_tracer_block(         &
            ciso_init_ecosys_init_file_fmt, &
            ciso_init_ecosys_init_file,     &
            tracer_d_module,                &
            ciso_ind_name_table,            &
            ciso_tracer_init_ext,           &
            tracer_module)

       if (n_topo_smooth > 0) then
          errorCode = POP_Success
          do n = 1, ecosys_ciso_tracer_cnt
             do k=1,km
                call fill_points(k,tracer_module(:,:,k,n,oldtime,:), errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(oldtime)')
                   return
                endif

                call fill_points(k,tracer_module(:,:,k,n,curtime,:), errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(newtime)')
                   return
                endif

             enddo
          enddo
       endif

       if (ciso_use_nml_surf_vals) then
          ciso_surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
          ciso_surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
       else
          call comp_surf_avg(&
               tracer_module(:,:,1,:,oldtime,:), &
               tracer_module(:,:,1,:,curtime,:), &
               ecosys_ciso_tracer_cnt,           &
               ciso_vflux_flag,ciso_surf_avg)
       endif

    case default
       call document(subname, 'ciso_init_ecosys_option', ciso_init_ecosys_option)
       call exit_POP(sigAbort, 'unknown ciso_init_ecosys_option')

    end select

    !-----------------------------------------------------------------------
    !  apply land mask to tracers - FIXME
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(iblock,n,k)
    do iblock=1,nblocks_clinic
       do n = 1,ecosys_ciso_tracer_cnt
          do k = 1,km
             where (.not. land_mask(:,:,iblock) .or. k > KMT(:,:,iblock))
                TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
                TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
             end where
          end do
       end do
    enddo
    !$OMP END PARALLEL DO

  end subroutine ecosys_driver_ciso_init_tracers

  !***********************************************************************

  subroutine ecosys_driver_ciso_init_atm_D13_D14

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Initialize surface flux computations for the ecosys_ciso tracer module.
    !  Includes reading CO2 and D13C and D14C data from file if option file is used
    !---------------------------------------------------------------------

    use communicate     , only : my_task, master_task
    use io_types        , only : stdout
    use exit_mod        , only : exit_POP, sigAbort
    use constants       , only : blank_fmt      
    use constants       , only : delim_fmt      
    use constants       , only : ndelim_fmt     
    use marbl_share_mod , only : ciso_atm_d13c_opt          
    use marbl_share_mod , only : ciso_atm_d14c_opt          
    use marbl_share_mod , only : ciso_atm_d13c_const
    use marbl_share_mod , only : ciso_atm_d14c_const

    implicit none

    !-------------------------------------------------------------------------
    !     Set D13C data source
    !-------------------------------------------------------------------------

    select case (ciso_atm_d13c_opt)
    case ('const')
       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Using constant D13C values of ',ciso_atm_d13c_const
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
    case('file')
       call ecosys_driver_ciso_read_atm_D13C_data ! READ in D13C data from file
    case default
       call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt in ecosys_ciso_init_sflux')
    end select

    !-------------------------------------------------------------------------
    !     Set D14C data source
    !-------------------------------------------------------------------------

    select case (ciso_atm_d14c_opt)
    case ('const')
       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Using constant D14C values of ',ciso_atm_d14c_const
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
    case('file')
       call ecosys_driver_ciso_read_atm_D14C_data ! READ in D14C data from files
    case default
       call exit_POP(sigAbort, 'unknown ciso_atm_d14c_opt in ecosys_ciso_init_sflux')
    end select

  end subroutine ecosys_driver_ciso_init_atm_D13_D14

  !***********************************************************************

  subroutine ecosys_driver_ciso_read_atm_D13C_data()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Read atmospheric D13C [permil] data from file
    !
    !  Have the master_task do the following :
    !     1) get length of data
    !     2) allocate memory for data
    !     3) read in data, checking for consistent lengths
    !  Then, outside master_task conditional
    !     1) broadcast length of data
    !     2) have non-mastertasks allocate memory for data
    !     3) broadcast data
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : ciso_atm_d13c_data_nbval                                                       
    use marbl_share_mod , only : ciso_atm_d13c_data                                                             
    use marbl_share_mod , only : ciso_atm_d13c_data_yr                                                          
    use marbl_share_mod , only : ciso_atm_d13c_const                                                            
    use marbl_share_mod , only : ciso_atm_d13c_opt                                                              
    use marbl_share_mod , only : ciso_atm_d13c_filename                                                         
    use broadcast       , only : broadcast_array
    use broadcast       , only : broadcast_scalar
    use communicate     , only : master_task
    use communicate     , only : my_task
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt
    use io_types        , only : nml_in
    use io_tools        , only : document

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: sub_name = 'ecosys_driver:ecosys_driver_ciso_read_atm_D13C_data'

    integer (int_kind) ::    &
         stat,                 &  ! i/o status code
         irec,                 &  ! counter for looping
         skiplines,            &  ! number of comment lines at beginning of ascii file
         il                       ! looping index

    character (char_len) :: &
         sglchr                   ! variable to read characters from file into

    !-----------------------------------------------------------------------
    !     READ in D13C data from file
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*)'ciso: Using varying D13C values from file ',trim(ciso_atm_d13c_filename)
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       open (nml_in, file=ciso_atm_d13c_filename, status='old',iostat=stat)
       if (stat /= 0) then
          write(stdout,fmt=*) 'open failed'
          go to 99
       endif
       read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d13c_data_nbval
       if (stat /= 0) then
          write(stdout,fmt=*) '1st line read failed'
          go to 99
       endif
       allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
       allocate(ciso_atm_d13c_data   (ciso_atm_d13c_data_nbval))
       do irec=1,skiplines
          read(nml_in,FMT=*,iostat=stat) sglchr
          if (stat /= 0) then
             write(stdout,fmt=*) 'skipline read failed'
             go to 99
          endif
       enddo
       do irec=1,ciso_atm_d13c_data_nbval
          read(nml_in,FMT=*,iostat=stat) ciso_atm_d13c_data_yr(irec), ciso_atm_d13c_data(irec)
          if (stat /= 0) then
             write(stdout,fmt=*) 'data read failed'
             go to 99
          endif
       enddo
       close(nml_in)
    endif

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
         &/ sub_name)

    !---------------------------------------------------------------------
    !     Need to allocate and broadcast the variables to other tasks beside master-task
    !---------------------------------------------------------------------

    call broadcast_scalar(ciso_atm_d13c_data_nbval,master_task)

    if (my_task /= master_task) then
       allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
       allocate(ciso_atm_d13c_data(ciso_atm_d13c_data_nbval))
    endif

    call broadcast_array(ciso_atm_d13c_data   , master_task)
    call broadcast_array(ciso_atm_d13c_data_yr, master_task)

  end subroutine ecosys_driver_ciso_read_atm_D13C_data

  !***********************************************************************

  subroutine ecosys_driver_ciso_read_atm_D14C_data

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Read atmospheric D14C data from file
    !
    !  Have the master_task do the following :
    !     1) get length of data
    !     2) allocate memory for data
    !     3) read in data, checking for consistent lengths
    !  Then, outside master_task conditional
    !     1) broadcast length of data
    !     2) have non-mastertasks allocate memory for data
    !     3) broadcast data
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : ciso_atm_d14c_data_nbval                                                       
    use marbl_share_mod , only : ciso_atm_d14c_data                                                             
    use marbl_share_mod , only : ciso_atm_d14c_data_yr                                                          
    use marbl_share_mod , only : ciso_atm_d14c_const                                                            
    use marbl_share_mod , only : ciso_atm_d14c_opt                                                              
    use marbl_share_mod , only : ciso_atm_d14c_filename                                                         
    use broadcast       , only : broadcast_array
    use broadcast       , only : broadcast_scalar
    use communicate     , only : master_task
    use communicate     , only : my_task
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt
    use io_types        , only : nml_in
    use io_tools        , only : document

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: sub_name = 'ciso_read_atm_D14C_data:ciso_read_atm_D14C_data'

    integer (int_kind) ::      &
         stat,                   &  ! i/o status code
         irec,                   &  ! counter for looping
         skiplines,              &  ! number of comment lines at beginning of ascii file
         il                         ! looping index

    character (char_len) ::  &
         sglchr                     ! variable to read characters from file into

    integer (int_kind) :: &
         ciso_atm_d14c_data_nbval_tmp

    logical (log_kind) :: &
         nbval_mismatch

    !-----------------------------------------------------------------------
    !     ensure that three datafiles have same number of entries
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       write(stdout,*)'ciso DIC14 calculation: Using varying C14 values from files'
       do il=1,3
          write(stdout,*) trim(ciso_atm_d14c_filename(il))
       enddo
       nbval_mismatch = .false.
       do il=1,3
          open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
          if (stat /= 0) then
             write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
          if (stat /= 0) then
             write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          close(nml_in)
          if (il == 1) then
             ciso_atm_d14c_data_nbval = ciso_atm_d14c_data_nbval_tmp
          else
             if (ciso_atm_d14c_data_nbval /= ciso_atm_d14c_data_nbval_tmp) nbval_mismatch = .true.
          endif
       enddo
    endif

    call broadcast_scalar(nbval_mismatch, master_task)
    if (nbval_mismatch) then
       call document(sub_name, 'D14C data files must all have the same number of values')
       call exit_POP(sigAbort, 'stopping in ' /&
            &/ sub_name)
    endif

    call broadcast_scalar(ciso_atm_d14c_data_nbval, master_task)
    allocate(ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,3))
    allocate(ciso_atm_d14c_data   (ciso_atm_d14c_data_nbval,3))

    !-----------------------------------------------------------------------
    !     READ in C14 data from files - three files, for SH, EQ, NH
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       do il=1,3
          open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
          if (stat /= 0) then
             write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
          if (stat /= 0) then
             write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          do irec=1,skiplines
             read(nml_in,FMT=*,iostat=stat) sglchr
             if (stat /= 0) then
                write(stdout,fmt=*) 'skipline read failed for ', trim(ciso_atm_d14c_filename(il))
                go to 99
             endif
          enddo
          do irec=1,ciso_atm_d14c_data_nbval
             read(nml_in,FMT=*,iostat=stat) ciso_atm_d14c_data_yr(irec,il), ciso_atm_d14c_data(irec,il)
             if (stat /= 0) then
                write(stdout,fmt=*) 'data read failed for ', trim(ciso_atm_d14c_filename(il))
                go to 99
             endif
          enddo
          close(nml_in)
       enddo
    endif

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
         &/ sub_name)

    !---------------------------------------------------------------------
    ! Broadcast the variables to other tasks beside master_task
    !---------------------------------------------------------------------

    call broadcast_array(ciso_atm_d14c_data   , master_task)
    call broadcast_array(ciso_atm_d14c_data_yr, master_task)

  end subroutine ecosys_driver_ciso_read_atm_D14C_data

  !***********************************************************************

  subroutine ecosys_driver_ciso_update_atm_D13C_D14C ()

    ! Updates module variables D13C and D14C (for atmospheric ratios)

    use marbl_share_mod   , only : ciso_atm_d13c_opt          
    use marbl_share_mod   , only : ciso_atm_d14c_opt          
    use marbl_share_mod   , only : ciso_atm_d13c_const          
    use marbl_share_mod   , only : ciso_atm_d14c_const          
    use grid              , only : TAREA
    use domain            , only : blocks_clinic
    use blocks            , only : get_block
    use global_reductions , only : global_sum

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    type (block) :: &
         this_block      ! block info for the current block

    real (r8), dimension(nx_block,ny_block) :: &
         work1, &! local work space
         tfact   ! factor for normalizing sums

    integer (int_kind) :: &
         ib,ie,jb,je, &
         iblock  ! index for looping over blocks

    real (r8), dimension(max_blocks_clinic) :: &
         d14c_local_sums, & ! array for holding block sums when calculating global D14C
         tarea_local_sums   ! array for holding block sums of TAREA when calculating global D14C

    real (r8) :: &
         d14c_sum_tmp,  & ! temp for local sum of D14C
         tarea_sum_tmp    ! temp for local sum of TAREA

    !-----------------------------------------------------------------------

    work1(:,:) = c0
    d14c_local_sums(:)  = c0
    tarea_local_sums(:) = c0
    
    !-----------------------------------------------------------------------
    ! Loop over blocks
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic

       !-----------------------------------------------------------------------
       !  Set D13C (constant or from files read in _init) and put on global grid
       !-----------------------------------------------------------------------

       select case (ciso_atm_d13c_opt)
       case ('const')
          D13C(:,:,:) = ciso_atm_d13c_const
       case ('file')
          call ecosys_driver_ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c(iblock), D13C(:,:,iblock))
       case default
          call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt in ecosys_ciso_set_sflux')
       end select

       !-----------------------------------------------------------------------
       !  Set D14C (constant or from files read in _init) and put on global grid
       !-----------------------------------------------------------------------

       select case (ciso_atm_d14c_opt)
       case ('const')
          D14C(:,:,:) = ciso_atm_d14c_const
       case ('file')
          call ecosys_driver_ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c(iblock), D14C(:,:,iblock))
       case default
          call exit_POP(sigAbort, 'unknown ciso_atm_d14c_opt in ecosys_ciso_set_sflux')
       end select

       !-----------------------------------------------------------------------
       ! Save local D14C field for making global mean after end of iblock loop
       !-----------------------------------------------------------------------

       this_block = get_block(blocks_clinic(iblock),iblock)
       ib = this_block%ib
       ie = this_block%ie
       jb = this_block%jb
       je = this_block%je

       where (land_mask(:,:,iblock))
          tfact(:,:) = TAREA(:,:,iblock)
       elsewhere
          tfact(:,:) = 0.0_r8
       endwhere

       work1(:,:) = D14C(:,:,iblock) * tfact(:,:)
       d14c_local_sums(iblock)  = sum(work1(ib:ie,jb:je))
       tarea_local_sums(iblock) = sum(tfact(ib:ie,jb:je))

    end do

    !-----------------------------------------------------------------------
    ! Compute D14C making global mean
    !-----------------------------------------------------------------------

    d14c_sum_tmp  = sum(d14c_local_sums)
    tarea_sum_tmp = sum(tarea_local_sums)

    d14c_glo_avg  = global_sum(d14c_sum_tmp ,distrb_clinic) / global_sum(tarea_sum_tmp,distrb_clinic)

  end subroutine ecosys_driver_ciso_update_atm_D13C_D14C

  !***********************************************************************

  subroutine ecosys_driver_ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c, D13C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of d13c when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate data values to current model time step
    !  2. Spatial patern of D13Cis the same everywhere (90 S - 90 N)
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : ciso_atm_model_year                                                            
    use marbl_share_mod , only : ciso_atm_data_year                                                             
    use marbl_share_mod , only : ciso_atm_d13c_data_nbval                                                       
    use marbl_share_mod , only : ciso_atm_d13c_data                                                             
    use marbl_share_mod , only : ciso_atm_d13c_data_yr                                                          
    use marbl_share_mod , only : ciso_atm_d13c_const                                                            
    use marbl_share_mod , only : ciso_atm_d13c_opt                                                              
    use marbl_share_mod , only : ciso_atm_d13c_filename                                                         
    use time_management , only : days_in_year
    use time_management , only : frac_day
    use time_management , only : iday_of_year
    use time_management , only : iyear
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt
    use communicate     , only : master_task
    use communicate     , only : my_task

    implicit none

    ! note that ciso_data_ind_d13c is always strictly less than the length
    ! of the data and is initialized to -1 before the first call

    integer (int_kind) , intent(in)  :: iblock                  ! block index
    integer (int_kind) , intent(out) :: ciso_data_ind_d13c      ! inex for the data for current timestep,
    real (r8)          , intent(out) :: D13C(nx_block,ny_block) ! atmospheric D13C (permil)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: &
         i, j              ! loop indices

    real (r8) :: &
         model_date,     & ! date of current model timestep
         mapped_date,    & ! model_date mapped to data timeline
         weight            ! weighting for temporal interpolation
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Generate mapped_date and check to see if it is too large.
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year

    if (mapped_date >= ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval)) then
       call exit_POP(sigAbort, 'ciso: Model date maps to date after end of D13C data in file.')
    endif

    !--------------------------------------------------------------------------------------------------------------
    !  Set atmospheric D13C to first value in record for years before record begins
    !--------------------------------------------------------------------------------------------------------------

    if (mapped_date < ciso_atm_d13c_data_yr(1)) then
       D13C = ciso_atm_d13c_data(1)
       ciso_data_ind_d13c = 1
       if(my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Mapped date less than start of D13C data --> using first value in D13C data file'
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
       return
    endif

    !-----------------------------------------------------------------------
    !  On first time step, perform linear search to find data_ind_d13c
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d13c == -1) then
       do ciso_data_ind_d13c = ciso_atm_d13c_data_nbval-1,1,-1
          if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) exit
       end do
    endif

    !-----------------------------------------------------------------------
    !  See if ciso_data_ind_d13c needs to be updated,
    !  but do not set it to atm_d13c_data_nbval.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d13c < ciso_atm_d13c_data_nbval-1) then
       if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1)) ciso_data_ind_d13c = ciso_data_ind_d13c + 1
    endif

    !-----------------------------------------------------------------------
    !  Generate hemisphere values for current time step.
    !-----------------------------------------------------------------------

    weight = (mapped_date - ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) &
         / (ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1) - ciso_atm_d13c_data_yr(ciso_data_ind_d13c))

    D13C = weight * ciso_atm_d13c_data(ciso_data_ind_d13c+1) + (c1-weight) * ciso_atm_d13c_data(ciso_data_ind_d13c)

  end subroutine ecosys_driver_ciso_comp_varying_D13C

  !***********************************************************************

  subroutine ecosys_driver_ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c, D14C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of CO2 when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate hemispheric values to current time step
    !  2. Make global field of D14C, determined by:
    !   -Northern Hemisphere value is used for 20N - 90 N
    !   -Southern Hemisphere value is used for 20 S - 90 S
    !   -Equator value is used for 20 S- 20 N

    use marbl_share_mod, only : ciso_atm_model_year                                                            
    use marbl_share_mod, only : ciso_atm_data_year                                                             
    use marbl_share_mod, only : ciso_atm_d14c_data_nbval                                                       
    use marbl_share_mod, only : ciso_atm_d14c_data                                                             
    use marbl_share_mod, only : ciso_atm_d14c_data_yr                                                          
    use marbl_share_mod, only : ciso_atm_d14c_const                                                            
    use marbl_share_mod, only : ciso_atm_d14c_opt                                                              
    use marbl_share_mod, only : ciso_atm_d14c_filename                                                         
    use time_management, only : days_in_year
    use time_management, only : frac_day
    use time_management, only : iday_of_year
    use time_management, only : iyear
    use grid           , only : TLATD
    use communicate    , only : master_task
    use communicate    , only : my_task
    
    implicit none

    !  note that data_ind is always strictly less than the length of D14C data
    !  and is initialized to -1 before the first call

    integer (int_kind) , intent(in)  :: iblock                  ! block index
    integer (int_kind) , intent(out) :: ciso_data_ind_d14c      ! data_ind_d14c is the index into data for current timestep,
    real (r8)          , intent(out) :: D14C(nx_block,ny_block) ! atmospheric delta C14 in permil on global grid

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: &
         i, j, il        ! loop indices

    real (r8) :: &
         model_date,      & ! date of current model timestep
         mapped_date,     & ! model_date mapped to data timeline
         weight,          & ! weighting for temporal interpolation
         d14c_curr_sh,    & ! current atmospheric D14C value for SH (interpolated from data to model date)
         d14c_curr_nh,    & ! current atmospheric D14C value for NH (interpolated from data to model date)
         d14c_curr_eq       ! current atmospheric D14C value for EQ (interpolated from data to model date)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Generate mapped_date and check to see if it is too large.
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year
    do il=1,3
       if (mapped_date >= ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,il)) then
          call exit_POP(sigAbort, 'ciso: model date maps to date after end of D14C data in files.')
       endif
    enddo

    !--------------------------------------------------------------------------------------------------------------
    !  Set atmospheric D14C concentrations to zero before D14C record begins
    !--------------------------------------------------------------------------------------------------------------

    if (mapped_date < ciso_atm_d14c_data_yr(1,1)) then
       D14C = c0
       ciso_data_ind_d14c = 1
       if(my_task == master_task) then
          write(stdout,*)'ciso: Model date less than start of D14C data --> D14C=0'
       endif
       return
    endif

    !-----------------------------------------------------------------------
    !  On first time step, perform linear search to find data_ind_d14c.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d14c == -1) then
       do ciso_data_ind_d14c = ciso_atm_d14c_data_nbval-1,1,-1
          if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) exit
       end do
    endif

    !-----------------------------------------------------------------------
    !  See if data_ind_d14c need to be updated,
    !  but do not set it to atm_co2_data_nbval.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d14c < ciso_atm_d14c_data_nbval-1) then
       if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1))  then
          ciso_data_ind_d14c = ciso_data_ind_d14c + 1
       endif
    endif
    !
    !-----------------------------------------------------------------------
    !  Generate hemisphere values for current time step.
    !-----------------------------------------------------------------------

    weight = (mapped_date - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) &
         / (ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1) - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1))

    d14c_curr_sh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,1) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,1)
    d14c_curr_eq = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,2) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,2)
    d14c_curr_nh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,3) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,3)

    !-----------------------------------------------------------------------
    !  Merge hemisphere values for D14C
    !      -Northern Hemisphere value is used for >20N - 90 N
    !      -Southern Hemisphere value is used for >20 S - 90 S
    !      -Equatorial value is used for 20 S to 20 N
    !-----------------------------------------------------------------------

    do j = 1, ny_block
       do i = 1, nx_block
          if (TLATD(i,j,iblock) < -20.0_r8) then
             D14C(i,j) = d14c_curr_sh
          else if (TLATD(i,j,iblock) > 20.0_r8) then
             D14C(i,j) = d14c_curr_nh
          else
             D14C(i,j) = d14c_curr_eq
          endif
       end do
    end do

  end subroutine ecosys_driver_ciso_comp_varying_D14C

  !*****************************************************************************

  subroutine ecosys_saved_state_constructor(this)

    class(ecosys_saved_state_type), intent(inout) :: this

    ! NOTE(bja, 2015-08) we are intentionally allocating memory here
    ! but not initializing because marbl ecosys won't know anything
    ! about the domain. initialization of elements will occur in marbl/ecosys_init
    allocate(this%dust_flux_in      (nx_block, ny_block,     max_blocks_clinic))
    allocate(this%ph_prev_3d        (nx_block, ny_block, km, max_blocks_clinic))
    allocate(this%ph_prev_alt_co2_3d(nx_block, ny_block, km, max_blocks_clinic))
    allocate(this%ph_surf           (nx_block, ny_block, max_blocks_clinic) )
    allocate(this%ph_surf_alt_co2   (nx_block, ny_block, max_blocks_clinic) )

  end subroutine ecosys_saved_state_constructor

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

