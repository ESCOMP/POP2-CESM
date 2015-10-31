! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_driver

  !BOP
  ! !MODULE: ecosys_driver

  ! !DESCRIPTION:
  !  This module provides support for the ecosystem module and dependend tracer modules
  !  The base model calls subroutines in passive_tracers, which then call
  !  this module if ecosys_on is true. Ecosys_driver then calls subroutines
  !  in individual ecosystem modules (so far ecosys_mod and ecosys_ciso_mod)
  !
  !  Written by: Alexandra Jahn, NCAR, Nov/Dec 2012


  ! !REVISION HISTORY:
  !  SVN:$Id: passive_tracers.F90 28439 2011-05-18 21:40:58Z njn01 $

  ! !USES:

  use POP_KindsMod
  use POP_ErrorMod
  use POP_IOUnitsMod

  use kinds_mod                 , only : r8, int_kind, log_kind, char_len, char_len_long
  use blocks                    , only : block, nx_block, ny_block
  use domain                    , only : nblocks_clinic
  use domain                    , only : distrb_clinic
  use domain_size               , only : max_blocks_clinic, km, nt
  use communicate               , only : my_task, master_task
  use prognostic                , only : TRACER, tracer_d, oldtime, curtime, newtime
  use forcing_shf               , only : SHF_QSW_RAW, SHF_QSW
  use io_types                  , only : stdout, nml_in, nml_filename, io_field_desc, datafile
  use exit_mod                  , only : sigAbort, exit_pop
  use io_tools                  , only : document
  use prognostic                , only : tracer_field
  use constants                 , only : c0, c1, p5, delim_fmt, char_blank, ndelim_fmt
  use passive_tracer_tools      , only : set_tracer_indices
  use passive_tracer_tools      , only : ind_name_pair
  use broadcast                 , only : broadcast_scalar

  use marbl_parms               , only : dic_ind
  use marbl_parms               , only : alk_ind
  use marbl_parms               , only : dic_alt_co2_ind

  use marbl_interface           , only : marbl_interface_class
  use marbl_interface           , only : marbl_sizes_type
  use marbl_interface           , only : marbl_driver_sizes_type

  use marbl_interface_constants , only : marbl_status_ok
  use marbl_interface_types     , only : marbl_status_type
  use marbl_interface_types     , only: marbl_diagnostics_type
  use marbl_interface_types     , only : marbl_saved_state_type

  ! NOTE(bja, 2014-12) all other uses of marbl/ecosys modules need to be removed!
  use marbl_share_mod           , only : autotroph_cnt, zooplankton_cnt

  use ecosys_constants, only : ecosys_tracer_cnt

  use ecosys_diagnostics_mod, only : ecosys_diag_cnt_2d
  use ecosys_diagnostics_mod, only : ecosys_diag_cnt_3d
  use ecosys_diagnostics_mod, only : auto_diag_cnt_2d
  use ecosys_diagnostics_mod, only : auto_diag_cnt_3d
  use ecosys_diagnostics_mod, only : zoo_diag_cnt_2d
  use ecosys_diagnostics_mod, only : zoo_diag_cnt_3d
  use ecosys_diagnostics_mod, only : part_diag_cnt_2d
  use ecosys_diagnostics_mod, only : part_diag_cnt_3d
  use ecosys_diagnostics_mod, only : forcing_diag_cnt

  use ecosys_tavg, only : ecosys_tavg_init
  use ecosys_tavg, only : ecosys_tavg_accumulate
  use ecosys_tavg, only : ecosys_tavg_accumulate_flux

  use ecosys_mod, only:            &
       ecosys_init_nml,            &
       ecosys_init_tracer_metadata, &
       ecosys_init_postnml,        &
       ecosys_tracer_ref_val,      &
       ecosys_set_sflux,           &
       ecosys_tavg_forcing,        &
       ecosys_set_interior

  use ecosys_ciso_mod, only:        &
       ecosys_ciso_tracer_cnt,      &
       ecosys_ciso_init,            &
       ecosys_ciso_tracer_ref_val,  &
       ecosys_ciso_set_sflux,       &
       ecosys_ciso_tavg_forcing,    &
       ecosys_ciso_set_interior,    &
       ecosys_ciso_write_restart

  use ecosys_restore_mod, only : ecosys_restore_type

  use timers, only : timer_start
  use timers, only : timer_stop
  use timers, only : get_timer

  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:

  public ::                          &
       ecosys_driver_tracer_cnt_init,  &
       ecosys_driver_init,             &
       ecosys_driver_set_interior,     &
       ecosys_driver_set_sflux,        &
       ecosys_driver_tracer_cnt,       &
       ecosys_driver_tracer_ref_val,   &
       ecosys_driver_tavg_forcing,     &
       ecosys_driver_write_restart,    &
       ecosys_driver_unpack_source_sink_terms

  private :: &
       ecosys_write_restart

  !EOP
  !BOC

  !-----------------------------------------------------------------------
  !  module variables required by forcing_passive_tracer
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       ecosys_driver_tracer_cnt

  !-----------------------------------------------------------------------
  !     index bounds of passive tracer module variables in TRACER
  !-----------------------------------------------------------------------

  integer (kind=int_kind) :: &
       ecosys_driver_ind_begin,  ecosys_driver_ind_end,  &
       ecosys_ind_begin,         ecosys_ind_end,         &
       ecosys_ciso_ind_begin,    ecosys_ciso_ind_end

  integer (int_kind) :: ecosys_interior_timer
  
  !--------------------------------------------------------------------
  ! removed from marbl_share because they are read from pop's
  ! namelist and passed into marbl (lmarginal_seas) or not used in
  ! marbl at all (tadvect).
  ! --------------------------------------------------------------------
  logical(log_kind) :: lmarginal_seas         ! Is ecosystem active in marginal seas ?
  character(char_len) :: ecosys_tadvect_ctype ! advection method for ecosys tracers


  !-----------------------------------------------------------------------
  !  data storage for interaction with the marbl bgc library
  !-----------------------------------------------------------------------
  type(marbl_interface_class) :: marbl
  type(marbl_sizes_type) :: marbl_sizes
  type(marbl_driver_sizes_type) :: marbl_driver_sizes
  type(marbl_status_type) :: marbl_status
  type(marbl_diagnostics_type), dimension(max_blocks_clinic) :: marbl_diagnostics
  ! FIXME(bja, 2015-08) this needs to go into marbl%private_data%saved_state !!!
  type(marbl_saved_state_type) :: marbl_saved_state


  type(ecosys_restore_type) :: ecosys_restore

  !-----------------------------------------------------------------------
  !  average surface tracer value related variables
  !  used as reference value for virtual flux computations
  !-----------------------------------------------------------------------
  ! NOTE (mvertens, 2015-10) the following has been moved from ecosys_mod module variables
  ! to ecosys_driver module variables as part of the marbilization of the ecosys initialization

  logical (log_kind) :: vflux_flag(ecosys_tracer_cnt) ! which tracers get virtual fluxes applied
  real (r8)          :: surf_avg(ecosys_tracer_cnt)   ! average surface tracer values
  integer (int_kind) :: comp_surf_avg_flag            ! time flag id for computing average surface tracer values

  real (r8), allocatable, target :: PH_PREV(:, :, :)         ! computed ph from previous time step
  real (r8), allocatable, target :: PH_PREV_ALT_CO2(:, :, :) ! computed ph from previous time step, alternative CO2

  type(ind_name_pair) :: ind_name_table(ecosys_tracer_cnt) !  derived type & parameter for tracer index lookup

  !EOC
  !***********************************************************************

contains

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_tracer_cnt_init
  ! !INTERFACE:

  subroutine ecosys_driver_tracer_cnt_init(ciso_on)

    ! !DESCRIPTION:
    !  Zero-level initialization of ecosys_driver,
    !  which involves setting the ecosys_driver_tracer_cnt
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:

    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    !EOP
    !BOC
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Determine ecosys_driver_tracer_cnt, depending on whether only ecosys
    !  or also other modules are on
    !-----------------------------------------------------------------------

    ecosys_driver_tracer_cnt = ecosys_tracer_cnt

    if (ciso_on) then
       ecosys_driver_tracer_cnt = ecosys_driver_tracer_cnt + ecosys_ciso_tracer_cnt
    end if


    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_tracer_cnt_init



  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_init
  ! !INTERFACE:

  subroutine ecosys_driver_init(ciso_on,init_ts_file_fmt,     &
       read_restart_filename,                         &
       tracer_d_module, TRACER_MODULE, tadvect_ctype, &
       errorCode)


    ! !DESCRIPTION:
    !  Initialize ecosys_driver passive tracers. This involves:
    !  1) setting ecosys and ecosys_ciso module index bounds
    !  2) calling ecosys and ecosys_ciso module init subroutine
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_interface_constants , only : marbl_nl_buffer_size
    use grid                      , only : n_topo_smooth
    use grid                      , only : fill_points
    use grid                      , only : REGION_MASK
    use grid                      , only : KMT
    use broadcast                 , only : broadcast_scalar
    use time_management           , only : eval_time_flag
    use passive_tracer_tools      , only : set_tracer_indices
    use passive_tracer_tools      , only : extract_surf_avg
    use passive_tracer_tools      , only : comp_surf_avg
    use passive_tracer_tools      , only : rest_read_tracer_block
    use passive_tracer_tools      , only : field_exists_in_file
    use passive_tracer_tools      , only : tracer_read
    use passive_tracer_tools      , only : file_read_tracer_block
    use passive_tracer_tools      , only : read_field

    ! !INPUT PARAMETERS:

    logical (kind=log_kind) , intent(in) :: ciso_on                ! ecosys_ciso on
    character (*)           , intent(in) :: init_ts_file_fmt       ! format (bin or nc) for input file
    character (*)           , intent(in) :: read_restart_filename  ! file name for restart file

    ! !INPUT/OUTPUT PARAMETERS:

    type (tracer_field) , intent(inout) :: tracer_d_module(:)   ! descriptors for each tracer
    real (r8)           , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)

    ! !OUTPUT PARAMETERS:

    character (char_len), intent(out) :: tadvect_ctype(:)     ! advection method for ecosys tracers
    integer (POP_i4)    , intent(out) :: errorCode

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter         :: subname = 'ecosys_driver:ecosys_driver_init'
    integer (int_kind)              :: cumulative_nt, n, bid, k
    integer (int_kind)              :: nml_error          ! error flag for nml read
    integer (int_kind)              :: iostat             ! io status flag
    character (char_len)            :: sname, lname, units, coordinates
    character (4)                   :: grid_loc
    character(marbl_nl_buffer_size) :: nl_buffer
    character(char_len_long)        :: ioerror_msg

    ! NOTE (mvertens, 2015-10) variables moved from ecosys_init
    character (char_len) ::  ecosys_restart_filename            ! modified file name for restart file

    ! NOTE (mvertens, 2015-10) namelist variables moved from ecosys_mod to ecosys_driver
    real (r8)           :: surf_avg_dic_const, surf_avg_alk_const
    logical (log_kind)  :: use_nml_surf_vals         ! do namelist surf values override values from restart file
    character(char_len) :: init_ecosys_option        ! namelist option for initialization of bgc
    character(char_len) :: init_ecosys_init_file     ! filename for option 'file'
    character(char_len) :: init_ecosys_init_file_fmt ! file format for option 'file'
    type(tracer_read)   :: tracer_init_ext(ecosys_tracer_cnt) ! namelist variable for initializing tracers

    !-----------------------------------------------------------------------
    !  read in ecosys_driver namelist, to set namelist parameters that
    !  should be the same for all ecosystem-related modules
    !-----------------------------------------------------------------------

    namelist /ecosys_driver_nml/ &
         lmarginal_seas, ecosys_tadvect_ctype

    errorCode = POP_Success

    lmarginal_seas        = .true.
    ecosys_tadvect_ctype  = 'base_model'

    nl_buffer = ''
    if (my_task == master_task) then
       ! read the namelist file into a buffer.
       open(unit=nml_in, file=nml_filename, action='read', access='stream', &
            form='unformatted', iostat=nml_error)
       if (nml_error == 0) then
          read(unit=nml_in, iostat=nml_error, iomsg=ioerror_msg) nl_buffer

          ! we should always reach the EOF to capture the entire file...
          if (.not. is_iostat_end(nml_error)) then
             write(stdout, '(a, a, i8)') subname, &
                  ": IO ERROR reading namelist file into buffer: ", nml_error
             write(stdout, '(a)') ioerror_msg
          else
             write(stdout, '(a, a, a)') "Read '", trim(nml_filename), "' until EOF."
          end if

          write(stdout, '(a, a, i7, a)') subname, ": Read buffer of ", &
               len_trim(nl_buffer), " characters."

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

    ! broadcast the namelist buffer
    call broadcast_scalar(nl_buffer, master_task)

    ! now every process reads the namelist from the buffer
    ioerror_msg=''
    read(nl_buffer, nml=ecosys_driver_nml, iostat=nml_error, iomsg=ioerror_msg)
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
    call get_timer(ecosys_interior_timer, 'ECOSYS_INTERIOR', &
         nblocks_clinic, distrb_clinic%nprocs)

    !--------------------------------------------------------------------
    !
    !  Initialize marbl interface.
    !
    !  FIXME(bja, 2014-12) should this go here or in
    !  ecosys_driver_tracer_cnt_init...?  It should eventually be
    !  called from tracer_cnt_init to return marbl_sizes back to pop
    !  correctly? But we can't put it there now because we are still
    !  passing big arrays into ecosys_init through marbl_init.
    !
    !--------------------------------------------------------------------
    marbl_status%status = marbl_status_ok
    marbl_status%message = ''

    ! tell marbl how big the domain is
    marbl_driver_sizes%nx_block = nx_block
    marbl_driver_sizes%ny_block = ny_block
    marbl_driver_sizes%max_blocks_clinic = max_blocks_clinic
    marbl_driver_sizes%km = km

    ! NOTE(bja, 2015-08) we are intentionally allocating memory here
    ! but not initializing because marbl ecosys won't know anything
    ! about the domain. initialization of elements will occur in
    ! marbl/ecosys_init
    allocate(marbl_saved_state%dust_flux_in(nx_block, ny_block, max_blocks_clinic))
    allocate(marbl_saved_state%par_out(nx_block, ny_block, max_blocks_clinic))
    allocate(marbl_saved_state%ph_prev_3d(nx_block, ny_block, km, max_blocks_clinic))
    allocate(marbl_saved_state%ph_prev_alt_co2_3d(nx_block, ny_block, km, max_blocks_clinic))
    allocate(marbl_saved_state%land_mask(nx_block, ny_block, nblocks_clinic) )
    
    call marbl%init(marbl_driver_sizes, marbl_sizes, nl_buffer, marbl_status)
    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: marbl_init returned status: "'//marbl_status%message//'"')
    end if

    ! initialize ecosys_diagnostics type
    do bid=1,nblocks_clinic
      call marbl_diagnostics(bid)%construct(km, ecosys_diag_cnt_2d,           &
                      ecosys_diag_cnt_3d, auto_diag_cnt_2d, auto_diag_cnt_3d, &
                      zoo_diag_cnt_2d, zoo_diag_cnt_3d, part_diag_cnt_2d,     &
                      part_diag_cnt_3d, ecosys_tracer_cnt, autotroph_cnt,     &
                      zooplankton_cnt)
    end do


    ! now we know how many tracers marbl has, we can verify that pop
    ! has the correctly sized data.

    !-----------------------------------------------------------------------
    !  set up indices for ecosys modules that are on
    !  These indices are relative to ecosys_driver_ind_begin/end
    !  This means they start at 1 and end at ecosys_driver_tracer_cnt
    !  ecosys_driver then passes the ecosys module tracers back to passive
    !  tracers within TRACER(:,:,:,ecosys_driver_ind_beg,ecosys_driver_ind_end)
    !-----------------------------------------------------------------------

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

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    !  initialize virtual flux flag array
    vflux_flag(:) = .false.
    vflux_flag(dic_ind) = .true.
    vflux_flag(alk_ind) = .true.
    vflux_flag(dic_alt_co2_ind) = .true.

    !  allocate various  allocatable module variables
    allocate( PH_PREV(nx_block, ny_block, max_blocks_clinic) )
    allocate( PH_PREV_ALT_CO2(nx_block, ny_block, max_blocks_clinic) )

    tadvect_ctype(ecosys_ind_begin:ecosys_ind_end) = ecosys_tadvect_ctype

    ! initialize marbl namelists
    call  ecosys_init_nml(nl_buffer, comp_surf_avg_flag,          &
       use_nml_surf_vals, surf_avg_dic_const, surf_avg_alk_const, &
       init_ecosys_option, init_ecosys_init_file,                 &
       init_ecosys_init_file_fmt,tracer_init_ext, errorCode, marbl_status)

    ! initialize marbl land mask
    if (lmarginal_seas) then
       marbl_saved_state%land_mask = REGION_MASK /= 0
    else
       marbl_saved_state%land_mask = REGION_MASK > 0
    endif

    ! initialize metadata for ecosys tracers
    call ecosys_init_tracer_metadata(tracer_d_module)

    ! initialize ecosys tracers
    call ecosys_driver_init_tracers(&
       init_ts_file_fmt, &
       init_ecosys_option, init_ecosys_init_file, init_ecosys_init_file_fmt, &
       read_restart_filename, vflux_flag, use_nml_surf_vals, &
       tracer_init_ext, &
       tracer_d_module, &
       TRACER_MODULE, &
       ecosys_restart_filename, PH_PREV, PH_PREV_ALT_CO2, marbl_saved_state, &
       comp_surf_avg_flag, surf_avg, surf_avg_dic_const, surf_avg_alk_const, &
       errorCode)       

    ! initialize remaining variables
    call ecosys_init_postnml(init_ts_file_fmt,                        &
       read_restart_filename, tracer_d_module, TRACER_MODULE,         &
       lmarginal_seas, ecosys_restore, marbl_saved_state, vflux_flag, &
       use_nml_surf_vals, comp_surf_avg_flag, surf_avg,               &
       surf_avg_dic_const, surf_avg_alk_const, init_ecosys_option,    &
       init_ecosys_init_file, init_ecosys_init_file_fmt,              &
       tracer_init_ext,                                               &
       PH_PREV, PH_PREV_ALT_CO2, errorCode, marbl_status)
    
    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: ecosys_init returned status: "'//marbl_status%message//'"')
    end if
    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_init')
       return
    endif

    call ecosys_tavg_init(ecosys_restore)

    !-----------------------------------------------------------------------
    !  ECOSYS CISO block
    !-----------------------------------------------------------------------
    if (ciso_on) then
       tadvect_ctype(ecosys_ciso_ind_begin:ecosys_ciso_ind_end) = ecosys_tadvect_ctype
       call ecosys_ciso_init(init_ts_file_fmt, read_restart_filename,                       &
            tracer_d_module(ecosys_ciso_ind_begin:ecosys_ciso_ind_end),         &
            TRACER_MODULE(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:,:), &
            lmarginal_seas,                               &
            errorCode)

       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, &
               'init_ecosys_driver: error in ecosys_ciso_init')
          return
       endif

    end if


    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_init

  !-----------------------------------------------------------------------

  subroutine ecosys_driver_init_tracers(&
       init_ts_file_fmt, &
       init_ecosys_option, init_ecosys_init_file, init_ecosys_init_file_fmt, &
       read_restart_filename, vflux_flag, use_nml_surf_vals, &
       tracer_init_ext, &
       tracer_d_module, &
       TRACER_MODULE, &
       ecosys_restart_filename, PH_PREV, PH_PREV_ALT_CO2, saved_state, &
       comp_surf_avg_flag, surf_avg, surf_avg_dic_const, surf_avg_alk_const, &
       errorCode)       

    use marbl_interface_types , only : marbl_saved_state_type
    use passive_tracer_tools  , only : rest_read_tracer_block
    use passive_tracer_tools  , only : file_read_tracer_block
    use passive_tracer_tools  , only : field_exists_in_file
    use passive_tracer_tools  , only : read_field
    use passive_tracer_tools  , only : extract_surf_avg
    use passive_tracer_tools  , only : tracer_read
    use passive_tracer_tools  , only : ind_name_pair
    use passive_tracer_tools  , only : comp_surf_avg
    use prognostic            , only : tracer_field
    use prognostic            , only : curtime
    use prognostic            , only : oldtime
    use time_management       , only : check_time_flag
    use time_management       , only : eval_time_flag
    use ecosys_constants      , only : ecosys_tracer_cnt
    use io_tools              , only : document
    use exit_mod              , only : exit_POP
    use grid                  , only : fill_points
    use grid                  , only : n_topo_smooth
    
    implicit none
    
    character (*)                , intent(in)    :: init_ts_file_fmt                   ! format (bin or nc) for input file
    character(char_len)          , intent(in)    :: init_ecosys_option                 ! namelist option for initialization of bgc
    character(char_len)          , intent(in)    :: init_ecosys_init_file              ! filename for option 'file'
    character (*)                , intent(in)    :: read_restart_filename              ! file name for restart file
    logical (log_kind)           , intent(in)    :: vflux_flag(:) 
    logical (log_kind)           , intent(in)    :: use_nml_surf_vals                  ! do namelist surf values override values from restart    
    type(tracer_read)            , intent(in)    :: tracer_init_ext(ecosys_tracer_cnt) ! namelist variable for initializing tracers
    type(tracer_field)           , intent(in)    :: tracer_d_module(:)                 ! descriptors for each tracer

    character(char_len)          , intent(inout) :: init_ecosys_init_file_fmt          ! file format for option 'file'
    real (r8)                    , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)
    type(marbl_saved_state_type) , intent(inout) :: saved_state

    character(char_len)          , intent(out)   :: ecosys_restart_filename            ! modified file name for restart file
    real (r8)                    , intent(out)   :: PH_PREV(:, :, :)            ! computed ph from previous time step
    real (r8)                    , intent(out)   :: PH_PREV_ALT_CO2(:, :, :)    ! computed ph from previous time step
    integer (int_kind)           , intent(out)   :: comp_surf_avg_flag          ! time flag id for computing average surface tracer
    real (r8)                    , intent(out)   :: surf_avg(ecosys_tracer_cnt) ! average surface tracer values
    real (r8)                    , intent(out)   :: surf_avg_dic_const
    real (r8)                    , intent(out)   :: surf_avg_alk_const
    integer (POP_i4)             , intent(out)   :: errorCode
    
    
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_driver_mod:ecosys_driver_init_tracers'
    integer :: n, k
    type(ind_name_pair) :: ind_name_table(ecosys_tracer_cnt)

    !-----------------------------------------------------------------------
    !  initialize tracers
    !-----------------------------------------------------------------------
    
    !  initialize ind_name_table
    do n = 1, ecosys_tracer_cnt
       ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
    end do

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
               'PH_SURF', PH_PREV)
       else
          call document(subname, 'PH_SURF does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV to 0')
          PH_PREV = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_SURF_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_SURF_ALT_CO2', PH_PREV_ALT_CO2)
       else
          call document(subname, 'PH_SURF_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2 to 0')
          PH_PREV_ALT_CO2 = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_3D')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D', saved_state%PH_PREV_3D)
       else
          call document(subname, 'PH_3D does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_3D to 0')
          saved_state%PH_PREV_3D  = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_3D_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D_ALT_CO2', saved_state%PH_PREV_ALT_CO2_3D)
       else
          call document(subname, 'PH_3D_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2_3D to 0')
          saved_state%PH_PREV_ALT_CO2_3D = c0
       endif

       if (use_nml_surf_vals) then
          surf_avg(:) = c0
          surf_avg(dic_ind) = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind) = surf_avg_alk_const
       else
          call extract_surf_avg(init_ecosys_init_file_fmt,     &
               ecosys_restart_filename,       &
               ecosys_tracer_cnt, vflux_flag, &
               ind_name_table, surf_avg)
       endif

       call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

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

       PH_PREV = c0
       PH_PREV_ALT_CO2 = c0
       saved_state%PH_PREV_3D = c0
       saved_state%PH_PREV_ALT_CO2_3D = c0

       if (use_nml_surf_vals) then
          surf_avg = c0
          surf_avg(dic_ind) = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind) = surf_avg_alk_const
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

  end subroutine ecosys_driver_init_tracers

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_set_interior
  ! !INTERFACE:

  subroutine ecosys_driver_set_interior(ciso_on, TEMP_OLD, &
       TEMP_CUR, SALT_OLD, SALT_CUR, &
       TRACER_MODULE_OLD, TRACER_MODULE_CUR, DTRACER_MODULE, &
       this_block)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute source-sink terms
    !  accumulate commnon tavg fields related to source-sink terms
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_interface_types , only : marbl_column_domain_type
    use marbl_interface_types , only : marbl_gcm_state_type
    use constants             , only : salt_to_ppt
    use grid                  , only : KMT
    use grid                  , only : DZT
    use grid                  , only : dz
    use grid                  , only : partial_bottom_cells
    
    ! !INPUT PARAMETERS:

    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

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
    !EOP
    !BOC

    ! local
    integer (int_kind) :: i ! nx_block loop index
    integer (int_kind) :: c ! ny_block / column loop index
    integer (int_kind) :: k ! vertical level index
    integer (int_kind) :: n, d
    integer (int_kind) :: bid ! local block address for this block

    type(marbl_column_domain_type) :: marbl_domain
    type(marbl_gcm_state_type)     :: marbl_gcm_state

    real (r8), dimension(nx_block, ny_block, km, ecosys_tracer_cnt) :: tracer_module_avg
    real (r8), dimension(nx_block, ny_block, km) :: &
         temperature,          & ! temperature (C)
         salinity                ! salinity(ppt)

    real(r8), dimension(ecosys_tracer_cnt, km) :: column_tracer_module
    real(r8), dimension(ecosys_tracer_cnt, km) :: column_dtracer

    call marbl%set_interior()

    bid = this_block%local_id
    call marbl_diagnostics(bid)%set_to_zero()

    ! FIXME(bja, 2015-07) one time copy of global marbl_domain
    ! related memory from slab to column ordering. move entire
    ! copy to ecosys_driver_set_interior.
    marbl_domain%km = km
    allocate(marbl_domain%dzt(marbl_domain%km))
    allocate(marbl_domain%dz(marbl_domain%km))
    allocate(marbl_gcm_state%temperature(marbl_domain%km))
    allocate(marbl_gcm_state%salinity(marbl_domain%km))

    !-----------------------------------------------------------------------
    !  ECOSYS computations
    !-----------------------------------------------------------------------

    ! gcm dependent quantities (i.e. time stepping). need to be
    ! averaged into a single quantity for marbl ecosys
    temperature  = p5*(TEMP_OLD + TEMP_CUR)
    salinity     = p5*(SALT_OLD + SALT_CUR)*salt_to_ppt
    tracer_module_avg = p5*(&
         TRACER_MODULE_OLD(:, :, :, ecosys_ind_begin:ecosys_ind_end) + &
         TRACER_MODULE_CUR(:, :, :, ecosys_ind_begin:ecosys_ind_end))

    call timer_start(ecosys_interior_timer, block_id=bid)

    do c = this_block%jb,this_block%je
       do i = this_block%ib,this_block%ie
       ! Copy data from slab to column for marbl
          marbl_domain%land_mask = marbl_saved_state%land_mask(i, c, bid)
          marbl_domain%kmt = KMT(i, c, bid)
          if (marbl_saved_state%land_mask(i,c,bid)) then
             do k = 1, marbl_domain%km
                marbl_gcm_state%temperature(k) = temperature(i,c,k)
                marbl_gcm_state%salinity(k) = salinity(i,c,k)
                ! FIXME(bja, 2015-07) marbl shouldn't know about partial
                ! bottom cells. just pass in delta_z, set to DZT(i, c,
                ! k, bid) or dz(k) and use the single value!
                marbl_domain%dz(k) = dz(k)
                if (partial_bottom_cells) then
                   marbl_domain%dzt(k) = DZT(i, c, k, bid)
                else
                   marbl_domain%dzt(k) = c0
                end if
                column_tracer_module(:, k) = tracer_module_avg(i, c, k, ecosys_ind_begin:ecosys_ind_end)
             end do ! do k
          end if ! land_mask

          
          if (marbl_saved_state%land_mask(i,c,bid) .and. KMT(i, c, bid).gt.0) then 

            call ecosys_set_interior(i, c, bid,                             &
                 marbl_domain, marbl_gcm_state,                             &
                 marbl_diagnostics(bid), marbl_saved_state, ecosys_restore, &
                 marbl%private_data%ecosys_interior_share,                  &
                 marbl%private_data%ecosys_zooplankton_share,               &
                 marbl%private_data%ecosys_autotroph_share,                 &
                 marbl%private_data%ecosys_particulate_share,               &
                 column_tracer_module,                                      &
                 column_dtracer,                                            &
                 ciso_on)

            ! copy marbl column data back to slab
            do k = 1, marbl_domain%km
               do n = ecosys_ind_begin, ecosys_ind_end
                  DTRACER_MODULE(i, c, k, n) = column_dtracer(n, k)
               end do ! do n
            end do ! do k
          end if ! KMT > 0

          call ecosys_tavg_accumulate(i, c, bid, marbl_diagnostics(bid), ecosys_restore)

       end do ! do i
    end do ! do c
          

    call timer_stop(ecosys_interior_timer, block_id=bid)

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO computations
    !-----------------------------------------------------------------------

    if (ciso_on) then
       do k = 1, km
          call ecosys_ciso_set_interior(k,                                      &
               marbl%private_data%ecosys_interior_share(k),                               &
               marbl%private_data%ecosys_zooplankton_share(k),                            &
               marbl%private_data%ecosys_autotroph_share(k),                            &
               marbl%private_data%ecosys_particulate_share(k),                            &
               TEMP_OLD(:, :, k), TEMP_CUR(:, :, k),                                                &
               TRACER_MODULE_OLD(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),&
               TRACER_MODULE_CUR(:,:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),&
               DTRACER_MODULE(:,:,k,ecosys_ciso_ind_begin:ecosys_ciso_ind_end),     &
               bid)
       end do
    end if

    deallocate(marbl_domain%dzt)
    deallocate(marbl_domain%dz)
    deallocate(marbl_gcm_state%temperature)
    deallocate(marbl_gcm_state%salinity)

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_set_interior

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_set_sflux
  ! !INTERFACE:

  subroutine ecosys_driver_set_sflux(ciso_on,SHF_QSW_RAW, SHF_QSW, &
       U10_SQR,IFRAC,PRESS,SST,SSS, &
       SURFACE_VALS_OLD,SURFACE_VALS_CUR,STF_MODULE)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute surface fluxes
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:

    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
         SHF_QSW_RAW,  &! penetrative solar heat flux, from coupler (degC*cm/s)
         SHF_QSW,      &! SHF_QSW used by physics, may have diurnal cylce imposed (degC*cm/s)
         U10_SQR,      &! 10m wind speed squared (cm/s)**2
         IFRAC,        &! sea ice fraction (non-dimensional)
         PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
         SST,          &! sea surface temperature (C)
         SSS            ! sea surface salinity (psu)

    real (r8), dimension(:,:,:,:), &
         intent(in) :: SURFACE_VALS_OLD, SURFACE_VALS_CUR ! module tracers

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8), dimension(:,:,:,:), &
         intent(inout) :: STF_MODULE

    !EOP
    !BOC

    call marbl%set_surface_flux()

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    call ecosys_set_sflux(                                        &
         marbl_saved_state,                                       &
         marbl%private_data%surface_share,                        &
         SHF_QSW_RAW, SHF_QSW,                                    &
         U10_SQR, IFRAC, PRESS,                                   &
         SST, SSS,                                                &
         SURFACE_VALS_OLD(:,:,ecosys_ind_begin:ecosys_ind_end,:), &
         SURFACE_VALS_CUR(:,:,ecosys_ind_begin:ecosys_ind_end,:), &
         STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:),       &
         ciso_on,                                                 &
         vflux_flag, PH_PREV, PH_PREV_ALT_CO2, comp_surf_avg_flag, surf_avg) !new args

    !-----------------------------------------------------------------------
    !  ECOSYSC_CISO block
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ecosys_ciso_set_sflux(                                         &
            marbl%private_data%surface_share,                                &
            SST, &
            SURFACE_VALS_OLD(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:), &
            SURFACE_VALS_CUR(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:), &
            STF_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:))
    end if


    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_set_sflux

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_write_restart
  ! !INTERFACE:

  subroutine ecosys_driver_write_restart(ciso_on,restart_file, action)

    ! !DESCRIPTION:
    !  call restart routines for each tracer module that
    !  write fields besides the tracers themselves
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:
    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    character(*), intent(in) :: action

    ! !INPUT/OUTPUT PARAMETERS:

    type (datafile), intent (inout)  :: restart_file

    !EOP
    !BOC

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    call ecosys_write_restart(marbl_saved_state, restart_file, action, &
         vflux_flag, PH_PREV, PH_PREV_ALT_CO2)

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO block
    !-----------------------------------------------------------------------
    if (ciso_on) then
       call ecosys_ciso_write_restart(restart_file, action)
    end if

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_write_restart


  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_driver_tavg_forcing
  ! !INTERFACE:

  subroutine ecosys_driver_tavg_forcing(ciso_on,STF_MODULE)

    ! !DESCRIPTION:
    !  accumulate common tavg fields for tracer surface fluxes
    !  call accumation subroutines for tracer modules that have additional
    !     tavg fields related to surface fluxes
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:
    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    real (r8), dimension(:,:,:,:), &
         intent(in) :: STF_MODULE

    !EOP
    !BOC

    real (r8), dimension(nx_block,ny_block, forcing_diag_cnt, nblocks_clinic) :: &
         FLUX_DIAGS                ! Computed diagnostics for surface fluxes

    !-----------------------------------------------------------------------
    !  call routines from modules that have additional sflux tavg fields
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    call ecosys_tavg_forcing(marbl_saved_state, STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:),&
         FLUX_DIAGS)
    call ecosys_tavg_accumulate_flux(FLUX_DIAGS)

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO block
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ecosys_ciso_tavg_forcing( &
            STF_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:))
    end if


    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_driver_tavg_forcing


  !***********************************************************************

  function ecosys_driver_tracer_ref_val(ciso_on,ind)
    !
    ! !DESCRIPTION:
    !  return reference value for tracer with global tracer index ind
    !  this is used in virtual flux computations
    !
    !INPUT PARAMETERS:
    logical (kind=log_kind) , intent(in) :: ciso_on                 ! ecosys_ciso on
    integer(int_kind)       , intent(in) :: ind

    !OUTPUT PARAMETERS:
    real(r8) :: ecosys_driver_tracer_ref_val

    !-----------------------------------------------------------------------
    !  default value for reference value is 0
    !-----------------------------------------------------------------------

    ecosys_driver_tracer_ref_val = c0

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    if (ind >= ecosys_ind_begin .and. ind <= ecosys_ind_end) then
       ecosys_driver_tracer_ref_val = ecosys_tracer_ref_val(ind-ecosys_ind_begin+1, vflux_flag, surf_avg)
    endif

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO block
    !-----------------------------------------------------------------------

    if (ciso_on) then
       if (ind >= ecosys_ciso_ind_begin .and. ind <= ecosys_ciso_ind_end) then
          ecosys_driver_tracer_ref_val = ecosys_ciso_tracer_ref_val(ind-ecosys_ciso_ind_begin+1)
       endif
    endif

  end function ecosys_driver_tracer_ref_val

  !***********************************************************************

  subroutine ecosys_driver_unpack_source_sink_terms(source, destination)

    ! input parameters
    real(r8), dimension(:, :, :), intent(in) :: &
         source

    ! output parameters
    real(r8), dimension(:, :, :), intent(out) :: &
         destination

    destination(:, :, :) = source(:, :, :)

  end subroutine ecosys_driver_unpack_source_sink_terms

  !*****************************************************************************

  subroutine ecosys_write_restart(saved_state, restart_file, action, &
       vflux_flag, PH_PREV, PH_PREV_ALT_CO2) 

    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files

    use marbl_interface_types , only : marbl_saved_state_type
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

    ! !INPUT PARAMETERS:
    type(marbl_saved_state_type) , intent(in) :: saved_state
    character(*)                 , intent(in) :: action
    logical (log_kind)           , intent(in) :: vflux_flag(:) 
    real (r8)                    , intent(in) :: PH_PREV(:, :, :)         ! computed ph from previous time step
    real (r8)                    , intent(in) :: PH_PREV_ALT_CO2(:, :, :) ! computed ph from previous time step

    ! !INPUT/OUTPUT PARAMETERS:
    type (datafile), intent (inout)  :: restart_file

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character (char_len)       :: short_name   ! tracer name temporaries
    type (io_dim)              :: i_dim, j_dim ! dimension descriptors
    type (io_dim)              :: k_dim        ! dimension descriptor for vertical levels
    integer (int_kind)         :: n
    type (io_field_desc), save :: PH_SURF
    type (io_field_desc), save :: PH_SURF_ALT_CO2
    type (io_field_desc), save :: PH_3D_ALT_CO2
    type (io_field_desc), save :: PH_3D

    !-----------------------------------------------------------------------

    if (trim(action) == 'add_attrib_file') then
       short_name = char_blank
       do n=1, ecosys_tracer_cnt
          if (vflux_flag(n)) then   
             short_name = 'surf_avg_' /&
                  &/ ind_name_table(n)%name
             call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
          endif
       end do
    endif

    if (trim(action) == 'define') then
       i_dim = construct_io_dim('i', nx_global)
       j_dim = construct_io_dim('j', ny_global)
       k_dim = construct_io_dim('k', km)

       PH_SURF = construct_io_field('PH_SURF', i_dim, j_dim,     &
            long_name='surface pH at current time',      &
            units='pH', grid_loc='2110',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d2d_array = PH_PREV(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_SURF)

       PH_SURF_ALT_CO2 = construct_io_field('PH_SURF_ALT_CO2', i_dim, j_dim, &
            long_name='surface pH, alternate CO2, at current time', &
            units='pH', grid_loc='2110',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d2d_array = PH_PREV_ALT_CO2(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_SURF_ALT_CO2)

       PH_3D_ALT_CO2 = construct_io_field('PH_3D_ALT_CO2', i_dim, j_dim, k_dim, &
            long_name='3D pH, alternate CO2, at current time', &
            units='pH', grid_loc='3111',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d3d_array = saved_state%PH_PREV_ALT_CO2_3D(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_3D_ALT_CO2)

       PH_3D = construct_io_field('PH_3D', i_dim, j_dim, k_dim, &
            long_name='3D pH at current time', &
            units='pH', grid_loc='3111',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d3d_array = saved_state%PH_PREV_3D(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_3D)

    endif

    if (trim(action) == 'write') then
       call data_set (restart_file, 'write', PH_SURF)
       call data_set (restart_file, 'write', PH_SURF_ALT_CO2)
       call data_set (restart_file, 'write', PH_3D)
       call data_set (restart_file, 'write', PH_3D_ALT_CO2)
    endif

  end subroutine ecosys_write_restart

  !*****************************************************************************

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

