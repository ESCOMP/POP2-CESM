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
  use prognostic                , only : TRACER
  use io_types                  , only : stdout, nml_in, nml_filename, datafile
  use exit_mod                  , only : sigAbort, exit_pop
  use constants                 , only : c0, c1, p5, delim_fmt, char_blank, ndelim_fmt
  use passive_tracer_tools      , only : ind_name_pair

  use marbl_interface           , only : marbl_interface_class
  use marbl_interface           , only : marbl_sizes_type
  use marbl_interface           , only : marbl_driver_sizes_type

  use marbl_interface_constants , only : marbl_status_ok
  use marbl_interface_types     , only : marbl_status_type
  use marbl_interface_types     , only : marbl_diagnostics_type
  use marbl_interface_types     , only : marbl_saved_state_type

  ! NOTE(bja, 2014-12) all other uses of marbl/ecosys modules need to be removed!
  use marbl_share_mod           , only : autotroph_cnt, zooplankton_cnt

  use ecosys_constants, only : ecosys_tracer_cnt

  use ecosys_tavg, only : ecosys_tavg_init
  use ecosys_tavg, only : ecosys_tavg_accumulate
  use ecosys_tavg, only : ecosys_tavg_accumulate_flux

  use ecosys_mod, only: ecosys_init_nml    
  use ecosys_mod, only: ecosys_init_tracer_metadata
  use ecosys_mod, only: ecosys_init_postnml ! TEMPORARY
  use ecosys_mod, only: ecosys_init_tavg    ! TEMPORARY
  use ecosys_mod, only: ecosys_tracer_ref_val
  use ecosys_mod, only: ecosys_set_sflux
  use ecosys_mod, only: ecosys_tavg_forcing
  use ecosys_mod, only: ecosys_set_interior

  use ecosys_ciso_mod, only: ecosys_ciso_tracer_cnt
  use ecosys_ciso_mod, only: ecosys_ciso_init
  use ecosys_ciso_mod, only: ecosys_ciso_tracer_ref_val
  use ecosys_ciso_mod, only: ecosys_ciso_set_sflux
  use ecosys_ciso_mod, only: ecosys_ciso_tavg_forcing
  use ecosys_ciso_mod, only: ecosys_ciso_set_interior
  use ecosys_ciso_mod, only: ecosys_ciso_write_restart

  use ecosys_restore_mod, only : ecosys_restore_type

  use timers, only : timer_start
  use timers, only : timer_stop
  use timers, only : get_timer

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_driver_tracer_cnt_init
  public :: ecosys_driver_init
  public :: ecosys_driver_set_interior
  public :: ecosys_driver_set_sflux
  public :: ecosys_driver_tracer_cnt
  public :: ecosys_driver_tracer_ref_val
  public :: ecosys_driver_tavg_forcing
  public :: ecosys_driver_write_restart
  public :: ecosys_driver_unpack_source_sink_terms

  private :: ecosys_write_restart
  private :: ecosys_driver_init_tracers
  private :: ecosys_driver_init_interior_restore 

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
  type(marbl_interface_class)   :: marbl
  type(marbl_sizes_type)        :: marbl_sizes
  type(marbl_driver_sizes_type) :: marbl_driver_sizes
  type(marbl_status_type)       :: marbl_status
  type(marbl_diagnostics_type)  :: marbl_diagnostics(max_blocks_clinic)
  ! FIXME(bja, 2015-08) this needs to go into marbl%private_data%saved_state !!!
  type(marbl_saved_state_type)  :: marbl_saved_state
  type(ecosys_restore_type)     :: ecosys_restore

  !-----------------------------------------------------------------------
  !  average surface tracer value related variables
  !  used as reference value for virtual flux computations
  !-----------------------------------------------------------------------

  ! NOTE (mvertens, 2015-10) the following has been moved from ecosys_mod module variables
  ! to ecosys_driver module variables as part of the marbilization of the ecosys initialization

  logical (log_kind) :: vflux_flag(ecosys_tracer_cnt) ! which tracers get virtual fluxes applied
  real (r8)          :: surf_avg(ecosys_tracer_cnt)   ! average surface tracer values

  type(ind_name_pair) :: ind_name_table(ecosys_tracer_cnt) !  derived type & parameter for tracer index lookup

  !-----------------------------------------------------------------------
  !  named field indices
  !-----------------------------------------------------------------------

!  integer (int_kind) :: totChl_surf_nf_ind = 0 ! total chlorophyll in surface layer (TEMPORARY)
!  integer (int_kind) :: sflux_co2_nf_ind   = 0 ! air-sea co2 gas flux (TEMPORARY)
  integer (int_kind) :: atm_co2_nf_ind     = 0 ! atmospheric co2

  !-----------------------------------------------------------------------
  !  module variables related to ph computations
  !-----------------------------------------------------------------------

  real (r8), allocatable, target :: PH_PREV(:, :, :)         ! computed ph from previous time step
  real (r8), allocatable, target :: PH_PREV_ALT_CO2(:, :, :) ! computed ph from previous time step, alternative CO2
  real (r8), allocatable, target :: IRON_PATCH_FLUX(:, :, :) ! localized iron patch flux

  !-----------------------------------------------------------------------
  !  restoring climatologies for nutrients
  !-----------------------------------------------------------------------

  real (r8), allocatable, target :: FESEDFLUX(:, :, :, :)      !  sedimentary Fe inputs

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
    !-----------------------------------------------------------------------

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

  subroutine ecosys_driver_init(ciso_on,init_ts_file_fmt,     &
       read_restart_filename,                         &
       tracer_d_module, TRACER_MODULE, tadvect_ctype, &
       errorCode)


    ! !DESCRIPTION:
    !  Initialize ecosys_driver passive tracers. This involves:
    !  1) setting ecosys and ecosys_ciso module index bounds
    !  2) calling ecosys and ecosys_ciso module init subroutine

    use marbl_interface_constants , only : marbl_nl_buffer_size
    use grid                      , only : n_topo_smooth
    use grid                      , only : fill_points
    use grid                      , only : REGION_MASK
    use grid                      , only : KMT
    use broadcast                 , only : broadcast_scalar
    use prognostic                , only : tracer_field
    use io_tools                  , only : document
    use time_management           , only : eval_time_flag
    use passive_tracer_tools      , only : set_tracer_indices
    use passive_tracer_tools      , only : extract_surf_avg
    use passive_tracer_tools      , only : comp_surf_avg
    use passive_tracer_tools      , only : rest_read_tracer_block
    use passive_tracer_tools      , only : field_exists_in_file
    use passive_tracer_tools      , only : tracer_read
    use passive_tracer_tools      , only : file_read_tracer_block
    use passive_tracer_tools      , only : read_field
    use ecosys_diagnostics_mod    , only : ecosys_diag_cnt_2d
    use ecosys_diagnostics_mod    , only : ecosys_diag_cnt_3d
    use ecosys_diagnostics_mod    , only : auto_diag_cnt_2d
    use ecosys_diagnostics_mod    , only : auto_diag_cnt_3d
    use ecosys_diagnostics_mod    , only : zoo_diag_cnt_2d
    use ecosys_diagnostics_mod    , only : zoo_diag_cnt_3d
    use ecosys_diagnostics_mod    , only : part_diag_cnt_2d
    use ecosys_diagnostics_mod    , only : part_diag_cnt_3d
    use ecosys_diagnostics_mod    , only : forcing_diag_cnt
    use communicate               , only : my_task, master_task
    use named_field_mod           , only : named_field_register
    use named_field_mod           , only : named_field_set
    use marbl_share_mod           , only : totChl_surf_nf_ind !TEMPORARY
    use marbl_share_mod           , only : sflux_co2_nf_ind   !TEMPORARY
    use marbl_share_mod           , only : autotrophs
    use prognostic                , only : curtime
    use prognostic                , only : oldtime

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use marbl_share_mod           , only : comp_surf_avg_flag

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

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter         :: subname = 'ecosys_driver:ecosys_driver_init'
    integer (int_kind)              :: cumulative_nt, n, bid, k
    integer (int_kind)              :: nml_error                ! error flag for nml read
    integer (int_kind)              :: iostat                   ! io status flag
    character (char_len)            :: sname, lname, units, coordinates
    character (4)                   :: grid_loc
    character(marbl_nl_buffer_size) :: nl_buffer
    character(char_len_long)        :: ioerror_msg
    integer (int_kind)              :: auto_ind                 ! autotroph functional group index
    integer (int_kind)              :: iblock                   ! index for looping over blocks
    real(r8)                        :: WORK(nx_block, ny_block) ! FIXME (mvertens, 2015-10) remove this

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

    tadvect_ctype(ecosys_ind_begin:ecosys_ind_end) = ecosys_tadvect_ctype

    ! initialize marbl namelists
    call  ecosys_init_nml(nl_buffer, errorCode, marbl_status)

    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: ecosys_init returned status: "'//marbl_status%message//'"')
    end if
    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_init')
       return
    endif

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
       init_ts_file_fmt, read_restart_filename, tracer_d_module, TRACER_MODULE, &
       ecosys_restart_filename, marbl_saved_state, surf_avg, errorCode)       

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_init')
       return
    endif

    ! initialize remaining variables
    call ecosys_init_postnml()                  !FIXME (mvertens, 2015-11) move routine to driver 
    call ecosys_init_tavg()                     !FIXME (mvertens, 2015-11) move routine to driver
    call ecosys_init_sflux(marbl_saved_state)   !FIXME (mvertens, 2015-11) move routine to driver

    ! initialize interior tavg restoring
    do n = 1, ecosys_tracer_cnt
       ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
    end do
    call ecosys_restore%init(nml_filename, nml_in, ind_name_table)
    call ecosys_driver_init_interior_restore(marbl_saved_state, ecosys_restore)
    call ecosys_tavg_init(ecosys_restore)

    !$OMP PARALLEL DO PRIVATE(iblock, n, k)
    do iblock=1, nblocks_clinic
       do n = 1, ecosys_tracer_cnt
          do k = 1, km
             where (.not. marbl_saved_state%land_mask(:, :, iblock) .or. k > KMT(:, :, iblock))
                TRACER_MODULE(:, :, k, n, curtime, iblock) = c0
                TRACER_MODULE(:, :, k, n, oldtime, iblock) = c0
             end where
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !  register Chl field for short-wave absorption
    !  apply land mask to tracers
    !  set Chl field for short-wave absorption

    call named_field_register('model_chlorophyll', totChl_surf_nf_ind)
    !$OMP PARALLEL DO PRIVATE(iblock, n, WORK)
    do iblock=1, nblocks_clinic
       WORK = c0
       do auto_ind = 1, autotroph_cnt
          n = autotrophs(auto_ind)%Chl_ind
          WORK = WORK + max(c0, p5*(TRACER_MODULE(:, :, 1, n, oldtime, iblock) + &
                                    TRACER_MODULE(:, :, 1, n, curtime, iblock)))
       end do
       call named_field_set(totChl_surf_nf_ind, iblock, WORK)
    enddo
    !$OMP END PARALLEL DO

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

  end subroutine ecosys_driver_init

  !-----------------------------------------------------------------------

  subroutine ecosys_driver_init_tracers(&
       init_ts_file_fmt, read_restart_filename, tracer_d_module, TRACER_MODULE, &
       ecosys_restart_filename, marbl_saved_state, surf_avg, errorCode)       

    use marbl_interface_types , only : marbl_saved_state_type
    use marbl_parms           , only : dic_ind
    use marbl_parms           , only : alk_ind
    use marbl_parms           , only : dic_alt_co2_ind
    use marbl_share_mod       , only : init_ecosys_option
    use marbl_share_mod       , only : init_ecosys_init_file
    use marbl_share_mod       , only : init_ecosys_init_file_fmt
    use marbl_share_mod       , only : use_nml_surf_vals         
    use marbl_share_mod       , only : surf_avg_dic_const
    use marbl_share_mod       , only : surf_avg_alk_const
    use marbl_share_mod       , only : tracer_init_ext
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
    use prognostic            , only : newtime
    use time_management       , only : check_time_flag
    use time_management       , only : eval_time_flag
    use ecosys_constants      , only : ecosys_tracer_cnt
    use io_tools              , only : document
    use grid                  , only : fill_points
    use grid                  , only : n_topo_smooth
    use exit_mod              , only : exit_POP

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use marbl_share_mod           , only : comp_surf_avg_flag

    implicit none
    
    ! !INPUT PARAMETERS:

    character (*)                , intent(in)    :: init_ts_file_fmt                   ! format (bin or nc) for input file
    character (*)                , intent(in)    :: read_restart_filename              ! file name for restart file
    type(tracer_field)           , intent(in)    :: tracer_d_module(:)                 ! descriptors for each tracer

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8)                    , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)
    type(marbl_saved_state_type) , intent(inout) :: marbl_saved_state

    ! !OUTPUT PARAMETERS:

    character(char_len)          , intent(out)   :: ecosys_restart_filename            ! modified file name for restart file
    real (r8)                    , intent(out)   :: surf_avg(ecosys_tracer_cnt) ! average surface tracer values
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
    
    ! allocate module variables PH_PREV an PH_PREV_ALT_CO2
    allocate( PH_PREV(nx_block, ny_block, max_blocks_clinic) )
    allocate( PH_PREV_ALT_CO2(nx_block, ny_block, max_blocks_clinic) )

    ! initialize module variable - virtual flux flag array
    vflux_flag(:) = .false.
    vflux_flag(dic_ind) = .true.
    vflux_flag(alk_ind) = .true.
    vflux_flag(dic_alt_co2_ind) = .true.

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
               'PH_3D', marbl_saved_state%PH_PREV_3D)
       else
          call document(subname, 'PH_3D does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_3D to 0')
          marbl_saved_state%PH_PREV_3D  = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, ecosys_restart_filename, &
            'PH_3D_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D_ALT_CO2', marbl_saved_state%PH_PREV_ALT_CO2_3D)
       else
          call document(subname, 'PH_3D_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2_3D to 0')
          marbl_saved_state%PH_PREV_ALT_CO2_3D = c0
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
       marbl_saved_state%PH_PREV_3D = c0
       marbl_saved_state%PH_PREV_ALT_CO2_3D = c0

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

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: i ! nx_block loop index
    integer (int_kind) :: c ! ny_block / column loop index
    integer (int_kind) :: k ! vertical level index
    integer (int_kind) :: n, d
    integer (int_kind) :: bid ! local block address for this block

    type(marbl_column_domain_type) :: marbl_domain
    type(marbl_gcm_state_type)     :: marbl_gcm_state

    real (r8), dimension(nx_block, ny_block, km, ecosys_tracer_cnt) :: tracer_module_avg
    real (r8), dimension(nx_block, ny_block, km)                    :: temperature ! temperature (C)
    real (r8), dimension(nx_block, ny_block, km)                    :: salinity    ! salinity(ppt)

    real(r8), dimension(ecosys_tracer_cnt, km) :: column_tracer_module
    real(r8), dimension(ecosys_tracer_cnt, km) :: column_dtracer

    ! FIXME(bja, 2014-10) size of restore_local should be
    ! (non_living_biomass_ecosys_tracer_cnt, km) but non-living-biomass-tracer-count isn't global and I'm reluctant
    ! to make it right now.
    real(r8) :: tracer_local(ecosys_tracer_cnt)
    real(r8) :: restore_local(ecosys_tracer_cnt, km) ! local restoring terms for nutrients (mmol ./m^3/sec)

    !-----------------------------------------------------------------------

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

          if (marbl_domain%land_mask .and. marbl_domain%kmt > 0) then 

             !  set tracer restore fields
             do k = 1, marbl_domain%km
                do n = 1, ecosys_tracer_cnt
                   tracer_local(n) = max(c0, column_tracer_module(n,k))
                end do

                call ecosys_restore%restore_tracers(ecosys_tracer_cnt, &
                     vert_level=k, x_index=i, y_index=c, block_id=bid, &
                     local_data=tracer_local(:), restore_data=restore_local(:, k))
             end do

             !  compute time derivatives for ecosystem state variables
             call ecosys_set_interior(i, c, bid,                             &
                  marbl_domain, marbl_gcm_state,                             &
                  marbl_diagnostics(bid), marbl_saved_state, restore_local,  &
                  marbl%private_data%ecosys_interior_share,                  &
                  marbl%private_data%ecosys_zooplankton_share,               &
                  marbl%private_data%ecosys_autotroph_share,                 &
                  marbl%private_data%ecosys_particulate_share,               &
                  column_tracer_module,                                      &
                  column_dtracer,                                            &
                  fesedflux(i,c,:,bid), &
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

  subroutine ecosys_driver_set_sflux(ciso_on,SHF_QSW_RAW, SHF_QSW, &
       U10_SQR,IFRAC,PRESS,SST,SSS, &
       SURFACE_VALS_OLD,SURFACE_VALS_CUR,STF_MODULE)

    ! !DESCRIPTION:
    !  call subroutines for each tracer module that compute surface fluxes

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use marbl_share_mod           , only : comp_surf_avg_flag

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

    !-----------------------------------------------------------------------

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
         vflux_flag, PH_PREV, PH_PREV_ALT_CO2, comp_surf_avg_flag, surf_avg, & 
         IRON_PATCH_FLUX)

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

  end subroutine ecosys_driver_set_sflux

  !***********************************************************************

  subroutine ecosys_driver_write_restart(ciso_on,restart_file, action)

    ! !DESCRIPTION:
    !  call restart routines for each tracer module that
    !  write fields besides the tracers themselves

    ! !INPUT PARAMETERS:

    logical (kind=log_kind) , intent(in) :: ciso_on  ! ecosys_ciso on
    character(*)            , intent(in) :: action

    ! !INPUT/OUTPUT PARAMETERS:

    type (datafile), intent (inout)  :: restart_file

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

  end subroutine ecosys_driver_write_restart


  !***********************************************************************

  subroutine ecosys_driver_tavg_forcing(ciso_on,STF_MODULE)

    ! !DESCRIPTION:
    !  accumulate common tavg fields for tracer surface fluxes
    !  call accumulation subroutines for tracer modules that have additional
    !     tavg fields related to surface fluxes

    use ecosys_diagnostics_mod , only : forcing_diag_cnt
    use domain                 , only : nblocks_clinic
    use blocks                 , only : nx_block, ny_block

    implicit none 

    ! !INPUT PARAMETERS:
    logical (kind=log_kind), intent(in) ::  ciso_on                 ! ecosys_ciso on
    real (r8), intent(in)               :: STF_MODULE(:,:,:,:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    real (r8) :: FLUX_DIAGS(nx_block,ny_block, forcing_diag_cnt, nblocks_clinic)  ! Computed diagnostics for surface fluxes

    !-----------------------------------------------------------------------
    !  call routines from modules that have additional sflux tavg fields
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    call ecosys_tavg_forcing(marbl_saved_state, STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:), FLUX_DIAGS)

    call ecosys_tavg_accumulate_flux(FLUX_DIAGS)

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO block
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ecosys_ciso_tavg_forcing(STF_MODULE(:,:,ecosys_ciso_ind_begin:ecosys_ciso_ind_end,:))
    end if

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
             short_name = 'surf_avg_' // ind_name_table(n)%name
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

  subroutine ecosys_driver_init_interior_restore(marbl_saved_state, ecosys_restore)

    ! !DESCRIPTION:
    !  Initialize interior restoring computations for ecosys tracer module.

    use ecosys_restore_mod    , only : ecosys_restore_type
    use marbl_interface_types , only : marbl_saved_state_type
    use grid                  , only : KMT
    use grid                  , only : zt
    use io_types              , only : nml_in, nml_filename
    use blocks                , only : nx_block, ny_block
    use domain_size           , only : max_blocks_clinic, km
    use marbl_share_mod       , only : fesedflux_input
    use passive_tracer_tools  , only : read_field

    ! !INPUT PARAMETERS:

    type(marbl_saved_state_type), intent(in) :: marbl_saved_state

    ! !INPUT/OUTPUT PARAMETERS:

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
    !  initialize restoring timescale (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%initialize_restoring_timescale(nml_filename, nml_in, zt)

    !-----------------------------------------------------------------------
    !  load restoring fields (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%read_restoring_fields(marbl_saved_state%land_mask)

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
          where (.not. marbl_saved_state%land_mask(:, :, iblock) .or. k > KMT(:, :, iblock)) &
               fesedflux(:, :, k, iblock) = c0
          fesedflux(:, :, k, iblock) = fesedflux(:, :, k, iblock) * fesedflux_input%scale_factor
       enddo
    end do

  end subroutine ecosys_driver_init_interior_restore

  !***********************************************************************

  subroutine ecosys_init_sflux(saved_state)

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
    use marbl_interface_types , only : marbl_saved_state_type
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

    use marbl_share_mod       , only : sflux_co2_nf_ind   !TEMPORARY

    ! !INPUT/OUTPUT PARAMETERS:

    type(marbl_saved_state_type), intent(inout) :: saved_state

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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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
             where (.not. saved_state%land_mask(:, :, iblock)) &
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

  end subroutine ecosys_init_sflux

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

