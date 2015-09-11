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
  use broadcast                 , only : broadcast_scalar

  use marbl_interface           , only : marbl_interface_class
  use marbl_interface           , only : marbl_sizes_type
  use marbl_interface           , only : marbl_driver_sizes_type

  use marbl_interface_constants , only : marbl_status_ok
  use marbl_interface_types     , only : marbl_status_type

  use marbl_interface_types     , only : marbl_saved_state_type

  ! NOTE(bja, 2014-12) all other uses of marbl/ecosys modules need to be removed!
  use marbl_share_mod           , only : autotroph_cnt, zooplankton_cnt

  use ecosys_constants          , only : ecosys_tracer_cnt
  use ecosys_constants          , only : ecosys_diag_cnt
  use ecosys_constants          , only : auto_diag_cnt
  use ecosys_constants          , only : zoo_diag_cnt
  use ecosys_constants          , only : part_diag_cnt
  use ecosys_constants          , only : forcing_diag_cnt

  use ecosys_tavg               , only : ecosys_tavg_init
  use ecosys_tavg               , only : ecosys_tavg_write
  use ecosys_tavg               , only : ecosys_tavg_write_flux

  use ecosys_mod, only:            &
       ecosys_init,                &
       ecosys_tracer_ref_val,      &
       ecosys_set_sflux,           &
       ecosys_tavg_forcing,        &
       ecosys_set_interior,        &
       ecosys_write_restart


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
  ! FIXME(bja, 2015-08) this needs to go into marbl%private_data%saved_state !!!
  type(marbl_saved_state_type) :: marbl_saved_state


  type(ecosys_restore_type) :: ecosys_restore

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
    use marbl_interface_constants, only : marbl_nl_buffer_size

    ! !INPUT PARAMETERS:

    character (*), intent(in) :: &
         init_ts_file_fmt,    & ! format (bin or nc) for input file
         read_restart_filename  ! file name for restart file


    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    ! !INPUT/OUTPUT PARAMETERS:

    type (tracer_field), dimension(:), intent(inout) :: &
         tracer_d_module   ! descriptors for each tracer

    real (r8), dimension(:,:,:,:,:,:), &
         intent(inout) :: TRACER_MODULE

    ! !OUTPUT PARAMETERS:

    character (char_len), dimension(:), intent(out) :: &
         tadvect_ctype     ! advection method for ecosys tracers

    integer (POP_i4), intent(out) :: &
         errorCode

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_driver:ecosys_driver_init'

    integer (int_kind) :: cumulative_nt, n, &
         nml_error,        &! error flag for nml read
         iostat             ! io status flag

    character (char_len) :: sname, lname, units, coordinates
    character (4) :: grid_loc
    character(marbl_nl_buffer_size) :: nl_buffer
    character(char_len_long) :: ioerror_msg


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
    
    call marbl%init(marbl_driver_sizes, marbl_sizes, &
         nl_buffer, marbl_status)
    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: marbl_init returned status: "'//marbl_status%message//'"')
    end if

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

    call ecosys_init(nl_buffer, init_ts_file_fmt, read_restart_filename,       &
         tracer_d_module(ecosys_ind_begin:ecosys_ind_end),         &
         TRACER_MODULE(:,:,:,ecosys_ind_begin:ecosys_ind_end,:,:), &
         lmarginal_seas, ecosys_restore, marbl_saved_state, &
         errorCode, marbl_status)
    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in ecosys_driver_init: ecosys_init returned status: "'//marbl_status%message//'"')
    end if

    call ecosys_tavg_init(ecosys_restore)

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'init_ecosys_driver: error in ecosys_init')
       return
    endif


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

    use marbl_interface_types, only: ecosys_diagnostics_type
    use marbl_interface_types, only: marbl_diagnostics_type
    use marbl_interface_types, only : marbl_column_domain_type

    use constants, only : salt_to_ppt

    use grid, only : KMT
    use grid, only : DZT
    use grid, only : dz
    use grid, only : partial_bottom_cells
    


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

    real (r8), dimension(:,:,:,:), intent(out) :: &
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
    type(ecosys_diagnostics_type) :: ecosys_diagnostics(km)
    type(marbl_diagnostics_type) :: marbl_diagnostics(km)

    type(marbl_column_domain_type) :: marbl_domain

    real (r8), dimension(nx_block, ny_block, km, ecosys_tracer_cnt) :: tracer_module_avg
    real (r8), dimension(nx_block, ny_block, km) :: &
         temperature,          & ! temperature (C)
         salinity                ! salinity(ppt)

    real(r8), dimension(ecosys_tracer_cnt, km) :: column_tracer_module
    real(r8), dimension(ecosys_tracer_cnt, km) :: column_dtracer
    ! NOTE(bja MNL, 2015-03) eventually marbl%set_interior will be
    ! called with marbl_diagnostics_type and we may convert to
    ! ecosys_diagnostics_type for tavg. But for now we just use
    ! ecosys_diagnostics_type == marbl_diagnostics_type; I'm keeping
    ! the marbl_diagnostics name because eventually this will be
    ! calling marbl%set_interior.

    call marbl%set_interior()

    do k = 1, km
       allocate(ecosys_diagnostics(k)%DIAGS(nx_block, ny_block, ecosys_diag_cnt))
       allocate(ecosys_diagnostics(k)%AUTO_DIAGS(nx_block, ny_block, autotroph_cnt, auto_diag_cnt))
       allocate(ecosys_diagnostics(k)%ZOO_DIAGS(nx_block, ny_block, zooplankton_cnt, zoo_diag_cnt))
       allocate(ecosys_diagnostics(k)%PART_DIAGS(nx_block, ny_block, part_diag_cnt))
       allocate(ecosys_diagnostics(k)%restore_diags(nx_block, ny_block, ecosys_tracer_cnt))

       ecosys_diagnostics(k)%DIAGS = c0
       ecosys_diagnostics(k)%AUTO_DIAGS = c0
       ecosys_diagnostics(k)%ZOO_DIAGS = c0
       ecosys_diagnostics(k)%PART_DIAGS = c0
       ecosys_diagnostics(k)%restore_diags = c0

       allocate(marbl_diagnostics(k)%DIAGS(ecosys_diag_cnt))
       allocate(marbl_diagnostics(k)%AUTO_DIAGS(auto_diag_cnt, autotroph_cnt))
       allocate(marbl_diagnostics(k)%ZOO_DIAGS(zoo_diag_cnt, zooplankton_cnt))
       allocate(marbl_diagnostics(k)%PART_DIAGS(part_diag_cnt))
       allocate(marbl_diagnostics(k)%restore_diags(ecosys_tracer_cnt))

    end do

    ! FIXME(bja, 2015-07) one time copy of global marbl_domain
    ! related memory from slab to column ordering. move entire
    ! copy to ecosys_driver_set_interior.
    marbl_domain%km = km
    allocate(marbl_domain%dzt(marbl_domain%km))
    allocate(marbl_domain%dz(marbl_domain%km))
    allocate(marbl_domain%temperature(marbl_domain%km))
    allocate(marbl_domain%salinity(marbl_domain%km))

    bid = this_block%local_id
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
                marbl_domain%temperature(k) = temperature(i,c,k)
                marbl_domain%salinity(k) = salinity(i,c,k)
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

          
          if (marbl_saved_state%land_mask(i,c,bid) .and. &
              KMT(i, c, bid).gt.0) then 
            call ecosys_set_interior(i, c, ny_block, marbl_domain,            &
                 marbl_diagnostics, marbl_saved_state, ecosys_restore,        &
                 marbl%private_data%ecosys_interior_share,                    &
                 marbl%private_data%ecosys_zooplankton_share,                 &
                 marbl%private_data%ecosys_autotroph_share,                   &
                 marbl%private_data%ecosys_particulate_share,                 &
                 column_tracer_module,                                        &
                 column_dtracer,                                              &
                 ciso_on,                                                     &
                 bid)

            ! copy marbl column data back to slab
            do k = 1, marbl_domain%km
               do n = ecosys_ind_begin, ecosys_ind_end
                  DTRACER_MODULE(i, c, k, n) = column_dtracer(n, k)
               end do ! do n

               do d = 1, ecosys_diag_cnt
                  ecosys_diagnostics(k)%diags(i, c, d) = marbl_diagnostics(k)%diags(d)
               end do
          
               do n = 1, autotroph_cnt
                  do d = 1, auto_diag_cnt
                     ecosys_diagnostics(k)%auto_diags(i, c, n, d) = marbl_diagnostics(k)%auto_diags(d, n)
                  end do ! do d
               end do ! do n

               do n = 1, zooplankton_cnt
                  do d = 1, zoo_diag_cnt
                     ecosys_diagnostics(k)%zoo_diags(i, c, n, d) = marbl_diagnostics(k)%zoo_diags(d, n)
                  end do ! do d
               end do ! do n

               do d = 1, part_diag_cnt
                  ecosys_diagnostics(k)%part_diags(i, c, d) = marbl_diagnostics(k)%part_diags(d)
               end do ! do d

               do d = 1, ecosys_tracer_cnt
                  ecosys_diagnostics(k)%restore_diags(i, c, d) = marbl_diagnostics(k)%restore_diags(d)
               end do ! do d
            end do ! do k
          end if ! KMT > 0

       end do ! do i
    end do ! do c
          

    call timer_stop(ecosys_interior_timer, block_id=bid)

    do k = 1, km
       call ecosys_tavg_write(k, bid, ecosys_diagnostics(k))
       do n = 1, ecosys_tracer_cnt
          call ecosys_restore%accumulate_tavg(tracer_index=n, &
               vert_level=k, block_id=bid, &
               restore_local=ecosys_diagnostics(k)%restore_diags(:, :, n))
       end do
    end do
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

    do k = 1, km
       deallocate(ecosys_diagnostics(k)%DIAGS)
       deallocate(ecosys_diagnostics(k)%AUTO_DIAGS)
       deallocate(ecosys_diagnostics(k)%ZOO_DIAGS)
       deallocate(ecosys_diagnostics(k)%PART_DIAGS)
       deallocate(ecosys_diagnostics(k)%restore_diags)

       deallocate(marbl_diagnostics(k)%DIAGS)
       deallocate(marbl_diagnostics(k)%AUTO_DIAGS)
       deallocate(marbl_diagnostics(k)%ZOO_DIAGS)
       deallocate(marbl_diagnostics(k)%PART_DIAGS)
       deallocate(marbl_diagnostics(k)%restore_diags)
    end do

    deallocate(marbl_domain%dzt)
    deallocate(marbl_domain%dz)
    deallocate(marbl_domain%temperature)
    deallocate(marbl_domain%salinity)

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

    call ecosys_set_sflux(                                       &
         marbl_saved_state, &
         marbl%private_data%surface_share,                         &
         SHF_QSW_RAW, SHF_QSW,                                     &
         U10_SQR, IFRAC, PRESS,                                    &
         SST, SSS,                                                 &
         SURFACE_VALS_OLD(:,:,ecosys_ind_begin:ecosys_ind_end,:),  &
         SURFACE_VALS_CUR(:,:,ecosys_ind_begin:ecosys_ind_end,:),  &
         STF_MODULE(:,:,ecosys_ind_begin:ecosys_ind_end,:),        &
         ciso_on)

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

    call ecosys_write_restart(marbl_saved_state, restart_file, action)

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
    call ecosys_tavg_write_flux(FLUX_DIAGS)

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
  !BOP
  !IROUTINE: ecosys_driver_tracer_ref_val
  ! !INTERFACE:
  !
  function ecosys_driver_tracer_ref_val(ciso_on,ind)
    !
    ! !DESCRIPTION:
    !  return reference value for tracer with global tracer index ind
    !  this is used in virtual flux computations
    !
    ! !REVISION HISTORY:
    !  same as module

    !INPUT PARAMETERS:
    logical (kind=log_kind), intent(in)  ::  &
         ciso_on                 ! ecosys_ciso on

    integer(int_kind), intent(in) :: ind

    !OUTPUT PARAMETERS:

    real(r8) :: ecosys_driver_tracer_ref_val

    !EOP
    !BOC

    !-----------------------------------------------------------------------
    !  default value for reference value is 0
    !-----------------------------------------------------------------------

    ecosys_driver_tracer_ref_val = c0

    !-----------------------------------------------------------------------
    !  ECOSYS block
    !-----------------------------------------------------------------------

    if (ind >= ecosys_ind_begin .and. ind <= ecosys_ind_end) then
       ecosys_driver_tracer_ref_val = ecosys_tracer_ref_val(ind-ecosys_ind_begin+1)
    endif

    !-----------------------------------------------------------------------
    !  ECOSYS_CISO block
    !-----------------------------------------------------------------------

    if (ciso_on) then
       if (ind >= ecosys_ciso_ind_begin .and. ind <= ecosys_ciso_ind_end) then
          ecosys_driver_tracer_ref_val = ecosys_ciso_tracer_ref_val(ind-ecosys_ciso_ind_begin+1)
       endif
    endif

    !-----------------------------------------------------------------------
    !EOC

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

end module ecosys_driver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

