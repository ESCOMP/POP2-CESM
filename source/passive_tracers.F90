!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module passive_tracers

!BOP
! !MODULE: passive_tracers

! !DESCRIPTION:
!  This module provides support for passive tracers.
!  The base model calls subroutines in this module which then call
!     subroutines in individual passive tracer modules.

! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod

   use kinds_mod, only: r8, int_kind, log_kind, char_len
   use blocks, only: block, nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km, nt
   use domain, only: nblocks_clinic
   use communicate, only: my_task, master_task
   use broadcast, only: broadcast_scalar
   use prognostic, only: TRACER, PSURF, tracer_d, oldtime, curtime, newtime
   use forcing_shf, only: SHF_QSW_RAW, SHF_QSW
   use forcing_fields, only : lhas_vflux, lhas_riv_flux
   use mcog, only: FRACR_BIN, QSW_RAW_BIN, QSW_BIN
   use io_types, only: stdout, nml_in, nml_filename, io_field_desc, &
       datafile
   use exit_mod, only: sigAbort, exit_pop
   use timers, only: timer_start, timer_stop
   use tavg, only: define_tavg_field, tavg_method_qflux,  &
       accumulate_tavg_field, accumulate_tavg_now
   use constants, only: c0, c1, p5, delim_fmt, char_blank, &
       grav, salt_to_ppt, ocn_ref_salinity, ppt_to_salt, sea_ice_salinity
   use time_management, only: mix_pass, c2dtt
   use grid, only: partial_bottom_cells, DZT, KMT, dz, zw, &
       sfc_layer_type, sfc_layer_varthick
   use registry, only: register_string, registry_match
   use io_tools, only: document
   use passive_tracer_tools, only: set_tracer_indices

   use ecosys_driver, only:               &
       ecosys_tracer_cnt,                 &
       ecosys_driver_init,                &
       ecosys_driver_set_sflux_forcing,   &
       ecosys_driver_set_sflux,           &
       ecosys_driver_post_set_sflux,      &
       ecosys_driver_tavg_forcing,        &
       ecosys_driver_set_interior_forcing,&
       ecosys_driver_set_interior,        &
       ecosys_driver_set_global_scalars,  &
       ecosys_driver_comp_global_averages,&
       ecosys_driver_write_restart,       &
       ecosys_driver_tracer_ref_val,      &
       ecosys_driver_unpack_source_sink_terms

   use cfc_mod, only:              &
       cfc_tracer_cnt,             &
       cfc_init,                   &
       cfc_set_sflux,              &
       cfc_tavg_forcing

   use sf6_mod, only:              &
       sf6_tracer_cnt,             &
       sf6_init,                   &
       sf6_set_sflux,              &
       sf6_tavg_forcing

   use iage_mod, only:             &
       iage_tracer_cnt,            &
       iage_init,                  &
       iage_set_interior,          &
       iage_reset

   use abio_dic_dic14_mod, only:      &
       abio_dic_dic14_tracer_cnt,     &
       abio_dic_dic14_init,           &
       abio_dic_dic14_tracer_ref_val, &
       abio_dic_dic14_set_sflux,      &
       abio_dic_dic14_tavg_forcing,   &
       abio_dic_dic14_set_interior,   &
       abio_dic_dic14_write_restart

   use IRF_mod, only: IRF_tracer_cnt
   use IRF_mod, only: IRF_init
   use IRF_mod, only: IRF_reset

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public ::                                  &
      init_passive_tracers,                   &
      set_interior_passive_tracers,           &
      set_interior_passive_tracers_3D,        &
      set_sflux_passive_tracers,              &
      reset_passive_tracers,                  &
      write_restart_passive_tracers,          &
      tavg_passive_tracers,                   &
      passive_tracers_tavg_sflux,             &
      passive_tracers_tavg_fvice,             &
      passive_tracers_timer_print_all,        &
      passive_tracers_send_time,              &
      tracer_ref_val,                         &
      tadvect_ctype_passive_tracers,          &
      ecosys_on

!EOP
!BOC

!-----------------------------------------------------------------------
!  tavg ids for automatically generated tavg passive-tracer fields
!-----------------------------------------------------------------------

   integer (int_kind), dimension(nt), public :: &
      tavg_var_tend,            & ! tavg id for tracer tendency
      tavg_var_tend_zint_100m,  & ! vertically integrated tracer tendency, 0-100m
      tavg_var_rf_tend            ! tavg id for Robert Filter tracer adjustment

   integer (int_kind), dimension (3:nt) ::  &
      tavg_var,                 & ! tracer
      tavg_var_2,               & ! tracer
      tavg_var_z_t_150m,        & ! tracer (top 150 m)
      tavg_var_sqr,             & ! tracer square
      tavg_var_surf,            & ! tracer surface value
      tavg_var_zint_100m,       & ! 0-100m integral of tracer
      tavg_var_J,               & ! tracer source sink term
      tavg_var_Jint,            & ! vertically integrated tracer source sink term
      tavg_var_Jint_100m,       & ! vertically integrated tracer source sink term, 0-100m
      tavg_var_stf,             & ! surface tracer flux
      tavg_var_stf_riv,         & ! riverine tracer flux
      tavg_var_resid,           & ! tracer residual surface flux
      tavg_var_fvper,           & ! virtual tracer flux from precip,evap,runoff
      tavg_var_fvice              ! virtual tracer flux from ice formation

!-----------------------------------------------------------------------
!  array containing advection type for each passive tracer
!-----------------------------------------------------------------------

   character (char_len), dimension(3:nt) :: &
      tadvect_ctype_passive_tracers

!-----------------------------------------------------------------------
!  PER virtual fluxes. The application of the flux happens in surface
!  forcing subroutines, before tavg flags are set, so the tavg accumulation
!  must be in a different subroutine than the application. The fluxes
!  are stored to avoid recomputing them when accumulating.
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: FvPER

!-----------------------------------------------------------------------
!  logical variables that denote if a passive tracer module is on
!-----------------------------------------------------------------------

   logical (log_kind) ::  &
      ecosys_on, cfc_on, sf6_on, iage_on,&
      abio_dic_dic14_on, IRF_on

   namelist /passive_tracers_on_nml/  &
      ecosys_on, cfc_on, sf6_on, iage_on, &
      abio_dic_dic14_on, IRF_on


!-----------------------------------------------------------------------
!     index bounds of passive tracer module variables in TRACER
!-----------------------------------------------------------------------

   integer (int_kind) ::                                 &
      ecosys_driver_ind_begin,   ecosys_driver_ind_end,  &
      iage_ind_begin,            iage_ind_end,           &
      cfc_ind_begin,             cfc_ind_end,            &
      sf6_ind_begin,             sf6_ind_end,            &
      abio_dic_dic14_ind_begin,  abio_dic_dic14_ind_end, &
      IRF_ind_begin,             IRF_ind_end

!-----------------------------------------------------------------------
!  filtered SST and SSS, if needed
!-----------------------------------------------------------------------

   logical (log_kind) :: filtered_SST_SSS_needed

   real (r8), dimension(:,:,:), allocatable :: &
      SST_FILT,      & ! SST with time filter applied, [degC]
      SSS_FILT         ! SSS with time filter applied, [psu]

   real (r8), dimension(:, :, :, :, :), pointer :: &
        ecosys_source_sink_3d ! (nx_block, ny_block, km, nt, nblocks_clinic)

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_passive_tracers
! !INTERFACE:

 subroutine init_passive_tracers(init_ts_file_fmt, &
                                 read_restart_filename, errorCode)

! !DESCRIPTION:
!  Initialize passive tracers. This involves:
!  1) reading passive_tracers_on_nml to see which module are on
!  2) setting tracer module index bounds
!  3) calling tracer module init subroutine
!  4) define common tavg fields
!  5) set up space for storing virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'passive_tracers:init_passive_tracers'

   integer (int_kind) :: cumulative_nt, n, &
      nml_error,        &! error flag for nml read
      iostat             ! io status flag

   character (char_len) :: sname, lname, units, coordinates
   character (4) :: grid_loc

!-----------------------------------------------------------------------

   if (.not. registry_match('init_ts')) then
      call exit_POP(sigAbort, 'init_ts not called ' /&
         &/ 'before init_passive_tracers. This is necessary for ' /&
         &/ 'init_passive_tracers to have correct read_restart_filename')
   end if

!-----------------------------------------------------------------------
!  register init_passive_tracers
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call register_string('init_passive_tracers')

   ecosys_on         = .false.
   cfc_on            = .false.
   sf6_on            = .false.
   iage_on           = .false.
   abio_dic_dic14_on = .false.
   IRF_on            = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
        nml_error = -1
      else
        nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
        read(nml_in, nml=passive_tracers_on_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading passive_tracers_on namelist')
   endif

   if (my_task == master_task) then
      write(stdout,*) ' '
      write(stdout,*) ' Document Namelist Parameters:'
      write(stdout,*) ' ============================ '
      write(stdout,*) ' '
      write(stdout, passive_tracers_on_nml)
      write(stdout,*) ' '
      call POP_IOUnitsFlush(POP_stdout)
   endif


   call broadcast_scalar(ecosys_on,         master_task)
   call broadcast_scalar(cfc_on,            master_task)
   call broadcast_scalar(sf6_on,            master_task)
   call broadcast_scalar(iage_on,           master_task)
   call broadcast_scalar(abio_dic_dic14_on, master_task)
   call broadcast_scalar(IRF_on,            master_task)

!-----------------------------------------------------------------------
!  check for modules that require the flux coupler
!-----------------------------------------------------------------------

   if (cfc_on .and. .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort,'cfc module requires the flux coupler')
   end if

   if (sf6_on .and. .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort,'sf6 module requires the flux coupler')
   end if

   if (abio_dic_dic14_on .and. .not. registry_match('lcoupled')) then
      call exit_POP(sigAbort,'Abiotic DIC_DIC14 module requires the flux coupler')
   end if

!-----------------------------------------------------------------------
!  default is for tracers to use same advection scheme as the base model
!-----------------------------------------------------------------------

   tadvect_ctype_passive_tracers(3:nt) = 'base_model'

!-----------------------------------------------------------------------
!  set up indices for passive tracer modules that are on
!-----------------------------------------------------------------------

   cumulative_nt = 2

   if (ecosys_on) then
      call set_tracer_indices('ECOSYS_DRIVER', ecosys_tracer_cnt, cumulative_nt, &
                              ecosys_driver_ind_begin, ecosys_driver_ind_end)
   end if

   if (cfc_on) then
      call set_tracer_indices('CFC', cfc_tracer_cnt, cumulative_nt,  &
                              cfc_ind_begin, cfc_ind_end)
   end if

   if (sf6_on) then
      call set_tracer_indices('SF6', sf6_tracer_cnt, cumulative_nt,  &
                              sf6_ind_begin, sf6_ind_end)
   end if

   if (iage_on) then
      call set_tracer_indices('IAGE', iage_tracer_cnt, cumulative_nt,  &
                              iage_ind_begin, iage_ind_end)
   end if

   if (abio_dic_dic14_on) then
      call set_tracer_indices('ABIO', abio_dic_dic14_tracer_cnt, cumulative_nt,  &
                              abio_dic_dic14_ind_begin, abio_dic_dic14_ind_end)
   end if

   if (IRF_on) then
      call set_tracer_indices('IRF', IRF_tracer_cnt, cumulative_nt,  &
                              IRF_ind_begin, IRF_ind_end)
   end if

   if (cumulative_nt /= nt) then
      call document(subname, 'nt', nt)
      call document(subname, 'cumulative_nt', cumulative_nt)
      call exit_POP(sigAbort, &
         'ERROR in init_passive_tracers: declared nt does not match cumulative nt')
   end if

!-----------------------------------------------------------------------
!  by default, all tracers are written to tavg as full depth
!-----------------------------------------------------------------------

   tracer_d(1:nt)%lfull_depth_tavg = .true.

!-----------------------------------------------------------------------
!  by default, all tracers have scale_factor equal to one
!-----------------------------------------------------------------------

   tracer_d(3:nt)%scale_factor = 1.0_POP_r8

   ! FIXME (mnl, mvertens) -- for completeness, Mariana wants tracer_d
   ! initialized in ecosys_driver for the ecosystem tracers, which would
   ! lead to the other tracer modules (iage, cfc, IRF, abio) also
   ! needing to initialize lfull_depth_tavg and scale_factor rather than
   ! counting on it being done in passive_tracers.

!-----------------------------------------------------------------------
!  ECOSYS  DRIVER block
!-----------------------------------------------------------------------

   if (ecosys_on) then
      call ecosys_driver_init(                                                           &
           ecosys_driver_ind_begin,                                                      &
           init_ts_file_fmt,                                                             &
           read_restart_filename,                                                        &
           tracer_d(ecosys_driver_ind_begin:ecosys_driver_ind_end),                      &
           TRACER(:,:,:,ecosys_driver_ind_begin:ecosys_driver_ind_end,:,:),              &
           tadvect_ctype_passive_tracers(ecosys_driver_ind_begin:ecosys_driver_ind_end), &
           lhas_riv_flux(ecosys_driver_ind_begin:ecosys_driver_ind_end),                 &
           errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_passive_tracers: error in ecosys_driver_init')
         return
      endif

      allocate(ecosys_source_sink_3d(nx_block, ny_block, km, nt, nblocks_clinic))
      ecosys_source_sink_3d = c0

   end if

!-----------------------------------------------------------------------
!  CFC block
!-----------------------------------------------------------------------

   if (cfc_on) then
      call cfc_init(cfc_ind_begin, init_ts_file_fmt, read_restart_filename, &
                    tracer_d(cfc_ind_begin:cfc_ind_end), &
                    TRACER(:,:,:,cfc_ind_begin:cfc_ind_end,:,:), &
                    errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_passive_tracers: error in cfc_init')
         return
      endif

   end if

!-----------------------------------------------------------------------
!  SF6 block
!-----------------------------------------------------------------------

   if (sf6_on) then
      call sf6_init(sf6_ind_begin, init_ts_file_fmt, read_restart_filename, &
                    tracer_d(sf6_ind_begin:sf6_ind_end), &
                    TRACER(:,:,:,sf6_ind_begin:sf6_ind_end,:,:), &
                    errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_passive_tracers: error in sf6_init')
         return
      endif

   end if

!-----------------------------------------------------------------------
!  Ideal Age (IAGE) block
!-----------------------------------------------------------------------

   if (iage_on) then
      call iage_init(iage_ind_begin, init_ts_file_fmt, read_restart_filename, &
                     tracer_d(iage_ind_begin:iage_ind_end), &
                     TRACER(:,:,:,iage_ind_begin:iage_ind_end,:,:), &
                     errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_passive_tracers: error in iage_init')
         return
      endif

   end if

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 block
!-----------------------------------------------------------------------

   if (abio_dic_dic14_on) then
      call abio_dic_dic14_init(abio_dic_dic14_ind_begin, init_ts_file_fmt, read_restart_filename, &
                    tracer_d(abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end), &
                    TRACER(:,:,:,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,:,:), &
                    errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'init_passive_tracers: error in abio_dic_dic14_init')
         return
      endif

   end if

!-----------------------------------------------------------------------
!  IRF (IRF) block
!-----------------------------------------------------------------------

   if (IRF_on) then
      call IRF_init(tracer_d(IRF_ind_begin:IRF_ind_end), &
                    TRACER(:,:,:,IRF_ind_begin:IRF_ind_end,:,:))
   end if

!-----------------------------------------------------------------------
!  print out tracer names from tracer modules that are on
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,*) 'TRACER INDEX    TRACER NAME'
      write(stdout,1010) 1, 'TEMP'
      write(stdout,1010) 2, 'SALT'
      call POP_IOUnitsFlush(POP_stdout)
      do n = 3, nt
         write(stdout,1010) n, TRIM(tracer_d(n)%long_name)
         call POP_IOUnitsFlush(POP_stdout)
      enddo
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout)
   end if

!-----------------------------------------------------------------------
!  set lhas_vflux for passive tracers
!-----------------------------------------------------------------------

   do n = 3, nt
      lhas_vflux(n) = tracer_ref_val(n) /= c0
   enddo

!-----------------------------------------------------------------------
!  generate common tavg fields for all tracers
!-----------------------------------------------------------------------

   do n = 3, nt
      sname = tracer_d(n)%short_name
      lname = tracer_d(n)%long_name
      units = tracer_d(n)%units
      if (tracer_d(n)%lfull_depth_tavg) then
         grid_loc = '3111'
         coordinates = 'TLONG TLAT z_t time'
      else
         grid_loc = '3114'
         coordinates = 'TLONG TLAT z_t_150m time'
      end if
      call define_tavg_field(tavg_var(n),                           &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates=coordinates)

      sname = trim(tracer_d(n)%short_name) /&
                                            &/ '_2'
      call define_tavg_field(tavg_var_2(n),                         &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates=coordinates)

      sname = trim(tracer_d(n)%short_name) /&
                                            &/ '_z_t_150m'
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ', top 150m'
      call define_tavg_field(tavg_var_z_t_150m(n),                  &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc='3114',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT z_t_150m time')

      sname = trim(tracer_d(n)%short_name) /&
                                            &/ '_SQR'
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Squared'
      units = '(' /&
                   &/ trim(tracer_d(n)%units) /&
                                               &/ ')^2'
      call define_tavg_field(tavg_var_sqr(n),                       &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor**2,&
                             coordinates=coordinates)

      sname = trim(tracer_d(n)%short_name) /&
                                            &/ '_SURF'
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Surface Value'
      units = tracer_d(n)%units
      call define_tavg_field(tavg_var_surf(n),                      &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      sname = trim(tracer_d(n)%short_name) /&
                                            &/ '_zint_100m'
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' 0-100m Vertical Integral'
      units = trim(tracer_d(n)%units) /&
                                       &/ ' cm'
      call define_tavg_field(tavg_var_zint_100m(n),                 &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      sname = 'J_' /&
                    &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Source Sink Term'
      units = tracer_d(n)%tend_units
      call define_tavg_field(tavg_var_J(n),                         &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates=coordinates)

      sname = 'Jint_' /&
                       &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Source Sink Term Vertical Integral'
      units = tracer_d(n)%flux_units
      call define_tavg_field(tavg_var_Jint(n),                      &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      sname = 'Jint_100m_' /&
                            &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Source Sink Term Vertical Integral, 0-100m'
      units = tracer_d(n)%flux_units
      call define_tavg_field(tavg_var_Jint_100m(n),                 &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      sname = 'STF_' /&
                      &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Surface Flux, excludes FvICE term'
      units = tracer_d(n)%flux_units
      call define_tavg_field(tavg_var_stf(n),                       &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      if (lhas_riv_flux(n)) then
         sname = trim(tracer_d(n)%short_name) /&
                                               &/ '_RIV_FLUX'
         lname = trim(tracer_d(n)%long_name) /&
                                              &/ ' Riverine Flux'
         units = tracer_d(n)%flux_units
         call define_tavg_field(tavg_var_stf_riv(n),                   &
                                sname, 2, long_name=lname,             &
                                units=units, grid_loc='2110',          &
                                scale_factor=tracer_d(n)%scale_factor, &
                                coordinates='TLONG TLAT time')
      endif

      sname = 'RESID_' /&
                        &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Residual Surface Flux'
      units = tracer_d(n)%flux_units
      call define_tavg_field(tavg_var_resid(n),                     &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      if (lhas_vflux(n)) then
         sname = 'FvPER_' /&
                           &/ trim(tracer_d(n)%short_name)
         lname = trim(tracer_d(n)%long_name) /&
                                              &/ ' Virtual Surface Flux, PER'
         units = tracer_d(n)%flux_units
         call define_tavg_field(tavg_var_fvper(n),                     &
                                sname, 2, long_name=lname,             &
                                units=units, grid_loc='2110',          &
                                scale_factor=tracer_d(n)%scale_factor, &
                                coordinates='TLONG TLAT time')

         sname = 'FvICE_' /&
                           &/ trim(tracer_d(n)%short_name)
         lname = trim(tracer_d(n)%long_name) /&
                                              &/ ' Virtual Surface Flux, ICE'
         units = tracer_d(n)%flux_units
         call define_tavg_field(tavg_var_fvice(n),                     &
                                sname, 2, long_name=lname,             &
                                units=units, grid_loc='2110',          &
                                scale_factor=tracer_d(n)%scale_factor, &
                                tavg_method=tavg_method_qflux,         &
                                coordinates='TLONG TLAT time')
      endif
   enddo

   do n=1,nt
      sname = 'TEND_' /&
           &/ trim(tracer_d(n)%short_name)
      lname = 'Tendency of Thickness Weighted '/&
           &/ trim(tracer_d(n)%short_name)
      units = tracer_d(n)%tend_units
      if (tracer_d(n)%lfull_depth_tavg) then
         grid_loc = '3111'
         coordinates = 'TLONG TLAT z_t time'
      else
         grid_loc = '3114'
         coordinates = 'TLONG TLAT z_t_150m time'
      end if

      call define_tavg_field(tavg_var_tend(n),                      &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates=coordinates)

      sname = 'tend_zint_100m_' /&
                            &/ trim(tracer_d(n)%short_name)
      lname = trim(tracer_d(n)%long_name) /&
                                           &/ ' Tendency Vertical Integral, 0-100m'
      units = tracer_d(n)%flux_units
      call define_tavg_field(tavg_var_tend_zint_100m(n),            &
                             sname, 2, long_name=lname,             &
                             units=units, grid_loc='2110',          &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates='TLONG TLAT time')

      sname = 'RF_TEND_' /&
           &/ trim(tracer_d(n)%short_name)
      lname = 'Robert Filter Tendency for '/&
           &/ trim(tracer_d(n)%short_name)
      units = tracer_d(n)%tend_units

      call define_tavg_field(tavg_var_rf_tend(n),                   &
                             sname, 3, long_name=lname,             &
                             units=units, grid_loc=grid_loc,        &
                             scale_factor=tracer_d(n)%scale_factor, &
                             coordinates=coordinates)

   end do

!-----------------------------------------------------------------------
!  allocate and initialize storage for virtual fluxes
!-----------------------------------------------------------------------

   allocate(FvPER(nx_block,ny_block,3:nt,nblocks_clinic))
   FvPER = c0

!-----------------------------------------------------------------------
!  allocate space for filtered SST and SSS, if needed
!-----------------------------------------------------------------------

   filtered_SST_SSS_needed = ecosys_on .or. cfc_on .or. sf6_on .or. &
                             abio_dic_dic14_on

   if (filtered_SST_SSS_needed) then
      allocate(SST_FILT(nx_block,ny_block,max_blocks_clinic), &
               SSS_FILT(nx_block,ny_block,max_blocks_clinic))
   endif

 1010 format(5X,I2,10X,A)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: set_interior_passive_tracers
! !INTERFACE:

 subroutine set_interior_passive_tracers(k, this_block, TRACER_SOURCE)

! !DESCRIPTION:
!  call subroutines for each tracer module that compute source-sink terms
!  accumulate commnon tavg fields related to source-sink terms
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block information for this block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      TRACER_SOURCE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      bid,                &! local block address for this block
      n                    ! tracer index

   real (r8) :: &
      ztop                 ! depth of top of cell

   real (r8), dimension(nx_block,ny_block) :: &
      WORK

!-----------------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------
!  ECOSYS DRIVER block is done as part of
!  set_interior_passive_tracers_3D. Here we are just unpacking the 3D
!  structure into the 2D
!  -----------------------------------------------------------------------
   if (ecosys_on) then
      call ecosys_driver_unpack_source_sink_terms( &
           ecosys_source_sink_3d(:, :, k, ecosys_driver_ind_begin:ecosys_driver_ind_end, bid), &
           TRACER_SOURCE(:, :, ecosys_driver_ind_begin:ecosys_driver_ind_end))
   endif

!-----------------------------------------------------------------------
!  CFC does not have source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not have source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Ideal Age (IAGE) block
!-----------------------------------------------------------------------

   if (iage_on) then
      call iage_set_interior(k,                                    &
         TRACER_SOURCE (:,:,iage_ind_begin:iage_ind_end) )
   end if

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 block
!-----------------------------------------------------------------------

   if (abio_dic_dic14_on) then
      call abio_dic_dic14_set_interior(k,                                  &
         TRACER(:,:,:,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,oldtime,bid),&
         TRACER(:,:,:,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,curtime,bid),&
         TRACER_SOURCE(:,:,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end),       &
         this_block)
    end if

!-----------------------------------------------------------------------
!  IRF does not have source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  accumulate time average if necessary
!-----------------------------------------------------------------------

   if (mix_pass /= 1) then
      do n = 3, nt
         call accumulate_tavg_field(TRACER_SOURCE(:,:,n),tavg_var_J(n),bid,k)

         if (accumulate_tavg_now(tavg_var_Jint(n))) then
               if (partial_bottom_cells) then
                  WORK = merge(DZT(:,:,k,bid) * TRACER_SOURCE(:,:,n), &
                               c0, k<=KMT(:,:,bid))
               else
                  WORK = merge(dz(k) * TRACER_SOURCE(:,:,n), &
                               c0, k<=KMT(:,:,bid))
               endif
               call accumulate_tavg_field(WORK,tavg_var_Jint(n),bid,k)
         endif
      enddo

      ztop = c0
      if (k > 1) ztop = zw(k-1)
      if (ztop < 100.0e2_r8) then
         do n = 3, nt
            if (accumulate_tavg_now(tavg_var_Jint_100m(n))) then
                  if (partial_bottom_cells) then
                     WORK = merge(min(100.0e2_r8 - ztop, DZT(:,:,k,bid)) &
                                  * TRACER_SOURCE(:,:,n), c0, k<=KMT(:,:,bid))
                  else
                     WORK = merge(min(100.0e2_r8 - ztop, dz(k)) &
                                  * TRACER_SOURCE(:,:,n), c0, k<=KMT(:,:,bid))
                  endif
                  call accumulate_tavg_field(WORK,tavg_var_Jint_100m(n),bid,k)
            endif
         enddo
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine set_interior_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: set_interior_passive_tracers_3D
! !INTERFACE:

 subroutine set_interior_passive_tracers_3D (TRACER_OLD, TRACER_CUR)

! !DESCRIPTION:
!  call subroutines for each tracer module that computes 3D source-sink terms
!  accumulate commnon tavg fields related to source-sink terms
!
! !REVISION HISTORY:
!  same as module

   use domain, only : blocks_clinic
   use blocks, only : get_block

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,km,nt,max_blocks_clinic), intent(in) :: &
        TRACER_OLD, & ! previous timestep tracer tendencies
        TRACER_CUR    ! new tracer tendencies

! !INPUT/OUTPUT PARAMETERS:

!EOP
!BOC

   integer (int_kind) :: iblock ! counter for block loops
   type (block) :: this_block   ! block information for this block
   integer (int_kind) :: bid ! local block address for this block

!-----------------------------------------------------------------------
!  ECOSYS DRIVER modules 3D source-sink terms
!-----------------------------------------------------------------------
   if (ecosys_on) then
      call ecosys_driver_set_global_scalars('interior_tendency')

      call ecosys_driver_set_interior_forcing(FRACR_BIN, QSW_RAW_BIN, QSW_BIN)

      !$OMP PARALLEL DO PRIVATE(iblock, this_block, bid)
      do iblock = 1, nblocks_clinic

         this_block = get_block(blocks_clinic(iblock), iblock)
         bid = this_block%local_id

         call ecosys_driver_set_interior(&
              TRACER(:, :, :, ecosys_driver_ind_begin:ecosys_driver_ind_end, oldtime, bid), &
              TRACER(:, :, :, ecosys_driver_ind_begin:ecosys_driver_ind_end, curtime, bid), &
              ecosys_source_sink_3d(:, :, :, ecosys_driver_ind_begin:ecosys_driver_ind_end, bid), &
              this_block)

      end do
      !$OMP END PARALLEL DO

      call ecosys_driver_comp_global_averages('interior_tendency')
   end if

!-----------------------------------------------------------------------
!  CFC does not compute and store 3D source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not compute and store 3D source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Ideal Age (IAGE) does not compute and store 3D source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC 14 does not compute and store 3D source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IRF does not compute and store 3D source-sink terms
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  accumulate time average if necessary
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!EOC

 end subroutine set_interior_passive_tracers_3D

!***********************************************************************
!BOP
! !IROUTINE: set_sflux_passive_tracers
! !INTERFACE:

 subroutine set_sflux_passive_tracers(U10_SQR,ICE_FRAC,PRESS,ATM_FINE_DUST_FLUX,ATM_COARSE_DUST_FLUX,SEAICE_DUST_FLUX, &
                                      ATM_BLACK_CARBON_FLUX,SEAICE_BLACK_CARBON_FLUX, &
                                      lvsf_river,MASK_ESTUARY,vsf_river_correction,STF,STF_RIV)

! !DESCRIPTION:
!  call subroutines for each tracer module that compute surface fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,                  & ! 10m wind speed squared
      ICE_FRAC,                 & ! sea ice fraction (non-dimensional)
      PRESS,                    & ! sea level atmospheric pressure (Pascals)
      ATM_FINE_DUST_FLUX,       & ! fine dust flux from atm (g/cm^2/s)
      ATM_COARSE_DUST_FLUX,     & ! coarse dust flux from atm (g/cm^2/s)
      SEAICE_DUST_FLUX,         & ! coarse dust flux from seaice (g/cm^2/s)
      ATM_BLACK_CARBON_FLUX,    & ! black carbon flux from atm (g/cm^2/s)
      SEAICE_BLACK_CARBON_FLUX, & ! black carbon flux from seaice (g/cm^2/s)
      MASK_ESTUARY                ! mask for estuary model, 1 where it runs and 0 elsewhere

  logical (log_kind), intent(in) :: &
      lvsf_river

   real (r8), dimension(nt), intent(in) :: &
      vsf_river_correction ! per-tracer correction for using local tracer concenctration in application of ROFF_F

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), intent(inout) :: &
      STF,               & ! surface tracer fluxes
      STF_RIV              ! riverine tracer fluxes

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   logical (log_kind) :: first_call = .true.
   real (r8)          :: ref_val
   integer (int_kind) :: iblock, n
   real (r8), dimension(nx_block,ny_block,nblocks_clinic) :: STF2_m_river_correction

!-----------------------------------------------------------------------

   if (first_call) then
      call register_string('set_sflux_passive_tracers')
   end if

!-----------------------------------------------------------------------
!  compute filtered SST and SSS, if needed
!-----------------------------------------------------------------------

   if (filtered_SST_SSS_needed) then
      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic
         SST_FILT(:,:,iblock) = p5*(TRACER(:,:,1,1,oldtime,iblock) + &
                                    TRACER(:,:,1,1,curtime,iblock))
         SSS_FILT(:,:,iblock) = p5*(TRACER(:,:,1,2,oldtime,iblock) + &
                                    TRACER(:,:,1,2,curtime,iblock)) * salt_to_ppt
      end do
      !$OMP END PARALLEL DO
   end if

!-----------------------------------------------------------------------
!  ECOSYS  DRIVER block
!-----------------------------------------------------------------------

   if (ecosys_on) then
      call ecosys_driver_set_sflux_forcing( &
         U10_SQR, ICE_FRAC, PRESS, ATM_FINE_DUST_FLUX, ATM_COARSE_DUST_FLUX, SEAICE_DUST_FLUX, &
         ATM_BLACK_CARBON_FLUX, SEAICE_BLACK_CARBON_FLUX, &
         SST_FILT, SSS_FILT)

      call ecosys_driver_set_global_scalars('surface_flux')

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic
         call ecosys_driver_set_sflux(                                                  &
            TRACER(:,:,1,ecosys_driver_ind_begin:ecosys_driver_ind_end,oldtime,iblock), &
            TRACER(:,:,1,ecosys_driver_ind_begin:ecosys_driver_ind_end,curtime,iblock), &
            STF(:,:,ecosys_driver_ind_begin:ecosys_driver_ind_end,iblock),              &
            STF_RIV(:,:,ecosys_driver_ind_begin:ecosys_driver_ind_end,iblock), iblock)
         end do
      !$OMP END PARALLEL DO

      call ecosys_driver_post_set_sflux()

      call ecosys_driver_comp_global_averages('surface_flux')
   end if

!-----------------------------------------------------------------------
!  CFC block
!-----------------------------------------------------------------------

   if (cfc_on) then
      call cfc_set_sflux(U10_SQR, ICE_FRAC, PRESS,                 &
         SST_FILT, SSS_FILT,                                       &
         TRACER(:,:,1,cfc_ind_begin:cfc_ind_end,oldtime,:),        &
         TRACER(:,:,1,cfc_ind_begin:cfc_ind_end,curtime,:),        &
         STF(:,:,cfc_ind_begin:cfc_ind_end,:))
   end if

!-----------------------------------------------------------------------
!  SF6 block
!-----------------------------------------------------------------------

   if (sf6_on) then
      call sf6_set_sflux(U10_SQR, ICE_FRAC, PRESS,                 &
         SST_FILT, SSS_FILT,                                       &
         TRACER(:,:,1,sf6_ind_begin:sf6_ind_end,oldtime,:),        &
         TRACER(:,:,1,sf6_ind_begin:sf6_ind_end,curtime,:),        &
         STF(:,:,sf6_ind_begin:sf6_ind_end,:))
   end if

!-----------------------------------------------------------------------
!  IAGE does not have surface fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 block
!-----------------------------------------------------------------------

   if (abio_dic_dic14_on) then
      call abio_dic_dic14_set_sflux(U10_SQR, ICE_FRAC, PRESS,                 &
         SST_FILT, SSS_FILT,                                       &
         TRACER(:,:,1,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,oldtime,:),        &
         TRACER(:,:,1,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,curtime,:),        &
         STF(:,:,abio_dic_dic14_ind_begin:abio_dic_dic14_ind_end,:))
   end if


!-----------------------------------------------------------------------
!  IRF does not have surface fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  add virtual fluxes for tracers that use them
!  add STF_RIV to STF where the ebm is not handling it
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,n,ref_val)
   do iblock = 1,nblocks_clinic
      ! If lvsf_river is enabled, do not include vsf_river_correction(2) when
      ! computing tracer virtual flux. vsf_river_correction is only present in STF(2)
      ! where MASK_ESTUARY=1.
      if (lvsf_river) then
         STF2_m_river_correction(:,:,iblock) = STF(:,:,2,iblock) - &
            MASK_ESTUARY(:,:,iblock)*vsf_river_correction(2)
      else
         STF2_m_river_correction(:,:,iblock) = STF(:,:,2,iblock)
      endif
      do n=3,nt
         if (lhas_vflux(n)) then
            ref_val = tracer_ref_val(n)
            FvPER(:,:,n,iblock) = &
               (ref_val/(ocn_ref_salinity*ppt_to_salt)) * STF2_m_river_correction(:,:,iblock)
            ! Add global correction for salt conservation, correcting for using
            ! local tracer concentration in application of ROFF_F.
            ! Correction is applied where MASK_ESTUARY=1.
            if (lvsf_river) then
               FvPER(:,:,n,iblock) = FvPER(:,:,n,iblock) + &
                  MASK_ESTUARY(:,:,iblock)*vsf_river_correction(n)
            endif
            STF(:,:,n,iblock) = STF(:,:,n,iblock) + FvPER(:,:,n,iblock)
         endif
         if (lhas_riv_flux(n)) then
            STF(:,:,n,iblock) = STF(:,:,n,iblock) + &
               (c1-MASK_ESTUARY(:,:,iblock))*STF_RIV(:,:,n,iblock)
         endif
      end do
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------

   first_call = .false.

!-----------------------------------------------------------------------
!EOC

 end subroutine set_sflux_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: write_restart_passive_tracers
! !INTERFACE:

 subroutine write_restart_passive_tracers(restart_file, action)

! !DESCRIPTION:
!  call restart routines for each tracer module that
!  write fields besides the tracers themselves
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC

!-----------------------------------------------------------------------
!  ECOSYS  DRIVER block
!-----------------------------------------------------------------------

   if (ecosys_on) then
      call ecosys_driver_write_restart(restart_file, action)
   end if

!-----------------------------------------------------------------------
!  CFC does not write additional restart fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not write additional restart fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IAGE does not write additional restart fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14  block
!-----------------------------------------------------------------------
 if (abio_dic_dic14_on) then
      call abio_dic_dic14_write_restart(restart_file, action)
   end if

!-----------------------------------------------------------------------
!  IRF does not write additional restart fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

 end subroutine write_restart_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: reset_passive_tracers
! !INTERFACE:

 subroutine reset_passive_tracers(TRACER_OLD, TRACER_RESET, bid)

! !DESCRIPTION:
!  call subroutines for each tracer module to reset tracer values
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: bid

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
      TRACER_OLD,   & ! all tracers at old time for a given block
      TRACER_RESET    ! all tracers at time level being reset for a given block

!EOP
!BOC

!-----------------------------------------------------------------------
!  ECOSYS DRIVER modules do not reset values
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  CFC does not reset values
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not reset values
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IAGE block
!-----------------------------------------------------------------------

   if (iage_on) then
      call iage_reset(  &
         TRACER_RESET(:,:,:,iage_ind_begin:iage_ind_end), bid)
   end if

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 does not reset values
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IRF block
!-----------------------------------------------------------------------

   if (IRF_on) then
      call IRF_reset( &
         TRACER_OLD(:,:,:,IRF_ind_begin:IRF_ind_end), &
         TRACER_RESET(:,:,:,IRF_ind_begin:IRF_ind_end) )
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: tavg_passive_tracers
! !INTERFACE:

 subroutine tavg_passive_tracers(bid, k)

! !DESCRIPTION:
!  accumulate common tavg fields for tracers
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k, bid  ! vertical level index

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n                    ! tracer index

   real (r8) :: &
      ztop                 ! depth of top of cell

   real (r8), dimension(nx_block,ny_block) :: &
      WORK

!-----------------------------------------------------------------------


   if (mix_pass /= 1) then
      do n = 3, nt
         call accumulate_tavg_field(TRACER(:,:,k,n,curtime,bid),tavg_var(n),bid,k)
         call accumulate_tavg_field(TRACER(:,:,k,n,curtime,bid),tavg_var_2(n),bid,k)
         call accumulate_tavg_field(TRACER(:,:,k,n,curtime,bid),tavg_var_z_t_150m(n),bid,k)

         if (accumulate_tavg_now(tavg_var_sqr(n))) then
            WORK = TRACER(:,:,k,n,curtime,bid) ** 2
            call accumulate_tavg_field(WORK,tavg_var_sqr(n),bid,k)
         endif
      enddo

      if (k == 1) then
         do n = 3, nt
            call accumulate_tavg_field(TRACER(:,:,k,n,curtime,bid), &
                                       tavg_var_surf(n),bid,k)
         enddo
      endif

      ztop = c0
      if (k > 1) ztop = zw(k-1)
      if (ztop < 100.0e2_r8) then
         do n = 3, nt
            if (accumulate_tavg_now(tavg_var_zint_100m(n))) then
                  if (sfc_layer_type == sfc_layer_varthick .and. k == 1) then
                     WORK = merge((dz(k)+PSURF(:,:,curtime,bid)/grav) &
                                  * TRACER(:,:,k,n,curtime,bid), c0, k<=KMT(:,:,bid))
                  else
                     if (partial_bottom_cells) then
                        WORK = merge(min(100.0e2_r8 - ztop, DZT(:,:,k,bid)) &
                                     * TRACER(:,:,k,n,curtime,bid), c0, k<=KMT(:,:,bid))
                     else
                        WORK = merge(min(100.0e2_r8 - ztop, dz(k)) &
                                     * TRACER(:,:,k,n,curtime,bid), c0, k<=KMT(:,:,bid))
                     endif
                  endif
                  call accumulate_tavg_field(WORK,tavg_var_zint_100m(n),bid,k)
            endif
         enddo
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_passive_tracers

!***********************************************************************
!BOP
! !IROUTINE: passive_tracers_tavg_sflux
! !INTERFACE:

 subroutine passive_tracers_tavg_sflux(STF, STF_RIV, lvsf_river, MASK_ESTUARY, FLUX_ROFF_VSF_SRF)

! !DESCRIPTION:
!  accumulate common tavg fields for surface tracer fluxes
!  call accumation subroutines for tracer modules that have additional
!     tavg fields related to surface fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), intent(in) :: &
      STF,               & ! surface tracer fluxes
      STF_RIV              ! riverine tracer fluxes

  logical(log_kind), intent(in) :: &
      lvsf_river

  real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      MASK_ESTUARY         ! mask for estuary model, 1 where it runs and 0 elsewhere

  real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), intent(in) :: &
      FLUX_ROFF_VSF_SRF    !

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock, n

   real (r8), dimension(nx_block,ny_block) :: WORK

!-----------------------------------------------------------------------
!  accumulate surface flux and FvPER flux for all tracers
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,n,WORK)
   do iblock = 1,nblocks_clinic
      do n = 3, nt
         call accumulate_tavg_field(STF(:,:,n,iblock),tavg_var_stf(n),iblock,1)
         if (lhas_riv_flux(n)) then
            call accumulate_tavg_field(STF_RIV(:,:,n,iblock),tavg_var_stf_riv(n),iblock,1)

            ! explicitly include STF_RIV in STF_var where ebm is handling it
            if (lvsf_river) then
              WORK = MASK_ESTUARY(:,:,iblock)*STF_RIV(:,:,n,iblock)
              call accumulate_tavg_field(WORK,tavg_var_stf(n),iblock,1)
            endif
         endif
         if (lhas_vflux(n)) then
            call accumulate_tavg_field(FvPER(:,:,n,iblock),tavg_var_fvper(n),iblock,1)

            ! explicitly include FLUX_ROFF_VSF_SRF in STF_var and FvPER_var
            ! note that FLUX_ROFF_VSF_SRF is positive up and that MASK_ESTUARY is already applied to it
            if (lvsf_river) then
              WORK = -FLUX_ROFF_VSF_SRF(:,:,n,iblock)
              call accumulate_tavg_field(WORK,tavg_var_stf(n),iblock,1)
              call accumulate_tavg_field(WORK,tavg_var_fvper(n),iblock,1)
            endif
         endif
      enddo
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  call routines from modules that have additional sflux tavg fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ECOSYS  DRIVER block
!-----------------------------------------------------------------------

   if (ecosys_on) then
     !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
       call ecosys_driver_tavg_forcing(STF(:,:,ecosys_driver_ind_begin:ecosys_driver_ind_end,iblock), iblock)
     end do
     !$OMP END PARALLEL DO
   end if

!-----------------------------------------------------------------------
!  CFC block
!-----------------------------------------------------------------------

   if (cfc_on) then
      call cfc_tavg_forcing
   end if

!-----------------------------------------------------------------------
!  SF6 block
!-----------------------------------------------------------------------

   if (sf6_on) then
      call sf6_tavg_forcing
   end if

!-----------------------------------------------------------------------
!  IAGE does not have additional sflux tavg fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 block
!-----------------------------------------------------------------------

   if (abio_dic_dic14_on) then
      call abio_dic_dic14_tavg_forcing
   end if

!-----------------------------------------------------------------------
!  IRF does not have additional sflux tavg fields
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

 end subroutine passive_tracers_tavg_sflux

!***********************************************************************
!BOP
! !IROUTINE: passive_tracers_tavg_FvICE
! !INTERFACE:

 subroutine passive_tracers_tavg_FvICE(cp_over_lhfusion, QICE)

! !DESCRIPTION:
!  accumulate FvICE fluxes passive tracers
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      QICE                ! tot column cooling from ice form (in C*cm)

   real (r8), intent(in) ::  &
      cp_over_lhfusion    ! cp_sw/latent_heat_fusion

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block) :: &
      WORK

   real (r8) ::        &
      ref_val   ! temporary work array

   integer (int_kind) :: iblock, n

!-----------------------------------------------------------------------


   !$OMP PARALLEL DO PRIVATE(iblock,n,ref_val,WORK)
   do iblock = 1,nblocks_clinic
      do n = 3, nt
         if (lhas_vflux(n)) then
           ref_val = tracer_ref_val(n)
           WORK = ref_val * (c1 - sea_ice_salinity / ocn_ref_salinity) * &
              cp_over_lhfusion * max(c0, QICE(:,:,iblock))
           call accumulate_tavg_field(WORK,tavg_var_fvice(n),iblock,1,c1)
         endif
      enddo
   enddo
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine passive_tracers_tavg_FvICE

!***********************************************************************
!BOP
! !IROUTINE: tracer_ref_val
! !INTERFACE:

 function tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracer with global tracer index ind
!  this is used in virtual flux computations
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: tracer_ref_val

!EOP
!BOC

!-----------------------------------------------------------------------
!  default value for reference value is 0
!-----------------------------------------------------------------------

   tracer_ref_val = c0

!-----------------------------------------------------------------------
!  ECOSYS  DRIVER block
!-----------------------------------------------------------------------

   if (ecosys_on) then
      if (ind >= ecosys_driver_ind_begin .and. ind <= ecosys_driver_ind_end) then
         tracer_ref_val = ecosys_driver_tracer_ref_val(ind - ecosys_driver_ind_begin + 1)
      endif
   endif

!-----------------------------------------------------------------------
!  CFC does not use virtual fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not use virtual fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IAGE does not use virtual fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 block
!-----------------------------------------------------------------------

   if (abio_dic_dic14_on) then
      if (ind >= abio_dic_dic14_ind_begin .and. ind <= abio_dic_dic14_ind_end) then
         tracer_ref_val = abio_dic_dic14_tracer_ref_val(ind-abio_dic_dic14_ind_begin+1)
      endif
   endif


!-----------------------------------------------------------------------
!EOC

 end function tracer_ref_val

!***********************************************************************
!BOP
! !IROUTINE: passive_tracers_timer_print_all
! !INTERFACE:

 subroutine passive_tracers_timer_print_all(stats)

! !DESCRIPTION:
!  If passive tracer packages provide their own timing numbers (rather
!  than using POP's timers), print them here

! !REVISION HISTORY:
!  same as module

   use ecosys_driver, only : ecosys_driver_print_marbl_timers

! !INPUT PARAMETERS:

   logical (log_kind), optional, intent(in) :: stats

!EOP
!BOC

!-----------------------------------------------------------------------
!  ECOSYS (and CISO) timers come from MARBL library
!-----------------------------------------------------------------------

  if (ecosys_on) then
    call ecosys_driver_print_marbl_timers(stats)
  end if

!-----------------------------------------------------------------------
!  CFC does not have separate timers
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  SF6 does not have separate timers
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IAGE does not have separate timers
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ABIO DIC & DIC14 do not have separate timers
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  IRF does not have separate timers
!-----------------------------------------------------------------------

!EOC

 end subroutine passive_tracers_timer_print_all

!***********************************************************************
!BOP
! !IROUTINE: passive_tracers_send_time
! !INTERFACE:

 subroutine passive_tracers_send_time

! !DESCRIPTION:
!  sends POP time information
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:


!EOP
!BOC

!-----------------------------------------------------------------------
!  IRF does not use virtual fluxes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!EOC

end subroutine passive_tracers_send_time

!***********************************************************************

 end module passive_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
