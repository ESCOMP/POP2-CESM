!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module sf6_mod

!BOP
! !MODULE: sf6_mod
!
!  Module for SF6
!
!  The units of concentration for these tracers are
!     fmol/cm^3 == nmol/m^3 == pmol/l ~= pmol/kg.
!  These units are chosen because ship measurements are typically
!  given in units of pmol/kg.
!
!  The units of surface fluxes for these tracers are
!     fmol/cm^3 * cm/s == fmol/cm^2/s.
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks,      only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km
   use domain,      only: nblocks_clinic, distrb_clinic
   use exit_mod,    only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants,   only: c0, c1
   use io_types,    only: stdout
   use io_tools,    only: document
   use tavg,        only: define_tavg_field, accumulate_tavg_field
   use time_management, only : iyear, iday_of_year, frac_day, days_in_year

   use passive_tracer_tools, only: forcing_monthly_every_ts, ind_name_pair
   use passive_tracer_tools, only : read_field, tracer_read
   use forcing_timeseries_mod, only: forcing_timeseries_dataset
   use broadcast

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
       sf6_tracer_cnt, &
       sf6_init, &
       sf6_set_sflux,  &
       sf6_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       sf6_tracer_cnt = 1

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      sf6_ind =  1 ! SF6


!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------
   
   type(ind_name_pair), dimension(sf6_tracer_cnt) :: &
        ind_name_table = (/ &
        ind_name_pair(sf6_ind, 'SF6')/)
   
!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) :: &
      sf6_formulation,     & ! how to calculate flux (ocmip or model)
      psf6_file              ! filename for netCDF time series of atm psf6

   integer (int_kind) ::  &
      model_year,             & ! arbitrary model year
      data_year,              & ! year in data that corresponds to model_year
      psf6_first_nonzero_year   ! first year of non-zero values in psf6_file

   type (forcing_timeseries_dataset) :: &
      psf6_atm_forcing_dataset  ! data structure for atm psf6 timeseries

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK            ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      fice_file,           & ! ice fraction, if read from file
      xkw_file,            & ! a * wind-speed ** 2, if read from file
      ap_file                ! atmoshperic pressure, if read from file

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      SF6_SFLUX_TAVG

   integer (int_kind) ::  &
      tavg_SF6_IFRAC,     & ! tavg id for ice fraction
      tavg_SF6_XKW,       & ! tavg id for xkw
      tavg_SF6_ATM_PRESS, & ! tavg id for atmospheric pressure
      tavg_pSF6,          & ! tavg id for sf6 partial pressure
      tavg_SF6_SCHMIDT,   & ! tavg id for sf6 Schmidt number
      tavg_SF6_PV,        & ! tavg id for sf6 piston velocity
      tavg_SF6_surf_sat    ! tavg id for sf6 surface saturation

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: sf6_sflux_timer

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: sf6_init
! !INTERFACE:

 subroutine sf6_init(sf6_ind_begin, init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize sf6 tracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants,  only: char_blank, delim_fmt
   use prognostic, only: curtime, oldtime, tracer_field
   use grid,       only: KMT, n_topo_smooth, fill_points
   use io_types,   only: nml_in, nml_filename
   use timers,     only: get_timer

   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block

   use io_read_fallback_mod, only: io_read_fallback_register_tracer

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      sf6_ind_begin          ! starting index of sf6 tracers in global tracer array
                             ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(sf6_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,sf6_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'sf6_mod:sf6_init'

   character(char_len) :: &
      init_sf6_option,        & ! option for initialization of bgc
      init_sf6_init_file,     & ! filename for option 'file'
      init_sf6_init_file_fmt    ! file format for option 'file'

   integer (int_kind) :: &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(sf6_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   type(tracer_read) :: &
      gas_flux_fice,          & ! ice fraction for gas fluxes
      gas_flux_ws,            & ! wind speed for gas fluxes
      gas_flux_ap               ! atmospheric pressure for gas fluxes

   namelist /sf6_nml/ &
      init_sf6_option, init_sf6_init_file, init_sf6_init_file_fmt, &
      tracer_init_ext, psf6_file, psf6_first_nonzero_year, model_year, data_year, &
      sf6_formulation, gas_flux_fice, gas_flux_ws, gas_flux_ap

   real (r8) :: &
      mapped_date               ! date of current model timestep mapped to data timeline

   character (char_len) ::  &
      sf6_restart_filename      ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   do n = 1, sf6_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'fmol/cm^3'
      tracer_d_module(n)%tend_units = 'fmol/cm^3/s'
      tracer_d_module(n)%flux_units = 'fmol/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_sf6_option        = 'unknown'
   init_sf6_init_file     = 'unknown'
   init_sf6_init_file_fmt = 'bin'

   do n = 1, sf6_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   psf6_file               = 'unknown'
   psf6_first_nonzero_year = 1953
   model_year              = 1
   data_year               = 1
   sf6_formulation         = 'model'

   gas_flux_fice%filename     = 'unknown'
   gas_flux_fice%file_varname = 'FICE'
   gas_flux_fice%scale_factor = c1
   gas_flux_fice%default_val  = c0
   gas_flux_fice%file_fmt     = 'bin'

   gas_flux_ws%filename     = 'unknown'
   gas_flux_ws%file_varname = 'XKW'
   gas_flux_ws%scale_factor = c1
   gas_flux_ws%default_val  = c0
   gas_flux_ws%file_fmt     = 'bin'

   gas_flux_ap%filename     = 'unknown'
   gas_flux_ap%file_varname = 'P'
   gas_flux_ap%scale_factor = c1
   gas_flux_ap%default_val  = c0
   gas_flux_ap%file_fmt     = 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=sf6_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'sf6_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_sf6_option, master_task)
   call broadcast_scalar(init_sf6_init_file, master_task)
   call broadcast_scalar(init_sf6_init_file_fmt, master_task)

   do n = 1, sf6_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(psf6_file, master_task)
   call broadcast_scalar(psf6_first_nonzero_year, master_task)
   call broadcast_scalar(model_year, master_task)
   call broadcast_scalar(data_year, master_task)
   call broadcast_scalar(sf6_formulation, master_task)

   call broadcast_scalar(gas_flux_fice%filename, master_task)
   call broadcast_scalar(gas_flux_fice%file_varname, master_task)
   call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
   call broadcast_scalar(gas_flux_fice%default_val, master_task)
   call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

   fice_file%input = gas_flux_fice

   call broadcast_scalar(gas_flux_ws%filename, master_task)
   call broadcast_scalar(gas_flux_ws%file_varname, master_task)
   call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ws%default_val, master_task)
   call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

   xkw_file%input = gas_flux_ws

   call broadcast_scalar(gas_flux_ap%filename, master_task)
   call broadcast_scalar(gas_flux_ap%file_varname, master_task)
   call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
   call broadcast_scalar(gas_flux_ap%default_val, master_task)
   call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

   ap_file%input = gas_flux_ap

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_sf6_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d SF6s set to all zeros'
          write(stdout,delim_fmt)
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      ! if mapped_date is less than psf6_first_nonzero_year then
      ! register c0 as an io_read_fallback option
      mapped_date = iyear + (iday_of_year-1+frac_day)/days_in_year &
                    - model_year + data_year
      if (mapped_date < psf6_first_nonzero_year) then
         call io_read_fallback_register_tracer(tracername='SF6', &
            fallback_opt='const', const_val=c0)
      endif

      sf6_restart_filename = char_blank

      if (init_sf6_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read SF6s from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         sf6_restart_filename = read_restart_filename
         init_sf6_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         sf6_restart_filename = trim(init_sf6_init_file)

      endif

      call rest_read_tracer_block(sf6_ind_begin,          &
                                  init_sf6_init_file_fmt, &
                                  sf6_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

   case ('file')

      call document(subname, 'SF6s being read from separate file')

      call file_read_tracer_block(init_sf6_init_file_fmt, &
                                  init_sf6_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, sf6_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'sf6_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(subname, 'init_sf6_option', init_sf6_option)
      call exit_POP(sigAbort, 'unknown init_sf6_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, sf6_tracer_cnt
      do k = 1, km
         where (k > KMT(:,:,iblock))
            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
         end where
      end do
   end do
   end do

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
   LAND_MASK = (KMT.gt.0)

   call get_timer(sf6_sflux_timer, 'SF6_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call sf6_init_tavg
   call sf6_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine sf6_init

!***********************************************************************
!BOP
! !IROUTINE: sf6_init_tavg
! !INTERFACE:

 subroutine sf6_init_tavg

! !DESCRIPTION:
!  Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_SF6_IFRAC,'SF6_IFRAC',2,           &
                          long_name='Ice Fraction for SF6 fluxes',&
                          units='fraction', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SF6_XKW,'SF6_XKW',2,               &
                          long_name='XKW for SF6 fluxes',         &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SF6_ATM_PRESS,'SF6_ATM_PRESS',2,   &
                          long_name='Atmospheric Pressure for SF6 fluxes',&
                          units='atmospheres', grid_loc='2110',   &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pSF6,'pSF6',2,                 &
                          long_name='SF6 atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SF6_SCHMIDT,'SF6_SCHMIDT',2,   &
                          long_name='SF6 Schmidt Number',       &
                          units='none', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SF6_PV,'SF6_PV',2,             &
                          long_name='SF6 piston velocity',      &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SF6_surf_sat,'SF6_surf_sat',2, &
                          long_name='SF6 Saturation',           &
                          units='fmol/cm^3', grid_loc='2110',     &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(SF6_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   SF6_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine sf6_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: sf6_init_sflux
! !INTERFACE:

 subroutine sf6_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times
   use forcing_timeseries_mod, only: forcing_timeseries_init_dataset

! !DESCRIPTION:
!  Initialize surface flux computations for sf6 tracer module.
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'sf6_mod:sf6_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   call forcing_timeseries_init_dataset(psf6_file, &
      varnames      = (/ 'SF6NH', 'SF6SH' /), &
      model_year    = model_year, &
      data_year     = data_year, &
      taxmode_start = 'endpoint', &
      taxmode_end   = 'extrapolate', &
      dataset       = psf6_atm_forcing_dataset)

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!  otherwise, use values passed in
!-----------------------------------------------------------------------

   select case (sf6_formulation)

   case ('ocmip')

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
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

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
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

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

   case ('model')

      if (my_task == master_task) then
         write(stdout,*)  &
            ' Using fields from model forcing for calculating SF6 flux'
      endif

   case default
      call document(subname, 'sf6_formulation', sf6_formulation)

      call exit_POP(sigAbort, &
                    'sf6_init_sflux: Unknown value for sf6_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine sf6_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: sf6_set_sflux
! !INTERFACE:

 subroutine sf6_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute SF6 surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: field_loc_center, field_type_scalar, p5, xkw_coeff
   use time_management, only: thour00
   use forcing_tools, only: update_forcing_data, interpolate_forcing
   use timers, only: timer_start, timer_stop
   use forcing_timeseries_mod, only: forcing_timeseries_dataset_update_data

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      PRESS,     & ! sea level atmospheric pressure (dyne/cm**2)
      SST,       & ! sea surface temperature (C)
      SSS          ! sea surface salinity (psu)

   real (r8), dimension(nx_block,ny_block,sf6_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,sf6_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock             ! block index

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,      & ! used ice fraction (non-dimensional)
      XKW_USED,        & ! part of piston velocity (cm/s)
      AP_USED            ! used atm pressure (converted from dyne/cm**2 to atm)

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS,       & ! filtered surface tracer values
      pSF6,          & ! atmospheric SF6 mole fraction (pmol/mol)
      SF6_SCHMIDT,   & ! SF6 Schmidt number
      SF6_SOL_0,     & ! solubility of SF6 at 1 atm (mol/l/atm)
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,              & ! piston velocity (cm/s)
      SF6_surf_sat       ! SF6 surface saturation (fmol/cm^3)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

!-----------------------------------------------------------------------

   call timer_start(sf6_sflux_timer)

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
      AP_USED(:,:,iblock) = c0
   end do

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_update_data(psf6_atm_forcing_dataset)

   if (sf6_formulation == 'ocmip') then
       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(          fice_file%data_time,   &
               fice_file%data_time_min_loc,  fice_file%interp_type, &
               fice_file%data_next,          fice_file%data_update, &
               fice_file%data_type,          fice_file%data_inc,    &
               fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm, &
               tracer_data_label,            tracer_data_names,     &
               tracer_bndy_loc,              tracer_bndy_type,      &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:,:,:,:,1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(         xkw_file%data_time,   &
               xkw_file%data_time_min_loc,  xkw_file%interp_type, &
               xkw_file%data_next,          xkw_file%data_update, &
               xkw_file%data_type,          xkw_file%data_inc,    &
               xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm, &
               tracer_data_label,           tracer_data_names,    &
               tracer_bndy_loc,             tracer_bndy_type,     &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            xkw_file%DATA(:,:,:,:,1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(        ap_file%data_time,   &
               ap_file%data_time_min_loc,  ap_file%interp_type, &
               ap_file%data_next,          ap_file%data_update, &
               ap_file%data_type,          ap_file%data_inc,    &
               ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm, &
               tracer_data_label,          tracer_data_names,   &
               tracer_bndy_loc,            tracer_bndy_type,    &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:,:,:,:,1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = INTERP_WORK(:,:,:,1)
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,pSF6,SF6_SCHMIDT, &
   !$OMP                     SF6_SOL_0,XKW_ICE,&
   !$OMP                     PV,SF6_surf_sat)
   do iblock = 1, nblocks_clinic

      if (sf6_formulation == 'ocmip') then
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < 0.2000_r8) &
            IFRAC_USED(:,:,iblock) = 0.2000_r8
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > 0.9999_r8) &
            IFRAC_USED(:,:,iblock) = 0.9999_r8
      endif

      if (sf6_formulation == 'model') then
         where (LAND_MASK(:,:,iblock))
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)

            XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)

            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
         endwhere
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
            IFRAC_USED(:,:,iblock) = c0
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
            IFRAC_USED(:,:,iblock) = c1
      endif

!-----------------------------------------------------------------------
!  assume PRESS is in cgs units (dyne/cm**2) since that is what is
!    required for pressure forcing in barotropic
!  want units to be atmospheres
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) * (c1 / 1013.25e+3_r8)

      call comp_psf6(iblock, LAND_MASK(:,:,iblock), pSF6)

      call comp_sf6_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                            SF6_SCHMIDT)

      call comp_sf6_sol_0(LAND_MASK(:,:,iblock), SST(:,:,iblock), SSS(:,:,iblock), &
                          SF6_SOL_0)

      where (LAND_MASK(:,:,iblock))
         SF6_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         SF6_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         SF6_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
         SF6_SFLUX_TAVG(:,:,4,iblock) = pSF6
         SF6_SFLUX_TAVG(:,:,5,iblock) = SF6_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)
         PV = XKW_ICE * sqrt(660.0_r8 / SF6_SCHMIDT)
         SF6_surf_sat = AP_USED(:,:,iblock) * SF6_SOL_0 * pSF6
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,sf6_ind,iblock) + &
                         SURF_VALS_CUR(:,:,sf6_ind,iblock))
         STF_MODULE(:,:,sf6_ind,iblock) = &
            PV * (SF6_surf_sat - SURF_VALS)

         SF6_SFLUX_TAVG(:,:,6,iblock) = PV
         SF6_SFLUX_TAVG(:,:,7,iblock) = SF6_surf_sat

      elsewhere
         STF_MODULE(:,:,sf6_ind,iblock) = c0
      endwhere

   end do
   !$OMP END PARALLEL DO

   call timer_stop(sf6_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine sf6_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: comp_psf6
! !INTERFACE:

 subroutine comp_psf6(iblock, LAND_MASK, pSF6)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of SF6s
!  Linearly interpolate hemispheric values to current time step
!  Spatial pattern is determined by :
!     Northern Hemisphere value is used North of 10N
!     Southern Hemisphere value is used North of 10S
!     Linear Interpolation (in latitude) is used between 10N & 10S

! !REVISION HISTORY:
!  same as module

! !USES:

   use grid, only : TLATD
   use constants, only : c10
   use forcing_timeseries_mod, only: forcing_timeseries_dataset_get_var

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   integer (int_kind) :: &
      iblock          ! block index

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      pSF6  ! atmospheric SF6 mole fraction (pmol/mol)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j              ! loop indices

   real (r8) :: &
      psf6_nh_curr,   & ! psf6_nh for current time step (pmol/mol)
      psf6_sh_curr      ! psf6_sh for current time step (pmol/mol)

!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!
!  varind in the following calls must match varname ordering in 
!  call to forcing_timeseries_init_dataset in subroutine sf6_init_sflux
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_get_var(psf6_atm_forcing_dataset, varind=1, data_1d=psf6_nh_curr)
   call forcing_timeseries_dataset_get_var(psf6_atm_forcing_dataset, varind=2, data_1d=psf6_sh_curr)

!-----------------------------------------------------------------------
!     Merge hemisphere values.
!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            if (TLATD(i,j,iblock) < -c10) then
               pSF6(i,j) = psf6_sh_curr
            else if (TLATD(i,j,iblock) > c10) then
               pSF6(i,j) = psf6_nh_curr
            else
               pSF6(i,j) = psf6_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (psf6_nh_curr - psf6_sh_curr)
            endif
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_psf6

!***********************************************************************
!BOP
! !IROUTINE: comp_sf6_schmidt
! !INTERFACE:

 subroutine comp_sf6_schmidt(LAND_MASK, SST_IN, SF6_SCHMIDT)

! !DESCRIPTION:
!  Compute Schmidt numbers of SF6s.
!
!  range of validity of fit is -2:40
!
!  Ref : Wanninkhof 2014, Relationship between wind speed 
!        and gas exchange over the ocean revisited,
!        Limnol. Oceanogr.: Methods, 12, 
!        doi:10.4319/lom.2014.12.351
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in)  :: LAND_MASK(nx_block,ny_block)    ! land mask for this block
   real (r8)         , intent(in)  :: SST_IN(nx_block,ny_block)       ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8)         , intent(out) :: SF6_SCHMIDT(nx_block,ny_block)  ! Schmidt number of SF6 (non-dimensional)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   integer(int_kind)    :: i, j
   real (r8)            :: SST(nx_block,ny_block)

   real (r8), parameter :: a = 3177.5_r8
   real (r8), parameter :: b = -200.57_r8
   real (r8), parameter :: c =    6.8865_r8
   real (r8), parameter :: d =   -0.13335_r8
   real (r8), parameter :: e =    0.0010877_r8

!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            SST(i,j) = max(-2.0_r8, min(40.0_r8, SST_IN(i,j)))
            SF6_SCHMIDT(i,j) = a + SST(i,j) * (b + SST(i,j) * (c + SST(i,j) * (d + SST(i,j) * e)))
         else
            SF6_SCHMIDT(i,j) = c0
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_sf6_schmidt

!***********************************************************************
!BOP
! !IROUTINE: comp_sf6_sol_0
! !INTERFACE:

 subroutine comp_sf6_sol_0(LAND_MASK, SST, SSS, SF6_SOL_0)

! !DESCRIPTION:
!  Compute solubilities of SF6s at 1 atm.
!  Ref: Bullister et al., 2002: The solubility of sulfur 
!       hexafluoride in water and seawater, DSR, 49(1),
!       doi:10.1016/S0967-0637(01)00051-6.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: T0_Kelvin

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block) :: &
      SST,             & ! sea surface temperature (C)
      SSS                ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      SF6_SOL_0  ! solubility of SF6 at 1 atm (mol/l/atm)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a1 = -96.5975_r8,    &
      a2 = 139.883_r8,     &
      a3 =  37.8193_r8,    &
      a4 =   0.00000_r8,   &
      b1 =   0.0310693_r8, &
      b2 =  -0.0356385_r8, &
      b3 =   0.00743254_r8

   real (r8), dimension(nx_block,ny_block) :: &
      SSTKp01  ! .01 * sea surface temperature (in Kelvin)

!-----------------------------------------------------------------------

   SSTKp01 = merge( ((SST + T0_Kelvin)* 0.01_r8), c1, LAND_MASK)

   where (LAND_MASK)
      SF6_SOL_0 = EXP(a1 + a2 / SSTKp01 &
                        + a3 * LOG(SSTKp01) + a4 * SSTKp01 ** 2 &
                        + SSS * (b1 + SSTKp01 * (b2 + b3 * SSTKp01)))
   elsewhere
      SF6_SOL_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_sf6_sol_0

!***********************************************************************
!BOP
! !IROUTINE: sf6_tavg_forcing
! !INTERFACE:

 subroutine sf6_tavg_forcing

! !DESCRIPTION:
!  Make accumulation calls for forcing related tavg fields. This is
!  necessary because the forcing routines are called before tavg flags
!  are set.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,1,iblock),tavg_SF6_IFRAC,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,2,iblock),tavg_SF6_XKW,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,3,iblock),tavg_SF6_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,4,iblock),tavg_pSF6,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,6,iblock),tavg_SF6_SCHMIDT,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,6,iblock),tavg_SF6_PV,iblock,1)
         call accumulate_tavg_field(SF6_SFLUX_TAVG(:,:,7,iblock),tavg_SF6_surf_sat,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine sf6_tavg_forcing

!***********************************************************************

end module sf6_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
