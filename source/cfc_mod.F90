!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module cfc_mod

!BOP
! !MODULE: cfc_mod
!
!  Module for Chlorofluorocarbons (CFCs)
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
   use blocks, only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic, distrb_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants, only: c0, c1
   use io_types, only: stdout
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use passive_tracer_tools, only: forcing_monthly_every_ts, ind_name_pair
   use passive_tracer_tools, only : read_field, tracer_read
   use time_management, only : iyear, iday_of_year, frac_day, days_in_year
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
       cfc_tracer_cnt, &
       cfc_init, &
       cfc_set_sflux,  &
       cfc_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       cfc_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      cfc11_ind =  1,  & ! CFC11
      cfc12_ind =  2     ! CFC12

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(cfc_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(cfc11_ind, 'CFC11'), &
      ind_name_pair(cfc12_ind, 'CFC12') /)

!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) :: &
      cfc_formulation,        & ! how to calculate flux (ocmip or model)
      pcfc_file                 ! filename for netCDF time series of atm pcfc

   integer (int_kind) ::  &
      model_year,             & ! arbitrary model year
      data_year,              & ! year in data that corresponds to model_year
      pcfc_first_nonzero_year   ! first year of non-zero values in pcfc_file

   type (forcing_timeseries_dataset) :: &
      pcfc_atm_forcing_dataset  ! data structure for atm pcfc timeseries

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK               ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      fice_file,              & ! ice fraction, if read from file
      xkw_file,               & ! a * wind-speed ** 2, if read from file
      ap_file                   ! atmoshperic pressure, if read from file

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      CFC_SFLUX_TAVG

   integer (int_kind) :: &
      tavg_CFC_IFRAC,      & ! tavg id for ice fraction
      tavg_CFC_XKW,        & ! tavg id for xkw
      tavg_CFC_ATM_PRESS,  & ! tavg id for atmospheric pressure
      tavg_pCFC11,         & ! tavg id for cfc11 partial pressure
      tavg_pCFC12,         & ! tavg id for cfc12 partial pressure
      tavg_CFC11_SCHMIDT,  & ! tavg id for cfc11 Schmidt number
      tavg_CFC12_SCHMIDT,  & ! tavg id for cfc12 Schmidt number
      tavg_CFC11_PV,       & ! tavg id for cfc11 piston velocity
      tavg_CFC11_surf_sat, & ! tavg id for cfc11 surface saturation
      tavg_CFC12_PV,       & ! tavg id for cfc12 piston velocity
      tavg_CFC12_surf_sat    ! tavg id for cfc12 surface saturation

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: cfc_sflux_timer

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: cfc_init
! !INTERFACE:

 subroutine cfc_init(cfc_ind_begin, init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize cfc tracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: char_blank, delim_fmt
   use prognostic, only: curtime, oldtime, tracer_field
   use grid, only: KMT, n_topo_smooth, fill_points
   use io_types, only: nml_in, nml_filename
   use timers, only: get_timer
   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block
    use io_read_fallback_mod, only: io_read_fallback_register_tracer

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      cfc_ind_begin          ! starting index of cfc tracers in global tracer array
                             ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(cfc_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,cfc_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'cfc_mod:cfc_init'

   character(char_len) :: &
      init_cfc_option,        & ! option for initialization of bgc
      init_cfc_init_file,     & ! filename for option 'file'
      init_cfc_init_file_fmt    ! file format for option 'file'

   integer (int_kind) :: &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(cfc_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   type(tracer_read) :: &
      gas_flux_fice,          & ! ice fraction for gas fluxes
      gas_flux_ws,            & ! wind speed for gas fluxes
      gas_flux_ap               ! atmospheric pressure for gas fluxes

   namelist /cfc_nml/ &
      init_cfc_option, init_cfc_init_file, init_cfc_init_file_fmt, &
      tracer_init_ext, pcfc_file, pcfc_first_nonzero_year, model_year, data_year, &
      cfc_formulation, gas_flux_fice, gas_flux_ws, gas_flux_ap

   real (r8) :: &
      mapped_date               ! date of current model timestep mapped to data timeline

   character (char_len) ::  &
      cfc_restart_filename      ! modified file name for restart file

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

   do n = 1, cfc_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'fmol/cm^3'
      tracer_d_module(n)%tend_units = 'fmol/cm^3/s'
      tracer_d_module(n)%flux_units = 'fmol/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_cfc_option        = 'unknown'
   init_cfc_init_file     = 'unknown'
   init_cfc_init_file_fmt = 'bin'

   do n = 1, cfc_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   pcfc_file               = 'unknown'
   pcfc_first_nonzero_year = 1936
   model_year              = 1
   data_year               = 1
   cfc_formulation         = 'model'

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
         read(nml_in, nml=cfc_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'cfc_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_cfc_option, master_task)
   call broadcast_scalar(init_cfc_init_file, master_task)
   call broadcast_scalar(init_cfc_init_file_fmt, master_task)

   do n = 1, cfc_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(pcfc_file, master_task)
   call broadcast_scalar(pcfc_first_nonzero_year, master_task)
   call broadcast_scalar(model_year, master_task)
   call broadcast_scalar(data_year, master_task)
   call broadcast_scalar(cfc_formulation, master_task)

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

   select case (init_cfc_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d CFCs set to all zeros'
          write(stdout,delim_fmt)
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      ! if mapped_date is less than pcfc_first_nonzero_year then
      ! register c0 as an io_read_fallback option
      mapped_date = iyear + (iday_of_year-1+frac_day)/days_in_year &
                    - model_year + data_year
      if (mapped_date < pcfc_first_nonzero_year) then
         call io_read_fallback_register_tracer(tracername='CFC11', &
            fallback_opt='const', const_val=c0)
         call io_read_fallback_register_tracer(tracername='CFC12', &
            fallback_opt='const', const_val=c0)
      endif

      cfc_restart_filename = char_blank

      if (init_cfc_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read CFCs from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         cfc_restart_filename = read_restart_filename
         init_cfc_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         cfc_restart_filename = trim(init_cfc_init_file)

      endif

      call rest_read_tracer_block(cfc_ind_begin, &
                                  init_cfc_init_file_fmt, &
                                  cfc_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

   case ('file')

      call document(subname, 'CFCs being read from separate file')

      call file_read_tracer_block(init_cfc_init_file_fmt, &
                                  init_cfc_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, cfc_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'cfc_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(subname, 'init_cfc_option', init_cfc_option)
      call exit_POP(sigAbort, 'unknown init_cfc_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, cfc_tracer_cnt
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

   call get_timer(cfc_sflux_timer, 'CFC_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call cfc_init_tavg
   call cfc_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init

!***********************************************************************
!BOP
! !IROUTINE: cfc_init_tavg
! !INTERFACE:

 subroutine cfc_init_tavg

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

   call define_tavg_field(tavg_CFC_IFRAC,'CFC_IFRAC',2,           &
                          long_name='Ice Fraction for CFC fluxes',&
                          units='fraction', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC_XKW,'CFC_XKW',2,               &
                          long_name='XKW for CFC fluxes',         &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC_ATM_PRESS,'CFC_ATM_PRESS',2,   &
                          long_name='Atmospheric Pressure for CFC fluxes',&
                          units='atmospheres', grid_loc='2110',   &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCFC11,'pCFC11',2,                 &
                          long_name='CFC11 atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCFC12,'pCFC12',2,                 &
                          long_name='CFC12 atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_SCHMIDT,'CFC11_SCHMIDT',2,   &
                          long_name='CFC11 Schmidt Number',       &
                          units='none', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_SCHMIDT,'CFC12_SCHMIDT',2,   &
                          long_name='CFC12 Schmidt Number',       &
                          units='none', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_PV,'CFC11_PV',2,             &
                          long_name='CFC11 piston velocity',      &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_surf_sat,'CFC11_surf_sat',2, &
                          long_name='CFC11 Saturation',           &
                          units='fmol/cm^3', grid_loc='2110',     &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_PV,'CFC12_PV',2,             &
                          long_name='CFC12 piston velocity',      &
                          units='cm/s', grid_loc='2110',          &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_surf_sat,'CFC12_surf_sat',2, &
                          long_name='CFC12 Saturation',           &
                          units='fmol/cm^3', grid_loc='2110',     &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(CFC_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   CFC_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: cfc_init_sflux
! !INTERFACE:

 subroutine cfc_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times
   use forcing_timeseries_mod, only: forcing_timeseries_init_dataset

! !DESCRIPTION:
!  Initialize surface flux computations for cfc tracer module.
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'cfc_mod:cfc_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   call forcing_timeseries_init_dataset(pcfc_file, &
      varnames      = (/ 'pcfc11_nh', 'pcfc11_sh', 'pcfc12_nh', 'pcfc12_sh' /), &
      model_year    = model_year, &
      data_year     = data_year, &
      taxmode_start = 'endpoint', &
      taxmode_end   = 'extrapolate', &
      dataset       = pcfc_atm_forcing_dataset)

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!  otherwise, use values passed in
!-----------------------------------------------------------------------

   select case (cfc_formulation)

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
            ' Using fields from model forcing for calculating CFC flux'
      endif

   case default
      call document(subname, 'cfc_formulation', cfc_formulation)

      call exit_POP(sigAbort, &
                    'cfc_init_sflux: Unknown value for cfc_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: cfc_set_sflux
! !INTERFACE:

 subroutine cfc_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute CFC11 surface flux and store related tavg fields for
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

   real (r8), dimension(nx_block,ny_block,cfc_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,cfc_tracer_cnt,max_blocks_clinic), &
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
      pCFC11,          & ! atmospheric CFC11 mole fraction (pmol/mol)
      pCFC12,          & ! atmospheric CFC11 mole fraction (pmol/mol)
      CFC11_SCHMIDT,   & ! CFC11 Schmidt number
      CFC12_SCHMIDT,   & ! CFC12 Schmidt number
      CFC11_SOL_0,     & ! solubility of CFC11 at 1 atm (mol/l/atm)
      CFC12_SOL_0,     & ! solubility of CFC12 at 1 atm (mol/l/atm)
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,              & ! piston velocity (cm/s)
      CFC_surf_sat       ! CFC surface saturation (either CFC11 or CFC12) (fmol/cm^3)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

!-----------------------------------------------------------------------

   call timer_start(cfc_sflux_timer)

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
      AP_USED(:,:,iblock) = c0
   end do

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_update_data(pcfc_atm_forcing_dataset)

   if (cfc_formulation == 'ocmip') then
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

   !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,pCFC11,pCFC12,CFC11_SCHMIDT, &
   !$OMP                     CFC12_SCHMIDT,CFC11_SOL_0,CFC12_SOL_0,XKW_ICE,&
   !$OMP                     PV,CFC_surf_sat)
   do iblock = 1, nblocks_clinic

      if (cfc_formulation == 'ocmip') then
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < 0.2000_r8) &
            IFRAC_USED(:,:,iblock) = 0.2000_r8
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > 0.9999_r8) &
            IFRAC_USED(:,:,iblock) = 0.9999_r8
      endif

      if (cfc_formulation == 'model') then
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

      call comp_pcfc(iblock, LAND_MASK(:,:,iblock), pCFC11, pCFC12)

      call comp_cfc_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                            CFC11_SCHMIDT, CFC12_SCHMIDT)

      call comp_cfc_sol_0(LAND_MASK(:,:,iblock), SST(:,:,iblock), SSS(:,:,iblock), &
                          CFC11_SOL_0, CFC12_SOL_0)

      where (LAND_MASK(:,:,iblock))
         CFC_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,4,iblock) = pCFC11
         CFC_SFLUX_TAVG(:,:,5,iblock) = pCFC12
         CFC_SFLUX_TAVG(:,:,6,iblock) = CFC11_SCHMIDT
         CFC_SFLUX_TAVG(:,:,7,iblock) = CFC12_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)

         PV = XKW_ICE * sqrt(660.0_r8 / CFC11_SCHMIDT)
         CFC_SFLUX_TAVG(:,:,8,iblock) = PV
         CFC_surf_sat = AP_USED(:,:,iblock) * CFC11_SOL_0 * pCFC11
         CFC_SFLUX_TAVG(:,:,9,iblock) = CFC_surf_sat
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,cfc11_ind,iblock) + &
                         SURF_VALS_CUR(:,:,cfc11_ind,iblock))
         STF_MODULE(:,:,cfc11_ind,iblock) = &
            PV * (CFC_surf_sat - SURF_VALS)

         PV = XKW_ICE * sqrt(660.0_r8 / CFC12_SCHMIDT)
         CFC_SFLUX_TAVG(:,:,10,iblock) = PV
         CFC_surf_sat = AP_USED(:,:,iblock) * CFC12_SOL_0 * pCFC12
         CFC_SFLUX_TAVG(:,:,11,iblock) = CFC_surf_sat
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,cfc12_ind,iblock) + &
                         SURF_VALS_CUR(:,:,cfc12_ind,iblock))
         STF_MODULE(:,:,cfc12_ind,iblock) = &
            PV * (CFC_surf_sat - SURF_VALS)
      elsewhere
         STF_MODULE(:,:,cfc11_ind,iblock) = c0
         STF_MODULE(:,:,cfc12_ind,iblock) = c0
      endwhere

   end do
   !$OMP END PARALLEL DO

   call timer_stop(cfc_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: comp_pcfc
! !INTERFACE:

 subroutine comp_pcfc(iblock, LAND_MASK, pCFC11, pCFC12)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of CFCs
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
      pCFC11,       & ! atmospheric CFC11 mole fraction (pmol/mol)
      pCFC12          ! atmospheric CFC11 mole fraction (pmol/mol)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j              ! loop indices

   real (r8) :: &
      pcfc11_nh_curr, & ! pcfc11_nh for current time step (pmol/mol)
      pcfc11_sh_curr, & ! pcfc11_sh for current time step (pmol/mol)
      pcfc12_nh_curr, & ! pcfc12_nh for current time step (pmol/mol)
      pcfc12_sh_curr    ! pcfc12_sh for current time step (pmol/mol)

!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!
!  varind in the following calls must match varname ordering in 
!  call to forcing_timeseries_init_dataset in subroutine cfc_init_sflux
!-----------------------------------------------------------------------

   call forcing_timeseries_dataset_get_var(pcfc_atm_forcing_dataset, varind=1, data_1d=pcfc11_nh_curr)
   call forcing_timeseries_dataset_get_var(pcfc_atm_forcing_dataset, varind=2, data_1d=pcfc11_sh_curr)
   call forcing_timeseries_dataset_get_var(pcfc_atm_forcing_dataset, varind=3, data_1d=pcfc12_nh_curr)
   call forcing_timeseries_dataset_get_var(pcfc_atm_forcing_dataset, varind=4, data_1d=pcfc12_sh_curr)

!-----------------------------------------------------------------------
!     Merge hemisphere values.
!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            if (TLATD(i,j,iblock) < -c10) then
               pCFC11(i,j) = pcfc11_sh_curr
               pCFC12(i,j) = pcfc12_sh_curr
            else if (TLATD(i,j,iblock) > c10) then
               pCFC11(i,j) = pcfc11_nh_curr
               pCFC12(i,j) = pcfc12_nh_curr
            else
               pCFC11(i,j) = pcfc11_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (pcfc11_nh_curr - pcfc11_sh_curr)
               pCFC12(i,j) = pcfc12_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (pcfc12_nh_curr - pcfc12_sh_curr)
            endif
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_pcfc

!***********************************************************************
!BOP
! !IROUTINE: comp_cfc_schmidt
! !INTERFACE:

 subroutine comp_cfc_schmidt(LAND_MASK, SST_IN, CFC11_SCHMIDT, CFC12_SCHMIDT)

! !DESCRIPTION:
!  Compute Schmidt numbers of CFCs.
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

   logical (log_kind), intent(in)  :: LAND_MASK(nx_block,ny_block)      ! land mask for this block
   real (r8)         , intent(in)  :: SST_IN(nx_block,ny_block)         ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8)         , intent(out) :: CFC11_SCHMIDT(nx_block,ny_block)  ! Schmidt number of CFC11 (non-dimensional)
   real (r8)         , intent(out) :: CFC12_SCHMIDT(nx_block,ny_block)  ! Schmidt number of CFC12 (non-dimensional)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   integer(int_kind)    :: i, j
   real (r8)            :: SST(nx_block,ny_block)

   real (r8), parameter :: a_11 = 3579.2_r8
   real (r8), parameter :: b_11 = -222.63_r8
   real (r8), parameter :: c_11 =    7.5749_r8
   real (r8), parameter :: d_11 =   -0.14595_r8
   real (r8), parameter :: e_11 =    0.0011874_r8

   real (r8), parameter :: a_12 = 3828.1_r8
   real (r8), parameter :: b_12 = -249.86_r8
   real (r8), parameter :: c_12 =    8.7603_r8
   real (r8), parameter :: d_12 =   -0.1716_r8
   real (r8), parameter :: e_12 =    0.001408_r8

!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            SST(i,j) = max(-2.0_r8, min(40.0_r8, SST_IN(i,j)))
            CFC11_SCHMIDT(i,j) = a_11 + SST(i,j) * (b_11 + SST(i,j) * (c_11 + SST(i,j) * (d_11 + SST(i,j) * e_11)))
            CFC12_SCHMIDT(i,j) = a_12 + SST(i,j) * (b_12 + SST(i,j) * (c_12 + SST(i,j) * (d_12 + SST(i,j) * e_12)))
         else
            CFC11_SCHMIDT(i,j) = c0
            CFC12_SCHMIDT(i,j) = c0
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_cfc_schmidt

!***********************************************************************
!BOP
! !IROUTINE: comp_cfc_sol_0
! !INTERFACE:

 subroutine comp_cfc_sol_0(LAND_MASK, SST, SSS, CFC11_SOL_0, CFC12_SOL_0)

! !DESCRIPTION:
!  Compute solubilities of CFCs at 1 atm.
!  Ref : Warner & Weiss (1985), Deep Sea Reasearch,
!        Vol 32, No. 12, pp. 1485-1497 (Table 5)
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
      CFC11_SOL_0,     & ! solubility of CFC11 at 1 atm (mol/l/atm)
      CFC12_SOL_0        ! solubility of CFC12 at 1 atm (mol/l/atm)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a1_11 = -229.9261_r8,     a1_12 = -218.0971_r8, &
      a2_11 =  319.6552_r8,     a2_12 =  298.9702_r8, &
      a3_11 =  119.4471_r8,     a3_12 =  113.8049_r8, &
      a4_11 =   -1.39165_r8,    a4_12 =   -1.39165_r8, &
      b1_11 =   -0.142382_r8,   b1_12 =   -0.143566_r8, &
      b2_11 =    0.091459_r8,   b2_12 =    0.091015_r8, &
      b3_11 =   -0.0157274_r8,  b3_12 =   -0.0153924_r8

   real (r8), dimension(nx_block,ny_block) :: &
      SSTKp01  ! .01 * sea surface temperature (in Kelvin)

!-----------------------------------------------------------------------

   SSTKp01 = merge( ((SST + T0_Kelvin)* 0.01_r8), c1, LAND_MASK)

   where (LAND_MASK)
      CFC11_SOL_0 = EXP(a1_11 + a2_11 / SSTKp01 &
                        + a3_11 * LOG(SSTKp01) + a4_11 * SSTKp01 ** 2 &
                        + SSS * (b1_11 + SSTKp01 * (b2_11 + b3_11 * SSTKp01)))
      CFC12_SOL_0 = EXP(a1_12 + a2_12 / SSTKp01 &
                        + a3_12 * LOG(SSTKp01) + a4_12 * SSTKp01 ** 2 &
                        + SSS * (b1_12 + SSTKp01 * (b2_12 + b3_12 * SSTKp01)))
   elsewhere
      CFC11_SOL_0 = c0
      CFC12_SOL_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_cfc_sol_0

!***********************************************************************
!BOP
! !IROUTINE: cfc_tavg_forcing
! !INTERFACE:

 subroutine cfc_tavg_forcing

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
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,1,iblock),tavg_CFC_IFRAC,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,2,iblock),tavg_CFC_XKW,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,3,iblock),tavg_CFC_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,4,iblock),tavg_pCFC11,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,5,iblock),tavg_pCFC12,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,6,iblock),tavg_CFC11_SCHMIDT,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,7,iblock),tavg_CFC12_SCHMIDT,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,8,iblock),tavg_CFC11_PV,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,9,iblock),tavg_CFC11_surf_sat,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,10,iblock),tavg_CFC12_PV,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,11,iblock),tavg_CFC12_surf_sat,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_tavg_forcing

!***********************************************************************

end module cfc_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
