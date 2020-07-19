!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module abio_dic_dic14_mod

!BOP
! !MODULE: abio_dic_dic14_mod
!
!  Module for Abiotic DIC & DIC14, essentially following OCMIP2 protocols
!  but using ice fracetion and winds from the model/data atmos instead of from files
!
!  There is an option to read CO2 and D14C from files, use the atmospheric CO2
!  from the coupler and the D14C from the file, or to use constant CO2 and D14C values
!
!  The units of concentration for these tracers are
!     nmol/cm^3 == mmol/m^3 == umol/l ~= umol/kg.
!
!  The units of surface fluxes for these tracers are
!     nmol/cm^3 * cm/s == nmol/cm^2/s.
!
! Written by Alexandra Jahn, NCAR
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:


   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size, only: max_blocks_clinic, km, nx_global, ny_global
   use domain, only: nblocks_clinic, distrb_clinic, blocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants, only: c0, c1, c1000, p5, c2,rho_sw
   use io_types
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use passive_tracer_tools, only: ind_name_pair, read_field, comp_surf_avg
   use passive_tracer_tools, only: extract_surf_avg, tracer_read
   use broadcast, only: broadcast_array, broadcast_scalar
   use netcdf
   use co2calc
   use time_management
   use constants, only: xkw_coeff
   use schmidt_number, only: SCHMIDT_CO2

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
      abio_dic_dic14_tracer_cnt,     &
      abio_dic_dic14_init,           &
      abio_dic_dic14_tracer_ref_val, &
      abio_dic_dic14_set_sflux,      &
      abio_dic_dic14_tavg_forcing,   &
      abio_dic_dic14_set_interior,   &
      abio_dic_dic14_write_restart


!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      abio_dic_dic14_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      abio_dic_ind   =  1,  & ! Abiotic DIC
      abio_dic14_ind =  2     ! Abiotic DIC14

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(abio_dic_dic14_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(abio_dic_ind, 'ABIO_DIC'), &
      ind_name_pair(abio_dic14_ind, 'ABIO_DIC14') /)

!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) ::    &
      abio_atm_co2_opt              ! option for CO2 varying or constant forcing

   real (r8) :: &
      abio_atm_co2_const            !  atmospheric CO2 constant [ppm]

   character(char_len) ::    &
      abio_atm_co2_filename         ! filename for varying atm CO2 data

   integer (int_kind) :: &
      abio_atm_co2_data_nbval       !  number of values in abio_atm_co2_filename

   real (r8), dimension(:), allocatable :: &
      abio_atm_co2_data_ppm,      & !  atmospheric pCO2 values in datafile [ppm]
      abio_atm_co2_data_yr          !  date of atmospheric pCO2 values in datafile

   character(char_len) :: &
      abio_atm_d14c_opt             ! option for D14C varying or constant forcing

   real (r8) :: &
      abio_atm_d14c_const           !  atmospheric 14CO2 constant [permil]

   real (r8), dimension(3) :: &
      abio_atm_d14c_lat_band_vals   !  atmospheric 14CO2 constant [permil]

   character(char_len) :: &
      abio_atm_d14c_filename        ! filename for varying atm D14C

   integer (int_kind) :: &
      abio_atm_model_year,        & !  arbitrary model year
      abio_atm_data_year            !  year in data that corresponds to abio_atm_model_year

!-------------------------------------------------------------------------
!  named field indices
!-------------------------------------------------------------------------
   integer (int_kind) :: &
      atm_co2_nf_ind     = 0      ! atmospheric co2

!-----------------------------------------------------------------------
!  module variables related to ph computations
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:), allocatable, target :: &
      ABIO_PH_PREV                     ! computed ph from previous time step

!-----------------------------------------------------------------------
!  define array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      ABIO_DIC_SFLUX_TAVG

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_ABIO_CO2_IFRAC,      & ! tavg id for ice fraction
      tavg_ABIO_CO2_XKW,        & ! tavg id for xkw
      tavg_ABIO_CO2_ATM_PRESS,  & ! tavg id for atmospheric pressure
      tavg_ABIO_pCO2,           & ! tavg id for atmospheric CO2 partial pressure
      tavg_ABIO_D14Catm,        & ! tavg id for atmospheric D14C in permil
      tavg_ABIO_CO2_SCHMIDT,    & ! tavg id for CO2 Schmidt number
      tavg_ABIO_CO2_PV,         & ! tavg id for CO2 piston velocity
      tavg_ABIO_pCO2SURF,       & ! tavg id for surface pco2
      tavg_ABIO_DCO2STAR,       & ! tavg id for surface DCO2STAR (CO2STAR_AIR-CO2STAR)
      tavg_ABIO_CO2STAR,        & ! tavg id for surface CO2Star
      tavg_ABIO_DpCO2,          & ! tavg id for delta pco2
      tavg_FG_ABIO_DIC14,       & ! tavg id for surface gas flux of C14
      tavg_FG_ABIO_DIC,         & ! tavg id for surface gas flux of CO2
      tavg_ABIO_ALK,            & ! tavg id for surface Alkalinity
      tavg_ABIO_PH,             & ! tavg id for surface PH
      tavg_d_SF_ABIO_DIC_d_ABIO_DIC,    & ! tavg id for derivative of SF_ABIO_DIC wrt ABIO_DIC
      tavg_d_SF_ABIO_DIC14_d_ABIO_DIC,  & ! tavg id for derivative of SF_ABIO_DIC14 wrt ABIO_DIC
      tavg_d_SF_ABIO_DIC14_d_ABIO_DIC14   ! tavg id for derivative of SF_ABIO_DIC14 wrt ABIO_DIC14

!-----------------------------------------------------------------------
!  define tavg id for 3d fields related to surface fluxes
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_ABIO_D14Cocn           ! tavg id for (DIC14/DIC -1 )* 1000

!-----------------------------------------------------------------------
!  data_ind_co2 is the index for the CO2 data for the current timestep
!  Note that data_ind_co2 is always less than abio_atm_co2_data_nbval
!  To enable OpenMP parallelism, duplicating data_ind for each block
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:), allocatable :: &
     data_ind_co2                 ! data index for CO2 data

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: abio_dic_dic14_sflux_timer

!-----------------------------------------------------------------------
!  scalar constants for decay calculation
!-----------------------------------------------------------------------

   real (r8), parameter :: c14_halflife_years = 5730.0_r8 !C14 half file
   real (r8) :: c14_lambda_inv_sec           ! Decay variable in seconds

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   logical (log_kind), dimension(abio_dic_dic14_tracer_cnt) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      comp_surf_avg_flag        ! time flag id for computing average
                                ! surface tracer values

   real (r8), dimension(abio_dic_dic14_tracer_cnt) :: &
      surf_avg                  ! average surface tracer values

!-----------------------------------------------------------------------
! Values for co2_calc
!-----------------------------------------------------------------------

 real (r8), parameter :: &
      abio_phlo_surf_init = 7.0_r8, &      ! low bound for surface ph for no prev soln
      abio_phhi_surf_init = 9.0_r8, &      ! high bound for surface ph for no prev soln
      abio_del_ph         = 0.20_r8        ! delta-ph for prev soln


!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_init
! !INTERFACE:


 subroutine abio_dic_dic14_init(abio_dic_dic14_ind_begin, init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize Abiotic DIC DIC14 tracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:
   use constants, only: char_blank
   use prognostic, only: curtime, oldtime, tracer_field
   use grid, only: KMT, n_topo_smooth, fill_points
   use timers, only: get_timer
   use passive_tracer_tools, only: field_exists_in_file
   use passive_tracer_tools, only: rest_read_tracer_block, file_read_tracer_block
   use io_read_fallback_mod, only: io_read_fallback_register_field


! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      abio_dic_dic14_ind_begin ! starting index of abio_dic_dic14 tracers in global tracer array
                               ! passed through to rest_read_tracer_block

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(abio_dic_dic14_tracer_cnt), intent(inout) :: &
      tracer_d_module        ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,abio_dic_dic14_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'abio_dic_dic14_mod:abio_dic_dic14_init'

   character(char_len) :: &
      init_abio_dic_dic14_option,        &  ! option for initialization of abio dic/dic14
      init_abio_dic_dic14_init_file,     &  ! filename for option 'file'
      init_abio_dic_dic14_init_file_fmt, &  ! file format for option 'file'
      abio_dic_dic14_restfile_fallback,  &  ! restfile to read tracers from, if they are not found in init_file
      fieldname,                         &  ! restfile field name
      abio_comp_surf_avg_freq_opt,       &  ! choice for freq of comp_surf_avg
      abio_dic_dic14_restart_filename       ! modified file name for restart file

   integer (int_kind) :: &
      n,                                  & ! index for looping over tracers
      k,                                  & ! index for looping over depth levels
      iblock,                             & ! index for looping over blocks
      nml_error                             ! namelist i/o error flag

   type(tracer_read), dimension(abio_dic_dic14_tracer_cnt) :: &
      abio_tracer_init_ext                  ! namelist variable for initializing tracers


   integer (int_kind) :: &
      abio_comp_surf_avg_freq_iopt,       & ! choice for freq of comp_surf_avg
      abio_comp_surf_avg_freq               ! choice for freq of comp_surf_avg

   real (r8) :: &
      abio_surf_avg_dic_const,            & ! Constant surface DIC
      abio_surf_avg_dic14_const             ! Constant surface DIC14

   logical (log_kind) :: &
      abio_use_nml_surf_vals,             & ! do namelist surf values override values
                                            ! from restart file?
      all_fields_exist_in_restfile          ! are all needed fields present in the restart file?

   namelist /abio_dic_dic14_nml/ &
      init_abio_dic_dic14_option, init_abio_dic_dic14_init_file, abio_tracer_init_ext, &
      init_abio_dic_dic14_init_file_fmt, abio_dic_dic14_restfile_fallback, &
      abio_surf_avg_dic_const, abio_use_nml_surf_vals, &
      abio_comp_surf_avg_freq_opt, abio_comp_surf_avg_freq, abio_surf_avg_dic14_const, &
      abio_atm_co2_opt, abio_atm_co2_const, abio_atm_co2_filename, &
      abio_atm_d14c_opt, abio_atm_d14c_const, abio_atm_d14c_lat_band_vals, &
      abio_atm_d14c_filename, abio_atm_model_year, abio_atm_data_year

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   do n = 1, abio_dic_dic14_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'nmol/cm^3'
      tracer_d_module(n)%tend_units = 'nmol/cm^3/s'
      tracer_d_module(n)%flux_units = 'nmol/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_abio_dic_dic14_option         = 'unknown'
   init_abio_dic_dic14_init_file      = 'unknown'
   init_abio_dic_dic14_init_file_fmt  = 'bin'
   abio_dic_dic14_restfile_fallback   = 'unknown'

   do n = 1, abio_dic_dic14_tracer_cnt
      abio_tracer_init_ext(n)%mod_varname  = 'unknown'
      abio_tracer_init_ext(n)%filename     = 'unknown'
      abio_tracer_init_ext(n)%file_varname = 'unknown'
      abio_tracer_init_ext(n)%scale_factor = c1
      abio_tracer_init_ext(n)%default_val  = c0
      abio_tracer_init_ext(n)%file_fmt     = 'bin'
   end do


   abio_use_nml_surf_vals            = .false.
   abio_surf_avg_dic_const           = 1944._r8
   abio_surf_avg_dic14_const         = 1944._r8
   abio_comp_surf_avg_freq_opt       = 'never'
   abio_comp_surf_avg_freq           = 1

   abio_atm_co2_opt                  = 'const'
   abio_atm_co2_const                = 280.0_r8
   abio_atm_co2_filename             = 'unknown'

   abio_atm_d14c_opt                 = 'lat_bands'
   abio_atm_d14c_const               = c0
   abio_atm_d14c_lat_band_vals(:)    = (/ -2.3_r8, -4.0_r8, -5.8_r8 /)
   abio_atm_d14c_filename            = 'unknown'

   abio_atm_model_year               = 1
   abio_atm_data_year                = 1

!-----------------------------------------------------------------------
!  read namelist settings from namelist
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=abio_dic_dic14_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'abio_dic_dic14_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' abio_dic_dic14:'
      write(stdout,blank_fmt)
      write(stdout,*) ' abio_dic_dic14_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,abio_dic_dic14_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_abio_dic_dic14_option, master_task)
   call broadcast_scalar(init_abio_dic_dic14_init_file, master_task)
   call broadcast_scalar(init_abio_dic_dic14_init_file_fmt, master_task)
   call broadcast_scalar(abio_dic_dic14_restfile_fallback, master_task)

   do n = 1, abio_dic_dic14_tracer_cnt
      call broadcast_scalar(abio_tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(abio_tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(abio_tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(abio_tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(abio_tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(abio_tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(abio_use_nml_surf_vals, master_task)
   call broadcast_scalar(abio_surf_avg_dic_const, master_task)
   call broadcast_scalar(abio_surf_avg_dic14_const, master_task)
   call broadcast_scalar(abio_comp_surf_avg_freq_opt, master_task)
   call broadcast_scalar(abio_comp_surf_avg_freq, master_task)

   call broadcast_scalar(abio_atm_co2_opt, master_task)
   call broadcast_scalar(abio_atm_co2_const, master_task)
   call broadcast_scalar(abio_atm_co2_filename, master_task)

   call broadcast_scalar(abio_atm_d14c_opt, master_task)
   call broadcast_scalar(abio_atm_d14c_const, master_task)
   call broadcast_array(abio_atm_d14c_lat_band_vals, master_task)
   call broadcast_scalar(abio_atm_d14c_filename, master_task)

   call broadcast_scalar(abio_atm_model_year, master_task)
   call broadcast_scalar(abio_atm_data_year, master_task)

!-----------------------------------------------------------------------
!  set variables immediately dependent on namelist variables
!-----------------------------------------------------------------------

   select case (abio_comp_surf_avg_freq_opt)
   case ('never')
      abio_comp_surf_avg_freq_iopt = freq_opt_never
   case ('nyear')
      abio_comp_surf_avg_freq_iopt = freq_opt_nyear
   case ('nmonth')
      abio_comp_surf_avg_freq_iopt = freq_opt_nmonth
   case default
      call document(subname, 'abio_comp_surf_avg_freq_opt', abio_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'unknown abio_comp_surf_avg_freq_opt')
   end select

   call init_time_flag('abio_dic_dic14_comp_surf_avg', comp_surf_avg_flag, &
     default=.false., freq_opt=abio_comp_surf_avg_freq_iopt,  &
     freq=abio_comp_surf_avg_freq, owner='abio_dic_dic14_init')

!-----------------------------------------------------------------------
!  namelist consistency checking
!-----------------------------------------------------------------------

   if (abio_use_nml_surf_vals .and. abio_comp_surf_avg_freq_iopt /= freq_opt_never) then
      call document(subname, 'abio_use_nml_surf_vals', abio_use_nml_surf_vals)
      call document(subname, 'abio_comp_surf_avg_freq_opt', abio_comp_surf_avg_freq_opt)
      call exit_POP(sigAbort, 'abio_use_nml_surf_vals can only be .true. if ' /&
                           &/ ' abio_comp_surf_avg_freq_opt is never')
   endif

!-----------------------------------------------------------------------
!  initialize virtual flux flag array
!-----------------------------------------------------------------------

   vflux_flag = .false.
   vflux_flag(abio_dic_ind) = .true.
   vflux_flag(abio_dic14_ind) = .true.

!-----------------------------------------------------------------------
!  Initialize PH
!-----------------------------------------------------------------------

   allocate( ABIO_PH_PREV(nx_block,ny_block,max_blocks_clinic) )

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_abio_dic_dic14_option)


   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      abio_dic_dic14_restart_filename = char_blank

      if (init_abio_dic_dic14_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read Abiotic DIC & DIC14 from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         abio_dic_dic14_restart_filename = read_restart_filename
         init_abio_dic_dic14_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         abio_dic_dic14_restart_filename = trim(init_abio_dic_dic14_init_file)

      endif

      all_fields_exist_in_restfile = .true.
      do n = 1, abio_dic_dic14_tracer_cnt
         fieldname = trim(tracer_d_module(n)%short_name) /&
              &/ '_CUR'
         if (.not. field_exists_in_file(init_abio_dic_dic14_init_file_fmt, &
              abio_dic_dic14_restart_filename, fieldname)) then
            call document(subname, trim(fieldname) /&
                 &/ ' not found in ' /&
                 &/ trim(abio_dic_dic14_restart_filename))
            all_fields_exist_in_restfile = .false.
         endif

         fieldname = trim(tracer_d_module(n)%short_name) /&
              &/ '_OLD'
         if (.not. field_exists_in_file(init_abio_dic_dic14_init_file_fmt, &
              abio_dic_dic14_restart_filename, fieldname)) then
            call document(subname, trim(fieldname) /&
                 &/ ' not found in ' /&
                 &/ trim(abio_dic_dic14_restart_filename))
            all_fields_exist_in_restfile = .false.
         endif
      end do

      if (.not. all_fields_exist_in_restfile) then
         if (trim(abio_dic_dic14_restfile_fallback) == 'unknown') then
            call document(subname, 'some abio_dic_dic14 fields missing from restfile, ' /&
                 &/ 'but abio_dic_dic14_restfile_fallback == unknown')
            call exit_POP(sigAbort, 'stopping in ' /&
                 &/ subname)
         endif
         call document(subname, 'some abio_dic_dic14 fields missing from restfile, ' /&
              &/ 'using abio_dic_dic14_restfile_fallback')
         abio_dic_dic14_restart_filename = abio_dic_dic14_restfile_fallback
      endif

      call io_read_fallback_register_field('ABIO_PH_SURF', &
           fallback_opt='const', const_val=c0)

      call rest_read_tracer_block(abio_dic_dic14_ind_begin,          &
                                  init_abio_dic_dic14_init_file_fmt, &
                                  abio_dic_dic14_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

      call read_field(init_abio_dic_dic14_init_file_fmt, &
                      abio_dic_dic14_restart_filename,   &
                      'ABIO_PH_SURF', ABIO_PH_PREV)


      if (abio_use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(abio_dic_ind)   = abio_surf_avg_dic_const
         surf_avg(abio_dic14_ind) = abio_surf_avg_dic14_const
      else
         call extract_surf_avg(init_abio_dic_dic14_init_file_fmt, &
                               abio_dic_dic14_restart_filename,   &
                               abio_dic_dic14_tracer_cnt,&
                               vflux_flag,ind_name_table,surf_avg)
      endif

      call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

      if (check_time_flag(comp_surf_avg_flag)) &
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:), &
                            abio_dic_dic14_tracer_cnt,        &
                            vflux_flag,surf_avg)


   case ('file', 'ccsm_startup', 'zero', 'ccsm_startup_spunup')

      call document(subname, 'Abiotic DIC and DIC14 being read from separate file')

      call file_read_tracer_block(init_abio_dic_dic14_init_file_fmt, &
                                  init_abio_dic_dic14_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  abio_tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, abio_dic_dic14_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,oldtime,:), errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'abio_dic_dic14_init: error in fill_points for tracers(oldtime)')
                  return
               endif
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'abio_dic_dic14_init: error in fill_points for tracers(newtime)')
                  return
               endif

            end do
         end do
      endif

      ABIO_PH_PREV = c0

      if (abio_use_nml_surf_vals) then
         surf_avg = c0
         surf_avg(abio_dic_ind)   = abio_surf_avg_dic_const
         surf_avg(abio_dic14_ind) = abio_surf_avg_dic14_const
      else
         call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
                            TRACER_MODULE(:,:,1,:,curtime,:), &
                            abio_dic_dic14_tracer_cnt,        &
                            vflux_flag,surf_avg)
      endif

   case default
      call document(subname, 'init_abio_dic_dic14_option', init_abio_dic_dic14_option)
      call exit_POP(sigAbort, 'unknown init_abio_dic_dic14_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, abio_dic_dic14_tracer_cnt
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
   LAND_MASK(:,:,:) = merge(.true., .false., KMT > 0)

   call get_timer(abio_dic_dic14_sflux_timer,&
                  'ABIO_DIC_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  Define decay variable for DIC14, using the earlier defined half life of c14
!-----------------------------------------------------------------------

   c14_lambda_inv_sec = log(c2) / (c14_halflife_years * seconds_in_year) !following OCCMIP2

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call abio_dic_dic14_init_tavg

   call abio_dic_dic14_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_init

!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_init_tavg
! !INTERFACE:

 subroutine abio_dic_dic14_init_tavg

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
      var_cnt             ! how many surface tavg variables are defined?

!-----------------------------------------------------------------------
! Surface tavg variables

   var_cnt = 0

   call define_tavg_field(tavg_ABIO_CO2_IFRAC,'ABIO_CO2_IFRAC',2,         &
                          long_name='Ice Fraction for Abiotic CO2 fluxes',&
                          units='fraction', grid_loc='2110',              &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_CO2_XKW,'ABIO_CO2_XKW',2,             &
                          long_name='XKW for Abiotic CO2 fluxes',         &
                          units='cm/s', grid_loc='2110',                  &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_CO2_ATM_PRESS,'ABIO_CO2_ATM_PRESS',2,         &
                          long_name='Atmospheric Pressure for Abiotic CO2 fluxes',&
                          units='atmospheres', grid_loc='2110',                   &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_pCO2,'ABIO_pCO2',2,                        &
                          long_name='Abiotic CO2 atmospheric partial pressure',&
                          units='ppm', grid_loc='2110',                        &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_D14Catm,'ABIO_D14Catm',2,                 &
                          long_name='Abiotic atmospheric Delta C14 in permil',&
                          units='permil', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1


   call define_tavg_field(tavg_ABIO_CO2_SCHMIDT,'ABIO_CO2_SCHMIDT',2,   &
                          long_name='Abiotic CO2 Schmidt Number',       &
                          units='none', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_CO2_PV,'ABIO_CO2_PV',2,             &
                          long_name='Abiotic CO2 piston velocity',      &
                          units='cm/s', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1


   call define_tavg_field(tavg_ABIO_CO2STAR,'ABIO_CO2STAR',2, &
                          long_name='ABIO_CO2STAR',           &
                          units='nmol/cm^3', grid_loc='2110', &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_DCO2STAR,'ABIO_DCO2STAR',2, &
                          long_name='ABIO_DCO2STAR',            &
                          units='nmol/cm^3', grid_loc='2110',   &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_pCO2SURF,'ABIO_pCO2SURF',2, &
                          long_name='ABIO_pCO2SURF',            &
                          units='ppmv', grid_loc='2110',        &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_DpCO2,'ABIO_DpCO2',2, &
                          long_name='ABIO_DpCO2',         &
                          units='ppmv', grid_loc='2110',  &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_PH,'ABIO_PH_SURF',2, &
                          long_name='Abiotic Surface PH',&
                          units='none', grid_loc='2110', &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ABIO_ALK,'ABIO_ALK_SURF',2,       &
                          long_name='Abiotic Surface Alkalinity',&
                          units='neq/cm3', grid_loc='2110',      &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1


   call define_tavg_field(tavg_FG_ABIO_DIC,'FG_ABIO_DIC',2,           &
                          long_name='Surface gas flux of abiotic DIC',&
                          units='nmol/cm^2/s', grid_loc='2110',       &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_FG_ABIO_DIC14,'FG_ABIO_DIC14',2,         &
                          long_name='Surface gas flux of abiotic DIC14',&
                          units='nmol/cm^2/s', grid_loc='2110',         &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_d_SF_ABIO_DIC_d_ABIO_DIC,'d_SF_ABIO_DIC_d_ABIO_DIC',2, &
                          long_name='derivative of SF_ABIO_DIC wrt ABIO_DIC', &
                          units='cm/s', grid_loc='2110', &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_d_SF_ABIO_DIC14_d_ABIO_DIC,'d_SF_ABIO_DIC14_d_ABIO_DIC',2, &
                          long_name='derivative of SF_ABIO_DIC14 wrt ABIO_DIC', &
                          units='cm/s', grid_loc='2110', &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_d_SF_ABIO_DIC14_d_ABIO_DIC14,'d_SF_ABIO_DIC14_d_ABIO_DIC14',2, &
                          long_name='derivative of SF_ABIO_DIC14 wrt ABIO_DIC14', &
                          units='cm/s', grid_loc='2110', &
                          coordinates='TLONG TLAT time')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------
! Allocate variable to save surface variables to in surf_flux,
! so tavg can be called later for them (in tavg_forcing)
!-----------------------------------------------------------------------

   allocate(ABIO_DIC_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   ABIO_DIC_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
! Additional interior tavg variables
!-----------------------------------------------------------------------

  call define_tavg_field(tavg_ABIO_D14Cocn,'ABIO_D14Cocn',3,              &
                          long_name='Abiotic oceanic Delta C14 in permil',&
                          units='permil', grid_loc='3111',                &
                          coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_init_tavg


!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_init_sflux
! !INTERFACE:

 subroutine abio_dic_dic14_init_sflux

! !USES:
   use named_field_mod, only: named_field_get_index
   use registry, only: registry_match
   use c14_atm_forcing_mod, only: c14_atm_forcing_init

! !DESCRIPTION:
!  Initialize surface flux computations for the abio_dic_dic14 tracer module.
!  Includes reading CO2 and D14C data from file if option file is used
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'abio_dic_dic14_mod:abio_dic_dic14_init_sflux'

   character(len=char_len) :: message

!-------------------------------------------------------------------------
!     READ in CO2 data from files for option file or get
!     CO2 from coupler for option drv_diag
!-------------------------------------------------------------------------

   call document(subname, 'abio_atm_co2_opt', abio_atm_co2_opt)

   select case (abio_atm_co2_opt)

   case ('const')

      call document(subname, 'abio_atm_co2_const', abio_atm_co2_const)

   case ('drv_diag')

!-----------------------------------------------------------------------
!     Verify running coupled if gas fluxes use drv_diag forcing
!-----------------------------------------------------------------------

      if (.not. registry_match('lcoupled')) then
         write(message, *) 'abio_dic_dic14_init: abio_dic_dic14 module requires ', &
                           ' the flux coupler when abio_atm_co2_opt=drv_diag'
         call exit_POP(sigAbort, message)
      endif

!-----------------------------------------------------------------------
!    Get co2 index from coupler for reading of CO2 data later
!    (each timesstep) and abort if it isn't there
!-----------------------------------------------------------------------

      call named_field_get_index('ATM_CO2_DIAG', atm_co2_nf_ind, &
                                 exit_on_err=.false.)

      if (atm_co2_nf_ind == 0) then
         write(message, *) 'abio_dic_dic14_init: abio_dic_dic14 module requires ', &
                           'atmopsheric CO2 from the flux coupler and it is not ', &
                           'present'
         call exit_POP(sigAbort, message)
      endif


   case('file')

!-----------------------------------------------------------------------
!     READ in CO2 data from file
!-----------------------------------------------------------------------

      call read_atm_CO2_data

  case default

      call exit_POP(sigAbort, 'unknown abio_atm_co2_opt')

  end select

!-------------------------------------------------------------------------

  call c14_atm_forcing_init('abio_dic_dic14', abio_atm_d14c_opt, &
      abio_atm_d14c_const, abio_atm_d14c_lat_band_vals, &
      abio_atm_d14c_filename, abio_atm_model_year, abio_atm_data_year)

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_init_sflux


!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_set_sflux
! !INTERFACE:

 subroutine abio_dic_dic14_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute abiotic CO2 and C14 surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants,           only: ocn_ref_salinity, rho_sw
   use timers,              only: timer_start, timer_stop
   use named_field_mod,     only: named_field_get
   use grid,                only: REGION_MASK
   use c14_atm_forcing_mod, only: c14_atm_forcing_update_data, c14_atm_forcing_get_data
   use tavg,                only: accumulate_tavg_now

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      PRESS,     & ! sea level atmospheric pressure (dyne/cm**2)
      SST,       & ! sea surface temperature (C)
      SSS          ! sea surface salinity (psu)

   real (r8), dimension(nx_block,ny_block,abio_dic_dic14_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,abio_dic_dic14_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock,            & ! block index
      j                    ! looping index

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,        & ! used ice fraction (non-dimensional)
      XKW_USED,          & ! part of piston velocity (cm/s)
      AP_USED              ! used atm pressure (converted from dyne/cm**2 to atm)

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS_DIC,      & ! filtered DIC surface tracer values
      SURF_VALS_DIC14,    & ! filtered DIC14 surface tracer values
      pCO2,               & ! atmospheric CO2 mole fraction (pmol/mol)
      D14C,               & ! atmospheric Delta C14 in permil
      CO2_SCHMIDT_USED,   & ! CO2 Schmidt number
      XKW_ICE,            & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,                 & ! piston velocity (cm/s)
      SiO2,               & ! Silicate (constant), in mol/m3 =mmol/cm3
      PO4,                & ! Phosphate (constant), in mol/m3 =mmol/cm3
      R14C_ocn,           & ! Rocn=DIC14/DIC
      R14C_atm,           & ! Ratm = 1+ D14C/1000
      GAS_FLUX_ABIO_DIC,  & ! Surface gas flux of DIC
      GAS_FLUX_ABIO_DIC14   ! Surface flux of DIC14


      real (r8), dimension(nx_block) :: &
      PHLO,             & ! lower bound for ph in solver
      PHHI,             & ! upper bound for ph in solver
      ABIO_DIC_ROW,     & ! row of DIC values for solver
      CO2STAR_ROW,      & ! CO2STAR from solver
      CO2STAR_ROW_TMP,  & ! CO2STAR from solver
      DCO2STAR_ROW,     & ! DCO2STAR from solver
      pCO2SURF_ROW,     & ! pCO2SURF from solver
      DpCO2_ROW,        & ! DpCO2 from solver
      PH_NEW,           & ! new PH
      ALK_ROW,          & ! Alkalinty in nano eq/cm3
      CO3_ROW             ! CO3 at the surface

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------
   real (r8), parameter :: &
      ALK_bar_global   = 2310._r8     ! Alkbar, in microeq/kg, based on OCMIP2

!-----------------------------------------------------------------------

   call timer_start(abio_dic_dic14_sflux_timer)

   if (check_time_flag(comp_surf_avg_flag))  &
      call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR,abio_dic_dic14_tracer_cnt, &
                            vflux_flag,surf_avg)

!-----------------------------------------------------------------------

   call c14_atm_forcing_update_data

!-----------------------------------------------------------------------

! Initilize fields to zero
!$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock)   = c0
      XKW_USED(:,:,iblock)     = c0
      AP_USED(:,:,iblock)      = c0
      STF_MODULE(:,:,:,iblock) = c0
   end do
!$OMP END PARALLEL DO
!-----------------------------------------------------------------------------------------------
!  Compute CO2 flux, computing disequilibrium one row at a time - modified from ecosy_mod.F90
!-----------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate and initialize IFRAC, XKW and AP (pressure)
!-----------------------------------------------------------------------
 !$OMP PARALLEL DO PRIVATE(iblock,j,XKW_ICE,CO2_SCHMIDT_USED,PV,SiO2, PO4,&
 !$OMP                     pCO2, D14C,SURF_VALS_DIC,SURF_VALS_DIC14,&
 !$OMP                     R14C_ocn,R14C_atm,PHLO,PHHI,ABIO_DIC_ROW,&
 !$OMP                     ALK_ROW,PH_NEW,CO2STAR_ROW,CO2STAR_ROW_TMP,DCO2STAR_ROW,&
 !$OMP                     pCO2SURF_ROW,DpCO2_ROW,GAS_FLUX_ABIO_DIC,&
 !$OMP                     GAS_FLUX_ABIO_DIC14,CO3_ROW)

   do iblock = 1, nblocks_clinic
      where (LAND_MASK(:,:,iblock))
         IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)
         XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)
         AP_USED(:,:,iblock) = PRESS(:,:,iblock)
      endwhere

      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
         IFRAC_USED(:,:,iblock) = c0
      where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
         IFRAC_USED(:,:,iblock) = c1

!-----------------------------------------------------------------------
!  assume PRESS is in cgs units (dyne/cm**2) since that is what is
!    required for pressure forcing in barotropic
!  want units to be atmospheres
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) * (c1 / 1013.25e+3_r8)

!-----------------------------------------------------------------------
!  Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too
!-----------------------------------------------------------------------

      XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)

!-----------------------------------------------------------------------
! Compute CO2 Schmidt number
!-----------------------------------------------------------------------

      CO2_SCHMIDT_USED = SCHMIDT_CO2(SST(:,:,iblock), LAND_MASK(:,:,iblock))

!-----------------------------------------------------------------------
!  compute piston velocity PV
!-----------------------------------------------------------------------

      where (LAND_MASK(:,:,iblock))
         PV = XKW_ICE * SQRT(660.0_r8 / CO2_SCHMIDT_USED)
      elsewhere
         PV = c0
      end where

!-----------------------------------------------------------------------
! Assign constant SiO2, PO4 data to variables
!-----------------------------------------------------------------------

   SiO2  = 7.5_r8 * rho_sw      ! Silicate (constant), in mmol/m3 = nmol/cm3 (original value is 7.5 micromol/kg)
   PO4   = 0.5_r8 * rho_sw      ! Phosphate (constant), in mmol/m3 = nmol/cm3 (orignal value is 0.5 micromol/kg)

!-----------------------------------------------------------------------
! Assign CO2 data to variables
!-----------------------------------------------------------------------
      select case (abio_atm_co2_opt)

      case ('const')
            pCO2 = abio_atm_co2_const

      case ('drv_diag')

         call named_field_get(atm_co2_nf_ind,iblock, pCO2)

      case('file')
         call comp_varying_CO2(data_ind_co2(iblock), pCO2)

      end select

!-----------------------------------------------------------------------

      call c14_atm_forcing_get_data(iblock, D14C)

!-----------------------------------------------------------------------
! Make filtered Surface DIC and DIC14 variables
!-----------------------------------------------------------------------

      SURF_VALS_DIC = p5*(SURF_VALS_OLD(:,:,abio_dic_ind,iblock) + &
                          SURF_VALS_CUR(:,:,abio_dic_ind,iblock))

      SURF_VALS_DIC14 = p5*(SURF_VALS_OLD(:,:,abio_dic14_ind,iblock) + &
                            SURF_VALS_CUR(:,:,abio_dic14_ind,iblock))

!-----------------------------------------------------------------------
! Calculate R14C_ocn and R14C_atm
!-----------------------------------------------------------------------

      where(SURF_VALS_DIC/=c0)
         R14C_ocn = SURF_VALS_DIC14/SURF_VALS_DIC
      elsewhere
         R14C_ocn = c0
      endwhere

      R14C_atm = c1 + D14C / c1000

!-----------------------------------------------------------------------
!  Use co2calc_row to calculate Csurf following OCMIP-2, using ecosystem
!  code and co2calc.F90, but adding alkalinity calculation
!-----------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------
! Set PH bracket
!--------------------------------------------------------------------------------------------------------

      do j = 1,ny_block
         where (ABIO_PH_PREV(:,j,iblock) /= c0)
            PHLO = ABIO_PH_PREV(:,j,iblock) - abio_del_ph
            PHHI = ABIO_PH_PREV(:,j,iblock) + abio_del_ph
         elsewhere
            PHLO = abio_phlo_surf_init
            PHHI = abio_phhi_surf_init
         end where

         ABIO_DIC_ROW = p5*(SURF_VALS_OLD(:,j,abio_dic_ind,iblock) + &
                        SURF_VALS_CUR(:,j,abio_dic_ind,iblock))

!--------------------------------------------------------------------------------------------------------
!  Calculate Alkalinity, with different values in the marginal seas,
!  following the same procedure as used for the creation of the ecosystem initial alkalinity conditions
!  This means the Baltic and Black Sea need to be treated differently than the rest of the global ocean
!  In the Baltic Sea the Surface Alkalinity depends on the salinity, in the Black Sea the Surface Alkalinity
!  does not depend on salinity
!--------------------------------------------------------------------------------------------------------

         where ( (REGION_MASK(:,j,iblock) /= -12) .AND. (REGION_MASK(:,j,iblock) /= -13) )
            ALK_ROW=ALK_bar_global * rho_sw * SSS(:,j,iblock) / ocn_ref_salinity        ! units= nano eq/cm3
         elsewhere ( (REGION_MASK(:,j,iblock) == -12) .AND. (SSS(:,j,iblock) <=7.3) )   ! Baltic Sea, fresh part
            ALK_ROW=119._r8 + 196._r8 * SSS(:,j,iblock)                                 ! units= nano eq/cm3
         elsewhere ( (REGION_MASK(:,j,iblock) == -12) .AND. (SSS(:,j,iblock) >7.3) )    ! Baltic Sea, saltier part
            ALK_ROW=1237._r8 + 43._r8 * SSS(:,j,iblock)                                 ! units= nano eq/cm3
         elsewhere (REGION_MASK(:,j,iblock) == -13)                                     ! Black Sea
            ALK_ROW= 3300 * rho_sw                                                      ! units= nano eq/cm3
         end where

!-----------------------------------------------------------------------
!  Call co2calc_row which calculated PH, CO2* etc
!-----------------------------------------------------------------------

         call co2calc_row(iblock, j, LAND_MASK(:,j,iblock), .true., &
                          .true., SST(:,j,iblock), SSS(:,j,iblock), &
                          ABIO_DIC_ROW, ALK_ROW, PO4(:,j), SiO2(:,j), &
                          PHLO, PHHI, PH_NEW, pCO2(:,j), &
                          AP_USED(:,j,iblock), CO2STAR_ROW, &
                          DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW,CO3_ROW)

         ABIO_PH_PREV(:,j,iblock) = PH_NEW

!-----------------------------------------------------------------------
! Calculate the surface gas flux for DIC and DIC14
!-----------------------------------------------------------------------
         GAS_FLUX_ABIO_DIC(:,j)   = PV(:,j) * DCO2STAR_ROW

         GAS_FLUX_ABIO_DIC14(:,j) = PV(:,j) * ((DCO2STAR_ROW + CO2STAR_ROW) * &
                                     R14C_atm(:,j) - CO2STAR_ROW * R14C_ocn(:,j))

 !-----------------------------------------------------------------------
 ! Assign calculated values to TAVG fields
 !-----------------------------------------------------------------------
         ABIO_DIC_SFLUX_TAVG(:,j,10,iblock)  = CO2STAR_ROW
         ABIO_DIC_SFLUX_TAVG(:,j,11,iblock)  = DCO2STAR_ROW
         ABIO_DIC_SFLUX_TAVG(:,j,12,iblock)  = pCO2SURF_ROW
         ABIO_DIC_SFLUX_TAVG(:,j,13,iblock)  = DpCO2_ROW
         ABIO_DIC_SFLUX_TAVG(:,j,14,iblock)  = ALK_ROW

!-----------------------------------------------------------------------
!  compute derivatives of surface flux terms, if requested
!-----------------------------------------------------------------------

         if (accumulate_tavg_now(tavg_d_SF_ABIO_DIC14_d_ABIO_DIC14)) then

            ! d_SF_ABIO_DIC14_d_ABIO_DIC14
            where (SURF_VALS_DIC(:,j) /= c0)
               ABIO_DIC_SFLUX_TAVG(:,j,18,iblock) = (-1) * PV(:,j) * CO2STAR_ROW / SURF_VALS_DIC(:,j)
            elsewhere
               ABIO_DIC_SFLUX_TAVG(:,j,18,iblock) = c0
            endwhere

         endif

         if (accumulate_tavg_now(tavg_d_SF_ABIO_DIC_d_ABIO_DIC) .or. &
             accumulate_tavg_now(tavg_d_SF_ABIO_DIC14_d_ABIO_DIC)) then

            call co2calc_row(iblock, j, LAND_MASK(:,j,iblock), .true., &
                             .false., SST(:,j,iblock), SSS(:,j,iblock), &
                             ABIO_DIC_ROW + c1, ALK_ROW, PO4(:,j), SiO2(:,j), &
                             PHLO, PHHI, PH_NEW, pCO2(:,j), &
                             AP_USED(:,j,iblock), CO2STAR_ROW_TMP, &
                             DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, CO3_ROW)

            call co2calc_row(iblock, j, LAND_MASK(:,j,iblock), .true., &
                             .false., SST(:,j,iblock), SSS(:,j,iblock), &
                             ABIO_DIC_ROW - c1, ALK_ROW, PO4(:,j), SiO2(:,j), &
                             PHLO, PHHI, PH_NEW, pCO2(:,j), &
                             AP_USED(:,j,iblock), CO2STAR_ROW, &
                             DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, CO3_ROW)

            ! d_SF_ABIO_DIC_d_ABIO_DIC
            ABIO_DIC_SFLUX_TAVG(:,j,16,iblock) = (-1) * PV(:,j) * &
               (CO2STAR_ROW_TMP - CO2STAR_ROW) * p5

            ! d_SF_ABIO_DIC14_d_ABIO_DIC
            where (SURF_VALS_DIC(:,j) /= c0)
               ABIO_DIC_SFLUX_TAVG(:,j,17,iblock) = (-1) * PV(:,j) * SURF_VALS_DIC14(:,j) * &
                  (CO2STAR_ROW_TMP / (SURF_VALS_DIC(:,j) + c1) - CO2STAR_ROW / (SURF_VALS_DIC(:,j) - c1)) * p5
            elsewhere
               ABIO_DIC_SFLUX_TAVG(:,j,17,iblock) = c0
            endwhere

         endif

      end do !j = 1,ny_block
!-----------------------------------------------------------------------
! Update tracer fields with the surface fluxes just calculated
!-----------------------------------------------------------------------

      STF_MODULE(:,:,abio_dic14_ind,iblock) = STF_MODULE(:,:,abio_dic14_ind,iblock) + &
                                              GAS_FLUX_ABIO_DIC14

      STF_MODULE(:,:,abio_dic_ind,iblock)   = STF_MODULE(:,:,abio_dic_ind,iblock) +   &
                                              GAS_FLUX_ABIO_DIC

 !-----------------------------------------------------------------------
 ! Write all remaining fields to TAVG fields, and set to zero over land
 !-----------------------------------------------------------------------

      where (LAND_MASK(:,:,iblock))
         ABIO_DIC_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         ABIO_DIC_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         ABIO_DIC_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
         ABIO_DIC_SFLUX_TAVG(:,:,4,iblock) = pCO2
         ABIO_DIC_SFLUX_TAVG(:,:,5,iblock) = PV
         ABIO_DIC_SFLUX_TAVG(:,:,6,iblock) = CO2_SCHMIDT_USED
         ABIO_DIC_SFLUX_TAVG(:,:,7,iblock) = ABIO_PH_PREV(:,:,iblock)
         ABIO_DIC_SFLUX_TAVG(:,:,8,iblock) = GAS_FLUX_ABIO_DIC
         ABIO_DIC_SFLUX_TAVG(:,:,9,iblock) = GAS_FLUX_ABIO_DIC14
         ABIO_DIC_SFLUX_TAVG(:,:,15,iblock)= D14C
      elsewhere
         ABIO_DIC_SFLUX_TAVG(:,:,1,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,2,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,3,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,4,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,5,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,6,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,7,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,8,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,9,iblock) = c0
         ABIO_DIC_SFLUX_TAVG(:,:,15,iblock)= c0
      endwhere

   end do !iblock = 1, nblocks_clinic
!$OMP END PARALLEL DO


 call timer_stop(abio_dic_dic14_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_set_sflux



!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_set_interior
! !INTERFACE:

 subroutine abio_dic_dic14_set_interior(k, TRACER_MODULE_OLD, &
    TRACER_MODULE_CUR, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Compute interior source-sink terms for tracers - here the decay of C14
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   real (r8), dimension(nx_block,ny_block,km,abio_dic_dic14_tracer_cnt), intent(in) :: &
      TRACER_MODULE_OLD, &! old tracer values
      TRACER_MODULE_CUR   ! current tracer values

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,abio_dic_dic14_tracer_cnt), intent(out) :: &
      DTRACER_MODULE      ! computed source/sink terms

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                 ! local_block id

   real (r8), dimension(nx_block,ny_block) :: &
      DIC14_loc, &        ! local copy of model DIC14
      DIC_loc, &          ! local copy of model DIC
      ABIO_D14Cocn        ! D14Cocn value


   bid = this_block%local_id

   DTRACER_MODULE(:,:,abio_dic_ind)   = c0

!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!  apply mask to local copies
!-----------------------------------------------------------------------

   DIC14_loc = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,abio_dic14_ind) + &
                           TRACER_MODULE_CUR(:,:,k,abio_dic14_ind)))

   DIC_loc = max(c0, p5*(TRACER_MODULE_OLD(:,:,k,abio_dic_ind) + &
                           TRACER_MODULE_CUR(:,:,k,abio_dic_ind)))

   where (.not. LAND_MASK(:,:,bid))
      DIC14_loc = c0
      DIC_loc   = c0
   end where

   DTRACER_MODULE(:,:,abio_dic14_ind) = -c14_lambda_inv_sec * DIC14_loc

   where (DIC_loc == 0)
      ABIO_D14Cocn = c0
   elsewhere
      ABIO_D14Cocn =  (DIC14_loc / DIC_loc -1) *1000
   end where

   call accumulate_tavg_field(ABIO_D14Cocn,tavg_ABIO_D14Cocn,bid,k)

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_set_interior

!***********************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_tavg_forcing
! !INTERFACE:

 subroutine abio_dic_dic14_tavg_forcing

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
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,1,iblock),tavg_ABIO_CO2_IFRAC,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,2,iblock),tavg_ABIO_CO2_XKW,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,3,iblock),tavg_ABIO_CO2_ATM_PRESS,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,4,iblock),tavg_ABIO_pCO2,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,5,iblock),tavg_ABIO_CO2_PV,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,6,iblock),tavg_ABIO_CO2_SCHMIDT,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,7,iblock),tavg_ABIO_PH,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,8,iblock),tavg_FG_ABIO_DIC,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,9,iblock),tavg_FG_ABIO_DIC14,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,10,iblock),tavg_ABIO_CO2STAR,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,11,iblock),tavg_ABIO_DCO2STAR,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,12,iblock),tavg_ABIO_pCO2SURF,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,13,iblock),tavg_ABIO_DpCO2,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,14,iblock),tavg_ABIO_ALK,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,15,iblock),tavg_ABIO_D14Catm,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,16,iblock),tavg_d_SF_ABIO_DIC_d_ABIO_DIC,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,17,iblock),tavg_d_SF_ABIO_DIC14_d_ABIO_DIC,iblock,1)
      call accumulate_tavg_field(ABIO_DIC_SFLUX_TAVG(:,:,18,iblock),tavg_d_SF_ABIO_DIC14_d_ABIO_DIC14,iblock,1)
   end do



   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC


end subroutine abio_dic_dic14_tavg_forcing


!***********************************************************************
! !IROUTINE: read_atm_CO2_data
! !INTERFACE:

 subroutine read_atm_CO2_data

! !DESCRIPTION:
!  Read atmospheric CO2 [ppm] data from file
!
!  Have the master_task do the following :
!     1) get length of data
!     2) allocate memory for data
!     3) read in data, checking for consistent lengths
!  Then, outside master_task conditional
!     1) broadcast length of data
!     2) have non-mastertasks allocate memory for data
!     3) broadcast data
!
! !REVISION HISTORY:
!  same as module
!
! !USES:
!
!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
   character(*), parameter :: subname = 'abio_dic_dci14_mod:read_atm_CO2_data'

   integer (int_kind) ::    &
      stat,                 &  ! i/o status code
      irec,                 &  ! counter for looping
      skiplines                ! number of comment lines at beginning of ascii file

   character (char_len) ::  &
      sglchr                   ! variable to read characters from file into


!-----------------------------------------------------------------------
!     READ in CO2 data from file
!-----------------------------------------------------------------------
   call document(subname, 'abio_atm_co2_filename', abio_atm_co2_filename)

   if (my_task == master_task) then
      open (nml_in, file=abio_atm_co2_filename, status='old',iostat=stat)
      if (stat /= 0) then
         write(stdout,fmt=*) 'open failed'
         go to 99
      endif
      read(nml_in,FMT=*,iostat=stat) skiplines,abio_atm_co2_data_nbval
      if (stat /= 0) then
         write(stdout,fmt=*) '1st line read failed'
         go to 99
      endif
      allocate(abio_atm_co2_data_yr(abio_atm_co2_data_nbval))
      allocate(abio_atm_co2_data_ppm(abio_atm_co2_data_nbval))
      do irec=1,skiplines
         read(nml_in,FMT=*,iostat=stat) sglchr
         if (stat /= 0) then
            write(stdout,fmt=*) 'skipline read failed'
            go to 99
         endif
      enddo
      do irec=1,abio_atm_co2_data_nbval
         read(nml_in,FMT=*,iostat=stat) abio_atm_co2_data_yr(irec), abio_atm_co2_data_ppm(irec)
         if (stat /= 0) then
            write(stdout,fmt=*) 'data read failed'
            go to 99
         endif
      enddo
      close(nml_in)
   endif

99 call broadcast_scalar(stat, master_task)
   if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                          &/ subname)

!---------------------------------------------------------------------
! Need to allocate and broadcast the variables to other tasks beside master_task
!---------------------------------------------------------------------

   call broadcast_scalar(abio_atm_co2_data_nbval,master_task)

   if (my_task /= master_task) then
      allocate(abio_atm_co2_data_yr(abio_atm_co2_data_nbval))
      allocate(abio_atm_co2_data_ppm(abio_atm_co2_data_nbval))
   endif

   call broadcast_array(abio_atm_co2_data_ppm, master_task)
   call broadcast_array(abio_atm_co2_data_yr, master_task)

   allocate( data_ind_co2(nblocks_clinic) )
   data_ind_co2 = -1

!-----------------------------------------------------------------------
!EOC

 end subroutine read_atm_CO2_data

!***********************************************************************
! !IROUTINE: comp_varying_CO2
! !INTERFACE:

 subroutine comp_varying_CO2(data_ind_co2, pCO2)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of CO2 when temporarily
!  varying data is read from files
!  1. Linearly interpolate data values to current model time step
!  2. Spatial patern of CO2 is the same everywhere (90 S - 90 N)
!
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind) :: &
      data_ind_co2     ! data_ind_co2 is the index for the data for current timestep,
                       ! note that data_ind_co2 is always strictly less than the length
                       ! of the data and is initialized to -1 during init

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      pCO2             ! atmospheric CO2 mole fraction (pmol/mol)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8) :: &
      model_date,     & ! date of current model timestep
      mapped_date,    & ! model_date mapped to data timeline
      weight            ! weighting for temporal interpolation

!-----------------------------------------------------------------------
!  Generate mapped_date and check to see if it is too large.
!-----------------------------------------------------------------------

   model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
   mapped_date = model_date - abio_atm_model_year + abio_atm_data_year

   if (mapped_date >= abio_atm_co2_data_yr(abio_atm_co2_data_nbval)) then
      call exit_POP(sigAbort, 'model date maps to date after end of CO2 data in file')
   endif

!--------------------------------------------------------------------------------------------------------------
!  Set atmospheric CO2 to first value in record for years before record begins
!--------------------------------------------------------------------------------------------------------------

   if (mapped_date < abio_atm_co2_data_yr(1)) then
      pCO2 = abio_atm_co2_data_ppm(1)
      data_ind_co2 = 1
      if(my_task == master_task) then
         write(stdout,*)'Mapped date less than start of CO2 data --> using first value in CO2 data file'
      endif
      return
   endif

!-----------------------------------------------------------------------
!  On first time step, perform linear search to find data_ind_co2
!-----------------------------------------------------------------------

   if (data_ind_co2 == -1) then
      do data_ind_co2 = abio_atm_co2_data_nbval-1,1,-1
         if (mapped_date >= abio_atm_co2_data_yr(data_ind_co2)) exit
      end do
   endif

!-----------------------------------------------------------------------
!  See if data_ind_co2 needs to be updated,
!  but do not set it to abio_atm_co2_data_nbval.
!-----------------------------------------------------------------------

  if (data_ind_co2 < abio_atm_co2_data_nbval-1) then
      if (mapped_date >= abio_atm_co2_data_yr(data_ind_co2+1)) data_ind_co2 = data_ind_co2 + 1
  endif


!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!-----------------------------------------------------------------------

   weight = (mapped_date - abio_atm_co2_data_yr(data_ind_co2)) &
            / (abio_atm_co2_data_yr(data_ind_co2+1) - abio_atm_co2_data_yr(data_ind_co2))

   pCO2 = weight * abio_atm_co2_data_ppm(data_ind_co2+1) + (c1-weight) * abio_atm_co2_data_ppm(data_ind_co2)

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_varying_CO2

!*****************************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_write_restart
! !INTERFACE:

 subroutine abio_dic_dic14_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module
  use constants, only: char_blank, field_loc_center, field_type_scalar

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (char_len) :: &
      short_name   ! tracer name temporaries

   type (io_dim) :: &
      i_dim, j_dim ! dimension descriptors

   integer (int_kind) :: n

   type (io_field_desc), save :: ABIO_PH_SURF

!-----------------------------------------------------------------------

   if (trim(action) == 'add_attrib_file') then
      short_name = char_blank
      do n=1,abio_dic_dic14_tracer_cnt
         if (vflux_flag(n)) then
            short_name = 'surf_avg_' /&
                      &/ trim(ind_name_table(n)%name)
            call add_attrib_file(restart_file,trim(short_name),surf_avg(n))
         endif
      end do
   endif

   if (trim(action) == 'define') then
      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)

      ABIO_PH_SURF = construct_io_field('ABIO_PH_SURF', i_dim, j_dim,&
                   long_name='Abiotic surface pH at current time',   &
                   units    ='pH', grid_loc ='2110',                 &
                   field_loc = field_loc_center,                     &
                   field_type = field_type_scalar,                   &
                   d2d_array = ABIO_PH_PREV)
      call data_set (restart_file, 'define', ABIO_PH_SURF)
   endif

   if (trim(action) == 'write') then
      call data_set (restart_file, 'write', ABIO_PH_SURF)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine abio_dic_dic14_write_restart


!*****************************************************************************
!BOP
! !IROUTINE: abio_dic_dic14_tracer_ref_val
! !INTERFACE:

 function abio_dic_dic14_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: abio_dic_dic14_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------


   if (vflux_flag(ind)) then
      abio_dic_dic14_tracer_ref_val = surf_avg(ind)
   else
      abio_dic_dic14_tracer_ref_val = c0
   endif

!-----------------------------------------------------------------------
!EOC

 end function abio_dic_dic14_tracer_ref_val


end module abio_dic_dic14_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
