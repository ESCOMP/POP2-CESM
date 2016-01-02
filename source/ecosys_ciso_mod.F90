! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_ciso_mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !
  !  Carbon 13 module and biotic 14C module
  !  13C code is based on code form G. Xavier, ETH, 2010, which
  !  was written for pop1 (CCSM3)
  !  This code needs the ecosystem model to run, as it uses several
  !  variables computed there. Data is shared using marbl_share_mod.F90
  !  This module adds 7 carbon pools for 13C and another 7 for 14C
  !
  !  Developer: Alexandra Jahn, NCAR, Started Nov 2012, last edited July 2014
  !  A first version of the 13C code for POP1, which was used to develop 
  !  the code for the current model version, was written by Xavier Giraud 
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !  variables/subroutines/function used from other modules
  !  The following are used extensively in this module, so they are used at
  !  the module level. The use statements for variables that are only needed
  !  locally are located at the module subprogram level.
  !-----------------------------------------------------------------------

  ! !USES:

  use POP_KindsMod, only : POP_i4
  use POP_ErrorMod, only : POP_ErrorSet
  use POP_ErrorMod, only : POP_Success

  use kinds_mod, only : r8
  use kinds_mod, only : int_kind
  use kinds_mod, only : log_kind
  use kinds_mod, only : char_len

  use constants, only : c0
  use constants, only : c1
  use constants, only : c2
  use constants, only : c1000
  use constants, only : p5
  use constants, only : char_blank
  use constants, only : blank_fmt
  use constants, only : delim_fmt
  use constants, only : ndelim_fmt
  use constants, only : mpercm

  use communicate, only : master_task
  use communicate, only : my_task

  use broadcast, only : broadcast_array
  use broadcast, only : broadcast_scalar

  use global_reductions, only : global_sum

  use blocks, only : block
  use blocks, only : nx_block
  use blocks, only : ny_block
  use blocks, only : get_block

  use domain_size, only : km
  use domain_size, only : max_blocks_clinic

  use domain, only : distrb_clinic
  use domain, only : nblocks_clinic
  use domain, only : blocks_clinic

  use exit_mod, only : exit_POP
  use exit_mod, only : sigAbort

  use marbl_interface_types, only : marbl_tracer_metadata_type
  use marbl_interface_types, only : marbl_tracer_read_type
  use prognostic, only : curtime
  use prognostic, only : oldtime

  use grid, only : KMT
  use grid, only : TAREA
  use grid, only : DZT
  use grid, only : n_topo_smooth
  use grid, only : REGION_MASK
  use grid, only : TLATD
  use grid, only : dz
  use grid, only : zw
  use grid, only : fill_points

  use io, only : datafile

  use io_types, only : stdout
  use io_types, only : nml_in
  use io_types, only : nml_filename
  use io_types, only : add_attrib_file

  use io_tools, only : document

  use tavg, only : accumulate_tavg_now
  use tavg, only : accumulate_tavg_field
  use tavg, only : define_tavg_field

  use timers, only : timer_start
  use timers, only : timer_stop
  use timers, only : get_timer

  use passive_tracer_tools, only : ind_name_pair
  use passive_tracer_tools, only : comp_surf_avg
  use passive_tracer_tools, only : extract_surf_avg
  use passive_tracer_tools, only : file_read_tracer_block
  use passive_tracer_tools, only : rest_read_tracer_block

  use time_management, only : days_in_year
  use time_management, only : frac_day
  use time_management, only : iday_of_year
  use time_management, only : iyear
  use time_management, only : seconds_in_day
  use time_management, only : seconds_in_year
  use time_management, only : freq_opt_never
  use time_management, only : freq_opt_nmonth
  use time_management, only : freq_opt_nyear
  use time_management, only : check_time_flag
  use time_management, only : init_time_flag
  use time_management, only : eval_time_flag

  use marbl_share_mod, only : ecosys_ciso_tracer_cnt
  use marbl_share_mod, only : autotrophs, autotroph_cnt
  use marbl_share_mod, only : ciso_lsource_sink
  use marbl_share_mod, only : ciso_fract_factors
  use marbl_share_mod, only : ciso_init_ecosys_option
  use marbl_share_mod, only : ciso_init_ecosys_init_file
  use marbl_share_mod, only : ciso_init_ecosys_init_file_fmt
  use marbl_share_mod, only : ciso_atm_model_year                                                            
  use marbl_share_mod, only : ciso_atm_data_year                                                             
  use marbl_share_mod, only : ciso_atm_d13c_data_nbval                                                       
  use marbl_share_mod, only : ciso_atm_d14c_data_nbval                                                       
  use marbl_share_mod, only : ciso_atm_d13c_data                                                             
  use marbl_share_mod, only : ciso_atm_d13c_data_yr                                                          
  use marbl_share_mod, only : ciso_atm_d14c_data                                                             
  use marbl_share_mod, only : ciso_atm_d14c_data_yr                                                          
  use marbl_share_mod, only : ciso_atm_d13c_const                                                            
  use marbl_share_mod, only : ciso_atm_d14c_const                                                            
  use marbl_share_mod, only : ciso_atm_d13c_opt                                                              
  use marbl_share_mod, only : ciso_atm_d14c_opt                                                              
  use marbl_share_mod, only : ciso_atm_d13c_filename                                                         
  use marbl_share_mod, only : ciso_atm_d14c_filename                                                         
  use marbl_share_mod, only : ciso_comp_surf_avg_freq_opt !?

  use marbl_parms, only : di13c_ind
  use marbl_parms, only : do13c_ind
  use marbl_parms, only : zoo13C_ind
  use marbl_parms, only : di14c_ind
  use marbl_parms, only : do14c_ind
  use marbl_parms, only : zoo14C_ind

  !-----------------------------------------------------------------------
  !  include ecosystem parameters
  !  all variables from this modules have a parm_ prefix
  !-----------------------------------------------------------------------

  implicit none
  private

  !-----------------------------------------------------------------------
  !  public/private declarations
  !-----------------------------------------------------------------------

  public :: &
       ecosys_ciso_init,                  &
       ecosys_ciso_tracer_ref_val,        &
       ecosys_ciso_set_sflux,             &
       ecosys_ciso_tavg_forcing,          &
       ecosys_ciso_write_restart

  !-----------------------------------------------------------------------
  !  module variables required by forcing_passive_tracer
  !-----------------------------------------------------------------------

  logical (log_kind), dimension(:,:,:), allocatable :: LAND_MASK

  !-----------------------------------------------------------------------
  !  derived type & parameter for tracer index lookup
  !-----------------------------------------------------------------------

  type(ind_name_pair), dimension(ecosys_ciso_tracer_cnt) :: ciso_ind_name_table

  !-----------------------------------------------------------------------
  !  tavg ids and buffer indices (into ECO_CISO_SFLUX_TAVG) for 2d fields
  !  related to surface fluxes. Suplicates, which are used for placing fields
  !  into multiple tavg streams, do not need separate buffer indices.
  !  fields that are recoverable from the STF field do not need separate
  !  buffer indices
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       tavg_CISO_DI13C_GAS_FLUX,       buf_ind_CISO_FG_13CO2,     & ! di13c flux
       tavg_CISO_DI13C_AS_GAS_FLUX,    buf_ind_CISO_FG_as_13CO2,  & ! air-sea di13c flux
       tavg_CISO_DI13C_SA_GAS_FLUX,    buf_ind_CISO_FG_sa_13CO2,  & ! sea-air di13c flux
       tavg_CISO_d13C_GAS_FLUX,        buf_ind_CISO_FG_d13C,      & ! surface ocean delta 13C
       tavg_CISO_R13C_DIC_surf,        buf_ind_CISO_R13C_DIC_surf,& ! 13C/12C ratio in total DIC
       tavg_CISO_R13C_atm,             buf_ind_CISO_R13C_atm,     & ! atmospheric ratio of 13C/12C
       tavg_CISO_D13C_atm,             buf_ind_CISO_D13C_atm,     & ! atmospheric delta13C in permil
       tavg_CISO_DI14C_GAS_FLUX,       buf_ind_CISO_FG_14CO2,     & ! di14c flux
       tavg_CISO_DI14C_AS_GAS_FLUX,    buf_ind_CISO_FG_as_14CO2,  & ! air-sea di14c flux
       tavg_CISO_DI14C_SA_GAS_FLUX,    buf_ind_CISO_FG_sa_14CO2,  & ! sea-air di14c flux
       tavg_CISO_d14C_GAS_FLUX,        buf_ind_CISO_FG_d14C,      & ! surface ocean delta 14C
       tavg_CISO_R14C_DIC_surf,        buf_ind_CISO_R14C_DIC_surf,& ! 14C/12C ratio in total DIC
       tavg_CISO_R14C_atm,             buf_ind_CISO_R14C_atm,     & ! atmospheric ratio of 14C/12C
       tavg_CISO_D14C_atm,             buf_ind_CISO_D14C_atm,     & ! atmospheric delta14C in permil
       tavg_CISO_DI13C_RIV_FLUX,       buf_ind_DI13C_RIV_FLUX,    & ! river input of DI13C
       tavg_CISO_DO13C_RIV_FLUX,       buf_ind_DO13C_RIV_FLUX,    & ! river input of DO13C
       tavg_CISO_DI14C_RIV_FLUX,       buf_ind_DI14C_RIV_FLUX,    & ! river input of DI14C
       tavg_CISO_DO14C_RIV_FLUX,       buf_ind_DO14C_RIV_FLUX,    & ! river input of DO14C
       tavg_CISO_eps_aq_g_surf,        buf_ind_CISO_eps_aq_g_surf,& ! tavg id for eps_aq_g_surf
       tavg_CISO_eps_dic_g_surf,       buf_ind_CISO_eps_dic_g_surf  ! tavg id for eps_dic_g_surf

  ! for debugging
  integer (int_kind) ::       &
       tavg_CISO_GLOBAL_D14C,   & ! tavg id for the global averaged atmos. D14C
       buf_ind_CISO_GLOBAL_D14C   ! buffer index for the global averaged atmos. D14C

  !-----------------------------------------------------------------------
  !  define array for holding flux-related quantities that need to be time-averaged
  !  this is necessary since the forcing routines are called before tavg flags
  !-----------------------------------------------------------------------

  real (r8), dimension(:,:,:,:), allocatable :: ECO_CISO_SFLUX_TAVG

  !-----------------------------------------------------------------------
  !  ciso_data_ind_d13c is the index for the D13C data for the
  !  current timestep
  !  Note that ciso_data_ind_d13c is always less than ciso_atm_d13c_data_nbval.
  !  To enable OpenMP parallelism, duplicating data_ind for each block
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(:), allocatable :: &
       ciso_data_ind_d13c, &      ! data index for D13C data
       ciso_data_ind_d14c         ! data index for D14C data

  !-----------------------------------------------------------------------
  !  average surface tracer value related variables
  !  used as reference value for virtual flux computations
  !-----------------------------------------------------------------------

  logical (log_kind), dimension(ecosys_ciso_tracer_cnt) :: &
       ciso_vflux_flag                ! which tracers get virtual fluxes applied

  integer (int_kind) :: &
       ciso_comp_surf_avg_flag        ! time flag id for computing average
                                      ! surface tracer values

  real (r8), dimension(ecosys_ciso_tracer_cnt) :: &
       ciso_surf_avg                  ! average surface tracer values

  !-----------------------------------------------------------------------
  !  timers
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       ecosys_ciso_interior_timer,                   &
       ecosys_ciso_sflux_timer

  !-----------------------------------------------------------------------
  !  scalar constants for 14C decay calculation
  !-----------------------------------------------------------------------

  real (r8), parameter :: c14_halflife_years = 5730.0_r8 !C14 half file
  real (r8) :: c14_lambda_inv_sec           ! Decay variable in seconds

  !---------------------------------------------------------------------
  !     Isotope standards
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  !     Isotope standards
  !---------------------------------------------------------------------

  real(r8) :: pi
  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine ecosys_ciso_init(init_ts_file_fmt, read_restart_filename, &
       tracer_d_module, TRACER_MODULE, lmarginal_seas, &
       errorCode)

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Initialize ecosys_ciso tracer module. This involves setting metadata, reading
    !  the module namelist, setting initial conditions, setting up forcing,
    !  and defining additional tavg variables.
    !---------------------------------------------------------------------

    implicit none

    character (*)                     , intent(in)    :: init_ts_file_fmt      ! format (bin or nc) for input file
    character (*)                     , intent(in)    :: read_restart_filename ! file name for restart file
    logical (kind=log_kind)           , intent(in)    :: lmarginal_seas               
    type (marbl_tracer_metadata_type) , intent(inout) :: tracer_d_module(:)    ! descriptors for each tracer
    real (r8)                         , intent(inout) :: TRACER_MODULE(:,:,:,:,:,:)
    integer (POP_i4)                  , intent(out)   :: errorCode

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ciso_mod:ecosys_ciso_init'

    type(marbl_tracer_read_type), dimension(ecosys_ciso_tracer_cnt) :: &
         ciso_tracer_init_ext              ! namelist variable for initializing tracers

    logical (log_kind) :: &
         default                      ! arg to init_time_flag

    integer (int_kind) :: &
         non_autotroph_ciso_tracer_cnt, & ! number of non-autotroph ecosystem tracers
         auto_ind,                        & ! autotroph functional group index
         n,                               & ! index for looping over tracers
         k,                               & ! index for looping over depth levels
         l,                               & ! index for looping over time levels
         ind,                             & ! tracer index for tracer name from namelist
         iblock,                          & ! index for looping over blocks
         nml_error                          ! namelist i/o error flag

    integer (int_kind) :: &
         freq_opt, freq,              & ! args for init_time_flag
         ciso_comp_surf_avg_freq_iopt,& ! choice for freq of comp_surf_avg
         ciso_comp_surf_avg_freq        ! choice for freq of comp_surf_avg

    logical (log_kind) :: &
         ciso_use_nml_surf_vals         ! do namelist surf values override values from restart file

    logical (log_kind) :: &
         ciso_lecovars_full_depth_tavg  ! should ecosystem vars be written full depth

    !  values to be used when comp_surf_avg_freq_opt==never
    real (r8) :: &
         ciso_surf_avg_di13c_const, &
         ciso_surf_avg_di14c_const

    ! ecosys_ciso_nml namelist
    namelist /ecosys_ciso_nml/ &
         ciso_init_ecosys_option, ciso_init_ecosys_init_file, &
         ciso_init_ecosys_init_file_fmt, ciso_tracer_init_ext, &
         ciso_comp_surf_avg_freq_opt, ciso_comp_surf_avg_freq,  &
         ciso_use_nml_surf_vals, ciso_surf_avg_di13c_const, &
         ciso_surf_avg_di14c_const, &
         ciso_lsource_sink, &
         ciso_lecovars_full_depth_tavg, &
         ciso_atm_d13c_opt, ciso_atm_d13c_const, ciso_atm_d13c_filename, &
         ciso_atm_d14c_opt, ciso_atm_d14c_const, ciso_atm_d14c_filename, &
         ciso_fract_factors, ciso_atm_model_year, ciso_atm_data_year

    character (char_len) :: &
         ecosys_ciso_restart_filename  ! modified file name for restart file
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  initialize name table
    !-----------------------------------------------------------------------

    errorCode = POP_Success

    !-----------------------------------------------------------------------
    !  initialize non-autotroph tracer_d values
    !  accumulate non_autotroph_ciso_tracer_cnt
    !-----------------------------------------------------------------------
    non_autotroph_ciso_tracer_cnt = 0

    tracer_d_module(di13c_ind)%short_name='DI13C'
    tracer_d_module(di13c_ind)%long_name='Dissolved Inorganic Carbon-13'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    tracer_d_module(do13c_ind)%short_name='DO13C'
    tracer_d_module(do13c_ind)%long_name='Dissolved Organic Carbon-13'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    tracer_d_module(zoo13C_ind)%short_name='zoo13C'
    tracer_d_module(zoo13C_ind)%long_name='Zooplankton Carbon-13'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    tracer_d_module(di14c_ind)%short_name='DI14C'
    tracer_d_module(di14c_ind)%long_name='Dissolved Inorganic Carbon-14'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    tracer_d_module(do14c_ind)%short_name='DO14C'
    tracer_d_module(do14c_ind)%long_name='Dissolved Organic Carbon-14'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    tracer_d_module(zoo14C_ind)%short_name='zoo14C'
    tracer_d_module(zoo14C_ind)%long_name='Zooplankton Carbon-14'
    non_autotroph_ciso_tracer_cnt = non_autotroph_ciso_tracer_cnt + 1

    do n = 1, non_autotroph_ciso_tracer_cnt
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
    end do

    !-----------------------------------------------------------------------
    !  confirm that ecosys_ciso_tracer_cnt is consistent with autotroph declarations
    !-----------------------------------------------------------------------

    n = non_autotroph_ciso_tracer_cnt
    do auto_ind = 1, autotroph_cnt
       n = n + 2 ! C13, C14 tracers
       if (autotrophs(auto_ind)%imp_calcifier .or. &
            autotrophs(auto_ind)%exp_calcifier) n = n + 2 ! Ca13CO3 & Ca14CO3 tracers
    end do

    if (ecosys_ciso_tracer_cnt /= n) then
       call document(subname, 'actual ecosys_ciso_tracer_cnt', ecosys_ciso_tracer_cnt)
       call document(subname, 'computed ecosys_ciso_tracer_cnt', n)
       call exit_POP(sigAbort, 'inconsistency between actual ecosys_ciso_tracer_cnt and computed ecosys_ciso_tracer_cnt')
    endif

    !-----------------------------------------------------------------------
    !  initialize autotroph tracer_d values and tracer indices
    !-----------------------------------------------------------------------

    n = non_autotroph_ciso_tracer_cnt + 1

    do auto_ind = 1, autotroph_cnt
       tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // '13C'
       tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon-13'
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%C13_ind = n
       n = n + 1

       tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // '14C'
       tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon-14'
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%C14_ind = n
       n = n + 1

       if (autotrophs(auto_ind)%imp_calcifier .or. autotrophs(auto_ind)%exp_calcifier) then
          tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Ca13CO3'
          tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Ca13CO3'
          tracer_d_module(n)%units      = 'mmol/m^3'
          tracer_d_module(n)%tend_units = 'mmol/m^3/s'
          tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%Ca13CO3_ind = n
          n = n + 1

          tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Ca14CO3'
          tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Ca14CO3'
          tracer_d_module(n)%units      = 'mmol/m^3'
          tracer_d_module(n)%tend_units = 'mmol/m^3/s'
          tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%Ca14CO3_ind = n
          n = n + 1
       else
          autotrophs(auto_ind)%Ca13CO3_ind = 0
          autotrophs(auto_ind)%Ca14CO3_ind = 0
       endif
    end do

    if (my_task == master_task) THEN
       write (stdout,*) '----- autotroph tracer indices -----'
       do auto_ind = 1, autotroph_cnt
          write (stdout,*) 'C13_ind('     , trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%C13_ind
          write (stdout,*) 'C14_ind('     , trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%C14_ind
          write (stdout,*) 'Ca13CO3_ind(' , trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Ca13CO3_ind
          write (stdout,*) 'Ca14CO3_ind(' , trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Ca14CO3_ind
          write (stdout,*) 'autotroph_cnt =',autotroph_cnt
       end do
       write (stdout,*) '------------------------------------'
    endif

    !-----------------------------------------------------------------------
    !  initialize ind_name_table
    !-----------------------------------------------------------------------

    do n = 1, ecosys_ciso_tracer_cnt
       ciso_ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
    end do

    !-----------------------------------------------------------------------
    !  default namelist settings
    !-----------------------------------------------------------------------

    ciso_init_ecosys_option                 = 'unknown'
    ciso_init_ecosys_init_file              = 'unknown'
    ciso_init_ecosys_init_file_fmt          = 'bin'
    do n = 1,ecosys_ciso_tracer_cnt
       ciso_tracer_init_ext(n)%mod_varname  = 'unknown'
       ciso_tracer_init_ext(n)%filename     = 'unknown'
       ciso_tracer_init_ext(n)%file_varname = 'unknown'
       ciso_tracer_init_ext(n)%scale_factor = c1
       ciso_tracer_init_ext(n)%default_val  = c0
       ciso_tracer_init_ext(n)%file_fmt     = 'bin'
    end do

    ciso_lsource_sink                       = .true.

    ciso_comp_surf_avg_freq_opt             = 'never'
    ciso_comp_surf_avg_freq                 = 1
    ciso_use_nml_surf_vals                  = .false.
    ciso_surf_avg_di13c_const               = 1944.0_r8
    ciso_surf_avg_di14c_const               = 1944.0_r8

    ciso_atm_d13c_opt                       = 'const'
    ciso_atm_d13c_const                     = -6.379_r8
    ciso_atm_d13c_filename                  = 'unknown'

    ciso_atm_d14c_opt                       = 'const'
    ciso_atm_d14c_const                     = 0.0_r8
    ciso_atm_d14c_filename(1)               = 'unknown'
    ciso_atm_d14c_filename(2)               = 'unknown'
    ciso_atm_d14c_filename(3)               = 'unknown'

    ciso_fract_factors                      = 'Rau'

    ciso_lecovars_full_depth_tavg           = .false.

    ciso_atm_model_year                     = 1
    ciso_atm_data_year                      = 1

    if (my_task == master_task) then
       open (nml_in, file=nml_filename, status='old',iostat=nml_error)
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       endif
       do while (nml_error > 0)
          read(nml_in, nml=ecosys_ciso_nml,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call document(subname, 'ecosys_ciso_nml not found')
       call exit_POP(sigAbort, 'ERROR : stopping in '/&
            &/ subname)
    endif

    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' ecosys_ciso:'
       write(stdout,blank_fmt)
       write(stdout,*) ' ecosys_ciso_nml namelist settings:'
       write(stdout,blank_fmt)
       write(stdout,ecosys_ciso_nml)
       write(stdout,blank_fmt)
       write(stdout,delim_fmt)
    endif

    !-----------------------------------------------------------------------
    !  broadcast all namelist variables
    !-----------------------------------------------------------------------

    call broadcast_scalar(ciso_init_ecosys_option, master_task)
    call broadcast_scalar(ciso_init_ecosys_init_file, master_task)
    call broadcast_scalar(ciso_init_ecosys_init_file_fmt, master_task)

    call broadcast_scalar(ciso_atm_d13c_opt, master_task)
    call broadcast_scalar(ciso_atm_d13c_const, master_task)
    call broadcast_scalar(ciso_atm_d13c_filename, master_task)

    call broadcast_scalar(ciso_atm_d14c_opt, master_task)
    call broadcast_scalar(ciso_atm_d14c_const, master_task)
    call broadcast_scalar(ciso_atm_d14c_filename(1), master_task)
    call broadcast_scalar(ciso_atm_d14c_filename(2), master_task)
    call broadcast_scalar(ciso_atm_d14c_filename(3), master_task)

    do n = 1,ecosys_ciso_tracer_cnt
       call broadcast_scalar(ciso_tracer_init_ext(n)%mod_varname, master_task)
       call broadcast_scalar(ciso_tracer_init_ext(n)%filename, master_task)
       call broadcast_scalar(ciso_tracer_init_ext(n)%file_varname, master_task)
       call broadcast_scalar(ciso_tracer_init_ext(n)%scale_factor, master_task)
       call broadcast_scalar(ciso_tracer_init_ext(n)%default_val, master_task)
       call broadcast_scalar(ciso_tracer_init_ext(n)%file_fmt, master_task)
    end do

    call broadcast_scalar(ciso_comp_surf_avg_freq_opt, master_task)
    call broadcast_scalar(ciso_comp_surf_avg_freq, master_task)
    call broadcast_scalar(ciso_use_nml_surf_vals, master_task)
    call broadcast_scalar(ciso_surf_avg_di13c_const, master_task)
    call broadcast_scalar(ciso_surf_avg_di14c_const, master_task)
    call broadcast_scalar(ciso_lsource_sink, master_task)
    call broadcast_scalar(ciso_lecovars_full_depth_tavg, master_task)
    call broadcast_scalar(ciso_fract_factors, master_task)
    call broadcast_scalar(ciso_atm_model_year, master_task)
    call broadcast_scalar(ciso_atm_data_year, master_task)

    !-----------------------------------------------------------------------
    !  set variables immediately dependent on namelist variables
    !-----------------------------------------------------------------------

    select case (ciso_comp_surf_avg_freq_opt)
    case ('never')
       ciso_comp_surf_avg_freq_iopt = freq_opt_never
    case ('nyear')
       ciso_comp_surf_avg_freq_iopt = freq_opt_nyear
    case ('nmonth')
       ciso_comp_surf_avg_freq_iopt = freq_opt_nmonth
    case default
       call document(subname, 'ciso_comp_surf_avg_freq_opt', ciso_comp_surf_avg_freq_opt)
       call exit_POP(sigAbort, 'unknown ciso_comp_surf_avg_freq_opt')
    end select

    call init_time_flag('ciso_ecosys_comp_surf_avg', ciso_comp_surf_avg_flag, &
         default=.false., freq_opt=ciso_comp_surf_avg_freq_iopt,  &
         freq=ciso_comp_surf_avg_freq, owner='ciso_ecosys_init')

    !-----------------------------------------------------------------------
    !  namelist consistency checking
    !-----------------------------------------------------------------------

    if (ciso_use_nml_surf_vals .and. ciso_comp_surf_avg_freq_iopt /= freq_opt_never) then
       call document(subname, 'ciso_use_nml_surf_vals', ciso_use_nml_surf_vals)
       call document(subname, 'ciso_comp_surf_avg_freq_opt', ciso_comp_surf_avg_freq_opt)
       call exit_POP(sigAbort, 'ciso_use_nml_surf_vals can only be .true. if ' /&
            &/ ' ciso_comp_surf_avg_freq_opt is never')
    endif

    !-----------------------------------------------------------------------
    !  initialize virtual flux flag array
    !-----------------------------------------------------------------------

    ciso_vflux_flag = .false.
    ciso_vflux_flag(di13c_ind) = .true.
    ciso_vflux_flag(di14c_ind) = .true.

    !-----------------------------------------------------------------------
    !  allocate and initialize LAND_MASK
    !-----------------------------------------------------------------------

    allocate( LAND_MASK(nx_block,ny_block,nblocks_clinic) )

    if (lmarginal_seas) then
       LAND_MASK = REGION_MASK /= c0
    else
       LAND_MASK = REGION_MASK > c0
    endif

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

       call rest_read_tracer_block(ciso_init_ecosys_init_file_fmt, &
            ecosys_ciso_restart_filename,   &
            tracer_d_module,           &
            TRACER_MODULE)

       if (ciso_use_nml_surf_vals) then

          ciso_surf_avg = c0
          ciso_surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
          ciso_surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
       else
          call extract_surf_avg(ciso_init_ecosys_init_file_fmt, &
               ecosys_ciso_restart_filename,   &
               ecosys_ciso_tracer_cnt,          &
               ciso_vflux_flag,                &
               ciso_ind_name_table,ciso_surf_avg)
       endif

       call eval_time_flag(ciso_comp_surf_avg_flag) ! evaluates time_flag(ciso_comp_surf_avg_flag)%value via time_to_do

       if (check_time_flag(ciso_comp_surf_avg_flag)) &
            call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
            TRACER_MODULE(:,:,1,:,curtime,:), &
            ecosys_ciso_tracer_cnt,     &
            ciso_vflux_flag,ciso_surf_avg)

    case ('file', 'ccsm_startup')
       call document(subname, 'ciso vars being read from separate files')

       call file_read_tracer_block(ciso_init_ecosys_init_file_fmt, &
            ciso_init_ecosys_init_file,     &
            tracer_d_module,           &
            ciso_ind_name_table,            &
            ciso_tracer_init_ext,           &
            TRACER_MODULE)

       if (n_topo_smooth > 0) then
          do n = 1, ecosys_ciso_tracer_cnt
             do k=1,km
                call fill_points(k,TRACER_MODULE(:,:,k,n,oldtime,:), &
                     errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(oldtime)')
                   return
                endif

                call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                     errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(newtime)')
                   return
                endif

             enddo
          enddo
       endif

       if (ciso_use_nml_surf_vals) then
          ciso_surf_avg = c0
          ciso_surf_avg(di13c_ind) = ciso_surf_avg_di13c_const
          ciso_surf_avg(di14c_ind) = ciso_surf_avg_di14c_const
       else
          call comp_surf_avg(TRACER_MODULE(:,:,1,:,oldtime,:), &
               TRACER_MODULE(:,:,1,:,curtime,:), &
               ecosys_ciso_tracer_cnt,           &
               ciso_vflux_flag,ciso_surf_avg)
       endif

    case default
       call document(subname, 'ciso_init_ecosys_option', ciso_init_ecosys_option)
       call exit_POP(sigAbort, 'unknown ciso_init_ecosys_option')

    end select

    !-----------------------------------------------------------------------
    !  apply land mask to tracers
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(iblock,n,k)
    do iblock=1,nblocks_clinic
       do n = 1,ecosys_ciso_tracer_cnt
          do k = 1,km
             where (.not. LAND_MASK(:,:,iblock) .or. k > KMT(:,:,iblock))
                TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
                TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
             end where
          end do
       end do
    enddo
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !  timer init
    !-----------------------------------------------------------------------

    call get_timer(ecosys_ciso_interior_timer, 'ECOSYS_CISO_INTERIOR', nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_ciso_sflux_timer, 'ECOSYS_CISO_SFLUX',1, distrb_clinic%nprocs)

    !-----------------------------------------------------------------------
    !  Define decay variable for DI14C, using earlier defined half-life of 14C
    !-----------------------------------------------------------------------

    c14_lambda_inv_sec = log(c2) / (c14_halflife_years * seconds_in_year)

    !-----------------------------------------------------------------------
    !  call other initialization subroutines
    !-----------------------------------------------------------------------

    call ecosys_ciso_init_sflux_tavg
    call ecosys_ciso_init_sflux

    !-----------------------------------------------------------------------
    !  set lfull_depth_tavg flag for short-lived ecosystem tracers
    !-----------------------------------------------------------------------

    tracer_d_module(zoo13C_ind   )%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
    tracer_d_module(zoo14C_ind   )%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%C13_ind
       tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%C14_ind
       tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%Ca13CO3_ind
       if (n > 0) then
          tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
       endif
       n = autotrophs(auto_ind)%Ca14CO3_ind
       if (n > 0) then
          tracer_d_module(n)%lfull_depth_tavg = ciso_lecovars_full_depth_tavg
       endif
    end do


    ! Set module variable
    pi  = 4.0_r8 * atan( 1.0_r8 )

  end subroutine ecosys_ciso_init

  !***********************************************************************

  subroutine ecosys_ciso_init_sflux_tavg

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  call define_tavg_field for nonstandard tavg fields
    !---------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: &
         subname = 'ecosys_ciso_mod:ecosys_ciso_init_tavg'

    integer (int_kind) :: &
         auto_ind,       & ! autotroph functional group index
         buf_len           ! how many surface flux fields are stored in ECO_CISO_SFLUX_TAVG

    character(char_len) :: &
         sname             ! short-name of tavg variable

    !-----------------------------------------------------------------------
    !  2D fields related to surface fluxes
    !-----------------------------------------------------------------------

    buf_len = 0

    !13C
    call define_tavg_field(tavg_CISO_DI13C_GAS_FLUX,'CISO_FG_13CO2',2,            &
         long_name='DI13C Surface Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_13CO2 = buf_len

    call define_tavg_field(tavg_CISO_DI13C_AS_GAS_FLUX,'CISO_FG_as_13CO2',2,            &
         long_name='DI13C Surface Air-Sea Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_as_13CO2 = buf_len

    call define_tavg_field(tavg_CISO_DI13C_SA_GAS_FLUX,'CISO_FG_sa_13CO2',2,            &
         long_name='DI13C Surface Sea-Air Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_sa_13CO2 = buf_len

    call define_tavg_field(tavg_CISO_d13C_GAS_FLUX,'CISO_FG_d13C',2,                 &
         long_name='D13C Surface GAS FLUX',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_d13C = buf_len

    call define_tavg_field(tavg_CISO_D13C_atm,'CISO_D13C_atm',2,                 &
         long_name='Atmospheric Delta 13C in permil',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_D13C_atm = buf_len

    call define_tavg_field(tavg_CISO_R13C_DIC_surf,'CISO_R13C_DIC_surf',2,                 &
         long_name='13C/12C ratio in total DIC',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_R13C_DIC_surf = buf_len

    call define_tavg_field(tavg_CISO_R13C_atm,'CISO_R13C_atm',2,                 &
         long_name='13C/12C ratio in atmosphere',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_R13C_atm = buf_len

    call define_tavg_field(tavg_CISO_DI13C_RIV_FLUX,'CISO_DI13C_RIV_FLUX',2,          &
         long_name='Flux of DI13C from rivers',         &
         units='nmol/cm^2/s', grid_loc='2110',        &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_DI13C_RIV_FLUX = buf_len

    call define_tavg_field(tavg_CISO_DO13C_RIV_FLUX,'CISO_DO13C_RIV_FLUX',2,          &
         long_name='Flux of DO13C from rivers',         &
         units='nmol/cm^2/s', grid_loc='2110',        &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_DO13C_RIV_FLUX = buf_len

    ! Fractionation (for 14C and 13C)

    call define_tavg_field(tavg_CISO_eps_aq_g_surf,'CISO_eps_aq_g_surf',2,              &
         long_name='Surface equilibrium fractionation (CO2_gaseous <-> CO2_aq)', &
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')

    buf_len = buf_len+1
    buf_ind_CISO_eps_aq_g_surf = buf_len

    call define_tavg_field(tavg_CISO_eps_dic_g_surf,'CISO_eps_dic_g_surf',2,              &
         long_name='Surface equilibrium fractionation between total DIC and gaseous CO2', &
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')

    buf_len = buf_len+1
    buf_ind_CISO_eps_dic_g_surf = buf_len

    !14C
    call define_tavg_field(tavg_CISO_DI14C_GAS_FLUX,'CISO_FG_14CO2',2,            &
         long_name='DI14C Surface Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_14CO2 = buf_len

    call define_tavg_field(tavg_CISO_DI14C_AS_GAS_FLUX,'CISO_FG_as_14CO2',2,            &
         long_name='DI14C Surface Air-Sea Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_as_14CO2 = buf_len

    call define_tavg_field(tavg_CISO_DI14C_SA_GAS_FLUX,'CISO_FG_sa_14CO2',2,            &
         long_name='DI14C Surface Sea-Air Gas Flux',            &
         units='mmol/m^3 cm/s', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_sa_14CO2 = buf_len

    call define_tavg_field(tavg_CISO_d14C_GAS_FLUX,'CISO_FG_d14C',2,                 &
         long_name='D14C Surface GAS FLUX',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_FG_d14C = buf_len

    call define_tavg_field(tavg_CISO_D14C_atm,'CISO_D14C_atm',2,                 &
         long_name='Atmospheric Delta 14C in permil',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_D14C_atm = buf_len

    call define_tavg_field(tavg_CISO_R14C_DIC_surf,'CISO_R14C_DIC_surf',2,                 &
         long_name='14C/12C ratio in total DIC',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_R14C_DIC_surf = buf_len

    call define_tavg_field(tavg_CISO_R14C_atm,'CISO_R14C_atm',2,                 &
         long_name='14C/12C ratio in atmosphere',&
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_R14C_atm = buf_len

    call define_tavg_field(tavg_CISO_DI14C_RIV_FLUX,'CISO_DI14C_RIV_FLUX',2,          &
         long_name='Flux of DI14C from rivers',         &
         units='nmol/cm^2/s', grid_loc='2110',        &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_DI14C_RIV_FLUX = buf_len

    call define_tavg_field(tavg_CISO_DO14C_RIV_FLUX,'CISO_DO14C_RIV_FLUX',2,          &
         long_name='Flux of DO14C from rivers',         &
         units='nmol/cm^2/s', grid_loc='2110',        &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_DO14C_RIV_FLUX = buf_len


    ! for debugging

    call define_tavg_field(tavg_CISO_GLOBAL_D14C,'CISO_GLOBAL_D14C',2,            &
         long_name='GLOBAL_D14C',            &
         units='permil', grid_loc='2110',      &
         coordinates='TLONG TLAT time')
    buf_len = buf_len+1
    buf_ind_CISO_GLOBAL_D14C = buf_len

    !-----------------------------------------------------------------------
    !  allocate array for holding flux-related quantities that need to be time-averaged
    !  this is necessary since the forcing routines are called before tavg flags
    !-----------------------------------------------------------------------

    allocate(ECO_CISO_SFLUX_TAVG(nx_block,ny_block,buf_len,max_blocks_clinic))
    ECO_CISO_SFLUX_TAVG = c0

  end subroutine ecosys_ciso_init_sflux_tavg

  !***********************************************************************

  subroutine ecosys_ciso_init_sflux

    ! !USES:
    use named_field_mod, only: named_field_get_index
    use registry, only: registry_match

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Initialize surface flux computations for the ecosys_ciso tracer module.
    !  Includes reading CO2 and D13C data from file if option file is used
    ! !REVISION HISTORY:
    !  same as module
    !---------------------------------------------------------------------

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

       !-----------------------------------------------------------------------
       !     READ in D13C data from file
       !-----------------------------------------------------------------------

    case('file')
       call ciso_read_atm_D13C_data

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

       !-----------------------------------------------------------------------
       !     READ in D14C data from files
       !-----------------------------------------------------------------------

    case('file')
       call ciso_read_atm_D14C_data

    case default
       call exit_POP(sigAbort, 'unknown ciso_atm_d14c_opt in ecosys_ciso_init_sflux')

    end select

  end subroutine ecosys_ciso_init_sflux

  !*****************************************************************************

  subroutine ecosys_ciso_set_sflux(&
       ecosys_surface_share, &
       SST,SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

    use marbl_share_mod, only: ecosys_surface_share_type
    use marbl_parms, only : R13c_std, R14c_std

    ! !INPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
         SST           ! sea surface temperature (C)

    real (r8), dimension(:,:,:,:), intent(in) :: &
         SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8), dimension(:,:,:,:), intent(inout) :: STF_MODULE

    type(ecosys_surface_share_type), intent(inout) :: ecosys_surface_share

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_ciso_mod:ecosys_ciso_set_sflux'

    logical (log_kind), save :: &
         first = .true.  ! Logical for first iteration test

    type (block) :: &
         this_block      ! block info for the current block

    integer (int_kind) :: &
         i,j,iblock,n, & ! loop indices
         auto_ind,     & ! autotroph functional group index
         mcdate,sec,   & ! date vals for shr_strdata_advance
         errorCode       ! errorCode from HaloUpdate call

    real (r8), dimension(nx_block) :: &
         DI13C_ROW,    & ! row of DI13C values for solver
         DI14C_ROW,    & ! row of DI13C values for solver
         DIC_ROW,      & ! row of DIC values for solver
         DCO2STAR_ROW, & ! DCO2STAR from solver
         CO2STAR_ROW,  & ! CO2STAR from solver
         PV_ROW


    real (r8), dimension(nx_block,ny_block) :: &
         D13C,          &   ! atm 13co2 value
         R13C_DIC_surf, &   ! 13C/12C ratio in DIC
         R13C_atm,      &   ! 13C/12C ratio in atmospheric CO2
         FLUX13,        &   ! gas flux of 13CO2 (nmol/cm^2/s)
         FLUX13_as,     &   ! air-to-sea gas flux of 13CO2 (nmol/cm^2/s)
         FLUX13_sa,     &   ! sea-to-air gas flux of 13CO2 (nmol/cm^2/s)
         D14C,          &   ! atm 14co2 value
         R14C_DIC_surf, &   ! 14C/12C ratio in total DIC
         R14C_atm,      &   ! 14C/12C ratio in atmospheric CO2
         FLUX14,        &   ! gas flux of 14CO2 (nmol/cm^2/s)
         FLUX14_as,     &   ! air-to-sea gas flux of 14CO2 (nmol/cm^2/s)
         FLUX14_sa,     &   ! sea-to-air gas flux of 14CO2 (nmol/cm^2/s)
         FLUX,          &   ! gas flux of CO2 (nmol/cm^2/s)
         FLUX_as,       &   ! air-to-sea gas flux of CO2 (nmol/cm^2/s)
         FLUX_sa            ! sea-to-air gas flux of CO2 (nmol/cm^2/s)


    real (r8), dimension(nx_block,ny_block) :: &
         eps_aq_g_surf,       & ! equilibrium fractionation (CO2_gaseous <-> CO2_aq)
         alpha_aq_g_surf,     & ! alpha_xxx_g_surf => eps = ( alpa -1 ) * 1000
         eps_dic_g_surf,      & ! equilibrium fractionation between total DIC and gaseous CO2
         alpha_dic_g_surf,    & ! alpha_xxx_g_surf => eps = ( alpa -1 ) * 1000
         eps_hco3_g_surf,     & ! equilibrium fractionation between bicarbonate and gaseous CO2
         eps_co3_g_surf,      & ! equilibrium fractionation between carbonate and gaseous CO2
         frac_co3,            & ! carbonate fraction fCO3 = [CO3--]/DIC
         alpha_aq_g_surf_14c, & ! for 14C, with fractionation being twice as large for 14C than for 13C
         alpha_dic_g_surf_14c   ! for 14C, with fractionation being twice as large for 14C than for 13C


    real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
         di13c_riv_flux,   & ! River input of DI13C
         do13c_riv_flux,   & ! River input of DO13C
         di14c_riv_flux,   & ! River input of DI14C
         do14c_riv_flux      ! River input of DO14C

    ! For making global average of atmospheric D14C
    real (r8), dimension(nx_block,ny_block) :: &
         WORK1, &! local work space
         TFACT   ! factor for normalizing sums

    integer (int_kind) :: &
         ib,ie,jb,je

    real (r8), dimension(max_blocks_clinic) :: &
         D14C_local_sums, & ! array for holding block sums when calculating global D14C
         TAREA_local_sums   ! array for holding block sums of TAREA when calculating global D14C

    real (r8) :: &
         D14C_sum_tmp,  & ! temp for local sum of D14C
         TAREA_sum_tmp    ! temp for local sum of TAREA

    real (r8) :: &
         D14C_glo_avg  ! global average D14C over the ocean, computed from current D14C field

    !-----------------------------------------------------------------------
    !     local parameters for 13C
    !     Zhang et al, 1995, Geochim. et Cosmochim. Acta, 59 (1), 107-114
    !-----------------------------------------------------------------------
    real(r8) :: &
         alpha_k, &             ! eps = ( alpa -1 ) * 1000
         alpha_k_14c            ! for 14C, with fractionation being twice as large for 14C than for 13C


    real(r8), parameter :: &
         eps_k     = -0.81_r8  ! kinetic fraction during gas
    ! transfert (per mil) (air-sea CO2
    ! exchange) at 21C, Zhang et al 1995,
    ! eps_k = -0.95 at 5C

    associate(                                                               &
         DIC_SURF_fields      => ecosys_surface_share%DIC_SURF_fields,       & ! IN/OUT
         CO3_SURF_fields      => ecosys_surface_share%CO3_SURF_fields,       & ! IN/OUT
         CO2STAR_SURF_fields  => ecosys_surface_share%CO2STAR_SURF_fields,   & ! IN/OUT
         DCO2STAR_SURF_fields => ecosys_surface_share%DCO2STAR_SURF_fields , & ! IN/OUT
         PV_SURF_fields       => ecosys_surface_share%PV_SURF_fields,        & ! IN/OUT
         dic_riv_flux_fields  => ecosys_surface_share%dic_riv_flux_fields,   & ! IN/OUT
         doc_riv_flux_fields  => ecosys_surface_share%doc_riv_flux_fields    & ! IN/OUT
         )
      !-----------------------------------------------------------------------

      call timer_start(ecosys_ciso_sflux_timer)
      !-----------------------------------------------------------------------

      if (check_time_flag(ciso_comp_surf_avg_flag))     &
           call comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR,&
           ecosys_ciso_tracer_cnt,     &
           ciso_vflux_flag,ciso_surf_avg)



      if (first) then
         allocate( ciso_data_ind_d13c(max_blocks_clinic) )
         allocate( ciso_data_ind_d14c(max_blocks_clinic) )
         ciso_data_ind_d13c = -1
         ciso_data_ind_d14c = -1
         first = .false.
      endif

      !-----------------------------------------------------------------------
      !  fluxes initially set to 0
      !-----------------------------------------------------------------------
      WORK1            = c0
      D14C_local_sums  = c0
      TAREA_local_sums = c0

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic
         STF_MODULE(:,:,:,iblock) = c0
      end do
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(iblock, j, this_block, ib, ie, jb, je, D13C, R13C_atm,   &
      !$OMP                     DI13C_ROW, R13C_DIC_surf, FLUX13, FLUX13_as, &
      !$OMP                     FLUX13_sa, D14C, R14C_atm, DI14C_ROW, &
      !$OMP                     R14C_DIC_surf, FLUX14, FLUX14_as, FLUX14_sa, &
      !$OMP                     DIC_ROW, CO2STAR_ROW, DCO2STAR_ROW, FLUX, &
      !$OMP                     FLUX_as, FLUX_sa, PV_ROW, eps_aq_g_surf, &
      !$OMP                     eps_dic_g_surf, frac_co3, alpha_k, &
      !$OMP                     alpha_aq_g_surf, alpha_dic_g_surf, &
      !$OMP                     alpha_k_14c, alpha_aq_g_surf_14c, &
      !$OMP                     alpha_dic_g_surf_14c,TFACT,WORK1)
      !-----------------------------------------------------------------------

      do iblock = 1, nblocks_clinic

         !-----------------------------------------------------------------------
         !  Set D13C and D14C (constant or from files read in _init) and put on global grid
         !-----------------------------------------------------------------------
         select case (ciso_atm_d13c_opt)

         case ('const')
            D13C = ciso_atm_d13c_const

         case ('file')
            call ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c(iblock),D13C)

         case default
            call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt in ecosys_ciso_set_sflux')

         end select

         select case (ciso_atm_d14c_opt)

         case ('const')
            D14C = ciso_atm_d14c_const

         case ('file')
            call ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c(iblock),D14C)

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

         where (LAND_MASK(:,:,iblock))
            TFACT = TAREA(:,:,iblock)
         elsewhere
            TFACT = 0.0_r8
         endwhere

         WORK1 = D14C * TFACT
         D14C_local_sums(iblock) = sum(WORK1(ib:ie,jb:je))
         TAREA_local_sums(iblock) = sum(TFACT(ib:ie,jb:je))


         !-----------------------------------------------------------------------
         !     initialize R13C_atm  and R14C_atm
         !-----------------------------------------------------------------------

         R13C_atm = R13C_std * ( c1 + D13C / c1000 )

         R14C_atm = R14C_std * ( c1 + D14C / c1000 )

         !-----------------------------------------------------------------------
         !     compute 13C02 flux, based on CO2 flux calculated in ecosystem model
         !     Zhang et al, 1995, Geochim. et Cosmochim. Acta, 59 (1), 107-114
         !-----------------------------------------------------------------------

         do j = 1,ny_block
            !-----------------------------------------------------------------------
            !     compute R13C_DIC in surface ocean (assuming that DIC_ROW is 12C)
            !-----------------------------------------------------------------------
            DI13C_ROW = p5*(SURF_VALS_OLD(:,j,di13c_ind,iblock) + &
                 SURF_VALS_CUR(:,j,di13c_ind,iblock))
            DI14C_ROW = p5*(SURF_VALS_OLD(:,j,di14c_ind,iblock) + &
                 SURF_VALS_CUR(:,j,di14c_ind,iblock))

            DIC_ROW   = DIC_SURF_fields(:,j,iblock)

            where ( DIC_ROW /= c0 )
               R13C_DIC_surf(:,j) = DI13C_ROW /DIC_ROW
            elsewhere
               R13C_DIC_surf(:,j) = c0
            endwhere

            where ( DIC_ROW /= c0 )
               R14C_DIC_surf(:,j) = DI14C_ROW / DIC_ROW
            elsewhere
               R14C_DIC_surf(:,j) = c0
            endwhere


            !-----------------------------------------------------------------------
            !     individal discrimination factor of each species with respect to
            !     gaseous CO2, temperature dependent, based on Zhang et al. 95
            !-----------------------------------------------------------------------
            eps_aq_g_surf(:,j)   = 0.0049_r8 * SST(:,j,iblock) - 1.31_r8
            !!         eps_hco3_g_surf(:,j) = -0.141_r8  * SST(:,j,iblock) + 10.78_r8
            !!         eps_co3_g_surf(:,j)  = -0.052_r8  * SST(:,j,iblock) + 7.22_r8


            !-----------------------------------------------------------------------
            !     compute the equilibrium discrimination factor between DIC and
            !     gaseous CO2
            !-----------------------------------------------------------------------
            !     solution 1 : from individual species.
            !     Not used as Zhang et al. 95
            !     concluded that eps_dic_g_surf can not be calculated from the sum of
            !     the three individual species
            !----------------------------------------------------------------------
            !         eps_dic_g_surf(:,j) = eps_aq_g_surf(:,j) + eps_hco3_g_surf(:,j) &
            !                                + eps_co3_g_surf(:,j)
            !-----------------------------------------------------------------------
            !     solution 2: function of T and carbonate fraction (frac_co3)
            !     Using this one, which is based on the empirical relationship from
            !     the measured e_dic_g_surf of Zhang et al. 1995
            !---------------------------------------------------------------------

            where (.not. LAND_MASK(:,j,iblock))
               frac_co3(:,j) = c0
            elsewhere
               frac_co3(:,j) = CO3_SURF_fields(:,j,iblock) / DIC_ROW
            end where

            eps_dic_g_surf(:,j)  = 0.014_r8 * SST(:,j,iblock) * frac_co3(:,j) - &
                 0.105_r8 * SST(:,j,iblock) + 10.53_r8

            !-----------------------------------------------------------------------
            !     compute alpha coefficients from eps :  eps = ( alpha -1 ) * 1000
            !     => alpha = 1 + eps / 1000
            !-----------------------------------------------------------------------

            alpha_k               = c1 + eps_k               / c1000
            alpha_aq_g_surf(:,j)  = c1 + eps_aq_g_surf(:,j)  / c1000
            alpha_dic_g_surf(:,j) = c1 + eps_dic_g_surf(:,j) / c1000

            ! Fractionation is twice as large for 14C than for 13C, so eps needs to be multiplied by 2 for 14C
            alpha_k_14c               = c1 + eps_k * 2.0_r8              / c1000
            alpha_aq_g_surf_14c(:,j)  = c1 + eps_aq_g_surf(:,j) *2.0_r8  / c1000
            alpha_dic_g_surf_14c(:,j) = c1 + eps_dic_g_surf(:,j) *2.0_r8 / c1000

            !-----------------------------------------------------------------------
            !     compute 13C flux and C flux
            !-----------------------------------------------------------------------
            CO2STAR_ROW  = CO2STAR_SURF_fields(:,j,iblock)
            DCO2STAR_ROW = DCO2STAR_SURF_fields(:,j,iblock)
            PV_ROW       = PV_SURF_fields(:,j,iblock)

            FLUX13(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                 (( CO2STAR_ROW + DCO2STAR_ROW ) * R13C_atm(:,j) - &
                 CO2STAR_ROW * R13C_DIC_surf(:,j) / alpha_dic_g_surf(:,j) )

            FLUX14(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                 (( CO2STAR_ROW + DCO2STAR_ROW ) * R14C_atm(:,j) - &
                 CO2STAR_ROW * R14C_DIC_surf(:,j) / alpha_dic_g_surf_14C(:,j) )

            FLUX(:,j)   = PV_ROW * DCO2STAR_ROW



            !-----------------------------------------------------------------------
            !     compute fluxes in and out
            !-----------------------------------------------------------------------

            FLUX_as(:,j)   = PV_ROW * ( DCO2STAR_ROW + CO2STAR_ROW )
            FLUX_sa(:,j)   = PV_ROW * CO2STAR_ROW

            FLUX13_as(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                 (( CO2STAR_ROW + DCO2STAR_ROW ) * R13C_atm(:,j))

            FLUX13_sa(:,j) = PV_ROW * alpha_k * alpha_aq_g_surf(:,j) * &
                 ( CO2STAR_ROW * R13C_DIC_surf(:,j) / alpha_dic_g_surf(:,j) )

            FLUX14_as(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                 (( CO2STAR_ROW + DCO2STAR_ROW ) * R14C_atm(:,j))

            FLUX14_sa(:,j) = PV_ROW * alpha_k_14c * alpha_aq_g_surf_14c(:,j) * &
                 ( CO2STAR_ROW * R14C_DIC_surf(:,j) / alpha_dic_g_surf_14c(:,j) )

            !-----------------------------------------------------------------------
            !     end of 13C computation for gass exchange
            !-----------------------------------------------------------------------

         end do !j loop

         !-----------------------------------------------------------------------
         !     Adding 13C FLux to total DI13C
         !-----------------------------------------------------------------------

         STF_MODULE(:,:,di13c_ind,iblock) = STF_MODULE(:,:,di13c_ind,iblock) + FLUX13
         STF_MODULE(:,:,di14c_ind,iblock) = STF_MODULE(:,:,di14c_ind,iblock) + FLUX14

         !-----------------------------------------------------------------------
         !    Tavg variables
         !-----------------------------------------------------------------------

         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_13CO2,iblock)    = FLUX13
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_13CO2,iblock) = FLUX13_as
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_13CO2,iblock) = FLUX13_sa

         where ( FLUX /= c0 )
            ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d13C,iblock) = ( FLUX13 / FLUX  &
                 / R13C_std - c1 ) * c1000
         elsewhere
            ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d13C,iblock) = c0
         endwhere

         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_DIC_surf,iblock)  = R13C_DIC_surf
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R13C_atm,iblock)  = R13C_atm
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D13C_atm,iblock) = D13C
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_aq_g_surf,iblock) = eps_aq_g_surf
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_eps_dic_g_surf,iblock) = eps_dic_g_surf
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_14CO2,iblock) = FLUX14
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_as_14CO2,iblock) = FLUX14_as
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_sa_14CO2,iblock) = FLUX14_sa
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_DIC_surf,iblock) = R14C_DIC_surf
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_R14C_atm,iblock) = R14C_atm
         ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_D14C_atm,iblock) = D14C

         where ( FLUX /= c0 )
            ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d14C,iblock) = ( FLUX14 / FLUX &
                 / R14C_std - c1 ) * c1000
         elsewhere
            ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_FG_d14C,iblock) = c0
         endwhere

      enddo ! end of i-lopp
      !-----------------------------------------------------------------------
      !$OMP END PARALLEL DO
      !-----------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! River input of isotopic DIC and DOC.
      ! River input of BGC tracers in ecosys_mod is currently constant and from file
      ! So the isotopic carbon input is also done very simplified with one value
      ! globally, even though data shows it should vary from river to river.
      !
      ! Using constant delta values of
      ! D13C=-10 permil for DIC (Mook 1986, Raymond et al 2004)
      ! D13C=-27.6 permil for DOC (Raymond et al 2004)
      ! D14C=-50 permil for DOC (Raymond et al 2004), Gruber et al
      ! D14C= atmos_D14C - 50 permil for DIC (based on very few data points and 
      !       discussion with N. Gruber)
      !-------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      !  Make average global DI14C
      !-----------------------------------------------------------------------------

      D14C_sum_tmp  = sum(D14C_local_sums)
      TAREA_sum_tmp = sum(TAREA_local_sums)


      D14C_glo_avg  = global_sum(D14C_sum_tmp,distrb_clinic) / &
           global_sum(TAREA_sum_tmp,distrb_clinic)


      di13c_riv_flux = dic_riv_flux_fields * (-10.0_r8/c1000 +c1) * R13C_std
      di14c_riv_flux = dic_riv_flux_fields * ((D14C_glo_avg - 50.0_r8)/c1000 +c1) * R14C_std

      do13c_riv_flux = doc_riv_flux_fields * (-27.6_r8/c1000 +c1) * R13C_std
      do14c_riv_flux = doc_riv_flux_fields * (-50.0_r8/c1000 +c1) * R14C_std

      STF_MODULE(:,:,di13c_ind,:) = STF_MODULE(:,:,di13c_ind,:) + di13c_riv_flux
      STF_MODULE(:,:,do13c_ind,:) = STF_MODULE(:,:,do13c_ind,:) + do13c_riv_flux

      STF_MODULE(:,:,di14c_ind,:) = STF_MODULE(:,:,di14c_ind,:) + di14c_riv_flux
      STF_MODULE(:,:,do14c_ind,:) = STF_MODULE(:,:,do14c_ind,:) + do14c_riv_flux

      ! write to tavg
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI13C_RIV_FLUX,:) = di13c_riv_flux
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO13C_RIV_FLUX,:) = do13c_riv_flux

      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DI14C_RIV_FLUX,:) = di14c_riv_flux
      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_DO14C_RIV_FLUX,:) = do14c_riv_flux

      ECO_CISO_SFLUX_TAVG(:,:,buf_ind_CISO_GLOBAL_D14C,:) = D14C_glo_avg

      call timer_stop(ecosys_ciso_sflux_timer)

    end associate

  end subroutine ecosys_ciso_set_sflux

  !***********************************************************************

  subroutine ciso_read_atm_D13C_data

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

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: sub_name = 'ecosys_ciso_mod:ciso_read_atm_D13C_data'

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
       allocate(ciso_atm_d13c_data(ciso_atm_d13c_data_nbval))
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


    call broadcast_array(ciso_atm_d13c_data, master_task)
    call broadcast_array(ciso_atm_d13c_data_yr, master_task)

  end subroutine ciso_read_atm_D13C_data

  !***********************************************************************

  subroutine ciso_read_atm_D14C_data

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
    allocate(ciso_atm_d14c_data(ciso_atm_d14c_data_nbval,3))

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

    call broadcast_array(ciso_atm_d14c_data, master_task)
    call broadcast_array(ciso_atm_d14c_data_yr, master_task)

  end subroutine ciso_read_atm_D14C_data

  !***********************************************************************

  subroutine ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c, D13C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of d13c when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate data values to current model time step
    !  2. Spatial patern of D13Cis the same everywhere (90 S - 90 N)
    !-----------------------------------------------------------------------

    implicit none

    ! note that ciso_data_ind_d13c is always strictly less than the length
    ! of the data and is initialized to -1 before the first call

    integer (int_kind) , intent(in) :: iblock                   ! block index
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

  end subroutine ciso_comp_varying_D13C

  !***********************************************************************

  subroutine ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c, D14C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of CO2 when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate hemispheric values to current time step
    !  2. Make global field of D14C, determined by:
    !   -Northern Hemisphere value is used for 20N - 90 N
    !   -Southern Hemisphere value is used for 20 S - 90 S
    !   -Equator value is used for 20 S- 20 N

    implicit none

    ! !INPUT PARAMETERS:

    integer (int_kind) :: &
         iblock          ! block index

    integer (int_kind) :: &
         ciso_data_ind_d14c   ! data_ind_d14c is the index into data for current timestep,
                              !  note that data_ind is always strictly less than the length of D14C data
                              !  and is initialized to -1 before the first call

    ! !OUTPUT PARAMETERS:

    real (r8), dimension(nx_block,ny_block), intent(out) :: &
         D14C            ! atmospheric delta C14 in permil on global grid

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

  end subroutine ciso_comp_varying_D14C

  !*****************************************************************************

  subroutine ecosys_ciso_tavg_forcing()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Accumulate non-standard forcing related tavg variables.
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: iblock  ! block loop index
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Do components of flux calculations saved in hack array
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(iblock)

    do iblock = 1,nblocks_clinic
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_13CO2,iblock)       , tavg_CISO_DI13C_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_as_13CO2,iblock)    , tavg_CISO_DI13C_AS_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_sa_13CO2,iblock)    , tavg_CISO_DI13C_SA_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_d13C,iblock)        , tavg_CISO_d13C_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_R13C_DIC_surf,iblock)  , tavg_CISO_R13C_DIC_surf,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_R13C_atm,iblock)       , tavg_CISO_R13C_atm,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_D13C_atm,iblock)       , tavg_CISO_D13C_atm,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_eps_aq_g_surf,iblock)  , tavg_CISO_eps_aq_g_surf,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_eps_dic_g_surf,iblock) , tavg_CISO_eps_dic_g_surf,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_14CO2,iblock)       , tavg_CISO_DI14C_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_as_14CO2,iblock)    , tavg_CISO_DI14C_AS_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_sa_14CO2,iblock)    , tavg_CISO_DI14C_SA_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_R14C_DIC_surf,iblock)  , tavg_CISO_R14C_DIC_surf,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_R14C_atm,iblock)       , tavg_CISO_R14C_atm,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_D14C_atm,iblock)       , tavg_CISO_D14C_atm,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_FG_d14C,iblock)        , tavg_CISO_d14C_GAS_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_DI13C_RIV_FLUX,iblock)      , tavg_CISO_DI13C_RIV_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_DO13C_RIV_FLUX,iblock)      , tavg_CISO_DO13C_RIV_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_DI14C_RIV_FLUX,iblock)      , tavg_CISO_DI14C_RIV_FLUX,iblock,1)
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_DO14C_RIV_FLUX,iblock)      , tavg_CISO_DO14C_RIV_FLUX,iblock,1)
       !debugging
       call accumulate_tavg_field(ECO_CISO_SFLUX_TAVG(: ,:,buf_ind_CISO_GLOBAL_D14C,iblock)    , tavg_CISO_GLOBAL_D14C,iblock,1)
    end do

    !$OMP END PARALLEL DO

  end subroutine ecosys_ciso_tavg_forcing

  !*****************************************************************************

  function ecosys_ciso_tracer_ref_val(ind)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  return reference value for tracers using virtual fluxes
    !-----------------------------------------------------------------------

    implicit none

    integer (int_kind), intent(in) :: ind

    real (r8) :: ecosys_ciso_tracer_ref_val
    !-----------------------------------------------------------------------

    if (ciso_vflux_flag(ind)) then
       ecosys_ciso_tracer_ref_val = ciso_surf_avg(ind)
    else
       ecosys_ciso_tracer_ref_val = c0
    endif
    
  end function ecosys_ciso_tracer_ref_val

  !*****************************************************************************

  subroutine ecosys_ciso_write_restart(restart_file, action)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files
    !-----------------------------------------------------------------------

    use constants, only : char_blank
    use constants, only : field_loc_center
    use constants, only : field_type_scalar

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
       do n=1,ecosys_ciso_tracer_cnt
          if (ciso_vflux_flag(n)) then
             short_name = 'surf_avg_' // ciso_ind_name_table(n)%name
             call add_attrib_file(restart_file, trim(short_name), ciso_surf_avg(n))
          endif
       end do
    endif
    
  end subroutine ecosys_ciso_write_restart

end module ecosys_ciso_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


