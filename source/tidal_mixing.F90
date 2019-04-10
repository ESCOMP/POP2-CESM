!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module tidal_mixing

!BOP
! !MODULE: tidal_mixing
!
! !DESCRIPTION:
! This module computes the time-independent part of the tidally 
! driven mixing coefficients. 
!
! !DESCRIPTION:
! This module contains tidal-mixing method and initialization routines,
! including multiple options for tidal energy formulations and tidal mixing
! parameterization schemes. 
! 
! The user can select from three main tidal-mixing methods and four
! tidal energy initialization files. Additionally, the user can opt
! to modulate the 3D Egbert & Ray and Green & Nycander tidal energy fields
! with an 18.6-year lunar cycle and use with either the Schmittner (3D)
! method or the Jayne (2D) method. 
!
! Note that the tidal energy field in the Polzin/Melet method evolves in
! time as a function of N(bottom); this form of tidal energy may optionally
! be selected for use in the Jayne method.
! 
!
! Tidal Mixing Methods:
! ====================
! 
!   1) tidal_mixing_method_jayne
!      _________________________
!
!      Jayne, S. R., and L. C. St. Laurent, 2001: Parameterizing tidal
!        dissipation over rough topography. Geophys. Res. Lett., 
!        v28, 811-814.
!
!      Simmons, H. L., S. R. Jayne, L. C. St. Laurent, and
!        A. J. Weaver, 2004: Tidally driven mixing in a numerical
!        model of the ocean general circulation. Ocean Modelling,
!        vol 6, 245-263.
!
!      Jayne, Steven R., 2009: The Impact of Abyssal Mixing
!        Parameterizations in an Ocean General Circulation Model.
!        JPO, vol 39, 1756-1775.
!
!   2) tidal_mixing_method_polzin (based on Melet 2013)
!      ________________________________________________
!
!      Polzin, K. L., 2009:  An abyssal recipe. Ocean Modelling,
!        vol 30, 298-309
!
!      Melet, A. et al, 2013: Sensitivity of the ocean state to the
!        vertical distribution of the internal-tide-driven mixing. 
!        J. Phys Oceanography, vol 43, 602-615
!
!   3) tidal_mixing_method_schmittner
!      ______________________________
!
!      Schmittner, A. and G.D. Egbert, 2014: An improved parameterization
!        of tidal mixing for ocean models. Geosci. Model Dev., 7, 211-224, 2014
!
! Tidal Energy Input Files:
! ========================
! 
!   1) tidal_energy_arbic (Arbic 2014)
!      __________________
!
!
!   2) tidal_energy_jayne (Jayne 2009)
!      __________________
!
!
!   3) tidal_energy_ERO3 (Egbert and Ray 2003)
!      _________________
!
!
!   4) tidal_energy_GN13 (Green and Nycander 2013)
!      _________________
!
!   5) tidal_energy_LGM0 (Wilmes)  !active research -- do not use
!      _________________
!
!   6) tidal_energy_LGMi5g21 (Wilmes)  !active research -- do not use
!      _____________________
!
!   7) tidal_energy_LGMi6g21 (Wilmes)  !active research -- do not use
!      _____________________
!
!
! The following combinations of tidal energy files and options are allowed:
! 
!                    MIXING METHOD -->  Jayne   Polzin   Schmittner
!      ENERGY FILE  /
!          |
!          V
!        Arbic 2014                      yes     yes       no
!    
!        Jayne 2009                      yes     yes       no
! 
!        Egbert & Ray 2003               yes(2D) yes(2D)   yes(3D)
! 
!        Green & Nycander 2013           yes(2D) yes(2D)   yes(3D)
! 
!        LGM0 Wilmes                     yes(2D) yes(2D)   yes(3D)
! 
!        LGMi5g21 Wilmes                 yes(2D) yes(2D)   yes(3D)
! 
!        LGMi6g21 Wilmes                 yes(2D) yes(2D)   yes(3D)
! 
! Notes on initialization flow control:
!    1) init_tidal_mixing1    (tidal_mixing.F90) read tidal-mixing namelist variables
!       ------------------
!    2) init_ts                (initial.F90)
!       -------
!       |
!       +-> init_restart       (restart.F90) needs ltidal_lunar_cycle for LNC restart
!           ------------      
!    3) init_vertical_mix      (vertical_mix.F90)
!       -----------------
!       |
!       +-> init_tidal_mixing2 (tidal_mixing.F90) complete initialization of tidally driven vert. mixing
!          ------------------
!       |
!       +-> init_tidal_ts      (tidal_mixing.F90) initialize LNC; on restart uses info from init_restart
!           -------------
!       |
!       +-> init_vmix_kpp      (vmix_kpp.F90) needs info from init_tidal_mixing1,init_tidal_mixing2,init_restart
!           -------------
!       |
!       +-> cvmix_compute_Schmittner_invariant (cvmix_tidal.F90)  needs nlev from init_vmix_kpp
!           ----------------------------------
!       |
!       +-> cvmix_compute_Schmittner           (cvmix_tidal.F90) needs q*E(x,y,z) from init_tidal_mixing2
!           ------------------------
! 
! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use blocks
   use broadcast
   use communicate
   use constants
   use cvmix_tidal
   use domain_size
   use domain
   use exit_mod
   use gather_scatter
   use global_reductions
   use grid
   use io
   use io_tools
   use io_types
   use passive_tracer_tools
   use registry
   use shr_sys_mod
   use tavg
   use time_management
   use timers

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
 
   public :: init_tidal_mixing1
   public :: init_tidal_mixing2
   public :: init_tidal_ts
   public :: tidal_compute_diff
   public :: tidal_compute_diff_polzin
   public :: tidal_N2_integral
   public :: tidal_zstarp_inv
   public :: tidal_accumulate_tavg
   public :: tidal_accumulate_tavg_once
   public :: tidal_ts_driver
   public :: tidal_min_regions_set

! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:,:,:,:), allocatable, public ::  &
      TIDAL_DIFF,                            &! diffusivity due to tidal mixing 
      TIDAL_N2,                              &! buoyancy with no floor
      TIDAL_N2_eps                            ! buoyancy with epsilon floor
   real (r8), dimension(:,:,:,:), allocatable, public ::  &
      TIDAL_COEF_3D            ! time-independent coefficients for tidal-mixing methods

   ! FIXME(mnl,2016-01) -- moved from init -> module level for use by cvmix
   !                       branch of vmix_kpp; refactoring CVMix will let
   !                       this go back to init (and get deallocated there)
   !  njn01: this comment now applies to TIDAL_ENERGY_FLUX_2D

   public :: TIDAL_ENERGY_FLUX_2D
   public :: TIDAL_QE_2D
   public :: TIDAL_QE_3D
   public :: tavg_TIDAL_DIFF
   public :: tavg_TIDAL_COEF_3D
   public :: tavg_REGION_BOX3D

   interface tidal_min_regions_set
     module procedure tidal_min_regions_set_cvmixT
     module procedure tidal_min_regions_set_cvmixF
   end interface

   interface tidal_N2_integral
     module procedure tidal_N2_integral_clmn
     module procedure tidal_N2_integral_3D
   end interface

   interface tidal_zstarp_inv
     module procedure tidal_zstarp_inv_clmn
     module procedure tidal_zstarp_inv_3D
   end interface

   interface tidal_compute_diff_polzin
     module procedure tidal_compute_diff_polzin_clmn
     module procedure tidal_compute_diff_polzin_2D
   end interface

!-----------------------------------------------------------------------
!
!  choices for tidal mixing method 
!
!-----------------------------------------------------------------------

   integer (int_kind), public, parameter ::   &! integer ids for tidal mixing method
      tidal_mixing_method_jayne          = 1, &! Jayne & St. Laurent 2001; Simmons et al 2004
      tidal_mixing_method_polzin         = 2, &! Polzin, K.L. 2009; Melet et al 2013
      tidal_mixing_method_schmittner     = 3   ! Schmittner & Egbert 2013 

   integer (int_kind), public :: tidal_mixing_method_itype

!-----------------------------------------------------------------------
!
!  choices for tidal energy input file
!
!-----------------------------------------------------------------------

   integer (int_kind), public, parameter ::  &! integer ids for tidal energy input file
      tidal_energy_arbic                = 1, &! Arbic 2014
      tidal_energy_jayne                = 2, &! Jayne 2009
      tidal_energy_ER03                 = 3, &! Egbert & Ray 2003 
      tidal_energy_GN13                 = 4, &! Green & Nycander 2013
      tidal_energy_LGM0                 = 5, &! LGM 0ybp Wilmes      active research do not use without author permission
      tidal_energy_LGMi5g21             = 6, &! LGM 21kybpi5g Wilmes active research do not use without author permission
      tidal_energy_LGMi6g21             = 7   ! LGM 21kybpig6 Wilmes active research do not use without author permission

   integer (int_kind), public             :: &
      tidal_energy_itype                      ! tidal energy input file selected


!-----------------------------------------------------------------------
!
!  tidal maxima and minima;  LAT and LON for regional boxes
!
!-----------------------------------------------------------------------

   real (r8),public :: tidal_mix_max                                      ! upper limit for TIDAL_DIFF 

   integer (int_kind), public, parameter     :: max_tidal_min_regions = 9 ! max number of regions
   integer (int_kind), public                :: num_tidal_min_regions     ! number of regions
   logical (log_kind), public                :: ltidal_min_regions        ! apply min TIDAL_DIFF in regions
   character (char_len), dimension(max_tidal_min_regions),public ::   &
                                                tidal_min_regions_name
   real (r8),dimension(max_tidal_min_regions),public ::  &
                                                tidal_min_values,      &  ! minimum tidal-mixing values in specified regions
                                                tidal_TLATmin_regions, &  ! lower TLAT limit in specified regions
                                                tidal_TLATmax_regions ,&  ! upper TLAT limit in specified regions
                                                tidal_TLONmin_regions, &  ! lower TLON limit in specified regions
                                                tidal_TLONmax_regions     ! upper TLON limit in specified regions
   integer (int_kind),dimension(max_tidal_min_regions),public ::  &
                                                tidal_min_regions_klevels ! apply min over lowest k levels
   integer (int_kind),dimension(:,:,:,:), allocatable  :: REGION_BOX3D

   logical (log_kind), public :: ltidal_max     ! cap TIDAL_DIFF at tidal_mix_max

   character (char_len) :: tidal_mixing_method_choice  ! choice of tidal mixing method option
   character (char_len) :: tidal_energy_choice         ! choice of tidal energy input file option

!-----------------------------------------------------------------------
!
!  module logical control
!
!-----------------------------------------------------------------------

   logical (log_kind)         :: lccsm_control_compatible           !backwards compatibility with ccsm4 control
   logical (log_kind), public :: ltidal_schmittner_socn             !apply Schmittner southern ocean mod
   logical (log_kind), public :: ltidal_stabc                       !apply tidal diffusion stability controls
   logical (log_kind), public :: ltidal_lunar_cycle                 !activate 18.6-year lunar nodal cycle (LNC)
   logical (log_kind), public :: ltidal_lunar_cycle_read_restart    !read LNC info from restart file
   logical (log_kind)         :: ltidal_lunar_cycle_print = .false. !print output diagnostics
   logical (log_kind), public :: ltidal_melet_plot  = .false.       !create melet figure plot
   logical (log_kind), public :: ltidal_mixing                      !activate tidal mixing parameterization


!-----------------------------------------------------------------------
!
!  support for reading tidal-energy input data files
!
!-----------------------------------------------------------------------

   character (char_len_long) :: tidal_energy_file    ! tidal energy input file (generic )
   character (char_len_long) :: tidal_vars_file_polz ! tidal vars input file polzin

   character (char_len) ::        &! 'bin' or 'nc'
      tidal_energy_file_fmt,      &! tidal energy input file format generic
      tidal_vars_file_fmt_polz     ! tidal vars input file format polzin

   type (datafile) ::             &
      tidal_mixing_file_in         ! io file descriptor

   type (io_dim) ::               &
      i_dim,                      &! dimension descriptors for horiz dims
      j_dim,                      &
      k_dim                        ! dimension descriptor for vertical dim

   type (io_field_desc) ::        &
      H2_P_D,                     &! field descriptor for input H2_P flux
      TIDAL_ENERGY_FLUX_D,        &! field descriptor for input energy flux
      TCK1_D,                     &!                  for diurnal lunar
      TCM2_D,                     &!                  for semidiurnal lunar 
      TCO1_D,                     &!                  for diurnal solar
      TCS2_D,                     &!                  for semidirunal solar 
      U_P_D                        !                  for input U_P

!-----------------------------------------------------------------------
!
!  general module variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::  &
      TIDAL_QE_3D,            &! q*E(x,y,z) g/s^3 
      VERTICAL_FUNC,          &
      VERTICAL_NUM        

   real (r8), dimension(:,:,:), allocatable ::  &
      TIDAL_QE_2D,            &! q*E(x,y) g/s^3 
      TIDAL_COEF_2D            ! time-independent coefficients for tidal-mixing methods 

   real (r8), dimension(:,:,:), target, allocatable ::  &
      TIDAL_ENERGY_FLUX_2D     ! E(x,y) = input tidal energy flux at T-grid points (W/m^2)

   real (r8) ::                   &
      tidal_eps_n2,               &! lower limit N**2
      tidal_gamma_rhor,           &! tidal_mixing_efficiency/rho_fw
      zetar                        ! reciprocal of vertical decay scale for turbulence (cm)

   real (r8), public ::           &
      tidal_local_mixing_fraction  ! fraction of energy available for mixing local to the generation region

!-----------------------------------------------------------------------
!
!  Polzin formulation variables. Note that while the suffixes _P and _polz
!  are chosen to reference Polzin, the actual implementation is ala Melet 2013
!
!-----------------------------------------------------------------------


   real (r8), dimension(:,:,:), target,allocatable ::  &
      HTINV,                  &! 1/HT
      H2_P,                   &! bottom roughness used in Polzin scheme 
      U_P                      ! sqrt(barotropic tide variance) used in Polzin scheme 

   real (r8) ::               &
      coef_polz,              &! coefficient used in Polzin option
      kappa_polz,             &! wavenumber scale for topographic roughness
      mu_polz,                &! mu, a nondimensional constant
      mixing_eff_coef_polz,   &! fraction of local internal-tide energy dissipation
      nb_ref_polz,            &! reference value for N at sea bottom
      omega2_polz,            &! omega**2
      zstar_inv_coeff_polz     ! coefficient of 1/zstar_polz

!-----------------------------------------------------------------------
!
!  Schmittner method variables
!
!-----------------------------------------------------------------------

   character (char_len) ::                   &
      tidal_vert_decay_option_schm ! suboption for schmittner method: 'SSJ02' or 'P09'

   logical (log_kind) ::  &
      ltidal_all_TC_coefs_eq_1,  &! if .true., weight all TC equally for plotting (*c1)
      ltidal_all_TC_coefs_eq_p33  ! if .true., weight all TC equally for plotting (*p33)

   real (r8), dimension(:,:,:,:), allocatable ::  &
      TANH_SCHM_3D

   real (r8), dimension(:,:,:), allocatable ::  &
      TANH_TLAT_SCHM

   real (r8), dimension(:), allocatable ::  &
      decay_fn,               &! vert. decay fn, in Schmittner method
      tanhzw_schm         

   real (r8) ::               &
      hab,                    &! height above bottom -- Schmittner method
      tidal_diss_lim_TC        ! dissipation limit for tidal-energy-constituent data

!-----------------------------------------------------------------------
!
!  ER03 and NG13 tidal constituent variables
!     Energy flux out of the barotropic tide as estimated by
!     by Egbert using satellite altimetry or Green and Nycander (W/m^2)
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), target, allocatable ::  &
      TCK1,                   &! diurnal lunar
      TCM2,                   &! semidiurnal lunar 
      TCO1,                   &! diurnal solar
      TCS2                     ! semidirunal solar 

!-----------------------------------------------------------------------------
!
!   support for 18.6-year lunar cycle modulation (Lunar Nodal Cycle -- LNC)
!
!-----------------------------------------------------------------------------

   integer (int_kind), parameter :: &
      ntidal_ts_TC         = 4,                   &! number of tidal constituents (m2,s2,k1,o1)
      tidal_ts_DATA_max    = 1000.0_r8*366.00_r8, &! 1000+ years daily data, any calendar
      tidal_ts_NSTEPS_max  = 1000                  ! 1000 steps per day, max

   character (char_len_long), dimension(ntidal_ts_TC) ::   &
      tidal_energy_ts_files                        ! LNC timeseries filename

   character (char_len)  ::                       &
      tidal_energy_ts_calendar,                   &! LNC timeseries calendar
      tidal_energy_ts_file_fmt                     ! LNC timeseries file format
                                          
   integer (int_kind)        ::                   &
      tidal_energy_ts_data_first_year,            &! year of first ts data year used; this year is
                                                   !  aligned with tidal_energy_ts_model_yr_align
                                                   !  eg, ...model_year_align = 1  
                                                   !         ...ts_year_first = 1948
      tidal_energy_ts_data_first_month,           &! first month of first ts data year used
      tidal_energy_ts_data_first_day,             &! first day   of first ts data year used
      tidal_energy_ts_data_final_year,            &! year of last ts data year used 
      tidal_energy_ts_data_final_month,           &! last month of last ts data year used 
      tidal_energy_ts_data_final_day,             &! last day   of last ts data year used 
      tidal_energy_ts_model_yr_align               ! model year that is aligned with tidal_energy_ts_data_first_year

!-----------------------------------------------------------------------
!
!  public LNC variables
!
!-----------------------------------------------------------------------

   integer (int_kind), public :: tidal_ts_data_ind_low
   integer (int_kind), public :: tidal_ts_data_ind_high
   integer (int_kind), public :: tidal_ts_data_day_of_year

   type :: tidal_ts_desc
     character(char_len_long), dimension(ntidal_ts_TC)::   &
        filenames                  ! tidal modulation input timeseries filename

     character(char_len)     ::  &
        fmt,                     &! tidal modulation input timeseries file format
        data_calendar             ! data calendar ('gregorian', '365')

     logical(log_kind)     ::    &
        calendar_mismatch         ! if data calendar is gregorian and model calendar is 365

     integer(kind=int_kind)  ::  &
        model_first_year,        &! first model year  aligned with data_first_year
        model_month_first,       &! first model month aligned with data_first_month
        model_day_first,         &! first model day   aligned with data_first_day
        data_first_year,         &! first data year   aligned with model_first_year
        data_first_month,        &! first data month  aligned with model_month_first
        data_first_day,          &! first data day    aligned with model_day_first
        data_final_year,         &! final data year  used 
        data_final_month,        &! final data month used 
        data_final_day,          &! final data day   used 
        data_day_of_year,        &! data day number now in use [1,366]
        data_day_of_year_first,  &! first data day number [1,366]
        data_ind_low,            &! data index lower bound
        data_ind_high,           &! data index upper bound
        data_ind_first,          &! data index at beginning of data cycle
        data_ind_final,          &! data index at end of data cycle
        data_npoints,            &! number of timeseries data point selected
        k1,                      &! index for K1 tidal constituent
        m2,                      &! index for M2 tidal constituent
        o1,                      &! index for O1 tidal constituent
        s2                        ! index for S2 tidal constituent
 
     integer(kind=int_kind), dimension(tidal_ts_DATA_max)  ::  &
        data_year,               &! data year  now in use
        data_month,              &! data month now in use
        data_day                  ! data day   now in use


     real(kind=r8), dimension(tidal_ts_NSTEPS_max) ::  &
        factor,                  &! factor by which tidal energy components are
                                  !  multiplied in 18.6-year lunar cycle
        w1,                      &! daily weight lower
        w2                        ! daily weight upper
 
     real(kind=r8), dimension(tidal_ts_DATA_max,ntidal_ts_TC) ::    &
        data                      ! timeseries data values
  end type

  type (tidal_ts_desc) ::       &
     tidal_ts

!-----------------------------------------------------------------------
!
!  ids for tidal-mixing-related time-averaged history fields
!
!-----------------------------------------------------------------------

   integer (int_kind)         ::   &
      tavg_H2_P,                   &! id for Polzin scheme H2_P
      tavg_POLZIN_EQ2,             &! id for Jayne, Polzin, or Schmittner 2D
      tavg_TCK1,                   &! id for Schmittner diurnal lunar
      tavg_TCM2,                   &! id for Schmittner semidiurnal lunar 
      tavg_TCO1,                   &! id for Schmittner diurnal solar
      tavg_TCS2,                   &! id for Schmittner semidirunal solar 
      tavg_TIDAL_COEF_3D,          &! id for TIDAL_COEF_3D
      tavg_TIDAL_DIFF,             &! id for TIDAL_DIFF
      tavg_TIDAL_ENERGY_FLUX_2D,   &! id for Jayne, Polzin, or Schmittner 2D
      tavg_TIDAL_GAMMA_EPS,        &! id for Gamma*epsilon
      tavg_TIDAL_N2,               &! id for N**2 (no floor)
      tavg_TIDAL_N2_eps,           &! id for N**2 (epsilon floor)
      tavg_TIDAL_QE_2D,            &! id for Jayne or Polzin q*E
      tavg_TIDAL_QE_3D,            &! id for TIDAL_QE_3D (q*E)
      tavg_REGION_BOX3D,           &! id for REGION_BOX3D (dimensionless)
      tavg_U_P,                    &! id for Polzin scheme U_P 
      tavg_VERTICAL_FUNC            ! id for time-independent vertical function

   integer (int_kind)         ::   &! fields to be used in Melet comparison:
      tavg_TEMP1_1p5km,            &! id for temperature in [1km,1.5km] in 62-level models only
      tavg_TEMP1_2km,              &! id for temperature in [1km,2km] in 62-level models only
      tavg_TEMP1_3km,              &! id for temperature in [1km,3km] in 62-level models only
      tavg_TEMP2_3km,              &! id for temperature in [2km,3km] in 62-level models only
      tavg_TEMP3p5_6km,            &! id for temperature in [3.5km,6km] in 62-level models only
      tavg_TEMP2_4km,              &! id for temperature in [2km,4km] in 62-level models only
      tavg_TEMP4_6km                ! id for temperature in [4km,6km] in 62-level models only

!-----------------------------------------------------------------------------
!
!   module timer
!
!-----------------------------------------------------------------------------
   integer(kind=int_kind) ::  timer_tidal_ts  ! timer for 18.6-year LNC overhead

!-----------------------------------------------------------------------------
!
!   misc
!
!-----------------------------------------------------------------------------
   character (char_len_long) :: string               ! generic string
   character (char_len_long) :: tidal_energy_string  ! string describing tidal energy field
   character (char_len)      :: subname              ! generic subroutine name

   !... development diagnostics -- comparison with UVIC results
   logical (log_kind) :: ltidal_compare_VIC 
   real (r8) horiz_TIDAL_QE_3D(1:km)

   real (kind=r8), dimension(18), public :: sum_vic_TE, sum_dz
   real (kind=r8), dimension(18), public :: zw_vic = (/5000., 13000., 24000., 38000., 55000.,   &
            75000., 98000., 124000., 153000., 185000., 220000., 258000., 299000.,   &
           343000., 390000., 440000., 493000., 549000. /)
!          343000., 390000., 440000., 493000., 549000., 608000./)

  real (r8),public ::             &
      tidal_mixing_efficiency,    &! mixing efficiency
      vertical_decay_scale         ! vertical decay scale for turbulence (cm)

!EOP

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_tidal_mixing1
! !INTERFACE:

 subroutine init_tidal_mixing1

! !DESCRIPTION:
!  Only reads tidal-mixing namelist. Must be read prior to calling read_restart.
!  init_tidal_mixing2 initializes everything else.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     input namelist
!
!-----------------------------------------------------------------------

   namelist /tidal_nml/ tidal_local_mixing_fraction,               &
                        ltidal_all_TC_coefs_eq_1,                  &
                        ltidal_all_TC_coefs_eq_p33,                &
                        ltidal_max,                                &
                        ltidal_schmittner_socn,                    &
                        ltidal_stabc,                              &
                        ltidal_lunar_cycle,                        &
                        ltidal_melet_plot,                         &
                        ltidal_mixing,                             &
                        tidal_diss_lim_TC,                         &
                        tidal_energy_choice,                       &
                        tidal_energy_file,                         &
                        tidal_energy_file_fmt,                     &
                        tidal_energy_ts_files,                     &
                        tidal_energy_ts_file_fmt,                  &
                        tidal_energy_ts_calendar,                  &
                        tidal_energy_ts_model_yr_align,            &
                        tidal_energy_ts_data_first_year,           &
                        tidal_energy_ts_data_first_month,          &
                        tidal_energy_ts_data_first_day ,           &
                        tidal_energy_ts_data_final_year,           &
                        tidal_energy_ts_data_final_month,          &
                        tidal_energy_ts_data_final_day ,           &
                        tidal_eps_n2,                              &
                        tidal_mix_max,                             &
                        tidal_mixing_efficiency,                   &
                        tidal_mixing_method_choice,                &
                        tidal_vars_file_polz,                      &
                        tidal_vars_file_fmt_polz,                  &
                        tidal_vert_decay_option_schm,              &
                        ltidal_min_regions,                        &
                        num_tidal_min_regions,                     &
                        tidal_min_regions_name,                    &
                        tidal_min_regions_klevels,                 & 
                        tidal_min_values,                          &
                        tidal_TLATmin_regions,                     &
                        tidal_TLATmax_regions,                     &
                        tidal_TLONmin_regions,                     &
                        tidal_TLONmax_regions,                     &
                        vertical_decay_scale

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::          &
      elapsed_days,               &! temp var for day calculation
      elapsed_days_jan1,          &! temp var for day calculation
      iblock,                     &! block index
      ii,                         &
      i,j,                        &! horizontal indices
      i_size = 1,                 &! default allocation x
      j_size = 1,                 &! default allocation y
      k_size = 1,                 &! default allocation z
      k,k1,                       &! vertical level indices
      nml_error,                  &! namelist i/o error flag
      nu                           ! i/o unit number

   integer (int_kind), dimension(:,:),     allocatable  :: REGION_BOX2D_G
   integer (int_kind), dimension(:,:,:),   allocatable  :: REGION_BOX2D
   logical (log_kind), dimension(max_tidal_min_regions) :: wrap

   real (r8), dimension(:,:),   allocatable :: TLON_G, TLAT_G
   real (r8), dimension(:,:,:), allocatable ::  &
      WORK                         ! WORK array

   type (block) ::                &
      this_block                   ! block information for current block

!-----------------------------------------------------------------------
!
!     set defaults for tidal parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   tidal_local_mixing_fraction      = 0.33_r8
   ltidal_all_TC_coefs_eq_1         = .false.
   ltidal_all_TC_coefs_eq_p33       = .false.
   ltidal_max                       = .true.
   ltidal_compare_VIC               = .false.
   ltidal_schmittner_socn           = .false.
   ltidal_stabc                     = .true.
   ltidal_melet_plot                = .false.
   ltidal_mixing                    = .false.
   tidal_diss_lim_TC                = 0.0e02_r8 !cm    or 225.0e02_r8
   tidal_energy_choice              = 'unknown_tidal_energy_choice'
   tidal_energy_file                = 'unknown_tidal_energy_file'
   tidal_energy_file_fmt            = 'unknown_tidal_energy_file_fmt'
   tidal_energy_ts_files            = 'unknown_tidal_ts_file'
   tidal_energy_ts_file_fmt         = 'unknown_tidal_ts_file_fmt'
   tidal_energy_ts_calendar         = 'unknown_tidal_ts_file_calendar'
   tidal_energy_ts_model_yr_align   = 1
   tidal_energy_ts_data_first_year  = 1948
   tidal_energy_ts_data_first_month = 1
   tidal_energy_ts_data_first_day   = 1
   tidal_energy_ts_data_final_year  = 2009 
   tidal_energy_ts_data_final_month = 1 
   tidal_energy_ts_data_final_day   = 1 
   tidal_eps_n2                     = 1.0e-8_r8
   tidal_mix_max                    = 100.0_r8
   tidal_mixing_efficiency          = 0.20_r8
   tidal_mixing_method_choice       = 'unknown_tidal_mixing_method_choice'
   tidal_vars_file_polz             = 'unknown_vars_energy_file_polz'
   tidal_vars_file_fmt_polz         = 'nc'
   tidal_vert_decay_option_schm     = 'SSJ02'
   vertical_decay_scale             = 500.0e02_r8

   tidal_min_regions_name           = ' '
   ltidal_min_regions               = .false. 
   num_tidal_min_regions            = 3
   tidal_min_values                 = 20.0_r8
   tidal_min_regions_klevels        = 6
 
!-----------------------------------------------------------------------
!
!  read namelist input and broadcast variables
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tidal_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar (nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP (SigAbort, 'ERROR reading tidal_nml')
   endif

!-----------------------------------------------------------------------
!
!  document namelist in output log file
!
!-----------------------------------------------------------------------
   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,ndelim_fmt)
     write(stdout,blank_fmt)
     write(stdout,*) ' Tidal mixing information'
     write(stdout,blank_fmt)

     !-----------------------------------------------------------------------
     !
     !  tidal mixing method
     !
     !-----------------------------------------------------------------------

     select case (tidal_mixing_method_choice(1:3))
     case ('jay')
       tidal_mixing_method_itype = tidal_mixing_method_jayne
       string = 'Jayne'
     case ('pol')
       tidal_mixing_method_itype = tidal_mixing_method_polzin
       string = 'Polzin/Melet'
     case ('sch','Sch')
       tidal_mixing_method_itype = tidal_mixing_method_schmittner
       string = 'Schmittner'
     case default
       tidal_mixing_method_itype = -1000
       string = 'Unknown'
     end select
     write(stdout,*) trim(string)// ' Tidal Mixing'

     !-----------------------------------------------------------------------
     !
     !  tidal energy file
     !
     !-----------------------------------------------------------------------

     select case (tidal_energy_choice(1:4))
     case ('arbi')
       tidal_energy_string = 'Arbic'
       tidal_energy_itype = tidal_energy_arbic
     case ('jayn')
       tidal_energy_string = 'Jayne' 
       tidal_energy_itype = tidal_energy_jayne
     case ('ER03','er03')
       tidal_energy_string = 'ER03' 
       tidal_energy_itype = tidal_energy_ER03
     case ('GN13','gn13')
       tidal_energy_string = 'GN13' 
       tidal_energy_itype = tidal_energy_GN13
     case ('LGM0','lgm0')
       tidal_energy_string = 'LGM0 -- research only do not use without author permission' 
       tidal_energy_itype = tidal_energy_LGM0
     case ('LGMi','lgmi')
       select case (tidal_energy_choice(1:8))
        case('LGMi5g21','lgmi5g21')
          tidal_energy_string = 'LGMi5g21 -- research only do not use without author permission' 
          tidal_energy_itype = tidal_energy_LGMi5g21
        case('LGMi6g21','lgmi6g21')
          tidal_energy_string = 'LGMi6g21 -- research only do not use without author permission' 
          tidal_energy_itype = tidal_energy_LGMi6g21
       end select
     case default
       tidal_energy_itype = -1000
       tidal_energy_string = 'Unknown'
     end select


     write(stdout,*) 'with '//trim(tidal_energy_string)// ' Tidal Energy file'

     write(stdout,blank_fmt)
     write(stdout,*) ' tidal_nml namelist:'
     write(stdout,blank_fmt)
     write(stdout,tidal_nml)
     write(stdout,blank_fmt)
   endif

   call broadcast_scalar (tidal_local_mixing_fraction,      master_task)
   call broadcast_scalar (ltidal_all_TC_coefs_eq_1,         master_task)
   call broadcast_scalar (ltidal_all_TC_coefs_eq_p33,       master_task)
   call broadcast_scalar (ltidal_max,                       master_task)
   call broadcast_scalar (ltidal_schmittner_socn,           master_task)
   call broadcast_scalar (ltidal_stabc,                     master_task)
   call broadcast_scalar (ltidal_lunar_cycle,               master_task)
   call broadcast_scalar (ltidal_mixing,                    master_task)
   call broadcast_scalar (tidal_diss_lim_TC,                master_task)
   call broadcast_scalar (tidal_energy_choice,              master_task)
   call broadcast_scalar (tidal_energy_itype,               master_task)
   call broadcast_scalar (tidal_energy_file,                master_task) 
   call broadcast_scalar (tidal_energy_file_fmt,            master_task) 
   call broadcast_scalar (tidal_energy_ts_file_fmt,         master_task) 
   call broadcast_scalar (tidal_energy_ts_calendar,         master_task) 
   call broadcast_scalar (tidal_energy_ts_model_yr_align,   master_task) 
   call broadcast_scalar (tidal_energy_ts_data_first_year,  master_task) 
   call broadcast_scalar (tidal_energy_ts_data_first_month, master_task) 
   call broadcast_scalar (tidal_energy_ts_data_first_day,   master_task) 
   call broadcast_scalar (tidal_energy_ts_data_final_year,  master_task) 
   call broadcast_scalar (tidal_energy_ts_data_final_month, master_task) 
   call broadcast_scalar (tidal_energy_ts_data_final_day,   master_task) 
   call broadcast_scalar (tidal_eps_n2,                     master_task) 
   call broadcast_scalar (tidal_mix_max,                    master_task) 
   call broadcast_scalar (tidal_mixing_efficiency,          master_task)
   call broadcast_scalar (tidal_mixing_method_choice,       master_task)
   call broadcast_scalar (tidal_mixing_method_itype,        master_task)
   call broadcast_scalar (tidal_vars_file_polz,             master_task) 
   call broadcast_scalar (tidal_vars_file_fmt_polz,         master_task) 
   call broadcast_scalar (tidal_vert_decay_option_schm,     master_task) 
   call broadcast_scalar (vertical_decay_scale,             master_task)
   call broadcast_scalar (tidal_energy_string,              master_task)
   call broadcast_scalar (num_tidal_min_regions,            master_task)
   call broadcast_scalar (ltidal_min_regions,               master_task)         

   !... loop over vectors one element at a time
   do ii=1,num_tidal_min_regions
   call broadcast_scalar (tidal_min_regions_name(ii),       master_task)         
   call broadcast_scalar (tidal_min_values(ii),             master_task)        
   call broadcast_scalar (tidal_TLATmin_regions(ii),        master_task)         
   call broadcast_scalar (tidal_TLATmax_regions(ii),        master_task)      
   call broadcast_scalar (tidal_TLONmin_regions(ii),        master_task)       
   call broadcast_scalar (tidal_TLONmax_regions(ii),        master_task)        
   enddo
   
   !... loop over vectors one element at a time
   do ii=1,ntidal_ts_TC
   call broadcast_scalar (tidal_energy_ts_files(ii),        master_task) 
   enddo

   call register_string('tidal_energy_'//trim(tidal_energy_string))


!-----------------------------------------------------------------------
!
!  if applying a tidal-mixing floor in specified regions, define the
!   mask now
!
!-----------------------------------------------------------------------

 if (ltidal_min_regions) then

 !... document regions in output log file
   if (my_task == master_task) then
   do ii=1,num_tidal_min_regions
     write(stdout,1100)  '(init_tidal_mixing1) num_region = ', ii
     write(stdout,1101)  '(init_tidal_mixing1) tidal_min_regions_name = ', trim(tidal_min_regions_name(ii))
     write(stdout,1102)  '(init_tidal_mixing1) tidal_TLATmin_regions  = ', tidal_TLATmin_regions(ii)
     write(stdout,1102)  '(init_tidal_mixing1) tidal_TLATmax_regions  = ', tidal_TLATmax_regions(ii)
     write(stdout,1102)  '(init_tidal_mixing1) tidal_TLONmin_regions  = ', tidal_TLONmin_regions(ii)
     write(stdout,1102)  '(init_tidal_mixing1) tidal_TLONmax_regions  = ', tidal_TLONmax_regions(ii)
     write(stdout,*) ' '
     call shr_sys_flush(stdout)
   enddo !num_region
   endif !my_task

 !... determine if longitude is wrapped

   do ii=1,num_tidal_min_regions
    if ( tidal_TLONmin_regions(ii) .le. tidal_TLONmax_regions(ii)) then
      !... no wrap
      wrap(ii) = .false.
    else 
      !... assume box wraps around 360, so add 360 to max
      wrap(ii) = .true.
      call document ('init_tidal_mixing1','box wraps around 360 for tidal region box ', ii)
    endif
   enddo ! ii

 !... create 2D global tidal-region box mask
   
   allocate (REGION_BOX2D_G(nx_global,ny_global)); REGION_BOX2D_G = 0

   !... create 2D global LAT,LON in degrees
   allocate (TLON_G(nx_global,ny_global),TLAT_G(nx_global,ny_global))
   call gather_global(TLON_G, TLON, master_task,distrb_clinic)
   call gather_global(TLAT_G, TLAT, master_task,distrb_clinic)
   TLON_G = TLON_G*radian
   TLAT_G = TLAT_G*radian

   if (my_task == master_task) then
     do ii=1,num_tidal_min_regions
           write(stdout,1104)'(init_tidal_mixing1) tidal-mixing region latitudes/longitudes = ', trim(tidal_min_regions_name(ii))
       do j=1,ny_global
       do i=1,nx_global

        if(TLAT_G(i,j) .ge. tidal_TLATmin_regions(ii)  .and.  &
           TLAT_G(i,j) .le. tidal_TLATmax_regions(ii) ) then
           if (.not. wrap(ii)) then
             if (TLON_G(i,j) .ge. tidal_TLONmin_regions(ii)  .and.  &
                 TLON_G(i,j) .le. tidal_TLONmax_regions(ii))  then
                 REGION_BOX2D_G(i,j) = ii
                 write(stdout,1105)  &
                   '(init_tidal_mixing1) REGION_BOX2D_G(i,j), TLATmin, TLAT_G(i,j), TLATmax, TLONmin, TLON_G(i,j), TLONmax ', &
                                         REGION_BOX2D_G(i,j), &
                                         tidal_TLATmin_regions(ii), TLAT_G(i,j),tidal_TLATmax_regions(ii), &
                                         tidal_TLONmin_regions(ii), TLON_G(i,j),tidal_TLONmax_regions(ii)
             endif
           else ! box longitudes wrap around 360
             if ((TLON_G(i,j) .ge. tidal_TLONmin_regions(ii)) .and. &
                 (TLON_G(i,j) .le. 360.0_r8                 ) .or.  & 
                 (TLON_G(i,j) .le. tidal_TLONmax_regions(ii))) then
                 REGION_BOX2D_G(i,j) = ii
                 write(stdout,1105)  &
                   '(init_tidal_mixing1) REGION_BOX2D_G(i,j), TLATmin, TLAT_G(i,j), TLATmax, TLONmin, TLON_G(i,j), TLONmax ', &
                                         REGION_BOX2D_G(i,j), &
                                         tidal_TLATmin_regions(ii), TLAT_G(i,j),tidal_TLATmax_regions(ii), &
                                         tidal_TLONmin_regions(ii), TLON_G(i,j),tidal_TLONmax_regions(ii)
             endif
           endif ! wrap
        endif ! lat 
       enddo ! i
       enddo ! j
     enddo ! ii
   endif ! master_task

 !... create 2D local tidal-region box mask
   
   allocate (REGION_BOX2D (nx_block,ny_block,   nblocks_clinic)); REGION_BOX2D = 0
   allocate (REGION_BOX3D (nx_block,ny_block,km,nblocks_clinic)); REGION_BOX3D = 0
   call scatter_global(REGION_BOX2D, REGION_BOX2D_G, master_task, distrb_clinic, &
                       field_loc_center, field_type_scalar)

 !... sanity test; note that non-zero REGION_BOX2D are not on master_task -- see the cesm log file 
   if (my_task == master_task) then
   if (any(REGION_BOX2D_G .ne. 0)) then
     write(stdout,*) '(init_tidal_mixing1) not all REGION_BOX2D_G values are zero'
   else
     write(stdout,*) '(init_tidal_mixing1)     ALL REGION_BOX2D_G values are zero'
     call exit_POP (SigAbort, 'ERROR  ALL REGION_BOX2D_G values are zero')
   endif
   endif ! master_task

   if (any(REGION_BOX2D .ne. 0)) then
     write(stdout,*) my_task, ': (init_tidal_mixing1) not all REGION_BOX2D values are zero'
   else
     write(stdout,*) my_task, ': (init_tidal_mixing1)     ALL REGION_BOX2D values are zero'
   endif

 !... create 3D local tidal-region box mask
   do iblock = 1,nblocks_clinic
      do ii=1,num_tidal_min_regions
      do k=1,km
       do j=1,ny_block
       do i=1,nx_block
         if (REGION_BOX2D(i,j,iblock) == ii ) then
          if ( tidal_min_regions_klevels(ii) == 2 .and. k > 2) then
            if (k == KMT(i,j,iblock)-1 .or. k == KMT(i,j,iblock)-2) then
                REGION_BOX3D(i,j,k,iblock) = ii
            endif
          else if (tidal_min_regions_klevels(ii) == 6 .and. k > 6) then
            if (k == KMT(i,j,iblock)-1 .or. k == KMT(i,j,iblock)-2  .or.  &
                k == KMT(i,j,iblock)-3 .or. k == KMT(i,j,iblock)-4  .or.  &
                k == KMT(i,j,iblock)-5 .or. k == KMT(i,j,iblock)-6 ) then
                REGION_BOX3D(i,j,k,iblock) = ii
            endif
          endif 
         endif ! REGION_BOX2D > 0
       enddo ! i
       enddo ! j
      enddo ! k
     enddo ! ii
   enddo ! iblock

   call define_tavg_field(tavg_REGION_BOX3D,'REGION_BOX3D',3,     &
                          long_name='Regions where minimum tidal diffusion is applied', &
                          tavg_method=tavg_method_constant,       &
                          units=' ',                              &
                          grid_loc='3111',                        &
                          coordinates  ='TLONG TLAT z_w time'     ) 

   do iblock = 1,nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)
     do ii=1,num_tidal_min_regions
     !write(stdout,*) trim(tidal_min_regions_name(ii))
      do k=1,km
       do j=1,ny_block
       do i=1,nx_block
         if (REGION_BOX3D(i,j,k,iblock) > 0) then
          write(stdout,1212) '(init_tidal_mixing1) reg num, REGION_BOX3D, i,j,k,iblock,lat,lon',  &
            trim(tidal_min_regions_name(ii)), ii, REGION_BOX3D(i,j,k,iblock), &
            this_block%i_glob(i), this_block%j_glob(j),k,iblock,  &
            TLATD(i,j,iblock), TLOND(i,j,iblock)
         endif ! REGION_BOX3D > 0
       enddo ! i
       enddo ! j
      enddo ! k
     enddo ! ii
   enddo ! iblock
1212 format(1x, a, a, 6i4, 1x, 2F10.2)
1213 format(1x, a, 6i3)
write(stdout,*) ' after REGION_BOX3D print test'

   deallocate (REGION_BOX2D_G, REGION_BOX2D)
   deallocate (TLAT_G, TLON_G)

1100 format(1x, a,i2)
1101 format(1x, a,1x,a)
1102 format(1x, a,1x,F9.2)
1104 format(1x,a, a)
1105 format(1x, a, i3, 3x, 3F9.2, 3x, 3F9.2)
 endif ! ltidal_min_regions

!-----------------------------------------------------------------------
!EOC

 end subroutine init_tidal_mixing1

!BOP
! !IROUTINE: init_tidal_mixing2
! !INTERFACE:

 subroutine init_tidal_mixing2

! !DESCRIPTION:
!  Reads tidal-mixing namelist. Must be read prior to calling read_restart.
!  Initializes parameters for tidally driven mixing, reads input file,
!  and computes the time-independent part of the mixing coefficients

! !REVISION HISTORY:
!  same as module

!EOP
!BOC


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::          &
      elapsed_days,               &! temp var for day calculation
      elapsed_days_jan1,          &! temp var for day calculation
      iblock,                     &! block index
      ii,                         &
      i,j,                        &! horizontal indices
      i_size = 1,                 &! default allocation x
      j_size = 1,                 &! default allocation y
      k_size = 1,                 &! default allocation z
      k,k1,                       &! vertical level indices
      nml_error,                  &! namelist i/o error flag
      nu                           ! i/o unit number

   real (r8), dimension(:,:,:), allocatable ::  &
      WORK                         ! WORK array

   type (block) ::                &
      this_block                   ! block information for current block

!-----------------------------------------------------------------------
!
!     set defaults for tidal parameters, then read them from namelist
!
!-----------------------------------------------------------------------

   if (ltidal_mixing) then
     i_size = nx_block
     j_size = ny_block
     k_size = km
   endif

!-----------------------------------------------------------------------
!
!  allocate and initialize some arrays and variables
!    minimal allocation if tidal_mixing is not enabled
!
!-----------------------------------------------------------------------

   allocate (TIDAL_DIFF   (i_size,j_size,k_size,nblocks_clinic)); TIDAL_DIFF    = c0
   allocate (TIDAL_N2     (i_size,j_size,k_size,nblocks_clinic)); TIDAL_N2      = c0
   allocate (TIDAL_N2_eps (i_size,j_size,k_size,nblocks_clinic)); TIDAL_N2_eps  = c0
   allocate (VERTICAL_FUNC(i_size,j_size,k_size,nblocks_clinic)); VERTICAL_FUNC = c0
   allocate (VERTICAL_NUM (i_size,j_size,k_size,nblocks_clinic)); VERTICAL_NUM  = c0

   allocate (TIDAL_ENERGY_FLUX_2D(i_size,j_size,nblocks_clinic)); TIDAL_ENERGY_FLUX_2D = c0
   allocate (TIDAL_QE_2D         (i_size,j_size,nblocks_clinic)); TIDAL_QE_2D = c0
 
!-----------------------------------------------------------------------
!
!  return after minimal allocation if tidal_mixing is not enabled
!
!-----------------------------------------------------------------------

   if (.not. ltidal_mixing) return

!-----------------------------------------------------------------------
!
!  internal consistency checks and exit conditions
!
!-----------------------------------------------------------------------

   call tidal_check

!-----------------------------------------------------------------------
!
!     Set up CVMix
!
!-----------------------------------------------------------------------

   select case (tidal_mixing_method_itype)

      case (tidal_mixing_method_jayne)
      !==============================

        call cvmix_init_tidal(mix_scheme='Simmons',                                &
                              efficiency=tidal_mixing_efficiency,                  &
                              vertical_decay_scale=vertical_decay_scale*1.0e-2_r8, &
                              max_coefficient=tidal_mix_max*1.0e-4_r8,             &
                              local_mixing_frac=tidal_local_mixing_fraction,       &
                              depth_cutoff=0.0_r8,                                 &
                              ltidal_Schmittner_socn=ltidal_schmittner_socn) 

      case (tidal_mixing_method_schmittner)
           !==============================
        call cvmix_init_tidal(mix_scheme='Schmittner',                             &
                              efficiency=tidal_mixing_efficiency,                  &
                              vertical_decay_scale=vertical_decay_scale*1.0e-2_r8, &
                              max_coefficient=tidal_mix_max*1.0e-4_r8,             &
                              local_mixing_frac=tidal_local_mixing_fraction,       &
                              depth_cutoff=0.0_r8,                                 &
                              ltidal_Schmittner_socn=ltidal_schmittner_socn) 
   end select ! tidal_mixing_method_itype


!-----------------------------------------------------------------------
!
!     misc definitions
!
!-----------------------------------------------------------------------

   tidal_gamma_rhor = tidal_mixing_efficiency/rho_fw
   zetar = c1/vertical_decay_scale

!-----------------------------------------------------------------------
!
!  read the selected tidal-energy file
!
!-----------------------------------------------------------------------

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', km)

   select case (tidal_energy_itype)

      !---------------------------------------------
      ! 2D tidal energy datasets:
      !    read E(x,y)
      !    form TIDAL_QE_2D = q*E(x,y)
      !
      ! 3D tidal energy datasets:
      !    read TC(x,y,z,ntc)
      !    form TIDAL_QE_3D = sum_tc(q_tc*TC(x,y,z)) = q*E(x,y,z)
      !    form TIDAL_QE_2D = sum_k(E(x,y,z))
      !      always form tidal_qE_2D, although it is not always used
      !---------------------------------------------

      case (tidal_energy_jayne)
           !==================

         call tidal_read_energy_jayne 
         call tidal_form_qE_2D (tidal_local_mixing_fraction)

      case (tidal_energy_arbic)
           !==================

         call tidal_read_energy_arbic  
         call tidal_form_qE_2D (tidal_local_mixing_fraction)

      case (tidal_energy_ER03)
           !=================

         call tidal_read_energy_TC ('ER03') 
         call tidal_form_qE_3D  (c1,c1,c1,c1)
         call tidal_form_qE_2D_from_3D

      case (tidal_energy_GN13)
           !=================

         call tidal_read_energy_TC ('GN13')
         call tidal_form_qE_3D  (c1,c1,c1,c1)
         call tidal_form_qE_2D_from_3D

      case (tidal_energy_LGM0)
           !=================

         !active research option; do not use
         call tidal_read_energy_TC ('LGM0')
         call tidal_form_qE_3D  (c1,c1,c1,c1)
         call tidal_form_qE_2D_from_3D

      case (tidal_energy_LGMi5g21)
           !=====================

         !active research option; do not use
         call tidal_read_energy_TC ('LGMi5g21')
         call tidal_form_qE_3D  (c1,c1,c1,c1)
         call tidal_form_qE_2D_from_3D

      case (tidal_energy_LGMi6g21)
           !=====================

         !active research option; do not use
         call tidal_read_energy_TC ('LGMi6g21')
         call tidal_form_qE_3D  (c1,c1,c1,c1)
         call tidal_form_qE_2D_from_3D

   end select ! tidal_energy_itype

!-----------------------------------------------------------------------
!
!  option-specific mixing-method initializations 
!
!-----------------------------------------------------------------------

   select case (tidal_mixing_method_itype)

      case (tidal_mixing_method_jayne)
      !==============================

        !---------------------------------------------------------------
        !
        !  compute time-independent parts of Jayne tidal mixing method 
        !
        !---------------------------------------------------------------

        allocate (TIDAL_COEF_3D(nx_block,ny_block,km,nblocks_clinic)) ; TIDAL_COEF_3D = c0
        allocate (TIDAL_COEF_2D(nx_block,ny_block,   nblocks_clinic)) ; TIDAL_COEF_2D = c0
        allocate (WORK         (nx_block,ny_block,   nblocks_clinic)) ; WORK          = c0

        do iblock = 1,nblocks_clinic
          ! see Gokhan notes 2014-09-23
          do k=1,km
            where ( k < KMT(:,:,iblock) ) 
              WORK(:,:,iblock) = WORK(:,:,iblock) + exp(-(HT(:,:,iblock)-zw(k))/vertical_decay_scale)*dzw(k)
            endwhere
          end do ! k

          do k=1,km
            !where ( k <= KMT(:,:,iblock) )
            !  VERTICAL_NUM  (:,:,k,iblock) = exp(-(HT(:,:,iblock) - zw(k))/vertical_decay_scale)
            !  VERTICAL_FUNC(:,:,k,iblock) = VERTICAL_NUM(:,:,k,iblock)/WORK(:,:,iblock)
!           !  !*** the following line was moved to tidal_form_coef_jayne
!           !  TIDAL_COEF_3D(:,:,k,iblock) = tidal_gamma_rhor*TIDAL_QE_2D(:,:,iblock)*VERTICAL_FUNC(:,:,k,iblock)
            !endwhere
            ! The above where-block is divided as follows for generalizing it for partial bottom cells.
            where ( k < KMT(:,:,iblock) )
              VERTICAL_NUM  (:,:,k,iblock) = exp(-(HT(:,:,iblock) - zw(k))/vertical_decay_scale)
              VERTICAL_FUNC(:,:,k,iblock) = VERTICAL_NUM(:,:,k,iblock)/WORK(:,:,iblock)
            endwhere
            where ( k == KMT(:,:,iblock) )
              VERTICAL_NUM  (:,:,k,iblock) = 1.0_r8
              VERTICAL_FUNC(:,:,k,iblock) = VERTICAL_NUM(:,:,k,iblock)/WORK(:,:,iblock)
            endwhere
          enddo ! k
        enddo ! iblock

        !*** form TIDAL_COEF_3D according to Jayne method
        call tidal_form_coef_jayne

        if (ltidal_compare_VIC) & 
        call tidal_compare_VIC(WORK, 'jayne')

        deallocate (WORK)

      case (tidal_mixing_method_polzin)
           !========================== 

        !--------------------------------------------------------------------
        !
        !  compute time-independent parts of Polzin/Melet tidal mixing method 
        !
        !--------------------------------------------------------------------

        allocate (HTINV         (nx_block,ny_block,nblocks_clinic)) ; HTINV = c0
        allocate (H2_P          (nx_block,ny_block,nblocks_clinic)) ; H2_P  = c0
        allocate (U_P           (nx_block,ny_block,nblocks_clinic)) ; U_P   = c0
        allocate (TIDAL_COEF_2D (nx_block,ny_block,nblocks_clinic)) ; TIDAL_COEF_2D = c0

        !***  Polzin/Melet contstants
        mu_polz             = 6.97e-2_r8
        nb_ref_polz         = 9.6e-4_r8
        kappa_polz          = pi2/125._r8             ! kappa:lengthscale = "kappa = 2*pi/125 km"
                                                      !   corresponds to Jayne vars file value
        kappa_polz          = kappa_polz*1.0e-5_r8    ! convert to 1/cm
        omega2_polz         = omega**2                ! omega=earth angular velocity, from init_constants
        zstar_inv_coeff_polz= kappa_polz**2/(mu_polz*nb_ref_polz**2)

        !***  read 3D vars used in Polzin/Melet formulation (roughness & urms)
        call tidal_read_roughness_RMS 

        where (HT /= c0) 
          HTINV = c1/HT
        elsewhere
          HTINV = 0.001_r8 ! assume minimum H is 10m*100cm/m
        endwhere
         
        !*** (gamma/rho)*q*E(x,y)

        TIDAL_COEF_2D = tidal_gamma_rhor*RCALCT(:,:,1:nblocks_clinic)*TIDAL_QE_2D

        if (any(TIDAL_COEF_2D < c0)) call shr_sys_abort ('Polzin/Melet method: negative TIDAL_COEF_2D terms')

      case (tidal_mixing_method_schmittner)
           !==============================

        !---------------------------------------------------------------
        !
        !  compute time-independent parts of Schmittner tidal mixing 
        !   subgridscale parameterization
        !
        !---------------------------------------------------------------

        allocate (TIDAL_COEF_3D(nx_block,ny_block,km,nblocks_clinic)); TIDAL_COEF_3D = c0
        allocate (decay_fn      (km)); decay_fn      = c0

        if (trim(tidal_vert_decay_option_schm) == 'SSJ02') then

          !*** in schmittner method, use vertical decay function as in St Laurent et al 2002
          do k=1,km
          decay_fn(k) = zetar/(c1-exp(-zetar*zw(k)))
          enddo !k

          call document ('init_tidal_mixing2','using SSJ02 tidal_vert_decay_option_schm ', &
              trim(tidal_vert_decay_option_schm) )

        else if (trim(tidal_vert_decay_option_schm) == 'P09') then
          !*** in polzin method, use vertical decay function as in Melet 2013
          do k=1,km
       !   decay_fn(k) = TBD --  UNCOMMENT WHEN decay_fn is determined
          enddo !k

          call document ('init_tidal_mixing2','using P09 tidal_vert_decay_option_schm ', &
              trim(tidal_vert_decay_option_schm) )
          call shr_sys_abort ('invalid tidal_vert_decay_option_schm')

        else
          call document ('init_tidal_mixing2','invalid tidal_vert_decay_option_schm ', &
              trim(tidal_vert_decay_option_schm) )
          call shr_sys_abort ('invalid tidal_vert_decay_option_schm')
        endif

        !*** form TIDAL_COEF_3D using Schmittner method
        call tidal_form_coef_schm

        if (ltidal_compare_VIC) then 
          if (.not.allocated(WORK)) &
            allocate (WORK (nx_block,ny_block,nblocks_clinic)) ; WORK = c0

          call tidal_compare_VIC(WORK, 'schmittner')
          deallocate (WORK)
        endif

      case default

   end select ! tidal_mixing_method_itype

!--------------------------------------------------------------------
!  time-independent terms in schmittner high mixing in Southern Ocean 
!  below 500m as in observations.  Original code:
!  zkappa = max(zkappa,tanh((zw(k)-50000.)/10000)*(1-tanh((tlat(i,j+joff)+40)/8))/2)
!--------------------------------------------------------------------

   if (ltidal_schmittner_socn)  then
     allocate (tanhzw_schm   (km)) ; tanhzw_schm   = c0
     allocate (TANH_TLAT_SCHM(nx_block,ny_block,   nblocks_clinic)) ; TANH_TLAT_SCHM = c0
     allocate (TANH_SCHM_3D  (nx_block,ny_block,km,nblocks_clinic)) ; TANH_SCHM_3D = c0

     TANH_TLAT_SCHM(:,:,:) = p5*(c1-tanh((TLATD(:,:,:)+40._r8)/8._r8))

     do k=1,km
      tanhzw_schm (k) = tanh((zw(k)-50000._r8)/10000._r8)
      TANH_SCHM_3D(:,:,k,:) = tanhzw_schm (k)*TANH_TLAT_SCHM(:,:,:)
     enddo !k

     deallocate (TANH_TLAT_SCHM, tanhzw_schm)

   endif !ltidal_schmittner_socn


!---------------------------------------------------------------
!
!   set up tidal-mixing tavg ids
!
!---------------------------------------------------------------

   call init_tidal_tavg


!---------------------------------------------------------------
!
!   DO NOT call init_tidal_ts here, because it needs info from
!     restart file which is not yet available -- call from 
!     initial.F90 instead.
!
!---------------------------------------------------------------


!-----------------------------------------------------------------------
!EOC

 end subroutine init_tidal_mixing2

!***********************************************************************
!BOP
! !IROUTINE: init_tidal_ts
! !INTERFACE:

     subroutine init_tidal_ts

! !DESCRIPTION: This subroutine reads the tidal-energy 18.6-year
!               lunar nodal cycle (LNC) modulation timeseries file and initializes
!               related variables
   
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character(80)   :: readline
   character (3)   ::  &
     cmonths(12) = (/'Jan','Feb','Mar','Apr','May','Jun',  &
                     'Jul','Aug','Sep','Oct','Nov','Dec'/) 

   character (3)      ::  &
     cmonth

   integer (int_kind) ::  &
     ii, ind, index,      &! local indices
     ioerr,               &! io error index
     nn,nnn,              &! local indices
     nu,                  &! temporary unit number for reading data
     yyyy,mm,dd,          &! temporary year,month,day data holders
     m2,s2,k1,o1           ! temporary tidal constituent indices

   integer (int_kind) :: tidal_elapsed_days
   integer (int_kind) :: tidal_elapsed_days_jan1

   real (r8), dimension(ntidal_ts_TC)    :: ts_data
   real (r8), dimension(:), allocatable  :: dayfrac

!-----------------------------------------------------------------------
!  if time-varying lunar cycle is not active, return
!-----------------------------------------------------------------------
   if (.not. ltidal_lunar_cycle) return

   allocate (dayfrac(0:nsteps_per_day))

   !*** initialize local indices. The indices must be sequential, match the
   !    input data indices, and start at the value one.
   m2 = 1
   s2 = 2
   k1 = 3
   o1 = 4

!-------------------------------------------------------------------------
!  initialize tidal timeseries variables with defaults and namelist values
!-------------------------------------------------------------------------
   do ii=1,ntidal_ts_TC
   tidal_ts%filenames(ii)    = trim(tidal_energy_ts_files(ii))
   enddo
   tidal_ts%fmt              = trim(tidal_energy_ts_file_fmt)
   tidal_ts%data_calendar    = trim(tidal_energy_ts_calendar)
   tidal_ts%model_first_year = tidal_energy_ts_model_yr_align
   tidal_ts%data_first_year  = tidal_energy_ts_data_first_year
   tidal_ts%data_first_month = tidal_energy_ts_data_first_month
   tidal_ts%data_first_day   = tidal_energy_ts_data_first_day
   tidal_ts%data_final_year  = tidal_energy_ts_data_final_year
   tidal_ts%data_final_month = tidal_energy_ts_data_final_month
   tidal_ts%data_final_day   = tidal_energy_ts_data_final_day
   tidal_ts%data_year        = tidal_ts%model_first_year
   tidal_ts%m2               = m2
   tidal_ts%s2               = s2
   tidal_ts%k1               = k1
   tidal_ts%o1               = o1
   tidal_ts%data_year        = 0
   tidal_ts%data_month       = 0
   tidal_ts%data_day         = 0
   tidal_ts%data             = c0

!-----------------------------------------------------------------------
!  read timeseries data from file
!-----------------------------------------------------------------------
   file_loop: do ii=1,ntidal_ts_TC
     !*** initialize counter to zero
     tidal_ts%data_npoints = 0

     !*** get unit and open files for each tidal constituent
     call get_unit(nu)

     !*** should add a test for this: does ii=tidal_ts%m2, etc?
     if (my_task == master_task) then
        open(nu,file=trim(tidal_ts%filenames(ii)), form='formatted',iostat=ioerr)
     endif

     write(stdout,*) '(init_tidal_ts) open ', trim(tidal_ts%filenames(ii))

     call broadcast_scalar(ioerr, master_task)
     if (ioerr /= 0) call exit_POP(sigAbort, 'Error opening tidal timeseries file')

     !*** should test if each file the same length


     !*** initialize to negative number for error checking
     tidal_ts%data_ind_final = -1

     if (my_task == master_task) then
        tidal_ts_loop: do nn=1,tidal_ts_DATA_max

        !*** read timeseries data
        read(nu,'(a)',iostat=ioerr,end=9999) readline
        read (readline,1100) dd,cmonth,yyyy,ts_data(ii)
        1100 format (i2,1x,a3,1x,i4,1x,8x,F25.6)
        if (ioerr .ne. 0) exit tidal_ts_loop

        month_loop: do mm=1,12
          if (trim(cmonth) == cmonths(mm)) then
             exit month_loop
          endif
          if (mm .eq. 12) then
              call exit_POP(sigAbort, 'Error reading tidal timeseries file')
          endif
        enddo month_loop

        !*** save index of first data date
        if (yyyy .eq. tidal_energy_ts_data_first_year  .and.  &
            mm   .eq. tidal_energy_ts_data_first_month .and.  &
            dd   .eq. tidal_energy_ts_data_first_day          ) then
            tidal_ts%data_ind_first = nn
        endif

        !*** store all data from the beginning of the record until the 
        !    final data date selected at run time
        tidal_ts%data_year(nn)  = yyyy
        tidal_ts%data_month(nn) = mm
        tidal_ts%data_day(nn)   = dd
        tidal_ts%data(nn,ii)    = ts_data(ii)

        !*** count the number of active timeseries data points
        if (nn .ge. tidal_ts%data_ind_first) then
          tidal_ts%data_npoints = tidal_ts%data_npoints + 1
        endif

        !*** at final data date, save index and exit loop
        if (yyyy .eq. tidal_energy_ts_data_final_year  .and.  &
            mm   .eq. tidal_energy_ts_data_final_month .and.  &
            dd   .eq. tidal_energy_ts_data_final_day          ) then
            tidal_ts%data_ind_final = nn
            exit tidal_ts_loop
        endif

        enddo tidal_ts_loop

        !*** close unit
   endif ! master_task
   close(nu)
   call release_unit(nu)
9999 continue

   call broadcast_scalar(ioerr, master_task)
   if (ioerr /= 0) call exit_POP(sigAbort, &
                                 'Error reading tidal_ts%filename')
   enddo file_loop

!-----------------------------------------------------------------------
!  share information with all processors
!-----------------------------------------------------------------------
   call broadcast_scalar(tidal_ts%data_ind_first,  master_task)
   call broadcast_scalar(tidal_ts%data_ind_final,  master_task)
   call broadcast_scalar(tidal_ts%data_npoints,    master_task)
   call broadcast_array (tidal_ts%data_year,       master_task)
   call broadcast_array (tidal_ts%data_month,      master_task)
   call broadcast_array (tidal_ts%data_day,        master_task)
   call broadcast_array (tidal_ts%data,            master_task)

   !*** error checking -- after broadcast
   if (tidal_ts%data_ind_final < 0) then
     call document ('init_tidal_ts', 'tidal_ts_DATA_max',tidal_ts_DATA_max)
     call document ('init_tidal_ts', 'tidal_ts%data_npoints',tidal_ts%data_npoints)
     call document ('init_tidal_ts', ' you probably need to increase tidal_ts_DATA_max')
     call exit_POP (SigAbort, 'ERROR reading tidal timeseries dataset')
   endif

!-----------------------------------------------------------------------
!  initialize tidal timeseries variable using ts data
!-----------------------------------------------------------------------

   call ymd2eday (tidal_ts%data_first_year,   &
                  tidal_ts%data_first_month,  &
                  tidal_ts%data_first_day,    &
                  tidal_elapsed_days)

   call ymd2eday (tidal_ts%data_first_year, 1,1,tidal_elapsed_days_jan1)
   tidal_ts%data_day_of_year_first = tidal_elapsed_days - tidal_elapsed_days_jan1 + 1

!-----------------------------------------------------------------------
!  necessary information for exact restarts
!-----------------------------------------------------------------------
   if (registry_match('ccsm_continue') .or. registry_match('restart')) then
       !*** from restart file
       tidal_ts%data_ind_low     = tidal_ts_data_ind_low  
       tidal_ts%data_ind_high    = tidal_ts_data_ind_high
       tidal_ts%data_day_of_year = tidal_ts_data_day_of_year
    else
       !*** initialization
       tidal_ts%data_ind_low     = tidal_ts%data_ind_first
                                                             !assume for now it is not Feb 28th or 29th.
       tidal_ts%data_ind_high    = tidal_ts%data_ind_low + 1 !*** unless leap year;
       tidal_ts%data_day_of_year = tidal_ts%data_day_of_year_first
   endif

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   !*** is data frequency ok? only daily frequency supported now
   !*** is data record too large? increase tidal_ts_DATA_max if necessary
   !*** is data calendar supported? only fixed 365-day or gregorian data
   !      calendars are supported

   !*** tidal ts limitation: for now, only daily coupling is supported
   if (fit_freq .ne. 1) then
     call document ('init_tidal_ts', &
                    'temporal tidal variations supported with daily coupling only')
    !call shr_sys_abort ('invalid coupling frequency')
   endif

   !*** is model timestep so small that the number of timesteps exceeds 
   !    the size/dimension of the daily wieghts?
   if (nsteps_per_day > tidal_ts_NSTEPS_max) then
     call document ('init_tidal_ts','nsteps_per_day',nsteps_per_day)
     call document ('init_tidal_ts','tidal_ts_NSTEPS_max',tidal_ts_NSTEPS_max)
     call shr_sys_abort ('must increase tidal_ts_NSTEPS_max')
   endif

   !***  is there a data/model calendar mismatch? Model can handle this:
   if (trim(tidal_ts%data_calendar) == 'gregorian' .and. .not. allow_leapyear) then
     tidal_ts%calendar_mismatch = .true.
   else
     tidal_ts%calendar_mismatch = .false.
   endif

   !***  is there a data/model calendar mismatch? Model CANNOT handle this:
   if (trim(tidal_ts%data_calendar) == '365' .and. allow_leapyear) then
     call document ('init_tidal_ts','FATAL data calendar mismatch ')
     call document ('init_tidal_ts','data calendar = ',trim(tidal_ts%data_calendar))
     call document ('init_tidal_ts','model calendar allow leapyear ', allow_leapyear)
     call shr_sys_abort ('must adjust calendars')
   endif

!-----------------------------------------------------------------------
!  pre-compute fixed daily weights
!-----------------------------------------------------------------------
   dayfrac(0) = c0

   do nn=1,nsteps_per_day
     dayfrac(nn)     = dayfrac(nn-1) + stepsize/seconds_in_day
     tidal_ts%w1(nn) = c1 - dayfrac(nn)
     tidal_ts%w2(nn) =      dayfrac(nn)
   enddo ! n
   tidal_ts%w1(nsteps_per_day) = c0 !eliminate rounding/accumulation error
   tidal_ts%w2(nsteps_per_day) = c1 !eliminate rounding/accumulation error

   write(stdout,*) ' nn  dayfrac  w1  w2'
   do nn=1,nsteps_per_day
     write(stdout,'(1x,i2.2,1x,3(1pe25.10,1x))') nn, dayfrac(nn), tidal_ts%w1(nn), tidal_ts%w2(nn)
   enddo ! n

!-----------------------------------------------------------------------
!  form modulation factors by overwriting
!      tidal_ts%data with 1.0 + tidal_ts%data/100.0  
!-----------------------------------------------------------------------

   tidal_ts%data = 1.0 + 0.01_r8*tidal_ts%data

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------
   if (ltidal_lunar_cycle)  &
      call get_timer(timer_tidal_ts,'TIDAL_TS_DRIVER', &
                                  nblocks_clinic, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  document tidal timeseries variables in output log file
!-----------------------------------------------------------------------
   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' Tidal mixing 18.6-year lunar tidal cycle information'
     write(stdout,blank_fmt)
     do ii=1,ntidal_ts_TC
     write(stdout,*) '     tidal ts file                 = ', trim(tidal_ts%filenames(ii))
     enddo
     write(stdout,*) '     tidal ts file fmt             = ', trim(tidal_ts%fmt)
     write(stdout,*) '     tidal ts calendar             = ', trim(tidal_ts%data_calendar)
     write(stdout,*) '     tidal ts calendar mismatch    = ', tidal_ts%calendar_mismatch
     write(stdout,*) '     tidal ts model year alignment = ', tidal_ts%model_first_year
     write(stdout,*) '     tidal ts model allow_leapyear = ', allow_leapyear
     write(stdout,*) '     tidal ts first data year      = ', tidal_ts%data_first_year
     write(stdout,*) '     tidal ts first data month     = ', tidal_ts%data_first_month
     write(stdout,*) '     tidal ts first data day       = ', tidal_ts%data_first_day
     write(stdout,*) '     tidal_ts data day of year     = ', tidal_ts%data_day_of_year
     write(stdout,*) '     tidal ts last  data year      = ', tidal_ts%data_final_year
     write(stdout,*) '     tidal ts last  data month     = ', tidal_ts%data_final_month
     write(stdout,*) '     tidal ts last  data day       = ', tidal_ts%data_final_day 
     write(stdout,*) '     tidal_ts%data_ind_first       = ', tidal_ts%data_ind_first
     write(stdout,*) '     tidal_ts%data_ind_final       = ', tidal_ts%data_ind_final
     write(stdout,*) '     tidal_ts%data_npoints         = ', tidal_ts%data_npoints
     write(stdout,blank_fmt)

     write(stdout,*) '     tidal ts data points'
     write(stdout,*) '     nn        yyyy-mm-dd        data'
     do nn=tidal_ts%data_ind_first,tidal_ts%data_ind_final
      !1200 format (5x,i7,2x,i4.4,a,i2.2,a,i2.2,2x,4(1pe25.15))
       1200 format (5x,i7,2x,i4.4,a,i2.2,a,i2.2,2x,4f25.15)
       write(stdout,1200) nn, tidal_ts%data_year(nn),'-',  &
         tidal_ts%data_month(nn), '-', tidal_ts%data_day(nn), &
         (tidal_ts%data(nn,nnn),nnn=1,ntidal_ts_TC)
     enddo
     write(stdout,blank_fmt)

     write(stdout,*) '     tidal ts daily weights  '
     write(stdout,*) '     nn             w1               w2              sum'
     do nn=1,nsteps_per_day
       write(stdout,'(5x,i4,3(1pe20.10))')  &
           nn,tidal_ts%w1(nn),tidal_ts%w2(nn),tidal_ts%w1(nn)+tidal_ts%w2(nn)
     enddo ! nn
     write(stdout,blank_fmt)
   endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine init_tidal_ts

!***********************************************************************
!BOP
! !IROUTINE: init_tidal_tavg
! !INTERFACE:

 subroutine init_tidal_tavg

! !DESCRIPTION: This subroutine initializes various ids for time-averaged
!               output
   

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len_long) ::  string   ! dummy string
   character (char_len_long) ::  string1  ! dummy string
   character (char_len_long) ::  string2  ! dummy string

   character (char_len) :: field_name  ! dummy field name
   character (char_len) :: substring   ! dummy sring

   integer (int_kind) :: tavg_method

!-----------------------------------------------------------------------
!
!  set up support for writing out time-independent TIDAL_QE_3D
!  via the tavg mechanism (put into "once" stream)
!
!-----------------------------------------------------------------------

   !*** pick one from column A... (tidal mixing method)
   select case (tidal_mixing_method_itype)
    case (tidal_mixing_method_jayne)
      string = 'Jayne'
    case (tidal_mixing_method_polzin)
      string = 'Polzin'
    case (tidal_mixing_method_schmittner)
      string = 'Schmittner'
    case default
      string = 'unknown'
   end select
   string = trim(string)//' Tidal Mixing Method with '

   !*** and one from column B... (tidal energy file)
   select case (tidal_energy_itype)
    case (tidal_energy_jayne)
      string = trim(string)//'Jayne'
    case (tidal_energy_arbic)
      string = trim(string)//'Arbic'
    case (tidal_energy_ER03)
      string = trim(string)//'Egbert & Ray'
    case (tidal_energy_GN13)
      string = trim(string)//'Green & Nycander'
    case (tidal_energy_LGM0)
      string = trim(string)//'LGM0 Wilmes active research do not use without author permission'
    case (tidal_energy_LGMi5g21)
      string = trim(string)//'LGMi5g21 Wilmes active research do not use without author permission' 
    case (tidal_energy_LGMi6g21)
      string = trim(string)//'LGMi6g21 Wilmes active research do not use without author permission' 
    case default
      string = trim(string)//'unknown'
   end select
   string = trim(string)//' Tidal Energy'

   field_name = 'TIDAL_QE_2D'
   if (ltidal_lunar_cycle) then
     tavg_method = tavg_method_avg
   else
     tavg_method = tavg_method_constant
   endif
   string = trim(field_name)//trim(string)
   call define_tavg_field(tavg_TIDAL_QE_2D,trim(field_name),2, &
                          tavg_method=tavg_method,             &
                          long_name=trim(string),              &
                          units='g/s^3',                       &
                          grid_loc='2110',                     &
                          coordinates  ='TLONG TLAT'           ) 

   field_name = 'TIDAL_QE_3D'
   if (ltidal_lunar_cycle) then
     tavg_method = tavg_method_avg
   else
     tavg_method = tavg_method_constant
   endif
   string = trim(field_name)//trim(string)
   call define_tavg_field(tavg_TIDAL_QE_3D,trim(field_name),3,    &
                          long_name=trim(string),                 &
                          tavg_method=tavg_method,                &
                          units='g/s^3',                          &
                          grid_loc='3111',                        &
                          coordinates  ='TLONG TLAT z_w time'     ) 

   field_name = 'TIDAL_COEF_3D'
   if (ltidal_lunar_cycle) then
     tavg_method = tavg_method_avg
   else
     tavg_method = tavg_method_constant
   endif
   string = trim(field_name)//' tidal coefficients'
   call define_tavg_field(tavg_TIDAL_COEF_3D,trim(field_name),3,  &
                          long_name=trim(string),                 &
                          tavg_method=tavg_method,                &
                          units='unknown',                        &
                          grid_loc='3111',                        &
                          coordinates  ='TLONG TLAT z_w time'     ) 
   field_name = 'VERTICAL_FUNC'
   string = trim(field_name)//' vertical function'
   call define_tavg_field(tavg_VERTICAL_FUNC,trim(field_name),3,  &
                          long_name=trim(string),                 &
                          tavg_method=tavg_method_constant,       &
                          units='unknown',                        &
                          grid_loc='3111',                        &
                          coordinates  ='TLONG TLAT z_w time'     ) 

!-----------------------------------------------------------------------
!
!  set up support for accumulating TIDAL_DIFF and TIDAL_GAMMA_EPS
!
!-----------------------------------------------------------------------

   !*** tidal mixing method
   select case (tidal_mixing_method_itype)
    case (tidal_mixing_method_jayne)
      string1 = 'Jayne'
    case (tidal_mixing_method_polzin)
      string1 = 'Polzin'
    case (tidal_mixing_method_schmittner)
      string1 = 'Schmittner'
    case default
      string1 = 'unknown'
   end select

   string  = trim(string1)//' Tidal Diffusion'
   field_name = 'TIDAL_DIFF'
   call define_tavg_field(tavg_TIDAL_DIFF,trim(field_name),3,     &
                          long_name=trim(string),                 &
                          units=' ',                              &
                          grid_loc='3113',                        &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   string = trim(string1)//' Tidal Gamma*epsilon'
   field_name = 'TIDAL_GAMMA_EPS'
   call define_tavg_field(tavg_TIDAL_GAMMA_EPS,trim(field_name),3,&
                          long_name=trim(string),                 &
                          units=' ',                              &
                          grid_loc='3113',                        &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   string = trim(string1)//' N^2'
   field_name = 'TIDAL_N2'
   call define_tavg_field(tavg_TIDAL_N2,trim(field_name),3,       &
                          long_name=trim(string),                 &
                          units=' ',                              &
                          grid_loc='3113',                        &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   string = trim(string1)//' max(N^2,tidal_eps_n2)'
   field_name = 'TIDAL_N2_eps'
   call define_tavg_field(tavg_TIDAL_N2_eps,trim(field_name),3,   &
                          long_name=trim(string),                 &
                          units=' ',                              &
                          grid_loc='3113',                        &
                          coordinates  ='TLONG TLAT z_w_bot time' ) 

   select case (tidal_energy_itype)

      case (tidal_energy_arbic,tidal_energy_jayne)
           !=====================================

        field_name = 'TIDAL_ENERGY_FLUX_2D'

        select case (tidal_energy_itype)
          case (tidal_energy_jayne)
            string = trim(field_name)//' Jayne Tidal Energy Flux'
          case (tidal_energy_arbic)
            string = trim(field_name)//' Arbic Tidal Energy Flux'
        end select ! tidal_energy_itype


        call define_tavg_field(tavg_TIDAL_ENERGY_FLUX_2D,trim(field_name),2,     &
                               tavg_method=tavg_method_constant,  &
                               long_name=trim(string),            &
                               units='g/s^3',                     &
                               grid_loc='2110',                   &
                               coordinates  ='TLONG TLAT'         ) 

           !===================================
      case (tidal_energy_ER03,tidal_energy_GN13,  &
            tidal_energy_LGM0,tidal_energy_LGMi5g21,tidal_energy_LGMi6g21)
           !============================================================

        select case (tidal_energy_itype)
          case (tidal_energy_ER03)
            substring = 'ER03'
          case (tidal_energy_GN13)
            substring = 'GN13'
          case (tidal_energy_LGM0)
            substring = 'LGM0'
          case (tidal_energy_LGMi5g21)
            substring = 'LGMi5g21'
          case (tidal_energy_LGMi6g21)
            substring = 'LGMi6g21'
        end select

        !---------------------------------------------------------------
        ! define ER03- GN13- or LGM-specific tavg output variables
        !---------------------------------------------------------------

        field_name = 'TCM2'
        string = trim(field_name)//' '//trim(substring)//' semidiurnal lunar'
        call define_tavg_field(tavg_TCM2,trim(field_name),3,      &
                               tavg_method=tavg_method_constant,  &
                               long_name=trim(string),            &
                               units='W/m^2',                     &
                               grid_loc='3111',                   &
                               coordinates  ='TLONG TLAT  z_t'    ) 

        field_name = 'TCS2'
        string = trim(field_name)//' '//trim(substring)//' semidiurnal solar'
        call define_tavg_field(tavg_TCS2,trim(field_name),3,      &
                               tavg_method=tavg_method_constant,  &
                               long_name=trim(string),            &
                               units='W/m^2',                     &
                               grid_loc='3111',                   &
                               coordinates  ='TLONG TLAT  z_t'    ) 

        field_name = 'TCK1'
        string = trim(field_name)//' '//trim(substring)//' diurnal lunar'
        call define_tavg_field(tavg_TCK1,trim(field_name),3,      &
                               tavg_method=tavg_method_constant,  &
                               long_name=trim(string),            &
                               units='W/m^2',                     &
                               grid_loc='3111',                   &
                               coordinates  ='TLONG TLAT  z_t'    ) 

        field_name = 'TCO1'
        string = trim(field_name)//' '//trim(substring)//' diurnal solar'
        call define_tavg_field(tavg_TCO1,trim(field_name),3,      &
                               tavg_method=tavg_method_constant,  &
                               long_name=trim(string),            &
                               units='W/m^2',                     &
                               grid_loc='3111',                   &
                               coordinates  ='TLONG TLAT  z_t'    ) 
   end select


   select case (tidal_mixing_method_itype)

      case (tidal_mixing_method_jayne,tidal_mixing_method_polzin)
           !====================================================

           !---------------------------------------------------------------
           !
           !  define option-specific tavg output variable
           !
           !---------------------------------------------------------------


           field_name = 'TIDAL_QE_2D'

           select case (tidal_mixing_method_itype)
             case (tidal_mixing_method_jayne)
               string = ': Jayne'
             case (tidal_mixing_method_polzin)
               string = ': Polzin'
           end select ! tidal_mixing_method_itype

           string = trim(field_name)//trim(string)//' Tidal Mixing Method with '

           select case (tidal_energy_itype)
             case (tidal_energy_jayne)
               string = trim(string)//'Jayne'
             case (tidal_energy_ER03)
               string = trim(string)//'ER03 (2D)'
             case (tidal_energy_GN13)
               string = trim(string)//'GN13 (2D)'
             case (tidal_energy_LGM0)
               string = trim(string)//'LGM0 (2D)'
             case (tidal_energy_LGMi5g21)
               string = trim(string)//'LGMi5g21 (2D)'
             case (tidal_energy_LGMi6g21)
               string = trim(string)//'LGMi6g21 (2D)'
             case (tidal_energy_arbic)
               string = trim(string)//'Arbic'
           end select ! tidal_energy_itype

           string = trim(string)//' Tidal Energy'

   end select


   select case (tidal_mixing_method_itype)
      case (tidal_mixing_method_polzin)
           !========================== 

        !---------------------------------------------------------------
        !
        !  define Polzin-scheme-specific tavg output variables
        !
        !---------------------------------------------------------------

        field_name = 'H2_P'
        string = trim(field_name)//' from Jayne initial conditions'
        call define_tavg_field(tavg_H2_P,trim(field_name),2,     &
                               tavg_method=tavg_method_constant, &
                               long_name=trim(string),           &
                               units='cm^2',                     &
                               grid_loc='2110',                  &
                               coordinates  ='TLONG TLAT') 

        field_name = 'U_P'
        string = trim(field_name)//' from Jayne initial conditions'
        call define_tavg_field(tavg_U_P,trim(field_name),2,      &
                               tavg_method=tavg_method_constant, &
                               long_name=trim(string),           &
                               units='cm/s',                     &
                               grid_loc='2110',                  &
                               coordinates  ='TLONG TLAT') 

        field_name = 'POLZIN_EQ2'
        string = trim(field_name)//' from Jayne initial conditions'
        call define_tavg_field(tavg_POLZIN_EQ2,trim(field_name),2,      &
                               tavg_method=tavg_method_constant, &
                               long_name=trim(string),           &
                               units='cm/s',                     &
                               grid_loc='2110',                  &
                               coordinates  ='TLONG TLAT') 

   end select

!-----------------------------------------------------------------------
!
!  set up support for accumulating TEMP fields used to compare with
!   Melet et al figure
!  WARNING: these sums are NOT weighted by thickness
!
!-----------------------------------------------------------------------

     field_name = 'TEMP1_1p5km'
     string = 'Potential Temperature [1km,1.5km]'
     call define_tavg_field(tavg_TEMP1_1p5km,trim(field_name),3,  &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
  
     field_name = 'TEMP1_2km'
     string = 'Potential Temperature [1km,2km]'
     call define_tavg_field(tavg_TEMP1_2km,trim(field_name),3,    &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
  
     field_name = 'TEMP1_3km'
     string = 'Potential Temperature [1km,3km]'
     call define_tavg_field(tavg_TEMP1_3km,trim(field_name),3,    &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
  
     field_name = 'TEMP2_3km'
     string = 'Potential Temperature [2km,3km]'
     call define_tavg_field(tavg_TEMP2_3km,trim(field_name),3,    &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
     field_name = 'TEMP2_4km'
     string = 'Potential Temperature [2km,4km]'
     call define_tavg_field(tavg_TEMP2_4km,trim(field_name),3,    &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
  
     field_name = 'TEMP3p5_6km'
     string = 'Potential Temperature [3.5km,6km]'
     call define_tavg_field(tavg_TEMP3p5_6km,trim(field_name),3,  &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')
  
     field_name = 'TEMP4_6km'
     string = 'Potential Temperature [4km,6km]'
     call define_tavg_field(tavg_TEMP4_6km,trim(field_name),3,    &
                            long_name=trim(string),               &
                            units='degC', grid_loc='3111',        &
                            coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!EOC
 
 end subroutine init_tidal_tavg

!***********************************************************************
!BOP
! !IROUTINE: tidal_read_energy_arbic
! !INTERFACE:

 subroutine tidal_read_energy_arbic 

! !DESCRIPTION: This subroutine reads the Arbic 20?? input tidal energy
!               flux dataset.  Presently, this is only a place-holder routine.
   
!EOP
!BOC

!---------------------------------------------------------------
!
!  read Arbic input tidal energy flux
!
!---------------------------------------------------------------

   subname = 'tidal_read_energy_arbic'

   string = 'read '//trim(tidal_energy_file)
   call document (trim(subname),trim(string))

   tidal_mixing_file_in = construct_file(         &
               tidal_energy_file_fmt,             &
               full_name=trim(tidal_energy_file), &
               record_length=rec_type_dbl,        &
               recl_words=nx_global*ny_global     )

   call data_set(tidal_mixing_file_in, 'open_read')

!*** tbd...
!  TIDAL_ENERGY_FLUX_D = construct_io_field(                        &
!              'TIDAL_ENERGY_FLUX', dim1=i_dim, dim2=j_dim,         &
!              long_name='Input Tidal Energy Flux at T-grid Points',&
!              units    ='W/m^2', grid_loc ='2110',                 &
!              field_loc = field_loc_center,                        &
!              field_type = field_type_scalar,                      &
!              d2d_array =TIDAL_ENERGY_FLUX_2D(:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TIDAL_ENERGY_FLUX_D)
   call data_set (tidal_mixing_file_in, 'read'  , TIDAL_ENERGY_FLUX_D)

!---------------------------------------------------------------
!     convert from to g/s^3 if necessary
!---------------------------------------------------------------
   TIDAL_ENERGY_FLUX_2D = c1000*TIDAL_ENERGY_FLUX_2D

   call destroy_io_field (TIDAL_ENERGY_FLUX_D)
   call data_set (tidal_mixing_file_in, 'close')

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_read_energy_arbic

!***********************************************************************
!BOP
! !IROUTINE: tidal_read_energy_jayne
! !INTERFACE:

 subroutine tidal_read_energy_jayne

! !DESCRIPTION: This subroutine reads the Jayne 2009 input tidal energy
!               flux dataset and converts from mks to cgs
   
!EOP
!BOC

!---------------------------------------------------------------
!
!  read Jayne 2009 input tidal energy flux 
!
!---------------------------------------------------------------

   subname = 'tidal_read_energy_jayne'

   string = 'read '//trim(tidal_energy_file)
   call document (trim(subname),trim(string))

   tidal_mixing_file_in = construct_file(         &
               tidal_energy_file_fmt,             &
               full_name=trim(tidal_energy_file), &
               record_length=rec_type_dbl,        &
               recl_words=nx_global*ny_global     )

   call data_set(tidal_mixing_file_in, 'open_read')

   TIDAL_ENERGY_FLUX_D = construct_io_field(                        &
               'TIDAL_ENERGY_FLUX', dim1=i_dim, dim2=j_dim,         &
               long_name='Input Tidal Energy Flux at T-grid Points',&
               units    ='W/m^2', grid_loc ='2110',                 &
               field_loc = field_loc_center,                        &
               field_type = field_type_scalar,                      &
               d2d_array =TIDAL_ENERGY_FLUX_2D(:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TIDAL_ENERGY_FLUX_D)
   call data_set (tidal_mixing_file_in, 'read'  , TIDAL_ENERGY_FLUX_D)

!---------------------------------------------------------------
!     convert from W/m^2 to g/s^3 
!---------------------------------------------------------------
   TIDAL_ENERGY_FLUX_2D = c1000*TIDAL_ENERGY_FLUX_2D

   call destroy_io_field (TIDAL_ENERGY_FLUX_D)
   call data_set (tidal_mixing_file_in, 'close')

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_read_energy_jayne

!***********************************************************************
!BOP
! !IROUTINE: tidal_read_energy_TC
! !INTERFACE:

 subroutine tidal_read_energy_TC (file_id_string)

! !DESCRIPTION: This subroutine reads the 3D tidal constituents from the
!  Egbert & Ray 2003 or the Nycander & Green 2013 input tidal datasets
!  and converts from mks to cgs
 
! !INPUT PARAMETERS:
   character(*), intent(in) :: file_id_string

!EOP
!BOC

!---------------------------------------------------------------
!
!  read ER03, GN13, or LGM input diurnal and semidiurnal solar and
!  lunar tidal energy files and compute time-independent E_TC
!
!  In POP2, the M2,S2,K1,O1 variables are stored in variables
!  named EDR{M2,S2,K1,O1}, even if the input dataset is GN13 or LGM.
!
!---------------------------------------------------------------

   subname = 'tidal_read_energy_TC'

   string = 'read '//trim(tidal_energy_file)
   call document (trim(subname),trim(string))

   allocate (TCM2(nx_block,ny_block,km,nblocks_clinic)) ; TCM2 = c0
   allocate (TCS2(nx_block,ny_block,km,nblocks_clinic)) ; TCS2 = c0
   allocate (TCK1(nx_block,ny_block,km,nblocks_clinic)) ; TCK1 = c0
   allocate (TCO1(nx_block,ny_block,km,nblocks_clinic)) ; TCO1 = c0

   tidal_mixing_file_in = construct_file(         &
               tidal_energy_file_fmt,             &
               full_name=trim(tidal_energy_file), &
               record_length=rec_type_dbl,        &
               recl_words=nx_global*ny_global     )

   call data_set(tidal_mixing_file_in, 'open_read')

   TCM2_D = construct_io_field(                                    &
             'M2', dim1=i_dim, dim2=j_dim,dim3=k_dim,              &
              long_name='Input Semidiurnal Lunar Tidal Energy Flux at T-grid Points',&
              units    ='W/m^2', grid_loc ='3111',                 &
              field_loc = field_loc_center,                        &
              field_type = field_type_scalar,                      &
              d3d_array =TCM2(:,:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TCM2_D)
   call data_set (tidal_mixing_file_in, 'read'  , TCM2_D)
   call destroy_io_field (TCM2_D)

   TCS2_D = construct_io_field(                                    &
             'S2', dim1=i_dim, dim2=j_dim,dim3=k_dim,              &
              long_name='Input Semidiurnal Solar Tidal Energy Flux at T-grid Points',&
              units    ='W/m^2', grid_loc ='3111',                 &
              field_loc = field_loc_center,                        &
              field_type = field_type_scalar,                      &
              d3d_array =TCS2(:,:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TCS2_D)
   call data_set (tidal_mixing_file_in, 'read'  , TCS2_D)
   call destroy_io_field (TCS2_D)

   TCK1_D = construct_io_field(                                    &
             'K1', dim1=i_dim, dim2=j_dim,dim3=k_dim,              &
              long_name='Input Diurnal Lunar Tidal Energy Flux at T-grid Points',&
              units    ='W/m^2', grid_loc ='3111',                 &
              field_loc = field_loc_center,                        &
              field_type = field_type_scalar,                      &
              d3d_array =TCK1(:,:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TCK1_D)
   call data_set (tidal_mixing_file_in, 'read'  , TCK1_D)
   call destroy_io_field (TCK1_D)

   TCO1_D = construct_io_field(                                    &
             'O1', dim1=i_dim, dim2=j_dim,dim3=k_dim,              &
              long_name='Input Diurnal Solar Tidal Energy Flux at T-grid Points',&
              units    ='W/m^2', grid_loc ='3111',                 &
              field_loc = field_loc_center,                        &
              field_type = field_type_scalar,                      &
              d3d_array =TCO1(:,:,:,:))

   call data_set (tidal_mixing_file_in, 'define', TCO1_D)
   call data_set (tidal_mixing_file_in, 'read'  , TCO1_D)
   call destroy_io_field (TCO1_D)

   call data_set (tidal_mixing_file_in, 'close')

!---------------------------------------------------------------
! E&R and N&G fields are in W/m^2. convert to g/s^3
!---------------------------------------------------------------

   TCM2 = c1000*TCM2
   TCS2 = c1000*TCS2
   TCK1 = c1000*TCK1
   TCO1 = c1000*TCO1

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_read_energy_TC

!***********************************************************************
!BOP
! !IROUTINE: tidal_read_roughness_RMS
! !INTERFACE:

 subroutine tidal_read_roughness_RMS

! !DESCRIPTION: This subroutine reads roughness and urms fields needed for
!               the Polzin method
   
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len_long) ::   &
      input_file_name,            &! generic tidal energy filename
      string                       ! dummy string

   character (char_len)      ::   &
      input_file_fmt,             &! generic tidal energy file format string
      subname                      ! subroutine name

   subname = 'tidal_read_roughness_RMS'
   input_file_name = tidal_vars_file_polz
   input_file_fmt  = tidal_vars_file_fmt_polz

   string = 'read '//trim(input_file_name)
   call document (trim(subname),trim(string))

!---------------------------------------------------------------
!
!  read roughness and urms from the input file
!  and compute time-independent parts of tidal mixing method (Polzin)
!
!---------------------------------------------------------------

   !***  read roughness (h)
        !roughness:units = "meters" ;

   tidal_mixing_file_in = construct_file(          &
               input_file_fmt,                     &
               full_name=trim(input_file_name),    &
               record_length=rec_type_dbl,         &
               recl_words=nx_global*ny_global      )

   call data_set(tidal_mixing_file_in, 'open_read')

   H2_P_D = construct_io_field(                             &
               'roughness', dim1=i_dim, dim2=j_dim,        &
               long_name='Input Roughness at T-grid Points',&
               units    ='m', grid_loc ='2110',             &
               field_loc = field_loc_center,                &
               field_type = field_type_scalar,              &
               d2d_array = H2_P(:,:,:))


   call data_set (tidal_mixing_file_in, 'define', H2_P_D)
   call data_set (tidal_mixing_file_in, 'read'  , H2_P_D)

   call destroy_io_field (H2_P_D)

   call document (trim(subname), 'min H2_P (m) ', minval(H2_P))
   call document (trim(subname), 'max H2_P (m) ', maxval(H2_P))

  !***  read U_P(urms)
        !urms:units = "meters/second" ;
   U_P_D = construct_io_field(                         &
               'urms', dim1=i_dim, dim2=j_dim,        &
               long_name='Input URMS at T-grid Points',&
               units    ='m/s', grid_loc ='2110',      &
               field_loc = field_loc_center,           &
               field_type = field_type_scalar,         &
               d2d_array =U_P(:,:,:))

   call data_set (tidal_mixing_file_in, 'define', U_P_D)
   call data_set (tidal_mixing_file_in, 'read'  , U_P_D)

   call destroy_io_field (U_P_D)
   call data_set (tidal_mixing_file_in, 'close')

   call document (trim(subname), 'min U_P (m/s)', minval(U_P))
   call document (trim(subname), 'max U_P (m/s)', maxval(U_P))

   !***  convert from m/s to cm/s 
   U_P = 100.0_r8*U_P*RCALCT(:,:,1:nblocks_clinic)

   !***  convert from m to cm
   H2_P = 100.0_r8*RCALCT(:,:,1:nblocks_clinic)*H2_P

   !***  create h**2  (roughness = h)
   H2_P = H2_P**2

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_read_roughness_RMS

!***********************************************************************
!BOP
! !IROUTINE: tidal_form_coef_jayne
! !INTERFACE:

 subroutine tidal_form_coef_jayne

! !DESCRIPTION:
!  Form TIDAL_COEF_3D when using Jayne method
!  Note: this routine is called once, unless ltidal_lunar_cycle=.true.
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  iblock  ! block index
   integer (int_kind) ::  k       ! vertical loop index

   TIDAL_COEF_2D = tidal_gamma_rhor*RCALCT(:,:,1:nblocks_clinic)*TIDAL_QE_2D
   TIDAL_COEF_3D = c0

   do iblock = 1,nblocks_clinic

     do k=1,km
       where ( k <= KMT(:,:,iblock) )
         TIDAL_COEF_3D(:,:,k,iblock) = TIDAL_COEF_2D(:,:,iblock)*VERTICAL_FUNC(:,:,k,iblock)
       endwhere
     enddo ! k

   enddo ! iblock

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_form_coef_jayne

!***********************************************************************
!BOP
! !IROUTINE: tidal_form_coef_schm
! !INTERFACE:

 subroutine tidal_form_coef_schm

! !DESCRIPTION:
!  Form TIDAL_COEF_3D using Schmittner method
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r8), dimension(:,:,:,:,:), allocatable :: EXP_HAB_ZETAR

   integer (int_kind) :: iblock    ! block index
   integer (int_kind) :: i,j,k,k1  ! loop indices

   logical (log_kind) :: first_call = .true.

   save EXP_HAB_ZETAR
   save first_call

   if (first_call) then
     allocate (EXP_HAB_ZETAR(nx_block,ny_block,km,km,nblocks_clinic)) 
     EXP_HAB_ZETAR = c0

     do iblock = 1,nblocks_clinic
       do j=1,ny_block
         do i=1,nx_block
           do k=1,KMT(i,j,iblock)-1
             do k1=k+1,KMT(i,j,iblock)
               hab = zw(k) - zw(k1) !height above bottom
               EXP_HAB_ZETAR(i,j,k,k1,iblock) = exp(hab*zetar)*decay_fn(k1)
             enddo ! k1
           enddo ! k
         enddo !i
       enddo !j
     enddo !iblock

     EXP_HAB_ZETAR = tidal_gamma_rhor*EXP_HAB_ZETAR

     first_call = .false.
   endif

   TIDAL_COEF_3D = c0
   do iblock = 1,nblocks_clinic
     do j=1,ny_block
       do i=1,nx_block
         do k=1,KMT(i,j,iblock)-1

           !***subgrid-scale scheme: sum over all levels below k
           do k1=k+1,KMT(i,j,iblock)
             TIDAL_COEF_3D(i,j,k,iblock) = TIDAL_COEF_3D(i,j,k,iblock) +  &
                  TIDAL_QE_3D(i,j,k1,iblock)*EXP_HAB_ZETAR(i,j,k,k1,iblock)
           enddo ! k1
         enddo ! k
       enddo !i
     enddo !j
   enddo !iblock

  !TIDAL_COEF_3D = tidal_gamma_rhor * TIDAL_COEF_3D  !included in EXP_HAB_ZETAR above

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_form_coef_schm

!***********************************************************************
!BOP
! !IROUTINE: tidal_form_qE_2D
! !INTERFACE:

 subroutine tidal_form_qE_2D (tidal_local_mixing_fraction)

! !DESCRIPTION:
!  Form q*TIDAL_ENERGY_FLUX_2D from tidal constituents 
!
! !REVISION HISTORY:
!  njn01 2016

! !INPUT:
   real (r8), intent(in) :: tidal_local_mixing_fraction  ! "q"

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: iblock ! block index
      
!---------------------------------------------------------------
!     Form q*TIDAL_ENERGY_FLUX_2D
!---------------------------------------------------------------

   TIDAL_QE_2D = tidal_local_mixing_fraction*TIDAL_ENERGY_FLUX_2D

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_form_qE_2D

!***********************************************************************
!BOP
! !IROUTINE: tidal_form_qE_2D_from_3D
! !INTERFACE:

 subroutine tidal_form_qE_2D_from_3D

! !DESCRIPTION:
!  Collapse 3D q*E to 2D q*E
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: iblock,k  ! block and vertical indices
      
   !-----------------------------------------------------------
   !  collapse 3D q*E to 2D q*E 
   !-----------------------------------------------------------

   TIDAL_QE_2D(:,:,:) = c0

   do iblock = 1,nblocks_clinic
     do k=1,km
       where (k <= KMT(:,:,iblock)) 
        TIDAL_QE_2D(:,:,iblock) = TIDAL_QE_2D(:,:,iblock) + TIDAL_QE_3D(:,:,k,iblock)
       endwhere
     enddo ! k
   enddo ! iblock

  !-----------------------------------------------------------------------
  !  Test: is qE positive? 
  !-----------------------------------------------------------------------
   if (any(TIDAL_QE_2D < c0)) call shr_sys_abort ('tidal_form_qE_2D_from_3D: negative TIDAL_QE_2D terms')
   !*** do not deallocate (TIDAL_QE_3D) here

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_form_qE_2D_from_3D

!***********************************************************************
!BOP
! !IROUTINE: tidal_form_qE_3D
! !INTERFACE:

 subroutine tidal_form_qE_3D (cm2,cs2,ck1,co1)

! !DESCRIPTION:
!  Form TIDAL_QE_3D from weighted tidal constituents 
!
! !REVISION HISTORY:
!  njn01 2016

! !INPUT PARAMETERS:
   real(r8), intent(in) :: cm2,cs2,ck1,co1

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: iblock  ! block index
   integer (int_kind) :: k       ! vertical level index

   logical (log_kind) :: first_call = .true.
      
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      WORK                         ! WORK array

   real (r8), dimension(nx_block,ny_block) ::  &
      TIDAL_QK1,                  &! qk1 coefficient used in Schmittner & Egbert
      TIDAL_QO1                    ! qo1 coefficient used in Schmittner & Egbert

   real (r8), dimension(nx_block,ny_block) ::  &
      ck1TIDAL_QK1,                  &! time-varying ck1 coefficient*TIDAL_QK1
      co1TIDAL_QO1                    ! time-varying qo1 coefficient*TIDAL_QO1

   real (r8) ::  &
    p33cm2,  &
    p33cs2

   save first_call
   save TIDAL_QK1
   save TIDAL_QO1
      
   if (.not. allocated (TIDAL_QE_3D)) &
     allocate (TIDAL_QE_3D(nx_block,ny_block,km,nblocks_clinic)) 

   TIDAL_QE_3D(:,:,:,:) = c0

   do iblock = 1,nblocks_clinic

    if (first_call) then
       where (abs(TLATD(:,:,iblock)) < 30._r8 )
           TIDAL_QK1(:,:) = p33
           TIDAL_QO1(:,:) = p33
       elsewhere
           TIDAL_QK1(:,:) = c1
           TIDAL_QO1(:,:) = c1
       endwhere
    endif

    ck1TIDAL_QK1 = ck1*TIDAL_QK1
    co1TIDAL_QO1 = co1*TIDAL_QO1

    p33cm2 = p33*cm2
    p33cs2 = p33*cs2

!-----------------------------------------------------------------------
!  Note: zero dissipation contribution above tidal_diss_lim_TC 
!-----------------------------------------------------------------------
     if (ltidal_all_TC_coefs_eq_1) then  ! for checking/plotting purposes
       do k=1,km
         where ( k <= KMT(:,:,iblock) .and. zw(k) > tidal_diss_lim_TC )
           TIDAL_QE_3D(:,:,k,iblock) =   &
              cm2*TCM2(:,:,k,iblock) + cs2*TCS2(:,:,k,iblock) +  &
              ck1*TCK1(:,:,k,iblock) + co1*TCO1(:,:,k,iblock)
         endwhere
       enddo ! k
     elseif (ltidal_all_TC_coefs_eq_p33) then  ! for checking/plotting purposes
       do k=1,km
         where ( k <= KMT(:,:,iblock) .and. zw(k) > tidal_diss_lim_TC )
           TIDAL_QE_3D(:,:,k,iblock) =  &
                  (cm2*TCM2(:,:,k,iblock) + cs2*TCS2(:,:,k,iblock) +  &
                   ck1*TCK1(:,:,k,iblock) + co1*TCO1(:,:,k,iblock))*p33
         endwhere
       enddo ! k
     else  ! standard production version
       do k=1,km
         where ( k <= KMT(:,:,iblock) .and. zw(k) > tidal_diss_lim_TC )
           TIDAL_QE_3D(:,:,k,iblock) =   &
                        p33cm2*TCM2(:,:,k,iblock) +            p33cs2*TCS2(:,:,k,iblock) +  &
             ck1TIDAL_QK1(:,:)*TCK1(:,:,k,iblock) + co1TIDAL_QO1(:,:)*TCO1(:,:,k,iblock)

          ! TIDAL_QE_3D(:,:,k,iblock) =   &
          !           p33*(cm2*TCM2(:,:,k,iblock) + cs2*TCS2(:,:,k,iblock)) +  &
          ! TIDAL_QK1(:,:)*ck1*TCK1(:,:,k,iblock) + co1*TIDAL_QO1(:,:)*TCO1(:,:,k,iblock)
         endwhere
       enddo ! k
     endif ! ltidal_all_TC_coefs_eq_1

   enddo ! iblock

!-----------------------------------------------------------------------
!  Document input tidal-constituent coefficients
!-----------------------------------------------------------------------
  if (first_call) then
    !*** this prints from all threads
    !*** note: coefficients vary in time when using LNC
    call document ('tidal_form_qE_3D', 'cm2', cm2)
    call document ('tidal_form_qE_3D', 'cs2', cs2)
    call document ('tidal_form_qE_3D', 'ck1', ck1)
    call document ('tidal_form_qE_3D', 'co1', co1)
  endif

!-----------------------------------------------------------------------
!  Test: is qE positive? 
!-----------------------------------------------------------------------
   if (any(TIDAL_QE_3D < c0)) call shr_sys_abort ('tidal_form_qE_3D: negative TIDAL_QE_3D terms')

  !*** do not deallocate here
  !if (.not. ltidal_lunar_cycle) deallocate (TIDAL_QK1, TIDAL_QO1)


   first_call = .false.

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_form_qE_3D

!***********************************************************************
!BOP
! !IROUTINE: tidal_check
! !INTERFACE:

 subroutine tidal_check

! !DESCRIPTION:
!  Perform tidal-mixing module consistency checking
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len_long) :: string  ! dummy string

   integer (int_kind)        :: number_of_fatal_errors ! counter for number of fatal error conditions


!-----------------------------------------------------------------------
!
!  exit and error conditions
!
!-----------------------------------------------------------------------

   !--------------------------------------------------------------------
   !  exit if tidal mixing is not enabled
   !--------------------------------------------------------------------
   if (.not. ltidal_mixing) return

   number_of_fatal_errors = 0
   call document ('tidal_check','Begin tidal-mixing error checking')

   !--------------------------------------------------------------------
   !  abort if partial_bottom_cells is enabled when polzin or schmittner
   !  scheme is to be used
   !--------------------------------------------------------------------
   if ( partial_bottom_cells .and.&
        (tidal_mixing_method_itype .eq. tidal_mixing_method_schmittner.or.&
         tidal_mixing_method_itype .eq. tidal_mixing_method_polzin) ) then
      string = 'ERROR: partial_bottom_cells not applicable for Schmittner'&
                //' or Polzin tidal_mixing option'
      call document ('tidal_check', trim(string))
      call document ('tidal_check', 'tidal_mixing_method_itype', tidal_mixing_method_itype)
      call document ('tidal_check', 'partial_bottom_cells', partial_bottom_cells)
      number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   !--------------------------------------------------------------------
   !  abort if unknown tidal mixing option
   !--------------------------------------------------------------------
   if (tidal_mixing_method_itype == -1000) then
      string = 'ERROR: unknown tidal mixing type option'
      call document ('tidal_check', trim(string))
      call document ('tidal_check', 'tidal_mixing_method_itype', tidal_mixing_method_itype)
      number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   !--------------------------------------------------------------------
   !  abort if unknown Schmittner vertical decay option
   !--------------------------------------------------------------------
   if (trim(tidal_vert_decay_option_schm) .ne. 'SSJ02' .and.           &
       trim(tidal_vert_decay_option_schm) .ne. 'P09' ) then
      string = 'ERROR: unknown tidal_vert_decay_option_schm option'
      call document ('tidal_check', trim(string))
      call document ('tidal_check', 'tidal_vert_decay_option_schm ',   &
                                trim(tidal_vert_decay_option_schm))
      number_of_fatal_errors = number_of_fatal_errors + 1
   endif

   !--------------------------------------------------------------------
   !  abort if tidal mixing option/tidal energy file mismatch
   !--------------------------------------------------------------------
     if (tidal_mixing_method_itype == tidal_mixing_method_jayne .or.  &
         tidal_mixing_method_itype == tidal_mixing_method_polzin      ) then
      !*** allowable energy file choices are: arbic, jayne, ER03, GN13, or LGM
      if (tidal_energy_itype == tidal_energy_arbic     .or.  &
          tidal_energy_itype == tidal_energy_jayne     .or.  &
          tidal_energy_itype == tidal_energy_ER03      .or.  &
          tidal_energy_itype == tidal_energy_GN13      .or.  &
          tidal_energy_itype == tidal_energy_LGM0      .or.  &
          tidal_energy_itype == tidal_energy_LGMi5g21  .or.  &
          tidal_energy_itype == tidal_energy_LGMi6g21) then
          !*** ok
      else
          string = 'ERROR: mismatch of tidal mixing method and energy type selections'
          call document ('tidal_check', trim(string))
          call document ('tidal_check','tidal_mixing_method_itype', tidal_mixing_method_itype)
          call document ('tidal_check','tidal_energy_itype', tidal_energy_itype)
          number_of_fatal_errors = number_of_fatal_errors + 1
      endif
     endif

     if (tidal_mixing_method_itype == tidal_mixing_method_schmittner) then
      !*** allowable energy file choices are schmittner, Green & Nycander or LGM
       if (tidal_energy_itype == tidal_energy_ER03      .or. &
           tidal_energy_itype == tidal_energy_GN13      .or. &
           tidal_energy_itype == tidal_energy_LGM0      .or. &
           tidal_energy_itype == tidal_energy_LGMi5g21  .or. &
           tidal_energy_itype == tidal_energy_LGMi6g21) then
           !*** ok
       else
          string = 'ERROR: mismatch of tidal mixing method and energy type selections'
          call document ('tidal_check', trim(string))
          call document ('tidal_check','tidal_mixing_method_itype', tidal_mixing_method_itype)
          call document ('tidal_check','tidal_energy_itype', tidal_energy_itype)
          number_of_fatal_errors = number_of_fatal_errors + 1
       endif
     endif

     if (ltidal_lunar_cycle) then
       !*** allowable energy file choices are: ER03, GN13, and LGM
       if (tidal_energy_itype == tidal_energy_ER03      .or. &
           tidal_energy_itype == tidal_energy_GN13      .or. &
           tidal_energy_itype == tidal_energy_LGM0      .or. &
           tidal_energy_itype == tidal_energy_LGMi5g21  .or. &
           tidal_energy_itype == tidal_energy_LGMi6g21) then
         !*** ok
     else
          string = 'ERROR: mismatch of tidal energy file and lunar cycle option'
          call document ('tidal_check', trim(string))
          call document ('tidal_check','tidal_energy_itype', tidal_energy_itype)
          call document ('tidal_check','ltidal_lunar_cycle', ltidal_lunar_cycle)
          number_of_fatal_errors = number_of_fatal_errors + 1
     endif
     endif

     
   !--------------------------------------------------------------------
   !  abort if tidal mixing option and suboptions mismatch
   !--------------------------------------------------------------------
!######## DEBUG disable for now, fix later. ltidal_schmittner_socn works with 'jayne' method, too
!    if (ltidal_schmittner_socn) then
!      if (tidal_mixing_method_itype == tidal_mixing_method_schmittner) then
!       !*** ok
!      else
!         string = 'ERROR: ltidal_schmittner_socn can only be used with schmittner method'
!         call document ('tidal_check', trim(string))
!         call document ('tidal_check','tidal_mixing_method_itype', tidal_mixing_method_itype)
!         call document ('tidal_check','ltidal_schmittner_socn',    ltidal_schmittner_socn)
!         number_of_fatal_errors = number_of_fatal_errors + 1
!     endif
!    endif
!######## DEBUG disable for now

   !--------------------------------------------------------------------
   !  abort if tidal energy file not yet available
   !--------------------------------------------------------------------
   if (tidal_energy_itype == tidal_energy_arbic) then
      string = 'ERROR: Arbic tidal energy file is not yet available'
      call document ('tidal_check', trim(string))
      number_of_fatal_errors = number_of_fatal_errors + 1
   endif

!-----------------------------------------------------------------------
!
!  Now that error messages have been written, stop if there are fatal errors
!
!-----------------------------------------------------------------------

 
   if (number_of_fatal_errors > 0 ) then
     string = 'correct the Tidal-Mixing error condition(s) listed above before continuing'
     call exit_POP (sigAbort, trim(string))
   else
     string = 'No fatal tidal-mixing error conditions detected'
     call document ('tidal_check', string)
   endif

!-----------------------------------------------------------------------
!
!  Now report non-fatal conceptual mismatch and continue execution
!
!-----------------------------------------------------------------------

     if (ltidal_melet_plot) then
      if (tidal_mixing_method_itype .ne. tidal_mixing_method_polzin) then
          string = 'WARNING: mismatch of tidal mixing method and tavg field selections; normally not used together'
          call document ('tidal_check', trim(string))
          call document ('tidal_check','tidal_mixing_method_itype', tidal_mixing_method_itype)
          call document ('tidal_check','ltidal_melet_plot', ltidal_melet_plot)
      endif
     endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_check

!***********************************************************************
!BOP
! !IROUTINE: tidal_compute_diff
! !INTERFACE:

 subroutine tidal_compute_diff (WORK1,k,bid,this_block)

! !DESCRIPTION:
!  Compute tidal diffusion
!
! !REVISION HISTORY:
!  Jayne, S.R.  original version
!  subroutinized and moved from vmix_kpp.F90 by njn01
!  generalized to support Schmittner method as well as Jayne's method

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block), intent(in) :: WORK1        

   integer (int_kind), intent(in) :: k    ! vertical index
   integer (int_kind), intent(in) :: bid  ! block id

   type (block), intent(in)  :: this_block ! block information for current block

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: i,j

!-----------------------------------------------------------------------
!
!  compute tidal diffusion at k
!
!-----------------------------------------------------------------------

   where (WORK1 > c0)
     TIDAL_DIFF(:,:,k,bid) = TIDAL_COEF_3D(:,:,k,bid)/WORK1
   endwhere
    
!-----------------------------------------------------------------------
!
! high mixing in Southern Ocean below 500m as in observations. Schmittner only
!    coded this way for comparability to columnized (cvmix) version
!
!-----------------------------------------------------------------------

   if (ltidal_schmittner_socn)  then

      where (k <= KMT(:,:,bid)-1)
        TIDAL_DIFF(:,:,k,bid) = max(TIDAL_DIFF(:,:,k,bid),TANH_SCHM_3D(:,:,k,bid))
      endwhere

   endif 

!-----------------------------------------------------------------------
!
!  impose maximum tidal diffusion in same sequence as is 
!    done in the columnized (cvmix) version
!
!-----------------------------------------------------------------------

   if (ltidal_max) &
     TIDAL_DIFF(:,:,k,bid) = min(TIDAL_DIFF(:,:,k,bid),tidal_mix_max)


!-----------------------------------------------------------------------
!
!  apply computational stability controls in the same sequence as is
!    done in the columnized (cvmix) version
!
!-----------------------------------------------------------------------

   if (ltidal_stabc .and. .not. lccsm_control_compatible .and. k.gt.2) then
     do j=1,ny_block
     do i=1,nx_block
       if (k == KMT(i,j,bid)-1 .or. k == KMT(i,j,bid)-2) &
         TIDAL_DIFF(i,j,k,bid) = max(TIDAL_DIFF(i,j,k,bid),TIDAL_DIFF(i,j,k-1,bid))
     end do ! i
     end do ! j
   end if ! ltidal_stabc

!-----------------------------------------------------------------------
!
!  impose minimum tidal diffusion at depth in the Denmark Strait to control
!   instability in the same sequence as is done in the columnized (cvmix) version
!
!-----------------------------------------------------------------------

   if (ltidal_min_regions) &
     call tidal_min_regions_set(1,nx_block,1,ny_block,k,bid,TIDAL_DIFF(:,:,k,bid),  &
                           this_block)

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_compute_diff

!***********************************************************************
!BOP
! !IROUTINE: tidal_compute_diff_polzin_2D
! !INTERFACE:

 subroutine tidal_compute_diff_polzin_2D (ZSTARP_INV,N2_AVG,DBLOC,k,bid)

! !DESCRIPTION:
!  Compute tidal diffusion based on Polzin, K. L., 2009; Melet et al 2013
!  2D (horizontal) version
!
! !REVISION HISTORY:
!  njn01 2016

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      ZSTARP_INV,          &! c1/z^*_sub p
      N2_AVG                ! (N**2(z))bar z

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      DBLOC                 ! buoyancy difference between adjacent levels

   integer, intent(in) :: k   ! vertical index
   integer, intent(in) :: bid ! block id

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind)  :: kk
   integer (int_kind)  :: i,j  ! horizontal indices

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1,               &! work array
      WORK2,               &! work array
      ZSTARZ                ! z*(z)

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      WORK3                 ! work array for tavg accumulation

!-----------------------------------------------------------------------
!
!  compute time-dependent dissipation terms (eq 13 Melet 2013)
!
!-----------------------------------------------------------------------

   call tidal_N2_integral(WORK1,DBLOC,k,bid)

   !*** zstarz and (1/H + 1/zstarp)/N2_AVG
   !    note: test N2_AVG = c0 in tidal_N2_integral
   ZSTARZ = WORK1/N2_AVG
   WORK2  = (HTINV(:,:,bid)+ZSTARP_INV(:,:))/N2_AVG(:,:)

   !*** WORK2*1/(1+zstarp(z)/zstarp)**2
   where ( (c1+ZSTARZ*ZSTARP_INV) /= c0)
     WORK1=WORK2/(c1+ZSTARZ*ZSTARP_INV)**2
   elsewhere
     WORK1 = c0  !DEBUG check on this
   endwhere

   !*** N**2, limited by tidal_eps_n2
   WORK2(:,:) = max(DBLOC(:,:,k)/dzw(k),tidal_eps_n2)

   !*** N**2/(N**2+omega**2)
   WORK2 = WORK2/(WORK2 + omega2_polz)
   
!-----------------------------------------------------------------------
!
!  tidal diffusion computed from eq 14 Melet 2013
!
!    note: TIDAL_COEF_2D = 0.2*q*E(x,y)/rho
!
!-----------------------------------------------------------------------

   TIDAL_DIFF(:,:,k,bid) = WORK2(:,:)*TIDAL_COEF_2D(:,:,bid)*WORK1(:,:)

   if (any(TIDAL_DIFF < c0))  then
     write(stdout,*)    '(tidal_compute_diff_polzin_2D) TIDAL_DIFF < c0'
     call shr_sys_abort ('tidal_compute_diff_polzin_2D) TIDAL_DIFF < c0')
   endif

!-----------------------------------------------------------------------
!
!  compute diagnostic quantities for comparison with Melet 2013 Figure 3
!    note: N^2 and TIDAL_DIFF are accumulated in module vmix_kpp.F90
!
!-----------------------------------------------------------------------
   WORK3(:,:,bid) = max(WORK1(:,:),tidal_eps_n2) ! N^2
   WORK3(:,:,bid) = WORK3(:,:,bid)*TIDAL_DIFF(:,:,k,bid)    ! Gamma*epsilon
   call accumulate_tavg_field(WORK3(:,:,bid),tavg_TIDAL_GAMMA_EPS,bid,k)

   WORK1 = c0
   WORK2 = c0
   do kk=2,km
     !*** find sea floor, starting from kk=2
     where (kk == KMT(:,:,bid) )
       WORK1(:,:) = max(DBLOC(:,:,kk-1)/dzw(kk-1),tidal_eps_n2) ! N^2(bottom)
       WORK2 = sqrt(WORK1) ! N(bottom)
     end where
   enddo ! kk
  
  !WORK3(:,:,bid) = p5*WORK2*kappa_polz*H2_P(:,:,bid)*U_P(:,:,bid)
   WORK3(:,:,bid) =    WORK2*kappa_polz*H2_P(:,:,bid)*U_P(:,:,bid)
   call accumulate_tavg_field(WORK3(:,:,bid),tavg_POLZIN_EQ2,bid,k) 

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_compute_diff_polzin_2D

!***********************************************************************
!BOP
! !IROUTINE: tidal_compute_diff_polzin_clmn
! !INTERFACE:

 subroutine tidal_compute_diff_polzin_clmn (ZSTARP_INV,N2_AVG,DBLOC,i,j,nlev,bid)

! !DESCRIPTION:
!  Compute tidal diffusion based on Polzin, K. L., 2009; Melet et al 2013
!  Columnized version
!
! !REVISION HISTORY:
!  njn01 2017

! !INPUT PARAMETERS:
   real (r8), intent(in) :: & 
      ZSTARP_INV,          &! c1/z^*_sub p
      N2_AVG                ! (N**2(z))bar z

   real (r8), dimension(km), intent(in) :: &
      DBLOC                 ! buoyancy difference between adjacent levels

   integer, intent(in) :: i,j  ! horizontal indices
   integer, intent(in) :: nlev ! vertical range index
   integer, intent(in) :: bid  ! block id

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind)  :: k,kk

   real (r8)  ::  &
      WORK1,      &! work array
      WORK2,      &! work array
      ZSTARZ       ! z*(z)

!-----------------------------------------------------------------------
!
!  compute time-dependent dissipation terms (eq 13 Melet 2013)
!
!-----------------------------------------------------------------------

 !*** zero out column
 TIDAL_DIFF(i,j,:,bid) = c0

 !*** then compute over defined k levels
 do k=1,nlev
   call tidal_N2_integral(WORK1,DBLOC,i,j,k,bid)

   !*** zstarz and (1/H + 1/zstarp)/N2_AVG
   !    note: test N2_AVG = c0 in tidal_N2_integral
   ZSTARZ = WORK1/N2_AVG
   WORK2  = (HTINV(i,j,bid)+ZSTARP_INV)/N2_AVG

   !*** WORK2*1/(1+zstarp(z)/zstarp)**2
   if ( (c1+ZSTARZ*ZSTARP_INV) /= c0) then
     WORK1=WORK2/(c1+ZSTARZ*ZSTARP_INV)**2
   else
     WORK1 = c0  !DEBUG check on this
   endif

   !*** N**2, limited by tidal_eps_n2
   WORK2 = max(DBLOC(k)/dzw(k),tidal_eps_n2)

   !*** N**2/(N**2+omega**2)
   WORK2 = WORK2/(WORK2 + omega2_polz)
   
!-----------------------------------------------------------------------
!
!  tidal diffusion computed from eq 14 Melet 2013
!
!    note: TIDAL_COEF_2D = 0.2*q*E(x,y)/rho
!
!-----------------------------------------------------------------------

   TIDAL_DIFF(i,j,k,bid) = WORK2*TIDAL_COEF_2D(i,j,bid)*WORK1

   if (TIDAL_DIFF(i,j,k,bid) < c0)  then
     write(stdout,*)    '(tidal_compute_diff_polzin_clmn) TIDAL_DIFF < c0'
     call shr_sys_abort ('tidal_compute_diff_polzin_clmn) TIDAL_DIFF < c0')
   endif
  
 enddo !k

!-----------------------------------------------------------------------
!  diagnostics -- cannot do here in column mode
!-----------------------------------------------------------------------
!  WORK1 = c0
!  WORK2 = c0
!  do kk=1,km
!    !*** find sea floor
!    where (kk == KMT(i,j,bid) )
!      WORK1 = max(DBLOC(kk-1)/dzw(kk-1),tidal_eps_n2) ! N^2(bottom)
!      WORK2 = sqrt(WORK1) ! N(bottom)
!    end where
!  enddo ! kk

!  WORK3(i,j,bid) =    WORK2*kappa_polz*H2_P(i,j,bid)*U_P(i,j,bid)
!  call accumulate_tavg_field(WORK3(:,:,bid),tavg_POLZIN_EQ2,bid,k)

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_compute_diff_polzin_clmn


!***********************************************************************
!BOP
! !IROUTINE: tidal_min_regions_set_cvmixF
! !INTERFACE:

 subroutine tidal_min_regions_set_cvmixF(i1,i2,j1,j2,k,bid,Tdiff,this_block)

! !DESCRIPTION:
!  Apply tidal_mixing minimum over specified regions, such as the Denmark Strait, 
!   for stability control.
!  This is a POP-specific override to the tidal diffusion computed via cvmix,
!   so it resides in the POP-specific tidal-mixing module. 
!
! !REVISION HISTORY:
!  njn01 February 2017

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) :: i1,i2,j1,j2  ! horizontal indices
   integer (int_kind), intent(in) :: k            ! vertical index
   integer (int_kind), intent(in) :: bid          ! block id

   real (r8), dimension(i1:i2,j1:j2),intent(inout) :: Tdiff ! tidal diffusion

   type (block), intent(in) :: this_block ! block information for current block

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: i,j          ! loop indices
   integer (int_kind) :: num_region   ! loop index
   logical (log_kind) :: first = .true.

   if (ltidal_min_regions) then

       do j=j1,j2
       do i=i1,i2
         do num_region=1,num_tidal_min_regions
           if (REGION_BOX3D(i,j,k,bid) == num_region) then
              Tdiff(i,j) = max(Tdiff(i,j), tidal_min_values(num_region))
              if (nsteps_total == 1) &
              write(stdout,1212) '(tidal_regions_set_min) num_region,k,lat,lon,tidal_diff = ', &
              num_region, this_block%i_glob(i), this_block%j_glob(j),k, &
              TLATD(i,j,bid),TLOND(i,j,bid), Tdiff(i,j)
           endif ! lat lon
         enddo !nn
       end do ! i
       end do ! j

   end if ! ltidal_regions_set_min

1212 format (a, 2i6, 2F10.2, 1pE26.15)

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_min_regions_set_cvmixF

!***********************************************************************

!BOP
! !IROUTINE: tidal_min_regions_set_cvmixT
! !INTERFACE:
 subroutine tidal_min_regions_set_cvmixT(i,j,k,k1,k2,bid,Tdiff,  &
            conversion_factor,this_block)

! !DESCRIPTION:
!  Apply tidal_mixing minimum over specified regions, such as the Denmark Strait, 
!   for stability control.
!  This is a POP-specific override to the tidal diffusion computed via cvmix,
!   so it resides in the POP-specific tidal-mixing module. 
!
! !REVISION HISTORY:
!  njn01 January 2017

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) :: i,j    ! horizontal indices
   integer (int_kind), intent(in) :: k,k1,k2! vertical indices
   integer (int_kind), intent(in) :: bid    ! block id

   real (r8),          intent(in) :: conversion_factor

   type (block),       intent(in) :: this_block ! block information for current block

! !INPUT/OUTPUT PARAMETERS:
   real (r8), dimension(k1:k2),intent(inout) :: Tdiff ! tidal diffusion

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   logical (log_kind) :: first = .true.
   integer (int_kind) :: num_region   ! horizontal index
   real (r8) :: tidal_min_value
   real (r8) :: Tdiff_save


   if (ltidal_min_regions) then
     do num_region=1,num_tidal_min_regions
       if (REGION_BOX3D(i,j,k,bid) == num_region) then
         !  note: tidal_min_value could be pre-computed and stored;
         !  no need to recompute every step, every point
         tidal_min_value = conversion_factor*tidal_min_values(num_region)
         Tdiff(k+1) = max(Tdiff(k+1), tidal_min_value)
!######## DEBUG/TEST -- remove after testing
          if (nsteps_total == 1) then
            write(stdout,1212) '(tidal_min_regions_set_cvmixT) num_region,k,lat,lon,tidal_diff = ', &
            num_region, k, TLATD(i,j,bid),TLOND(i,j,bid), Tdiff(k+1)/conversion_factor
          endif
!######## DEBUG/TEST
       end if ! REGION_BOX3D
     enddo ! num_region
   end if ! ltidal_min_regions

1212 format (1x, a, 2i6, 2F10.2, F25.15)

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_min_regions_set_cvmixT

!***********************************************************************
!BOP
! !IROUTINE: tidal_N2_integral_3D
! !INTERFACE:

 subroutine tidal_N2_integral_3D (N2_INT,DBLOC,ktop,bid)

! !DESCRIPTION:
!  Compute vertical integral of N2 for use in Polzin scheme
!  Full 3D version
!
! !REVISION HISTORY:
!  njn01 2016

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,km), intent(in) :: & 
      DBLOC                !buoyancy difference between adjacent levels

   integer (int_kind), intent(in) ::  ktop  ! index for vertical levels
   integer (int_kind), intent(in) ::  bid   ! block index

! !OUTPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      N2_INT               !vertical integral of N2 from sea bottom to ktop

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: k      ! index for vertical levels

!-----------------------------------------------------------------------
!  compute vertical integral in three parts: bottom layer (half cell),
!    interior (where dzw terms cancel), and the top layer.
!    If computing the full-depth integral, input is ztop = 0, and the
!    formula is modified to account for extrapolating to the surface.
!-----------------------------------------------------------------------
   
    N2_INT(:,:) = c0

   !-----------------------
   ! begin at the sea floor
   !-----------------------
   do k=1,km-1   ! reduce upper bound for intel compiler guidence in debug mode
     where (k == KMT(:,:,bid)-1) 
       N2_INT(:,:) = dzw(k+1)*max(DBLOC(:,:,k)/dzw(k),tidal_eps_n2)
     endwhere
   enddo

   !--------------------------------------------
   ! continue integrating over vertical interior
   !--------------------------------------------
   do k=1,km-1
     where (k <= KMT(:,:,bid)-1 .and. k >= ktop+1)
       N2_INT(:,:) = N2_INT(:,:) + dzw(k)*max(DBLOC(:,:,k)/dzw(k),tidal_eps_n2)
     endwhere
   enddo

   !-----------------------------
   ! end at the top half interval
   !-----------------------------
   if (ktop == 0) then
     N2_INT(:,:) = HTINV(:,:,bid)*(N2_INT(:,:) + dzw(0)*max(DBLOC(:,:,1)/dzw(1),tidal_eps_n2))
   else if (ktop < km) then
     N2_INT(:,:) = N2_INT(:,:) + p5*dz(ktop+1)*max(DBLOC(:,:,ktop)/dzw(ktop),tidal_eps_n2)
   endif

   N2_INT(:,:) = max(RCALCT(:,:,bid)*N2_INT(:,:),tidal_eps_n2)
  
!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_N2_integral_3D

!***********************************************************************
!BOP
! !IROUTINE: tidal_N2_integral_clmn
! !INTERFACE:

 subroutine tidal_N2_integral_clmn (N2_INT,DBLOC,i,j,ktop,bid)

! !DESCRIPTION:
!  Compute vertical integral of N2 for use in Polzin scheme
!  Columnized version
!
! !REVISION HISTORY:
!  njn01 2017

! !INPUT PARAMETERS:
   real (r8), dimension(km), intent(in) :: & 
      DBLOC                !buoyancy difference between adjacent levels

   integer (int_kind), intent(in) :: ktop ! index for vertical levels
   integer (int_kind), intent(in) :: i,j  ! horizontal indices
   integer (int_kind), intent(in) :: bid  ! block id

! !OUTPUT PARAMETERS:
   real (r8), intent(out)         :: N2_INT !vertical integral of N2 from sea bottom to ktop

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: k      ! index for vertical levels

!-----------------------------------------------------------------------
!  compute vertical integral in three parts: bottom layer (half cell),
!    interior (where dzw terms cancel), and the top layer.
!    If computing the full-depth integral, input is ztop = 0, and the
!    formula is modified to account for extrapolating to the surface.
!-----------------------------------------------------------------------
   
    N2_INT = c0

   !-----------------------
   ! begin at the sea floor
   !-----------------------
   do k=1,km-1
     if (k == KMT(i,j,bid)-1)  then
       N2_INT = dzw(k+1)*max(DBLOC(k)/dzw(k),tidal_eps_n2)
     endif
   enddo

   !--------------------------------------------
   ! continue integrating over vertical interior
   !--------------------------------------------
   do k=1,km
     if (k <= KMT(i,j,bid)-1 .and. k >= ktop+1) then
       N2_INT = N2_INT + dzw(k)*max(DBLOC(k)/dzw(k),tidal_eps_n2)
     endif
   enddo

   !-----------------------------
   ! end at the top half interval
   !-----------------------------
   if (ktop == 0) then
     !*** integral is over entire ocean depth; normalize by HT
     N2_INT = HTINV(i,j,bid)*(N2_INT + dzw(0)*max(DBLOC(1)/dzw(1),tidal_eps_n2))
   elseif (ktop < km) then 
     !*** integrate over entire ocean depth
     N2_INT = N2_INT + p5*dz(ktop+1)*max(DBLOC(ktop)/dzw(ktop),tidal_eps_n2)
   endif

   N2_INT = max(RCALCT(i,j,bid)*N2_INT,tidal_eps_n2)
  
!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_N2_integral_clmn

!***********************************************************************
!BOP
! !IROUTINE: tidal_ts_driver
! !INTERFACE:

     subroutine tidal_ts_driver

! !DESCRIPTION: This subroutine modulates the 3D (ER03 or NG13) tidal energy
!               field with the 18.6-year lunar cycle 
   
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len) :: docdate
   character (22      ) :: sub_name = 'tidal_ts_driver'

   integer (int_kind) :: index,ind1,ind2,nn
   integer (int_kind) :: bid, k
   integer (int_kind) :: inx,jny,nlev

   real (r8) ::  cm2,cs2,ck1,co1

!-----------------------------------------------------------------------
!  only recompute qE if using lunar cycle; otherwise, return now
!-----------------------------------------------------------------------
   if (.not. ltidal_lunar_cycle) return

   call timer_start(timer_tidal_ts)

!-----------------------------------------------------------------------
!  compute tidal timeseries factor for today
!    assumption: tidal-constituent indices are 1,2,3,...,ntidal_ts_TC
!    assumption: tidal_ts%data contains 1.0 + original data in percent/100
!-----------------------------------------------------------------------
   ind1 = tidal_ts%data_ind_low
   ind2 = tidal_ts%data_ind_high
   do nn=1,ntidal_ts_TC
   tidal_ts%factor(nn) = tidal_ts%w1(nsteps_this_interval)*tidal_ts%data(ind1,nn) +  &
                         tidal_ts%w2(nsteps_this_interval)*tidal_ts%data(ind2,nn)
   enddo

!-----------------------------------------------------------------------
!  apply factors to tidal energy constituents and reform 3D q*E
!-----------------------------------------------------------------------
   cm2 = tidal_ts%factor(tidal_ts%m2)
   cs2 = tidal_ts%factor(tidal_ts%s2)
   ck1 = tidal_ts%factor(tidal_ts%k1)
   co1 = tidal_ts%factor(tidal_ts%o1)

   call tidal_form_qE_3D (cm2,cs2,ck1,co1)
 
   if (ltidal_lunar_cycle_print) then
     if (my_task == master_task) then
       write(stdout,1410) sub_name, iyear, cmonth3, iday, seconds_this_day, tidal_ts%data_ind_low,  &
                          tidal_ts%data_ind_high, tidal_ts%data_day_of_year, &
                          cm2,cs2,ck1,co1
1410 format(3x,a, i4.4,'-',a3,'-',i2.2,' ',1pe13.6,' sec',3x,3i7,1x,4F25.15)
     endif
   endif

!-----------------------------------------------------------------------
!  reform TIDAL_COEF_3D
!-----------------------------------------------------------------------
   select case (tidal_mixing_method_itype)

      case (tidal_mixing_method_jayne)
      !==============================

        !*** reform 2D q*E before reforming TIDAL_COEF_3D
        call tidal_form_qE_2D_from_3D
        if (.not. registry_match('lcvmix')) then
          call tidal_form_coef_jayne
        else
          !reform the coefficients in module vmix_kpp.F90
        endif

        do bid = 1,nblocks_clinic
          call accumulate_tavg_field(TIDAL_QE_2D(:,:,bid),tavg_TIDAL_QE_2D,bid,1)
        enddo

      case (tidal_mixing_method_schmittner)
      !===================================

        if (.not. registry_match('lcvmix')) then
          call tidal_form_coef_schm
        else
          !reform the coefficients in module vmix_kpp.F90
        endif

     case default

        call exit_POP (sigAbort, '(tidal_ts_driver) only valid for jayne or schmittner method')

   end select

!-----------------------------------------------------------------------
!  set index for next timestep
!-----------------------------------------------------------------------
   call tidal_ts_index_advance

!-----------------------------------------------------------------------
!  after advancing indices, store them for restart purposes
!-----------------------------------------------------------------------
   tidal_ts_data_ind_low     = tidal_ts%data_ind_low
   tidal_ts_data_ind_high    = tidal_ts%data_ind_high
   tidal_ts_data_day_of_year = tidal_ts%data_day_of_year

   if (ltidal_lunar_cycle_print) then
     write(stdout,*) sub_name, ' tidal_ts%data_ind_low = ',     tidal_ts%data_ind_low
     write(stdout,*) sub_name, ' tidal_ts%data_ind_high = ',    tidal_ts%data_ind_high
     write(stdout,*) sub_name, ' tidal_ts%data_day_of_year = ', tidal_ts%data_day_of_year
     call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!  accumulate time-averaged TIDAL_QE
!-----------------------------------------------------------------------
   do bid = 1,nblocks_clinic
     do k=1,km
       call accumulate_tavg_field(TIDAL_QE_3D(:,:,k,bid),tavg_TIDAL_QE_3D,bid,k)
     enddo
   enddo

  call timer_stop(timer_tidal_ts)

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_ts_driver

!***********************************************************************
!BOP
! !IROUTINE: tidal_ts_index_advance
! !INTERFACE:

     subroutine tidal_ts_index_advance

! !DESCRIPTION: This subroutine sets the tidal luna-cycle modulation
!               timeseries indices for the NEXT model timestep
   
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  index, iyear

!-----------------------------------------------------------------------
!  set/advance tidal data timeseries counters as needed
!-----------------------------------------------------------------------

   !--------------------------------------------------------------------
   !  first, test on MODEL current-day boundary. This model-date test is valid 
   !  for daily, monthly, or annual data interpolation.
   !--------------------------------------------------------------------
   if (eod) then

     !*** increment DATA counters, regardless of DATA day, for NEXT model timestep
     tidal_ts%data_day_of_year  = tidal_ts%data_day_of_year + 1
     tidal_ts%data_ind_low      = tidal_ts%data_ind_low     + 1
     tidal_ts%data_ind_high     = tidal_ts%data_ind_high    + 1

     index = tidal_ts%data_ind_low 

     !*** test: is DATA Leap Year adjustment needed?
     if (tidal_ts%calendar_mismatch) then

       iyear = tidal_ts%data_year(index)
       if (mod(iyear,4) == 0 .and. (mod(iyear,100) /= 0 .or. mod(iyear,400) == 0 ) ) then
         if (ltidal_lunar_cycle_print) then
           if (tidal_ts%data_day_of_year == 59 .or.  &
               tidal_ts%data_day_of_year == 60)  then
                 write(stdout,*) 'DATA YEAR is LEAP YEAR', ' iyear = ', iyear
           endif
         endif
         if (tidal_ts%data_day_of_year == 59) then 
           !*** increment upper data index by one more (skip Feb 29 data point)
           tidal_ts%data_ind_high = tidal_ts%data_ind_high + 1
         else if (tidal_ts%data_day_of_year == 60) then 
           !*** increment lower data index by one more (skip Feb 29 data point)
           tidal_ts%data_ind_low = tidal_ts%data_ind_low + 1
           tidal_ts%data_day_of_year  = tidal_ts%data_day_of_year + 1
         endif
       endif
     endif ! calendar_mismatch

     !*** test: is it time to reset DATA day_of_year counter?
     if (eoy) then
         tidal_ts%data_day_of_year  = 1
     endif

     !*** increment DATA counters, regardless of DATA day, for NEXT model timestep
     !*** test: is it time to cycle DATA dates?
     if (index == tidal_ts%data_ind_final) then
       if (ltidal_lunar_cycle_print) then
          write(stdout,*) '(tidal_ts_index_advance) CYCLE 1'
       endif
       tidal_ts%data_ind_high     = tidal_ts%data_ind_first
     else if (index  > tidal_ts%data_ind_final) then
       if (ltidal_lunar_cycle_print) then
          write(stdout,*) '(tidal_ts_index_advance) CYCLE 2'
       endif
       tidal_ts%data_ind_low      = tidal_ts%data_ind_first
       tidal_ts%data_ind_high     = tidal_ts%data_ind_first + 1
       tidal_ts%data_day_of_year  = tidal_ts%data_day_of_year_first
     endif

   else
     !*** do not increment counters
   endif  !eoy/eod/non

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_ts_index_advance


!***********************************************************************
!BOP
! !IROUTINE: tidal_zstarp_inv_3D
! !INTERFACE:

 subroutine tidal_zstarp_inv_3D(ZSTARP_INV,N2_AVG,DBLOC,bid)

! !DESCRIPTION:
!  Compute  1/zstarp(z) for use in Polzin scheme
!  Full 3D version
!
! !REVISION HISTORY:
!  njn01 2016

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,km), intent(in) :: & 
      DBLOC                !buoyancy difference between adjacent levels

   real (r8), dimension(nx_block,ny_block   ), intent(in) :: & 
      N2_AVG               !average N**(z) over full ocn depth

   integer (int_kind), intent(in) :: bid  ! block index

! !OUTPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block), intent(out) :: & 
      ZSTARP_INV           !c1/zstar_p(z)

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  k   ! vertical grid indices

   real (r8), dimension(nx_block,ny_block) ::  &
      NB,                &! N at sea floor
      WORK                ! work array

!-----------------------------------------------------------------------
!
!  compute NB = N at the sea floor
!
!-----------------------------------------------------------------------

   NB = c0
   do k=2,km
     !*** find sea floor (ok to start from k=2 for intel compiler bounds checking work-around
     where (k == KMT(:,:,bid) )
       WORK(:,:) = max(DBLOC(:,:,k-1)/dzw(k-1),tidal_eps_n2) ! N^2(bottom)
       NB = sqrt(WORK) ! N(bottom)
     end where
   enddo ! k
  
!-----------------------------------------------------------------------
!
!  compute 1/zstarp
!
!-----------------------------------------------------------------------

   where (CALCT(:,:,bid) .and. U_P(:,:,bid) /= c0)
     ZSTARP_INV(:,:) = zstar_inv_coeff_polz*H2_P(:,:,bid)*  &
                       NB(:,:)*N2_AVG(:,:)/U_P(:,:,bid)
   elsewhere
     ZSTARP_INV(:,:) = c0
   endwhere

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_zstarp_inv_3D

!***********************************************************************
!BOP
! !IROUTINE: tidal_zstarp_inv_clmn
! !INTERFACE:

 subroutine tidal_zstarp_inv_clmn(ZSTARP_INV,N2_AVG,DBLOC,i,j,bid)

! !DESCRIPTION:
!  Compute  1/zstarp(z) for use in Polzin scheme
!  Columnized version
!
! !REVISION HISTORY:
!  njn01 2017

! !INPUT PARAMETERS:
   real (r8), dimension(km), intent(in) :: &
      DBLOC                !buoyancy difference between adjacent levels
   real (r8), intent(in) :: N2_AVG        ! average N**(z) over full ocn depth

   integer (int_kind), intent(in) :: i,j  ! horizontal indices
   integer (int_kind), intent(in) :: bid  ! block id


! !OUTPUT PARAMETERS:
   real (r8), intent(out) :: ZSTARP_INV   ! c1/zstar_p(z)

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  k  ! index for vertical levels

   real (r8) ::  &
      NB,                &! N at sea floor
      WORK                ! work array

!-----------------------------------------------------------------------
!
!  compute NB = N at the sea floor
!
!-----------------------------------------------------------------------

   NB = c0
   do k=1,km
     !*** find sea floor
     if (k == KMT(i,j,bid) ) then
       WORK = max(DBLOC(k-1)/dzw(k-1),tidal_eps_n2) ! N^2(bottom)
       NB = sqrt(WORK) ! N(bottom)
     endif
   enddo ! k
  
!-----------------------------------------------------------------------
!
!  compute 1/zstarp
!
!-----------------------------------------------------------------------

   if (CALCT(i,j,bid) .and. U_P(i,j,bid) /= c0) then
     ZSTARP_INV = zstar_inv_coeff_polz*H2_P(i,j,bid)*NB*N2_AVG/U_P(i,j,bid)
   else
     ZSTARP_INV = c0
   endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_zstarp_inv_clmn

!***********************************************************************
!BOP
! !IROUTINE: tidal_accumulate_tavg
! !INTERFACE:

 subroutine tidal_accumulate_tavg (WORK,iblock,k,field)

! !DESCRIPTION:
!  Call accumulate_tavg for time-dependent tidal tavg fields
!
! !REVISION HISTORY:
!  njn01 December 2016

! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block), intent(in) :: & 
      WORK                !input field at k,curtime

   integer, intent(in) :: iblock      ! block index
   integer, intent(in) :: k           ! k     index

   character(*), intent(in) :: field  ! field name

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r8), dimension(nx_block,ny_block) ::  WORK1,WORK2  !temporary work arrays

   if (.not. ltidal_mixing) return

  select case (trim(field))

    case ('temp','TEMP','Temp')
!-----------------------------------------------------------------------
!
!  accumulate potential temperatures over pre-defined thicknesses
!    used in creating Melet et al 2013 figure
!    Note: build-namelist scripts activate these only when using gx1v6
!
!-----------------------------------------------------------------------

      if (k >= 40 .and. k <= 44) call accumulate_tavg_field(WORK(:,:),tavg_TEMP1_1p5km,iblock,k)
      if (k >= 40 .and. k <= 46) call accumulate_tavg_field(WORK(:,:),tavg_TEMP1_2km  ,iblock,k)
      if (k >= 40 .and. k <= 51) call accumulate_tavg_field(WORK(:,:),tavg_TEMP1_3km  ,iblock,k)
      if (k >= 46 .and. k <= 51) call accumulate_tavg_field(WORK(:,:),tavg_TEMP2_3km  ,iblock,k)
      if (k >  46 .and. k <= 54) call accumulate_tavg_field(WORK(:,:),tavg_TEMP2_4km  ,iblock,k)
      if (k >  53 .and. k <= km) call accumulate_tavg_field(WORK(:,:),tavg_TEMP3p5_6km,iblock,k)
      if (k >  54 .and. k <= km) call accumulate_tavg_field(WORK(:,:),tavg_TEMP4_6km  ,iblock,k)

    case ('N2')
      call accumulate_tavg_field(WORK(:,:),tavg_TIDAL_N2,iblock,k)

      WORK1(:,:)=max(WORK(:,:),tidal_eps_n2)
      call accumulate_tavg_field(WORK1(:,:),tavg_TIDAL_N2_eps,iblock,k)

      WORK2(:,:) = WORK1(:,:)*TIDAL_DIFF(:,:,k,iblock)    ! Gamma*epsilon
      call accumulate_tavg_field(WORK2(:,:),tavg_TIDAL_GAMMA_EPS,iblock,k) 

    case ('TIDAL_DIFF')
      call accumulate_tavg_field(WORK(:,:),tavg_TIDAL_DIFF,iblock,k)
       
    case ('TIDAL_COEF_3D')
      call accumulate_tavg_field(WORK(:,:),tavg_TIDAL_COEF_3D,iblock,k)
       
    case default
      !*** ignore and return

  end select

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_accumulate_tavg

!***********************************************************************
!BOP
! !IROUTINE: tidal_accumulate_tavg_once
! !INTERFACE:

 subroutine tidal_accumulate_tavg_once

! !DESCRIPTION:
!  Call accumulate_tavg for one-time fields
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r8),dimension(nx_block,ny_block) :: WORK        
   integer (int_kind) :: iblock  ! block index
   integer (int_kind) :: k       ! vertical index

!-----------------------------------------------------------------------
!
!  call accumulate_tavg_field for time-invariant tidal mixing fields,
!   based upon which tidal mixing option is active
!
!------------------------------------------------------------------------

   if (.not. ltidal_mixing) return

   !----------------------------------
   ! mixing-method-specific fields
   !----------------------------------
   select case (tidal_mixing_method_itype)

    case (tidal_mixing_method_polzin)
      do iblock = 1,nblocks_clinic
        call accumulate_tavg_field (H2_P(:,:,iblock),tavg_H2_P,iblock,1)
        call accumulate_tavg_field ( U_P(:,:,iblock),tavg_U_P, iblock,1)
      enddo

    case (tidal_mixing_method_jayne,tidal_mixing_method_schmittner)
      if (.not. ltidal_lunar_cycle) then
!######### DEBUG
      !*** for now, call this from vmix_kpp... will need to figure out "once" vs  time-dependent later
!     do iblock = 1,nblocks_clinic
!       do k=1,km
!       call accumulate_tavg_field(TIDAL_COEF_3D(:,:,k,iblock),tavg_TIDAL_COEF_3D,iblock,k)
!       enddo 
!     enddo 
!######### DEBUG
      endif

      !*** DO NOT deallocate; these fields are needed for computations
   end select

   !----------------------------------
   ! more mixing-method-specific fields
   !----------------------------------
   select case (tidal_mixing_method_itype)
    case (tidal_mixing_method_jayne)
      do iblock = 1,nblocks_clinic
       do k=1,km
        call accumulate_tavg_field(VERTICAL_FUNC(:,:,k,iblock),tavg_VERTICAL_FUNC,iblock,k)
       enddo 
      enddo 
   end select

   !----------------------------------
   !    energy-file-specific fields
   !----------------------------------

   select case (tidal_energy_itype)

    case (tidal_energy_arbic,tidal_energy_jayne)

     do iblock = 1,nblocks_clinic
      call accumulate_tavg_field(TIDAL_ENERGY_FLUX_2D(:,:,iblock),tavg_TIDAL_ENERGY_FLUX_2D,iblock,1)
     enddo ! iblock

    case (tidal_energy_ER03,tidal_energy_GN13,  &
          tidal_energy_LGM0,tidal_energy_LGMi5g21,tidal_energy_LGMi6g21)

     do iblock = 1,nblocks_clinic
       do k=1,km
        call accumulate_tavg_field(TCM2(:,:,k,iblock),tavg_TCM2,iblock,k)
        call accumulate_tavg_field(TCS2(:,:,k,iblock),tavg_TCS2,iblock,k)
        call accumulate_tavg_field(TCK1(:,:,k,iblock),tavg_TCK1,iblock,k)
        call accumulate_tavg_field(TCO1(:,:,k,iblock),tavg_TCO1,iblock,k)
       enddo ! k
       if (.not. ltidal_lunar_cycle) then
         do k=1,km
          call accumulate_tavg_field(TIDAL_QE_3D (:,:,k,iblock),tavg_TIDAL_QE_3D, iblock,k)
         enddo ! k
       endif
      enddo ! iblock

      !*** deallocate after writing to one-time tavg output file
      if (.not. ltidal_lunar_cycle) deallocate (TCM2, TCS2, TCK1, TCO1)

   end select

   !----------------------------------
   !    all energy files
   !----------------------------------
   if (.not. ltidal_lunar_cycle) then
     do iblock = 1,nblocks_clinic
       call accumulate_tavg_field(TIDAL_QE_2D(:,:,iblock),tavg_TIDAL_QE_2D,iblock,1)
     enddo ! iblock
   endif

   !----------------------------------
   !    tidal-mixing floor regions
   !----------------------------------
   if (ltidal_min_regions) then
   do iblock = 1,nblocks_clinic
   do k=1,km
     WORK(:,:) = REGION_BOX3D(:,:,k,iblock)
     call accumulate_tavg_field(WORK,tavg_REGION_BOX3D,iblock,k)
   enddo ! k
   enddo ! iblock
   endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_accumulate_tavg_once

!***********************************************************************
!BOP
! !IROUTINE: tidal_compare_VIC
! !INTERFACE:

 subroutine tidal_compare_VIC (WORK,method)

! !DESCRIPTION: Compare POP to UVic 19-level model

! !INPUT PARAMETERS:
   real (r8),dimension(nx_block,ny_block,nblocks_clinic),intent(inout) :: WORK        
   character (*), intent(in) :: method

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r8) :: factor
   real (r8) :: sum_horiz_TIDAL_QE_3D

   integer (int_kind) :: k  ! vertical index
   integer (int_kind) :: kk ! vertical index

   WORK = c0
   sum_horiz_TIDAL_QE_3D = c0
   kk = 1
   sum_vic_TE = c0 

   if (trim(method) == 'jayne') then
     factor = c1/tidal_gamma_rhor
   else if (trim(method) == 'schmittner') then
     factor = c1
   endif

   do k=1,km
      !*** TIDAL_COEF_3D = (gamma/rho)*q*E*F(z1,z2)
      !*** TIDAL_QE_nD = q*E 2D or 3D
      if (trim(method) == 'jayne') then
        WORK(:,:,:) = factor*TIDAL_COEF_3D(:,:,k,:)*RCALCT(:,:,1:nblocks_clinic)*TAREA(:,:,1:nblocks_clinic)
      else if (trim(method) == 'schmittner') then
         WORK(:,:,:) = factor*TIDAL_QE_3D(:,:,k,:)*RCALCT(:,:,1:nblocks_clinic)*TAREA(:,:,1:nblocks_clinic)
      endif

      horiz_TIDAL_QE_3D(k) = global_sum(WORK,distrb_clinic,field_loc_center)*1.0e-19
     !                 10^-3 W * 10^-4 m^2 --> cm^2 * 10^-12 terra   ==> TW *1.0e-19
      write(stdout,*) 'global_sum(TIDAL_QE_3D(',k,') TW = ', horiz_TIDAL_QE_3D(k)

      sum_horiz_TIDAL_QE_3D = sum_horiz_TIDAL_QE_3D + horiz_TIDAL_QE_3D(k)
      if (kk .lt. 19) then
        sum_vic_TE(kk) = sum_vic_TE(kk) + horiz_TIDAL_QE_3D(k)
        if (zt(k) .le. zw_vic(kk)) then
          ! continue accumulating 
        else
          kk = kk+1
        endif
      endif !kk 
   enddo ! k

   write(stdout,*) 'total global_sum TW = ',  sum_horiz_TIDAL_QE_3D

   do kk=1,18
     write(stdout,*) 'sum over dz_vic TIDAL_QE_3D(',kk,') TW = ', sum_vic_TE(kk)
   enddo

!-----------------------------------------------------------------------
!EOC
 
 end subroutine tidal_compare_VIC

!***********************************************************************

end module tidal_mixing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

