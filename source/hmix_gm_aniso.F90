!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      module hmix_gm_aniso

!BOP
! !MODULE: hmix_gm_aniso

! !DESCRIPTION:
!  This module contains routines for computing horizontal mixing
!  using the Gent-McWilliams eddy transport parameterization
!  and isopycnal diffusion.

! !REVISION HISTORY:
!  SVN:$Id: hmix_gm.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

      use kinds_mod
      use blocks
      use distribution
      use domain
      use constants
      use broadcast
      use grid
      use io
      use vertical_mix
      use vmix_kpp
      use state_mod
      use time_management
      use tavg
      use diagnostics
      use exit_mod
      use registry
      use hmix_gm_submeso_share
      use operators, only: zcurl, grad

#ifdef CCSMCOUPLED
   use shr_sys_mod
#endif
      use timers, only: timer_start, timer_stop, get_timer !, timer_print


      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: init_gm_aniso,   &
                hdifft_gm_aniso

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     variables to save from one call to next
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ktp = 1, kbt = 2     ! refer to the top and bottom halves of a
                              !  grid cell, respectively

      real (r8), dimension(:), allocatable :: &
         kappa_depth          ! depth dependence for KAPPA 

      real (r8), dimension(:,:,:), allocatable :: &
         HYXW, HXYS, &        ! west and south-shifted values of above
         RB,         &        ! Rossby radius
         RBR,        &        ! inverse of Rossby radius
         BTP,        &        ! beta plane approximation
         BL_DEPTH,   &        ! boundary layer depth
         UIT, VIT,   &        ! work arrays for isopycnal mixing velocities
         WTOP_ISOP, WBOT_ISOP ! vertical component of isopycnal velocities

      real (r8), dimension(:,:,:,:,:,:), allocatable :: &
         SF_SLXX, SF_SLYY,  & !SJR MOD ! components of the merged streamfunction
         SF_SLXY, SF_SLYX,  & !SJR ADD
         DD_SLX,  DD_SLY      !SJR ADD ! isopycnal slopes with grid spacings 

      real (r8), dimension(:,:,:,:,:,:,:), allocatable :: & !SJR ADD
         SUBCELLV, &  !SJR ADD ! eight subcell volume
         SL_MAG       !SJR ADD ! isopycnal slope magnitude

      real (r8), dimension(:,:,:,:,:), allocatable :: &
         SLA_SAVE             ! isopycnal slopes

      real (r8), dimension(:,:,:,:), allocatable :: &
         DZTW,   &  !SJR ADD
         FZTOP                ! vertical flux

      logical (log_kind), dimension(:), allocatable :: &
         compute_kappa        ! compute spatially varying coefficients
                              !  this time step?

      logical (log_kind) ::   &
         diff_tapering,       &   ! different tapering for two diffusivities
         cancellation_occurs, &   ! specified choices for the isopycnal and
                                  !  thickness diffusion coefficients result in 
                                  !  cancellation of some tensor elements
         read_n2_data             ! if .true., use climatological N^2 data 
                                  !  instead of model dependent N^2

!-----------------------------------------------------------------------
!
!     KAPPA_LATERAL and KAPPA_VERTICAL are 2D and 3D arrays, respectively,
!     containing the spatial variations of the isopycnal (KAPPA_ISOP)
!     and thickness (KAPPA_THIC) diffusion coefficients. Except in kappa_type_eg,
!     KAPPA_LATERAL has the actual diffusion coefficient values in cm^2/s,
!     whereas KAPPA_VERTICAL is unitless. So, the total coefficients are
!     constructed using
!
!      KAPPA_ISOP or KAPPA_THIC (:,:,:,k,bid) ~ KAPPA_LATERAL(:,:,bid)
!                                             * KAPPA_VERTICAL(:,:,k,bid)
!
!     When kappa_type_eg, KAPPA_VERTICAL contains the isopycnal diffusivity
!     coefficients in cm^2/s and KAPPA_LATERAL is not used!
!
!-----------------------------------------------------------------------

      real (r8), dimension(:,:,:,:,:), allocatable :: &
         KAPPA_ISOP, &      ! 3D isopycnal diffusion coefficient
                            !  for top and bottom half of a grid cell
         KAPPA_THIC, &      ! 3D thickness diffusion coefficient
                            !  for top and bottom half of a grid cell
         KXX_ISOP, &      !SJR ADD
         KXY_ISOP, &      !SJR ADD
         KYY_ISOP, &      !SJR ADD
         KXX_THIC, &      !SJR ADD
         KXY_THIC, &      !SJR ADD
         KYY_THIC, &      !SJR ADD
         HXX_DIFF, &      !SJR ADD
         HXY_DIFF, &      !SJR ADD
         HYY_DIFF, &      !SJR ADD
         K_EIGENVAL_RAT, &  !SJR ADD
         KRAT_SHRD, &  !SJR ADD
         HOR_DIFF           ! 3D horizontal diffusion coefficient
                            !  for top and bottom half of a grid cell

      real (r8), dimension(:,:,:,:), allocatable :: & !SJR ADD
         KXX_IN, &  !SJR ADD !KXX read in from file
         KXY_IN, &  !SJR ADD
         KYY_IN, &  !SJR ADD
         MAJOR_EIGENVAL, &  !SJR ADD
         MINOR_EIGENVAL, &  !SJR ADD
         NX_ISOP, NY_ISOP, &        !SJR ADD  !principal axis of major eigenvalue
         NX_SHRD, NY_SHRD        !SJR ADD  !principal axis of major eigenvalue


      real (r8), dimension(:,:,:), allocatable :: &
         KAPPA_LATERAL      ! horizontal variation of KAPPA in cm^2/s

      real (r8), dimension(:,:,:,:), allocatable :: &
         KAPPA_VERTICAL     ! vertical variation of KAPPA (unitless),
                            !  e.g. normalized buoyancy frequency dependent 
                            !  profiles at the tracer grid points
                            !  ( = N^2 / N_ref^2 ) OR a time-independent
                            !  user-specified function

      real (r8), dimension(:,:,:,:), allocatable :: &
         BUOY_FREQ_SQ,    & ! N^2 defined at level interfaces
         SIGMA_TOPO_MASK    ! bottom topography mask used with kappa_type_eg

!-----------------------------------------------------------------------
!
!     GM specific options
!
!     kappa_freq = how often spatial variations of the diffusion 
!                  coefficients are computed. Same frequency is 
!                  used for both coefficients.
!     slope_control = tanh function (Danabasoglu and McWilliams 1995) or
!                     DM95 with replacement function to tanh or
!                     slope clipping or
!                     method of Gerdes et al (1991)
!     diag_gm_bolus = .true. for diagnostic bolus velocity computation.
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter ::   &
         kdir_type_shrd           = 0,   & !SJR ADD !using shear dispersion parameterization
         kdir_type_east           = 1,   & !SJR ADD !using grad(f)
         kdir_type_zonl           = 2,   & !SJR ADD !using angle
         kdir_type_flow           = 3,   & !SJR ADD !align with flow velcoity
         kdir_type_apvg           = 4,   & !SJR ADD !across pv gradient
         kdir_type_read           = 5,   & !SJR ADD !read from files for KXX, KXY, KYY
         kmin_type_simp           = 1,   & !SJR ADD
         kmin_type_read           = 2,   & !SJR ADD
         krat_type_shrd           = 0,   & !SJR ADD
         krat_type_simp           = 1,   & !SJR ADD
         krat_type_read           = 2,   & !SJR ADD
         kappa_type_const         = 1,   &
         kappa_type_depth         = 2,   &
         kappa_type_vmhs          = 3,   &
         kappa_type_hdgr          = 4,   &
         kappa_type_dradius       = 5,   &
         kappa_type_bfreq         = 6,   &
         kappa_type_bfreq_vmhs    = 7,   &
         kappa_type_bfreq_hdgr    = 8,   &
         kappa_type_bfreq_dradius = 9,   &
         kappa_type_eg            = 10,  &
         slope_control_tanh   = 1,       &
         slope_control_notanh = 2,       &
         slope_control_clip   = 3,       &
         slope_control_Gerd   = 4,       &
         kappa_freq_never           = 1, &
         kappa_freq_every_time_step = 2, &
         kappa_freq_once_a_day      = 3  

      integer (int_kind) :: &
         kappa_isop_type,   &   ! choice of KAPPA_ISOP
         kappa_thic_type,   &   ! choice of KAPPA_THIC
         kappa_freq,        &   ! frequency of KAPPA computations
         slope_control          ! choice for slope control

      logical (log_kind) ::  &
         diag_gm_bolus          ! true for diagnostic bolus velocity computation 

!-----------------------------------------------------------------------
!
!     if use_const_ah_bkg_srfbl = .true., then the specified constant
!     value of ah_bkg_srfbl is used as the "maximum" background horizontal 
!     diffusivity within the surface boundary layer. Otherwise,
!     KAPPA_ISOP is utilized as this "maximum".
!
!-----------------------------------------------------------------------
  
      logical (log_kind) ::      &
         use_const_ah_bkg_srfbl, & ! see above 
         transition_layer_on       ! control for transition layer parameterization
                               
      real (r8) ::      &
         ah,            &       ! isopycnal diffusivity
         ah_bolus,      &       ! thickness (GM bolus) diffusivity
         ah_bkg_bottom, &       ! backgroud horizontal diffusivity at k = KMT
         ah_bkg_srfbl,  &       ! backgroud horizontal diffusivity within the
                                !  surface boundary layer
         slm_r,         &       ! max. slope allowed for redi diffusion
         slm_b                  ! max. slope allowed for bolus transport

!-----------------------------------------------------------------------
!
!     the following set of variables are used in Eden and Greatbatch 
!     (2008) KAPPA formulation. They are in the input namelist.
!
!-----------------------------------------------------------------------

      real (r8) ::      &
         const_eg,      &  ! tuning parameter (unitless)
         gamma_eg,      &  ! (> 0) effective upper limit for inverse eddy 
                           !  time scale (unitless)
         kappa_min_eg,  &  ! minimum KAPPA (cm^2/s)
         kappa_max_eg      ! maximum KAPPA (cm^2/s)

!-----------------------------------------------------------------------
!
!     transition layer type variables
!
!----------------------------------------------------------------------- 

      type tlt_info
        real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
           DIABATIC_DEPTH,  &   ! depth of the diabatic region at the
                                !  surface, i.e. mean mixed or boundary layer
                                !  depth
           THICKNESS,       &   ! transition layer thickness
           INTERIOR_DEPTH       ! depth at which the interior, adiabatic
                                !  region starts, i.e.
                                !   = TLT%DIABATIC_DEPTH + TLT%THICKNESS
        integer (int_kind), &
              dimension(nx_block,ny_block,max_blocks_clinic) :: &
           K_LEVEL,  &          ! k level at or below which the interior,
                                !  adiabatic region starts
           ZTW                  ! designates if the interior region
                                !  starts below depth zt or zw.
                                !  ( = 1 for zt, = 2 for zw )
      end type tlt_info

      type (tlt_info) :: &
         TLT                    ! transition layer thickness related fields 

!-----------------------------------------------------------------------
!
!     tavg ids for tavg diagnostics related to diffusivities and
!     isopycnal velocities. Zonal and meridional refer here to logical 
!     space only.
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         tavg_UISOP,        &   ! zonal      isopycnal velocity
         tavg_VISOP,        &   ! meridional isopycnal velocity
         tavg_WISOP,        &   ! vertical   isopycnal velocity
         tavg_KAPPA_ISOP,   &   ! isopycnal  diffusion coefficient 
         tavg_KXX_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_KXY_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_KYY_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_MIN_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_MAJ_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_C2T_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_S2T_ISOP,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_KAPPA_THIC,   &   ! thickness  diffusion coefficient
         tavg_KXX_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_KXY_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_KYY_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_MIN_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_MAJ_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_C2T_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_S2T_THIC,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_HOR_DIFF,     &   ! horizontal diffusion coefficient
         tavg_HXX_DIFF,   &   ! horizontal diffusion coefficient !SJR ADD
         tavg_HXY_DIFF,   &   ! horizontal diffusion coefficient !SJR ADD
         tavg_HYY_DIFF,   &   ! horizontal diffusion coefficient !SJR ADD
         tavg_MIN_DIFF,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_MAJ_DIFF,   &   ! isopycnal  diffusion coefficient  !SJR ADD
         tavg_C2T_DIFF,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_S2T_DIFF,   &   ! isopycnal  diffusion coefficient  !SJR ADD 
         tavg_DIA_DEPTH,    &   ! depth of the diabatic region at the surface
         tavg_TLT,          &   ! transition layer thickness
         tavg_INT_DEPTH,    &   ! depth at which the interior region starts
         tavg_ADVT_ISOP,    &   ! vertically-integrated T eddy-induced
                                !  advection tendency
         tavg_ADVS_ISOP,    &   ! vertically-integrated S eddy-induced
                                !  advection tendency
         tavg_VNT_ISOP,     &   ! heat flux tendency in grid-y direction
                                !  due to eddy-induced velocity
         tavg_VNS_ISOP          ! salt flux tendency in grid-y direction
                                !  due to eddy-induced velocity

!-----------------------------------------------------------------------
!
!     input namelist for setting GM ANISO options
!
!-----------------------------------------------------------------------
      logical (log_kind) ::   &
         addrandfluc,  & !SJR ADD ! T to add random fluctuation to orientation
         cflmajoronly, & !SJR ADD ! T to reduce major only for cfl violations, F to reduce entire tensor
         vertdiffhere, & !SJR ADD ! do VDC here, F to do in vertical_mix 
         simpsubcells, & !SJR ADD ! T to use simple subcell volume = 1/8 T-cell volume, F to use HTN & HTE
         isoonly,      & !SJR ADD ! T to do isotropic using diagnosis
         isominoronly, & !SJR ADD ! T to set isotropic diffusivity with minor, F to use avg of major and minor 
         savenewtavgs    !SJR ADD ! T to set isotropic diffusivity with minor, F to use avg of major and minor 

      real (r8) ::      &
         cflmult,       &  !SJR ADD ! multiplication factor for cfl check
         erat_const,    &  !SJR ADD ! constant eigenvalue ratio 
         minorfactor,   &  !SJR ADD ! minor eigenvalue multiplicative factor 
         erat_factor,   &  !SJR ADD ! max negative factor of major to set minor, set to 0 to force minor to be positive 
         shrdispfac        !SJR ADD ! multiplicative coefficient for shear dispersion term: MAJOR = MINOR + shrdispfac/MINOR*<(U*dy)^2+(V*dx)^2> 

      integer (int_kind) :: &
         kdir_type,         &   !SJR ADD - choice for major axis direction
         kmin_type,         &   !SJR ADD - choice for minor diffusivity
         krat_type              !SJR ADD - choice for diffusivity ratio (major/minor)

      character (char_len) :: &
         kdir_type_choice,    & 
         kmin_type_choice,    & 
         krat_type_choice

!-----------------------------------------------------------------------
!
!  timers
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
      timer_gm0, &         ! main n loop
      timer_gm1, &         ! main n loop
      timer_gm2, &         ! main n loop
      timer_gm3, &         ! main n loop
      timer_gm4, &         ! main n loop
      timer_gm5, &         ! main n loop
      timer_gm6, &         ! main n loop
      timer_gm7, &         ! main n loop
      timer_gm7a, &         ! main n loop
      timer_gm7b, &         ! main n loop
      timer_gm8, &         ! main n loop
      timer_gm9, &         ! main n loop
      timer_nloop         ! main n loop

!-----------------------------------------------------------------------
!
!     GM specific options
!
!     kappa_freq = how often spatial variations of the diffusion 
!                  coefficients are computed. Same frequency is 
!                  used for both coefficients.
!     slope_control = tanh function (Danabasoglu and McWilliams 1995) or
!                     DM95 with replacement function to tanh or
!                     slope clipping or
!                     method of Gerdes et al (1991)
!     diag_gm_bolus = .true. for diagnostic bolus velocity computation.
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

      contains

!***********************************************************************
!BOP
! !IROUTINE: init_gm_aniso
! !INTERFACE:

      subroutine init_gm_aniso

! !DESCRIPTION:
!  Initializes various choices and allocates necessary space.
!  Also computes some time-independent arrays.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   type (block) :: &
      this_block         ! block information for current block

   character (char_len) ::    &
      KXX_filename, &! file name for the time-independent KXX !SJR ADD 
      KXY_filename, &! file name for the time-independent KXY !SJR ADD 
      KYY_filename, &! file name for the time-independent KYY !SJR ADD 
      buoyancy_freq_filename, &! file name for the time-independent
                               !  buoyancy frequency squared
      buoyancy_freq_fmt,      &! format (bin or netcdf) of buoyancy file
      message                  ! string to hold error message

   type (datafile) :: &
      KXX_data_file, &  ! data file descriptor for KXX data !SJR ADD
      KXY_data_file, &  ! data file descriptor for KXY data !SJR ADD
      KYY_data_file, &  ! data file descriptor for KYY data !SJR ADD
      buoyancy_freq_data_file  ! data file descriptor for buoyancy freq data

   type (io_field_desc) :: &
      KXX_data_in, &  ! io field descriptor for input KXX data !SJR ADD
      KXY_data_in, &  ! io field descriptor for input KXY data !SJR ADD
      KYY_data_in, &  ! io field descriptor for input KYY data !SJR ADD
      buoyancy_freq_data_in  ! io field descriptor for input buoyancy freq data

   type (io_dim) :: &
      i_dim, j_dim, &  ! dimension descriptors for horiz dims
      k_dim            ! dimension descriptor  for depth

!-----------------------------------------------------------------------
!
!     input namelist for setting various GM options
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error,         &  ! error flag for namelist
      k,                 &  ! vertical level index
      iblock,            &  ! block index
      istrt, jstrt,      &  ! SJR ADD ! indices for subcell volumes   
      istop, jstop,      &  ! SJR ADD ! indices for subcell volumes   
      i_sub, j_sub,      &  ! SJR ADD ! indices for subcell volumes   
      i, j                  ! lateral indices

   real (r8) ::          &
      kappa_depth_1,     &  ! parameters for variation of KAPPA 
      kappa_depth_2,     &  !  with depth used with the kappa_type_depth
      kappa_depth_scale     !  option

   character (char_len) :: &
      kappa_choice,        &  ! (generic) choice for KAPPA
      kappa_freq_choice,   &  ! frequency choice for KAPPA computation
      kappa_isop_choice,   &  ! choice for KAPPA_ISOP
      kappa_thic_choice,   &  ! choice for KAPPA_THIC
      slope_control_choice    ! choice for slope control


      real (r8), dimension(nx_block,ny_block) :: &
         WORK3            ! local work space

   namelist /hmix_gm_nml/ kappa_isop_choice,                     &
                          kappa_thic_choice,                     &
                          kappa_freq_choice,                     &
                          slope_control_choice,                  &
                          kappa_depth_1, kappa_depth_2,          &
                          kappa_depth_scale, ah, ah_bolus,       &
                          use_const_ah_bkg_srfbl, ah_bkg_srfbl,  &
                          ah_bkg_bottom,                         &
                          slm_r, slm_b, diag_gm_bolus,           &
                          transition_layer_on, read_n2_data,     &
                          buoyancy_freq_filename,                &
                          buoyancy_freq_fmt,                     &
                          const_eg, gamma_eg,                    &
                          kappa_min_eg, kappa_max_eg

!-----------------------------------------------------------------------
!
!     read input namelist for additional GM options
!
!     DEFAULT SETUP: 
!     KAPPA_ISOP      : constant (= ah)
!     KAPPA_THIC      : constant (= ah_bolus)
!     kappa_freq      : never
!     slope control   : method by DM95 with replacing tanh by polynomial
!                        (= slope_control_notanh)
!
!     variation of kappa with depth used with kappa_type_depth is
!       kappa_depth_1 + kappa_depth_2*exp(-z/kappa_depth_scale)
!
!      with kappa_depth_1 = 1, kappa_depth_2 = 0, and 
!      kappa_depth_scale = 150000 cm, i.e. no depth variation
!
!     ah_bolus        : thickness diffusion coefficient (= 0.8e07 cm^2/s)
!     ah_bkg_bottom   : background horizontal diffusion at k=KMT (= 0)
!     use_const_ah_bkg_srfbl : use ah_bkg_srfbl value
!     ah_bkg_srfbl    : background horizontal diffusion within the
!                        surface boundary layer (= 0.8e07 cm^2/s)
!     slm_r           : max. slope allowed for isopycnal (redi) diffusion (= 0.3)
!     slm_b           : max. slope allowed for thickness (bolus) diffusion (= 0.3)
!     diag_gm_bolus   : .true.
!     transition_layer_on: .false.
!     read_n2_data    : .false.
!     const_eg        : 1.
!     gamma_eg        : 300.
!     kappa_min_eg    : 0.35e07 cm^2/s
!     kappa_max_eg    : 5.00e07 cm^2/s
!
!-----------------------------------------------------------------------


   namelist /hmix_gm_aniso_nml/ addrandfluc,                     &
                                cflmajoronly,                    &
                                cflmult,                         &
                                erat_const,                      &
                                erat_factor,                     &
                                isominoronly,                    &
                                isoonly,                         &
                                kdir_type_choice,                &
                                kmin_type_choice,                &
                                krat_type_choice,                &
                                minorfactor,                     &
                                savenewtavgs,                    &
                                shrdispfac,                      &
                                simpsubcells,                    &
                                vertdiffhere

!-----------------------------------------------------------------------
!
!     read input namelist for additional GM ANISO options
!

!-----------------------------------------------------------------------
!
!     register init_gm_aniso
!
!-----------------------------------------------------------------------

   call register_string('init_gm_aniso')
 
   kappa_isop_type = kappa_type_const
   kappa_thic_type = kappa_type_const
   kappa_freq    = kappa_freq_never
   slope_control = slope_control_notanh
   kappa_depth_1 = c1
   kappa_depth_2 = c0
   kappa_depth_scale = 150000.0_r8
   ah            = 0.8e7_r8
   ah_bolus      = 0.8e7_r8
   ah_bkg_bottom = c0
   ah_bkg_srfbl  = 0.8e7_r8
   slm_r = 0.3_r8
   slm_b = 0.3_r8
   diag_gm_bolus          = .true.
   use_const_ah_bkg_srfbl = .true.
   transition_layer_on    = .false.
   read_n2_data           = .false.
   buoyancy_freq_filename = 'unknown-buoyancy'
   buoyancy_freq_fmt      = 'nc'
   const_eg               = c1 
   gamma_eg               = 300.0_r8 
   kappa_min_eg           = 0.35e7_r8
   kappa_max_eg           = 5.0e7_r8

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error =  1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=hmix_gm_nml,iostat=nml_error)
     end do
     if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading hmix_gm_nml')
   endif

!-----------------------------------------------------------------------
!
!     Defaults for aniso namelist
!
!-----------------------------------------------------------------------
   
   addrandfluc = .false.
   cflmajoronly = .true.
   cflmult = 0.175_r8
   erat_const = c5
   erat_factor = c0
   isominoronly = .true.
   isoonly = .false.
   kdir_type = kdir_type_shrd !flow
   kmin_type = kmin_type_simp
   krat_type = krat_type_shrd !simp
   minorfactor = c1 !sqrt(c5)
   savenewtavgs = .true.
   shrdispfac = p5 / pi**2
   simpsubcells = .false.
   vertdiffhere = .false.						 

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error =  1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=hmix_gm_aniso_nml,iostat=nml_error)
     end do
     if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading hmix_gm_aniso_nml')
   endif

!-------------v--SJR ADD--v--------------------------------------------- 
   if (kmin_type == kmin_type_simp) erat_factor=c0
   if (isoonly) then
      erat_const=c1
      kdir_type=kdir_type_flow
      krat_type=krat_type_simp
      savenewtavgs = .false.
   endif
!-------------^--SJR ADD--^--------------------------------------------- 

   if (my_task == master_task) then

     write(stdout,*) ' '
     write(stdout,*) ' Document Namelist Parameters:'
     write(stdout,*) ' ============================ '
     write(stdout,*) ' '
     write(stdout,  hmix_gm_nml)
     write(stdout,*) ' '

     write(stdout,*) '  Gent-McWilliams options:'
     write(stdout,*) '    kappa_isop choice is ',  &
                     trim(kappa_isop_choice)
     write(stdout,*) '    kappa_thic choice is ',  &
                     trim(kappa_thic_choice)
     write(stdout,*) '    kappa computation frequency choice is ',  &
                     trim(kappa_freq_choice)
     write(stdout,*) '    slope control choice is ',  &
                                    trim(slope_control_choice)
     write(stdout,'(a28,1pe13.6)') ' isopycnal diffusion      = ',  &
                                     ah
     write(stdout,'(a28,1pe13.6)') ' thickness diffusion      = ',  &
                                     ah_bolus
     write(stdout,'(a28,1pe13.6)') ' backgroud diff. (bottom) = ',  &
                                     ah_bkg_bottom
     write(stdout,'(a28,1pe13.6)') ' backgroud diff. (srfbl)  = ',  &
                                     ah_bkg_srfbl
     write(stdout,'(a28,1pe13.6)') ' max. slope for redi      = ',  &
                                     slm_r
     write(stdout,'(a28,1pe13.6)') ' max. slope for bolus     = ',  &
                                     slm_b
     if ( diag_gm_bolus ) then
       write(stdout,'(a47)')  &
           '    diagnostic bolus velocity computation is on'
     endif

     if ( use_const_ah_bkg_srfbl ) then
       write(stdout,1001)
     else
       write(stdout,1002)
     endif
 1001   format(/,' Maximum horizontal background diffusion ', &
                 'coefficient', /,                            &
           '  within the boundary layer is constant at ah_bkg_srfbl.')
 1002   format(/,' Maximum horizontal background diffusion ', &
                 'coefficient', /,                            &
           '  within the boundary layer depends on KAPPA_ISOP.')

     if ( transition_layer_on ) then
       write(stdout,'(a33)') '  transition layer scheme is on. '
     endif

     write(stdout,*) ' '
     write(stdout,  hmix_gm_aniso_nml)
     write(stdout,*) ' '

     write(stdout,*) '  Gent-McWilliams Anisotropic options:'
     write(stdout,*) '     major axis direction             ',  &
                     kdir_type_choice
     write(stdout,*) '     minor diffusivity                ',  &
                     kmin_type_choice
     write(stdout,*) '     diffusivity ratio (major/minor)  ',  &
                     krat_type_choice

     if (addrandfluc) &
     write(stdout,*) '     add random fluctuation to orientation '

     if (cflmajoronly) then
        write(stdout,*) '     reduce major only for cfl violations'
     else
        write(stdout,*) '     reduce entire tensor for cfl violations'
     end if

     if (isominoronly) then
        write(stdout,*) '     set isotropic diffusivity with minor'
     else
        write(stdout,*) '     set isotropic diffusivity with avg of major and minor'
     end if

     if (isoonly) then
        write(stdout,*) '     do isotropic using diagnosis'
     else
        write(stdout,*) '     dont do isotropic using diagnosis'
     end if

     if (simpsubcells) then
        write(stdout,*) '     subcell volume = 1/8 T-cell volume (simple)'
     else
        write(stdout,*) '     subcell volume determined using HTN and HTE'
     end if

     if (savenewtavgs) then
        write(stdout,*) '     save aniso time averaged diagnostics'
     end if

     if (vertdiffhere) then
        write(stdout,*) '     calculate effective vertical diffusion coefficient here'
     else
        write(stdout,*) '     calculate effective vertical diffusion coefficient in vertical_mix'
     end if

     write(stdout,'(a28,1pe13.6)') '     mult factor for cfl check is ', &
                                   cflmult
     write(stdout,'(a28,1pe13.6)') '     constant eigenvalue ratio is ', &
                                   erat_const
     write(stdout,'(a28,1pe13.6)') '     max neg factor of major to set minor ', &
                                   erat_factor
     write(stdout,'(a28,1pe13.6)') '     minor eigenvalue mult factor ', &
                                   minorfactor
     write(stdout,'(a28,1pe13.6)') '     mult coef for shear dispersion', &
                                   shrdispfac

     select case (kdir_type_choice(1:4))
     case ('shea')
        kdir_type=kdir_type_shrd
     case ('east')
        kdir_type=kdir_type_east
     case ('zona')
        kdir_type=kdir_type_zonl
     case ('flow')
        kdir_type=kdir_type_flow
     case ('pvgr')
        kdir_type=kdir_type_apvg
     case ('read')
        kdir_type=kdir_type_read
     case default
        kdir_type = -1000
     end select

     select case (kmin_type_choice(1:4))
     case ('simp')
        kmin_type=kmin_type_simp
     case ('read')
        kmin_type=kmin_type_read
     case default
        kmin_type = -1000
     end select

     select case (krat_type_choice(1:4))
     case ('shea')
        krat_type=krat_type_shrd
     case ('simp')
        krat_type=krat_type_simp
     case ('read')
        krat_type=krat_type_read
     case default
        krat_type = -1000
     end select

     select case (kappa_isop_choice(1:4))
     case ('cons')
       kappa_isop_type = kappa_type_const
     case ('dept')
       kappa_isop_type = kappa_type_depth
     case ('vmhs')
       kappa_isop_type = kappa_type_vmhs
     case ('hdgr')
       kappa_isop_type = kappa_type_hdgr
     case ('drad')
       kappa_isop_type = kappa_type_dradius 
     case ('bfre')
       kappa_isop_type = kappa_type_bfreq 
     case ('bfvm')
       kappa_isop_type = kappa_type_bfreq_vmhs
     case ('bfhd')
       kappa_isop_type = kappa_type_bfreq_hdgr
     case ('bfdr')
       kappa_isop_type = kappa_type_bfreq_dradius
     case ('edgr')
       kappa_isop_type = kappa_type_eg
     case default
       kappa_isop_type = -1000
     end select

     select case (kappa_thic_choice(1:4))
     case ('cons')
       kappa_thic_type = kappa_type_const
     case ('dept')
       kappa_thic_type = kappa_type_depth
     case ('vmhs')
       kappa_thic_type = kappa_type_vmhs
     case ('hdgr')
       kappa_thic_type = kappa_type_hdgr
     case ('drad')
       kappa_thic_type = kappa_type_dradius
     case ('bfre')
       kappa_thic_type = kappa_type_bfreq
     case ('bfvm')
       kappa_thic_type = kappa_type_bfreq_vmhs
     case ('bfhd')
       kappa_thic_type = kappa_type_bfreq_hdgr
     case ('bfdr')
       kappa_thic_type = kappa_type_bfreq_dradius
     case ('edgr')
       kappa_thic_type = kappa_type_eg
     case default
       kappa_thic_type = -1000
     end select

     select case (kappa_freq_choice(1:3))
     case ('nev')
       kappa_freq = kappa_freq_never
     case ('eve')
       kappa_freq = kappa_freq_every_time_step
     case ('onc')
       kappa_freq = kappa_freq_once_a_day
     case default
       kappa_freq = -1000
     end select

     select case (slope_control_choice(1:4))
     case ('tanh')
       slope_control = slope_control_tanh
     case ('nota')
       slope_control = slope_control_notanh
     case ('clip')
       slope_control = slope_control_clip
     case ('Gerd')
       slope_control = slope_control_Gerd
     case default
       slope_control = -1000
     end select

     if ( kappa_isop_type == kappa_type_bfreq          .or.  &
          kappa_isop_type == kappa_type_bfreq_vmhs     .or.  &
          kappa_isop_type == kappa_type_bfreq_hdgr     .or.  &
          kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
          kappa_thic_type == kappa_type_bfreq          .or.  &
          kappa_thic_type == kappa_type_bfreq_vmhs     .or.  &
          kappa_thic_type == kappa_type_bfreq_hdgr     .or.  &
          kappa_thic_type == kappa_type_bfreq_dradius ) then
       if ( read_n2_data ) then
         write(stdout,'(a32)') ' using climatological N^2 data. '
       else
         write(stdout,'(a34)') ' using model data to compute N^2. '
       endif
      endif

     if ( kappa_isop_type == kappa_type_depth  .or.  &
          kappa_thic_type == kappa_type_depth ) then

       if ( kappa_isop_type == kappa_type_depth )  &
         write(stdout,'(a30)') ' KAPPA_ISOP varies with depth.'

       if ( kappa_thic_type == kappa_type_depth )  &
         write(stdout,'(a30)') ' KAPPA_THIC varies with depth.'

       write(stdout,'(a20,1pe13.6)') '    kappa_depth_1 = ',  &
                                          kappa_depth_1
       write(stdout,'(a20,1pe13.6)') '    kappa_depth_2 = ',  &
                                          kappa_depth_2
       write(stdout,'(a24,1pe13.6)') '    kappa_depth_scale = ',  &
                                          kappa_depth_scale

     endif

     if ( kappa_isop_type == kappa_type_eg  .or.  &
          kappa_thic_type == kappa_type_eg ) then

       write(stdout,'(a16,1pe13.6)') '     const_eg = ',  &
                                           const_eg 
       write(stdout,'(a16,1pe13.6)') '     gamma_eg = ',  &
                                           gamma_eg
       write(stdout,'(a16,1pe13.6)') ' kappa_min_eg = ',  &
                                       kappa_min_eg
       write(stdout,'(a16,1pe13.6)') ' kappa_max_eg = ',  &
                                       kappa_max_eg

     endif


   endif

   call broadcast_scalar(kappa_isop_type,        master_task)
   call broadcast_scalar(kappa_thic_type,        master_task)
   call broadcast_scalar(kappa_freq,             master_task)
   call broadcast_scalar(slope_control,          master_task)
   call broadcast_scalar(kappa_depth_1,          master_task)
   call broadcast_scalar(kappa_depth_2,          master_task)
   call broadcast_scalar(kappa_depth_scale,      master_task)
   call broadcast_scalar(ah,                     master_task)
   call broadcast_scalar(ah_bolus,               master_task)
   call broadcast_scalar(ah_bkg_bottom,          master_task)
   call broadcast_scalar(ah_bkg_srfbl,           master_task)
   call broadcast_scalar(slm_r,                  master_task)
   call broadcast_scalar(slm_b,                  master_task)
   call broadcast_scalar(diag_gm_bolus,          master_task)
   call broadcast_scalar(use_const_ah_bkg_srfbl, master_task)
   call broadcast_scalar(transition_layer_on,    master_task)
   call broadcast_scalar(read_n2_data,           master_task)
   call broadcast_scalar(buoyancy_freq_filename, master_task)
   call broadcast_scalar(buoyancy_freq_fmt,      master_task)
   call broadcast_scalar(const_eg,               master_task)
   call broadcast_scalar(gamma_eg,               master_task)
   call broadcast_scalar(kappa_min_eg,           master_task)
   call broadcast_scalar(kappa_max_eg,           master_task)

   call broadcast_scalar(addrandfluc,            master_task)
   call broadcast_scalar(cflmajoronly,           master_task)
   call broadcast_scalar(cflmult,                master_task)
   call broadcast_scalar(erat_const,             master_task)
   call broadcast_scalar(erat_factor,            master_task)
   call broadcast_scalar(isominoronly,           master_task)
   call broadcast_scalar(isoonly,                master_task)
   call broadcast_scalar(kdir_type,              master_task)
   call broadcast_scalar(kmin_type,              master_task)
   call broadcast_scalar(krat_type,              master_task)
   call broadcast_scalar(minorfactor,            master_task)
   call broadcast_scalar(savenewtavgs,           master_task)
   call broadcast_scalar(shrdispfac,             master_task)
   call broadcast_scalar(simpsubcells,           master_task)
   call broadcast_scalar(vertdiffhere,           master_task)

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   if ( kappa_isop_type == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown type for KAPPA_ISOP in GM setup')
   endif
   if ( kappa_thic_type == -1000 ) then
    call exit_POP(sigAbort,  &
                  'unknown type for KAPPA_THIC in GM setup')
   endif
   if ( kappa_freq == -1000 ) then
     call exit_POP(sigAbort,                                       &
                  'unknown type for KAPPA computation frequency in GM setup')
   endif
   if ( slope_control == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown slope control method in GM setup')
   endif

   if ( kappa_isop_type == kappa_type_depth  .and.      &
        ( kappa_thic_type /= kappa_type_const  .and.     &
          kappa_thic_type /= kappa_type_depth ) ) then   
     message = 'when kappa_isop_type is kappa_type_depth, kappa_thic_type '/&
           &/  'should be either kappa_type_const or kappa_type_depth.'
     call exit_POP(sigAbort, message)
   endif

   if ( kappa_thic_type == kappa_type_depth  .and.      &
        ( kappa_isop_type /= kappa_type_const  .and.     &
          kappa_isop_type /= kappa_type_depth ) ) then   
     message = 'when kappa_thic_type is kappa_type_depth, kappa_isop_type ' /&
            &/ 'should be either kappa_type_const or kappa_type_depth.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( ( kappa_isop_type /= kappa_type_const   .and.     &
            kappa_isop_type /= kappa_type_depth ) .and.     &
          ( kappa_thic_type /= kappa_type_const   .and.     &
            kappa_thic_type /= kappa_type_depth  ) ) .and.  &
          ( kappa_isop_type /= kappa_thic_type ) ) then
     message = 'kappa_isop_type and kappa_thic_type should be the same ' /&
           &/  'if they are BOTH model fields dependent.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( ( kappa_isop_type == kappa_type_vmhs           .and.     &
            kappa_thic_type == kappa_type_vmhs )          .or.     &
          ( kappa_isop_type == kappa_type_bfreq_vmhs     .and.     &
            kappa_thic_type == kappa_type_bfreq_vmhs )    .or.     &
          ( kappa_isop_type == kappa_type_dradius        .and.     &
            kappa_thic_type == kappa_type_dradius )       .or.     &
          ( kappa_isop_type == kappa_type_bfreq_dradius  .and.     &
            kappa_thic_type == kappa_type_bfreq_dradius ) ) .and.  &
         ah /= ah_bolus ) then
     message = 'ah and ah_bolus should be equal when vmhs or dradius '  /&
          &/ ' related dependencies are used.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_depth   .or.  &
          kappa_thic_type == kappa_type_depth ) .and. & 
        kappa_depth_2 == c0 ) then
     message = 'kappa_depth_2 should be non-zero if a time-independent '  /&
          &/   'depth variation is requested.'
     call exit_POP(sigAbort, message)
      endif

   if ( ( kappa_isop_type == kappa_type_vmhs         .or.   &
          kappa_thic_type == kappa_type_vmhs         .or.   &
          kappa_isop_type == kappa_type_bfreq_vmhs   .or.   &
          kappa_thic_type == kappa_type_bfreq_vmhs   .or.   &
          kappa_isop_type == kappa_type_hdgr         .or.   &
          kappa_thic_type == kappa_type_hdgr         .or.   &
          kappa_isop_type == kappa_type_bfreq_hdgr   .or.   &
          kappa_thic_type == kappa_type_bfreq_hdgr ) .and.  &
                    zt(km) <= 2.0e5_r8 ) then
     message = 'the max bottom depth should be > 2000.0 m' /&
           &/ ' when vmhs or hdgr dependency is chosen.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_eg    .or.  &
          kappa_thic_type == kappa_type_eg )  .and. &
          use_const_ah_bkg_srfbl ) then
     message = 'use_const_ah_bkg_srfbl should be false when ' /&
            &/ 'kappa_type_eg is used.'
     call exit_POP(sigAbort, message)
   endif

   if ( ( kappa_isop_type == kappa_type_vmhs           .or.    &
          kappa_thic_type == kappa_type_vmhs           .or.    &
          kappa_isop_type == kappa_type_hdgr           .or.    &
          kappa_thic_type == kappa_type_hdgr           .or.    &
          kappa_isop_type == kappa_type_dradius        .or.    &
          kappa_thic_type == kappa_type_dradius        .or.    &
          kappa_isop_type == kappa_type_bfreq_vmhs     .or.    &
          kappa_thic_type == kappa_type_bfreq_vmhs     .or.    &
          kappa_isop_type == kappa_type_bfreq_hdgr     .or.    &
          kappa_thic_type == kappa_type_bfreq_hdgr     .or.    &
          kappa_isop_type == kappa_type_bfreq_dradius  .or.    &
          kappa_thic_type == kappa_type_bfreq_dradius  .or.    &
          kappa_isop_type == kappa_type_eg             .or.    &
          kappa_thic_type == kappa_type_eg             .or.    &
          ( ( kappa_isop_type == kappa_type_bfreq  .or.        &
              kappa_thic_type == kappa_type_bfreq ) .and.      &
            .not. read_n2_data ) )                    .and.    &
          kappa_freq == kappa_freq_never ) then

     message = 'kappa_freq should not be set to never when model ' /&
            &/ 'fields dependent kappa types are chosen.'
     call exit_POP(sigAbort, message)
   endif

!  SJR DEL  if ( partial_bottom_cells ) then
!  SJR DEL    call exit_POP(sigAbort, &
!  SJR DEL     'hmix_gm currently incompatible with partial bottom cells')
!  SJR DEL  endif

   if (read_n2_data .and. trim(buoyancy_freq_filename) == 'unknown-buoyancy' ) then
     call exit_POP(sigAbort,  &
                   'Must define buoyancy_freq_filename if read_n2_data is .true.')
   endif

   if ( kmin_type == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown type for minor diffusivity KMIN_TYPE in GM Aniso setup')
   endif
   if ( krat_type == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown type for diffusivity ratio KRAT_TYPE in GM Aniso setup')
   endif
   if ( kdir_type == -1000 ) then
     call exit_POP(sigAbort,  &
                   'unknown type for major axis direction KDIR_TYPE in GM Aniso setup')
   endif

!-----------------------------------------------------------------------
!
!  allocate GM arrays
!
!-----------------------------------------------------------------------


    allocate (HYXW(nx_block,ny_block,nblocks_clinic),    &
             HXYS(nx_block,ny_block,nblocks_clinic),    &
             RBR (nx_block,ny_block,nblocks_clinic),    &
             BTP (nx_block,ny_block,nblocks_clinic),    &
             BL_DEPTH(nx_block,ny_block,nblocks_clinic))

    allocate (SF_SLXX(nx_block,ny_block,2,2,km,nblocks_clinic),  & !SJR MOD
              SF_SLYY(nx_block,ny_block,2,2,km,nblocks_clinic),  & !SJR MOD
              SF_SLXY(nx_block,ny_block,2,2,km,nblocks_clinic),  & !SJR ADD
              SF_SLYX(nx_block,ny_block,2,2,km,nblocks_clinic),  & !SJR ADD
              DD_SLX (nx_block,ny_block,2,2,km,nblocks_clinic),  & !SJR MOD
              DD_SLY (nx_block,ny_block,2,2,km,nblocks_clinic))    !SJR ADD
    
    allocate (FZTOP(nx_block,ny_block,nt,nblocks_clinic))
    if ( partial_bottom_cells ) allocate (DZTW(nx_block,ny_block,0:km,nblocks_clinic)) !SJR ADD

    allocate (kappa_depth(km))

    allocate (KAPPA_ISOP(nx_block,ny_block,2,km,nblocks_clinic),  &
              KAPPA_THIC(nx_block,ny_block,2,km,nblocks_clinic),  &
              KXX_ISOP(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              KXY_ISOP(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              KYY_ISOP(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              KXX_THIC(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              KXY_THIC(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              KYY_THIC(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              HXX_DIFF(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              HXY_DIFF(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              HYY_DIFF(nx_block,ny_block,2,km,nblocks_clinic),  &   !SJR ADD
              K_EIGENVAL_RAT(nx_block,ny_block,2,km,nblocks_clinic), & !SJR ADD
              HOR_DIFF(nx_block,ny_block,2,km,nblocks_clinic))

    allocate (KXX_IN(nx_block,ny_block,km,nblocks_clinic), & !SJR ADD
              KXY_IN(nx_block,ny_block,km,nblocks_clinic), & !SJR ADD
              KYY_IN(nx_block,ny_block,km,nblocks_clinic), & !SJR ADD
              MAJOR_EIGENVAL(nx_block,ny_block,km,nblocks_clinic), & !SJR ADD
              MINOR_EIGENVAL(nx_block,ny_block,km,nblocks_clinic), & !SJR ADD
              NX_ISOP(nx_block,ny_block,km,nblocks_clinic),        & !SJR ADD
              NY_ISOP(nx_block,ny_block,km,nblocks_clinic) )         !SJR ADD

    if (kdir_type == kdir_type_shrd .or. krat_type == krat_type_shrd) then
       allocate (NX_SHRD(nx_block,ny_block,km,nblocks_clinic),        & !SJR ADD
                 NY_SHRD(nx_block,ny_block,km,nblocks_clinic) )         !SJR ADD
       allocate (KRAT_SHRD(nx_block,ny_block,2,km,nblocks_clinic))  !SJR ADD

    endif

    allocate (SUBCELLV(nx_block,ny_block,2,2,2,km,nblocks_clinic)) !SJR ADD
    allocate (SL_MAG(nx_block,ny_block,2,2,2,km,nblocks_clinic)) !SJR ADD

    allocate (KAPPA_LATERAL (nx_block,ny_block,nblocks_clinic),  &
             KAPPA_VERTICAL(nx_block,ny_block,km,nblocks_clinic))

    allocate (BUOY_FREQ_SQ(nx_block,ny_block,km,nblocks_clinic))

    allocate (VDC_GM(nx_block,ny_block,km,nblocks_clinic))

    allocate (compute_kappa(nblocks_clinic))

   HYXW     = c0
   HXYS     = c0
   RBR      = c0
   BTP      = c0
   BL_DEPTH = c0
   SF_SLXX   = c0 !SJR MOD
   SF_SLYY   = c0 !SJR MOD
   SF_SLXY   = c0 !SJR ADD
   SF_SLYX   = c0 !SJR ADD
   DD_SLX    = c0 !SJR ADD
   DD_SLY    = c0 !SJR ADD
   FZTOP    = c0
   VDC_GM   = c0
   SUBCELLV = c0 !SJR ADD
   SL_MAG   = c0 !SJR ADD
   
   if ( partial_bottom_cells ) DZTW(:,:,0:km,:) = p5*(DZT(:,:,0:km,:) + DZT(:,:,1:km+1,:)) !SJR ADD

   if ( transition_layer_on ) then
     allocate (SLA_SAVE(nx_block,ny_block,2,km,nblocks_clinic))
     allocate (RB(nx_block,ny_block,nblocks_clinic))
     SLA_SAVE = c0
     RB = c0
   endif

!-----------------------------------------------------------------------
!
!  initialize various time-independent arrays
!
!-----------------------------------------------------------------------

   do k=1,km
     kappa_depth(k) = kappa_depth_1  &
                    + kappa_depth_2  &
                     *exp(-zt(k)/kappa_depth_scale)
   enddo

   do iblock = 1,nblocks_clinic

     this_block = get_block(blocks_clinic(iblock),iblock)

     KAPPA_LATERAL(:,:,iblock)    = ah
     KAPPA_VERTICAL(:,:,:,iblock) = c1

     KAPPA_ISOP(:,:,:,:,iblock) = ah
     KAPPA_THIC(:,:,:,:,iblock) = ah_bolus
     HOR_DIFF(:,:,:,:,iblock)   = ah

     if ( kappa_isop_type == kappa_type_depth  .or.  &
          kappa_thic_type == kappa_type_depth ) then

       do k=1,km 
         KAPPA_VERTICAL(:,:,k,iblock) = kappa_depth(k)
       enddo

     endif


     HYXW(:,:,iblock) = eoshift(HYX(:,:,iblock), dim=1, shift=-1)
     HXYS(:,:,iblock) = eoshift(HXY(:,:,iblock), dim=2, shift=-1)

!-----------------------------------------------------------------------
!  compute the Rossby radius which will be used to
!  contol KAPPA [Large et al (1997), JPO, 27,
!  pp 2418-2447]. Rossby radius = c/(2*omg*sin(latitude))
!  where c=200cm/s is the first baroclinic wave speed.
!  15km < Rossby radius < 100km
!-----------------------------------------------------------------------

     !*** Inverse of Rossby radius

     RBR(:,:,iblock) = abs(FCORT(:,:,iblock))  &
                       / 200.0_r8             ! |f|/Cg, Cg = 2 m/s = 200 cm/s
     RBR(:,:,iblock) = min(RBR(:,:,iblock),    &
                           c1/1.5e+6_r8)      ! Cg/|f| .ge. 15 km = 1.5e+6 cm
     RBR(:,:,iblock) = max(RBR(:,:,iblock),    &
                           1.e-7_r8)          ! Cg/|f| .le. 100 km = 1.e+7 cm

     if ( transition_layer_on ) then
       RB(:,:,iblock) = c1 / RBR(:,:,iblock)
     endif

     !*** beta at t-points

     call ugrid_to_tgrid(BTP(:,:,iblock),ULAT(:,:,iblock),iblock)
      
     BTP(:,:,iblock) = c2*omega*cos(BTP(:,:,iblock))/radius

!-------------v--SJR ADD--v--------------------------------------------- 
!-------------v--SJR ADD--v--------------------------------------------- 
     if (simpsubcells) then
        do k=1,km
           if ( partial_bottom_cells ) then
              SUBCELLV(:,:,1,1,1,k,iblock) = p125 * DXT(:,:,iblock) * DYT(:,:,iblock) * DZT(:,:,k,iblock)
           else
              SUBCELLV(:,:,1,1,1,k,iblock) = p125 * DXT(:,:,iblock) * DYT(:,:,iblock) * dz(k)
           endif
           SUBCELLV(:,:,2,1,1,k,iblock) = SUBCELLV(:,:,1,1,1,k,iblock) 
           SUBCELLV(:,:,:,2,1,k,iblock) = SUBCELLV(:,:,:,1,1,k,iblock) 
           SUBCELLV(:,:,:,:,2,k,iblock) = SUBCELLV(:,:,:,:,1,k,iblock) 
        enddo
     else
        do k=1,km
           do j_sub=1,2
              jstrt = 3-j_sub
              jstop = ny_block-j_sub+1
              do i_sub=1,2
                 istrt = 3-i_sub
                 istop = nx_block-i_sub+1
                 SUBCELLV(2:nx_block,2:ny_block,i_sub,j_sub,1,k,iblock) = &
                              p125 * HTE(istrt:istop,2:ny_block,iblock) * &
                                     HTN(2:nx_block,jstrt:jstop,iblock) * dz(k)
                 SUBCELLV(1,2:ny_block,i_sub,j_sub,1,k,iblock) = &
                              p125 * DXT(1,2:ny_block,iblock) * &
                                     DYT(1,2:ny_block,iblock) * dz(k)
                 SUBCELLV(:,1,i_sub,j_sub,1,k,iblock) = &
                              p125 * DXT(:,1,iblock) * &
                                     DYT(:,1,iblock) * dz(k)
                 if ( partial_bottom_cells ) &
                    SUBCELLV(:,:,i_sub,j_sub,1,k,iblock) = SUBCELLV(:,:,i_sub,j_sub,1,k,iblock) * DZT(:,:,k,iblock) / dz(k)
              enddo
           enddo
        enddo
        SUBCELLV(:,:,:,:,2,:,iblock) = SUBCELLV(:,:,:,:,1,:,iblock)
     endif
!-------------^--SJR ADD--^--------------------------------------------- 

   enddo

!-------------v--SJR ADD--v--------------------------------------------- 
   if (kmin_type == kmin_type_read .OR. krat_type == krat_type_read .OR. kdir_type == kdir_type_read) then

      KXX_filename = '/glade/scratch/scotreck/Kfiles/KXX_sm.nc'
      KXY_filename = '/glade/scratch/scotreck/Kfiles/KXY_sm.nc'
      KYY_filename = '/glade/scratch/scotreck/Kfiles/KYY_sm.nc'

      KXX_data_file = construct_file ( 'nc', full_name=trim(KXX_filename) )

      i_dim = construct_io_dim ( 'nlon', nx_global )
      j_dim = construct_io_dim ( 'nlat', ny_global )
      k_dim = construct_io_dim ( 'nz', km )

      call data_set(KXX_data_file, 'open_read')
      KXX_data_in = construct_io_field('KXX',             &
                     dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                     field_loc = field_loc_center,        &
                     field_type = field_type_scalar,      &
                     d3d_array = KXX_IN)
      call data_set (KXX_data_file, 'define', KXX_data_in)
      call data_set (KXX_data_file, 'read',   KXX_data_in)
      call data_set (KXX_data_file, 'close')
      call destroy_io_field (KXX_data_in)
      call destroy_file (KXX_data_file)
 
      KXY_data_file = construct_file ( 'nc', full_name=trim(KXY_filename) )
      call data_set(KXY_data_file, 'open_read')
      KXY_data_in = construct_io_field('KXY',             &
                     dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                     field_loc = field_loc_center,        &
                     field_type = field_type_scalar,      &
                     d3d_array = KXY_IN)
      call data_set (KXY_data_file, 'define', KXY_data_in)
      call data_set (KXY_data_file, 'read',   KXY_data_in)
      call data_set (KXY_data_file, 'close')
      call destroy_io_field (KXY_data_in)
      call destroy_file (KXY_data_file)

      KYY_data_file = construct_file ( 'nc', full_name=trim(KYY_filename) )
      call data_set(KYY_data_file, 'open_read')
      KYY_data_in = construct_io_field('KYY',             &
                     dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                     field_loc = field_loc_center,        &
                     field_type = field_type_scalar,      &
                     d3d_array = KYY_IN)
      call data_set (KYY_data_file, 'define', KYY_data_in)
      call data_set (KYY_data_file, 'read',   KYY_data_in)
      call data_set (KYY_data_file, 'close')
      call destroy_io_field (KYY_data_in)
      call destroy_file (KYY_data_file)

      NX_ISOP = sqrt((KXX_IN+KYY_IN)**2-c4*(KXX_IN*KYY_IN-KXY_IN**2))/c2 !EIGENVALUES are +/- this
      MAJOR_EIGENVAL = (KXX_IN+KYY_IN)/c2+NX_ISOP
      MINOR_EIGENVAL = (KXX_IN+KYY_IN)/c2-NX_ISOP
      KYY_IN = KXY_IN**2 + (MAJOR_EIGENVAL - KXX_IN)**2 !NORMALIZATION FACTOR
      where (KYY_IN < eps2)
         NY_ISOP = c0
         NX_ISOP = c1
      elsewhere
         NY_ISOP = (MAJOR_EIGENVAL - KXX_IN)/sqrt(KYY_IN)
         NX_ISOP = KXY_IN/sqrt(KYY_IN)
      endwhere

!---------------------------------------------------
      MAJOR_EIGENVAL = MAX(eps,MAJOR_EIGENVAL) !FORCE MAJOR TO BE POSITIVE
      if (erat_factor<eps) then
         MINOR_EIGENVAL = MAX(eps,MINOR_EIGENVAL) !FORCE MINOR TO BE POSITIVE
      else
         MINOR_EIGENVAL = MAX(-MAJOR_EIGENVAL/erat_factor,MINOR_EIGENVAL) !ALLOW NEGATIVE MINOR, BUT MAGNITUDE LESS THAN MAJOR/erat_factor
         where (ABS(MINOR_EIGENVAL)<eps )
            MINOR_EIGENVAL = eps
         endwhere
      endif
      KXX_IN = MAJOR_EIGENVAL*NX_ISOP**2 + MINOR_EIGENVAL*NY_ISOP**2 
      KYY_IN = MINOR_EIGENVAL*NX_ISOP**2 + MAJOR_EIGENVAL*NY_ISOP**2 
      KXY_IN = (MAJOR_EIGENVAL-MINOR_EIGENVAL)*NX_ISOP*NY_ISOP 

      K_EIGENVAL_RAT(:,:,1,:,:) = MAJOR_EIGENVAL/MINOR_EIGENVAL
      K_EIGENVAL_RAT(:,:,2,:,:) = K_EIGENVAL_RAT(:,:,1,:,:)  

!---------------------------------------------------

      if (isoonly) then
         if (isominoronly) then
            MINOR_EIGENVAL = MAX(eps,MINOR_EIGENVAL)
         else
            MINOR_EIGENVAL = MAX(eps,sqrt(MINOR_EIGENVAL*MAJOR_EIGENVAL))
         endif
         KXX_IN = MINOR_EIGENVAL
         KYY_IN = MINOR_EIGENVAL
         KXY_IN = c0
         K_EIGENVAL_RAT = c1
         MAJOR_EIGENVAL = MINOR_EIGENVAL
      endif
   endif
   if ( krat_type == krat_type_simp ) then
      K_EIGENVAL_RAT = erat_const
   endif
!-------------^--SJR ADD--^--------------------------------------------- 

!-----------------------------------------------------------------------
!  HYXW, HXYS only needed in physical domain and should
!  be defined correctly there.  BTP invalid on south
!  and westernmost ghost cells, but not needed there
!  as long as number of ghost cells is >1
!-----------------------------------------------------------------------
 
!  call update_ghost_cells(HYXW, bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)
!  call update_ghost_cells(HXYS, bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)
!  call update_ghost_cells(BTP , bndy_clinic, field_loc_t,     &
!                                             field_type_scalar)

   BUOY_FREQ_SQ = c0

   if ( read_n2_data ) then

     buoyancy_freq_filename = '/home/bluesky/gokhan/buoyancy_freq.nc'

     buoyancy_freq_data_file = construct_file ( 'nc',         &
                          full_name=trim(buoyancy_freq_filename) )

     call data_set(buoyancy_freq_data_file, 'open_read')

     i_dim = construct_io_dim ( 'i', nx_global )
     j_dim = construct_io_dim ( 'j', ny_global )
     k_dim = construct_io_dim ( 'k', km )

     buoyancy_freq_data_in = construct_io_field                &
                         ('BUOYANCY_FREQUENCY',                &
                          dim1=i_dim, dim2=j_dim, dim3=k_dim,  &
                          field_loc = field_loc_center,        &
                          field_type = field_type_scalar,      &
                          d3d_array = BUOY_FREQ_SQ)

     call data_set (buoyancy_freq_data_file, 'define', &
                    buoyancy_freq_data_in)
     call data_set (buoyancy_freq_data_file, 'read',   &
                    buoyancy_freq_data_in)
     call data_set (buoyancy_freq_data_file, 'close')
     call destroy_io_field (buoyancy_freq_data_in)
     call destroy_file (buoyancy_freq_data_file)

     if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,'(a43,a)')                                 &
               ' Buoyancy frequency climatology file read: ',  &
                              trim(buoyancy_freq_filename)
     endif

   endif

   compute_kappa = .false.

   if (slm_r /= slm_b) then
     diff_tapering = .true.
   else
     diff_tapering = .false.
   endif

   cancellation_occurs = .true.

   if ( ( kappa_isop_type /= kappa_thic_type )  .or.     &
       ( diff_tapering )  .or.                           &
       ( ( kappa_isop_type == kappa_type_const  .and.    &
           kappa_thic_type == kappa_type_const )  .and.  &
         ah /= ah_bolus )  .or.                          &
       ( ( kappa_isop_type == kappa_type_depth  .and.    &
           kappa_thic_type == kappa_type_depth )  .and.  &
         ah /= ah_bolus )  .or.                          &
       ( ( kappa_isop_type == kappa_type_bfreq  .and.    &
           kappa_thic_type == kappa_type_bfreq )  .and.  &
         ah /= ah_bolus ) )                              &
    cancellation_occurs = .false. 

    !*** for transition layer cases, the following will always be true!!

   if ( transition_layer_on )  cancellation_occurs = .false.

!-----------------------------------------------------------------------
!
!  initialize topography mask used with kappa_type_eg. This mask eliminates
!  excessively large values of KAPPA near sloping topography.
!
!-----------------------------------------------------------------------

   if ( kappa_isop_type == kappa_type_eg  .or.  &
        kappa_thic_type == kappa_type_eg ) then 

     allocate (SIGMA_TOPO_MASK(nx_block,ny_block,km,nblocks_clinic))

     do iblock=1,nblocks_clinic

       do k=1,km
         where ( k < KMT(:,:,iblock) ) 
           SIGMA_TOPO_MASK(:,:,k,iblock) = c1
         elsewhere
           SIGMA_TOPO_MASK(:,:,k,iblock) = c0
         endwhere
       enddo

       do k=1,km-1
         do j=2,ny_block-1
           do i=2,nx_block-1 
             if ( k < KMT(i,j,iblock) ) then
               if ( k == KMT(i-1,j+1,iblock)  .or.  &
                    k == KMT(i  ,j+1,iblock)  .or.  &
                    k == KMT(i+1,j+1,iblock)  .or.  &
                    k == KMT(i-1,j  ,iblock)  .or.  &
                    k == KMT(i+1,j  ,iblock)  .or.  &
                    k == KMT(i-1,j-1,iblock)  .or.  &
                    k == KMT(i  ,j-1,iblock)  .or.  &
                    k == KMT(i+1,j-1,iblock) )      &
                 SIGMA_TOPO_MASK(i,j,k,iblock) = c0 
             endif
           enddo
         enddo
       enddo

     enddo

   endif

!-----------------------------------------------------------------------
!
!  define tavg fields related to diffusivities 
!
!-----------------------------------------------------------------------

   call define_tavg_field (tavg_KAPPA_ISOP, 'KAPPA_ISOP', 3,     &
                   long_name='Isopycnal diffusion coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

  if (savenewtavgs) then
   call define_tavg_field (tavg_KXX_ISOP, 'KXX_ISOP', 3,     &
                   long_name='Isopycnal diffusion XX coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_KXY_ISOP, 'KXY_ISOP', 3,     &
                   long_name='Isopycnal diffusion XY coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_KYY_ISOP, 'KYY_ISOP', 3,     &
                   long_name='Isopycnal diffusion YY coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MIN_ISOP, 'MINOR_ISOP', 3,     &
                   long_name='Isopycnal minor diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MAJ_ISOP, 'MAJOR_ISOP', 3,     &
                   long_name='Isopycnal major diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_C2T_ISOP, 'COS2T_ISOP', 3,     &
                   long_name='Isopycnal cos(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_S2T_ISOP, 'SIN2T_ISOP', 3,     &
                   long_name='Isopycnal sin(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')
  endif

   call define_tavg_field (tavg_KAPPA_THIC, 'KAPPA_THIC', 3,     &
                   long_name='Thickness diffusion coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

  if (savenewtavgs) then
   call define_tavg_field (tavg_KXX_THIC, 'KXX_THIC', 3,     &
                   long_name='Thickness diffusion XX coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_KXY_THIC, 'KXY_THIC', 3,     &
                   long_name='Thickness diffusion XY coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_KYY_THIC, 'KYY_THIC', 3,     &
                   long_name='Thickness diffusion YY coefficient',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MIN_THIC, 'MINOR_THIC', 3,     &
                   long_name='Thickness minor diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MAJ_THIC, 'MAJOR_THIC', 3,     &
                   long_name='Thickness major diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_C2T_THIC, 'COS2T_THIC', 3,     &
                   long_name='Thickness cos(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_S2T_THIC, 'SIN2T_THIC', 3,     &
                   long_name='Thickness sin(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')
  endif

   call define_tavg_field (tavg_HOR_DIFF, 'HOR_DIFF', 3,         &
                   long_name='Horizontal diffusion coefficient', &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

  if (savenewtavgs) then
   call define_tavg_field (tavg_HXX_DIFF, 'HXX_DIFF', 3,         &
                   long_name='Horizontal diffusion XX coefficient', &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_HXY_DIFF, 'HXY_DIFF', 3,         &
                   long_name='Horizontal diffusion XY coefficient', &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_HYY_DIFF, 'HYY_DIFF', 3,         &
                   long_name='Horizontal diffusion YY coefficient', &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MIN_DIFF, 'MINOR_DIFF', 3,     &
                   long_name='Horizontal minor diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_MAJ_DIFF, 'MAJOR_DIFF', 3,     &
                   long_name='Horizontal major diffusivity',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_C2T_DIFF, 'COS2T_DIFF', 3,     &
                   long_name='Thickness cos(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')

   call define_tavg_field (tavg_S2T_DIFF, 'SIN2T_DIFF', 3,     &
                   long_name='Thickness sin(2*theta) major axis angle',  &
                   units='cm^2/s', grid_loc='3111',              &
                   coordinates='TLONG TLAT z_t time')
  endif

!-----------------------------------------------------------------------
!
!  define tavg fields related to the near-surface eddy flux
!  parameterization (i.e. transition layer is on). 
!
!-----------------------------------------------------------------------

   if ( transition_layer_on ) then

     call define_tavg_field (tavg_DIA_DEPTH, 'DIA_DEPTH', 2,       &
         long_name='Depth of the Diabatic Region at the Surface',  &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_TLT, 'TLT', 2,                   &
         long_name='Transition Layer Thickness',                   &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_INT_DEPTH, 'INT_DEPTH', 2,       &
         long_name='Depth at which the Interior Region Starts',    &
                     units='cm', grid_loc='2110',                  &
                     coordinates='TLONG TLAT time')

   endif

!-----------------------------------------------------------------------
!
!  allocate and initialize bolus velocity arrays if requested
!  define tavg fields related to bolus velocity
!
!-----------------------------------------------------------------------

   if ( diag_gm_bolus ) then

     call register_string ('diag_gm_bolus')

     allocate(WTOP_ISOP(nx_block,ny_block,nblocks_clinic), &
              WBOT_ISOP(nx_block,ny_block,nblocks_clinic), &
                    UIT(nx_block,ny_block,nblocks_clinic), &
                    VIT(nx_block,ny_block,nblocks_clinic))

     WTOP_ISOP = c0
     WBOT_ISOP = c0
     UIT       = c0
     VIT       = c0

     call define_tavg_field (tavg_UISOP, 'UISOP', 3,                &
      long_name='Bolus Velocity in grid-x direction (diagnostic)',  &
                   units='cm/s', grid_loc='3211',                   &
                   coordinates='ULONG TLAT z_t time')

     call define_tavg_field (tavg_VISOP, 'VISOP', 3,                &
      long_name='Bolus Velocity in grid-y direction (diagnostic)',  &
                   units='cm/s', grid_loc='3121',                   &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_WISOP, 'WISOP', 3,                &
      long_name='Vertical Bolus Velocity (diagnostic)',             &
                   units='cm/s', grid_loc='3112',                   &
                   coordinates='TLONG TLAT z_w time')

     call define_tavg_field (tavg_ADVT_ISOP, 'ADVT_ISOP', 2,                            &
      long_name='Vertically-Integrated T Eddy-Induced Advection Tendency (diagnostic)', &
                   units='cm degC/s', grid_loc='2110',                                  &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_ADVS_ISOP, 'ADVS_ISOP', 2,                            &
      long_name='Vertically-Integrated S Eddy-Induced Advection Tendency (diagnostic)', &
                   scale_factor=1000.0_r8,                                              &
                   units='cm gram/kilogram/s', grid_loc='2110',                         &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_VNT_ISOP, 'VNT_ISOP', 3,                               &
      long_name='Heat Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   units='degC/s', grid_loc='3121',                                      &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_VNS_ISOP, 'VNS_ISOP', 3,                               &
      long_name='Salt Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   scale_factor=1000.0_r8,                                               &
                   units='gram/kilogram/s', grid_loc='3121',                             &
                   coordinates='TLONG ULAT z_t time')

   endif

      call get_timer(timer_gm0,'GM0',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm1,'GM1',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm2,'GM2',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm3,'GM3',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm4,'GM4',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm5,'GM5',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm6,'GM6',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm7,'GM7',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm7a,'GM7a',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm7b,'GM7b',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm8,'GM8',nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(timer_gm9,'GM9',nblocks_clinic, distrb_clinic%nprocs)

      call get_timer(timer_nloop,'HMIX_TRACER_GM_NLOOP', &
                                  nblocks_clinic, distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

   end subroutine init_gm_aniso

!***********************************************************************
!BOP
! !IROUTINE: hdifft_gm_aniso
! !INTERFACE:

      subroutine hdifft_gm_aniso (k, GTK, TMIX, UMIX, VMIX, tavg_HDIFE_TRACER, &
                            tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Gent-McWilliams eddy transport parameterization
!  and isopycnal diffusion.
!
!  This routine must be called successively with k = 1,2,3,...
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      integer (int_kind), intent(in) :: k  ! depth level index

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX            ! U,V  at all vertical levels
                               !   at mixing time level

      integer (int_kind), dimension(nt), intent(in) :: &
         tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
         tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
         tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

! !OUTPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,nt), intent(out) :: &
         GTK     ! diffusion+bolus advection for nth tracer at level k

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind), parameter :: &
         ieast  = 1, iwest  = 2,       &
         jnorth = 1, jsouth = 2

      integer (int_kind) :: &
         i,j,n,kk,          &! dummy loop counters
         kid,ktmp,          &! array indices
         kk_sub, kp1,       & 
         kk_shft, kz_shft,  &  ! SJR ADD 
         istrt, istop,      &  ! SJR ADD 
         jstrt, jstop,      &  ! SJR ADD
         astrt, astop,      &  ! SJR ADD   
         bstrt, bstop,      &  ! SJR ADD   
         kstrt, kstop,      &  ! SJR ADD   
         i_sub, j_sub,      &  ! SJR ADD
         k_sub, k_shft,     &  ! SJR ADD   
         bid                 ! local block address for this sub block

      real (r8) :: &
         fz, dz_bottom, factor, globval

      real (r8), dimension(nx_block,ny_block) :: &
         CX, CY,                  &
         RZ,                      &! Dz(rho)
         SLA,                     &! absolute value of slope
         WORK1, WORK2,            &! local work space
         WORK3, WORK4,            &! local work space
         WORK5, WORK6,            &! local work space
         KMASK,                   &! ocean mask
         KMASKE, KMASKN,          &! ocean mask !SJR ADD
         TAPER1, TAPER2, TAPER3,  &! tapering factors
         UIB, VIB,                &! work arrays for isopycnal mixing velocities
         U_ISOP, V_ISOP            ! horizontal components of isopycnal velocities

      real (r8), dimension(nx_block,ny_block,nt) :: &
         FZBOT, &  !SJR ADD   ! vertical flux
         FX, FY                     ! fluxes across east, north faces

      real (r8), dimension(2) :: &
         reference_depth

      real (r8), dimension(nx_block,ny_block) :: &            !SJR ADD
         GRADX_FCOR, GRADY_FCOR                               !SJR ADD

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      call timer_start(timer_gm0, block_id=this_block%local_id)

      bid = this_block%local_id

      U_ISOP = c0
      V_ISOP = c0
      WORK1  = c0
      WORK2  = c0
      WORK3  = c0
      WORK4  = c0

      if ( .not. implicit_vertical_mix )  &
        call exit_POP (sigAbort, &
         'implicit vertical mixing must be used with GM horiz mixing')

      if ( k == 1 ) then
      call timer_start(timer_gm1, block_id=this_block%local_id)
      call timer_start(timer_gm2, block_id=this_block%local_id)

        if ( diag_gm_bolus ) then
          UIB = c0
          VIB = c0
          UIT(:,:,bid) = c0
          VIT(:,:,bid) = c0
          WBOT_ISOP(:,:,bid) = c0
        endif

        HOR_DIFF(:,:,ktp,k,bid) = ah_bkg_srfbl

        BL_DEPTH(:,:,bid) = zw(k)
        if ( vmix_itype == vmix_type_kpp )  &
                    BL_DEPTH(:,:,bid) = KPP_HBLT(:,:,bid)

      endif

      if (diag_gm_bolus)  WTOP_ISOP(:,:,bid) = WBOT_ISOP(:,:,bid)

      CX = merge(HYX(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTE(:,:,bid)))
      CY = merge(HXY(:,:,bid)*p25, c0, (k <= KMT (:,:,bid))   &
                                 .and. (k <= KMTN(:,:,bid)))

      if ( k == 1 ) then

        if ( transition_layer_on ) then

          if ( vmix_itype == vmix_type_kpp ) then

            call smooth_hblt ( .false., .true., bid,  &
                               SMOOTH_OUT=TLT%DIABATIC_DEPTH(:,:,bid) )
          else
            TLT%DIABATIC_DEPTH(:,:,bid) = zw(k)
          endif

          do kk=1,km
            do kk_sub = ktp,kbt
              kid = kk + kk_sub - 2
              if ( partial_bottom_cells ) then
                 SLA_SAVE(:,:,kk_sub,kk,bid) = DZTW(:,:,kid,bid)*sqrt(p5*(       &
                        (SLX(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                      + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                      + eps
              else
                 SLA_SAVE(:,:,kk_sub,kk,bid) = dzw(kid)*sqrt(p5*(       &
                        (SLX(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                      + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                      + eps
              endif
            enddo
          enddo

          call transition_layer ( this_block )

        endif

	

!-----------------------------------------------------------------------
!
!     compute isopycnal and thickness diffusion coefficients if
!     they depend on the model fields
!
!-----------------------------------------------------------------------
        
	if ( ( kappa_isop_type == kappa_type_vmhs           .or.    &
               kappa_thic_type == kappa_type_vmhs           .or.    &
               kappa_isop_type == kappa_type_hdgr           .or.    &
               kappa_thic_type == kappa_type_hdgr           .or.    &
               kappa_isop_type == kappa_type_dradius        .or.    &
               kappa_thic_type == kappa_type_dradius        .or.    &
               kappa_isop_type == kappa_type_bfreq          .or.    &
               kappa_thic_type == kappa_type_bfreq          .or.    &
               kappa_isop_type == kappa_type_bfreq_vmhs     .or.    &
               kappa_thic_type == kappa_type_bfreq_vmhs     .or.    &
               kappa_isop_type == kappa_type_bfreq_hdgr     .or.    &
               kappa_thic_type == kappa_type_bfreq_hdgr     .or.    &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.    &
               kappa_thic_type == kappa_type_bfreq_dradius  .or.    &
               kappa_isop_type == kappa_type_eg             .or.    &
               kappa_thic_type == kappa_type_eg )            .and.  &
           ( ( kappa_freq == kappa_freq_every_time_step )           &
        .or. ( kappa_freq == kappa_freq_once_a_day .and. eod_last ) &
        .or. ( nsteps_total == 1 ) ) )  compute_kappa(bid) = .true.

        if ( compute_kappa(bid) ) then

          if ( kappa_isop_type == kappa_type_vmhs        .or.  &
               kappa_thic_type == kappa_type_vmhs        .or.  &
               kappa_isop_type == kappa_type_bfreq_vmhs  .or.  & 
               kappa_thic_type == kappa_type_bfreq_vmhs ) then

            if ( nsteps_total == 1 ) then
              KAPPA_LATERAL(:,:,bid) = ah
              if ( kappa_isop_type == kappa_type_const )  &
                KAPPA_LATERAL(:,:,bid) = ah_bolus
            else
              call kappa_lon_lat_vmhs (TMIX, UMIX, VMIX, this_block)
            endif

          endif

          if ( kappa_isop_type == kappa_type_hdgr        .or.  &
               kappa_thic_type == kappa_type_hdgr        .or.  &
               kappa_isop_type == kappa_type_bfreq_hdgr  .or.  &
               kappa_thic_type == kappa_type_bfreq_hdgr )      &
            call kappa_lon_lat_hdgr (TMIX, this_block)

          if ( kappa_isop_type == kappa_type_dradius        .or.  &
               kappa_thic_type == kappa_type_dradius        .or.  &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
               kappa_thic_type == kappa_type_bfreq_dradius )      &
            call kappa_lon_lat_dradius (this_block)

          if ( kappa_isop_type == kappa_type_bfreq          .or.  &
               kappa_thic_type == kappa_type_bfreq          .or.  &
               kappa_isop_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_thic_type == kappa_type_bfreq_vmhs     .or.  &
               kappa_isop_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_thic_type == kappa_type_bfreq_hdgr     .or.  &
               kappa_isop_type == kappa_type_bfreq_dradius  .or.  &
               kappa_thic_type == kappa_type_bfreq_dradius )      &
            call buoyancy_frequency_dependent_profile (TMIX, this_block)

          if ( kappa_isop_type == kappa_type_eg  .or.  &
               kappa_thic_type == kappa_type_eg ) &
            call kappa_eg (TMIX, UMIX, VMIX, this_block) 

          compute_kappa(bid) = .false.

        endif  ! end of ( compute_kappa ) if statement


!-----------------------------------------------------------------------
!
!     reinitialize the diffusivity coefficients 
!
!-----------------------------------------------------------------------
	
        if ( kappa_isop_type == kappa_type_const ) then
          KAPPA_ISOP(:,:,:,:,bid) = ah
        elseif ( kappa_isop_type == kappa_type_eg ) then
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_ISOP(:,:,kk_sub,kk,bid) = KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        else
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_ISOP(:,:,kk_sub,kk,bid) =  KAPPA_LATERAL(:,:,bid)  &
                                         * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        endif

        if ( .not. use_const_ah_bkg_srfbl )  &
          HOR_DIFF(:,:,ktp,k,bid) = KAPPA_ISOP(:,:,ktp,k,bid) 

        if ( kappa_thic_type == kappa_type_const ) then
          KAPPA_THIC(:,:,:,:,bid) = ah_bolus
        else if ( kappa_thic_type == kappa_type_depth  .or.  &
                  kappa_thic_type == kappa_type_bfreq ) then
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_THIC(:,:,kk_sub,kk,bid) =  ah_bolus  &
                                        * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo 
        else if ( kappa_thic_type == kappa_type_eg ) then
          KAPPA_THIC(:,:,:,:,bid) = KAPPA_ISOP(:,:,:,:,bid)
        else
          do kk_sub=ktp,kbt
            do kk=1,km
              KAPPA_THIC(:,:,kk_sub,kk,bid) = KAPPA_LATERAL(:,:,bid)  &
                                        * KAPPA_VERTICAL(:,:,kk,bid)
            enddo
          enddo
        endif


!-----------------------------------------------------------------------
!
!     control slope of isopycnal surfaces or KAPPA
!
!-----------------------------------------------------------------------

        if (kmin_type == kmin_type_read ) then !store MINOR diffusivity in the XX components 
           KXX_ISOP(:,:,:,:,bid) = KAPPA_ISOP(:,:,:,:,bid) !SJR ADD !storage for KAPPA_ISOP for overall tapering factor
           KXX_THIC(:,:,:,:,bid) = KAPPA_THIC(:,:,:,:,bid) !SJR ADD !storage for KAPPA_THIC for overall tapering factor
           HXX_DIFF(:,:,:,:,bid) =   HOR_DIFF(:,:,:,:,bid) !SJR ADD !storage for   HOR_DIFF for overall tapering factor
        endif

        do kk=1,km

          kp1 = min(kk+1,km)
          reference_depth(ktp) = zt(kp1)
          reference_depth(kbt) = zw(kp1)
          if ( kk == km )  reference_depth(ktp) = zw(kp1)

          do kk_sub = ktp,kbt 

            kid = kk + kk_sub - 2

!-----------------------------------------------------------------------
!
!     control KAPPA to reduce the isopycnal mixing near the
!     ocean surface Large et al (1997), JPO, 27, pp 2418-2447.
!     WORK1 = ratio between the depth of water parcel and
!     the vertical displacement of isopycnal surfaces
!     where the vertical displacement =
!     Rossby radius * slope of isopycnal surfaces
!
!-----------------------------------------------------------------------

            if ( transition_layer_on ) then
              SLA = SLA_SAVE(:,:,kk_sub,kk,bid)
            else
              if ( partial_bottom_cells ) then
                 SLA = DZTW(:,:,kid,bid)*sqrt(p5*(                               &
                        (SLX(:,:,1,kk_sub,kk,bid)**2                    & 
                       + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                      + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                      + eps
              else
                 SLA = dzw(kid)*sqrt(p5*(                               &
                        (SLX(:,:,1,kk_sub,kk,bid)**2                    & 
                       + SLX(:,:,2,kk_sub,kk,bid)**2)/DXT(:,:,bid)**2   &
                      + (SLY(:,:,1,kk_sub,kk,bid)**2                    &
                       + SLY(:,:,2,kk_sub,kk,bid)**2)/DYT(:,:,bid)**2)) &
                      + eps
              endif
            endif

            TAPER1 = c1 
            if ( .not. transition_layer_on ) then

              if ( kk == 1 ) then
                dz_bottom = c0
              else
                dz_bottom = zt(kk-1)
              endif

              if (slope_control == slope_control_tanh) then

                WORK1 = min(c1,zt(kk)*RBR(:,:,bid)/SLA)
                TAPER1 = p5*(c1+sin(pi*(WORK1-p5)))

!     use the Rossby deformation radius tapering
!     only within the boundary layer

                TAPER1 = merge(TAPER1, c1,  &
                               dz_bottom <= BL_DEPTH(:,:,bid))

              else

!     sine function is replaced by
!     function = 4.*x*(1.-abs(x)) for |x|<0.5

                WORK1 = min(c1,zt(kk)*RBR(:,:,bid)/SLA)
                TAPER1 = (p5+c2*(WORK1-p5)*(c1-abs(WORK1-p5)))

                TAPER1 = merge(TAPER1, c1,  &
                               dz_bottom <= BL_DEPTH(:,:,bid))

              endif

            endif

!-----------------------------------------------------------------------
!
!     control KAPPA for numerical stability
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     methods to control slope
!
!-----------------------------------------------------------------------

            TAPER2 = c1 
            TAPER3 = c1 

            select case (slope_control)
            case (slope_control_tanh)

!     method by Danabasoglu & Mcwilliams (1995)

              TAPER2 = merge(p5*  &
                 (c1-tanh(c10*SLA/slm_r-c4)), c0, SLA < slm_r)

              if ( diff_tapering ) then
                TAPER3 = merge(p5*  &
                 (c1-tanh(c10*SLA/slm_b-c4)), c0, SLA < slm_b)
              else
                TAPER3 = TAPER2
              endif

            case (slope_control_notanh)

!     similar to DM95 except replacing tanh by
!     function = x*(1.-0.25*abs(x)) for |x|<2
!              = sign(x)            for |x|>2
!     (faster than DM95)

              do j=1,ny_block
                do i=1,nx_block
                  if (SLA(i,j) > 0.2_r8*slm_r .and. &
                      SLA(i,j) < 0.6_r8*slm_r) then
                    TAPER2(i,j) = &
                       p5*(c1-(2.5_r8*SLA(i,j)/slm_r-c1)*  &
                          (c4-abs(c10*SLA(i,j)/slm_r-c4)))
                  else if (SLA(i,j) >= 0.6_r8*slm_r) then
                    TAPER2(i,j) = c0
                  endif
                enddo
              enddo

              if ( diff_tapering ) then
                do j=1,ny_block
                  do i=1,nx_block
                    if (SLA(i,j) > 0.2_r8*slm_b .and. &
                        SLA(i,j) < 0.6_r8*slm_b) then
                      TAPER3(i,j) = &
                         p5*(c1-(2.5_r8*SLA(i,j)/slm_b-c1)* &
                            (c4-abs(c10*SLA(i,j)/slm_b-c4)))
                    else if (SLA(i,j) >= 0.6_r8*slm_b) then
                      TAPER3(i,j) = c0
                    endif
                  end do
                end do
              else
                TAPER3 = TAPER2
              endif

            case (slope_control_clip)

!     slope clipping

             if ( partial_bottom_cells ) then
              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLX(i,j,n,kk_sub,kk,bid)  &
                          * DZTW(i,j,kid,bid) / HUS(i,j,bid)) > slm_r) then
                      SLX(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUS(i,j,bid)  &
                                     / DZTW(i,j,kid,bid),            &
                                  SLX(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo

              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLY(i,j,n,kk_sub,kk,bid)  &
                          * DZTW(i,j,kid,bid) / HUW(i,j,bid)) > slm_r) then  
                      SLY(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUW(i,j,bid)  &
                                     / DZTW(i,j,kid,bid),            &
                                  SLY(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo
             else
              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLX(i,j,n,kk_sub,kk,bid)  &
                          * dzw(kid) / HUS(i,j,bid)) > slm_r) then
                      SLX(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUS(i,j,bid)  &
                                     * dzwr(kid),            &
                                  SLX(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo

              do n=1,2
                do j=1,ny_block
                  do i=1,nx_block
                    if (abs(SLY(i,j,n,kk_sub,kk,bid)  &
                          * dzw(kid) / HUW(i,j,bid)) > slm_r) then  
                      SLY(i,j,n,kk_sub,kk,bid) =             &
                                  sign(slm_r * HUW(i,j,bid)  &
                                     * dzwr(kid),            &
                                  SLY(i,j,n,kk_sub,kk,bid))
                    endif
                  enddo
                enddo
              enddo

             endif

            case (slope_control_Gerd)

!     method by Gerdes et al (1991)

              do j=1,ny_block
                do i=1,nx_block
                  if (SLA(i,j) > slm_r)  &
                    TAPER2(i,j) = (slm_r/SLA(i,j))**2
                enddo
              enddo

              if (diff_tapering) then
                do j=1,ny_block
                  do i=1,nx_block
                    if (SLA(i,j) > slm_b)  &
                      TAPER3(i,j) = (slm_b/SLA(i,j))**2
                  enddo
                enddo
              else
                TAPER3 = TAPER2
              endif

            end select

            if ( transition_layer_on ) then
              TAPER2 = merge(c1, TAPER2, reference_depth(kk_sub) &
                                         <= TLT%DIABATIC_DEPTH(:,:,bid))
              TAPER3 = merge(c1, TAPER3, reference_depth(kk_sub) &
                                         <= TLT%DIABATIC_DEPTH(:,:,bid))
            endif

            if ( transition_layer_on  .and.  use_const_ah_bkg_srfbl ) then

              HOR_DIFF(:,:,kk_sub,kk,bid) = ah_bkg_srfbl

            else if ( transition_layer_on .and.               &
                    ( .not. use_const_ah_bkg_srfbl      .or.  &
                      kappa_isop_type == kappa_type_eg  .or.  &  !SJR COM - THE EG CHECKS DO NOTHING HERE
                      kappa_thic_type == kappa_type_eg ) ) then  !SJR COM - IF transition_layer_on=T, 
                                                                 !SJR COM - the two checks on use_const_ah_bkg_srfl are sufficient
                                                                 !SJR COM - IF transition_layer_on=F, it will pass both to 'else' 

              HOR_DIFF(:,:,kk_sub,kk,bid) = KAPPA_ISOP(:,:,kk_sub,kk,bid) 

            else

              if ( .not. ( kk == 1 .and. kk_sub == ktp ) ) then 

                if ( use_const_ah_bkg_srfbl ) then 
                  HOR_DIFF(:,:,kk_sub,kk,bid) =                       &
                       merge( ah_bkg_srfbl * (c1 - TAPER1 * TAPER2)   &
                              * KAPPA_VERTICAL(:,:,kk,bid),           &
                              c0, dz_bottom <= BL_DEPTH(:,:,bid) )
                else
                  HOR_DIFF(:,:,kk_sub,kk,bid) =                &
                         merge( KAPPA_ISOP(:,:,kk_sub,kk,bid)  &
                               * (c1 - TAPER1 * TAPER2),       &
                              c0, dz_bottom <= BL_DEPTH(:,:,bid) )
                endif

              endif

            endif

            KAPPA_ISOP(:,:,kk_sub,kk,bid) =  &
                  TAPER1 * TAPER2 * KAPPA_ISOP(:,:,kk_sub,kk,bid)
            KAPPA_THIC(:,:,kk_sub,kk,bid) =  &
                  TAPER1 * TAPER3 * KAPPA_THIC(:,:,kk_sub,kk,bid)

          end do  ! end of kk_sub loop


!-----------------------------------------------------------------------
!
!     impose the boundary conditions by setting KAPPA=0
!     in the quarter cells adjacent to rigid boundaries.
!     bottom B.C.s are considered within the kk-loop.
!
!-----------------------------------------------------------------------

!     B.C. at the bottom

          where (kk == KMT(:,:,bid))
            KAPPA_ISOP(:,:,kbt,kk,bid) = c0
            KAPPA_THIC(:,:,kbt,kk,bid) = c0
          end where

        enddo              ! end of kk-loop

!     B.C. at the top

        KAPPA_ISOP(:,:,ktp,1,bid) = c0
        KAPPA_THIC(:,:,ktp,1,bid) = c0

        if ( transition_layer_on ) call apply_vertical_profile_to_isop_hor_diff ( this_block ) 

        FZTOP(:,:,:,bid) = c0        ! zero flux B.C. at the surface
     
      call timer_stop(timer_gm2, block_id=this_block%local_id)
      call timer_start(timer_gm3, block_id=this_block%local_id)

!-------------v--SJR ADD--v--------------------------------------------- 
      if (kdir_type == kdir_type_shrd .or. krat_type == krat_type_shrd) then
         call shear_dispersion (UMIX, VMIX, TMIX, this_block) !sets NX_SHRD, NY_SHRD, KRAT_SHRD
      endif
!_________________________________________________
           if (kdir_type == kdir_type_shrd) then

              NX_ISOP(:,:,:,bid) = NX_SHRD(:,:,:,bid)
              NY_ISOP(:,:,:,bid) = NY_SHRD(:,:,:,bid)

           elseif (kdir_type == kdir_type_east) then

              GRADX_FCOR = c0
              GRADY_FCOR = c0
              GRADX_FCOR(2:nx_block,2:ny_block) = DXTR(2:nx_block  ,2:ny_block  ,bid) * & !x-gradient of coriolis parameter at T-pts
                                           p5 * ( FCOR(2:nx_block  ,2:ny_block  ,bid) + &
                                                  FCOR(2:nx_block  ,1:ny_block-1,bid) - & 
                                                  FCOR(1:nx_block-1,2:ny_block  ,bid) - &
                                                  FCOR(1:nx_block-1,1:ny_block-1,bid) ) 
              GRADY_FCOR(2:nx_block,2:ny_block) = DYTR(2:nx_block  ,2:ny_block  ,bid) * & !y-gradient of coriolis parameter at T-pts
                                           p5 * ( FCOR(2:nx_block  ,2:ny_block  ,bid) + &
                                                  FCOR(1:nx_block-1,2:ny_block  ,bid) - & 
                                                  FCOR(2:nx_block  ,1:ny_block-1,bid) - &
                                                  FCOR(1:nx_block-1,1:ny_block-1,bid) ) 

              call ugrid_to_tgrid(WORK1, GRADY_FCOR,bid)
              call ugrid_to_tgrid(WORK2,-GRADX_FCOR,bid)
              do kk=1,km
                 NX_ISOP(:,:,kk,bid) = WORK1 
                 NY_ISOP(:,:,kk,bid) = WORK2 
              enddo

!_________________________________________________
           elseif (kdir_type == kdir_type_zonl) then
              
              call ugrid_to_tgrid(WORK1,cos(ANGLE(:,:,bid)),bid)
              call ugrid_to_tgrid(WORK2,sin(ANGLE(:,:,bid)),bid)
              do kk=1,km
                 NX_ISOP(:,:,kk,bid) = WORK1 
                 NY_ISOP(:,:,kk,bid) = WORK2 
              enddo

!_________________________________________________
           elseif (kdir_type == kdir_type_flow) then

              do kk=1,km
                 call ugrid_to_tgrid(NX_ISOP(:,:,kk,bid),UMIX(:,:,kk),bid)
                 call ugrid_to_tgrid(NY_ISOP(:,:,kk,bid),VMIX(:,:,kk),bid) !Velocity... cross stream/along stream
              enddo

!_________________________________________________
           elseif (kdir_type == kdir_type_apvg) then

              do kk=1,km

                 call state (kk,kk,TMIX(:,:,kk,1),TMIX(:,:,kk,2),this_block, RHOOUT=WORK1) !RHO at kk

                 if (kk == 1 ) then
                    WORK3 = WORK1
                 else
                    call state(kk-1,kk,TMIX(:,:,kk-1,1),TMIX(:,:,kk-1,2),this_block, RHOOUT=WORK3) !RHO at kk-1
                    WORK3 = p5*(WORK3 + WORK1)
                 endif

                 if (kk == km ) then
                    WORK4 = WORK1
                 else
                    call state(kk+1,kk,TMIX(:,:,kk+1,1),TMIX(:,:,kk+1,2),this_block, RHOOUT=WORK4) !RHO at kk+1
                    where ( kk /= KMT(:,:,bid) ) 
                       WORK4 = p5*(WORK4 + WORK1)
                    elsewhere
                       WORK4 = WORK1
                    endwhere
                 endif

                 kp1=min(kk+1,km)
                 if ( partial_bottom_cells ) then
                    WORK1 = (   WORK3*DZTW(:,:,kp1,bid) / DZTW(:,:,kk,bid) &
                              - WORK4*DZTW(:,:,kk,bid)  / DZTW(:,:,kp1,bid) ) / &
                                  (   DZTW(:,:,kp1,bid) + DZTW(:,:,kk,bid) )  !Dz(rho) at T-pt
                 else
                    WORK1 = (   WORK3*dzw(kp1) / dzw(kk) &
                              - WORK4*dzw(kk)  / dzw(kp1) ) / &
                                  (   dzw(kp1) + dzw(kk) )  !Dz(rho) at T-pt
                 endif

                 where ( kk > KMT(:,:,bid) ) 
                    WORK1 = c0
                 endwhere

                 call zcurl(kk,WORK3,UMIX(:,:,kk),VMIX(:,:,kk),this_block) !vorticity at T-pt
                 WORK2 = WORK1*(WORK3*TAREA_R(:,:,bid) + FCORT(:,:,bid)) ! PV = pot vorticity at T-pt
                 WORK2 = WORK2 / max(maxval(abs(WORK2)),eps2)

                 call grad(kk,WORK3,WORK4,WORK2,this_block) !grad(PV) at U-pt

                 globval = maxval(sqrt(WORK3**2+WORK4**2))
                 WORK3 = WORK3 / max(globval,eps2)
                 WORK4 = WORK4 / max(globval,eps2)

                 call ugrid_to_tgrid(NX_ISOP(:,:,kk,bid), WORK4,bid)
                 call ugrid_to_tgrid(NY_ISOP(:,:,kk,bid),-WORK3,bid) !perpendicular to grad(PV)

              enddo

           endif
!_________________________________________________

        do kk=1,km
           WORK1=sqrt(NX_ISOP(:,:,kk,bid)**2+NY_ISOP(:,:,kk,bid)**2)
           where (WORK1 < eps2)
              NY_ISOP(:,:,kk,bid) = c0
              NX_ISOP(:,:,kk,bid) = c1
           elsewhere
              NY_ISOP(:,:,kk,bid)=NY_ISOP(:,:,kk,bid)/WORK1
              NX_ISOP(:,:,kk,bid)=NX_ISOP(:,:,kk,bid)/WORK1
           endwhere 
        enddo

        if (addrandfluc) then
           do kk=1,km
              call random_seed
              call random_number(WORK1)
              call random_seed
              call random_number(WORK2)
              where (abs(NX_ISOP(:,:,kk,bid)) < eps2)
                 WORK3 = pi/c2 + pi/c4*sqrt(-c2*log(WORK1))*cos(pi2*WORK2) 
              elsewhere
                 WORK3 = atan(NY_ISOP(:,:,kk,bid) / NX_ISOP(:,:,kk,bid)) + pi/c4*sqrt(-c2*log(WORK1))*cos(pi2*WORK2) 
              endwhere 
              NX_ISOP(:,:,kk,bid) = cos(WORK3)
              NY_ISOP(:,:,kk,bid) = sin(WORK3)
           enddo
        end if

        if (kmin_type == kmin_type_read ) then !store MINOR diffusivity in the XX components 
           do kk_sub=ktp,kbt
              where (ABS(KXX_ISOP(:,:,kk_sub,:,bid)) > eps)
                 KXX_ISOP(:,:,kk_sub,:,bid) = KAPPA_ISOP(:,:,kk_sub,:,bid) / &
                                                KXX_ISOP(:,:,kk_sub,:,bid) * &
                                          MINOR_EIGENVAL(:,:,       :,bid) * minorfactor !KAPPA_ISOP overall tapering factor
              endwhere
              where (ABS(KXX_THIC(:,:,kk_sub,:,bid)) > eps)
                 KXX_THIC(:,:,kk_sub,:,bid) = KAPPA_THIC(:,:,kk_sub,:,bid) / &
                                                KXX_THIC(:,:,kk_sub,:,bid) * &
                                          MINOR_EIGENVAL(:,:,       :,bid) * minorfactor !KAPPA_THIC overall tapering factor
              endwhere
              where (ABS(HXX_DIFF(:,:,kk_sub,:,bid)) > eps)
                 HXX_DIFF(:,:,kk_sub,:,bid) =   HOR_DIFF(:,:,kk_sub,:,bid) / &
                                                HXX_DIFF(:,:,kk_sub,:,bid) * &
                                          MINOR_EIGENVAL(:,:,       :,bid) * minorfactor !  HOR_DIFF overall tapering factor
              endwhere
           enddo
        else
           KXX_ISOP(:,:,:,:,bid) = KAPPA_ISOP(:,:,:,:,bid) * minorfactor
           KXX_THIC(:,:,:,:,bid) = KAPPA_THIC(:,:,:,:,bid) * minorfactor
           HXX_DIFF(:,:,:,:,bid) =   HOR_DIFF(:,:,:,:,bid) * minorfactor
        endif

        if ( krat_type == krat_type_read ) then
           K_EIGENVAL_RAT(:,:,1,:,bid) = MAJOR_EIGENVAL(:,:,:,bid) / MINOR_EIGENVAL(:,:,:,bid)
           K_EIGENVAL_RAT(:,:,2,:,bid) = K_EIGENVAL_RAT(:,:,1,:,bid) 
        elseif ( krat_type == krat_type_shrd ) then
           K_EIGENVAL_RAT(:,:,:,:,bid) = KRAT_SHRD(:,:,:,:,bid)
        else
           K_EIGENVAL_RAT(:,:,:,:,bid) = erat_const
        endif

        if (cflmult > c0) then
         do kk=1,km
          do kk_sub=1,2
           WORK1 = ABS(KXX_ISOP(:,:,kk_sub,kk,bid)) + ABS(HXX_DIFF(:,:,kk_sub,kk,bid))
           if (cflmajoronly) then
              WORK3 = ( cflmult / dtt - &
                           ( ABS( WORK1 * NY_ISOP(:,:,kk,bid) * DXTR(:,:,bid)**2 ) + &
                             ABS( WORK1 * NX_ISOP(:,:,kk,bid) * DYTR(:,:,bid)**2 ) ) ) / &
                           ( ABS( WORK1 * NX_ISOP(:,:,kk,bid) * DXTR(:,:,bid)**2 ) + &
                             ABS( WORK1 * NY_ISOP(:,:,kk,bid) * DYTR(:,:,bid)**2 ) ) !r*
              K_EIGENVAL_RAT(:,:,kk_sub,kk,bid) = MIN(K_EIGENVAL_RAT(:,:,kk_sub,kk,bid),MAX(WORK3,c1))
           endif
           WORK2 = WORK1 * K_EIGENVAL_RAT(:,:,kk_sub,kk,bid) 
           WORK4 = dtt / cflmult * &
                         ( ( ABS( WORK2 * NX_ISOP(:,:,kk,bid) ) + &
                             ABS( WORK1 * NY_ISOP(:,:,kk,bid) ) ) * DXTR(:,:,bid)**2 + &
                           ( ABS( WORK2 * NY_ISOP(:,:,kk,bid) ) + &
                             ABS( WORK1 * NX_ISOP(:,:,kk,bid) ) ) * DYTR(:,:,bid)**2 )
           if (cflmajoronly) then
              where ( WORK4 > c1 ) 
                 K_EIGENVAL_RAT(:,:,kk_sub,kk,bid) = c1
                 KXX_ISOP(:,:,kk_sub,kk,bid) =  KXX_ISOP(:,:,kk_sub,kk,bid) / WORK4
                 HXX_DIFF(:,:,kk_sub,kk,bid) =  HXX_DIFF(:,:,kk_sub,kk,bid) / WORK4
                 KXX_THIC(:,:,kk_sub,kk,bid) =  KXX_THIC(:,:,kk_sub,kk,bid) / WORK4
              endwhere
           else
              where ( WORK4 > c1 ) 
                 KXX_ISOP(:,:,kk_sub,kk,bid) =  KXX_ISOP(:,:,kk_sub,kk,bid) / WORK4
                 HXX_DIFF(:,:,kk_sub,kk,bid) =  HXX_DIFF(:,:,kk_sub,kk,bid) / WORK4
                 KXX_THIC(:,:,kk_sub,kk,bid) =  KXX_THIC(:,:,kk_sub,kk,bid) / WORK4
              endwhere
           endif
          enddo
         enddo
        endif

        do kk=1,km
           do kk_sub=ktp,kbt
              WORK1 = NY_ISOP(:,:,kk,bid)**2 + NX_ISOP(:,:,kk,bid)**2 * K_EIGENVAL_RAT(:,:,kk_sub,kk,bid)
              WORK2 = NX_ISOP(:,:,kk,bid)**2 + NY_ISOP(:,:,kk,bid)**2 * K_EIGENVAL_RAT(:,:,kk_sub,kk,bid)
              WORK3 = NX_ISOP(:,:,kk,bid)    * NY_ISOP(:,:,kk,bid)    *(K_EIGENVAL_RAT(:,:,kk_sub,kk,bid)-c1)   

              KYY_THIC(:,:,kk_sub,kk,bid)=KXX_THIC(:,:,kk_sub,kk,bid)*WORK2
              KXY_THIC(:,:,kk_sub,kk,bid)=KXX_THIC(:,:,kk_sub,kk,bid)*WORK3
              KXX_THIC(:,:,kk_sub,kk,bid)=KXX_THIC(:,:,kk_sub,kk,bid)*WORK1
              KYY_ISOP(:,:,kk_sub,kk,bid)=KXX_ISOP(:,:,kk_sub,kk,bid)*WORK2
              KXY_ISOP(:,:,kk_sub,kk,bid)=KXX_ISOP(:,:,kk_sub,kk,bid)*WORK3
              KXX_ISOP(:,:,kk_sub,kk,bid)=KXX_ISOP(:,:,kk_sub,kk,bid)*WORK1
              HYY_DIFF(:,:,kk_sub,kk,bid)=HXX_DIFF(:,:,kk_sub,kk,bid)*WORK2
              HXY_DIFF(:,:,kk_sub,kk,bid)=HXX_DIFF(:,:,kk_sub,kk,bid)*WORK3
              HXX_DIFF(:,:,kk_sub,kk,bid)=HXX_DIFF(:,:,kk_sub,kk,bid)*WORK1
           enddo        
        enddo
        if ( ah_bkg_bottom /= c0 ) then
          where ( k == KMT(:,:,bid) ) 
!-------------v--SJR MOD--v--------------------------------------------- 
            HXX_DIFF(:,:,kbt,k,bid) = ah_bkg_bottom *                         &
                            ( NY_ISOP(:,:,k,bid)**2 + NX_ISOP(:,:,k,bid)**2 * &
                              K_EIGENVAL_RAT(:,:,kbt,k,bid) )
            HYY_DIFF(:,:,kbt,k,bid) = ah_bkg_bottom *                         &
                            ( NX_ISOP(:,:,k,bid)**2 + NY_ISOP(:,:,k,bid)**2 * &
                              K_EIGENVAL_RAT(:,:,kbt,k,bid) )
            HXY_DIFF(:,:,kbt,k,bid) = ah_bkg_bottom *                         &
                            NX_ISOP(:,:,k,bid)    * NY_ISOP(:,:,k,bid)    *   &
                            (K_EIGENVAL_RAT(:,:,kbt,k,bid)-c1)
!-------------^--SJR MOD--^--------------------------------------------- 
          endwhere
        endif
      call timer_stop(timer_gm3, block_id=this_block%local_id)
      call timer_start(timer_gm4, block_id=this_block%local_id)

        do kk=1,km
           do kk_sub=1,2
              kk_shft = kk-2+kk_sub
              do i_sub=1,2
                 istrt = 3-i_sub
                 jstrt = 3-i_sub
                 istop = nx_block-i_sub+1
                 jstop = ny_block-i_sub+1
                 if ( partial_bottom_cells ) then
                    DD_SLX(2:nx_block,:,i_sub,kk_sub,kk,bid) = &
                      -SLX(2:nx_block,:,i_sub,kk_sub,kk,bid) * &
                       DZTW(2:nx_block,:,kk_shft,bid) / HUS(istrt:istop,:,bid)
                    DD_SLY(:,2:ny_block,i_sub,kk_sub,kk,bid) = &
                      -SLY(:,2:ny_block,i_sub,kk_sub,kk,bid) * &
                       DZTW(:,2:ny_block,kk_shft,bid) / HUW(:,jstrt:jstop,bid)
                 else
                    DD_SLX(2:nx_block,:,i_sub,kk_sub,kk,bid) = &
                      -SLX(2:nx_block,:,i_sub,kk_sub,kk,bid) * &
                       dzw(kk_shft) / HUS(istrt:istop,:,bid)
                    DD_SLY(:,2:ny_block,i_sub,kk_sub,kk,bid) = &
                      -SLY(:,2:ny_block,i_sub,kk_sub,kk,bid) * &
                       dzw(kk_shft) / HUW(:,jstrt:jstop,bid)
                 endif
              enddo
           enddo
        enddo
        DD_SLX(1,:,:,:,:,bid) = DD_SLX(2,:,:,:,:,bid)
        DD_SLY(:,1,:,:,:,bid) = DD_SLY(:,2,:,:,:,bid)

        do j_sub=1,2
           do i_sub=1,2
              SL_MAG(:,:,i_sub,j_sub,:,:,bid) = sqrt( &
              DD_SLX(:,:,i_sub,      :,:,bid)**2 + &
              DD_SLY(:,:,      j_sub,:,:,bid)**2     )
           enddo
        enddo

!-------------^--SJR ADD--^--------------------------------------------- 

        if ( transition_layer_on ) then

          call merged_streamfunction ( this_block )

        else

          TLT%DIABATIC_DEPTH(:,:,bid) = c0
          TLT%THICKNESS(:,:,bid)      = c0
          TLT%INTERIOR_DEPTH(:,:,bid) = c0

           do kk=1,km
            do kk_sub=ktp,kbt
!-------------v--SJR MOD--v--------------------------------------------- 
              do i_sub=1,2
                where ( kk <= KMT(:,:,bid) ) 
                   SF_SLXX(:,:,i_sub,kk_sub,kk,bid) =          &
                           KXX_THIC(:,:,kk_sub,kk,bid)     &
                         * DD_SLX(:,:,i_sub,kk_sub,kk,bid)

                   SF_SLYY(:,:,i_sub,kk_sub,kk,bid) =          &
                           KYY_THIC(:,:,kk_sub,kk,bid)     &
                         * DD_SLY(:,:,i_sub,kk_sub,kk,bid)
                endwhere
              enddo
!-------------^--SJR MOD--^--------------------------------------------- 
            enddo  ! end of kk_sub-loop
           enddo    ! end of kk-loop
          if (.not. isoonly) then
           do kk=1,km
            do kk_sub=ktp,kbt
!-------------v--SJR MOD--v--------------------------------------------- 
              do i_sub=1,2
                where ( kk <= KMT(:,:,bid) ) 
                   SF_SLXY(:,:,i_sub,kk_sub,kk,bid) =          &
                           KXY_THIC(:,:,kk_sub,kk,bid)     &
                         * DD_SLX(:,:,i_sub,kk_sub,kk,bid)

                   SF_SLYX(:,:,i_sub,kk_sub,kk,bid) =          &
                           KXY_THIC(:,:,kk_sub,kk,bid)     &
                         * DD_SLY(:,:,i_sub,kk_sub,kk,bid)
                endwhere
              enddo
!-------------^--SJR MOD--^--------------------------------------------- 
            enddo  ! end of kk_sub-loop
           enddo    ! end of kk-loop
          endif
        endif

      call timer_stop(timer_gm4, block_id=this_block%local_id)
      call timer_stop(timer_gm1, block_id=this_block%local_id)
      endif  ! end of k==1 if statement

      KMASK = merge(c1, c0, k < KMT(:,:,bid))
      

!-----------------------------------------------------------------------
!
!     calculate effective vertical diffusion coefficient
!     NOTE: it is assumed that VDC has been set before this
!           in vmix_coeffs or something similar.
!
!     Dz(VDC * Dz(T)) where D is derivative rather than difference
!     VDC = (Az(dz*Ax(KAPPA*HYX*SLX**2)) + Az(dz*Ay(KAPPA*HXY*SLY**2)))*
!           dzw/TAREA
!
!-----------------------------------------------------------------------

      
      FZBOT    = c0 !SJR ADD

      if ( k < km ) then

!-------------v--SJR MOD--v--------------------------------------------- 
      call timer_start(timer_gm5, block_id=this_block%local_id)
        WORK1 = c0
        if (isoonly) then
           do kk=1,2  !ktp(k+1),kbt(k)
              kk_shft = k+2-kk
              WORK2 = c0
              do j_sub=1,2  !jnorth,jsouth  
                 do i_sub=1,2  !ieast,iwest
                    WORK2 = WORK2 +  SUBCELLV(:,:,i_sub,j_sub,kk,kk_shft,bid) * &
                                     KXX_ISOP(:,:,            kk,kk_shft,bid) * &
                                     ( DD_SLX(:,:,i_sub,      kk,kk_shft,bid)**2 + &
                                       DD_SLY(:,:,      j_sub,kk,kk_shft,bid)**2 )
                 enddo
              enddo
              if ( partial_bottom_cells ) WORK2 = WORK2 * DZT(:,:,kk_shft,bid) / dz(kk_shft)
              WORK1 = WORK1 + WORK2
           enddo
        else
           do kk=1,2  !ktp(k+1),kbt(k)
              kk_shft = k+2-kk
              WORK2 = c0
              do j_sub=1,2  !jnorth,jsouth  
                 do i_sub=1,2  !ieast,iwest
                    WORK2 = WORK2 +  SUBCELLV(:,:,i_sub,j_sub,kk,kk_shft,bid) * &
                                (    KXX_ISOP(:,:,            kk,kk_shft,bid) * &
                                       DD_SLX(:,:,i_sub,      kk,kk_shft,bid)**2 + &
                                     KYY_ISOP(:,:,            kk,kk_shft,bid) * &
                                       DD_SLY(:,:,      j_sub,kk,kk_shft,bid)**2 + &
                                c2 * KXY_ISOP(:,:,            kk,kk_shft,bid) * &
                                       DD_SLX(:,:,i_sub,      kk,kk_shft,bid) * &
                                       DD_SLY(:,:,      j_sub,kk,kk_shft,bid) )
                 enddo
              enddo
              if ( partial_bottom_cells ) WORK2 = WORK2 * DZT(:,:,kk_shft,bid) / dz(kk_shft)
              WORK1 = WORK1 + WORK2
           enddo
        endif
        if ( partial_bottom_cells ) then
           WORK1=WORK1*KMASK*TAREA_R(:,:,bid)/DZTW(:,:,k,bid)  
        else
           WORK1=WORK1*KMASK*TAREA_R(:,:,bid)*dzwr(k)  
        endif

        if ( .not. vertdiffhere ) then
           VDC_GM(:,:,k,bid) = WORK1
           do n=1,size(VDC,DIM=4)
              VDC(:,:,k,n,bid) = VDC(:,:,k,n,bid) + WORK1
           end do
        end if
!-------------^--SJR MOD--^--------------------------------------------- 
      call timer_stop(timer_gm5, block_id=this_block%local_id)

      end if

!-----------------------------------------------------------------------
!
!     check if some horizontal diffusion needs to be added to the
!     bottom half of the bottom cell
!
!-----------------------------------------------------------------------

      call timer_start(timer_gm6, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!
!     combine isopycnal and horizontal diffusion coefficients
!
!-----------------------------------------------------------------------
      
!-------------v--SJR MOD--v--------------------------------------------- 

            istop = nx_block-1
            jstop = ny_block-1

            WORK1(1:nx_block-1,:) = &
                   ( SUBCELLV(1:istop,   :,1,1,1,k,bid) + &
                     SUBCELLV(1:istop,   :,1,2,1,k,bid) ) * &
                   ( KXX_ISOP(1:istop,   :,    1,k,bid) + &
                     HXX_DIFF(1:istop,   :,    1,k,bid) ) + &
                   ( SUBCELLV(1:istop,   :,1,1,2,k,bid) + &
                     SUBCELLV(1:istop,   :,1,2,2,k,bid) ) * &
                   ( KXX_ISOP(1:istop,   :,    2,k,bid) + &
                     HXX_DIFF(1:istop,   :,    2,k,bid) ) + &
                   ( SUBCELLV(2:nx_block,:,2,1,1,k,bid) + &
                     SUBCELLV(2:nx_block,:,2,2,1,k,bid) ) * &
                   ( KXX_ISOP(2:nx_block,:,    1,k,bid) + &
                     HXX_DIFF(2:nx_block,:,    1,k,bid) ) + &
                   ( SUBCELLV(2:nx_block,:,2,1,2,k,bid) + &
                     SUBCELLV(2:nx_block,:,2,2,2,k,bid) ) * &
                   ( KXX_ISOP(2:nx_block,:,    2,k,bid) + &
                     HXX_DIFF(2:nx_block,:,    2,k,bid) ) 

            WORK2(:, 1:ny_block-1) = &
                   ( SUBCELLV(:,1:jstop,   1,1,1,k,bid) + &
                     SUBCELLV(:,1:jstop,   2,1,1,k,bid) ) * &
                   ( KYY_ISOP(:,1:jstop,       1,k,bid) + &
                     HYY_DIFF(:,1:jstop,       1,k,bid) ) + &
                   ( SUBCELLV(:,1:jstop,   1,1,2,k,bid) + &
                     SUBCELLV(:,1:jstop,   2,1,2,k,bid) ) * &
                   ( KYY_ISOP(:,1:jstop,       2,k,bid) + &
                     HYY_DIFF(:,1:jstop,       2,k,bid) ) + &
                   ( SUBCELLV(:,2:ny_block,1,2,1,k,bid) + &
                     SUBCELLV(:,2:ny_block,2,2,1,k,bid) ) * &
                   ( KYY_ISOP(:,2:ny_block,    1,k,bid) + &
                     HYY_DIFF(:,2:ny_block,    1,k,bid) ) + &
                   ( SUBCELLV(:,2:ny_block,1,2,2,k,bid) + &
                     SUBCELLV(:,2:ny_block,2,2,2,k,bid) ) * &
                   ( KYY_ISOP(:,2:ny_block,    2,k,bid) + &
                     HYY_DIFF(:,2:ny_block,    2,k,bid) ) 

!-------------^--SJR MOD--^--------------------------------------------- 

!-----------------------------------------------------------------------
!
!     start loop over tracers
!
!-----------------------------------------------------------------------

      call timer_start(timer_nloop, block_id=this_block%local_id)

      do n = 1,nt

!-----------------------------------------------------------------------
!
!     calculate horizontal fluxes thru vertical faces of T-cell
!     FX = dz*HYX*Ax(Az(KAPPA))*Dx(T) : flux in x-direction
!     FY = dz*HXY*Ay(Az(KAPPA))*Dy(T) : flux in y-direction
!
!-----------------------------------------------------------------------

!-------------v--SJR MOD--v--------------------------------------------- 
         FX(:,:,n) = TX(:,:,k,n,bid) / HUS(:,:,bid) * WORK1
         FY(:,:,n) = TY(:,:,k,n,bid) / HUW(:,:,bid) * WORK2 
!-------------^--SJR MOD--^--------------------------------------------- 

      end do
      call timer_stop(timer_gm6, block_id=this_block%local_id)
      call timer_start(timer_gm7, block_id=this_block%local_id)

!-------------v--SJR ADD--v--------------------------------------------- 
      if (.not. isoonly) then
         do j_sub=1,2  
            jstrt = j_sub
            jstop = ny_block-2+j_sub
            astrt = 3-j_sub
            astop = ny_block+1-j_sub
            do i_sub=1,2  
               istrt = i_sub
               istop = nx_block-2+i_sub
               bstrt = 3-i_sub
               bstop = nx_block+1-i_sub
   
               WORK3(1:nx_block-1,2:ny_block) = & !FX KXY term
                 (   SUBCELLV(istrt:istop,2:ny_block,i_sub,j_sub,1,k,bid) * &
                   ( KXY_ISOP(istrt:istop,2:ny_block,            1,k,bid) + &
                     HXY_DIFF(istrt:istop,2:ny_block,            1,k,bid) ) + &
                     SUBCELLV(istrt:istop,2:ny_block,i_sub,j_sub,2,k,bid) * &
                   ( KXY_ISOP(istrt:istop,2:ny_block,            2,k,bid) + &
                     HXY_DIFF(istrt:istop,2:ny_block,            2,k,bid) ) ) / &
                             HUW(istrt:istop,astrt:astop,bid)

               WORK4(2:nx_block,1:ny_block-1) = & !FY KXY term
                (   SUBCELLV(2:nx_block,jstrt:jstop,i_sub,j_sub,1,k,bid) * &
                  ( KXY_ISOP(2:nx_block,jstrt:jstop,            1,k,bid) + &
                    HXY_DIFF(2:nx_block,jstrt:jstop,            1,k,bid) ) + &
                    SUBCELLV(2:nx_block,jstrt:jstop,i_sub,j_sub,2,k,bid) * &
                  ( KXY_ISOP(2:nx_block,jstrt:jstop,            2,k,bid) + &
                    HXY_DIFF(2:nx_block,jstrt:jstop,            2,k,bid) ) ) / &
                            HUS(bstrt:bstop,jstrt:jstop,bid)

               do n=1,nt
                  FX(1:nx_block-1,2:ny_block,n) = FX(1:nx_block-1,2:ny_block,n) + &
                        WORK3(1:nx_block-1,2:ny_block) * &
                        TY(istrt:istop,astrt:astop,k,n,bid)
                  FY(2:nx_block,1:ny_block-1,n) = FY(2:nx_block,1:ny_block-1,n) + &
                        WORK4(2:nx_block,1:ny_block-1) * &
                        TX(bstrt:bstop,jstrt:jstop,k,n,bid)
               enddo
            enddo
         enddo
      endif
!-------------^--SJR ADD--^--------------------------------------------- 
      call timer_stop(timer_gm7, block_id=this_block%local_id)
      call timer_start(timer_gm7a, block_id=this_block%local_id)

      if ( .not. cancellation_occurs ) then

!-------------v--SJR MOD--v--------------------------------------------- 
         do n = 3,nt !TZ=0 for k==1... other TZ's calculated in init_meso_mixing 
             if (k < km) TZ(:,:,k+1,n,bid) = TMIX(:,:,k  ,n) - TMIX(:,:,k+1,n)
         enddo

         do k_sub=1,2
            k_shft = min(k-2+k_sub,km-1) !for dzw(k_shft)
            kk_shft = k_shft+1
            do i_sub=1,2  
               istrt = i_sub
               istop = nx_block-2+i_sub
               jstop = ny_block-2+i_sub

               if (isoonly) then
                  WORK3(1:nx_block-1,:) = & !FX \partial_z term
                           ( SUBCELLV(istrt:istop,:,i_sub,1,k_sub,k,bid) + &
                             SUBCELLV(istrt:istop,:,i_sub,2,k_sub,k,bid) ) * &
                         (   KXX_ISOP(istrt:istop,:,        k_sub,k,bid) * &
                               DD_SLX(istrt:istop,:,i_sub,  k_sub,k,bid) - & 
                              SF_SLXX(istrt:istop,:,i_sub,  k_sub,k,bid) ) / &
                                  dzw(k_shft)

                  WORK4(:,1:ny_block-1) = & !FY \partial_z term
                           ( SUBCELLV(:,istrt:jstop,1,i_sub,k_sub,k,bid) + &
                             SUBCELLV(:,istrt:jstop,2,i_sub,k_sub,k,bid) ) * &
                           ( KYY_ISOP(:,istrt:jstop,        k_sub,k,bid) * &
                               DD_SLY(:,istrt:jstop,  i_sub,k_sub,k,bid) - & 
                              SF_SLYY(:,istrt:jstop,  i_sub,k_sub,k,bid) ) / &
                                  dzw(k_shft)
               else
                  WORK3(1:nx_block-1,:) = & !FX \partial_z term
                     (       SUBCELLV(istrt:istop,:,i_sub,1,k_sub,k,bid) * &
                         ( ( KXX_ISOP(istrt:istop,:,        k_sub,k,bid) * &
                               DD_SLX(istrt:istop,:,i_sub,  k_sub,k,bid) - & 
                              SF_SLXX(istrt:istop,:,i_sub,  k_sub,k,bid) ) + &
                           ( KXY_ISOP(istrt:istop,:,        k_sub,k,bid) * &
                               DD_SLY(istrt:istop,:,      1,k_sub,k,bid) - & 
                              SF_SLYX(istrt:istop,:,      1,k_sub,k,bid) ) ) + &
                             SUBCELLV(istrt:istop,:,i_sub,2,k_sub,k,bid) * &
                         ( ( KXX_ISOP(istrt:istop,:,        k_sub,k,bid) * &
                               DD_SLX(istrt:istop,:,i_sub,  k_sub,k,bid) - & 
                              SF_SLXX(istrt:istop,:,i_sub,  k_sub,k,bid) ) + &
                           ( KXY_ISOP(istrt:istop,:,        k_sub,k,bid) * &
                               DD_SLY(istrt:istop,:,      2,k_sub,k,bid) - & 
                              SF_SLYX(istrt:istop,:,      2,k_sub,k,bid) ) )  ) / &
                                  dzw(k_shft)

                  WORK4(:,1:ny_block-1) = & !FY \partial_z term
                       (     SUBCELLV(:,istrt:jstop,1,i_sub,k_sub,k,bid) * &
                         ( ( KXY_ISOP(:,istrt:jstop,        k_sub,k,bid) * &
                               DD_SLX(:,istrt:jstop,1,      k_sub,k,bid) - & 
                              SF_SLXY(:,istrt:jstop,1,      k_sub,k,bid) ) + &
                           ( KYY_ISOP(:,istrt:jstop,        k_sub,k,bid) * &
                               DD_SLY(:,istrt:jstop,  i_sub,k_sub,k,bid) - & 
                              SF_SLYY(:,istrt:jstop,  i_sub,k_sub,k,bid) ) ) + &
                             SUBCELLV(:,istrt:jstop,2,i_sub,k_sub,k,bid) * &
                         ( ( KXY_ISOP(:,istrt:jstop,        k_sub,k,bid) * &
                               DD_SLX(:,istrt:jstop,2,      k_sub,k,bid) - & 
                              SF_SLXY(:,istrt:jstop,2,      k_sub,k,bid) ) + &
                           ( KYY_ISOP(:,istrt:jstop,        k_sub,k,bid) * &
                               DD_SLY(:,istrt:jstop,  i_sub,k_sub,k,bid) - & 
                              SF_SLYY(:,istrt:jstop,  i_sub,k_sub,k,bid) ) )  ) / &
                                  dzw(k_shft)

               endif

               if (partial_bottom_cells) then
                  WORK3 = WORK3 /  DZTW(:,:,k_shft,bid)*dzw(k_shft)
                  WORK4 = WORK4 /  DZTW(:,:,k_shft,bid)*dzw(k_shft)
               endif

               do n=1,nt
                  FX(1:nx_block-1,:,n) = FX(1:nx_block-1,:,n) + &
                           WORK3(1:nx_block-1,:) * TZ(istrt:istop,:,kk_shft,n,bid)
                  FY(:,1:ny_block-1,n) = FY(:,1:ny_block-1,n) + &
                           WORK4(:,1:ny_block-1) * TZ(:,istrt:jstop,kk_shft,n,bid)
               enddo
            enddo
         enddo

         if ( k < km ) then
            do k_sub=1,2
               kz_shft = k+2-k_sub
               do i_sub=1,2  
                  bstrt = 3-i_sub
                  bstop = nx_block+1-i_sub
                  astop = ny_block+1-i_sub

                  if (isoonly) then
                     WORK1(2:nx_block,:) = & !FZ \partial_x term
                              ( SUBCELLV(2:nx_block,:,i_sub,1,k_sub,kz_shft,bid) + &
                                SUBCELLV(2:nx_block,:,i_sub,2,k_sub,kz_shft,bid) ) * &
                              ( KXX_ISOP(2:nx_block,:,        k_sub,kz_shft,bid) * &
                                  DD_SLX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) + & 
                                 SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) ) / &
                                     HUS(bstrt:bstop,:,bid)

                     WORK2(:,2:ny_block) = & !FZ \partial_y term
                              ( SUBCELLV(:,2:ny_block,1,i_sub,k_sub,kz_shft,bid) + &
                                SUBCELLV(:,2:ny_block,2,i_sub,k_sub,kz_shft,bid) ) * &
                              ( KYY_ISOP(:,2:ny_block,        k_sub,kz_shft,bid) * &
                                  DD_SLY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) + & 
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) ) / &
                                     HUW(:,bstrt:astop,bid)
                  else
                     WORK1(2:nx_block,:) = & !FZ \partial_x term
                          (     SUBCELLV(2:nx_block,:,i_sub,1,k_sub,kz_shft,bid) * &
                            ( ( KXX_ISOP(2:nx_block,:,        k_sub,kz_shft,bid) * &
                                  DD_SLX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) + & 
                                 SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) ) + &
                              ( KXY_ISOP(2:nx_block,:,        k_sub,kz_shft,bid) * &
                                  DD_SLY(2:nx_block,:,      1,k_sub,kz_shft,bid) + & 
                                 SF_SLYX(2:nx_block,:,      1,k_sub,kz_shft,bid) ) ) + &
                                SUBCELLV(2:nx_block,:,i_sub,2,k_sub,kz_shft,bid) * &
                            ( ( KXX_ISOP(2:nx_block,:,        k_sub,kz_shft,bid) * &
                                  DD_SLX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) + & 
                                 SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) ) + &
                              ( KXY_ISOP(2:nx_block,:,        k_sub,kz_shft,bid) * &
                                  DD_SLY(2:nx_block,:,      2,k_sub,kz_shft,bid) + & 
                                 SF_SLYX(2:nx_block,:,      2,k_sub,kz_shft,bid) ) ) ) / &
                                  HUS(bstrt:bstop,:,bid)

                     WORK2(:,2:ny_block) = & !FZ \partial_y term
                          (     SUBCELLV(:,2:ny_block,1,i_sub,k_sub,kz_shft,bid) * &
                            ( ( KXY_ISOP(:,2:ny_block,        k_sub,kz_shft,bid) * &
                                  DD_SLX(:,2:ny_block,1,      k_sub,kz_shft,bid) + & 
                                 SF_SLXY(:,2:ny_block,1,      k_sub,kz_shft,bid) ) + &
                              ( KYY_ISOP(:,2:ny_block,        k_sub,kz_shft,bid) * &
                                  DD_SLY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) + & 
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) ) ) + &
                                SUBCELLV(:,2:ny_block,2,i_sub,k_sub,kz_shft,bid) * &
                            ( ( KXY_ISOP(:,2:ny_block,        k_sub,kz_shft,bid) * &
                                  DD_SLX(:,2:ny_block,2,      k_sub,kz_shft,bid) + & 
                                 SF_SLXY(:,2:ny_block,2,      k_sub,kz_shft,bid) ) + &
                              ( KYY_ISOP(:,2:ny_block,        k_sub,kz_shft,bid) * &
                                  DD_SLY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) + & 
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) ) ) ) / &
                                     HUW(:,bstrt:astop,bid)
                  endif

                  do n=1,nt
                     FZBOT(2:nx_block,2:ny_block,n) = FZBOT(2:nx_block,2:ny_block,n) + &
                           WORK1(2:nx_block,2:ny_block) * TX(bstrt:bstop,2:ny_block,kz_shft,n,bid) + &
                           WORK2(2:nx_block,2:ny_block) * TY(2:nx_block,bstrt:astop,kz_shft,n,bid)
                  enddo
               enddo
            enddo
         end if
      elseif ( k < km ) then
         do k_sub=1,2
            kz_shft = k+2-k_sub
            do i_sub=1,2  
               bstrt = 3-i_sub
               bstop = nx_block+1-i_sub
               astop = ny_block+1-i_sub

               if (isoonly) then
                  WORK1(2:nx_block,:) = c2 * & !FZ \partial_x term
                              ( SUBCELLV(2:nx_block,:,i_sub,1,k_sub,kz_shft,bid) + &
                                SUBCELLV(2:nx_block,:,i_sub,2,k_sub,kz_shft,bid) ) * &
                                 SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) / &
                                  HUS(bstrt:bstop,:,bid)

                  WORK2(:,2:ny_block) = c2 * & !FZ \partial_y term
                              ( SUBCELLV(:,2:ny_block,1,i_sub,k_sub,kz_shft,bid) + &
                                SUBCELLV(:,2:ny_block,2,i_sub,k_sub,kz_shft,bid) ) * &
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid) / &
                                     HUW(:,bstrt:astop,bid)
               else
                  WORK1(2:nx_block,:) = c2 * & !FZ \partial_x term
                          (     SUBCELLV(2:nx_block,:,i_sub,1,k_sub,kz_shft,bid) * &
                                (SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) + &
                                 SF_SLYX(2:nx_block,:,      1,k_sub,kz_shft,bid)) + &
                                SUBCELLV(2:nx_block,:,i_sub,2,k_sub,kz_shft,bid) * &
                                (SF_SLXX(2:nx_block,:,i_sub,  k_sub,kz_shft,bid) + &
                                 SF_SLYX(2:nx_block,:,      2,k_sub,kz_shft,bid)) ) / &
                                  HUS(bstrt:bstop,:,bid)

                  WORK2(:,2:ny_block) = c2 * & !FZ \partial_y term
                          (     SUBCELLV(:,2:ny_block,1,i_sub,k_sub,kz_shft,bid) * &
                                (SF_SLXY(:,2:ny_block,1,      k_sub,kz_shft,bid) + &
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid)) + &
                                SUBCELLV(:,2:ny_block,2,i_sub,k_sub,kz_shft,bid) * &
                                (SF_SLXY(:,2:ny_block,2,      k_sub,kz_shft,bid) + &
                                 SF_SLYY(:,2:ny_block,  i_sub,k_sub,kz_shft,bid)) ) / &
                                     HUW(:,bstrt:astop,bid)
               endif

               do n=1,nt
                  FZBOT(2:nx_block,2:ny_block,n) = FZBOT(2:nx_block,2:ny_block,n) + &
                        WORK1(2:nx_block,2:ny_block) * TX(bstrt:bstop,2:ny_block,kz_shft,n,bid) + &
                        WORK2(2:nx_block,2:ny_block) * TY(2:nx_block,bstrt:astop,kz_shft,n,bid)
               enddo
            enddo
         enddo
      end if
!-------------^--SJR MOD--^--------------------------------------------- 
      call timer_stop(timer_gm7a, block_id=this_block%local_id)
      call timer_start(timer_gm7b, block_id=this_block%local_id)

!-------------v--SJR ADD--v--------------------------------------------- 
      if ( vertdiffhere ) then
         do n=1,nt
               FZBOT(:,:,n) = FZBOT(:,:,n) + VDC_GM(:,:,k,bid) * TZ(:,:,k+1,n,bid) * TAREA(:,:,bid) !Remove volume normalization... dz cancels for PBCs or not 
         enddo
      end if
!-------------^--SJR ADD--^--------------------------------------------- 

!-------------v--SJR ADD--v--------------------------------------------- 
      KMASKE = merge(c1, c0, k <= KMT(:,:,bid) .and. k <= KMTE(:,:,bid))
      KMASKN = merge(c1, c0, k <= KMT(:,:,bid) .and. k <= KMTN(:,:,bid))
      do n=1,nt
         FX(:,:,n)     =    FX(:,:,n)     / HUS(:,:,bid) * KMASKE
         FY(:,:,n)     =    FY(:,:,n)     / HUW(:,:,bid) * KMASKN
         if (partial_bottom_cells) then
            FZBOT(:,:,n) = FZBOT(:,:,n) * KMASK / DZTW(:,:,k,bid) 
         else
            FZBOT(:,:,n) = FZBOT(:,:,n) * KMASK / dzw(k) 
         endif
      enddo
!-------------^--SJR ADD--^--------------------------------------------- 

!-------------v--SJR MOD--v--------------------------------------------- 
      GTK = c0

      if (partial_bottom_cells) then
         do n=1,nt
            GTK(2:nx_block-1,2:ny_block-1,n) =                                &
                ( ( FX(2:nx_block-1,2:ny_block-1,n) *                         &
                    MIN(DZT(3:nx_block  ,2:ny_block-1,k,bid),                 &
                        DZT(2:nx_block-1,2:ny_block-1,k,bid)) -               &
                    FX(1:nx_block-2,2:ny_block-1,n) *                         &
                    MIN(DZT(1:nx_block-2,2:ny_block-1,k,bid),                 &
                        DZT(2:nx_block-1,2:ny_block-1,k,bid)) +               &
                    FY(2:nx_block-1,2:ny_block-1,n) *                         &
                    MIN(DZT(2:nx_block-1,3:ny_block  ,k,bid),                 &
                        DZT(2:nx_block-1,2:ny_block-1,k,bid)) -               &
                    FY(2:nx_block-1,1:ny_block-2,n) *                         &
                    MIN(DZT(2:nx_block-1,1:ny_block-2,k,bid),                 &
                        DZT(2:nx_block-1,2:ny_block-1,k,bid)) ) /             &
                  DZT(2:nx_block-1,2:ny_block-1,k,bid) +                      &
                  FZTOP(2:nx_block-1,2:ny_block-1,n,bid) -                    &
                  FZBOT(2:nx_block-1,2:ny_block-1,n) ) /                      &
                DZT(2:nx_block-1,2:ny_block-1,k,bid)
         enddo
      else
         GTK(2:nx_block-1,2:ny_block-1,:) = dzr(k) * &
             (       FX(2:nx_block-1,2:ny_block-1,:) - &
                     FX(1:nx_block-2,2:ny_block-1,:) + &
                     FY(2:nx_block-1,2:ny_block-1,:) - &
                     FY(2:nx_block-1,1:ny_block-2,:) + &
                  FZTOP(2:nx_block-1,2:ny_block-1,:,bid) - &
                  FZBOT(2:nx_block-1,2:ny_block-1,:) )
      endif
      do n=1,nt
         GTK(:,:,n) = GTK(:,:,n) * TAREA_R(:,:,bid)
      enddo

      FZTOP(:,:,:,bid) = FZBOT(:,:,:)
!-------------^--SJR MOD--^--------------------------------------------- 

      call timer_stop(timer_gm7b, block_id=this_block%local_id)
      call timer_stop(timer_nloop, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!
!     diagnostic computation of the bolus velocities 
!
!-----------------------------------------------------------------------
      call timer_start(timer_gm8, block_id=this_block%local_id)
      
      if ( diag_gm_bolus ) then

!-------------v--SJR MOD--v--------------------------------------------- 
        WORK1 = c0 ! \psi_y at i+1/2,j,k+1/2 (bottom of east  face, center in y)
        WORK2 = c0 ! \psi_x at i,j+1/2,k+1/2 (bottom of north face, center in x)
        do k_sub=1,2  !ktp(k+1),kbt(k)
           kk_shft = min(k+2-k_sub,km)
           do j_sub=1,2  !jnorth(j),jsouth(j) for WORK1 and jnorth(j),jsouth(j+1) for WORK2 
                         !ieast(i),iwest(i+1) for WORK1 and ieast(i),iwest(i) for WORK2
              jstrt = j_sub
              jstop = ny_block-2+j_sub
              istop = nx_block-2+j_sub

              if (isoonly) then
                 WORK1(1:nx_block-1,:) = WORK1(1:nx_block-1,:) + &
                                    ( SUBCELLV(jstrt:istop ,:,j_sub,1,k_sub,kk_shft,bid) + &
                                      SUBCELLV(jstrt:istop ,:,j_sub,2,k_sub,kk_shft,bid) ) * &
                                       SF_SLXX(jstrt:istop ,:,j_sub,  k_sub,kk_shft,bid) 
                 WORK2(:,1:ny_block-1) = WORK2(:,1:ny_block-1) + &
                                    ( SUBCELLV(:, jstrt:jstop,1,j_sub,k_sub,kk_shft,bid) + &
                                      SUBCELLV(:, jstrt:jstop,2,j_sub,k_sub,kk_shft,bid) ) * &
                                       SF_SLYY(:, jstrt:jstop,  j_sub,k_sub,kk_shft,bid) 
              else
                 WORK1(1:nx_block-1,:) = WORK1(1:nx_block-1,:) + &
                                      SUBCELLV(jstrt:istop ,:,j_sub,1,k_sub,kk_shft,bid) * &
                                     ( SF_SLXX(jstrt:istop ,:,j_sub,  k_sub,kk_shft,bid) + &
                                       SF_SLYX(jstrt:istop ,:,      1,k_sub,kk_shft,bid) ) + & 
                                      SUBCELLV(jstrt:istop ,:,j_sub,2,k_sub,kk_shft,bid) * &
                                     ( SF_SLXX(jstrt:istop ,:,j_sub,  k_sub,kk_shft,bid) + &
                                       SF_SLYX(jstrt:istop ,:,      2,k_sub,kk_shft,bid) ) 
                 WORK2(:,1:ny_block-1) = WORK2(:,1:ny_block-1) + &
                                      SUBCELLV(:, jstrt:jstop,1,j_sub,k_sub,kk_shft,bid) * &
                                     ( SF_SLXY(:, jstrt:jstop,1,      k_sub,kk_shft,bid) + &
                                       SF_SLYY(:, jstrt:jstop,  j_sub,k_sub,kk_shft,bid) ) + & 
                                      SUBCELLV(:, jstrt:jstop,2,j_sub,k_sub,kk_shft,bid) * &
                                     ( SF_SLXY(:, jstrt:jstop,2,      k_sub,kk_shft,bid) + &
                                       SF_SLYY(:, jstrt:jstop,  j_sub,k_sub,kk_shft,bid) ) 
              endif
           enddo
        enddo
        if (partial_bottom_cells) then
           WORK1(:,:) =  WORK1(:,:) / HUS(:,:,bid) / HTE(:,:,bid) / DZTW(:,:,k,bid)  
           WORK2(:,:) = -WORK2(:,:) / HUW(:,:,bid) / HTN(:,:,bid) / DZTW(:,:,k,bid)  
        else
           WORK1(:,:) =  WORK1(:,:) / HUS(:,:,bid) / HTE(:,:,bid) / dzw(k)  
           WORK2(:,:) = -WORK2(:,:) / HUW(:,:,bid) / HTN(:,:,bid) / dzw(k)  
        endif
        if (k == km) then
           WORK1 = p5 * WORK1 !instead of "factor"
           WORK2 = p5 * WORK2
        end if
!-------------^--SJR MOD--^--------------------------------------------- 

        UIB = merge( WORK1, c0, k < KMT(:,:,bid)  &
                          .and. k < KMTE(:,:,bid) )
        VIB = merge( WORK2, c0, k < KMT(:,:,bid)  &
                          .and. k < KMTN(:,:,bid) )

!-------------v--SJR MOD--v--------------------------------------------- 
        if (partial_bottom_cells) then
           U_ISOP = merge( -( UIT(:,:,bid) - UIB ) / DZT(:,:,k,bid) , c0, &       !u_b = -\partial_z(\psi_y)
                            k <= KMT(:,:,bid) .and. k <= KMTE(:,:,bid) )
           V_ISOP = merge(  ( VIT(:,:,bid) - VIB ) / DZT(:,:,k,bid) , c0, &       !v_b =  \partial_z(\psi_x)
                            k <= KMT(:,:,bid) .and. k <= KMTN(:,:,bid) )
           !WBOT_ISOP(2:nx_block,2:ny_block,bid) =  &   !w_b = \partial_x(\psi_y) - \partial_y(\psi_x)
           !         merge(  ( WORK1(2:nx_block,2:ny_block) - WORK1(1:nx_block-1,2:ny_block) ) / &
           !                     DXT(2:nx_block,2:ny_block,bid) - &
           !                 ( WORK2(2:nx_block,2:ny_block) - WORK2(2:nx_block,1:ny_block-1) ) / &
           !                     DYT(2:nx_block,2:ny_block,bid) , c0, k < KMT(2:nx_block,2:ny_block,bid) )  
              !w_b = \partial_x(\psi_y) - \partial_y(\psi_x) = w_b^top + dz*(dxu+dyv)
           WBOT_ISOP(2:nx_block,2:ny_block,bid) = &
                  merge( WTOP_ISOP(2:nx_block,2:ny_block,  bid) + &
                               DZT(2:nx_block,2:ny_block,k,bid) * &   
                          ( U_ISOP(2:nx_block,2:ny_block) - U_ISOP(1:nx_block-1,2:ny_block) ) / &
                               DXT(2:nx_block,2:ny_block,  bid) + &
                          ( V_ISOP(2:nx_block,2:ny_block) - V_ISOP(2:nx_block,1:ny_block-1) ) / &
                               DYT(2:nx_block,2:ny_block,  bid) , c0, k < KMT(2:nx_block,2:ny_block,bid) )  
        else
           U_ISOP = merge( -( UIT(:,:,bid) - UIB ) * dzr(k) , c0, &       !u_b = -\partial_z(\psi_y)
                            k <= KMT(:,:,bid) .and. k <= KMTE(:,:,bid) )
           V_ISOP = merge(  ( VIT(:,:,bid) - VIB ) * dzr(k) , c0, &       !v_b =  \partial_z(\psi_x)
                            k <= KMT(:,:,bid) .and. k <= KMTN(:,:,bid) )
           !WBOT_ISOP(2:nx_block,2:ny_block,bid) =  &   !w_b = \partial_x(\psi_y) - \partial_y(\psi_x)
           !         merge(  ( WORK1(2:nx_block,2:ny_block) - WORK1(1:nx_block-1,2:ny_block) ) / &
           !                     DXT(2:nx_block,2:ny_block,bid) - &
           !                 ( WORK2(2:nx_block,2:ny_block) - WORK2(2:nx_block,1:ny_block-1) ) / &
           !                     DYT(2:nx_block,2:ny_block,bid) , c0, k < KMT(2:nx_block,2:ny_block,bid) )  
              !w_b = \partial_x(\psi_y) - \partial_y(\psi_x) = w_b^top + dz*(dxu+dyv)
           WBOT_ISOP(2:nx_block,2:ny_block,bid) = &
                  merge( WTOP_ISOP(2:nx_block,2:ny_block,bid) + dz(k) * &   
                          ( U_ISOP(2:nx_block,2:ny_block) - U_ISOP(1:nx_block-1,2:ny_block) ) / &
                               DXT(2:nx_block,2:ny_block,bid) + &
                          ( V_ISOP(2:nx_block,2:ny_block) - V_ISOP(2:nx_block,1:ny_block-1) ) / &
                               DYT(2:nx_block,2:ny_block,bid) , c0, k < KMT(2:nx_block,2:ny_block,bid) )  
        endif
!-------------^--SJR MOD--^--------------------------------------------- 

        if ( linertial .and. k == 2 ) then
          BOLUS_SP(:,:,bid) = 50.0_r8 * sqrt(U_ISOP**2 + V_ISOP**2) !SJR MOD - integer exponent for efficiency  
        endif

      endif

!-----------------------------------------------------------------------
!
!     update remaining bottom-face fields to top-face fields for next
!     pass
!
!-----------------------------------------------------------------------

      
      if ( diag_gm_bolus ) then
        UIT(:,:,bid) = UIB
        VIT(:,:,bid) = VIB
      endif

!-----------------------------------------------------------------------
!
!     compute isopycnal diffusion cfl diagnostics if required
!
!-----------------------------------------------------------------------

      if (ldiag_cfl) then

        WORK2 = p5 * (KAPPA_ISOP(:,:,ktp,k,bid)  &
                    + KAPPA_ISOP(:,:,kbt,k,bid)) 

        WORK1 = merge(c4*WORK2*(DXTR(:,:,bid)**2 + DYTR(:,:,bid)**2),  &
                      c0, KMT(:,:,bid) > k)
        WORK2 = abs(WORK1)
        call cfl_hdiff (k,bid,WORK2,1,this_block)

      endif
      call timer_stop(timer_gm8, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!
!     accumulate time average if necessary; testing is internal to
!       accumulate_tavg_field
!
!-----------------------------------------------------------------------

      if ( mix_pass /= 1 ) then
      call timer_start(timer_gm9, block_id=this_block%local_id)

          call accumulate_tavg_field                      &
                   (p5 * (KAPPA_ISOP(:,:,ktp,k,bid)    &
                       +  KAPPA_ISOP(:,:,kbt,k,bid)),  &
                          tavg_KAPPA_ISOP, bid, k)

         if (savenewtavgs) then
          WORK1 = p5 * ( KXX_ISOP(:,:,ktp,k,bid) + KXX_ISOP(:,:,kbt,k,bid) )
          WORK2 = p5 * ( KXY_ISOP(:,:,ktp,k,bid) + KXY_ISOP(:,:,kbt,k,bid) )
          WORK3 = p5 * ( KYY_ISOP(:,:,ktp,k,bid) + KYY_ISOP(:,:,kbt,k,bid) )

          call accumulate_tavg_field(WORK1,tavg_KXX_ISOP, bid, k)
          call accumulate_tavg_field(WORK2,tavg_KXY_ISOP, bid, k)
          call accumulate_tavg_field(WORK3,tavg_KYY_ISOP, bid, k)

          WORK4 = sqrt( (WORK1+WORK3)**2 - c4*(WORK1*WORK3-WORK2**2) )/c2 
          WORK5 = (WORK1+WORK3)/c2+WORK4
          WORK6 = (WORK1+WORK3)/c2-WORK4
          WORK4 = WORK2**2+(WORK5-WORK1)**2
          where (WORK4 < eps2)
             WORK3 = c0
             WORK2 = c1
          elsewhere
             WORK3 = c2*WORK2*(WORK5-WORK1) / WORK4 !sin(2*theta)
             WORK2 = ( WORK2**2-(WORK5-WORK1)**2 ) / WORK4 !cos(2*theta)
          endwhere

          call accumulate_tavg_field(WORK5,tavg_MAJ_ISOP, bid, k)
          call accumulate_tavg_field(WORK6,tavg_MIN_ISOP, bid, k)
          call accumulate_tavg_field(WORK2,tavg_C2T_ISOP, bid, k)
          call accumulate_tavg_field(WORK3,tavg_S2T_ISOP, bid, k)
         endif

          call accumulate_tavg_field                      &
                   (p5 * (KAPPA_THIC(:,:,ktp,k,bid)    &
                       +  KAPPA_THIC(:,:,kbt,k,bid)),  &
                          tavg_KAPPA_THIC, bid, k)

         if (savenewtavgs) then
          WORK1 = p5 * ( KXX_THIC(:,:,ktp,k,bid) + KXX_THIC(:,:,kbt,k,bid) )
          WORK2 = p5 * ( KXY_THIC(:,:,ktp,k,bid) + KXY_THIC(:,:,kbt,k,bid) )
          WORK3 = p5 * ( KYY_THIC(:,:,ktp,k,bid) + KYY_THIC(:,:,kbt,k,bid) )

          call accumulate_tavg_field(WORK1,tavg_KXX_THIC, bid, k)
          call accumulate_tavg_field(WORK2,tavg_KXY_THIC, bid, k)
          call accumulate_tavg_field(WORK3,tavg_KYY_THIC, bid, k)

          WORK4 = sqrt( (WORK1+WORK3)**2 - c4*(WORK1*WORK3-WORK2**2) )/c2 
          WORK5 = (WORK1+WORK3)/c2+WORK4
          WORK6 = (WORK1+WORK3)/c2-WORK4
          WORK4 = WORK2**2+(WORK5-WORK1)**2
          where (WORK4 < eps2)
             WORK3 = c0
             WORK2 = c1
          elsewhere
             WORK3 = c2*WORK2*(WORK5-WORK1) / WORK4 !sin(2*theta)
             WORK2 = ( WORK2**2-(WORK5-WORK1)**2 ) / WORK4 !cos(2*theta)
          endwhere

          call accumulate_tavg_field(WORK5,tavg_MAJ_THIC, bid, k)
          call accumulate_tavg_field(WORK6,tavg_MIN_THIC, bid, k)
          call accumulate_tavg_field(WORK2,tavg_C2T_THIC, bid, k)
          call accumulate_tavg_field(WORK3,tavg_S2T_THIC, bid, k)
         endif

          call accumulate_tavg_field                      &
                   (p5 * (HOR_DIFF(:,:,ktp,k,bid)      &
                       +  HOR_DIFF(:,:,kbt,k,bid)),    &
                          tavg_HOR_DIFF, bid, k)

         if (savenewtavgs) then
          WORK1 = p5 * ( HXX_DIFF(:,:,ktp,k,bid) + HXX_DIFF(:,:,kbt,k,bid) )
          WORK2 = p5 * ( HXY_DIFF(:,:,ktp,k,bid) + HXY_DIFF(:,:,kbt,k,bid) )
          WORK3 = p5 * ( HYY_DIFF(:,:,ktp,k,bid) + HYY_DIFF(:,:,kbt,k,bid) )

          call accumulate_tavg_field(WORK1,tavg_HXX_DIFF, bid, k)
          call accumulate_tavg_field(WORK2,tavg_HXY_DIFF, bid, k)
          call accumulate_tavg_field(WORK3,tavg_HYY_DIFF, bid, k)

          WORK4 = sqrt( (WORK1+WORK3)**2 - c4*(WORK1*WORK3-WORK2**2) )/c2 
          WORK5 = (WORK1+WORK3)/c2+WORK4
          WORK6 = (WORK1+WORK3)/c2-WORK4
          WORK4 = WORK2**2+(WORK5-WORK1)**2
          where (WORK4 < eps2)
             WORK3 = c0
             WORK2 = c1
          elsewhere
             WORK3 = c2*WORK2*(WORK5-WORK1) / WORK4 !sin(2*theta)
             WORK2 = ( WORK2**2-(WORK5-WORK1)**2 ) / WORK4 !cos(2*theta)
          endwhere

          call accumulate_tavg_field(WORK5,tavg_MAJ_DIFF, bid, k)
          call accumulate_tavg_field(WORK6,tavg_MIN_DIFF, bid, k)
          call accumulate_tavg_field(WORK2,tavg_C2T_DIFF, bid, k)
          call accumulate_tavg_field(WORK3,tavg_S2T_DIFF, bid, k)
         endif

        if ( transition_layer_on  .and.  k == 1 ) then

            call accumulate_tavg_field (TLT%DIABATIC_DEPTH(:,:,bid),  &
                                        tavg_DIA_DEPTH, bid, 1)  

            call accumulate_tavg_field (TLT%THICKNESS(:,:,bid),       &
                                        tavg_TLT, bid, 1)  

            call accumulate_tavg_field (TLT%INTERIOR_DEPTH(:,:,bid),  &
                                        tavg_INT_DEPTH, bid, 1)  

        endif

        if ( diag_gm_bolus ) then

            call accumulate_tavg_field (U_ISOP, tavg_UISOP, bid, k) 
            call accumulate_tavg_field (V_ISOP, tavg_VISOP, bid, k) 
            call accumulate_tavg_field (WTOP_ISOP(:,:,bid), tavg_WISOP,bid, k)

          if (accumulate_tavg_now(tavg_ADVT_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )  
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                 if (partial_bottom_cells) then
                  WORK1(i,j) = - DZT(i,j,k,bid) * TAREA_R(i,j,bid) * WORK3(i,j)
                 else
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                 endif
                endif
              enddo
            enddo

            call accumulate_tavg_field (WORK1, tavg_ADVT_ISOP, bid, k)

          endif

           if (accumulate_tavg_now(tavg_ADVS_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            if (partial_bottom_cells) then
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - DZT(i,j,k,bid) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            else
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            endif

            call accumulate_tavg_field (WORK1, tavg_ADVS_ISOP, bid, k)

          endif

          if ( accumulate_tavg_now(tavg_VNT_ISOP)  .or.  &
               accumulate_tavg_now(tavg_VNS_ISOP) ) then

            WORK1 = p5 * V_ISOP * HTN(:,:,bid) * TAREA_R(:,:,bid) 

            if (accumulate_tavg_now(tavg_VNT_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,1)  &
                         + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNT_ISOP, bid, k) 
            endif

            if (accumulate_tavg_now(tavg_VNS_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,2)  &
                         + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNS_ISOP, bid, k)
            endif

          endif

        endif ! bolus velocity option on

        do n = 1,nt
          if (accumulate_tavg_now(tavg_HDIFE_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FX(i,j,n)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FX(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FY(i,j,n)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FY(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FZTOP(i,j,n,bid)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FZTOP(i,j,n,bid)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFB_TRACER(n),bid,k)
          endif
        enddo

      call timer_stop(timer_gm9, block_id=this_block%local_id)
      endif   ! mix_pass ne 1

      call timer_stop(timer_gm0, block_id=this_block%local_id)

!-----------------------------------------------------------------------
!EOC

      end subroutine hdifft_gm_aniso

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_vmhs 
! !INTERFACE:

      subroutine kappa_lon_lat_vmhs (TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization by Visbeck et al (1997):
!  \begin{equation}
!     KAPPA_LATERAL = C {{l^2}\over{T}},
!  \end{equation}
!  where $C$ is a constant, $T$ is the (Eady) time scale
!  $f/\surd\overline{Ri}$ and $l$ is the Rossby radius.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX               ! tracers at all vertical levels
                            !  at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX         ! U,V  at all vertical levels
                            !  at mixing time level

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         kp1,i,j,kk,k,      &
         k1,k2,             &
         bid                ! local block address for this sub block

      real (r8) :: &
         const, zmin1, zmin2

      real (r8), dimension(nx_block,ny_block) :: &
         UKT, VKT,                               &
         UKPT, VKPT, RNUM,                       &
         WORK1, WORK2, WORK3, WORK4,             &
         TK, TKP,                                & 
         RHOT, RHOS,                             &
         GRATE, LSC                              

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

      bid = this_block%local_id

      GRATE = c0
      LSC = c0
      k1 = 0
      k2 = 0

      vert_loop: do k=1,km

!     -2000m < z < -100m

        if ( zt(k) >= 1.0e4_r8 .and. zt(k) <= 2.0e5_r8 ) then

!-----------------------------------------------------------------------
!
!     start computing the Richardson number and the baroclinic wave speed 
!     average velocities at T points of level k
!
!-----------------------------------------------------------------------

          if ( k1 == 0 ) then

            k1 = k              ! k1 is the upper limit for integration

            call ugrid_to_tgrid(UKT,UMIX(:,:,k),bid)
            call ugrid_to_tgrid(VKT,VMIX(:,:,k),bid)

            TK = max(-c2, TMIX(:,:,k,1))

          endif

!-----------------------------------------------------------------------
!
!     compute RZ=Dz(rho) at k=k
!
!-----------------------------------------------------------------------

          if ( k < km ) then
            kp1 = k + 1
          else
            kp1 = km
          endif

          call state (k,kp1,TMIX(:,:,k,1),TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS)

          TKP = max(-c2, TMIX(:,:,kp1,1))

          WORK1 = TK - TKP
          WORK2 = TMIX(:,:,k,2) - TMIX(:,:,kp1,2)

          WORK3 = RHOT * WORK1 + RHOS * WORK2
          WORK3 = min(WORK3,-eps2)

!-----------------------------------------------------------------------
!
!     average velocities at T points of level k+1
!
!-----------------------------------------------------------------------

          call ugrid_to_tgrid(UKPT,UMIX(:,:,kp1),bid)
          call ugrid_to_tgrid(VKPT,VMIX(:,:,kp1),bid)

!-----------------------------------------------------------------------
!
!     local Richardson number at the bottom of T box, level k
!     = -g*Dz(RHO)/RHO_0/(Dz(u)**2+Dz(v)**2)
!
!-----------------------------------------------------------------------

          if (partial_bottom_cells) then

            where (k < KMT(:,:,bid))

!     RHO_0 = 1 in denominator (approximately) in cgs unit.

              WORK4 = p5*(dz(k) + DZT(:,:,k+1,bid))  ! dzw equivalent

              RNUM = -WORK4/((UKT - UKPT)**2 + (VKT - VKPT)**2 + eps)
              GRATE = GRATE + grav*RNUM*WORK4*WORK3
              LSC = LSC - grav*WORK3

            end where

          else ! no partial bottom cells

            where (k < KMT(:,:,bid))

!     RHO_0 = 1 in denominator (approximately) in cgs unit.

              RNUM = -dzw(k)/((UKT - UKPT)**2 + (VKT - VKPT)**2 + eps)
              GRATE = GRATE + grav*RNUM*dzw(k)*WORK3
              LSC = LSC - grav*WORK3

            end where

          endif ! partial bottom cells

!-----------------------------------------------------------------------
!
!     bottom values become top values for next pass
!
!-----------------------------------------------------------------------

          UKT = UKPT
          VKT = VKPT
          TK  = TKP

        else

          if ( k1 /= 0 .and. k2 == 0 ) then
            k2 = k                ! k2 is the lower limit for integration
            exit vert_loop
          endif

        endif                     ! if(zt(k).ge.1.0e4.and.zt(k).lt.2.0e5)

      end do vert_loop

      do j=1,ny_block
        do i=1,nx_block
          kk = KMT(i,j,bid)
          if ( kk > 0 ) then
            zmin1 = min(zt(k1),zt(kk))
            zmin2 = min(zt(k2),zt(kk))
            GRATE(i,j) = GRATE(i,j)/(zmin2-zmin1+eps) ! Ri
            LSC(i,j) = LSC(i,j)*(zmin2-zmin1)         ! c_g^2=N^2*H^2
          else
            GRATE(i,j) = c0
            LSC(i,j)   = c0
          endif
        end do
      end do

!     equatorial inverse of time scale and square of length scale

      WORK1 = sqrt(c2*sqrt(LSC)*BTP(:,:,bid)) ! sqrt(c_g*2*beta)
      WORK2 = sqrt(LSC)/(c2*BTP(:,:,bid))     ! c_g/(2 beta)

!     inverse of time scale

      WORK1 = max(abs(FCORT(:,:,bid)),WORK1)
      GRATE = WORK1/sqrt(GRATE+eps)           ! 1/T = f/sqrt(Ri)

!     square of length scale

      LSC = LSC/(FCORT(:,:,bid)+eps)**2       ! L^2 = c_g^2/f^2
      LSC = min(LSC,WORK2)

!     grid size = lower limit of length scale

      WORK1 = min(DXT(:,:,bid)**2,DYT(:,:,bid)**2)
      LSC = max(LSC,WORK1)

!-----------------------------------------------------------------------
!
!     compute KAPPA_LATERAL
!
!-----------------------------------------------------------------------

!     const = 0.015_r8  ! constant taken from Visbeck et al (1997)
      const = 0.13_r8   
      WORK1 = const*GRATE*LSC                 ! const*L**2/T

!     KAPPA_LATERAL is bounded by 3.0e6 < KAPPA_LATERAL < 4.0e7

      KAPPA_LATERAL(:,:,bid) = min(4.0e7_r8,WORK1)
      KAPPA_LATERAL(:,:,bid) = max(3.0e6_r8,KAPPA_LATERAL(:,:,bid))

      where (KMT(:,:,bid) <= k1)     ! KAPPA_LATERAL=3.0e6 when depth < 100m
        KAPPA_LATERAL(:,:,bid) = 3.0e6_r8
      end where

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_vmhs 

!***********************************************************************
!BOP
! !IROUTINE: kappa_eg
! !INTERFACE:

      subroutine kappa_eg (TMIX, UMIX, VMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization by Eden and Greatbatch 
!  (2008, Ocean Modelling, v. 20, 223-239):
!  \begin{equation}
!     KAPPA_ISOP = KAPPA_THIC = c {{L^2}*{\sigma}},
!  \end{equation}
!  where $c$ is a tuning parameter (=const_eg), $\sigma$ is the inverse eddy 
!  time scale, and $L$ is the minimum of the Rossby deformation radius and
!  Rhines scale.
!
!  This subroutine returns KAPPA_ISOP in KAPPA_VERTICAL. Here, KAPPA_VERTICAL
!  serves as a temporary array because this subroutine may be called less
!  frequently than every time step and KAPPA_ISOP is modified in each time step.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX               ! tracers at all vertical levels
                            !  at mixing time level

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX         ! U,V at all vertical levels
                            !  at mixing time level

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kp1,            &  ! vertical level index
         bid                   ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         TK, TKP,                                &  ! level k and k+1 TMIX
         WORK1, WORK2, WORK3,                    &  ! work arrays
         C_ROSSBY,                               &  ! first baroclinic Rossby wave speed
         L_ROSSBY,                               &  ! Rossby deformation radius
         RHOT, RHOS                                 ! dRHOdT and dRHOdS


      real (r8), dimension(nx_block,ny_block,km) :: &
         SIGMA,         &  ! inverse eddy time scale
         RI_NO             ! local Richardson number 

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

      bid = this_block%local_id

      SIGMA = c0
      RI_NO = c0

      BUOY_FREQ_SQ(:,:,:,bid) = c0

!-----------------------------------------------------------------------
!
!     compute buoyancy frequency and Richardson number at the bottom of 
!     T box:
!             Ri = -g*Dz(RHO)/RHO_0/(Dz(u)**2+Dz(v)**2)
!
!     RHO_0 ~ 1 in cgs units.
!
!-----------------------------------------------------------------------

      do k=1,km-1

        if ( k == 1 ) then
          TK = max(-c2, TMIX(:,:,k,1))
        endif

        kp1 = k+1

        call state (k, kp1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                    this_block, DRHODT=RHOT, DRHODS=RHOS)

        TKP = max(-c2, TMIX(:,:,kp1,1))

        WORK1 = TK - TKP
        WORK2 = TMIX(:,:,k,2) - TMIX(:,:,kp1,2)
        WORK3 = RHOT * WORK1 + RHOS * WORK2
        WORK3 = min(WORK3,-eps2)

        WORK1 =  (UMIX(:,:,k) - UMIX(:,:,kp1))**2  &
               + (VMIX(:,:,k) - VMIX(:,:,kp1))**2

        call ugrid_to_tgrid (WORK2, WORK1, bid)

        if (partial_bottom_cells) then
         where (k < KMT(:,:,bid))
          BUOY_FREQ_SQ(:,:,k,bid) = - grav * WORK3 / DZTW(:,:,k,bid)
          WORK1 = DZTW(:,:,k,bid)**2 / ( WORK2 + eps2 )
          RI_NO(:,:,k) = WORK1 * BUOY_FREQ_SQ(:,:,k,bid) 
         end where
        else
         where (k < KMT(:,:,bid))
          BUOY_FREQ_SQ(:,:,k,bid) = - grav * WORK3 * dzwr(k)
          WORK1 = dzw(k)**2 / ( WORK2 + eps2 )
          RI_NO(:,:,k) = WORK1 * BUOY_FREQ_SQ(:,:,k,bid) 
         end where
        endif

        TK  = TKP

      end do

!-----------------------------------------------------------------------
!
!     compute the first baroclinic gravity-wave phase speed.
!     Computation of Rossby deformation radius follows Chelton et al.(1998)
!
!-----------------------------------------------------------------------

      C_ROSSBY = c0

      k = 1
      if (partial_bottom_cells) then
       where ( k < KMT(:,:,bid) )
        C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k-1,bid)
       endwhere

       do k=1,km
        where ( k < KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k,bid)
        endwhere
        where ( k > 1  .and.  k == KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k-1,bid)) * DZTW(:,:,k,bid)
        endwhere
       enddo
      else
       where ( k < KMT(:,:,bid) )
        C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k-1)
       endwhere

       do k=1,km
        where ( k < KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k)
        endwhere
        where ( k > 1  .and.  k == KMT(:,:,bid) )
          C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k-1,bid)) * dzw(k)
        endwhere
       enddo
      endif

      C_ROSSBY = C_ROSSBY / pi

      L_ROSSBY = min( C_ROSSBY / (abs(FCORT(:,:,bid))+eps), &
                      sqrt( C_ROSSBY / (c2*BTP(:,:,bid)) ) )

!-----------------------------------------------------------------------
!
!     compute the inverse time scale
!
!-----------------------------------------------------------------------

      do k=1,km-1
        WORK1 = max( abs( FCORT(:,:,bid) ), sqrt(C_ROSSBY * c2 * BTP(:,:,bid)) )
        where (k < KMT(:,:,bid))
          SIGMA(:,:,k) = SIGMA_TOPO_MASK(:,:,k,bid) * WORK1  &
                        / sqrt( RI_NO(:,:,k) + gamma_eg )
        end where
      enddo

!----------------------------------------------------------------------
!
!     compute KAPPA_ISOP = c_kappa * sigma * min( Rossby_radius, Rhines scales)^2.
!     note that KAPPA_ISOP is stored in KAPPA_VERTICAL. KAPPA_VERTICAL is
!     located at the T-grid points. in the following, KAPPA_ISOP computed at
!     the lower interface of a T-grid box is assigned to the T-grid point 
!     just above this interface for simplicity.
!
!-----------------------------------------------------------------------

      do k=1,km

        WORK1 = min(L_ROSSBY, SIGMA(:,:,k)/BTP(:,:,bid))

        KAPPA_VERTICAL(:,:,k,bid) = const_eg * SIGMA(:,:,k) * WORK1**2 

      enddo

!----------------------------------------------------------------------
!
!     use below-diabatic-layer values of KAPPA_ISOP within the surface 
!     diabatic layer
!
!----------------------------------------------------------------------

      WORK1 = BL_DEPTH(:,:,bid)
      if ( transition_layer_on )  &
        WORK1 = TLT%DIABATIC_DEPTH(:,:,bid)

      do k=km-1,1,-1
        where ( zw(k) <= WORK1 )
          KAPPA_VERTICAL(:,:,k,bid) = KAPPA_VERTICAL(:,:,k+1,bid)
        endwhere
      enddo

!----------------------------------------------------------------------
!
!     impose lower and upper limits on KAPPA
!
!----------------------------------------------------------------------

      KAPPA_VERTICAL(:,:,:,bid) = max(kappa_min_eg, KAPPA_VERTICAL(:,:,:,bid))
      KAPPA_VERTICAL(:,:,:,bid) = min(kappa_max_eg, KAPPA_VERTICAL(:,:,:,bid))

      end subroutine kappa_eg

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_hdgr 
! !INTERFACE:

      subroutine kappa_lon_lat_hdgr (TMIX, this_block)

! !DESCRIPTION:
!  Variable kappa parameterization based on the vertically averaged
!  horizontal density gradients as a measure of baroclinicity. This
!  approach is from GFDL. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k_min, k_max,  & ! min and max vertical indices for the vertical average
         k,             & ! vertical loop index
         i, j,          & ! horizontal loop indices
         kk,            & ! temporary value of KMT
         bid              ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         KMASKE, KMASKN, &  ! ocean masks for the east and north faces of a T-cell
         WORK,           &  ! work array
         RHOT,           &  ! dRHO/dT
         RHOS,           &  ! dRHO/dS
         RHO_X,          &  ! grid x-direction density difference
         RHO_Y,          &  ! grid y-direction density difference
         RHO_XY_AVG         ! vertically averaged horizontal density gradient

      real (r8), dimension(nx_block,ny_block,2) :: &
         TRACER_X,       &  ! grid x-direction T and S differences
         TRACER_Y           ! grid y-direction T and S differences

      real (r8) :: &
         depth_min, depth_max, & ! actual min and max depths (reflecting model
                                 !  discretization) for computing vertical
                                 !  averages at every horizontal grid point
         depth_w_min,    &       ! protection against k_min = k_max = 0 in
         depth_w_max,    &       !  array zw
         scaling_coefficient     ! multiplicative factor

!-----------------------------------------------------------------------
!
!     parameter values for the formulation 
!
!-----------------------------------------------------------------------

      real (r8), parameter :: &
         tuning_const      = 2.0_r8,    & ! dimensionless tuning constant 
         length_scale_ref  = 50.0e5_r8, & ! constant length scale in cm
         buoyancy_freq_ref = 4.0e-3_r8, & ! constant buoyancy frequency in 1/s
         vert_average_min  = 1.0e4_r8,  & ! min and max depth range for the
         vert_average_max  = 2.0e5_r8,  & !  vertical average in cm
         kappa_min         = 3.0e6_r8,  & ! min KAPPA_LATERAL value in cm^2/s
         kappa_max         = 4.0e7_r8     ! max KAPPA_LATERAL value in cm^2/s

!-----------------------------------------------------------------------
!
!     initialization
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      k_min = -1 
      k_max =  0

      TRACER_X   = c0
      TRACER_Y   = c0
      RHO_X      = c0
      RHO_Y      = c0
      RHO_XY_AVG = c0

      scaling_coefficient = tuning_const * (length_scale_ref**2) * grav &
                           / ( rho_sw * buoyancy_freq_ref ) 

      vert_average_loop: do k=1,km

        if ( zt(k) >= vert_average_min  .and. &
             zt(k) <= vert_average_max ) then

!-----------------------------------------------------------------------
!
!     start computing the vertical average of the horizontal density 
!     gradients 
!
!-----------------------------------------------------------------------

          if ( k_min == -1 )  k_min = k-1 ! k_min is the upper limit for averaging 

          KMASKE = merge(c1, c0, k <= KMT(:,:,bid) .and. &
                         k <= KMTE(:,:,bid))
          KMASKN = merge(c1, c0, k <= KMT(:,:,bid) .and. &
                         k <= KMTN(:,:,bid))

          WORK = max(-c2, TMIX(:,:,k,1))

!-----------------------------------------------------------------------
!
!     compute tracer horizontal differences 
!
!-----------------------------------------------------------------------

          do j=1,ny_block
            do i=1,nx_block-1
              TRACER_X(i,j,1) = KMASKE(i,j) * ( WORK(i+1,j) &
                                              - WORK(i,j) )
              TRACER_X(i,j,2) = KMASKE(i,j)                 &
                          * ( TMIX(i+1,j,k,2) - TMIX(i,j,k,2) )
            enddo
          enddo

          do j=1,ny_block-1
            do i=1,nx_block
              TRACER_Y(i,j,1) = KMASKN(i,j) * ( WORK(i,j+1) &
                                              - WORK(i,j) )
              TRACER_Y(i,j,2) = KMASKN(i,j)                 &
                          * ( TMIX(i,j+1,k,2) - TMIX(i,j,k,2) )
            enddo
          enddo

!-----------------------------------------------------------------------
!
!     now compute horizontal density differences
!
!-----------------------------------------------------------------------

          call state( k, k, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS )

          RHO_X = RHOT * TRACER_X(:,:,1) + RHOS * TRACER_X(:,:,2)
          RHO_Y = RHOT * TRACER_Y(:,:,1) + RHOS * TRACER_Y(:,:,2)

!-----------------------------------------------------------------------
!
!     compute horizontal density gradients and perform their vertical 
!     integral 
!
!-----------------------------------------------------------------------

          if (partial_bottom_cells) then
            do j=2,ny_block-1
              do i=2,nx_block-1
                RHO_XY_AVG(i,j) = RHO_XY_AVG(i,j) + DZT(i,j,k,bid) * sqrt( p5 &
                         * ( ( RHO_X(i,j)**2 + RHO_X(i-1,j)**2 )     &
                            * DXTR(i,j,bid)**2 )                         & 
                         + ( ( RHO_Y(i,j)**2 + RHO_Y(i,j-1)**2 )     &  
                            * DYTR(i,j,bid)**2 ) ) 
              enddo
            enddo
          else
            do j=2,ny_block-1
              do i=2,nx_block-1
                RHO_XY_AVG(i,j) = RHO_XY_AVG(i,j) + dz(k) * sqrt( p5 &
                         * ( ( RHO_X(i,j)**2 + RHO_X(i-1,j)**2 )     &
                            * DXTR(i,j,bid)**2 )                         & 
                         + ( ( RHO_Y(i,j)**2 + RHO_Y(i,j-1)**2 )     &  
                            * DYTR(i,j,bid)**2 ) ) 
              enddo
            enddo
          endif
        else

          if ( k_min >= 0  .and.  k_max == 0 ) then
            k_max = k-1            ! k_max is the lower limit for averaging 
            exit vert_average_loop 
          endif

        endif 

      enddo vert_average_loop 

!-----------------------------------------------------------------------
!
!     compute the vertical average of horizontal density gradients for 
!     the specified vertical range
!
!-----------------------------------------------------------------------

      depth_w_min = c0
      depth_w_max = c0
      if ( k_min  > 0 )  depth_w_min = zw(k_min) 
      if ( k_max /= 0 )  depth_w_max = zw(k_max)

      do j=2,ny_block-1
        do i=2,nx_block-1
          kk = KMT(i,j,bid)
          if ( kk > 0 ) then
            depth_min = min( depth_w_min, zw(kk) )
            depth_max = min( depth_w_max, zw(kk) )
            RHO_XY_AVG(i,j) = RHO_XY_AVG(i,j)        &
                             / ( depth_max - depth_min + eps )
          else
            RHO_XY_AVG(i,j) = c0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!
!     set KAPPA_LATERAL and enforce the limits 
!
!-----------------------------------------------------------------------

      KAPPA_LATERAL(:,:,bid) = scaling_coefficient * RHO_XY_AVG

      KAPPA_LATERAL(:,:,bid) = min( kappa_max, KAPPA_LATERAL(:,:,bid) )
      KAPPA_LATERAL(:,:,bid) = max( kappa_min, KAPPA_LATERAL(:,:,bid) )

      where ( KMT(:,:,bid) <= k_min )        ! KAPPA_LATERAL = kappa_min 
        KAPPA_LATERAL(:,:,bid) = kappa_min   !  when depth < vert_average_min
      endwhere

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_hdgr

!***********************************************************************
!BOP
! !IROUTINE: kappa_lon_lat_dradius 
! !INTERFACE:

      subroutine kappa_lon_lat_dradius (this_block)

! !DESCRIPTION:
!  Variable kappa parameterization based on
!  \begin{equation}
!      KAPPA_LATERAL = KAPPA_REF {LENGTH_SCALE \over LENGTH_SCALE_REF}
!  \end{equation}
!  where $LENGTH_SCALE=$ min (Rossby deformation radius, sqrt(TAREA)).
!  Computation of Rossby deformation radius follows Chelton et al.(1998)
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,            &    ! vertical loop index
         bid                ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block) :: &
         WORK1,        &    ! work array
         LENGTH_SCALE       ! mostly used as a temporary array,
                            !  but the final content is a length scale

      real (r8), parameter :: &
         length_scale_ref = 15.0e+5_r8, &  ! reference length scale in cm
         kappa_lateral_min = 0.1e+7_r8, &  ! minimum kappa in cm^2/s
         kappa_lateral_max = 8.0e+7_r8     ! maximum kappa in cm^2/s

!-----------------------------------------------------------------------
!
!     compute the first baroclinic gravity-wave phase speed, considering
!     the full-depth integral of N (consider pi later)
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      LENGTH_SCALE = c0

      k = 1
      if (partial_bottom_cells) then
       where ( k < KMT(:,:,bid) )
         LENGTH_SCALE = LENGTH_SCALE  &
                       + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k-1,bid)
       endwhere

       do k=1,km
         where ( k < KMT(:,:,bid) )
           LENGTH_SCALE = LENGTH_SCALE  &
                         + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k,bid)
         endwhere
         where ( k == KMT(:,:,bid) )
           LENGTH_SCALE = LENGTH_SCALE  &
                         + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * p5 * DZT(:,:,k,bid)
         endwhere
       enddo
      else
       where ( k < KMT(:,:,bid) )
         LENGTH_SCALE = LENGTH_SCALE  &
                       + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k-1)
       endwhere

       do k=1,km
         where ( k < KMT(:,:,bid) )
           LENGTH_SCALE = LENGTH_SCALE  &
                         + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k)
         endwhere
         where ( k == KMT(:,:,bid) )
           LENGTH_SCALE = LENGTH_SCALE  &
                         + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * p5 * dz(k)
         endwhere
       enddo
      endif

!-----------------------------------------------------------------------
!
!     now compute the Rossby radius of deformation
!
!-----------------------------------------------------------------------

      where ( abs(TLAT(:,:,bid)) > c5 / radian )
        LENGTH_SCALE = LENGTH_SCALE / ( pi * abs(FCORT(:,:,bid)) )
      elsewhere
        LENGTH_SCALE = sqrt( LENGTH_SCALE / ( c2 * pi * BTP(:,:,bid) ) )
      endwhere

!-----------------------------------------------------------------------
!
!     consider grid size as a possible lower limit of length scale.
!     note that if the vertical integral of BUOY_FREQ_SQ is zero, then
!     LENGTH_SCALE = 0. we choose to enforce limits on KAPPA later, rather
!     then on LENGTH_SCALE below.
!
!-----------------------------------------------------------------------

      WORK1 = min( DXT(:,:,bid), DYT(:,:,bid) )
      LENGTH_SCALE = min( LENGTH_SCALE, WORK1 )

!-----------------------------------------------------------------------
!
!     now compute KAPPA_LATERAL
!
!-----------------------------------------------------------------------

      KAPPA_LATERAL(:,:,bid) = ah * LENGTH_SCALE / length_scale_ref

      if ( kappa_isop_type == kappa_type_const )          &
        KAPPA_LATERAL(:,:,bid) = ah_bolus * LENGTH_SCALE  &
                                / length_scale_ref

      KAPPA_LATERAL(:,:,bid) = min(KAPPA_LATERAL(:,:,bid), &
                                   kappa_lateral_max)
      KAPPA_LATERAL(:,:,bid) = max(KAPPA_LATERAL(:,:,bid), &
                                   kappa_lateral_min)

!-----------------------------------------------------------------------
!EOC

      end subroutine kappa_lon_lat_dradius

!***********************************************************************
!BOP
! !IROUTINE: buoyancy_frequency_dependent_profile 
! !INTERFACE:

      subroutine buoyancy_frequency_dependent_profile (TMIX, this_block)

! !DESCRIPTION:
!  Computes normalized buoyancy frequency ($N$) dependent profiles that 
!  are used to represent the vertical variation of the isopycnal and/or
!  thickness diffusion coefficients. We use 
!  \begin{equation}
!   KAPPA_VERTICAL = {{N^2}\over {N_ref^2}}, 
!  \end{equation}
!  where $N_ref$ is the value of $N$ "just" below the surface diabatic
!  layer (SDL) provided that $N^2 > 0$ there. Otherwise, $N_ref$ is the
!  first stable $N^2$ below SDL. KAPPA_VERTICAL is set to unity for shallower
!  depths. Also, 0.1 <= KAPPA_VERTICAL <= 1.0 is enforced. Note that
!  this implies KAPPA_VERTICAL = 0.1 if the local $N^2 \le 0$.
!  KAPPA_VERTICAL is located at the model tracer points.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX                  ! tracers at all vertical levels
                               !   at mixing time level

      type (block), intent(in) :: &
         this_block            ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,        &          ! vertical loop index
         bid                  ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_MIN                ! k index below SDL 

      real (r8), dimension(nx_block,ny_block) :: &
         TEMP_K,           &  ! temperature at level k
         TEMP_KP1,         &  ! temperature at level k+1
         RHOT,             &  ! dRHO/dT
         RHOS,             &  ! dRHO/dS
         BUOY_FREQ_SQ_REF, &  ! reference (normalization) value
                              !  of N^2
         SDL                  ! surface diabatic layer (see below)

      real (r8), dimension(nx_block,ny_block,km) :: &
         BUOY_FREQ_SQ_NORM    ! normalized N^2 defined at level interfaces

!-----------------------------------------------------------------------
!
!     initialize variables
!
!     note that "km+1" value in K_MIN is used as a flag, i.e. it indicates
!     that a minimum value has not been found yet or that it does not
!     exist.
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      BUOY_FREQ_SQ_NORM         = c0
      BUOY_FREQ_SQ_REF          = c0
      KAPPA_VERTICAL(:,:,:,bid) = c1

      K_MIN = merge( km+1, 0, KMT(:,:,bid) /= 0 )  

      SDL = zw(1)
      if ( vmix_itype == vmix_type_kpp )  SDL = KPP_HBLT(:,:,bid) 
      if ( transition_layer_on )    SDL = TLT%INTERIOR_DEPTH(:,:,bid)

      if ( .not. read_n2_data ) then

!-----------------------------------------------------------------------
!
!     compute N^2
!
!-----------------------------------------------------------------------

        do k=1,km-1

          if ( k == 1 ) &
            TEMP_K = max( -c2, TMIX(:,:,k,1) )

          TEMP_KP1 = max( -c2, TMIX(:,:,k+1,1) )

          call state( k, k+1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                      this_block, DRHODT=RHOT, DRHODS=RHOS )

          if (partial_bottom_cells) then
           where ( k < KMT(:,:,bid) ) 
             BUOY_FREQ_SQ(:,:,k,bid) = max( c0, - grav / DZTW(:,:,k,bid) &
                 * ( RHOT * ( TEMP_K - TEMP_KP1 )                &
                   + RHOS * ( TMIX(:,:,k,  2) - TMIX(:,:,k+1,2) ) ) )
           endwhere
          else
           where ( k < KMT(:,:,bid) ) 
             BUOY_FREQ_SQ(:,:,k,bid) = max( c0, - grav * dzwr(k) &
                 * ( RHOT * ( TEMP_K - TEMP_KP1 )                &
                   + RHOS * ( TMIX(:,:,k,  2) - TMIX(:,:,k+1,2) ) ) )
           endwhere
          endif

          TEMP_K = TEMP_KP1

        enddo

      endif

!-----------------------------------------------------------------------
!
!     determine the reference buoyancy frequency and the associated 
!     level index (most likely just below SDL)
!
!-----------------------------------------------------------------------

      do k=1,km-1
        where ( ( K_MIN == km+1 ) .and. ( zw(k) > SDL ) .and.  &
                ( k <= KMT(:,:,bid) )  .and.                   &
                ( BUOY_FREQ_SQ(:,:,k,bid) > c0 ) ) 
          BUOY_FREQ_SQ_REF = BUOY_FREQ_SQ(:,:,k,bid)
          K_MIN = k
        endwhere
      enddo
        
!-----------------------------------------------------------------------
!
!     now compute the normalized profiles at the interfaces 
!
!-----------------------------------------------------------------------

      do k=1,km-1
        where ( ( k >= K_MIN ) .and. ( k < KMT(:,:,bid) ) .and. &
               ( BUOY_FREQ_SQ_REF /= c0 ) )
          BUOY_FREQ_SQ_NORM(:,:,k) =  &
              max( BUOY_FREQ_SQ(:,:,k,bid) / BUOY_FREQ_SQ_REF, 0.1_r8 )
          BUOY_FREQ_SQ_NORM(:,:,k) =  &
              min( BUOY_FREQ_SQ_NORM(:,:,k), c1 ) 
        elsewhere
          BUOY_FREQ_SQ_NORM(:,:,k) = c1
        endwhere
      enddo

      do k=1,km-1
        where ( k == KMT(:,:,bid)-1 ) 
          BUOY_FREQ_SQ_NORM(:,:,k+1) = BUOY_FREQ_SQ_NORM(:,:,k)
        endwhere
      enddo

!-----------------------------------------------------------------------
!
!     finally, do not average interface values of BUOY_FREQ_SQ_NORM to
!     the tracer points, but instead copy from above to preserve the 
!     extrema. 
!
!-----------------------------------------------------------------------

      do k=2,km
        where ( ( k > K_MIN ) .and. ( k <= KMT(:,:,bid) ) )
          KAPPA_VERTICAL(:,:,k,bid) = BUOY_FREQ_SQ_NORM(:,:,k-1)
        endwhere
      enddo

!-----------------------------------------------------------------------
!EOC

      end subroutine buoyancy_frequency_dependent_profile

!***********************************************************************
!BOP
! !IROUTINE: transition_layer 
! !INTERFACE:

      subroutine transition_layer ( this_block )

! !DESCRIPTION:
!  Compute transition layer related fields. the related algorithms
!  should work even with zero transition layer depth.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block

      integer (int_kind), dimension(nx_block,ny_block) :: &
         K_START,   &        ! work arrays for TLT%K_LEVEL and 
         K_SUB               !  TLT%ZTW, respectively

      logical (log_kind), dimension(nx_block,ny_block) :: &
         COMPUTE_TLT         ! flag

      real (r8), dimension(nx_block,ny_block) :: &
         WORK                ! work space for TLT%THICKNESS

      real (r8), dimension(2) :: &
         reference_depth     ! zt or zw

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      K_START = 0
      K_SUB   = 0

      TLT%THICKNESS(:,:,bid)      = c0
      TLT%INTERIOR_DEPTH(:,:,bid) = c0
      TLT%K_LEVEL(:,:,bid)        = 0
      TLT%ZTW(:,:,bid)            = 0

      COMPUTE_TLT = merge(.true., .false., KMT(:,:,bid) /= 0)

!-----------------------------------------------------------------------
!
!     initial pass to determine the minimum transition layer thickness.
!     the vertical level at or below which the interior region starts
!     is also computed.
!
!-----------------------------------------------------------------------

      do k=1,km

        where ( COMPUTE_TLT  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) < zw(k) )

          K_START = k+1
          K_SUB   = ktp

          TLT%THICKNESS(:,:,bid) = zw(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 2

          COMPUTE_TLT = .false.

        endwhere

        where ( k /= 1  .and.  K_START == k+1  .and. &
                TLT%DIABATIC_DEPTH(:,:,bid) < zt(k) )

          K_START = k
          K_SUB   = kbt

          TLT%THICKNESS(:,:,bid) = zt(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 1

        endwhere

      enddo

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect DIABATIC_DEPTH value in TLT'  &
                        //  ' computation')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!
!     using R * |S| as the vertical displacement length scale associated
!     with a horizontal displacement equivalent to the Rossby deformation
!     radius, R, determine if the transition layer thickness needs to be
!     modified. |S| is the absolute value of the slope of the isopycnal
!     surfaces. First, start with columns where K_SUB = kbt.
!
!-----------------------------------------------------------------------

      where ( KMT(:,:,bid) == 0  .or.  K_START > KMT(:,:,bid)  .or.  &
              ( K_START == KMT(:,:,bid)  .and.  K_SUB == kbt ) )
        COMPUTE_TLT = .false.
      elsewhere
        COMPUTE_TLT = .true.
      endwhere

      do k=1,km-1

        WORK = c0

        where ( COMPUTE_TLT  .and.  K_SUB == kbt  .and.  &
                K_START < KMT(:,:,bid)  .and.  K_START == k )
          WORK = max(SLA_SAVE(:,:,kbt,k,bid), &
                     SLA_SAVE(:,:,ktp,k+1,bid)) * RB(:,:,bid)
        endwhere

        where ( WORK /= c0  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) <  (zw(k) - WORK) )
          COMPUTE_TLT = .false.
        endwhere

        where ( WORK /= c0  .and.  &
                TLT%DIABATIC_DEPTH(:,:,bid) >= (zw(k) - WORK) )

          K_START = K_START + 1
          K_SUB   = ktp

          TLT%THICKNESS(:,:,bid) = zw(k) - TLT%DIABATIC_DEPTH(:,:,bid)
          TLT%K_LEVEL(:,:,bid)   = k
          TLT%ZTW(:,:,bid)       = 2

        endwhere

      enddo

!-----------------------------------------------------------------------
!
!     now consider the deeper levels
!
!-----------------------------------------------------------------------

      do k=2,km

        reference_depth(ktp) = zt(k)
        reference_depth(kbt) = zw(k)

        do kk=ktp,kbt

          WORK = c0

          if (kk == ktp) then
            where ( COMPUTE_TLT  .and.  K_START <= KMT(:,:,bid)  .and. &
                    K_START == k )
              WORK = max(SLA_SAVE(:,:,ktp,k,bid), &
                         SLA_SAVE(:,:,kbt,k,bid)) * RB(:,:,bid)
            endwhere
          else
            where ( COMPUTE_TLT  .and.  K_START < KMT(:,:,bid)  .and. &
                    K_START == k )
              WORK = max(SLA_SAVE(:,:,kbt,k,bid), &
                         SLA_SAVE(:,:,ktp,k+1,bid)) * RB(:,:,bid)
            endwhere
            where ( COMPUTE_TLT  .and.  K_START == KMT(:,:,bid)  .and. &
                    K_START == k )
              WORK = SLA_SAVE(:,:,kbt,k,bid) * RB(:,:,bid)
            endwhere
          endif

          where ( WORK /= c0  .and.  &
           TLT%DIABATIC_DEPTH(:,:,bid) <  (reference_depth(kk) - WORK) )
            COMPUTE_TLT = .false.
          endwhere

          where ( WORK /= c0  .and.  &
           TLT%DIABATIC_DEPTH(:,:,bid) >= (reference_depth(kk) - WORK) )
            TLT%THICKNESS(:,:,bid) = reference_depth(kk)  &
                                    - TLT%DIABATIC_DEPTH(:,:,bid)
            TLT%K_LEVEL(:,:,bid)   = k
            TLT%ZTW(:,:,bid)       = kk
          endwhere

        enddo

        where ( COMPUTE_TLT  .and.  K_START == k )
          K_START = K_START + 1
        endwhere

      enddo

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect TLT computations')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!
!     compute the depth at which the interior, adiabatic region starts
!
!-----------------------------------------------------------------------

      do k=1,km
        where ( TLT%K_LEVEL(:,:,bid) == k  .and.  &
                TLT%ZTW(:,:,bid) == 1 )
          TLT%INTERIOR_DEPTH(:,:,bid) = zt(k)
        endwhere
        where ( TLT%K_LEVEL(:,:,bid) == k  .and.  &
                TLT%ZTW(:,:,bid) == 2 )
          TLT%INTERIOR_DEPTH(:,:,bid) = zw(k)
        endwhere
      enddo

      COMPUTE_TLT = .false.
      where ( KMT(:,:,bid) /= 0  .and.  &
              TLT%INTERIOR_DEPTH(:,:,bid) == c0 )  &
        COMPUTE_TLT = .true.
      where ( KMT(:,:,bid) == 0  .and.  &
              TLT%INTERIOR_DEPTH(:,:,bid) /= c0 )  &
        COMPUTE_TLT = .true.

#ifdef CCSMCOUPLED
#ifndef _HIRES
      if ( any(COMPUTE_TLT) ) then
        call shr_sys_abort ('Incorrect TLT%INTERIOR_DEPTH computation')
      endif
#endif
#endif

!-----------------------------------------------------------------------
!EOC

      end subroutine transition_layer

!***********************************************************************
!BOP
! !IROUTINE: merged_streamfunction 
! !INTERFACE:


      subroutine merged_streamfunction ( this_block )

! !DESCRIPTION:
!  Construct a merged streamfunction that has the appropriate
!  behavior in the surface diabatic region, transition layer, and
!  adiabatic interior
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block
 
      real (r8), dimension(nx_block,ny_block,2) :: &
         WORKD1, WORKD2, WORKD3, WORKD4, & !SJR MOD   ! work arrays for diagonal terms
         WORKF1, WORKF2, WORKF3, WORKF4    !SJR ADD ! work arrays for off-diadgonal terms

      real (r8), dimension(nx_block,ny_block) :: &
         WORKD2_NEXT, WORKD4_NEXT, & !SJR ADD   ! work arrays for diagonal terms
         WORKF2_NEXT, WORKF4_NEXT    !SJR ADD ! work arrays for off-diadgonal terms

      real (r8), dimension(nx_block,ny_block) :: &
         WORK5, WORK6, WORK7          ! more work arrays

      logical (log_kind), dimension(nx_block,ny_block) :: &
         LMASK               ! flag

      real (r8), dimension(2) :: &
         reference_depth              ! zt or zw

!-----------------------------------------------------------------------
!
!     initialize various quantities
!
!-----------------------------------------------------------------------

      bid = this_block%local_id

      SF_SLXX(:,:,:,:,:,bid) = c0 !SJR MOD
      SF_SLYY(:,:,:,:,:,bid) = c0 !SJR MOD
      SF_SLXY(:,:,:,:,:,bid) = c0 !SJR ADD
      SF_SLYX(:,:,:,:,:,bid) = c0 !SJR ADD

      WORKD1 = c0 !SJR MOD
      WORKD2 = c0 !SJR MOD
      WORKD3 = c0 !SJR MOD
      WORKD4 = c0 !SJR MOD
      WORKF1 = c0 !SJR ADD
      WORKF2 = c0 !SJR ADD
      WORKF3 = c0 !SJR ADD
      WORKF4 = c0 !SJR ADD

      WORK5 = c0
      WORK6 = c0
      WORK7 = c0

      WORKD2_NEXT = c0
      WORKD4_NEXT = c0
      WORKF2_NEXT = c0
      WORKF4_NEXT = c0

!-----------------------------------------------------------------------
!
!     compute the interior streamfunction and its first derivative at the
!     INTERIOR_DEPTH level. WORK1 and WORK2 contain the streamfunction
!     and its first derivative, respectively, for the zonal component
!     for the east and west sides of a grid cell. WORK3 and WORK4 are
!     the corresponding fields for the meridional component for the
!     north and south sides of a grid cell. Note that these definitions
!     include a "dz". Also, the first derivative computations assume 
!     that the streamfunctions are located in the middle of the top or
!     bottom half of a grid cell, hence a factor of two in WORK2 and
!     WORK4 calculations. 
!
!-----------------------------------------------------------------------

      do k=1,km-1
!-------------v--SJR MOD--v--------------------------------------------- 
        do kk=1,2

          LMASK = TLT%K_LEVEL(:,:,bid) == k  .and.            &
                  TLT%K_LEVEL(:,:,bid) < KMT(:,:,bid)  .and.  &
                  TLT%ZTW(:,:,bid) == 1
       
          if (partial_bottom_cells) then
           where ( LMASK )

            WORKD1(:,:,kk) =  KXX_THIC(:,:,kbt,k,bid) * DD_SLX(:,:,kk,kbt,k,bid)
            WORKD2(:,:,kk) = c2 / DZTW(:,:,k,bid) * ( WORKD1(:,:,kk)            &
                            - KXX_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) ) 

            WORKD2_NEXT = c2 / DZT(:,:,k+1,bid) * ( &
              KXX_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) - &
              KXX_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) )

            WORKD3(:,:,kk) =  KYY_THIC(:,:,kbt,k,bid) * DD_SLY(:,:,kk,kbt,k,bid)
            WORKD4(:,:,kk) = c2 / DZTW(:,:,k,bid) * ( WORKD3(:,:,kk)            &
                            - KYY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) )

            WORKD4_NEXT = c2 / DZT(:,:,k+1,bid) * ( &
              KYY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) - &
              KYY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF1(:,:,kk) =  KXY_THIC(:,:,kbt,k,bid) * DD_SLX(:,:,kk,kbt,k,bid)
            WORKF2(:,:,kk) = c2 / DZTW(:,:,k,bid) * ( WORKF1(:,:,kk)            &
                            - KXY_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) )

            WORKF2_NEXT = c2 / DZT(:,:,k+1,bid) * ( &
              KXY_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) - &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) )

            WORKF3(:,:,kk) =  KXY_THIC(:,:,kbt,k,bid) * DD_SLY(:,:,kk,kbt,k,bid)
            WORKF4(:,:,kk) = c2 / DZTW(:,:,k,bid) * ( WORKF3(:,:,kk)            &
                            - KXY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) )

            WORKF4_NEXT = c2 / DZT(:,:,k+1,bid) * ( &
              KXY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) - &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) )

           endwhere
           endif
          else   
           where ( LMASK )

            WORKD1(:,:,kk) =  KXX_THIC(:,:,kbt,k,bid) * DD_SLX(:,:,kk,kbt,k,bid)
            WORKD2(:,:,kk) = c2 * dzwr(k) * ( WORKD1(:,:,kk)            &
                            - KXX_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) )

            WORKD2_NEXT = c2 * dzr(k+1) * ( &
              KXX_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) - &
              KXX_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) )

            WORKD3(:,:,kk) =  KYY_THIC(:,:,kbt,k,bid) * DD_SLY(:,:,kk,kbt,k,bid)
            WORKD4(:,:,kk) = c2 * dzwr(k) * ( WORKD3(:,:,kk)            &
                            - KYY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) )

            WORKD4_NEXT = c2 * dzr(k+1) * ( &
              KYY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) - &
              KYY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF1(:,:,kk) =  KXY_THIC(:,:,kbt,k,bid) * DD_SLX(:,:,kk,kbt,k,bid)
            WORKF2(:,:,kk) = c2 * dzwr(k) * ( WORKF1(:,:,kk)            &
                            - KXY_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) ) 

            WORKF2_NEXT = c2 * dzr(k+1) * ( &
              KXY_THIC(:,:,ktp,k+1,bid) * DD_SLX(:,:,kk,ktp,k+1,bid) - &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) )

            WORKF3(:,:,kk) =  KXY_THIC(:,:,kbt,k,bid) * DD_SLY(:,:,kk,kbt,k,bid)
            WORKF4(:,:,kk) = c2 * dzwr(k) * ( WORKF3(:,:,kk)            &
                            - KXY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) ) 

            WORKF4_NEXT = c2 * dzr(k+1) * ( &
              KXY_THIC(:,:,ktp,k+1,bid) * DD_SLY(:,:,kk,ktp,k+1,bid) - &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) )

           endwhere
           endif
          endif

           where ( LMASK .and. abs(WORKD2_NEXT) < abs(WORKD2(:,:,kk)) ) &
            WORKD2(:,:,kk) = WORKD2_NEXT

           where ( LMASK .and. abs(WORKD4_NEXT) < abs(WORKD4(:,:,kk)) ) &
            WORKD4(:,:,kk) = WORKD4_NEXT

           if (.not. isoonly) then
              where ( LMASK .and. abs(WORKF2_NEXT) < abs(WORKF2(:,:,kk)) ) &
               WORKF2(:,:,kk) = WORKF2_NEXT

              where ( LMASK .and. abs(WORKF4_NEXT) < abs(WORKF4(:,:,kk)) ) &
               WORKF4(:,:,kk) = WORKF4_NEXT
           endif


          LMASK = TLT%K_LEVEL(:,:,bid) == k  .and.           &
                  TLT%K_LEVEL(:,:,bid) < KMT(:,:,bid)  .and. &
                  TLT%ZTW(:,:,bid) == 2

          if (partial_bottom_cells) then
           where ( LMASK )

            WORKD1(:,:,kk) =  KXX_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLX(:,:,kk,ktp,k+1,bid)
            WORKD2(:,:,kk) =  c2 / DZT(:,:,k+1,bid) * ( WORKD1(:,:,kk)     &
                           - ( KXX_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLX(:,:,kk,kbt,k+1,bid) ) )

            WORKD3(:,:,kk) =  KYY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLY(:,:,kk,ktp,k+1,bid)
            WORKD4(:,:,kk) =  c2 / DZT(:,:,k+1,bid) * ( WORKD3(:,:,kk)     &
                           - ( KYY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLY(:,:,kk,kbt,k+1,bid) ) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF1(:,:,kk) =  KXY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLX(:,:,kk,ktp,k+1,bid)
            WORKF2(:,:,kk) =  c2 / DZT(:,:,k+1,bid) * ( WORKF1(:,:,kk)     &
                           - ( KXY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLX(:,:,kk,kbt,k+1,bid) ) )

            WORKF3(:,:,kk) =  KXY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLY(:,:,kk,ktp,k+1,bid)
            WORKF4(:,:,kk) =  c2 / DZT(:,:,k+1,bid) * ( WORKF3(:,:,kk)     &
                           - ( KXY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLY(:,:,kk,kbt,k+1,bid) ) )

           endwhere
           endif
          else
           where ( LMASK )

            WORKD1(:,:,kk) =  KXX_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLX(:,:,kk,ktp,k+1,bid)
            WORKD2(:,:,kk) =  c2 * dzr(k+1) * ( WORKD1(:,:,kk)     &
                           - ( KXX_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLX(:,:,kk,kbt,k+1,bid) ) )

            WORKD3(:,:,kk) =  KYY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLY(:,:,kk,ktp,k+1,bid)
            WORKD4(:,:,kk) =  c2 * dzr(k+1) * ( WORKD3(:,:,kk)     &
                           - ( KYY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLY(:,:,kk,kbt,k+1,bid) ) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF1(:,:,kk) =  KXY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLX(:,:,kk,ktp,k+1,bid)
            WORKF2(:,:,kk) =  c2 * dzr(k+1) * ( WORKF1(:,:,kk)     &
                           - ( KXY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLX(:,:,kk,kbt,k+1,bid) ) )

            WORKF3(:,:,kk) =  KXY_THIC(:,:,   ktp,k+1,bid)     &
                              * DD_SLY(:,:,kk,ktp,k+1,bid)
            WORKF4(:,:,kk) =  c2 * dzr(k+1) * ( WORKF3(:,:,kk)     &
                           - ( KXY_THIC(:,:,   kbt,k+1,bid)        &
                               * DD_SLY(:,:,kk,kbt,k+1,bid) ) )

           endwhere
           endif
          endif

          LMASK = LMASK .and. TLT%K_LEVEL(:,:,bid) + 1 < KMT(:,:,bid)

          if (partial_bottom_cells) then
           where ( LMASK )

            WORKD2_NEXT = c2 / DZTW(:,:,k+1,bid) * ( &
              KXX_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) - &
              KXX_THIC(:,:,ktp,k+2,bid) * DD_SLX(:,:,kk,ktp,k+2,bid) )

            WORKD4_NEXT = c2 / DZTW(:,:,k+1,bid) * ( &
              KYY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) - &
              KYY_THIC(:,:,ktp,k+2,bid) * DD_SLY(:,:,kk,ktp,k+2,bid) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF2_NEXT = c2 / DZTW(:,:,k+1,bid) * ( &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) - &
              KXY_THIC(:,:,ktp,k+2,bid) * DD_SLX(:,:,kk,ktp,k+2,bid) )

            WORKF4_NEXT = c2 / DZTW(:,:,k+1,bid) * ( &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) - &
              KXY_THIC(:,:,ktp,k+2,bid) * DD_SLY(:,:,kk,ktp,k+2,bid) )

           endwhere
           endif
          else
           where ( LMASK )

            WORKD2_NEXT = c2 * dzwr(k+1) * ( &
              KXX_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) - &
              KXX_THIC(:,:,ktp,k+2,bid) * DD_SLX(:,:,kk,ktp,k+2,bid) )

            WORKD4_NEXT = c2 * dzwr(k+1) * ( &
              KYY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) - &
              KYY_THIC(:,:,ktp,k+2,bid) * DD_SLY(:,:,kk,ktp,k+2,bid) )

           endwhere

           if (.not. isoonly) then
           where ( LMASK )

            WORKF2_NEXT = c2 * dzwr(k+1) * ( &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLX(:,:,kk,kbt,k+1,bid) - &
              KXY_THIC(:,:,ktp,k+2,bid) * DD_SLX(:,:,kk,ktp,k+2,bid) )

            WORKF4_NEXT = c2 * dzwr(k+1) * ( &
              KXY_THIC(:,:,kbt,k+1,bid) * DD_SLY(:,:,kk,kbt,k+1,bid) - &
              KXY_THIC(:,:,ktp,k+2,bid) * DD_SLY(:,:,kk,ktp,k+2,bid) )

           endwhere
           endif
          endif

           where ( LMASK .and. abs(WORKD2_NEXT) < abs(WORKD2(:,:,kk)) ) &
            WORKD2(:,:,kk) = WORKD2_NEXT

           where ( LMASK .and. abs(WORKD4_NEXT) < abs(WORKD4(:,:,kk)) ) &
            WORKD4(:,:,kk) = WORKD4_NEXT

           if (.not. isoonly) then
              where ( LMASK .and. abs(WORKF2_NEXT) < abs(WORKF2(:,:,kk)) ) &
               WORKF2(:,:,kk) = WORKF2_NEXT

              where ( LMASK .and. abs(WORKF4_NEXT) < abs(WORKF4(:,:,kk)) ) &
               WORKF4(:,:,kk) = WORKF4_NEXT
           endif
        enddo
!-------------^--SJR MOD--^--------------------------------------------- 
      enddo

!-----------------------------------------------------------------------
!
!     compute the depth independent interpolation factors used in the 
!     linear and quadratic interpolations within the diabatic and 
!     transition regions, respectively.
!
!-----------------------------------------------------------------------

      WORK5 = merge(c0,c1 / &
              ( c2 * TLT%DIABATIC_DEPTH(:,:,bid) + TLT%THICKNESS(:,:,bid) ), KMT(:,:,bid) == 0)
      WORK6 = merge(c0,WORK5 / TLT%THICKNESS(:,:,bid), TLT%THICKNESS(:,:,bid)<eps)

!-----------------------------------------------------------------------
!
!     start of interpolation to construct the merged streamfunction
!
!-----------------------------------------------------------------------

      do k=1,km

        reference_depth(ktp) = zt(k) - p25 * dz(k)
        reference_depth(kbt) = zt(k) + p25 * dz(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: use linear interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------

!-------------v--SJR MOD--v--------------------------------------------- 
            where ( reference_depth(kk) <= TLT%DIABATIC_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) ) 

               SF_SLXX(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKD1(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                        * WORKD2(:,:,1) )

               SF_SLXX(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKD1(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                        * WORKD2(:,:,2) )

               SF_SLYY(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKD3(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                        * WORKD4(:,:,1) )

               SF_SLYY(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKD3(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                        * WORKD4(:,:,2) )

            endwhere
            if (.not. isoonly) then
            where ( reference_depth(kk) <= TLT%DIABATIC_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) ) 

               SF_SLXY(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKF1(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                        * WORKF2(:,:,1) )

               SF_SLXY(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKF1(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                        * WORKF2(:,:,2) )

               SF_SLYX(:,:,1,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKF3(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                        * WORKF4(:,:,1) )

               SF_SLYX(:,:,2,kk,k,bid) = reference_depth(kk) * WORK5  &
                     * ( c2 * WORKF3(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                        * WORKF4(:,:,2) )

            endwhere
            endif
!-------------^--SJR MOD--^--------------------------------------------- 

!-----------------------------------------------------------------------
!
!     transition layer: use quadratic interpolation (in streamfunction) 
!
!-----------------------------------------------------------------------

!-------------v--SJR MOD--v--------------------------------------------- 
            where ( reference_depth(kk) > TLT%DIABATIC_DEPTH(:,:,bid)   &
               .and.  reference_depth(kk) <= TLT%INTERIOR_DEPTH(:,:,bid) &
               .and.  k <= KMT(:,:,bid) )

               WORK7 = (TLT%DIABATIC_DEPTH(:,:,bid)  &
                        - reference_depth(kk))**2

               SF_SLXX(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKD1(:,:,1) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKD2(:,:,1) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKD1(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                      * WORKD2(:,:,1) )

               SF_SLXX(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKD1(:,:,2) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKD2(:,:,2) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKD1(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                      * WORKD2(:,:,2) )

               SF_SLYY(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKD3(:,:,1) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKD4(:,:,1) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKD3(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                      * WORKD4(:,:,1) )

               SF_SLYY(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKD3(:,:,2) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKD4(:,:,2) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKD3(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                      * WORKD4(:,:,2) )

            endwhere
            if (.not. isoonly) then
            where ( reference_depth(kk) > TLT%DIABATIC_DEPTH(:,:,bid)   &
               .and.  reference_depth(kk) <= TLT%INTERIOR_DEPTH(:,:,bid) &
               .and.  k <= KMT(:,:,bid) )

               WORK7 = (TLT%DIABATIC_DEPTH(:,:,bid)  &
                        - reference_depth(kk))**2

               SF_SLXY(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKF1(:,:,1) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKF2(:,:,1) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKF1(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                      * WORKF2(:,:,1) )

               SF_SLXY(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKF1(:,:,2) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKF2(:,:,2) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKF1(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                      * WORKF2(:,:,2) )

               SF_SLYX(:,:,1,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKF3(:,:,1) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKF4(:,:,1) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKF3(:,:,1) - TLT%THICKNESS(:,:,bid)  &
                      * WORKF4(:,:,1) )

               SF_SLYX(:,:,2,kk,k,bid) = - WORK7 * WORK6            &
                   * ( WORKF3(:,:,2) - TLT%INTERIOR_DEPTH(:,:,bid)  &
                      * WORKF4(:,:,2) )                             &
                  + reference_depth(kk) * WORK5                     &
                   * ( c2 * WORKF3(:,:,2) - TLT%THICKNESS(:,:,bid)  &
                      * WORKF4(:,:,2) )

            endwhere
            endif
!-------------^--SJR MOD--^--------------------------------------------- 

!-----------------------------------------------------------------------
!
!     interior, adiabatic region: no interpolation is needed. note that
!     "dzw" is introduced here, too, for consistency. 
!
!-----------------------------------------------------------------------

!-------------v--SJR MOD--v--------------------------------------------- 
            where ( reference_depth(kk) > TLT%INTERIOR_DEPTH(:,:,bid)  & 
                  .and.  k <= KMT(:,:,bid) )

               SF_SLXX(:,:,1,kk,k,bid) =  KXX_THIC(:,:,kk,k,bid) * DD_SLX(:,:,1,kk,k,bid)
               SF_SLXX(:,:,2,kk,k,bid) =  KXX_THIC(:,:,kk,k,bid) * DD_SLX(:,:,2,kk,k,bid)

               SF_SLYY(:,:,1,kk,k,bid) =  KYY_THIC(:,:,kk,k,bid) * DD_SLY(:,:,1,kk,k,bid)
               SF_SLYY(:,:,2,kk,k,bid) =  KYY_THIC(:,:,kk,k,bid) * DD_SLY(:,:,2,kk,k,bid)

            endwhere
            if (.not. isoonly) then
            where ( reference_depth(kk) > TLT%INTERIOR_DEPTH(:,:,bid)  & 
                  .and.  k <= KMT(:,:,bid) )

               SF_SLXY(:,:,1,kk,k,bid) =  KXY_THIC(:,:,kk,k,bid) * DD_SLX(:,:,1,kk,k,bid)
               SF_SLXY(:,:,2,kk,k,bid) =  KXY_THIC(:,:,kk,k,bid) * DD_SLX(:,:,2,kk,k,bid)

               SF_SLYX(:,:,1,kk,k,bid) =  KXY_THIC(:,:,kk,k,bid) * DD_SLY(:,:,1,kk,k,bid)
               SF_SLYX(:,:,2,kk,k,bid) =  KXY_THIC(:,:,kk,k,bid) * DD_SLY(:,:,2,kk,k,bid)

            endwhere
            endif
!-------------^--SJR MOD--^--------------------------------------------- 
        enddo  ! end of kk-loop
      enddo    ! end of k-loop

!-----------------------------------------------------------------------
!EOC

      end subroutine merged_streamfunction

!***********************************************************************
!BOP
! !IROUTINE: apply_vertical_profile_to_isop_hor_diff 
! !INTERFACE:

      subroutine apply_vertical_profile_to_isop_hor_diff ( this_block ) 

! !DESCRIPTION:
!  Apply vertical tapers to KAPPA_ISOP and HOR_DIFF based on their
!  vertical location with respect to the diabatic, transition, and
!  adiabatic regions.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

      type (block), intent(in) :: &
         this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, kk,     &        ! loop indices
         bid                 ! local block address for this sub block

      real (r8), dimension(2) :: &
         reference_depth

      bid = this_block%local_id

!-----------------------------------------------------------------------
!
!     start of tapering
!
!-----------------------------------------------------------------------

      do k=1,km

        reference_depth(ktp) = zt(k) - p25 * dz(k)
        reference_depth(kbt) = zt(k) + p25 * dz(k)

        do kk=ktp,kbt

!-----------------------------------------------------------------------
!
!     diabatic region: no isopycnal diffusion 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) <= TLT%DIABATIC_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) ) 
            KAPPA_ISOP(:,:,kk,k,bid) = c0
          endwhere

!-----------------------------------------------------------------------
!
!      transition layer: a linear combination of isopcynal and horizontal
!      diffusion coefficients 
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%DIABATIC_DEPTH(:,:,bid)   &
            .and.  reference_depth(kk) <= TLT%INTERIOR_DEPTH(:,:,bid) &
            .and.  k <= KMT(:,:,bid)  .and.                           &
                   TLT%THICKNESS(:,:,bid) > eps )

            HOR_DIFF(:,:,kk,k,bid) = ( TLT%INTERIOR_DEPTH(:,:,bid)  &
                  - reference_depth(kk) ) * HOR_DIFF(:,:,kk,k,bid)  &
                  / TLT%THICKNESS(:,:,bid)
            KAPPA_ISOP(:,:,kk,k,bid) = ( reference_depth(kk)        &
                                   - TLT%DIABATIC_DEPTH(:,:,bid) )  &
                 * KAPPA_ISOP(:,:,kk,k,bid) / TLT%THICKNESS(:,:,bid)

          endwhere

!-----------------------------------------------------------------------
!
!     interior region: no horizontal diffusion
!
!-----------------------------------------------------------------------

          where ( reference_depth(kk) > TLT%INTERIOR_DEPTH(:,:,bid)  &
                  .and.  k <= KMT(:,:,bid) )
            HOR_DIFF(:,:,kk,k,bid) = c0
          endwhere

        enddo  ! end of kk-loop

      enddo    ! end of k-loop

!-----------------------------------------------------------------------
!EOC

      end subroutine apply_vertical_profile_to_isop_hor_diff

!***********************************************************************
!BOP
! !IROUTINE: shear_dispersion 
! !INTERFACE:

      subroutine shear_dispersion (UMIX, VMIX, TMIX, this_block)

! !DESCRIPTION:
!
! !REVISION HISTORY:

! !INPUT PARAMETERS:


      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         UMIX, VMIX         ! U,V  at all vertical levels
                            !  at mixing time level

      real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
         TMIX   ! T  at all vertical levels
                            !  at mixing time level

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!     shear_method: (0) MAJOR = kappa+1/kappa*<(u*dy)^2+(v*dx)^2>, flow aligned
!                   (1) MAJOR = kappa+1/kappa*[<u*dy>^2+<v*dx>^2], flow aligned
!                   (2) MAJOR = kappa+1/kappa*[<(u*dy)^2>+<(v*dx)^2>], flow aligned
!                   (3) MAJOR = kappa+1/kappa*<(|u*dy|+|v*dx|)^2>, flow aligned 
!                   (4) K_tensor = kappa*[1 0; 0 1]+1/kappa*[<u*dy>^2 <u*dy><v*dx>; <u*dy><v*dx> <v*dx>^2]
!                   (5) K_tensor = kappa*[1 0; 0 1]+1/kappa*[<(u*dy)^2> <u*dy*v*dx>; <u*dy*v*dx> <(v*dx)^2>]
!
!     grid_type: treatment for dx and dy for flux calculation across flow (Dx and Dy are grid spacings)
!                   (0) dx=Dx, dy=Dy
!                   (1) dx=dy=R_1 (first baroclinic Rossby deformation radius)
!                   (2) dx=dy=Delta/sqrt(2) (where Delta is actual distance cutting across T-cell perpendicular to velocity)
!                   (3) dx=dy=sqrt(Dx*Dy)
!                   (4) dx=Dx or dy=Dy and the other is reduced to match actual distance cut across T-cell, across flow
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k, k_sub, kp1,    &
         shear_method, &
         area_type, &  ! type of area handling - 0=just do ugrid_to_tgrid, 1=UAREA, 2=SUBCELLV
         grid_type, &
         bid                ! local block address for this sub block


      integer (int_kind), dimension(4) :: &
         stencil_type      ! type of stencil added to ME for averaging - 
                           ! index 1=+E+W, 2=+N+S, 3=+NE+SW, 4=+NW+SE

      real (r8), dimension(nx_block,ny_block,km) :: &
         UAVG, VAVG, UFLX, VFLX, &
         NXLOC, NYLOC, &
         RNLOC, WORK, KMAJ, WORKU, WORKV, WORKW, &
         LDX, LDY

      real (r8), dimension(nx_block,ny_block) :: &
         TK, TKP,                                &  ! level k and k+1 TMIX
         WORK1, WORK2, WORK3,                    &  ! work arrays
         C_ROSSBY,                               &  ! first baroclinic Rossby wave speed
         L_ROSSBY,                               &  ! Rossby deformation radius
         RHOT, RHOS,                             &   ! dRHOdT and dRHOdS
         RDLOC


!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

      bid = this_block%local_id
      
      stencil_type(1) = 0
      stencil_type(2) = 0
      stencil_type(3) = 0
      stencil_type(4) = 0

      shear_method = 0
      grid_type = 4
      area_type = 0

      if (grid_type == 1) then
         BUOY_FREQ_SQ(:,:,:,bid) = c0

         !-----------------------------------------------------------------------
         !
         !     compute buoyancy frequency and Richardson number at the bottom of 
         !     T box:
         !             Ri = -g*Dz(RHO)/RHO_0/(Dz(u)**2+Dz(v)**2)
         !
         !     RHO_0 ~ 1 in cgs units.
         !
         !-----------------------------------------------------------------------

         do k=1,km-1

            if ( k == 1 ) then
               TK = max(-c2, TMIX(:,:,k,1))
            endif

            kp1 = k+1

            call state (k, kp1, TMIX(:,:,k,1), TMIX(:,:,k,2), &
                        this_block, DRHODT=RHOT, DRHODS=RHOS)

            TKP = max(-c2, TMIX(:,:,kp1,1))

            WORK1 = TK - TKP
            WORK2 = TMIX(:,:,k,2) - TMIX(:,:,kp1,2)
            WORK3 = RHOT * WORK1 + RHOS * WORK2
            WORK3 = min(WORK3,-eps2)

            if (partial_bottom_cells) then
             where (k < KMT(:,:,bid))
              BUOY_FREQ_SQ(:,:,k,bid) = - grav * WORK3 / DZTW(:,:,k,bid)
             end where
            else
             where (k < KMT(:,:,bid))
              BUOY_FREQ_SQ(:,:,k,bid) = - grav * WORK3 * dzwr(k)
             end where
            endif

            TK  = TKP

         end do

         !-----------------------------------------------------------------------
         !
         !     compute the first baroclinic gravity-wave phase speed.
         !     Computation of Rossby deformation radius follows Chelton et al.(1998)
         !
         !-----------------------------------------------------------------------

         C_ROSSBY = c0

         k = 1
         if (partial_bottom_cells) then
          where ( k < KMT(:,:,bid) )
           C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k-1,bid)
          endwhere

          do k=1,km
           where ( k < KMT(:,:,bid) )
             C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * DZTW(:,:,k,bid)
           endwhere
           where ( k > 1  .and.  k == KMT(:,:,bid) )
             C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k-1,bid)) * DZTW(:,:,k,bid)
           endwhere
          enddo
         else
          where ( k < KMT(:,:,bid) )
           C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k-1)
          endwhere
   
          do k=1,km
           where ( k < KMT(:,:,bid) )
             C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k,bid)) * dzw(k)
           endwhere
           where ( k > 1  .and.  k == KMT(:,:,bid) )
             C_ROSSBY = C_ROSSBY + sqrt(BUOY_FREQ_SQ(:,:,k-1,bid)) * dzw(k)
           endwhere
          enddo
         endif

         C_ROSSBY = C_ROSSBY / pi

         L_ROSSBY = min( C_ROSSBY / (abs(FCORT(:,:,bid))+eps), &
                         sqrt( C_ROSSBY / (c2*BTP(:,:,bid)) ) )

         LDX(:,:,1) = L_ROSSBY !R1 - first Rossby radius of deformation
         do k=2,km
            LDX(:,:,k) = LDX(:,:,1) 
         enddo
         LDY = LDX

      elseif (grid_type == 2) then
         RNLOC = abs(VMIX) / (abs(UMIX)+eps)
         RDLOC = DXU(:,:,bid) / DYU(:,:,bid)
         do k=1,km
            where ( RNLOC(:,:,k) > RDLOC )
               LDX(:,:,k) = DXU(:,:,bid) * sqrt(c1 + c1/RNLOC(:,:,k)**2)
            elsewhere
               LDX(:,:,k) = DYU(:,:,bid) * sqrt(c1 +    RNLOC(:,:,k)**2)
            endwhere
         enddo
         LDY = LDX
      elseif (grid_type == 3) then
         LDX(:,:,1) = sqrt( DXU(:,:,bid)*DYU(:,:,bid) )
         do k=2,km
            LDX(:,:,k) = LDX(:,:,1) 
         enddo
         LDY = LDX
      elseif (grid_type == 4) then
         do k=1,km
            LDX(:,:,k) = DXU(:,:,bid)
            LDY(:,:,k) = DYU(:,:,bid)
         enddo

         RNLOC = abs(VMIX) / (abs(UMIX)+eps)
         RDLOC = DXU(:,:,bid) / DYU(:,:,bid)
         do k=1,km
            where ( RNLOC(:,:,k) > RDLOC )
               LDY(:,:,k) = DXU(:,:,bid) / RNLOC(:,:,k)
            elsewhere
               LDX(:,:,k) = DYU(:,:,bid) * RNLOC(:,:,k)
            endwhere
         enddo
      else
         do k=1,km
            LDX(:,:,k) = DXU(:,:,bid)
            LDY(:,:,k) = DYU(:,:,bid)
         enddo
      endif

      UFLX = UMIX * LDY
      VFLX = VMIX * LDX

      if (shear_method < 4) then
         call local_avg_at_t (UMIX, UAVG, stencil_type, area_type, this_block)
         call local_avg_at_t (VMIX, VAVG, stencil_type, area_type, this_block)
         NXLOC = UAVG
         NYLOC = VAVG
         WORK = sqrt(NXLOC**2 + NYLOC**2)
         where (WORK < eps2)
            NYLOC = c0
            NXLOC = c1
         elsewhere
            NYLOC = NYLOC / WORK
            NXLOC = NXLOC / WORK
         endwhere 
         if (shear_method == 0) then ! (0) MAJOR = kappa+1/kappa*<(u*dy)^2+(v*dx)^2>, flow aligned
            WORK = UFLX**2 + VFLX**2
            call local_avg_at_t (WORK, KMAJ, stencil_type, area_type, this_block)
         elseif (shear_method == 1) then ! (1) MAJOR = kappa+1/kappa*[<u*dy>^2+<v*dx>^2], flow aligned
            call local_avg_at_t (UFLX, WORKU, stencil_type, area_type, this_block)
            call local_avg_at_t (VFLX, WORKV, stencil_type, area_type, this_block)
            KMAJ = WORKU**2 + WORKV**2
         elseif (shear_method == 2) then ! (2) MAJOR = kappa+1/kappa*[<(u*dy)^2>+<(v*dx)^2>], flow aligned
            WORKU = UFLX**2 
            WORKV = VFLX**2 
            call local_avg_at_t (WORKU, KMAJ, stencil_type, area_type, this_block)
            call local_avg_at_t (WORKV, WORK, stencil_type, area_type, this_block)
            KMAJ = KMAJ + WORK
         elseif (shear_method == 3) then ! (3) MAJOR = kappa+1/kappa*<(|u*dy|+|v*dx|)^2>, flow aligned 
            WORK = ( abs(UFLX) + abs(VFLX) )**2
            call local_avg_at_t (WORK, KMAJ, stencil_type, area_type, this_block)
         endif

         NX_SHRD(:,:,:,bid) = NXLOC
         NY_SHRD(:,:,:,bid) = NYLOC
         do k_sub=1,2 
            KRAT_SHRD(:,:,k_sub,:,bid) = c1 + shrdispfac * KMAJ /             &
            ( abs( KAPPA_ISOP(:,:,k_sub,:,bid)+HOR_DIFF(:,:,k_sub,:,bid) ) + eps )**2
         enddo

         !KXX_ISOP(:,:,1,:,bid) = KMAJ * NXLOC**2
         !KXY_ISOP(:,:,1,:,bid) = KMAJ * NXLOC * NYLOC
         !KYY_ISOP(:,:,1,:,bid) = KMAJ * NYLOC**2
      elseif (shear_method == 4) then ! (4) K_tensor = kappa*[1 0; 0 1]+1/kappa*[<u*dy>^2 <u*dy><v*dx>; <u*dy><v*dx> <v*dx>^2]
         call local_avg_at_t (UFLX, WORKU, stencil_type, area_type, this_block)
         call local_avg_at_t (VFLX, WORKV, stencil_type, area_type, this_block)
         WORK  = WORKU*WORKV !KXY
         WORKV = WORKV**2 !KYY
         WORKU = WORKU**2 !KXX
      elseif (shear_method == 5) then ! (5) K_tensor = kappa*[1 0; 0 1]+1/kappa*[<(u*dy)^2> <u*dy*v*dx>; <u*dy*v*dx> <(v*dx)^2>]
         WORKU = UFLX*VFLX 
         call local_avg_at_t (WORKU, KMAJ, stencil_type, area_type, this_block)
         WORK = KMAJ !KXY
         
         WORKV = VFLX**2 
         call local_avg_at_t (WORKV, KMAJ, stencil_type, area_type, this_block)
         WORKV = KMAJ !KYY

         WORKU = UFLX**2 
         call local_avg_at_t (WORKU, KMAJ, stencil_type, area_type, this_block)
         WORKU = KMAJ !KXX
      endif

      if (shear_method > 3) then
         do k_sub=1,2
            KMAJ = abs( KAPPA_ISOP(:,:,k_sub,:,bid)+HOR_DIFF(:,:,k_sub,:,bid) ) !ISOTROPIC KAPPA
            WORKU = KMAJ + shrdispfac * WORKU / (KMAJ+eps) !KXX
            WORKV = KMAJ + shrdispfac * WORKV / (KMAJ+eps) !KYY
            WORK  =        shrdispfac * WORK  / (KMAJ+eps) !KXY

            WORKW = sqrt((WORKU+WORKV)**2-c4*(WORKU*WORKV-WORK**2))/c2 !EIGENVALUES are +/- this
            KMAJ  = (WORKU+WORKV)/c2+WORKW !MAJOR=(KXX+KYY)/2+sqrt((KXX+KYY)^2-4*(KXX*KYY-KXY^2))/2
            WORKW = (WORKU+WORKV)/c2-WORKW !MINOR=(KXX+KYY)/2-sqrt((KXX+KYY)^2-4*(KXX*KYY-KXY^2))/2
            WORKV = WORK**2 + (KMAJ - WORKU)**2 !NORMALIZATION FACTOR=KXY^2+(MAJOR-KXX)^2
            where (WORKV < eps2)
               NY_SHRD(:,:,:,bid) = c0
               NX_SHRD(:,:,:,bid) = c1
            elsewhere
               NY_SHRD(:,:,:,bid) = (KMAJ - WORKU)/sqrt(WORKV) !\vec{n}=(KXY,MAJOR-KXX)... normalized here
               NX_SHRD(:,:,:,bid) = WORK/sqrt(WORKV) 
            endwhere 

            KRAT_SHRD(:,:,k_sub,:,bid) = KMAJ / WORKW !RAT=MAJOR/MINOR
         enddo
      endif

!-----------------------------------------------------------------------
!EOC

      end subroutine shear_dispersion 


!***********************************************************************
!BOP
! !IROUTINE: local_avg_at_t 
! !INTERFACE:

      subroutine local_avg_at_t (FGEN, FAVG, stencil_type, area_type, this_block)

! !DESCRIPTION: local average of U-type variable neighbors at T-cell centers
!
! !REVISION HISTORY:

! !INPUT PARAMETERS:

      real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
         FGEN               ! general variable defined at U-pts

      integer (int_kind), intent(in) :: &
         area_type  ! type of area handling - 0=just do ugrid_to_tgrid, 1=UAREA, 2=SUBCELLV

      integer (int_kind), dimension(4), intent(in) :: &
         stencil_type      ! type of stencil added to ME for averaging - 
                           ! index 1=+E+W, 2=+N+S, 3=+NE+SW, 4=+NW+SE

      type (block), intent(in) :: &
         this_block         ! block info for this sub block

      real (r8), dimension(nx_block,ny_block,km), intent(out) :: &
         FAVG               ! general variable local average

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      integer (int_kind) :: &
         k,      &
         bid                ! local block address for this sub block

      real (r8), dimension(nx_block,ny_block,km) :: &
         SUMF, TOTVOL, TUMF, TTTVOL

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------

     bid = this_block%local_id
     if (area_type == 0) then

      do k=1,km
         call ugrid_to_tgrid(FAVG(:,:,k),FGEN(:,:,k),bid)
      enddo

     else
      
      SUMF = c0
      TOTVOL = c0
      TUMF = c0
      TTTVOL = c0
      FAVG = c0

      if (area_type == 1) then
         if (partial_bottom_cells) then
            do k=1,km
               TUMF(2:nx_block,2:ny_block,k) = &
                                             FGEN(1:nx_block-1,2:ny_block  ,k) * &
                                            UAREA(1:nx_block-1,2:ny_block  ,bid) * &
                                              DZT(1:nx_block-1,2:ny_block  ,k,bid) +  &
                                             FGEN(2:nx_block  ,2:ny_block  ,k) * &
                                            UAREA(2:nx_block  ,2:ny_block  ,bid) * &
                                              DZT(2:nx_block  ,2:ny_block  ,k,bid) +  &
                                             FGEN(1:nx_block-1,1:ny_block-1,k) * &
                                            UAREA(1:nx_block-1,1:ny_block-1,bid) * &
                                              DZT(1:nx_block-1,1:ny_block-1,k,bid) +  &
                                             FGEN(2:nx_block  ,1:ny_block-1,k) * &
                                            UAREA(2:nx_block  ,1:ny_block-1,bid) * &
                                              DZT(2:nx_block  ,1:ny_block-1,k,bid) 
               TTTVOL(2:nx_block,2:ny_block,k) = &
                                            UAREA(1:nx_block-1,2:ny_block  ,bid) * &
                                              DZT(1:nx_block-1,2:ny_block  ,k,bid) +  &
                                            UAREA(2:nx_block  ,2:ny_block  ,bid) * &
                                              DZT(2:nx_block  ,2:ny_block  ,k,bid) +  &
                                            UAREA(1:nx_block-1,1:ny_block-1,bid) * &
                                              DZT(1:nx_block-1,1:ny_block-1,k,bid) +  &
                                            UAREA(2:nx_block  ,1:ny_block-1,bid) * &
                                              DZT(2:nx_block  ,1:ny_block-1,k,bid) 
            enddo
         else
            do k=1,km
               TUMF(2:nx_block,2:ny_block,k) = dz(k) * ( &
                                             FGEN(1:nx_block-1,2:ny_block  ,k) * &
                                            UAREA(1:nx_block-1,2:ny_block  ,bid) +  &
                                             FGEN(2:nx_block  ,2:ny_block  ,k) * &
                                            UAREA(2:nx_block  ,2:ny_block  ,bid) +  &
                                             FGEN(1:nx_block-1,1:ny_block-1,k) * &
                                            UAREA(1:nx_block-1,1:ny_block-1,bid) +  &
                                             FGEN(2:nx_block  ,1:ny_block-1,k) * &
                                            UAREA(2:nx_block  ,1:ny_block-1,bid) )
               TTTVOL(2:nx_block,2:ny_block,k) = dz(k) * ( &
                                            UAREA(1:nx_block-1,2:ny_block  ,bid) +  &
                                            UAREA(2:nx_block  ,2:ny_block  ,bid) +  &
                                            UAREA(1:nx_block-1,1:ny_block-1,bid) +  &
                                            UAREA(2:nx_block  ,1:ny_block-1,bid) )
            enddo
         endif
      else
         do k=1,km
            TUMF(2:nx_block,2:ny_block,k) = &
                                             FGEN(1:nx_block-1,2:ny_block  ,k) * &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,2,1,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,2,1,2,k,bid) ) + &
                                             FGEN(2:nx_block  ,2:ny_block  ,k) * &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,1,1,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,1,1,2,k,bid) ) + &
                                             FGEN(1:nx_block-1,1:ny_block-1,k) * &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,2,2,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,2,2,2,k,bid) ) + &
                                             FGEN(2:nx_block  ,1:ny_block-1,k) * &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,1,2,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,1,2,2,k,bid) ) 
            TTTVOL(2:nx_block,2:ny_block,k) = &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,2,1,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,2,1,2,k,bid) ) + &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,1,1,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,1,1,2,k,bid) ) + &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,2,2,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,2,2,2,k,bid) ) + &
                                     (   SUBCELLV(2:nx_block  ,2:ny_block  ,1,2,1,k,bid) + &
                                         SUBCELLV(2:nx_block  ,2:ny_block  ,1,2,2,k,bid) ) 
         enddo
      endif

      if (stencil_type(1) == 1) then
            SUMF(2:nx_block-1,:,:) = SUMF(2:nx_block-1,:,:) + & !E+W
                                     TUMF(1:nx_block-2,:,:) + TUMF(3:nx_block,:,:)
            TOTVOL(2:nx_block-1,:,:) = TOTVOL(2:nx_block-1,:,:) + & !E+W
                                       TTTVOL(1:nx_block-2,:,:) + TTTVOL(3:nx_block,:,:)
      endif

      if (stencil_type(2) == 1) then
            SUMF(:,2:ny_block-1,:) = SUMF(:,2:ny_block-1,:) + & !N+S
                                     TUMF(:,1:ny_block-2,:) + TUMF(:,3:ny_block,:)
            TOTVOL(:,2:ny_block-1,:) = TOTVOL(:,2:ny_block-1,:) + & !N+S
                                       TTTVOL(:,1:ny_block-2,:) + TTTVOL(:,3:ny_block,:)
      endif

      if (stencil_type(3) == 1) then
            SUMF(2:nx_block-1,2:ny_block-1,:) = SUMF(2:nx_block-1,2:ny_block-1,:) + & !NE+SW
                                                TUMF(1:nx_block-2,1:ny_block-2,:) + TUMF(3:nx_block,3:ny_block,:)
            TOTVOL(2:nx_block-1,2:ny_block-1,:) = TOTVOL(2:nx_block-1,2:ny_block-1,:) + & !NE+SW
                                                  TTTVOL(1:nx_block-2,1:ny_block-2,:) + TTTVOL(3:nx_block,3:ny_block,:)
      endif

      if (stencil_type(4) == 1) then
            SUMF(2:nx_block-1,2:ny_block-1,:) = SUMF(2:nx_block-1,2:ny_block-1,:) + & !NW+SE
                                                TUMF(3:nx_block  ,1:ny_block-2,:) + TUMF(1:nx_block-2,3:ny_block,:)
            TOTVOL(2:nx_block-1,2:ny_block-1,:) = TOTVOL(2:nx_block-1,2:ny_block-1,:) + & !NE+SW
                                                  TTTVOL(3:nx_block  ,1:ny_block-2,:) + TTTVOL(1:nx_block-2,3:ny_block,:)
      endif

      SUMF = SUMF + TUMF
      TOTVOL = TOTVOL + TTTVOL

      FAVG(2:nx_block-1,2:ny_block-1,:) = SUMF(2:nx_block-1,2:ny_block-1,:) / ( TOTVOL(2:nx_block-1,2:ny_block-1,:) + eps )
     endif
!-----------------------------------------------------------------------
!EOC

      end subroutine local_avg_at_t

!***********************************************************************

      end module hmix_gm_aniso

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
