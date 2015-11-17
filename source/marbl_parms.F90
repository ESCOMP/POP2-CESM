module marbl_parms

  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module ecosys_mod.
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist marbl_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !   CVS:$Id: marbl_parms.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------
  !   Modified to include parameters for diazotrophs, JKM  4/2002
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this ecosys, so are used at
  !   the module level. The use statements for variables that are only needed
  !   locally are located at the module subprogram level.
  !-----------------------------------------------------------------------------

  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : log_kind

  use marbl_kinds_mod, only : c1

  use marbl_share_mod, only : zooplankton
  use marbl_share_mod, only : autotrophs
  use marbl_share_mod, only : grazing
  use marbl_share_mod, only : sp_ind
  use marbl_share_mod, only : diaz_ind
  use marbl_share_mod, only : diat_ind
  use marbl_share_mod, only : zooplankton_cnt
  use marbl_share_mod, only : autotroph_cnt
  use marbl_share_mod, only : grazer_prey_cnt

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   all module variables are public and should have their values preserved
  !-----------------------------------------------------------------------------

  PUBLIC
  SAVE

  !-----------------------------------------------------------------------
  !  MARBL indices for surface fluxes
  !-----------------------------------------------------------------------

  integer (int_kind), parameter :: &
       ind_nox_flux     =  1,  &
       ind_nhy_flux     =  2,  &
       ind_no3_flux     =  3,  &
       ind_nh4_flux     =  4,  &
       ind_din_riv_flux =  5,  &
       ind_dip_riv_flux =  6,  &
       ind_don_riv_flux =  7,  &
       ind_dop_riv_flux =  8,  &
       ind_dsi_riv_flux =  9,  &
       ind_dfe_riv_flux = 10,  &
       ind_dic_riv_flux = 11,  &
       ind_alk_riv_flux = 12,  &
       ind_doc_riv_flux = 13

  !-----------------------------------------------------------------------
  !  non-autotroph relative tracer indices
  !  autotroph relative tracer indices are in autotroph derived type and are determined at run time
  !-----------------------------------------------------------------------

  integer (int_kind), parameter :: &
       po4_ind         =  1,  & ! dissolved inorganic phosphate
       no3_ind         =  2,  & ! dissolved inorganic nitrate
       sio3_ind        =  3,  & ! dissolved inorganic silicate
       nh4_ind         =  4,  & ! dissolved ammonia
       fe_ind          =  5,  & ! dissolved inorganic iron
       o2_ind          =  6,  & ! dissolved oxygen
       dic_ind         =  7,  & ! dissolved inorganic carbon
       dic_alt_co2_ind =  8,  & ! dissolved inorganic carbon with alternative CO2
       alk_ind         =  9,  & ! alkalinity
       doc_ind         = 10,  & ! dissolved organic carbon
       don_ind         = 11,  & ! dissolved organic nitrogen
       dofe_ind        = 12,  & ! dissolved organic iron
       dop_ind         = 13,  & ! dissolved organic phosphorus
       dopr_ind        = 14,  & ! refractory DOP
       donr_ind        = 15     ! refractory DON

  !-----------------------------------------------------------------------------
  !   epsilon values
  !-----------------------------------------------------------------------------

   real(kind=r8), parameter :: &
      epsC      = 1.00e-8, & ! small C concentration (mmol C/m^3)
      epsTinv   = 3.17e-8    ! small inverse time scale (1/year) (1/sec)

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       spd = 86400.0_r8,        & ! number of seconds in a day
       dps = c1 / spd,          & ! number of days in a second
       yps = c1 / (365.0_r8*spd)  ! number of years in a second


  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       parm_Red_D_C_P  = 117.0_r8,                 & ! carbon:phosphorus
       parm_Red_D_N_P  =  16.0_r8,                 & ! nitrogen:phosphorus
       parm_Red_D_O2_P = 170.0_r8,                 & ! oxygen:phosphorus
       parm_Remin_D_O2_P = 138.0_r8,               & ! oxygen:phosphorus
       parm_Red_P_C_P  = parm_Red_D_C_P,                 & ! carbon:phosphorus
       parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,  & ! carbon:nitrogen
       parm_Red_P_C_N  = parm_Red_D_C_N,                 & ! carbon:nitrogen
       parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P, & ! carbon:oxygen
       parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P, & ! carbon:oxygen
       parm_Red_P_C_O2 = parm_Red_D_C_O2,                & ! carbon:oxygen
       parm_Red_Fe_C   = 3.0e-6_r8,                & ! iron:carbon
       parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_r8! carbon:oxygen
                                                           ! for diazotrophs

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

  REAL(KIND=r8) :: &
       parm_Fe_bioavail,      & ! fraction of Fe flux that is bioavailable
       parm_o2_min,           & ! min O2 needed for prod & consump. (nmol/cm^3)
       parm_o2_min_delta,     & ! width of min O2 range (nmol/cm^3)
       parm_kappa_nitrif,     & ! nitrification inverse time constant (1/sec)
       parm_nitrif_par_lim,   & ! PAR limit for nitrif. (W/m^2)
       parm_labile_ratio,     & ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
       parm_POMbury,          & ! scale factor for burial of POC, PON, and POP
       parm_BSIbury,          & ! scale factor burial of bSi
       parm_fe_scavenge_rate0,& ! base scavenging rate
       parm_f_prod_sp_CaCO3,  & !fraction of sp prod. as CaCO3 prod.
       parm_POC_diss,         & ! base POC diss len scale
       parm_SiO2_diss,        & ! base SiO2 diss len scale
       parm_CaCO3_diss          ! base CaCO3 diss len scale

  REAL(KIND=r8), DIMENSION(4) :: &
       parm_scalelen_z,       & ! depths of prescribed scalelen values
       parm_scalelen_vals       ! prescribed scalelen values

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       fe_scavenge_thres1 = 0.8e-3_r8,  & !upper thres. for Fe scavenging
       dust_fescav_scale  = 1.0e9,      & !dust scavenging scale factor
       fe_max_scale2      = 1200.0_r8     !unitless scaling coeff.

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out.
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !			    gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (nmol Fe/g Dust) 
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       dust_to_Fe=0.035_r8/55.847_r8*1.0e9_r8
 
  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      caco3_poc_min    = 0.4_r8,  & !minimum proportionality between 
                                          !   QCaCO3 and grazing losses to POC 
                                          !   (mmol C/mmol CaCO3)
      spc_poc_fac      = 0.11_r8, & !small phyto grazing factor (1/mmolC)
      f_graze_sp_poc_lim = 0.3_r8, & 
      f_photosp_CaCO3  = 0.4_r8,  & !proportionality between small phyto 
                                          !    production and CaCO3 production
      f_graze_CaCO3_remin = 0.33_r8, & !fraction of spCaCO3 grazing 
                                             !          which is remin
      f_graze_si_remin    = 0.35_r8      !fraction of diatom Si grazing 
                                             !          which is remin

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       r_Nfix_photo=1.25_r8         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      Q             = 0.137_r8,  & !N/C ratio (mmol/mmol) of phyto & zoo
      Qp_zoo_pom    = 0.00855_r8,& !P/C ratio (mmol/mmol) zoo & pom
      Qfe_zoo       = 3.0e-6_r8, & !zooplankton fe/C ratio
      gQsi_0        = 0.137_r8,  & !initial Si/C ratio
      gQsi_max      = 0.685_r8,  & !max Si/C ratio
      gQsi_min      = 0.0457_r8, & !min Si/C ratio
      QCaCO3_max    = 0.4_r8,    & !max QCaCO3
      ! carbon:nitrogen ratio for denitrification
      denitrif_C_N  = parm_Red_D_C_P/136.0_r8

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
      thres_z1          = 100.0e2_r8, & !threshold = C_loss_thres for z shallower than this (cm)
      thres_z2          = 150.0e2_r8, & !threshold = 0 for z deeper than this (cm)
      CaCO3_temp_thres1 = 6.0_r8,   & !upper temp threshold for CaCO3 prod
      CaCO3_temp_thres2 = -2.0_r8,  & !lower temp threshold
      CaCO3_sp_thres    = 4.0_r8      ! bloom condition thres (mmolC/m3)

  !---------------------------------------------------------------------
  !     grazing functions
  !---------------------------------------------------------------------

  INTEGER (INT_KIND), PARAMETER ::   &
         grz_fnc_michaelis_menten = 1,       &
         grz_fnc_sigmoidal        = 2

  !---------------------------------------------------------------------
  !     fraction of incoming shortwave assumed to be PAR
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       f_qsw_par = 0.45_r8   ! PAR fraction

        
  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       Tref = 30.0_r8, & ! reference temperature (C)
       Q_10 = 1.5_r8     ! factor for temperature dependence (non-dim)

  !---------------------------------------------------------------------
  !  DOM parameters for refractory components and DOP uptake
  !---------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: &
       DOC_reminR  = (c1/250.0_r8) * dps,          & ! rate for semi-labile DOC 1/250days
       DON_reminR  = (c1/160.0_r8) * dps,          & ! rate for semi-labile DON 1/160days
       DOFe_reminR = (c1/160.0_r8) * dps,          & ! rate for semi-labile DOFe 1/160days
       DOP_reminR  = (c1/160.0_r8) * dps,          & ! rate for semi-labile DOP 1/160days  
       DONr_reminR = (c1/(365.0_r8*2.5_r8)) * dps, & ! timescale for refrac DON 1/2.5yrs
       DOPr_reminR = (c1/(365.0_r8*2.5_r8)) * dps, & ! timescale for refrac DOP 1/2.5yrs
       DONrefract = 0.08_r8,                       & ! fraction of DON to refractory pool
       DOPrefract = 0.03_r8                          ! fraction of DOP to refractory pool

  !*****************************************************************************

  public :: &
       marbl_params_init, &
       marbl_params_print

  private :: &
       marbl_params_set_defaults
  
contains

  !*****************************************************************************

  subroutine marbl_params_set_defaults()
    ! assign default parameter values

    implicit none

    integer :: auto_ind, zoo_ind, prey_ind
    
    parm_Fe_bioavail       = 1.0_r8
    parm_o2_min            = 4.0_r8
    parm_o2_min_delta      = 2.0_r8
    parm_kappa_nitrif      = 0.06_r8 * dps  ! (= 1/( days))
    parm_nitrif_par_lim    = 1.0_r8
    parm_labile_ratio      = 0.85_r8
    parm_POMbury           = 1.4_r8         ! x1 default
    parm_BSIbury           = 0.65_r8        ! x1 default
    parm_fe_scavenge_rate0 = 3.0_r8         ! x1 default
    parm_f_prod_sp_CaCO3   = 0.055_r8       ! x1 default
    parm_POC_diss          = 88.0e2_r8
    parm_SiO2_diss         = 250.0e2_r8
    parm_CaCO3_diss        = 150.0e2_r8

    parm_scalelen_z    = (/ 130.0e2_r8, 290.0e2_r8, 670.0e2_r8, 1700.0e2_r8 /) ! x1 default
    parm_scalelen_vals = (/     1.0_r8,     3.0_r8,     5.0_r8,      9.0_r8 /) ! x1 default

    zoo_ind = 1
    zooplankton(zoo_ind)%sname          ='zoo'
    zooplankton(zoo_ind)%lname          = 'Zooplankton'
    zooplankton(zoo_ind)%z_mort_0       = 0.1_r8 * dps
    zooplankton(zoo_ind)%z_mort2_0      = 0.4_r8 * dps
    zooplankton(zoo_ind)%loss_thres     = 0.005_r8     !zoo conc. where losses go to zero

    auto_ind = sp_ind
    autotrophs(auto_ind)%sname         = 'sp'
    autotrophs(auto_ind)%lname         = 'Small Phyto'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .true.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%kFe           = 0.04e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.01_r8
    autotrophs(auto_ind)%kDOP          = 0.26_r8
    autotrophs(auto_ind)%kNO3          = 0.1_r8
    autotrophs(auto_ind)%kNH4          = 0.01_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_r8
    autotrophs(auto_ind)%Qp            = 0.00855_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.6_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_r8
    autotrophs(auto_ind)%temp_thres    = -10.0_r8
    autotrophs(auto_ind)%mort          = 0.12_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.01_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8

    auto_ind = diat_ind
    autotrophs(auto_ind)%sname         = 'diat'
    autotrophs(auto_ind)%lname         = 'Diatom'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%kFe           = 0.06e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.05_r8
    autotrophs(auto_ind)%kDOP          = 0.9_r8
    autotrophs(auto_ind)%kNO3          = 0.5_r8
    autotrophs(auto_ind)%kNH4          = 0.05_r8
    autotrophs(auto_ind)%kSiO3         = 0.8_r8
    autotrophs(auto_ind)%Qp            = 0.00855_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.465_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 4.0_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_r8
    autotrophs(auto_ind)%temp_thres    = -10.0_r8
    autotrophs(auto_ind)%mort          = 0.12_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.02_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8
   

    auto_ind = diaz_ind
    autotrophs(auto_ind)%sname         = 'diaz'
    autotrophs(auto_ind)%lname         = 'Diazotroph'
    autotrophs(auto_ind)%Nfixer        = .true.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%kFe           = 0.04e-3_r8
    autotrophs(auto_ind)%kPO4          = 0.02_r8
    autotrophs(auto_ind)%kDOP          = 0.09_r8
    autotrophs(auto_ind)%kNO3          = 1.0_r8
    autotrophs(auto_ind)%kNH4          = 0.15_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_r8
    autotrophs(auto_ind)%Qp            = 0.002735_r8
    autotrophs(auto_ind)%gQfe_0        = 60.0e-6_r8
    autotrophs(auto_ind)%gQfe_min      = 12.0e-6_r8
    autotrophs(auto_ind)%alphaPI       = 0.4_r8 * dps
    autotrophs(auto_ind)%PCref         = 0.7_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_r8
    autotrophs(auto_ind)%loss_thres    = 0.022_r8
    autotrophs(auto_ind)%loss_thres2   = 0.001_r8
    autotrophs(auto_ind)%temp_thres    = 14.0_r8
    autotrophs(auto_ind)%mort          = 0.15_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.0_r8
    autotrophs(auto_ind)%agg_rate_max  = 0.0_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.0_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_r8
 

    !---------------------------------------------------------------------------
    ! predator-prey relationships
    !---------------------------------------------------------------------------
    zoo_ind = 1
    prey_ind = sp_ind
    grazing(prey_ind,zoo_ind)%sname            = 'grz_' // autotrophs(prey_ind)%sname // '_' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%lname            = 'Grazing of ' // autotrophs(prey_ind)%sname // ' by ' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%auto_ind(1)      = prey_ind
    grazing(prey_ind,zoo_ind)%auto_ind_cnt       = 1
    grazing(prey_ind,zoo_ind)%zoo_ind          = -1
    grazing(prey_ind,zoo_ind)%zoo_ind_cnt        = 0
    grazing(prey_ind,zoo_ind)%z_umax_0         = 3.3_r8 * dps ! x1 default
    grazing(prey_ind,zoo_ind)%z_grz            = 1.05_r8              
    grazing(prey_ind,zoo_ind)%graze_zoo        = 0.3_r8
    grazing(prey_ind,zoo_ind)%graze_poc        = 0.0_r8
    grazing(prey_ind,zoo_ind)%graze_doc        = 0.15_r8
    grazing(prey_ind,zoo_ind)%f_zoo_detr       = 0.15_r8
    grazing(prey_ind,zoo_ind)%grazing_function = grz_fnc_michaelis_menten

    prey_ind = diat_ind
    grazing(prey_ind,zoo_ind)%sname            = 'grz_' // autotrophs(prey_ind)%sname // '_' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%lname            = 'Grazing of ' // autotrophs(prey_ind)%sname // ' by ' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%auto_ind(1)      = prey_ind
    grazing(prey_ind,zoo_ind)%auto_ind_cnt       = 1
    grazing(prey_ind,zoo_ind)%zoo_ind          = -1
    grazing(prey_ind,zoo_ind)%zoo_ind_cnt        = 0
    grazing(prey_ind,zoo_ind)%z_umax_0         = 3.08_r8 * dps ! x1 default
    grazing(prey_ind,zoo_ind)%z_grz            = 1.0_r8              
    grazing(prey_ind,zoo_ind)%graze_zoo        = 0.3_r8
    grazing(prey_ind,zoo_ind)%graze_poc        = 0.42_r8
    grazing(prey_ind,zoo_ind)%graze_doc        = 0.15_r8
    grazing(prey_ind,zoo_ind)%f_zoo_detr       = 0.2_r8
    grazing(prey_ind,zoo_ind)%grazing_function = grz_fnc_michaelis_menten

    prey_ind = diaz_ind
    grazing(prey_ind,zoo_ind)%sname            = 'grz_' // autotrophs(prey_ind)%sname // '_' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%lname            = 'Grazing of ' // autotrophs(prey_ind)%sname // ' by ' // zooplankton(zoo_ind)%sname
    grazing(prey_ind,zoo_ind)%auto_ind(1)      = prey_ind
    grazing(prey_ind,zoo_ind)%auto_ind_cnt       = 1
    grazing(prey_ind,zoo_ind)%zoo_ind          = -1
    grazing(prey_ind,zoo_ind)%zoo_ind_cnt        = 0
    grazing(prey_ind,zoo_ind)%z_umax_0         = 0.6_r8 * dps
    grazing(prey_ind,zoo_ind)%z_grz            = 1.2_r8              
    grazing(prey_ind,zoo_ind)%graze_zoo        = 0.3_r8
    grazing(prey_ind,zoo_ind)%graze_poc        = 0.05_r8
    grazing(prey_ind,zoo_ind)%graze_doc        = 0.15_r8
    grazing(prey_ind,zoo_ind)%f_zoo_detr       = 0.15_r8
    grazing(prey_ind,zoo_ind)%grazing_function = grz_fnc_michaelis_menten

    
  end subroutine marbl_params_set_defaults

  !*****************************************************************************

  subroutine marbl_params_init(nl_buffer, marbl_status)

    use marbl_interface_constants, only: marbl_nl_buffer_size
    use marbl_interface_constants, only: marbl_status_ok, marbl_status_could_not_read_namelist
    use marbl_interface_types, only: marbl_status_type

    implicit none

    character(marbl_nl_buffer_size), intent(in) :: nl_buffer
    type(marbl_status_type), intent(inout) :: marbl_status

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: subname = 'marbl_parms:marbl_parms_init'

    integer(kind=int_kind) :: io_error

    NAMELIST /ecosys_parms_nml/ &
         parm_Fe_bioavail, &
         parm_o2_min, &
         parm_o2_min_delta, &
         parm_kappa_nitrif, &
         parm_nitrif_par_lim, &
         parm_labile_ratio, &
         parm_POMbury, &
         parm_BSIbury, &
         parm_fe_scavenge_rate0, &
         parm_f_prod_sp_CaCO3, &
         parm_POC_diss, &
         parm_SiO2_diss, &
         parm_CaCO3_diss, &
         parm_scalelen_z, &
         parm_scalelen_vals, &
         autotrophs, & 
         zooplankton, &
         grazing

    marbl_status%status = marbl_status_ok
    marbl_status%message = ''

    call marbl_params_set_defaults()

    !---------------------------------------------------------------------------
    ! read in namelist to override some defaults
    !---------------------------------------------------------------------------
    read(nl_buffer, nml=ecosys_parms_nml, iostat=io_error)
    if (io_error /= 0) then
       marbl_status%status = marbl_status_could_not_read_namelist
       marbl_status%message = "ERROR: marbl_parmams_read_namelist(): could not read namelist 'ecosys_parms_nml'."
       return
    end if

  end subroutine marbl_params_init

  !*****************************************************************************

  subroutine marbl_params_print(stdout)
    ! echo all parameters to the specified output file

    implicit none

    integer, intent(in) :: stdout

    integer :: zoo_ind, auto_ind, prey_ind
    
    !---------------------------------------------------------------------------

    write(stdout, *) '----------------------------------------'
    write(stdout, *) '----- marbl_parms namelist values -----'
    write(stdout, *) 'parm_Fe_bioavail       = ', parm_Fe_bioavail
    write(stdout, *) 'parm_o2_min            = ', parm_o2_min
    write(stdout, *) 'parm_o2_min_delta      = ', parm_o2_min_delta
    write(stdout, *) 'parm_kappa_nitrif      = ', parm_kappa_nitrif
    write(stdout, *) 'parm_nitrif_par_lim    = ', parm_nitrif_par_lim
    write(stdout, *) 'parm_labile_ratio      = ', parm_labile_ratio
    write(stdout, *) 'parm_POMbury           = ', parm_POMbury
    write(stdout, *) 'parm_BSIbury           = ', parm_BSIbury
    write(stdout, *) 'parm_fe_scavenge_rate0 = ', parm_fe_scavenge_rate0
    write(stdout, *) 'parm_f_prod_sp_CaCO3   = ', parm_f_prod_sp_CaCO3
    write(stdout, *) 'parm_POC_diss          = ', parm_POC_diss
    write(stdout, *) 'parm_SiO2_diss         = ', parm_SiO2_diss
    write(stdout, *) 'parm_CaCO3_diss        = ', parm_CaCO3_diss
    write(stdout, *) 'parm_scalelen_z        = ', parm_scalelen_z
    write(stdout, *) 'parm_scalelen_vals     = ', parm_scalelen_vals

    do zoo_ind = 1, zooplankton_cnt
       write(stdout, *) 'lname(', trim(zooplankton(zoo_ind)%sname), ') = ', zooplankton(zoo_ind)%lname
       write(stdout, *) 'z_mort_0(', trim(zooplankton(zoo_ind)%sname), ') = ', zooplankton(zoo_ind)%z_mort_0
       write(stdout, *) 'z_mort2_0(', trim(zooplankton(zoo_ind)%sname), ') = ', zooplankton(zoo_ind)%z_mort2_0
       write(stdout, *) 'loss_thres(', trim(zooplankton(zoo_ind)%sname), ') = ', zooplankton(zoo_ind)%loss_thres
    end do
    
    do auto_ind = 1, autotroph_cnt
       write(stdout, *) 'lname(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%lname
       write(stdout, *) 'Nfixer(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Nfixer
       write(stdout, *) 'imp_calcifier(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%imp_calcifier
       write(stdout, *) 'exp_calcifier(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%exp_calcifier
       write(stdout, *) 'kFe(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kFe
       write(stdout, *) 'kPO4(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kPO4
       write(stdout, *) 'kdoP(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kdoP
       write(stdout, *) 'kNO3(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kNO3
       write(stdout, *) 'kNH4(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kNH4
       write(stdout, *) 'kSiO3(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%kSiO3
       write(stdout, *) 'Qp(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Qp
       write(stdout, *) 'gQfe_0(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%gQfe_0
       write(stdout, *) 'gQfe_min(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%gQfe_min
       write(stdout, *) 'alphaPI(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%alphaPI
       write(stdout, *) 'PCref(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%PCref
       write(stdout, *) 'thetaN_max(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%thetaN_max
       write(stdout, *) 'loss_thres(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_thres
       write(stdout, *) 'loss_thres2(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_thres2
       write(stdout, *) 'temp_thres(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%temp_thres
       write(stdout, *) 'mort(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%mort
       write(stdout, *) 'mort2(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%mort2
       write(stdout, *) 'agg_rate_max(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%agg_rate_max
       write(stdout, *) 'agg_rate_min(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%agg_rate_min
       write(stdout, *) 'loss_poc(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%loss_poc
    end do

    
    do prey_ind = 1, grazer_prey_cnt
       do zoo_ind = 1, zooplankton_cnt
          write(stdout, *) 'lname(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%lname
          write(stdout, *) 'auto_ind(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%auto_ind
          write(stdout, *) 'auto_ind_cnt(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%auto_ind_cnt
          write(stdout, *) 'zoo_ind(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%zoo_ind
          write(stdout, *) 'zoo_ind_cnt(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%zoo_ind_cnt
          write(stdout, *) 'z_umax_0(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%z_umax_0
          write(stdout, *) 'z_grz(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%z_grz
          write(stdout, *) 'graze_zoo(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%graze_zoo
          write(stdout, *) 'graze_poc(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%graze_poc
          write(stdout, *) 'graze_doc(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%graze_doc
          write(stdout, *) 'f_zoo_detr(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%f_zoo_detr
          write(stdout, *) 'grazing_function(', trim(grazing(prey_ind,zoo_ind)%sname), ') = ', grazing(prey_ind,zoo_ind)%grazing_function
       end do
    end do

    write(stdout, *) '----------------------------------------'

  end subroutine marbl_params_print

  !*****************************************************************************

end module marbl_parms
