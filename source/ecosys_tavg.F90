!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ecosys_tavg

!BOP
! !MODULE: ecosys_tavg

! !DESCRIPTION:
!  This module provides support for writing MARBL diagnostics into POP tavg
!  files
!
!  Written by: Michael Levy, NCAR, Dec 2014


! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

  use kinds_mod, only: r8
  use kinds_mod, only: int_kind
  use kinds_mod, only: char_len

  use blocks, only: nx_block
  use blocks, only: ny_block

  use domain, only : nblocks_clinic
  use domain_size, only : km

  use tavg, only : define_tavg_field
  use tavg, only : accumulate_tavg_field

  use ecosys_constants, only : ecosys_tracer_cnt
  ! standard diagnostics
  ! 2d
  use ecosys_diagnostics_mod, only : ecosys_diag_cnt_2d
  use ecosys_diagnostics_mod, only : zsatcalc_diag_ind
  use ecosys_diagnostics_mod, only : zsatarag_diag_ind
  use ecosys_diagnostics_mod, only : O2_ZMIN_diag_ind
  use ecosys_diagnostics_mod, only : O2_ZMIN_DEPTH_diag_ind
  use ecosys_diagnostics_mod, only : photoC_TOT_zint_diag_ind
  use ecosys_diagnostics_mod, only : photoC_NO3_TOT_zint_diag_ind
  use ecosys_diagnostics_mod, only : Jint_Ctot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_100m_Ctot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_Ntot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_100m_Ntot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_Ptot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_100m_Ptot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_Sitot_diag_ind
  use ecosys_diagnostics_mod, only : Jint_100m_Sitot_diag_ind
  ! 3d
  use ecosys_diagnostics_mod, only : ecosys_diag_cnt_3d
  use ecosys_diagnostics_mod, only : CO3_diag_ind
  use ecosys_diagnostics_mod, only : HCO3_diag_ind
  use ecosys_diagnostics_mod, only : H2CO3_diag_ind
  use ecosys_diagnostics_mod, only : pH_3D_diag_ind
  use ecosys_diagnostics_mod, only : CO3_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : HCO3_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : H2CO3_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : ph_3D_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : co3_sat_calc_diag_ind
  use ecosys_diagnostics_mod, only : co3_sat_arag_diag_ind
  use ecosys_diagnostics_mod, only : NITRIF_diag_ind
  use ecosys_diagnostics_mod, only : DENITRIF_diag_ind
  use ecosys_diagnostics_mod, only : O2_PRODUCTION_diag_ind
  use ecosys_diagnostics_mod, only : O2_CONSUMPTION_diag_ind
  use ecosys_diagnostics_mod, only : AOU_diag_ind
  use ecosys_diagnostics_mod, only : PAR_avg_diag_ind
  use ecosys_diagnostics_mod, only : auto_graze_TOT_diag_ind
  use ecosys_diagnostics_mod, only : photoC_TOT_diag_ind
  use ecosys_diagnostics_mod, only : photoC_NO3_TOT_diag_ind
  use ecosys_diagnostics_mod, only : DOC_prod_diag_ind
  use ecosys_diagnostics_mod, only : DOC_remin_diag_ind
  use ecosys_diagnostics_mod, only : DON_prod_diag_ind
  use ecosys_diagnostics_mod, only : DON_remin_diag_ind
  use ecosys_diagnostics_mod, only : DOP_prod_diag_ind
  use ecosys_diagnostics_mod, only : DOP_remin_diag_ind
  use ecosys_diagnostics_mod, only : DOFe_prod_diag_ind
  use ecosys_diagnostics_mod, only : DOFe_remin_diag_ind
  use ecosys_diagnostics_mod, only : Fe_scavenge_diag_ind
  use ecosys_diagnostics_mod, only : Fe_scavenge_rate_diag_ind

  ! autotroph diagnostics
  ! 2d
  use ecosys_diagnostics_mod, only : auto_diag_cnt_2d
  use ecosys_diagnostics_mod, only : photoC_zint_diag_ind
  use ecosys_diagnostics_mod, only : photoC_NO3_zint_diag_ind
  use ecosys_diagnostics_mod, only : CaCO3_form_zint_diag_ind
  ! 3d
  use ecosys_diagnostics_mod, only : auto_diag_cnt_3d
  use ecosys_diagnostics_mod, only : N_lim_diag_ind
  use ecosys_diagnostics_mod, only : P_lim_diag_ind
  use ecosys_diagnostics_mod, only : Fe_lim_diag_ind
  use ecosys_diagnostics_mod, only : SiO3_lim_diag_ind
  use ecosys_diagnostics_mod, only : light_lim_diag_ind
  use ecosys_diagnostics_mod, only : photoC_diag_ind
  use ecosys_diagnostics_mod, only : photoC_NO3_diag_ind
  use ecosys_diagnostics_mod, only : photoFe_diag_ind
  use ecosys_diagnostics_mod, only : photoNO3_diag_ind
  use ecosys_diagnostics_mod, only : photoNH4_diag_ind
  use ecosys_diagnostics_mod, only : DOP_uptake_diag_ind
  use ecosys_diagnostics_mod, only : PO4_uptake_diag_ind
  use ecosys_diagnostics_mod, only : auto_graze_diag_ind
  use ecosys_diagnostics_mod, only : auto_graze_poc_diag_ind
  use ecosys_diagnostics_mod, only : auto_graze_doc_diag_ind
  use ecosys_diagnostics_mod, only : auto_graze_zoo_diag_ind
  use ecosys_diagnostics_mod, only : auto_loss_diag_ind
  use ecosys_diagnostics_mod, only : auto_loss_poc_diag_ind
  use ecosys_diagnostics_mod, only : auto_loss_doc_diag_ind
  use ecosys_diagnostics_mod, only : auto_agg_diag_ind
  use ecosys_diagnostics_mod, only : bSi_form_diag_ind
  use ecosys_diagnostics_mod, only : CaCO3_form_diag_ind
  use ecosys_diagnostics_mod, only : Nfix_diag_ind

  ! zooplankton diagnostics
  ! 2d
  use ecosys_diagnostics_mod, only : zoo_diag_cnt_2d
  ! 3d
  use ecosys_diagnostics_mod, only : zoo_diag_cnt_3d
  use ecosys_diagnostics_mod, only : zoo_loss_diag_ind
  use ecosys_diagnostics_mod, only : zoo_loss_poc_diag_ind
  use ecosys_diagnostics_mod, only : zoo_loss_doc_diag_ind
  use ecosys_diagnostics_mod, only : zoo_graze_diag_ind
  use ecosys_diagnostics_mod, only : zoo_graze_poc_diag_ind
  use ecosys_diagnostics_mod, only : zoo_graze_doc_diag_ind
  use ecosys_diagnostics_mod, only : zoo_graze_zoo_diag_ind
  use ecosys_diagnostics_mod, only : x_graze_zoo_diag_ind

  ! particulate diagnostics
  ! 2D
  use ecosys_diagnostics_mod, only : part_diag_cnt_2d
  use ecosys_diagnostics_mod, only : calcToSed_diag_ind
  use ecosys_diagnostics_mod, only : bsiToSed_diag_ind
  use ecosys_diagnostics_mod, only : pocToSed_diag_ind
  use ecosys_diagnostics_mod, only : SedDenitrif_diag_ind
  use ecosys_diagnostics_mod, only : OtherRemin_diag_ind
  use ecosys_diagnostics_mod, only : ponToSed_diag_ind
  use ecosys_diagnostics_mod, only : popToSed_diag_ind
  use ecosys_diagnostics_mod, only : dustToSed_diag_ind
  use ecosys_diagnostics_mod, only : pfeToSed_diag_ind
  ! 3D
  use ecosys_diagnostics_mod, only : part_diag_cnt_3d
  use ecosys_diagnostics_mod, only : POC_FLUX_IN_diag_ind
  use ecosys_diagnostics_mod, only : POC_PROD_diag_ind
  use ecosys_diagnostics_mod, only : POC_REMIN_diag_ind
  use ecosys_diagnostics_mod, only : CaCO3_FLUX_IN_diag_ind
  use ecosys_diagnostics_mod, only : CaCO3_PROD_diag_ind
  use ecosys_diagnostics_mod, only : CaCO3_REMIN_diag_ind
  use ecosys_diagnostics_mod, only : SiO2_FLUX_IN_diag_ind
  use ecosys_diagnostics_mod, only : SiO2_PROD_diag_ind
  use ecosys_diagnostics_mod, only : SiO2_REMIN_diag_ind
  use ecosys_diagnostics_mod, only : dust_FLUX_IN_diag_ind
  use ecosys_diagnostics_mod, only : dust_REMIN_diag_ind
  use ecosys_diagnostics_mod, only : P_iron_FLUX_IN_diag_ind
  use ecosys_diagnostics_mod, only : P_iron_PROD_diag_ind
  use ecosys_diagnostics_mod, only : P_iron_REMIN_diag_ind

  ! tavg_forcing diagnostics
  use ecosys_diagnostics_mod, only : forcing_diag_cnt
  use ecosys_diagnostics_mod, only : ECOSYS_IFRAC_diag_ind
  use ecosys_diagnostics_mod, only : ECOSYS_XKW_diag_ind
  use ecosys_diagnostics_mod, only : ECOSYS_ATM_PRESS_diag_ind
  use ecosys_diagnostics_mod, only : PV_O2_diag_ind
  use ecosys_diagnostics_mod, only : SCHMIDT_O2_diag_ind
  use ecosys_diagnostics_mod, only : O2SAT_diag_ind
  use ecosys_diagnostics_mod, only : O2_GAS_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : CO2STAR_diag_ind
  use ecosys_diagnostics_mod, only : DCO2STAR_diag_ind
  use ecosys_diagnostics_mod, only : pCO2SURF_diag_ind
  use ecosys_diagnostics_mod, only : DpCO2_diag_ind
  use ecosys_diagnostics_mod, only : PV_CO2_diag_ind
  use ecosys_diagnostics_mod, only : SCHMIDT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DIC_GAS_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : PH_diag_ind
  use ecosys_diagnostics_mod, only : ATM_CO2_diag_ind
  use ecosys_diagnostics_mod, only : CO2STAR_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DCO2STAR_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : pCO2SURF_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DpCO2_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DIC_GAS_FLUX_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : PH_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : ATM_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : IRON_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DUST_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : NOx_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : NHy_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIN_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIP_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DoN_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DoNr_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOP_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOPr_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DSI_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DFE_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIC_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : ALK_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOC_RIV_FLUX_diag_ind

  use marbl_share_mod, only : autotrophs
  use marbl_share_mod, only : zooplankton
  use marbl_share_mod, only : autotroph_cnt
  use marbl_share_mod, only : zooplankton_cnt

  use marbl_interface_types, only : marbl_diagnostics_type

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_tavg_init
  public :: ecosys_tavg_accumulate
  public :: ecosys_tavg_accumulate_flux

  !-----------------------------------------------------------------------
  !  define tavg id for everything but zooplankton, autotrophs, surface
  !  flux terms and particulate terms
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(ecosys_diag_cnt_2d) :: tavg_ecosys_2d
  integer (int_kind), dimension(ecosys_diag_cnt_3d) :: tavg_ecosys_3d

  !-----------------------------------------------------------------------
  !  tavg ids for 2d fields related to surface fluxes
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(forcing_diag_cnt) :: tavg_forcing
  integer (int_kind) ::      &
      tavg_ECOSYS_IFRAC_2,   &! ice fraction duplicate
      tavg_ECOSYS_XKW_2,     &! xkw duplicate
      tavg_O2_GAS_FLUX_2,    &! O2 flux duplicate
      tavg_DpCO2_2,          &! delta pco2 duplicate
      tavg_DIC_GAS_FLUX_2     ! dic flux duplicate

  !-----------------------------------------------------------------------
  !  define tavg id for MORE nonstandard 3d fields
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
      tavg_DONr_remin,  &! tavg id for DONrefractory remin
      tavg_DOPr_remin    ! tavg id for DOPrefractory remin

  !-----------------------------------------------------------------------
  !  tavg ids for particulate terms
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(part_diag_cnt_2d) :: tavg_part_2d
  integer (int_kind), dimension(part_diag_cnt_3d) :: tavg_part_3d
  integer (int_kind) :: tavg_POC_ACCUM      ! tavg id for poc accumulation

  !-----------------------------------------------------------------------
  !  tavg ids for zooplankton fields
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(zoo_diag_cnt_2d, zooplankton_cnt) :: tavg_zoo_2d
  integer (int_kind), dimension(zoo_diag_cnt_3d, zooplankton_cnt) :: tavg_zoo_3d

  !-----------------------------------------------------------------------
  !  tavg ids for autotroph fields
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(auto_diag_cnt_2d, autotroph_cnt) :: tavg_auto_2d
  integer (int_kind), dimension(auto_diag_cnt_3d, autotroph_cnt) :: tavg_auto_3d

  integer (int_kind) ::         &
      tavg_tot_bSi_form,        &! tavg id for Si uptake
      tavg_tot_CaCO3_form,      &! tavg id for CaCO3 formation
      tavg_tot_CaCO3_form_zint, &! tavg id for CaCO3 formation vertical integral
      tavg_tot_Nfix              ! tavg id for N fixation


contains

  subroutine ecosys_tavg_init(marbl_diags, ecosys_restore)
! !DESCRIPTION:
!  call define_tavg_field for all tavg fields
!
! !REVISION HISTORY:
!  same as module
!
    use ecosys_restore_mod, only : ecosys_restore_type

    type(marbl_diagnostics_type), intent(inout) :: marbl_diags
    type(ecosys_restore_type), intent(inout) :: ecosys_restore
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'

    integer (int_kind) :: &
      auto_ind,       & ! autotroph functional group index
      zoo_ind           ! zooplankton functional group index

    character(char_len) :: sname             ! short-name of tavg variable
    integer(int_kind) :: n
    logical :: define_field

!-----------------------------------------------------------------------
!   Define tavg fields for surface flux terms
!-----------------------------------------------------------------------

    associate(                                                                &
              diags_2d      => marbl_diags%diags_2d(:),                       &
              diags_3d      => marbl_diags%diags_3d(:),                       &
              auto_diags_2d => marbl_diags%auto_diags_2d(:,:),                &
              auto_diags_3d => marbl_diags%auto_diags_3d(:,:),                &
              zoo_diags_2d  => marbl_diags%zoo_diags_2d(:,:),                 &
              zoo_diags_3d  => marbl_diags%zoo_diags_3d(:,:),                 &
              part_diags_2d => marbl_diags%part_diags_2d(:),                  &
              part_diags_3d => marbl_diags%part_diags_3d(:),                  &
              restore_diags => marbl_diags%restore_diags(:)                   &
             )
    call define_tavg_field(tavg_forcing(ECOSYS_IFRAC_diag_ind),         &
                           'ECOSYS_IFRAC',2,                            &
                           long_name='Ice Fraction for ecosys fluxes',  &
                           units='fraction', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ECOSYS_IFRAC_2,'ECOSYS_IFRAC_2',2,      &
                           long_name='Ice Fraction for ecosys fluxes',  &
                           units='fraction', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(ECOSYS_XKW_diag_ind),           &
                           'ECOSYS_XKW',2,                              &
                           long_name='XKW for ecosys fluxes',           &
                           units='cm/s', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ECOSYS_XKW_2,'ECOSYS_XKW_2',2,          &
                           long_name='XKW for ecosys fluxes',           &
                           units='cm/s', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(ECOSYS_ATM_PRESS_diag_ind),     &
                           'ECOSYS_ATM_PRESS',2,                        &
                           long_name='Atmospheric Pressure for ecosys fluxes', &
                           units='atmospheres', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(PV_O2_diag_ind),'PV_O2',2,      &
                           long_name='PV_O2',                           &
                           units='cm/s', grid_loc='2110',               &
                          coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(SCHMIDT_O2_diag_ind),           &
                           'SCHMIDT_O2',2,                              &
                           long_name='O2 Schmidt Number',               &
                           units='none', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(O2SAT_diag_ind),'O2SAT',2,      &
                           long_name='O2 Saturation',                   &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_O2_GAS_FLUX_2,'STF_O2_2',2,             &
                           long_name='Dissolved Oxygen Surface Flux',   &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(CO2STAR_diag_ind),'CO2STAR',2,  &
                           long_name='CO2 Star',                        &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DCO2STAR_diag_ind),'DCO2STAR',2,&
                           long_name='D CO2 Star',                      &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(pCO2SURF_diag_ind),'pCO2SURF',2,&
                           long_name='surface pCO2',                    &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DpCO2_diag_ind),'DpCO2',2,      &
                           long_name='D pCO2',                          &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_DpCO2_2,'DpCO2_2',2,                    &
                           long_name='D pCO2',                          &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(PV_CO2_diag_ind),'PV_CO2',2,    &
                           long_name='CO2 Piston Velocity',             &
                           units='cm/s', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(SCHMIDT_CO2_diag_ind),          &
                           'SCHMIDT_CO2',2,                             &
                           long_name='CO2 Schmidt Number',              &
                           units='none', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DIC_GAS_FLUX_diag_ind),         &
                           'FG_CO2',2,                                  &
                           long_name='DIC Surface Gas Flux',            &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_DIC_GAS_FLUX_2,'FG_CO2_2',2,            &
                           long_name='DIC Surface Gas Flux',            &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(PH_diag_ind),'PH',2,            &
                           long_name='Surface pH',                      &
                           units='none', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(ATM_CO2_diag_ind),'ATM_CO2',2,  &
                           long_name='Atmospheric CO2',                 &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(CO2STAR_ALT_CO2_diag_ind),      &
                           'CO2STAR_ALT_CO2',2,                         &
                           long_name='CO2 Star, Alternative CO2',       &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DCO2STAR_ALT_CO2_diag_ind),     &
                           'DCO2STAR_ALT_CO2',2,                        &
                           long_name='D CO2 Star, Alternative CO2',     &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(pCO2SURF_ALT_CO2_diag_ind),     &
                           'pCO2SURF_ALT_CO2',2,                        &
                           long_name='surface pCO2, Alternative CO2',   &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DpCO2_ALT_CO2_diag_ind),        &
                           'DpCO2_ALT_CO2',2,                           &
                           long_name='D pCO2, Alternative CO2',         &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DIC_GAS_FLUX_ALT_CO2_diag_ind), &
                           'FG_ALT_CO2',2,                              &
                           long_name='DIC Surface Gas Flux, Alternative CO2', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(PH_ALT_CO2_diag_ind),           &
                           'PH_ALT_CO2',2,                              &
                           long_name='Surface pH, Alternative CO2',     &
                           units='none', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(ATM_ALT_CO2_diag_ind),          &
                           'ATM_ALT_CO2',2,                             &
                           long_name='Atmospheric Alternative CO2',     &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(IRON_FLUX_diag_ind),            &
                           'IRON_FLUX',2,                               &
                           long_name='Atmospheric Iron Flux',           &
                           units='mmol/m^2/s', grid_loc='2110',         &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DUST_FLUX_diag_ind),            &
                           'DUST_FLUX',2,                               &
                           long_name='Dust Flux',                       &
                           units='g/cm^2/s', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(NOx_FLUX_diag_ind),'NOx_FLUX',2,&
                           long_name='Flux of NOx from Atmosphere',     &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(NHy_FLUX_diag_ind),'NHy_FLUX',2,&
                           long_name='Flux of NHy from Atmosphere',     &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DIN_RIV_FLUX_diag_ind),         &
                           'DIN_RIV_FLUX',2,                            &
                           long_name='Flux of DIN from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DIP_RIV_FLUX_diag_ind),         &
                           'DIP_RIV_FLUX',2,                            &
                           long_name='Flux of DIP from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DON_RIV_FLUX_diag_ind),         &
                           'DON_RIV_FLUX',2,                            &
                           long_name='Flux of DON from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DONr_RIV_FLUX_diag_ind),        &
                           'DONr_RIV_FLUX',2,                           &
                           long_name='Flux of DONr from rivers',        &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DOP_RIV_FLUX_diag_ind),         &
                           'DOP_RIV_FLUX',2,                            &
                           long_name='Flux of DOP from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DOPr_RIV_FLUX_diag_ind),        &
                           'DOPr_RIV_FLUX',2,                           &
                           long_name='Flux of DOPr from rivers',        &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DSI_RIV_FLUX_diag_ind),         &
                           'DSI_RIV_FLUX',2,                            &
                           long_name='Flux of DSI from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DFE_RIV_FLUX_diag_ind),         &
                           'DFE_RIV_FLUX',2,                            &
                           long_name='Flux of DFE from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DIC_RIV_FLUX_diag_ind),         &
                           'DIC_RIV_FLUX',2,                            &
                           long_name='Flux of DIC from rivers',         &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(ALK_RIV_FLUX_diag_ind),         &
                           'ALK_RIV_FLUX',2,                            &
                           long_name='Flux of ALK from rivers',         &
                           units='alk/cm^2/s', grid_loc='2110',         &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DOC_RIV_FLUX_diag_ind),         &
                           'DOC_RIV_FLUX',2,                            &
                           long_name='Flux of DOC from rivers',         &
                           units='alk/cm^2/s', grid_loc='2110',         &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!   Define 2D tavg fields for everything but zooplankton, autotrophs,
!   surface flux terms, and particulate terms
!-----------------------------------------------------------------------

    do n=1,ecosys_diag_cnt_2d
      call define_tavg_field(tavg_ecosys_2d(n), trim(diags_2d(n)%short_name), &
                            2, long_name=trim(diags_2d(n)%long_name),         &
                            units=trim(diags_2d(n)%units),                    &
                            grid_loc='2110', coordinates='TLONG TLAT time')
    end do

!-----------------------------------------------------------------------
!   Define 3D tavg fields for everything but zooplankton, autotrophs,
!   surface flux terms, and particulate terms
!-----------------------------------------------------------------------

    call define_tavg_field(tavg_ecosys_3d(CO3_diag_ind),'CO3',3,           &
                           long_name='Carbonate Ion Concentration',     &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(HCO3_diag_ind),'HCO3',3,         &
                           long_name='Bicarbonate Ion Concentration',   &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(H2CO3_diag_ind),'H2CO3',3,       &
                           long_name='Carbonic Acid Concentration',     &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(pH_3D_diag_ind),'pH_3D',3,       &
                           long_name='pH',                              &
                           units='none', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(CO3_ALT_CO2_diag_ind),'CO3_ALT_CO2',3, &
                           long_name='Carbonate Ion Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(HCO3_ALT_CO2_diag_ind),'HCO3_ALT_CO2',3, &
                           long_name='Bicarbonate Ion Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(H2CO3_ALT_CO2_diag_ind),'H2CO3_ALT_CO2',3, &
                           long_name='Carbonic Acid Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                     &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(pH_3D_ALT_CO2_diag_ind),'pH_3D_ALT_CO2',3, &
                           long_name='pH, Alternative CO2',                       &
                           units='none', grid_loc='3111',                         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(co3_sat_calc_diag_ind),'co3_sat_calc',3, &
                           long_name='CO3 concentration at calcite saturation', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(co3_sat_arag_diag_ind),'co3_sat_arag',3, &
                           long_name='CO3 concentration at aragonite saturation', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(NITRIF_diag_ind),'NITRIF',3,     &
                           long_name='Nitrification',                   &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(DENITRIF_diag_ind),'DENITRIF',3, &
                           long_name='Denitrification',                 &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(O2_PRODUCTION_diag_ind),'O2_PRODUCTION',3, &
                           long_name='O2 Production',                             &
                           units='mmol/m^3/s', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(O2_CONSUMPTION_diag_ind),'O2_CONSUMPTION',3, &
                           long_name='O2 Consumption',                              &
                           units='mmol/m^3/s', grid_loc='3111',                     &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(AOU_diag_ind),'AOU',3,           &
                           long_name='Apparent O2 Utilization ',        &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(PAR_avg_diag_ind),'PAR_avg',3,   &
                           long_name='PAR Average over Model Cell',     &
                           units='w/m^2', grid_loc='3114',              &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(auto_graze_TOT_diag_ind),'graze_auto_TOT',3, &
                           long_name='Total Autotroph Grazing',                     &
                           units='mmol/m^3/s', grid_loc='3114',                     &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(photoC_TOT_diag_ind),'photoC_TOT',3, &
                           long_name='Total C Fixation',                    &
                           units='mmol/m^3/s', grid_loc='3114',             &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(photoC_NO3_TOT_diag_ind),'photoC_NO3_TOT',3, &
                           long_name='Total C Fixation from NO3',                   &
                           units='mmol/m^3/s', grid_loc='3114',                     &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DOC_prod_diag_ind),'DOC_prod',3, &
                           long_name='DOC Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(DOC_remin_diag_ind),'DOC_remin',3, &
                           long_name='DOC Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(DON_prod_diag_ind),'DON_prod',3, &
                           long_name='DON Production',                  &
                           units='mmol/m^3/s', grid_loc='3114',         &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DON_remin_diag_ind),'DON_remin',3, &
                           long_name='DON Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3114',           &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DOP_prod_diag_ind),'DOP_prod',3, &
                           long_name='DOP Production',                  &
                           units='mmol/m^3/s', grid_loc='3114',         &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DOP_remin_diag_ind),'DOP_remin',3, &
                           long_name='DOP Remineralization',            &
                           units='mmol/m^3/s', grid_loc='3114',         &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DOFe_prod_diag_ind),'DOFe_prod',3, &
                           long_name='DOFe Production',                   &
                           units='mmol/m^3/s', grid_loc='3114',           &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(DOFe_remin_diag_ind),               &
                           'DOFe_remin',3,                                    &
                           long_name='DOFe Remineralization',                 &
                           units='mmol/m^3/s', grid_loc='3114',               &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys_3d(Fe_scavenge_diag_ind),              &
                           'Fe_scavenge',3,                                   &
                           long_name='Iron Scavenging',                       &
                           units='mmol/m^3/s', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys_3d(Fe_scavenge_rate_diag_ind),         &
                           'Fe_scavenge_rate',3,                              &
                           long_name='Iron Scavenging Rate',                  &
                           units='1/y', grid_loc='3111',                      &
                           coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!  Define 2D tavg fields for particulate terms
!-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Vars to sum up burial in sediments and sed Denitrif N losses
    !-----------------------------------------------------------------------
    call define_tavg_field(tavg_part_2d(calcToSed_diag_ind),'calcToSed',2, &
                           long_name='CaCO3 Flux to Sediments',         &
                           units='nmolC/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(pocToSed_diag_ind),'pocToSed',2,   &
                           long_name='POC Flux to Sediments',           &
                           units='nmolC/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(ponToSed_diag_ind),'ponToSed',2,   &
                           long_name='nitrogen burial Flux to Sediments',&
                           units='nmolN/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(SedDenitrif_diag_ind),'SedDenitrif',2, &
                           long_name='nitrogen loss in Sediments',          &
                           units='nmolN/cm^2/s', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(OtherRemin_diag_ind),'OtherRemin',2, &
                           long_name='non-oxic,non-dentr remin in Sediments', &
                           units='nmolC/cm^2/s', grid_loc='2110',         &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(popToSed_diag_ind),'popToSed',2, &
                           long_name='phosporus Flux to Sediments',   &
                           units='nmolP/cm^2/s', grid_loc='2110',     &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(bsiToSed_diag_ind),'bsiToSed',2,   &
                           long_name='biogenic Si Flux to Sediments',   &
                           units='nmolSi/cm^2/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(dustToSed_diag_ind),'dustToSed',2, &
                           long_name='dust Flux to Sediments',          &
                           units='g/cm^2/s', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part_2d(pfeToSed_diag_ind),'pfeToSed',2,   &
                           long_name='pFe Flux to Sediments',           &
                           units='nmolFe/cm^2/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  Define 3D tavg fields for particulate terms
!-----------------------------------------------------------------------

    call define_tavg_field(tavg_part_3d(POC_FLUX_IN_diag_ind),'POC_FLUX_IN',3, &
                           long_name='POC Flux into Cell',                  &
                           units='mmol/m^3 cm/s', grid_loc='3111',          &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(POC_PROD_diag_ind),'POC_PROD',3,   &
                           long_name='POC Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(POC_REMIN_diag_ind),'POC_REMIN',3, &
                           long_name='POC Remineralization',            &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(CaCO3_FLUX_IN_diag_ind),'CaCO3_FLUX_IN',3, &
                           long_name='CaCO3 flux into cell',                    &
                           units='mmol/m^3 cm/s', grid_loc='3111',              &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(CaCO3_PROD_diag_ind),'CaCO3_PROD',3, &
                           long_name='CaCO3 Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(CaCO3_REMIN_diag_ind),'CaCO3_REMIN',3, &
                           long_name='CaCO3 Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(SiO2_FLUX_IN_diag_ind),'SiO2_FLUX_IN',3, &
                           long_name='SiO2 Flux into Cell',                   &
                           units='mmol/m^3 cm/s', grid_loc='3111',            &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(SiO2_PROD_diag_ind),'SiO2_PROD',3, &
                           long_name='SiO2 Production',                 &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(SiO2_REMIN_diag_ind),'SiO2_REMIN',3, &
                           long_name='SiO2 Remineralization',             &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(dust_FLUX_IN_diag_ind),'dust_FLUX_IN',3, &
                           long_name='Dust Flux into Cell',                   &
                           units='ng/s/m^2', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(dust_REMIN_diag_ind),'dust_REMIN',3, &
                           long_name='Dust Remineralization',             &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(P_iron_FLUX_IN_diag_ind),'P_iron_FLUX_IN',3, &
                           long_name='P_iron Flux into Cell',                     &
                           units='mmol/m^3 cm/s', grid_loc='3111',                &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(P_iron_PROD_diag_ind),'P_iron_PROD',3, &
                           long_name='P_iron Production',                   &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part_3d(P_iron_REMIN_diag_ind),'P_iron_REMIN',3, &
                           long_name='P_iron Remineralization',               &
                           units='mmol/m^3/s', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_POC_ACCUM,'POC_ACCUM',3,                &
                           long_name='POC Accumulation',                &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!  Define 2D tavg fields for zooplankton
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Define 3D tavg fields for zooplankton
!-----------------------------------------------------------------------

    do zoo_ind = 1, zooplankton_cnt
       call define_tavg_field(tavg_zoo_3d(zoo_loss_diag_ind,zoo_ind),      &
            trim(zooplankton(zoo_ind)%sname) // '_loss', 3,             &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss',      &
            units='mmol/m^3/s', grid_loc='3114',                        &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_loss_poc_diag_ind,zoo_ind),     &
            trim(zooplankton(zoo_ind)%sname) // '_loss_poc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss to POC',  &
            units='mmol/m^3/s', grid_loc='3114',                           &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_loss_doc_diag_ind,zoo_ind),     &
            trim(zooplankton(zoo_ind)%sname) // '_loss_doc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss to DOC',  &
            units='mmol/m^3/s', grid_loc='3114',                           &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_graze_diag_ind,zoo_ind),          &
            'graze_' // trim(zooplankton(zoo_ind)%sname), 3,                 &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss',   &
            units='mmol/m^3/s', grid_loc='3114',                             &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_graze_poc_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_poc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to POC', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_graze_doc_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_doc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to DOC', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(zoo_graze_zoo_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_zoo', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to ZOO', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo_3d(x_graze_zoo_diag_ind,zoo_ind),       &
            'x_graze_' // trim(zooplankton(zoo_ind)%sname), 3,              &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing gain',  &
            units='mmol/m^3/s', grid_loc='3114',                            &
            coordinates='TLONG TLAT z_t_150m time')
    end do


!-----------------------------------------------------------------------
!  Define 2D tavg fields for autotrophs
!-----------------------------------------------------------------------

    do auto_ind = 1, autotroph_cnt
      do n=1,auto_diag_cnt_2d
        define_field = .true.
        if (n.eq.CaCO3_form_zint_diag_ind) then
          define_field = (autotrophs(auto_ind)%CaCO3_ind > 0)
        end if

        if (define_field) then
          call define_tavg_field(tavg_auto_2d(n, auto_ind),                   &
                        trim(auto_diags_2d(n, auto_ind)%short_name), 2,       &
                        long_name=trim(auto_diags_2d(n, auto_ind)%long_name), &
                        units=trim(auto_diags_2d(n, auto_ind)%units),         &
                        grid_loc='2110', coordinates='TLONG TLAT time')
        end if
      end do
    end do

    call define_tavg_field(tavg_tot_CaCO3_form_zint, 'CaCO3_form_zint', 2, &
                           long_name='Total CaCO3 Formation Vertical Integral', &
                           units='mmol/m^3 cm/s', grid_loc='2110', &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  Define 3D tavg fields for autotrophs
!-----------------------------------------------------------------------

    do auto_ind = 1, autotroph_cnt
       call define_tavg_field(tavg_auto_3d(N_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_N_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' N Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(P_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_P_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' P Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(Fe_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_Fe_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Fe Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(SiO3_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_SiO3_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' SiO3 Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(light_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_light_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Light Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(photoC_diag_ind,auto_ind), &
                              'photoC_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(photoC_NO3_diag_ind,auto_ind), &
                              'photoC_NO3_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation from NO3', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(photoFe_diag_ind,auto_ind), &
                              'photoFe_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Fe Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(photoNO3_diag_ind,auto_ind), &
                              'photoNO3_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' NO3 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(photoNH4_diag_ind,auto_ind), &
                              'photoNH4_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' NH4 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(DOP_uptake_diag_ind,auto_ind), &
                              'DOP_' // trim(autotrophs(auto_ind)%sname) // '_uptake', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' DOP Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(PO4_uptake_diag_ind,auto_ind), &
                              'PO4_' // trim(autotrophs(auto_ind)%sname) // '_uptake', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' PO4 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_graze_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_graze_poc_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_poc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to POC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_graze_doc_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_doc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to DOC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_graze_zoo_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_zoo', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to ZOO', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_loss_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_loss_poc_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss_poc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss to POC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_loss_doc_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss_doc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss to DOC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto_3d(auto_agg_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_agg', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Aggregate', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       if (autotrophs(auto_ind)%Si_ind > 0) then
          sname = trim(autotrophs(auto_ind)%sname) // 'bSi_form'
          call define_tavg_field(tavg_auto_3d(bSi_form_diag_ind,auto_ind), sname, 3, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' Si Uptake', &
                                 units='mmol/m^3/s', grid_loc='3114', &
                                 coordinates='TLONG TLAT z_t_150m time')
       endif

       if (autotrophs(auto_ind)%CaCO3_ind > 0) then
          sname = trim(autotrophs(auto_ind)%sname) // '_CaCO3_form'
          call define_tavg_field(tavg_auto_3d(CaCO3_form_diag_ind,auto_ind), sname, 3, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' CaCO3 Formation', &
                                 units='mmol/m^3/s', grid_loc='3114', &
                                 coordinates='TLONG TLAT z_t_150m time')
       endif

       if (autotrophs(auto_ind)%Nfixer) then
          call define_tavg_field(tavg_auto_3d(Nfix_diag_ind,auto_ind), &
                                 trim(autotrophs(auto_ind)%sname) // '_Nfix', 3, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' N Fixation', &
                                 units='mmol/m^3/s', grid_loc='3114',   &
                                 coordinates='TLONG TLAT z_t_150m time')
       endif
    end do

    !-----------------------------------------------------------------------
    !  Define tavg fields for sum over all autotrophs
    !-----------------------------------------------------------------------

    call define_tavg_field(tavg_tot_bSi_form, 'bSi_form', 3, &
                           long_name='Total Si Uptake', &
                           units='mmol/m^3/s', grid_loc='3114', &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_tot_CaCO3_form, 'CaCO3_form', 3, &
                           long_name='Total CaCO3 Formation', &
                           units='mmol/m^3/s', grid_loc='3114', &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_tot_Nfix, 'Nfix', 3, &
                           long_name='Total N Fixation', &
                           units='mmol/m^3/s', grid_loc='3114',   &
                           coordinates='TLONG TLAT z_t_150m time')

!  DEFINED BUT NEVER WRITTEN
    call define_tavg_field(tavg_DONr_REMIN,'DONr_REMIN',3,              &
                           long_name='DONr Remineralization',           &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_DOPr_REMIN,'DOPr_REMIN',3,              &
                           long_name='DOPr Remineralization',           &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call ecosys_restore%define_tavg_fields()
    end associate

!-----------------------------------------------------------------------
!EOC

  end subroutine ecosys_tavg_init

  subroutine ecosys_tavg_accumulate(i, c,bid, marbl_diagnostics,              &
                                    ecosys_restore)

    use ecosys_restore_mod, only : ecosys_restore_type

    integer, intent(in) :: i, c, bid ! column indices and block index
    type(marbl_diagnostics_type), intent(in) :: marbl_diagnostics
    type(ecosys_restore_type), intent(in) :: ecosys_restore

    integer :: n, auto_ind, zoo_ind, k
    logical :: accumulate

    associate(                                                                &
              DIAGS_2D      => marbl_diagnostics%diags_2d(:)%field,           &
              DIAGS_3D      => marbl_diagnostics%diags_3d(:),                 &
              AUTO_DIAGS_2D => marbl_diagnostics%auto_diags_2d(:,:)%field,    &
              AUTO_DIAGS_3D => marbl_diagnostics%auto_diags_3d(:,:),          &
              ZOO_DIAGS_2D  => marbl_diagnostics%zoo_diags_2d(:,:)%field,     &
              ZOO_DIAGS_3D  => marbl_diagnostics%zoo_diags_3d(:,:),           &
              PART_DIAGS_2D => marbl_diagnostics%part_diags_2d(:)%field,      &
              PART_DIAGS_3D => marbl_diagnostics%part_diags_3d(:),            &
              restore_diags => marbl_diagnostics%restore_diags(:)             &
             )

    ! Accumulate general diagnostics
    ! 2D
    do n=1,ecosys_diag_cnt_2d
      call accumulate_tavg_field(DIAGS_2D(n),tavg_ecosys_2d(n),bid,i,c)
    end do

    ! 3D
    do n=1,ecosys_diag_cnt_3d
      call accumulate_tavg_field(DIAGS_3D(n)%field(:),tavg_ecosys_3d(n),bid,i,c)
    end do

    ! Accumulate autotroph terms
    ! 2D autotrophs
    do n=1,auto_diag_cnt_2d
      do auto_ind=1,autotroph_cnt
        accumulate = .true.
        if (n.eq.CaCO3_form_zint_diag_ind) then
          accumulate = (autotrophs(auto_ind)%imp_calcifier)
          if ( accumulate) then
            call accumulate_tavg_field(AUTO_DIAGS_2D(n,auto_ind),             &
                                       tavg_tot_CaCO3_form_zint, bid, i, c)
          end if
        end if

        if (accumulate) then
          call accumulate_tavg_field(AUTO_DIAGS_2D(n,auto_ind),               &
                                     tavg_auto_2d(n,auto_ind), bid, i, c)
        end if
      end do
    end do

    ! 3D autotrophs
    do n=1,auto_diag_cnt_3d
      do auto_ind=1,autotroph_cnt
        accumulate = .true.
        ! Some autotrophs are only accumulated under specific conditions
        select case (n)
          case (bSi_form_diag_ind)
            accumulate = (autotrophs(auto_ind)%Si_ind.gt.0)
            if ( accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS_3D(n,auto_ind)%field(:),  &
                                     tavg_tot_bSi_form, bid, i, c)
          case (CaCO3_form_diag_ind) 
            accumulate = (autotrophs(auto_ind)%imp_calcifier)
            if (accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS_3D(n,auto_ind)%field(:),  &
                                     tavg_tot_CaCO3_form, bid, i, c)
          case (Nfix_diag_ind)
            accumulate = (autotrophs(auto_ind)%Nfixer)
            if (accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS_3D(n,auto_ind)%field(:),  &
                                     tavg_tot_Nfix, bid, i, c)
        end select

        if (accumulate) &
            call accumulate_tavg_field(AUTO_DIAGS_3D(n,auto_ind)%field(:),    &
                                       tavg_auto_3d(n,auto_ind), bid, i, c)
      end do
     end do

    ! Accumulate zooplankton terms
    ! 3D
    do n=1,zoo_diag_cnt_3d
      do zoo_ind=1,zooplankton_cnt
        call accumulate_tavg_field(ZOO_DIAGS_3d(n,zoo_ind)%field(:),           &
                                   tavg_zoo_3d(n,zoo_ind), bid, i, c)
      end do
    end do

    ! Accumulate particulate terms
    ! 2D
    do n=1,part_diag_cnt_2d
      call accumulate_tavg_field(PART_DIAGS_2D(n), tavg_part_2d(n), bid, i, c)
    end do

    ! 3D
    do n=1,part_diag_cnt_3d
      call accumulate_tavg_field(PART_DIAGS_3D(n)%field(:), tavg_part_3d(n),  &
                                 bid, i, c)
    end do

!    call ecosys_restore%accumulate_tavg(restore_diags, bid, i, c)

    end associate

  end subroutine ecosys_tavg_accumulate

  subroutine ecosys_tavg_accumulate_flux(FLUX_DIAGS)

    real (r8), dimension(nx_block,ny_block, forcing_diag_cnt, nblocks_clinic),   &
                                                                   intent(in) :: &
               FLUX_DIAGS                ! Computed diagnostics for surface fluxes

    integer :: i, iblock

    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock=1,nblocks_clinic
      do i=1,forcing_diag_cnt
        associate(FLUX_FIELD => FLUX_DIAGS(:,:,i,iblock))
          ! Note: for some reason, there is no tavg_O2_GAS_FLUX variable...
          if (i.ne.O2_GAS_FLUX_diag_ind) then
            call accumulate_tavg_field(FLUX_FIELD, tavg_forcing(i), iblock, i)
          else
            call accumulate_tavg_field(FLUX_FIELD, tavg_O2_GAS_FLUX_2, iblock, i)
          end if

          select case (i)
            case (ECOSYS_IFRAC_diag_ind)
              call accumulate_tavg_field(FLUX_FIELD, tavg_ECOSYS_IFRAC_2,     &
                                         iblock, i)
            case (ECOSYS_XKW_diag_ind)
              call accumulate_tavg_field(FLUX_FIELD, tavg_ECOSYS_XKW_2,       &
                                         iblock, i)
            case (DpCO2_diag_ind)
              call accumulate_tavg_field(FLUX_FIELD, tavg_DpCO2_2, iblock, i)
            case (DIC_GAS_FLUX_diag_ind)
              call accumulate_tavg_field(FLUX_FIELD, tavg_DIC_GAS_FLUX_2,     &
                                         iblock, i)
          end select
        end associate
      end do
    end do
   !$OMP END PARALLEL DO

  end subroutine ecosys_tavg_accumulate_flux

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


