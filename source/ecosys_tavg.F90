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
  use ecosys_constants, only : ecosys_diag_cnt
  use ecosys_constants, only : CO3_diag_ind
  use ecosys_constants, only : HCO3_diag_ind
  use ecosys_constants, only : H2CO3_diag_ind
  use ecosys_constants, only : pH_3D_diag_ind
  use ecosys_constants, only : CO3_ALT_CO2_diag_ind
  use ecosys_constants, only : HCO3_ALT_CO2_diag_ind
  use ecosys_constants, only : H2CO3_ALT_CO2_diag_ind
  use ecosys_constants, only : ph_3D_ALT_CO2_diag_ind
  use ecosys_constants, only : co3_sat_calc_diag_ind
  use ecosys_constants, only : co3_sat_arag_diag_ind
  use ecosys_constants, only : zsatcalc_diag_ind
  use ecosys_constants, only : zsatarag_diag_ind
  use ecosys_constants, only : NITRIF_diag_ind
  use ecosys_constants, only : DENITRIF_diag_ind
  use ecosys_constants, only : O2_ZMIN_diag_ind
  use ecosys_constants, only : O2_ZMIN_DEPTH_diag_ind
  use ecosys_constants, only : O2_PRODUCTION_diag_ind
  use ecosys_constants, only : O2_CONSUMPTION_diag_ind
  use ecosys_constants, only : AOU_diag_ind
  use ecosys_constants, only : PAR_avg_diag_ind
  use ecosys_constants, only : auto_graze_TOT_diag_ind
  use ecosys_constants, only : photoC_TOT_diag_ind
  use ecosys_constants, only : photoC_TOT_zint_diag_ind
  use ecosys_constants, only : photoC_NO3_TOT_diag_ind
  use ecosys_constants, only : photoC_NO3_TOT_zint_diag_ind
  use ecosys_constants, only : DOC_prod_diag_ind
  use ecosys_constants, only : DOC_remin_diag_ind
  use ecosys_constants, only : DOCr_remin_diag_ind
  use ecosys_constants, only : DON_prod_diag_ind
  use ecosys_constants, only : DON_remin_diag_ind
  use ecosys_constants, only : DONr_remin_diag_ind
  use ecosys_constants, only : DOP_prod_diag_ind
  use ecosys_constants, only : DOP_remin_diag_ind
  use ecosys_constants, only : DOPr_remin_diag_ind
  use ecosys_constants, only : Fe_scavenge_diag_ind
  use ecosys_constants, only : Fe_scavenge_rate_diag_ind
  use ecosys_constants, only : Jint_Ctot_diag_ind
  use ecosys_constants, only : Jint_100m_Ctot_diag_ind
  use ecosys_constants, only : Jint_Ntot_diag_ind
  use ecosys_constants, only : Jint_100m_Ntot_diag_ind
  use ecosys_constants, only : Jint_Ptot_diag_ind
  use ecosys_constants, only : Jint_100m_Ptot_diag_ind
  use ecosys_constants, only : Jint_Sitot_diag_ind
  use ecosys_constants, only : Jint_100m_Sitot_diag_ind
  use ecosys_constants, only : Jint_Fetot_diag_ind
  use ecosys_constants, only : Jint_100m_Fetot_diag_ind
  ! autotroph diagnostics
  use ecosys_constants, only : auto_diag_cnt
  use ecosys_constants, only : N_lim_diag_ind
  use ecosys_constants, only : P_lim_diag_ind
  use ecosys_constants, only : Fe_lim_diag_ind
  use ecosys_constants, only : SiO3_lim_diag_ind
  use ecosys_constants, only : light_lim_diag_ind
  use ecosys_constants, only : photoC_diag_ind
  use ecosys_constants, only : photoC_zint_diag_ind
  use ecosys_constants, only : photoC_NO3_diag_ind
  use ecosys_constants, only : photoC_NO3_zint_diag_ind
  use ecosys_constants, only : photoFe_diag_ind
  use ecosys_constants, only : photoNO3_diag_ind
  use ecosys_constants, only : photoNH4_diag_ind
  use ecosys_constants, only : DOP_uptake_diag_ind
  use ecosys_constants, only : PO4_uptake_diag_ind
  use ecosys_constants, only : auto_graze_diag_ind
  use ecosys_constants, only : auto_graze_poc_diag_ind
  use ecosys_constants, only : auto_graze_doc_diag_ind
  use ecosys_constants, only : auto_graze_zoo_diag_ind
  use ecosys_constants, only : auto_loss_diag_ind
  use ecosys_constants, only : auto_loss_poc_diag_ind
  use ecosys_constants, only : auto_loss_doc_diag_ind
  use ecosys_constants, only : auto_agg_diag_ind
  use ecosys_constants, only : bSi_form_diag_ind
  use ecosys_constants, only : CaCO3_form_diag_ind
  use ecosys_constants, only : CaCO3_form_zint_diag_ind
  use ecosys_constants, only : Nfix_diag_ind
  ! zooplankton diagnostics
  use ecosys_constants, only : zoo_diag_cnt
  use ecosys_constants, only : zoo_loss_diag_ind
  use ecosys_constants, only : zoo_loss_poc_diag_ind
  use ecosys_constants, only : zoo_loss_doc_diag_ind
  use ecosys_constants, only : zoo_graze_diag_ind
  use ecosys_constants, only : zoo_graze_poc_diag_ind
  use ecosys_constants, only : zoo_graze_doc_diag_ind
  use ecosys_constants, only : zoo_graze_zoo_diag_ind
  use ecosys_constants, only : x_graze_zoo_diag_ind
  ! particulate diagnostics
  use ecosys_constants, only : part_diag_cnt
  use ecosys_constants, only : POC_FLUX_IN_diag_ind
  use ecosys_constants, only : POC_PROD_diag_ind
  use ecosys_constants, only : POC_REMIN_diag_ind
  use ecosys_constants, only : POC_REMIN_DIC_diag_ind
  use ecosys_constants, only : PON_REMIN_NH4_diag_ind
  use ecosys_constants, only : POP_REMIN_PO4_diag_ind
  use ecosys_constants, only : CaCO3_FLUX_IN_diag_ind
  use ecosys_constants, only : CaCO3_PROD_diag_ind
  use ecosys_constants, only : CaCO3_REMIN_diag_ind
  use ecosys_constants, only : SiO2_FLUX_IN_diag_ind
  use ecosys_constants, only : SiO2_PROD_diag_ind
  use ecosys_constants, only : SiO2_REMIN_diag_ind
  use ecosys_constants, only : dust_FLUX_IN_diag_ind
  use ecosys_constants, only : dust_REMIN_diag_ind
  use ecosys_constants, only : P_iron_FLUX_IN_diag_ind
  use ecosys_constants, only : P_iron_PROD_diag_ind
  use ecosys_constants, only : P_iron_REMIN_diag_ind
  use ecosys_constants, only : calcToSed_diag_ind
  use ecosys_constants, only : bsiToSed_diag_ind
  use ecosys_constants, only : pocToSed_diag_ind
  use ecosys_constants, only : SedDenitrif_diag_ind
  use ecosys_constants, only : OtherRemin_diag_ind
  use ecosys_constants, only : ponToSed_diag_ind
  use ecosys_constants, only : popToSed_diag_ind
  use ecosys_constants, only : dustToSed_diag_ind
  use ecosys_constants, only : pfeToSed_diag_ind
  ! tavg_forcing diagnostics
  use ecosys_constants, only : forcing_diag_cnt
  use ecosys_constants, only : ECOSYS_IFRAC_diag_ind
  use ecosys_constants, only : ECOSYS_XKW_diag_ind
  use ecosys_constants, only : ECOSYS_ATM_PRESS_diag_ind
  use ecosys_constants, only : PV_O2_diag_ind
  use ecosys_constants, only : SCHMIDT_O2_diag_ind
  use ecosys_constants, only : O2SAT_diag_ind
  use ecosys_constants, only : O2_GAS_FLUX_diag_ind
  use ecosys_constants, only : CO2STAR_diag_ind
  use ecosys_constants, only : DCO2STAR_diag_ind
  use ecosys_constants, only : pCO2SURF_diag_ind
  use ecosys_constants, only : DpCO2_diag_ind
  use ecosys_constants, only : PV_CO2_diag_ind
  use ecosys_constants, only : SCHMIDT_CO2_diag_ind
  use ecosys_constants, only : DIC_GAS_FLUX_diag_ind
  use ecosys_constants, only : PH_diag_ind
  use ecosys_constants, only : ATM_CO2_diag_ind
  use ecosys_constants, only : CO2STAR_ALT_CO2_diag_ind
  use ecosys_constants, only : DCO2STAR_ALT_CO2_diag_ind
  use ecosys_constants, only : pCO2SURF_ALT_CO2_diag_ind
  use ecosys_constants, only : DpCO2_ALT_CO2_diag_ind
  use ecosys_constants, only : DIC_GAS_FLUX_ALT_CO2_diag_ind
  use ecosys_constants, only : PH_ALT_CO2_diag_ind
  use ecosys_constants, only : ATM_ALT_CO2_diag_ind
  use ecosys_constants, only : IRON_FLUX_diag_ind
  use ecosys_constants, only : DUST_FLUX_diag_ind
  use ecosys_constants, only : NOx_FLUX_diag_ind
  use ecosys_constants, only : NHy_FLUX_diag_ind
  use ecosys_constants, only : DIN_RIV_FLUX_diag_ind
  use ecosys_constants, only : DIP_RIV_FLUX_diag_ind
  use ecosys_constants, only : DON_RIV_FLUX_diag_ind
  use ecosys_constants, only : DONr_RIV_FLUX_diag_ind
  use ecosys_constants, only : DOP_RIV_FLUX_diag_ind
  use ecosys_constants, only : DOPr_RIV_FLUX_diag_ind
  use ecosys_constants, only : DSI_RIV_FLUX_diag_ind
  use ecosys_constants, only : DFE_RIV_FLUX_diag_ind
  use ecosys_constants, only : DIC_RIV_FLUX_diag_ind
  use ecosys_constants, only : ALK_RIV_FLUX_diag_ind
  use ecosys_constants, only : DOC_RIV_FLUX_diag_ind
  use ecosys_constants, only : DOCr_RIV_FLUX_diag_ind

  use marbl_share_mod, only : autotrophs
  use marbl_share_mod, only : zooplankton
  use marbl_share_mod, only : autotroph_cnt
  use marbl_share_mod, only : zooplankton_cnt
  use marbl_interface_types, only : ecosys_diagnostics_type

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_tavg_init
  public :: ecosys_tavg_write
  public :: ecosys_tavg_write_flux

  !-----------------------------------------------------------------------
  !  define tavg id for everything but zooplankton, autotrophs, surface
  !  flux terms and particulate terms
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(ecosys_diag_cnt) :: tavg_ecosys

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
  !  tavg ids for particulate terms
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(part_diag_cnt) :: tavg_part

  !-----------------------------------------------------------------------
  !  tavg ids for zooplankton fields
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(zoo_diag_cnt, zooplankton_cnt) :: tavg_zoo

  !-----------------------------------------------------------------------
  !  tavg ids for autotroph fields
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(auto_diag_cnt, autotroph_cnt) :: tavg_auto

  integer (int_kind) ::         &
      tavg_tot_bSi_form,        &! tavg id for Si uptake
      tavg_tot_CaCO3_form,      &! tavg id for CaCO3 formation
      tavg_tot_CaCO3_form_zint, &! tavg id for CaCO3 formation vertical integral
      tavg_tot_Nfix              ! tavg id for N fixation


contains

  subroutine ecosys_tavg_init(ecosys_restore)
! !DESCRIPTION:
!  call define_tavg_field for all tavg fields
!
! !REVISION HISTORY:
!  same as module
!
    use ecosys_restore_mod, only : ecosys_restore_type

    type(ecosys_restore_type), intent(inout) :: ecosys_restore
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'

    integer (int_kind) :: &
      auto_ind,       & ! autotroph functional group index
      zoo_ind           ! zooplankton functional group index

    character(char_len) :: sname             ! short-name of tavg variable

!-----------------------------------------------------------------------
!   Define tavg fields for surface flux terms
!-----------------------------------------------------------------------

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
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_forcing(DOCr_RIV_FLUX_diag_ind),        &
                           'DOCr_RIV_FLUX',2,                           &
                           long_name='Flux of DOCr from rivers',        &
                           units='nmol/cm^2/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!   Define tavg fields for everything but zooplankton, autotrophs,
!   surface flux terms, and particulate terms
!-----------------------------------------------------------------------

    call define_tavg_field(tavg_ecosys(CO3_diag_ind),'CO3',3,           &
                           long_name='Carbonate Ion Concentration',     &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(HCO3_diag_ind),'HCO3',3,         &
                           long_name='Bicarbonate Ion Concentration',   &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(H2CO3_diag_ind),'H2CO3',3,       &
                           long_name='Carbonic Acid Concentration',     &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(pH_3D_diag_ind),'pH_3D',3,       &
                           long_name='pH',                              &
                           units='none', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(CO3_ALT_CO2_diag_ind),'CO3_ALT_CO2',3, &
                           long_name='Carbonate Ion Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(HCO3_ALT_CO2_diag_ind),'HCO3_ALT_CO2',3, &
                           long_name='Bicarbonate Ion Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(H2CO3_ALT_CO2_diag_ind),'H2CO3_ALT_CO2',3, &
                           long_name='Carbonic Acid Concentration, Alternative CO2', &
                           units='mmol/m^3', grid_loc='3111',                     &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(pH_3D_ALT_CO2_diag_ind),'pH_3D_ALT_CO2',3, &
                           long_name='pH, Alternative CO2',                       &
                           units='none', grid_loc='3111',                         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(co3_sat_calc_diag_ind),'co3_sat_calc',3, &
                           long_name='CO3 concentration at calcite saturation', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(zsatcalc_diag_ind),'zsatcalc',2, &
                           long_name='Calcite Saturation Depth',        &
                           units='cm', grid_loc='2110',                 &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(co3_sat_arag_diag_ind),'co3_sat_arag',3, &
                           long_name='CO3 concentration at aragonite saturation', &
                           units='mmol/m^3', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(zsatarag_diag_ind),'zsatarag',2, &
                           long_name='Aragonite Saturation Depth',      &
                           units='cm', grid_loc='2110',                 &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(NITRIF_diag_ind),'NITRIF',3,     &
                           long_name='Nitrification',                   &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DENITRIF_diag_ind),'DENITRIF',3, &
                           long_name='Denitrification',                 &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(O2_ZMIN_diag_ind),'O2_ZMIN',2,   &
                           long_name='Vertical Minimum of O2',          &
                           units='mmol/m^3', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(O2_ZMIN_DEPTH_diag_ind),'O2_ZMIN_DEPTH',2, &
                           long_name='Depth of Vertical Minimum of O2',           &
                           units='cm', grid_loc='2110',                           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(O2_PRODUCTION_diag_ind),'O2_PRODUCTION',3, &
                           long_name='O2 Production',                             &
                           units='mmol/m^3/s', grid_loc='3111',                   &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(O2_CONSUMPTION_diag_ind),'O2_CONSUMPTION',3, &
                           long_name='O2 Consumption',                              &
                           units='mmol/m^3/s', grid_loc='3111',                     &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(AOU_diag_ind),'AOU',3,           &
                           long_name='Apparent O2 Utilization ',        &
                           units='mmol/m^3', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(PAR_avg_diag_ind),'PAR_avg',3,   &
                           long_name='PAR Average over Model Cell',     &
                           units='w/m^2', grid_loc='3114',              &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys(auto_graze_TOT_diag_ind),'graze_auto_TOT',3, &
                           long_name='Total Autotroph Grazing',                     &
                           units='mmol/m^3/s', grid_loc='3114',                     &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys(photoC_TOT_diag_ind),'photoC_TOT',3, &
                           long_name='Total C Fixation',                    &
                           units='mmol/m^3/s', grid_loc='3114',             &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys(photoC_TOT_zint_diag_ind),'photoC_TOT_zint',2, &
                           long_name='Total C Fixation Vertical Integral',            &
                           units='mmol/m^3 cm/s', grid_loc='2110',                    &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(photoC_NO3_TOT_diag_ind),'photoC_NO3_TOT',3, &
                           long_name='Total C Fixation from NO3',                   &
                           units='mmol/m^3/s', grid_loc='3114',                     &
                           coordinates='TLONG TLAT z_t_150m time')

    call define_tavg_field(tavg_ecosys(photoC_NO3_TOT_zint_diag_ind),'photoC_NO3_TOT_zint',2,&
                           long_name='Total C Fixation from NO3 Vertical Integral',&
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(DOC_prod_diag_ind),'DOC_prod',3, &
                           long_name='DOC Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DOC_remin_diag_ind),'DOC_remin',3, &
                           long_name='DOC Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DOCr_remin_diag_ind),'DOCr_remin',3, &
                           long_name='DOCr Remineralization',               &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DON_prod_diag_ind),'DON_prod',3, &
                           long_name='DON Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DON_remin_diag_ind),'DON_remin',3, &
                           long_name='DON Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DONr_remin_diag_ind),'DONr_remin',3, &
                           long_name='DONr Remineralization',               &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DOP_prod_diag_ind),'DOP_prod',3, &
                           long_name='DOP Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DOP_remin_diag_ind),'DOP_remin',3, &
                           long_name='DOP Remineralization',            &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(DOPr_remin_diag_ind),'DOPr_remin',3, &
                           long_name='DOPr Remineralization',               &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(Fe_scavenge_diag_ind),'Fe_scavenge',3, &
                           long_name='Iron Scavenging',                       &
                           units='mmol/m^3/s', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_ecosys(Fe_scavenge_rate_diag_ind),'Fe_scavenge_rate',3, &
                           long_name='Iron Scavenging Rate',                            &
                           units='1/y', grid_loc='3111',                                &
                           coordinates='TLONG TLAT z_t time')

    !-----------------------------------------------------------------------
    !  Define tavg for fields related to conservation of total C, N, P, Si
    !-----------------------------------------------------------------------
    call define_tavg_field(tavg_ecosys(Jint_Ctot_diag_ind),'Jint_Ctot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ctot', &
                           units='mmol/m^3 cm/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_100m_Ctot_diag_ind),'Jint_100m_Ctot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ctot, 0-100m', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_Ntot_diag_ind),'Jint_Ntot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ntot', &
                           units='mmol/m^3 cm/s', grid_loc='2110',        &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_100m_Ntot_diag_ind),'Jint_100m_Ntot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ntot, 0-100m', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_Ptot_diag_ind),'Jint_Ptot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ptot', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_100m_Ptot_diag_ind),'Jint_100m_Ptot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Ptot, 0-100m', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_Sitot_diag_ind),'Jint_Sitot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Sitot', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_100m_Sitot_diag_ind),'Jint_100m_Sitot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Sitot, 0-100m', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_Fetot_diag_ind),'Jint_Fetot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Fetot', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ecosys(Jint_100m_Fetot_diag_ind),'Jint_100m_Fetot',2, &
                           long_name='Vertical Integral of Conservative Subterms of Source Sink Term for Fetot, 0-100m', &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  Define tavg fields for particulate terms
!-----------------------------------------------------------------------

    call define_tavg_field(tavg_part(POC_FLUX_IN_diag_ind),'POC_FLUX_IN',3, &
                           long_name='POC Flux into Cell',                  &
                           units='mmol/m^3 cm/s', grid_loc='3111',          &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(POC_PROD_diag_ind),'POC_PROD',3,   &
                           long_name='POC Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(POC_REMIN_diag_ind),'POC_REMIN',3, &
                           long_name='POC Remineralization',            &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(POC_REMIN_DIC_diag_ind),'POC_REMIN_DIC',3, &
                           long_name='POC Remineralization routed to DIC',      &
                           units='mmol/m^3/s', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(PON_REMIN_NH4_diag_ind),'PON_REMIN_NH4',3, &
                           long_name='PON Remineralization routed to NH4',      &
                           units='mmol/m^3/s', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(POP_REMIN_PO4_diag_ind),'POP_REMIN_PO4',3, &
                           long_name='POP Remineralization routed to PO4',      &
                           units='mmol/m^3/s', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(CaCO3_FLUX_IN_diag_ind),'CaCO3_FLUX_IN',3, &
                           long_name='CaCO3 flux into cell',                    &
                           units='mmol/m^3 cm/s', grid_loc='3111',              &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(CaCO3_PROD_diag_ind),'CaCO3_PROD',3, &
                           long_name='CaCO3 Production',                  &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(CaCO3_REMIN_diag_ind),'CaCO3_REMIN',3, &
                           long_name='CaCO3 Remineralization',              &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(SiO2_FLUX_IN_diag_ind),'SiO2_FLUX_IN',3, &
                           long_name='SiO2 Flux into Cell',                   &
                           units='mmol/m^3 cm/s', grid_loc='3111',            &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(SiO2_PROD_diag_ind),'SiO2_PROD',3, &
                           long_name='SiO2 Production',                 &
                           units='mmol/m^3/s', grid_loc='3111',         &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(SiO2_REMIN_diag_ind),'SiO2_REMIN',3, &
                           long_name='SiO2 Remineralization',             &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(dust_FLUX_IN_diag_ind),'dust_FLUX_IN',3, &
                           long_name='Dust Flux into Cell',                   &
                           units='ng/s/m^2', grid_loc='3111',                 &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(dust_REMIN_diag_ind),'dust_REMIN',3, &
                           long_name='Dust Remineralization',             &
                           units='mmol/m^3/s', grid_loc='3111',           &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(P_iron_FLUX_IN_diag_ind),'P_iron_FLUX_IN',3, &
                           long_name='P_iron Flux into Cell',                     &
                           units='mmol/m^3 cm/s', grid_loc='3111',                &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(P_iron_PROD_diag_ind),'P_iron_PROD',3, &
                           long_name='P_iron Production',                   &
                           units='mmol/m^3/s', grid_loc='3111',             &
                           coordinates='TLONG TLAT z_t time')

    call define_tavg_field(tavg_part(P_iron_REMIN_diag_ind),'P_iron_REMIN',3, &
                           long_name='P_iron Remineralization',               &
                           units='mmol/m^3/s', grid_loc='3111',               &
                           coordinates='TLONG TLAT z_t time')

    !-----------------------------------------------------------------------
    !  Vars to sum up burial in sediments and sed Denitrif N losses
    !-----------------------------------------------------------------------

    call define_tavg_field(tavg_part(calcToSed_diag_ind),'calcToSed',2, &
                           long_name='CaCO3 Flux to Sediments',         &
                           units='nmolC/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(pocToSed_diag_ind),'pocToSed',2,   &
                           long_name='POC Flux to Sediments',           &
                           units='nmolC/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(ponToSed_diag_ind),'ponToSed',2,   &
                           long_name='nitrogen burial Flux to Sediments',&
                           units='nmolN/cm^2/s', grid_loc='2110',       &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(SedDenitrif_diag_ind),'SedDenitrif',2, &
                           long_name='nitrogen loss in Sediments',          &
                           units='nmolN/cm^2/s', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(OtherRemin_diag_ind),'OtherRemin',2, &
                           long_name='non-oxic,non-dentr remin in Sediments', &
                           units='nmolC/cm^2/s', grid_loc='2110',         &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(popToSed_diag_ind),'popToSed',2, &
                           long_name='phosphorus Flux to Sediments',  &
                           units='nmolP/cm^2/s', grid_loc='2110',     &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(bsiToSed_diag_ind),'bsiToSed',2,   &
                           long_name='biogenic Si Flux to Sediments',   &
                           units='nmolSi/cm^2/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(dustToSed_diag_ind),'dustToSed',2, &
                           long_name='dust Flux to Sediments',          &
                           units='g/cm^2/s', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_part(pfeToSed_diag_ind),'pfeToSed',2,   &
                           long_name='pFe Flux to Sediments',           &
                           units='nmolFe/cm^2/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  Define tavg fields for zooplankton
!-----------------------------------------------------------------------

    do zoo_ind = 1, zooplankton_cnt
       call define_tavg_field(tavg_zoo(zoo_loss_diag_ind,zoo_ind),      &
            trim(zooplankton(zoo_ind)%sname) // '_loss', 3,             &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss',      &
            units='mmol/m^3/s', grid_loc='3114',                        &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_loss_poc_diag_ind,zoo_ind),     &
            trim(zooplankton(zoo_ind)%sname) // '_loss_poc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss to POC',  &
            units='mmol/m^3/s', grid_loc='3114',                           &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_loss_doc_diag_ind,zoo_ind),     &
            trim(zooplankton(zoo_ind)%sname) // '_loss_doc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' Loss to DOC',  &
            units='mmol/m^3/s', grid_loc='3114',                           &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_graze_diag_ind,zoo_ind),          &
            'graze_' // trim(zooplankton(zoo_ind)%sname), 3,                 &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss',   &
            units='mmol/m^3/s', grid_loc='3114',                             &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_graze_poc_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_poc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to POC', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_graze_doc_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_doc', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to DOC', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(zoo_graze_zoo_diag_ind,zoo_ind),           &
            'graze_' // trim(zooplankton(zoo_ind)%sname) // '_zoo', 3,            &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing loss to ZOO', &
            units='mmol/m^3/s', grid_loc='3114',                                  &
            coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_zoo(x_graze_zoo_diag_ind,zoo_ind),       &
            'x_graze_' // trim(zooplankton(zoo_ind)%sname), 3,              &
            long_name=trim(zooplankton(zoo_ind)%lname) // ' grazing gain',  &
            units='mmol/m^3/s', grid_loc='3114',                            &
            coordinates='TLONG TLAT z_t_150m time')
    end do


!-----------------------------------------------------------------------
!  Define tavg fields for autotrophs
!-----------------------------------------------------------------------

    do auto_ind = 1, autotroph_cnt
       call define_tavg_field(tavg_auto(N_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_N_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' N Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(P_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_P_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' P Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(Fe_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_Fe_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Fe Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(SiO3_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_SiO3_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' SiO3 Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(light_lim_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_light_lim', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Light Limitation', &
                              units='none', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(photoC_diag_ind,auto_ind), &
                              'photoC_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(photoC_zint_diag_ind,auto_ind), &
                              'photoC_' // trim(autotrophs(auto_ind)%sname) // '_zint', 2, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation Vertical Integral', &
                              units='mmol/m^3 cm/s', grid_loc='2110', &
                              coordinates='TLONG TLAT time')

       call define_tavg_field(tavg_auto(photoC_NO3_diag_ind,auto_ind), &
                              'photoC_NO3_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation from NO3', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(photoC_NO3_zint_diag_ind,auto_ind), &
                              'photoC_NO3_' // trim(autotrophs(auto_ind)%sname) // '_zint',2, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' C Fixation from NO3 Vertical Integral', &
                              units='mmol/m^3 cm/s', grid_loc='2110', &
                              coordinates='TLONG TLAT time')

       call define_tavg_field(tavg_auto(photoFe_diag_ind,auto_ind), &
                              'photoFe_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Fe Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(photoNO3_diag_ind,auto_ind), &
                              'photoNO3_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' NO3 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(photoNH4_diag_ind,auto_ind), &
                              'photoNH4_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' NH4 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(DOP_uptake_diag_ind,auto_ind), &
                              'DOP_' // trim(autotrophs(auto_ind)%sname) // '_uptake', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' DOP Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(PO4_uptake_diag_ind,auto_ind), &
                              'PO4_' // trim(autotrophs(auto_ind)%sname) // '_uptake', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' PO4 Uptake', &
                              units='mmol/m^3/s', grid_loc='3114', &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_graze_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname), 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_graze_poc_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_poc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to POC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_graze_doc_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_doc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to DOC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_graze_zoo_diag_ind,auto_ind), &
                              'graze_' // trim(autotrophs(auto_ind)%sname) // '_zoo', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Grazing to ZOO', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_loss_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_loss_poc_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss_poc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss to POC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_loss_doc_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_loss_doc', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Loss to DOC', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       call define_tavg_field(tavg_auto(auto_agg_diag_ind,auto_ind), &
                              trim(autotrophs(auto_ind)%sname) // '_agg', 3, &
                              long_name=trim(autotrophs(auto_ind)%lname) // ' Aggregate', &
                              units='mmol/m^3/s', grid_loc='3114',         &
                              coordinates='TLONG TLAT z_t_150m time')

       if (autotrophs(auto_ind)%Si_ind > 0) then
          sname = trim(autotrophs(auto_ind)%sname) // 'bSi_form'
          call define_tavg_field(tavg_auto(bSi_form_diag_ind,auto_ind), sname, 3, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' Si Uptake', &
                                 units='mmol/m^3/s', grid_loc='3114', &
                                 coordinates='TLONG TLAT z_t_150m time')
       endif

       if (autotrophs(auto_ind)%CaCO3_ind > 0) then
          sname = trim(autotrophs(auto_ind)%sname) // '_CaCO3_form'
          call define_tavg_field(tavg_auto(CaCO3_form_diag_ind,auto_ind), sname, 3, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' CaCO3 Formation', &
                                 units='mmol/m^3/s', grid_loc='3114', &
                                 coordinates='TLONG TLAT z_t_150m time')

          sname = trim(sname) // '_zint'
          call define_tavg_field(tavg_auto(CaCO3_form_zint_diag_ind,auto_ind), sname, 2, &
                                 long_name=trim(autotrophs(auto_ind)%lname) // ' CaCO3 Formation Vertical Integral', &
                                 units='mmol/m^3 cm/s', grid_loc='2110', &
                                 coordinates='TLONG TLAT time')
       endif

       if (autotrophs(auto_ind)%Nfixer) then
          call define_tavg_field(tavg_auto(Nfix_diag_ind,auto_ind), &
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

    call define_tavg_field(tavg_tot_CaCO3_form_zint, 'CaCO3_form_zint', 2, &
                           long_name='Total CaCO3 Formation Vertical Integral', &
                           units='mmol/m^3 cm/s', grid_loc='2110', &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_tot_Nfix, 'Nfix', 3, &
                           long_name='Total N Fixation', &
                           units='mmol/m^3/s', grid_loc='3114',   &
                           coordinates='TLONG TLAT z_t_150m time')

    call ecosys_restore%define_tavg_fields()

!-----------------------------------------------------------------------
!EOC

  end subroutine ecosys_tavg_init

  subroutine ecosys_tavg_write(k,bid, ecosys_diagnostics)

    integer, intent(in) :: k, bid ! level index and block index
    ! MNL: eventually this will be type(marbl_diagnostics_type) so I am leaving
    !      the variable name as marbl_diagnostics
    type(ecosys_diagnostics_type), intent(in) :: ecosys_diagnostics
    !type(ecosys_diagnostics_type), intent(in) :: marbl_diagnostics

    integer :: i, auto_ind, zoo_ind
    logical :: accumulate

   associate(                                                                          &
        DIAGS                     => ecosys_diagnostics%DIAGS,                          &
        AUTO_DIAGS                => ecosys_diagnostics%AUTO_DIAGS,                     &
        ZOO_DIAGS                 => ecosys_diagnostics%ZOO_DIAGS,                      &
        PART_DIAGS                => ecosys_diagnostics%PART_DIAGS)
    do i=1,ecosys_diag_cnt
      accumulate = .true.
      select case (i)
        case (zsatcalc_diag_ind, zsatarag_diag_ind)
          accumulate = (k.eq.km)
        case (O2_ZMIN_diag_ind, O2_ZMIN_DEPTH_diag_ind)
          accumulate = (k.eq.1)
      end select

      if (accumulate) &
        call accumulate_tavg_field(DIAGS(:,:,i),tavg_ecosys(i),bid,k)
    end do

    ! Accumulate autotroph terms
    do i=1,auto_diag_cnt
      do auto_ind=1,autotroph_cnt
        accumulate = .true.
        ! Some autotrophs are only accumulated under specific conditions
        select case (i)
          case (bSi_form_diag_ind)
            accumulate = (autotrophs(auto_ind)%Si_ind.gt.0)
            if ( accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS(:,:,auto_ind,i),  &
                                     tavg_tot_bSi_form, bid, k)
          case (CaCO3_form_diag_ind) 
            accumulate = (autotrophs(auto_ind)%imp_calcifier)
            if (accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS(:,:,auto_ind,i),  &
                                     tavg_tot_CaCO3_form, bid, k)
          case (CaCO3_form_zint_diag_ind)
            accumulate = (autotrophs(auto_ind)%imp_calcifier)
            if (accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS(:,:,auto_ind,i),  &
                                     tavg_tot_CaCO3_form_zint, bid, k)
          case (Nfix_diag_ind)
            accumulate = (autotrophs(auto_ind)%Nfixer)
            if (accumulate) &
              call accumulate_tavg_field(AUTO_DIAGS(:,:,auto_ind,i),  &
                                     tavg_tot_Nfix, bid, k)
        end select

        if (accumulate) &
          call accumulate_tavg_field(AUTO_DIAGS(:,:,auto_ind,i),      &
                                     tavg_auto(i,auto_ind), bid, k)

       end do
     end do

    ! Accumulate zooplankton terms
    do i=1,zoo_diag_cnt
      do zoo_ind=1,zooplankton_cnt
        call accumulate_tavg_field(ZOO_DIAGS(:,:,zoo_ind,i),                  &
                                   tavg_zoo(i,zoo_ind), bid, k)
      end do
    end do

    ! Accumulate particulate terms
    do i=1,part_diag_cnt
      call accumulate_tavg_field(PART_DIAGS(:,:,i), tavg_part(i), bid, k)
    end do

    end associate

  end subroutine ecosys_tavg_write

  subroutine ecosys_tavg_write_flux(FLUX_DIAGS)

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

  end subroutine ecosys_tavg_write_flux

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


