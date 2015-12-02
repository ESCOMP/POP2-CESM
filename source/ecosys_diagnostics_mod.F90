! Will be part of MARBL library
module ecosys_diagnostics_mod

  use grid, only : partial_bottom_cells
  use domain_size, only : km

  use marbl_kinds_mod, only : c0
  use marbl_kinds_mod, only : c1
  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : log_kind
  use marbl_kinds_mod, only : char_len

  use marbl_share_mod, only : autotrophs
  use marbl_share_mod, only : autotroph_cnt
  use marbl_share_mod, only : zooplankton
  use marbl_share_mod, only : zooplankton_cnt
  use marbl_share_mod, only : ecosys_tracer_cnt

  use marbl_interface_types, only : max_interior_diags
  use marbl_interface_types, only : max_restore_diags

  use marbl_interface_types, only : carbonate_type
  use marbl_interface_types, only : zooplankton_secondary_species_type
  use marbl_interface_types, only : autotroph_secondary_species_type
  use marbl_interface_types, only : photosynthetically_available_radiation_type
  use marbl_interface_types, only : dissolved_organic_matter_type
  use marbl_interface_types, only : marbl_diagnostics_type
  use marbl_interface_types, only : marbl_column_domain_type
  use marbl_interface_types, only : marbl_gcm_state_type

  use marbl_parms, only : po4_ind 
  use marbl_parms, only : no3_ind         
  use marbl_parms, only : sio3_ind        
  use marbl_parms, only : nh4_ind         
  use marbl_parms, only : dic_ind         
  use marbl_parms, only : doc_ind         
  use marbl_parms, only : don_ind         
  use marbl_parms, only : dop_ind         
  use marbl_parms, only : dopr_ind        
  use marbl_parms, only : donr_ind        

  Implicit None
  Public

  !-----------------------------------------------------------------------
  !  indices for diagnostic values written to tavg files
  !-----------------------------------------------------------------------

  type, private :: marbl_interior_diagnostics_indexing_type
    ! General 2D diags
    integer(int_kind) :: zsatcalc
    integer(int_kind) :: zsatarag
    integer(int_kind) :: O2_ZMIN
    integer(int_kind) :: O2_ZMIN_DEPTH
    integer(int_kind) :: photoC_TOT_zint
    integer(int_kind) :: photoC_NO3_TOT_zint
    integer(int_kind) :: Jint_Ctot
    integer(int_kind) :: Jint_100m_Ctot
    integer(int_kind) :: Jint_Ntot
    integer(int_kind) :: Jint_100m_Ntot
    integer(int_kind) :: Jint_Ptot
    integer(int_kind) :: Jint_100m_Ptot
    integer(int_kind) :: Jint_Sitot
    integer(int_kind) :: Jint_100m_Sitot

    ! Particulate 2D diags
    integer(int_kind) :: calcToSed
    integer(int_kind) :: pocToSed
    integer(int_kind) :: ponToSed
    integer(int_kind) :: SedDenitrif
    integer(int_kind) :: OtherRemin
    integer(int_kind) :: popToSed
    integer(int_kind) :: bsiToSed
    integer(int_kind) :: dustToSed
    integer(int_kind) :: pfeToSed

    ! Autotroph 2D diags
    integer(int_kind), dimension(autotroph_cnt) :: photoC_zint
    integer(int_kind), dimension(autotroph_cnt) :: photoC_NO3_zint
    integer(int_kind), dimension(autotroph_cnt) :: CaCO3_form_zint
    integer(int_kind) :: tot_CaCO3_form_zint

    ! General 3D diags
    integer(int_kind) :: CO3
    integer(int_kind) :: HCO3
    integer(int_kind) :: H2CO3
    integer(int_kind) :: pH_3D
    integer(int_kind) :: CO3_ALT_CO2
    integer(int_kind) :: HCO3_ALT_CO2
    integer(int_kind) :: H2CO3_ALT_CO2
    integer(int_kind) :: pH_3D_ALT_CO2
    integer(int_kind) :: co3_sat_calc
    integer(int_kind) :: co3_sat_arag
    integer(int_kind) :: NITRIF
    integer(int_kind) :: DENITRIF
    integer(int_kind) :: O2_PRODUCTION
    integer(int_kind) :: O2_CONSUMPTION
    integer(int_kind) :: AOU
    integer(int_kind) :: PAR_avg
    integer(int_kind) :: auto_graze_TOT
    integer(int_kind) :: photoC_TOT
    integer(int_kind) :: photoC_NO3_TOT
    integer(int_kind) :: DOC_prod
    integer(int_kind) :: DOC_remin
    integer(int_kind) :: DON_prod
    integer(int_kind) :: DON_remin
    integer(int_kind) :: DOP_prod
    integer(int_kind) :: DOP_remin
    integer(int_kind) :: DOFe_prod
    integer(int_kind) :: DOFe_remin
    integer(int_kind) :: DONr_remin
    integer(int_kind) :: DOPr_remin
    integer(int_kind) :: Fe_scavenge
    integer(int_kind) :: Fe_scavenge_rate

    ! Particulate 3D diags
    integer(int_kind) :: POC_FLUX_IN
    integer(int_kind) :: POC_PROD
    integer(int_kind) :: POC_REMIN
    integer(int_kind) :: CaCO3_FLUX_IN
    integer(int_kind) :: CaCO3_PROD
    integer(int_kind) :: CaCO3_REMIN
    integer(int_kind) :: SiO2_FLUX_IN
    integer(int_kind) :: SiO2_PROD
    integer(int_kind) :: SiO2_REMIN
    integer(int_kind) :: dust_FLUX_IN
    integer(int_kind) :: dust_REMIN
    integer(int_kind) :: P_iron_FLUX_IN
    integer(int_kind) :: P_iron_PROD
    integer(int_kind) :: P_iron_REMIN

    ! Autotroph 3D diags
    integer(int_kind), dimension(autotroph_cnt) :: N_lim
    integer(int_kind), dimension(autotroph_cnt) :: P_lim
    integer(int_kind), dimension(autotroph_cnt) :: Fe_lim
    integer(int_kind), dimension(autotroph_cnt) :: SiO3_lim
    integer(int_kind), dimension(autotroph_cnt) :: light_lim
    integer(int_kind), dimension(autotroph_cnt) :: photoC
    integer(int_kind), dimension(autotroph_cnt) :: photoC_NO3
    integer(int_kind), dimension(autotroph_cnt) :: photoFe
    integer(int_kind), dimension(autotroph_cnt) :: photoNO3
    integer(int_kind), dimension(autotroph_cnt) :: photoNH4
    integer(int_kind), dimension(autotroph_cnt) :: DOP_uptake
    integer(int_kind), dimension(autotroph_cnt) :: PO4_uptake
    integer(int_kind), dimension(autotroph_cnt) :: auto_graze
    integer(int_kind), dimension(autotroph_cnt) :: auto_graze_poc
    integer(int_kind), dimension(autotroph_cnt) :: auto_graze_doc
    integer(int_kind), dimension(autotroph_cnt) :: auto_graze_zoo
    integer(int_kind), dimension(autotroph_cnt) :: auto_loss
    integer(int_kind), dimension(autotroph_cnt) :: auto_loss_poc
    integer(int_kind), dimension(autotroph_cnt) :: auto_loss_doc
    integer(int_kind), dimension(autotroph_cnt) :: auto_agg
    integer(int_kind), dimension(autotroph_cnt) :: bSi_form
    integer(int_kind), dimension(autotroph_cnt) :: CaCO3_form
    integer(int_kind), dimension(autotroph_cnt) :: Nfix
    integer(int_kind) :: tot_bSi_form
    integer(int_kind) :: tot_CaCO3_form
    integer(int_kind) :: tot_Nfix

    ! zooplankton 3D diags
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_loss
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_loss_poc
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_loss_doc
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_graze
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_graze_poc
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_graze_doc
    integer(int_kind), dimension(zooplankton_cnt) :: zoo_graze_zoo
    integer(int_kind), dimension(zooplankton_cnt) :: x_graze_zoo
  end type marbl_interior_diagnostics_indexing_type
  type(marbl_interior_diagnostics_indexing_type), public :: marbl_interior_diag_ind

  type marbl_forcing_diagnostics_indexing_type
     integer(int_kind) :: ECOSYS_IFRAC
     integer(int_kind) :: ECOSYS_IFRAC_2   
     integer(int_kind) :: ECOSYS_XKW
     integer(int_kind) :: ECOSYS_XKW_2     
     integer(int_kind) :: ECOSYS_ATM_PRESS
     integer(int_kind) :: PV_O2
     integer(int_kind) :: SCHMIDT_O2
     integer(int_kind) :: O2SAT
     integer(int_kind) :: O2_GAS_FLUX
     integer(int_kind) :: O2_GAS_FLUX_2    
     integer(int_kind) :: CO2STAR
     integer(int_kind) :: DCO2STAR
     integer(int_kind) :: pCO2SURF
     integer(int_kind) :: DpCO2
     integer(int_kind) :: DpCO2_2          
     integer(int_kind) :: PV_CO2
     integer(int_kind) :: SCHMIDT_CO2
     integer(int_kind) :: DIC_GAS_FLUX
     integer(int_kind) :: DIC_GAS_FLUX_2   
     integer(int_kind) :: PH
     integer(int_kind) :: ATM_CO2
     integer(int_kind) :: CO2STAR_ALT_CO2
     integer(int_kind) :: DCO2STAR_ALT_CO2
     integer(int_kind) :: pCO2SURF_ALT_CO2
     integer(int_kind) :: DpCO2_ALT_CO2
     integer(int_kind) :: DIC_GAS_FLUX_ALT_CO2
     integer(int_kind) :: PH_ALT_CO2
     integer(int_kind) :: ATM_ALT_CO2
     integer(int_kind) :: IRON_FLUX
     integer(int_kind) :: DUST_FLUX
     integer(int_kind) :: NOx_FLUX
     integer(int_kind) :: NHy_FLUX
     integer(int_kind) :: DIN_RIV_FLUX
     integer(int_kind) :: DIP_RIV_FLUX
     integer(int_kind) :: DoN_RIV_FLUX
     integer(int_kind) :: DoNr_RIV_FLUX
     integer(int_kind) :: DOP_RIV_FLUX
     integer(int_kind) :: DOPr_RIV_FLUX
     integer(int_kind) :: DSI_RIV_FLUX
     integer(int_kind) :: DFE_RIV_FLUX
     integer(int_kind) :: DIC_RIV_FLUX
     integer(int_kind) :: ALK_RIV_FLUX
     integer(int_kind) :: DOC_RIV_FLUX
  end type marbl_forcing_diagnostics_indexing_type
  type(marbl_forcing_diagnostics_indexing_type), public :: marbl_forcing_diag_ind

  !***********************************************************************

contains

  !***********************************************************************

  subroutine marbl_ecosys_diagnostics_init(marbl_interior_diags, &
       marbl_restore_diags, marbl_forcing_diags, tracer_d_module)

    use prognostic           , only : tracer_field
    use marbl_interface_types, only : max_forcing_diags
    use marbl_interface_types, only : max_interior_diags
    use marbl_interface_types, only : max_restore_diags

    type(marbl_diagnostics_type), intent(inout) :: marbl_interior_diags
    type(marbl_diagnostics_type), intent(inout) :: marbl_restore_diags
    type(marbl_diagnostics_type), intent(inout) :: marbl_forcing_diags
    type (tracer_field)         , intent(in)    :: tracer_d_module(:)   ! descriptors for each tracer

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer :: n, tmp_id
    character(len=char_len) :: lname, sname, units, vgrid
    logical :: truncate
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Surface forcing diagnostics
    !-----------------------------------------------------------------

    ! Allocate memory for forcing diagnostics
    call marbl_forcing_diags%construct(max_forcing_diags)

    lname='Ice Fraction for ecosys fluxes'
    sname='ECOSYS_IFRAC'
    units='fraction'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ECOSYS_IFRAC)

    lname='Ice Fraction for ecosys fluxes'
    sname='ECOSYS_IFRAC_2'
    units='fraction'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ECOSYS_IFRAC_2)

    lname='XKW for ecosys fluxes'
    sname='ECOSYS_XKW'
    units='cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ECOSYS_XKW)

    lname='XKW for ecosys fluxes'
    sname='ECOSYS_XKW_2'
    units='cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ECOSYS_XKW_2)

    lname='Atmospheric Pressure for ecosys fluxes'
    sname='ECOSYS_ATM_PRESS'
    units='atmospheres'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ECOSYS_ATM_PRESS)

    lname='PV_O2'
    sname='PV_O2'
    units='cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%PV_O2)

    lname='O2 Schmidt Number'
    sname='SCHMIDT_O2'
    units='none'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%SCHMIDT_O2)

    lname='O2 Saturation'
    sname='O2SAT'
    units='mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%O2SAT)

    lname='Dissolved Oxygen Surface Flux'
    sname='STF_O2'
    units='mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%O2_GAS_FLUX)

    lname='Dissolved Oxygen Surface Flux'
    sname='STF_O2_2'
    units='mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%O2_GAS_FLUX_2)

    lname='CO2 Star'
    sname='CO2STAR'
    units='mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%CO2STAR)

    lname='D CO2 Star'
    sname='DCO2STAR'
    units='mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DCO2STAR)

    lname='surface pCO2'
    sname='pCO2SURF'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%pCO2SURF)

    lname='D pCO2'
    sname='DpCO2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DpCO2)

    lname='D pCO2'
    sname='DpCO2_2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DpCO2_2)

    lname='CO2 Piston Velocity'
    sname='PV_CO2'
    units='cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%PV_CO2)

    lname='CO2 Schmidt Number'
    sname='SCHMIDT_CO2'
    units='none'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%SCHMIDT_CO2)

    lname='DIC Surface Gas Flux'
    sname='FG_CO2'
    units='mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIC_GAS_FLUX)

    lname='DIC Surface Gas Flux'
    sname='FG_CO2_2'
    units='mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIC_GAS_FLUX_2)

    lname='Surface pH'
    sname='PH'
    units='none'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%PH)

    lname='Atmospheric CO2'
    sname='ATM_CO2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ATM_CO2)

    lname='CO2 Star, Alternative CO2'
    sname='CO2STAR_ALT_CO2'
    units='mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%CO2STAR_ALT_CO2)

    lname='D CO2 Star, Alternative CO2'
    sname='DCO2STAR_ALT_CO2'
    units='mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DCO2STAR_ALT_CO2)

    lname='surface pCO2, Alternative CO2'
    sname='pCO2SURF_ALT_CO2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%pCO2SURF_ALT_CO2)

    lname='D pCO2, Alternative CO2'
    sname='DpCO2_ALT_CO2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DpCO2_ALT_CO2)

    lname='DIC Surface Gas Flux, Alternative CO2'
    sname='FG_ALT_CO2'
    units='mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIC_GAS_FLUX_ALT_CO2)

    lname='Surface pH, Alternative CO2'
    sname='PH_ALT_CO2'
    units='none'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%PH_ALT_CO2)

    lname='Atmospheric Alternative CO2'
    sname='ATM_ALT_CO2'
    units='ppmv'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ATM_ALT_CO2)

    lname='Atmospheric Iron Flux'
    sname='IRON_FLUX'
    units='mmol/m^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%IRON_FLUX)

    lname='Dust Flux'
    sname='DUST_FLUX'
    units='g/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DUST_FLUX)

    lname='Flux of NOx from Atmosphere'
    sname='NOx_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%NOx_FLUX)

    lname='Flux of NHy from Atmosphere'
    sname='NHy_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%NHy_FLUX)

    lname='Flux of DIN from rivers'
    sname='DIN_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIN_RIV_FLUX)

    lname='Flux of DIP from rivers'
    sname='DIP_RIV_FLUX'
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIP_RIV_FLUX)
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.

    lname='Flux of DON from rivers'
    sname='DON_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DON_RIV_FLUX)

    lname='Flux of DONr from rivers'
    sname='DONr_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DONr_RIV_FLUX)

    lname='Flux of DOP from rivers'
    sname='DOP_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DOP_RIV_FLUX)

    lname='Flux of DOPr from rivers'
    sname='DOPr_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DOPr_RIV_FLUX)

    lname='Flux of DSI from rivers'
    sname='DSI_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DSI_RIV_FLUX)

    lname='Flux of DFE from rivers'
    sname='DFE_RIV_FLUX'
    units='nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DFE_RIV_FLUX)

    lname = 'Flux of DIC from rivers'
    sname = 'DIC_RIV_FLUX'
    units = 'nmol/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DIC_RIV_FLUX)

    lname='Flux of ALK from rivers'
    sname='ALK_RIV_FLUX'
    units='alk/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%ALK_RIV_FLUX)

    lname='Flux of DOC from rivers'
    sname='DOC_RIV_FLUX'
    units='alk/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_forcing_diags%add_diagnostic(&
         lname, sname, units, vgrid, truncate, marbl_forcing_diag_ind%DOC_RIV_FLUX)

    !-----------------------------------------------------------------
    ! Interior diagnostics
    !-----------------------------------------------------------------

    ! Allocate memory for interior diagnostics
    call marbl_interior_diags%construct(max_interior_diags)

    ! General 2D diags
    lname = 'Calcite Saturation Depth'
    sname = 'zsatcalc'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%zsatcalc)

    lname = 'Aragonite Saturation Depth'
    sname = 'zsatarag'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%zsatarag)

    lname = 'Vertical Minimum of O2'
    sname = 'O2_ZMIN'
    units = 'mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%O2_ZMIN)

    lname = 'Depth of Vertical Minimum of O2'
    sname = 'O2_ZMIN_DEPTH'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%O2_ZMIN_DEPTH)

    lname = 'Total C Fixation Vertical Integral'
    sname = 'photoC_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%photoC_TOT_zint)

    lname = 'Total C Fixation from NO3 Vertical Integral'
    sname = 'photoC_NO3_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%photoC_NO3_TOT_zint)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot'
    sname = 'Jint_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot, 0-100m'
    sname = 'Jint_100m_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_100m_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot'
    sname = 'Jint_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot, 0-100m'
    sname = 'Jint_100m_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_100m_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot'
    sname = 'Jint_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot, 0-100m'
    sname = 'Jint_100m_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_100m_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot'
    sname = 'Jint_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_Sitot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot, 0-100m'
    sname = 'Jint_100m_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Jint_100m_Sitot)

    ! Particulate 2D diags
    lname = 'CaCO3 Flux to Sediments'
    sname = 'calcToSed'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%calcToSed)

    lname = 'POC Flux to Sediments'
    sname = 'pocToSed'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%pocToSed)

    lname = 'nitrogen burial Flux to Sediments'
    sname = 'ponToSed'
    units = 'nmolN/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%ponToSed)

    lname = 'nitrogen loss in Sediments'
    sname = 'SedDenitrif'
    units = 'nmolN/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%SedDenitrif)

    lname = 'non-oxic,non-dentr remin in Sediments'
    sname = 'OtherRemin'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%OtherRemin)

    lname = 'phosphorus Flux to Sediments'
    sname = 'popToSed'
    units = 'nmolP/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%popToSed)

    lname = 'biogenic Si Flux to Sediments'
    sname = 'bsiToSed'
    units = 'nmolSi/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%bsiToSed)

    lname = 'dust Flux to Sediments'
    sname = 'dustToSed'
    units = 'g/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%dustToSed)

    lname = 'pFe Flux to Sediments'
    sname = 'pfeToSed'
    units = 'nmolFe/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%pfeToSed)

    ! Autotroph 2D diags
    do n=1,autotroph_cnt
      lname = trim(autotrophs(n)%lname) // ' C Fixation Vertical Integral'
      sname = 'photoC_' // trim(autotrophs(n)%sname) // '_zint'
      units = 'mmol/m^3 cm/s'
      vgrid = 'none'
      truncate = .false.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoC_zint(n))

      lname = trim(autotrophs(n)%lname) // ' C Fixation from NO3 Vertical Integral'
      sname = 'photoC_NO3_' // trim(autotrophs(n)%sname) // '_zint'
      units = 'mmol/m^3 cm/s'
      vgrid = 'none'
      truncate = .false.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoC_NO3_zint(n))

      if (autotrophs(n)%CaCO3_ind.gt.0) then
        lname = trim(autotrophs(n)%lname) //                                  &
                ' CaCO3 Formation Vertical Integral'
        sname = trim(autotrophs(n)%sname) // '_CaCO3_form_zint'
        units = 'mmol/m^3 cm/s'
        vgrid = 'none'
        truncate = .false.
        call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                        marbl_interior_diag_ind%CaCO3_form_zint(n))
      else
        marbl_interior_diag_ind%CaCO3_form_zint(n) = -1
      end if
    end do

    lname = 'Total CaCO3 Formation Vertical Integral'
    sname = 'CaCO3_form_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                      marbl_interior_diag_ind%tot_CaCO3_form_zint)

    ! General 3D diags
    lname = 'Carbonate Ion Concentration'
    sname = 'CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%CO3)

    lname = 'Bicarbonate Ion Concentration'
    sname = 'HCO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%HCO3)

    lname = 'Carbonic Acid Concentration'
    sname = 'H2CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%H2CO3)

    lname = 'pH'
    sname = 'pH_3D'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%ph_3D)

    lname = 'Carbonate Ion Concentration, Alternative CO2'
    sname = 'CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%CO3_ALT_CO2)

    lname = 'Bicarbonate Ion Concentration, Alternative CO2'
    sname = 'HCO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%HCO3_ALT_CO2)

    lname = 'Carbonic Acid Concentration, Alternative CO2'
    sname = 'H2CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%H2CO3_ALT_CO2)

    lname = 'pH, Alternative CO2'
    sname = 'pH_3D_ALT_CO2'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%ph_3D_ALT_CO2)

    lname = 'CO3 concentration at calcite saturation'
    sname = 'co3_sat_calc'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%co3_sat_calc)

    lname = 'CO3 concentration at aragonite saturation'
    sname = 'co3_sat_arag'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%co3_sat_arag)

    lname = 'Nitrification'
    sname = 'NITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%NITRIF)

    lname = 'Denitrification'
    sname = 'DENITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DENITRIF)

    lname = 'O2 Production'
    sname = 'O2_PRODUCTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%O2_PRODUCTION)

    lname = 'O2 Consumption'
    sname = 'O2_CONSUMPTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%O2_CONSUMPTION)

    lname = 'Apparent O2 Utilization'
    sname = 'AOU'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%AOU)

    lname = 'PAR Average over Model Cell'
    sname = 'PAR_avg'
    units = 'W/m^2'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%PAR_avg)

    lname = 'Total Autotroph Grazing'
    sname = 'graze_auto_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%auto_graze_TOT)

    lname = 'Total C Fixation'
    sname = 'photoC_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%photoC_TOT)

    lname = 'Total C Fixation from NO3'
    sname = 'photoC_NO3_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%photoC_NO3_TOT)

    lname = 'DOC Production'
    sname = 'DOC_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOC_prod)

    lname = 'DOC Remineralization'
    sname = 'DOC_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOC_remin)

    lname = 'DON Production'
    sname = 'DON_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DON_prod)

    lname = 'DON Remineralization'
    sname = 'DON_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DON_remin)

    lname = 'DOP Production'
    sname = 'DOP_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOP_prod)

    lname = 'DOP Remineralization'
    sname = 'DOP_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOP_remin)

    lname = 'DOFe Production'
    sname = 'DOFe_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOFe_prod)

    lname = 'DOFe Remineralization'
    sname = 'DOFe_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOFe_remin)

    lname = 'DONr Remineralization'
    sname = 'DONr_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DONr_remin)

    lname = 'DOPr Remineralization'
    sname = 'DOPr_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%DOPr_remin)

    lname = 'Iron Scavenging'
    sname = 'Fe_scavenge'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Fe_scavenge)

    lname = 'Iron Scavenging Rate'
    sname = 'Fe_scavenge_rate'
    units = '1/y'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%Fe_scavenge_rate)

    ! Particulate 3D diags
    lname = 'POC Flux into Cell'
    sname = 'POC_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%POC_FLUX_IN)

    lname = 'POC Production'
    sname = 'POC_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%POC_PROD)

    lname = 'POC Remineralization'
    sname = 'POC_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%POC_REMIN)

    lname = 'CaCO3 Flux into Cell'
    sname = 'CaCO3_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%CaCO3_FLUX_IN)

    lname = 'CaCO3 Production'
    sname = 'CaCO3_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%CaCO3_PROD)

    lname = 'CaCO3 Remineralization'
    sname = 'CaCO3_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%CaCO3_REMIN)

    lname = 'SiO2 Flux into Cell'
    sname = 'SiO2_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%SiO2_FLUX_IN)

    lname = 'SiO2 Production'
    sname = 'SiO2_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%SiO2_PROD)

    lname = 'SiO2 Remineralization'
    sname = 'SiO2_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%SiO2_REMIN)

    lname = 'Dust Flux into Cell'
    sname = 'dust_FLUX_IN'
    units = 'ng/s/m^2'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%dust_FLUX_IN)

    lname = 'Dust Remineralization'
    sname = 'dust_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%dust_REMIN)

    lname = 'P_iron Flux into Cell'
    sname = 'P_iron_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%P_iron_FLUX_IN)

    lname = 'P_iron Production'
    sname = 'P_iron_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%P_iron_PROD)

    lname = 'P_iron Remineralization'
    sname = 'P_iron_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_interior_diag_ind%P_iron_REMIN)

    ! Autotroph 3D diags
    do n=1,autotroph_cnt
      lname = trim(autotrophs(n)%lname) // ' N Limitation'
      sname = trim(autotrophs(n)%sname) // '_N_lim'
      units = 'none'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%N_lim(n))

      lname = trim(autotrophs(n)%lname) // ' P Limitation'
      sname = trim(autotrophs(n)%sname) // '_P_lim'
      units = 'none'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%P_lim(n))

      lname = trim(autotrophs(n)%lname) // ' Fe Limitation'
      sname = trim(autotrophs(n)%sname) // '_Fe_lim'
      units = 'none'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%Fe_lim(n))

      if (autotrophs(n)%kSiO3 > c0) then
        lname = trim(autotrophs(n)%lname) // ' SiO3 Limitation'
        sname = trim(autotrophs(n)%sname) // '_SiO3_lim'
        units = 'none'
        vgrid = 'layer_avg'
        truncate = .true.
        call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                        marbl_interior_diag_ind%SiO3_lim(n))
      else
        marbl_interior_diag_ind%SiO3_lim(n) = -1
      end if

      lname = trim(autotrophs(n)%lname) // ' Light Limitation'
      sname = trim(autotrophs(n)%sname) // '_light_lim'
      units = 'none'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%light_lim(n))

      lname = trim(autotrophs(n)%lname) // ' C Fixation'
      sname = 'photoC_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoC(n))

      lname = trim(autotrophs(n)%lname) // ' C Fixation from NO3'
      sname = 'photoC_NO3_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoC_NO3(n))

      lname = trim(autotrophs(n)%lname) // ' Fe Uptake'
      sname = 'photoFe_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoFe(n))

      lname = trim(autotrophs(n)%lname) // ' NO3 Uptake'
      sname = 'photoNO3_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%photoNO3(n))

      lname = trim(autotrophs(n)%lname) // ' NH4 Uptake'
      sname = 'photoNH4_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%photoNH4(n))

      lname = trim(autotrophs(n)%lname) // ' DOP Uptake'
      sname = 'DOP_' // trim(autotrophs(n)%sname) // '_uptake'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%DOP_uptake(n))

      lname = trim(autotrophs(n)%lname) // ' PO4 Uptake'
      sname = 'PO4_' // trim(autotrophs(n)%sname) // '_uptake'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%PO4_uptake(n))

      lname = trim(autotrophs(n)%lname) // ' Grazing'
      sname = 'graze_' // trim(autotrophs(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_graze(n))

      lname = trim(autotrophs(n)%lname) // ' Grazing to POC'
      sname = 'graze_' // trim(autotrophs(n)%sname) // '_poc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_graze_poc(n))

      lname = trim(autotrophs(n)%lname) // ' Grazing to DOC'
      sname = 'graze_' // trim(autotrophs(n)%sname) // '_doc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_graze_doc(n))

      lname = trim(autotrophs(n)%lname) // ' Grazing to ZOO'
      sname = 'graze_' // trim(autotrophs(n)%sname) // '_zoo'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_graze_zoo(n))

      lname = trim(autotrophs(n)%lname) // ' Loss'
      sname = trim(autotrophs(n)%sname) // '_loss'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_loss(n))

      lname = trim(autotrophs(n)%lname) // ' Loss to POC'
      sname = trim(autotrophs(n)%sname) // '_loss_poc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_loss_poc(n))

      lname = trim(autotrophs(n)%lname) // ' Loss to DOC'
      sname = trim(autotrophs(n)%sname) // '_loss_doc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_loss_doc(n))

      lname = trim(autotrophs(n)%lname) // ' Aggregate'
      sname = trim(autotrophs(n)%sname) // '_agg'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
           marbl_interior_diag_ind%auto_agg(n))

      if (autotrophs(n)%Si_ind.gt.0) then
        lname=trim(autotrophs(n)%lname) // ' Si Uptake'
        !sname = trim(autotrophs(n)%sname) // '_bSi_form'
        sname = trim(autotrophs(n)%sname) // 'bSi_form' ! FIXME: eventually add _
        units = 'mmol/m^3/s'
        vgrid = 'layer_avg'
        truncate = .true.
        call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
             marbl_interior_diag_ind%bSi_form(n))
      else
        marbl_interior_diag_ind%bSi_form(n) = -1
      end if

      if (autotrophs(n)%CaCO3_ind.gt.0) then
        lname=trim(autotrophs(n)%lname) // ' CaCO3 Formation'
        sname = trim(autotrophs(n)%sname) // '_CaCO3_form'
        units = 'mmol/m^3/s'
        vgrid = 'layer_avg'
        truncate = .true.
        call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
             marbl_interior_diag_ind%CaCO3_form(n))
      else
        marbl_interior_diag_ind%CaCO3_form(n) = -1
      end if

      if (autotrophs(n)%Nfixer) then
        lname=trim(autotrophs(n)%lname) // ' N Fixation'
        sname = trim(autotrophs(n)%sname) // '_Nfix'
        units = 'mmol/m^3/s'
        vgrid = 'layer_avg'
        truncate = .true.
        call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                        marbl_interior_diag_ind%Nfix(n))
      else
        marbl_interior_diag_ind%Nfix(n) = -1
      end if
    end do

    lname= 'Total Si Uptake'
    sname = 'bSi_form'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                             marbl_interior_diag_ind%tot_bSi_form)

    lname = 'Total CaCO3 Formation'
    sname = 'CaCO3_form'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                           marbl_interior_diag_ind%tot_CaCO3_form)

    lname = 'Total N Fixation'
    sname = 'Nfix'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate, &
                                        marbl_interior_diag_ind%tot_Nfix)

    ! Zooplankton 3D diags
    do n=1,zooplankton_cnt
      lname = trim(zooplankton(n)%lname) // ' Loss'
      sname = trim(zooplankton(n)%sname) // '_loss'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_loss(n))

      lname = trim(zooplankton(n)%lname) // ' Loss to POC'
      sname = trim(zooplankton(n)%sname) // '_loss_poc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_loss_poc(n))

      lname = trim(zooplankton(n)%lname) // ' Loss to DOC'
      sname = trim(zooplankton(n)%sname) // '_loss_doc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_loss_doc(n))

      lname = trim(zooplankton(n)%lname) // ' grazing loss'
      sname = 'graze_' // trim(zooplankton(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_graze(n))

      lname = trim(zooplankton(n)%lname) // ' grazing loss to POC'
      sname = 'graze_' // trim(zooplankton(n)%sname) // '_poc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_graze_poc(n))

      lname = trim(zooplankton(n)%lname) // ' grazing loss to DOC'
      sname = 'graze_' // trim(zooplankton(n)%sname) // '_doc'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_graze_doc(n))

      lname = trim(zooplankton(n)%lname) // ' grazing loss to ZOO'
      sname = 'graze_' // trim(zooplankton(n)%sname) // '_zoo'
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%zoo_graze_zoo(n))

      lname = trim(zooplankton(n)%lname) // ' grazing gain'
      sname = 'x_graze_' // trim(zooplankton(n)%sname)
      units = 'mmol/m^3/s'
      vgrid = 'layer_avg'
      truncate = .true.
      call marbl_interior_diags%add_diagnostic(lname, sname, units, vgrid, truncate,   &
                                      marbl_interior_diag_ind%x_graze_zoo(n))

    end do

    !-----------------------------------------------------------------
    ! Restoring diagnostics
    !-----------------------------------------------------------------

    ! Allocate memory for restore diagnostics
    call marbl_restore_diags%construct(max_restore_diags)

    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    do n=1,ecosys_tracer_cnt
      lname = trim(tracer_d_module(n)%long_name) // " Restoring"
      sname = trim(tracer_d_module(n)%short_name) // "_RESTORE"

      ! Note that tmp_id is a temp variable because restoring diagnostics
      ! have same indexing as the ecosys tracers
      call marbl_restore_diags%add_diagnostic(lname, sname, units, vgrid, .false., tmp_id)
    end do

  end subroutine marbl_ecosys_diagnostics_init

  !***********************************************************************

  subroutine store_diagnostics_carbonate(marbl_domain, carbonate, zt, zw, marbl_diags)

    type(carbonate_type), dimension(:), intent(in) :: carbonate
    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zt, zw
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    real(r8) :: zsat_calcite, zsat_aragonite

    associate(                                                                &
              diags => marbl_diags%diags(:),                                  &
              CO3 => carbonate(:)%CO3,                                        &
              CO3_sat_calcite => carbonate(:)%CO3_sat_calcite,                &
              CO3_sat_aragonite => carbonate(:)%CO3_sat_aragonite             &
              )
      ! Find depth where CO3 = CO3_sat_calcite or CO3_sat_argonite
      diags(marbl_interior_diag_ind%zsatcalc)%field_2d =                               &
           marbl_compute_saturation_depth(marbl_domain, zt, zw, CO3, CO3_sat_calcite)
      diags(marbl_interior_diag_ind%zsatarag)%field_2d =                               &
           marbl_compute_saturation_depth(marbl_domain, zt, zw, CO3, CO3_sat_aragonite)

      diags(marbl_interior_diag_ind%CO3)%field_3d(:)           = carbonate%CO3
      diags(marbl_interior_diag_ind%HCO3)%field_3d(:)          = carbonate%HCO3
      diags(marbl_interior_diag_ind%H2CO3)%field_3d(:)         = carbonate%H2CO3
      diags(marbl_interior_diag_ind%pH_3D)%field_3d(:)         = carbonate%pH
      diags(marbl_interior_diag_ind%CO3_ALT_CO2)%field_3d(:)   = carbonate%CO3_ALT_CO2
      diags(marbl_interior_diag_ind%HCO3_ALT_CO2)%field_3d(:)  = carbonate%HCO3_ALT_CO2
      diags(marbl_interior_diag_ind%H2CO3_ALT_CO2)%field_3d(:) = carbonate%H2CO3_ALT_CO2
      diags(marbl_interior_diag_ind%pH_3D_ALT_CO2)%field_3d(:) = carbonate%pH_ALT_CO2
      diags(marbl_interior_diag_ind%co3_sat_calc)%field_3d(:)  = carbonate%CO3_sat_calcite
      diags(marbl_interior_diag_ind%co3_sat_arag)%field_3d(:)  = carbonate%CO3_sat_aragonite

    end associate

  end subroutine store_diagnostics_carbonate

  !***********************************************************************

  function marbl_compute_saturation_depth(marbl_domain, zt, zw, CO3, sat_val)

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zt, zw, CO3, sat_val
    real(r8) :: marbl_compute_saturation_depth

    real(r8) :: anomaly(km) ! CO3 concentration above saturation at level
    integer :: k

    associate(column_kmt => marbl_domain%kmt)
      anomaly   = CO3 - sat_val

      if (all(anomaly(1:column_kmt).gt.c0)) then
         marbl_compute_saturation_depth = zw(column_kmt)
      elseif (anomaly(1).le.c0) then
         marbl_compute_saturation_depth = c0
      else
        do k=2,column_kmt
          if (anomaly(k).le.c0) exit
        end do
        ! saturation depth is location of root of anomaly
        marbl_compute_saturation_depth = linear_root(zt(k-1:k), anomaly(k-1:k))
      end if
    end associate

  end function marbl_compute_saturation_depth

  !***********************************************************************

  function linear_root(x,y)
  ! TO-DO (MNL): if we end up with a marbl_math_mod, this can be generalized
  !              to a better root-finding routine; otherwise maybe we compute
  !              the root inside compute_saturation_depth rather than as a
  !              separate function?

    real(kind=r8), dimension(2), intent(in) :: x,y
    real(kind=r8) :: linear_root

    real(kind=r8) :: m_inv

    if (y(1)*y(2).gt.c0) then
      ! TO-DO (MNL): do we have a marbl_abort() routine? How do I exit if we
      !              hit this error?
      print*, "MNL MNL MNL: can not find root, y-values are same sign!"
    end if
    if (y(2).eq.c0) then
      linear_root = x(2)
    else
      m_inv = (x(2)-x(1))/(y(2)-y(1))
      linear_root = x(1)-m_inv*y(1)
    end if

  end function linear_root

  !***********************************************************************

  subroutine store_diagnostics_nitrification(nitrif, denitrif, marbl_diags)

    real(r8), dimension(:), intent(in) :: nitrif
    real(r8), dimension(:), intent(in) :: denitrif
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_interior_diag_ind%NITRIF)%field_3d(:) = nitrif
      diags(marbl_interior_diag_ind%DENITRIF)%field_3d(:) = denitrif
    end associate

  end subroutine store_diagnostics_nitrification

  !***********************************************************************

  subroutine store_diagnostics_autotrophs(marbl_domain, autotroph_secondary_species,  &
       marbl_diags)

    use marbl_share_mod, only : autotroph_type

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    type(autotroph_secondary_species_type), dimension(:,:), intent(in) :: autotroph_secondary_species ! autotroph_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: n, k
    real(r8), dimension(km) :: delta_z

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if

    associate(diags => marbl_diags%diags(:))
      diags(marbl_interior_diag_ind%tot_CaCO3_form_zint)%field_2d = c0
      diags(marbl_interior_diag_ind%tot_bSi_form)%field_3d(:) = c0
      diags(marbl_interior_diag_ind%tot_Nfix)%field_3d(:) = c0
      diags(marbl_interior_diag_ind%tot_CaCO3_form)%field_3d(:) = c0
      do n = 1, autotroph_cnt
         diags(marbl_interior_diag_ind%N_lim(n))%field_3d(:) =                         &
                                       autotroph_secondary_species(n,:)%VNtot
         diags(marbl_interior_diag_ind%Fe_lim(n))%field_3d(:) =                        &
                                         autotroph_secondary_species(n,:)%VFe
         diags(marbl_interior_diag_ind%P_lim(n))%field_3d(:) =                         &
                                       autotroph_secondary_species(n,:)%VPtot

         if (marbl_interior_diag_ind%SiO3_lim(n).ne.-1) then
            diags(marbl_interior_diag_ind%SiO3_lim(n))%field_3d(:) =                   &
                                       autotroph_secondary_species(n,:)%VSiO3
         end if

         diags(marbl_interior_diag_ind%light_lim(n))%field_3d(:) =                     &
                                   autotroph_secondary_species(n,:)%light_lim
         diags(marbl_interior_diag_ind%photoNO3(n))%field_3d(:) =                      &
                                       autotroph_secondary_species(n,:)%NO3_V
         diags(marbl_interior_diag_ind%photoNH4(n))%field_3d(:) =                      &
                                       autotroph_secondary_species(n,:)%NH4_V
         diags(marbl_interior_diag_ind%PO4_uptake(n))%field_3d(:) =                    &
                                       autotroph_secondary_species(n,:)%PO4_V
         diags(marbl_interior_diag_ind%DOP_uptake(n))%field_3d(:) =                    &
                                       autotroph_secondary_species(n,:)%DOP_V
         diags(marbl_interior_diag_ind%photoFE(n))%field_3d(:) =                       &
                                     autotroph_secondary_species(n,:)%photoFe

         if (marbl_interior_diag_ind%bSi_form(n).ne.-1) then
           diags(marbl_interior_diag_ind%bSi_form(n))%field_3d(:) =                    &
                                     autotroph_secondary_species(n,:)%photoSi
           diags(marbl_interior_diag_ind%tot_bSi_form)%field_3d(:) =          &
                    diags(marbl_interior_diag_ind%tot_bSi_form)%field_3d(:) + &
                    diags(marbl_interior_diag_ind%bSi_form(n))%field_3d(:)
         endif

         if (marbl_interior_diag_ind%CaCO3_form(n).ne.-1) then
           diags(marbl_interior_diag_ind%CaCO3_form(n))%field_3d(:) =                  &
                                  autotroph_secondary_species(n,:)%CaCO3_PROD
           diags(marbl_interior_diag_ind%tot_CaCO3_form)%field_3d(:) =        &
                  diags(marbl_interior_diag_ind%tot_CaCO3_form)%field_3d(:) + &
                  diags(marbl_interior_diag_ind%CaCO3_form(n))%field_3d(:)
         end if

         if (marbl_interior_diag_ind%Nfix(n).ne.-1) then
           diags(marbl_interior_diag_ind%Nfix(n))%field_3d(:) =                        &
                                        autotroph_secondary_species(n,:)%Nfix
           diags(marbl_interior_diag_ind%tot_Nfix)%field_3d(:) =              &
                        diags(marbl_interior_diag_ind%tot_Nfix)%field_3d(:) + &
                        diags(marbl_interior_diag_ind%Nfix(n))%field_3d(:)
         end if

         diags(marbl_interior_diag_ind%auto_graze(n))%field_3d(:) =                    &
                                  autotroph_secondary_species(n,:)%auto_graze
         diags(marbl_interior_diag_ind%auto_graze_poc(n))%field_3d(:) =                &
                              autotroph_secondary_species(n,:)%auto_graze_poc
         diags(marbl_interior_diag_ind%auto_graze_doc(n))%field_3d(:) =                &
                              autotroph_secondary_species(n,:)%auto_graze_doc
         diags(marbl_interior_diag_ind%auto_graze_zoo(n))%field_3d(:) =                &
                              autotroph_secondary_species(n,:)%auto_graze_zoo
         diags(marbl_interior_diag_ind%auto_loss(n))%field_3d(:) =                     &
                                   autotroph_secondary_species(n,:)%auto_loss
         diags(marbl_interior_diag_ind%auto_loss_poc(n))%field_3d(:) =                 &
                               autotroph_secondary_species(n,:)%auto_loss_poc
         diags(marbl_interior_diag_ind%auto_loss_doc(n))%field_3d(:) =                 &
                               autotroph_secondary_species(n,:)%auto_loss_doc
         diags(marbl_interior_diag_ind%auto_agg(n))%field_3d(:) =                      &
                                    autotroph_secondary_species(n,:)%auto_agg
         diags(marbl_interior_diag_ind%photoC(n))%field_3d(:) =                        &
                                      autotroph_secondary_species(n,:)%photoC

         diags(marbl_interior_diag_ind%photoC_NO3(n))%field_3d(:) = c0
         where (autotroph_secondary_species(n,:)%VNtot > c0)
           diags(marbl_interior_diag_ind%photoC_NO3(n))%field_3d(:) =                  &
                autotroph_secondary_species(n,:)%photoC *                              &
                (autotroph_secondary_species(n,:)%VNO3 / autotroph_secondary_species(n,:)%VNtot)
         end where

         ! vertical integrals
         if (marbl_interior_diag_ind%CaCO3_form_zint(n).ne.-1) then
           diags(marbl_interior_diag_ind%CaCO3_form_zint(n))%field_2d = c0
           do k = 1,marbl_domain%kmt
             diags(marbl_interior_diag_ind%CaCO3_form_zint(n))%field_2d =              &
                    diags(marbl_interior_diag_ind%CaCO3_form_zint(n))%field_2d +       &
                    delta_z(k) * autotroph_secondary_species(n, k)%CaCO3_PROD
           end do
           diags(marbl_interior_diag_ind%tot_CaCO3_form_zint)%field_2d =      &
                diags(marbl_interior_diag_ind%tot_CaCO3_form_zint)%field_2d + &
                diags(marbl_interior_diag_ind%CaCO3_form_zint(n))%field_2d
         end if

         diags(marbl_interior_diag_ind%photoC_zint(n))%field_2d = c0
         diags(marbl_interior_diag_ind%photoC_NO3_zint(n))%field_2d = c0
         do k = 1,marbl_domain%kmt
           diags(marbl_interior_diag_ind%photoC_zint(n))%field_2d =                    &
                        diags(marbl_interior_diag_ind%photoC_zint(n))%field_2d +       &
                        delta_z(k) * autotroph_secondary_species(n, k)%photoC
           diags(marbl_interior_diag_ind%photoC_NO3_zint(n))%field_2d =                &
                 diags(marbl_interior_diag_ind%photoC_NO3_zint(n))%field_2d +          &
                 delta_z(k) * diags(marbl_interior_diag_ind%photoC_NO3(n))%field_3d(k)
         end do ! do k
      end do ! do n
    end associate

  end subroutine store_diagnostics_autotrophs

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_autotroph_sums(marbl_domain, autotroph_secondary_species,    &
       marbl_diags)

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    type(autotroph_secondary_species_type), dimension(:,:), intent(in) :: autotroph_secondary_species ! autotroph_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags


    integer(int_kind) :: n
    real(r8), dimension(km) :: delta_z

    if (partial_bottom_cells) then
          delta_z = marbl_domain%dzt
    else
          delta_z = marbl_domain%dz
    endif

    associate(diags => marbl_diags%diags)
      diags(marbl_interior_diag_ind%auto_graze_TOT)%field_3d(:) =                      &
                           sum(autotroph_secondary_species%auto_graze, dim=1)
      diags(marbl_interior_diag_ind%photoC_TOT)%field_3d(:) =                          &
                           sum(autotroph_secondary_species%photoC, dim=1)

      diags(marbl_interior_diag_ind%photoC_NO3_TOT)%field_3d(:) = c0
      do n = 1, autotroph_cnt
        where (autotroph_secondary_species(n,:)%VNtot > c0)
          diags(marbl_interior_diag_ind%photoC_NO3_TOT)%field_3d(:) =                  &
                           diags(marbl_interior_diag_ind%photoC_NO3_TOT)%field_3d(:) + &
                           (autotroph_secondary_species(n,:)%VNO3 /           &
                           autotroph_secondary_species(n,:)%VNtot) *          &
                           autotroph_secondary_species(n,:)%photoC
        end where
      end do

      diags(marbl_interior_diag_ind%photoC_TOT_zint)%field_2d = sum(                   &
                    delta_z * sum(autotroph_secondary_species%photoC, dim=1))
      diags(marbl_interior_diag_ind%photoC_NO3_TOT_zint)%field_2d = sum(               &
                    delta_z * diags(marbl_interior_diag_ind%photoC_NO3_TOT)%field_3d)
    end associate

  end subroutine store_diagnostics_autotroph_sums

  !***********************************************************************

  subroutine store_diagnostics_particulates(marbl_domain, POC, P_CaCO3,       &
                                            P_SiO2, dust, P_iron,             &
                                            sed_denitrif, other_remin,        &
                                            marbl_diags)
    !-----------------------------------------------------------------------
    ! - Set tavg variables.
    ! - Accumulte losses of BGC tracers to sediments
    !-----------------------------------------------------------------------
    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_parms, only : Q
    use marbl_parms, only : Qp_zoo_pom

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    type(column_sinking_particle_type), intent(in) :: POC
    type(column_sinking_particle_type), intent(in) :: P_CaCO3
    type(column_sinking_particle_type), intent(in) :: P_SiO2
    type(column_sinking_particle_type), intent(in) :: dust
    type(column_sinking_particle_type), intent(in) :: P_iron
    real(r8), dimension(:), intent(in) :: sed_denitrif ! km
    real(r8), dimension(:), intent(in) :: other_remin  ! km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    real(r8), dimension(km) :: delta_z

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if
    associate(diags => marbl_diags%diags)
      diags(marbl_interior_diag_ind%POC_FLUX_IN)%field_3d(:) = POC%sflux_in +          &
                                                      POC%hflux_in
      diags(marbl_interior_diag_ind%POC_PROD)%field_3d(:) = POC%prod
      diags(marbl_interior_diag_ind%POC_REMIN)%field_3d(:) = POC%remin

      diags(marbl_interior_diag_ind%CaCO3_FLUX_IN)%field_3d(:) = P_CaCO3%sflux_in +    &
                                                        P_CaCO3%hflux_in
      diags(marbl_interior_diag_ind%CaCO3_PROD)%field_3d(:) = P_CaCO3%prod
      diags(marbl_interior_diag_ind%CaCO3_REMIN)%field_3d(:) = P_CaCO3%remin

      diags(marbl_interior_diag_ind%SiO2_FLUX_IN)%field_3d(:) = P_SiO2%sflux_in +      &
                                                       P_SiO2%hflux_in
      diags(marbl_interior_diag_ind%SiO2_PROD)%field_3d(:) = P_SiO2%prod
      diags(marbl_interior_diag_ind%SiO2_REMIN)%field_3d(:) = P_SiO2%remin

      diags(marbl_interior_diag_ind%dust_FLUX_IN)%field_3d(:) = dust%sflux_in +        &
                                                       dust%hflux_in
      diags(marbl_interior_diag_ind%dust_REMIN)%field_3d(:) = P_SiO2%remin

      diags(marbl_interior_diag_ind%P_iron_FLUX_IN)%field_3d(:) = P_iron%sflux_in +    &
                                                         P_iron%hflux_in
      diags(marbl_interior_diag_ind%P_iron_PROD)%field_3d(:) = P_iron%prod
      diags(marbl_interior_diag_ind%P_iron_REMIN)%field_3d(:) = P_iron%remin

      diags(marbl_interior_diag_ind%calcToSed)%field_2d = sum(P_CaCO3%sed_loss)
      diags(marbl_interior_diag_ind%bsiToSed)%field_2d = sum(P_SiO2%sed_loss)
      diags(marbl_interior_diag_ind%pocToSed)%field_2d = sum(POC%sed_loss)
      diags(marbl_interior_diag_ind%SedDenitrif)%field_2d = sum(sed_denitrif * delta_z)
      diags(marbl_interior_diag_ind%OtherRemin)%field_2d = sum(other_remin * delta_z)
      diags(marbl_interior_diag_ind%ponToSed)%field_2d = sum(POC%sed_loss * Q)
      diags(marbl_interior_diag_ind%popToSed)%field_2d = sum(POC%sed_loss * Qp_zoo_pom)
      diags(marbl_interior_diag_ind%dustToSed)%field_2d = sum(dust%sed_loss)
      diags(marbl_interior_diag_ind%pfeToSed)%field_2d = sum(P_iron%sed_loss)
    end associate

  end subroutine store_diagnostics_particulates

   !***********************************************************************

  subroutine store_diagnostics_oxygen(marbl_domain, marbl_gcm_state, column_zt, &
                                      column_o2, o2_production, o2_consumption,&
                                      marbl_diags)

    use marbl_oxygen, only : o2sat_scalar

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    type(marbl_gcm_state_type), intent(in) :: marbl_gcm_state
    real(r8), dimension(:), intent(in) :: column_zt
    real(r8), dimension(:), intent(in) :: column_o2
    real(r8), dimension(:), intent(in) :: o2_production
    real(r8), dimension(:), intent(in) :: o2_consumption
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: k, min_ind

    ! Find min_o2 and depth_min_o2
    associate(diags => marbl_diags%diags(:))
      min_ind = minloc(column_o2(1:marbl_domain%kmt), dim=1)
      diags(marbl_interior_diag_ind%O2_ZMIN)%field_2d = column_o2(min_ind)
      diags(marbl_interior_diag_ind%O2_ZMIN_DEPTH)%field_2d = column_zt(min_ind)

      diags(marbl_interior_diag_ind%O2_PRODUCTION)%field_3d(:) = o2_production
      diags(marbl_interior_diag_ind%O2_CONSUMPTION)%field_3d(:) = o2_consumption

      diags(marbl_interior_diag_ind%AOU)%field_3d(:) = -column_o2
      do k=1,marbl_domain%kmt
        if (marbl_domain%land_mask) then
          diags(marbl_interior_diag_ind%AOU)%field_3d(k) =                             &
        O2SAT_scalar(marbl_gcm_state%temperature(k), marbl_gcm_state%salinity(k)) - &
        column_o2(k)
        end if
      end do
    end associate

  end subroutine store_diagnostics_oxygen

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_photosynthetically_available_radiation(PAR, marbl_diags)

    type(photosynthetically_available_radiation_type), dimension(:), intent(in) :: PAR
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_interior_diag_ind%PAR_avg)%field_3d(:) = PAR%avg
    end associate

  end subroutine store_diagnostics_photosynthetically_available_radiation

  !***********************************************************************

  subroutine store_diagnostics_misc(marbl_diags)
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags
  end subroutine store_diagnostics_misc

  !***********************************************************************

  subroutine store_diagnostics_zooplankton(zooplankton_secondary_species, marbl_diags)

    type(zooplankton_secondary_species_type), dimension(:,:), intent(in) ::   &
                                              zooplankton_secondary_species
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: n

    associate(diags => marbl_diags%diags(:))
      do n = 1, zooplankton_cnt
        diags(marbl_interior_diag_ind%zoo_loss(n))%field_3d(:) =                       &
                                  zooplankton_secondary_species(n,:)%zoo_loss
        diags(marbl_interior_diag_ind%zoo_loss_poc(n))%field_3d(:) =                   &
                              zooplankton_secondary_species(n,:)%zoo_loss_poc
        diags(marbl_interior_diag_ind%zoo_loss_doc(n))%field_3d(:) =                   &
                              zooplankton_secondary_species(n,:)%zoo_loss_doc
        diags(marbl_interior_diag_ind%zoo_graze(n))%field_3d(:) =                      &
                                 zooplankton_secondary_species(n,:)%zoo_graze
        diags(marbl_interior_diag_ind%zoo_graze_poc(n))%field_3d(:) =                  &
                             zooplankton_secondary_species(n,:)%zoo_graze_poc
        diags(marbl_interior_diag_ind%zoo_graze_doc(n))%field_3d(:) =                  &
                             zooplankton_secondary_species(n,:)%zoo_graze_doc
        diags(marbl_interior_diag_ind%zoo_graze_zoo(n))%field_3d(:) =                  &
                             zooplankton_secondary_species(n,:)%zoo_graze_zoo
        diags(marbl_interior_diag_ind%x_graze_zoo(n))%field_3d(:) =                    &
                               zooplankton_secondary_species(n,:)%x_graze_zoo
      end do
    end associate

  end subroutine store_diagnostics_zooplankton

  !***********************************************************************

  subroutine store_diagnostics_dissolved_organic_matter(&
       dissolved_organic_matter, fe_scavenge, fe_scavenge_rate, marbl_diags)

    type(dissolved_organic_matter_type), dimension(:), intent(in) ::          &
                                            dissolved_organic_matter
    real(r8), dimension(:), intent(in) :: fe_scavenge, fe_scavenge_rate
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_interior_diag_ind%DOC_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DOC_prod
      diags(marbl_interior_diag_ind%DOC_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOC_remin
      diags(marbl_interior_diag_ind%DON_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DON_prod
      diags(marbl_interior_diag_ind%DON_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DON_remin
      diags(marbl_interior_diag_ind%DOP_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DOP_prod
      diags(marbl_interior_diag_ind%DOP_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOP_remin
      diags(marbl_interior_diag_ind%DONr_remin)%field_3d(:) =                          &
                                          dissolved_organic_matter%DONr_remin
      diags(marbl_interior_diag_ind%DOPr_remin)%field_3d(:) =                          &
                                          dissolved_organic_matter%DOPr_remin
      diags(marbl_interior_diag_ind%DOFe_prod)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOFe_prod
      diags(marbl_interior_diag_ind%DOFe_remin)%field_3d(:) =                          &
                                          dissolved_organic_matter%DOFe_remin

      diags(marbl_interior_diag_ind%Fe_scavenge)%field_3d(:) = Fe_scavenge
      diags(marbl_interior_diag_ind%Fe_scavenge_rate)%field_3d(:) = Fe_scavenge_rate
    end associate

  end subroutine store_diagnostics_dissolved_organic_matter

  !***********************************************************************

  subroutine store_diagnostics_carbon_fluxes(marbl_domain, zw, POC, P_CaCO3,  &
                                             dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_share_mod, only : autotroph_type, zooplankton_type

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zw  ! km
    type(column_sinking_particle_type), intent(in) :: POC
    type(column_sinking_particle_type), intent(in) :: P_CaCO3
    real(r8), dimension(:,:), intent(in) :: dtracer ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: k, n, auto_ind
    real(r8), dimension(km) :: delta_z
    real(r8) :: ztop, work1

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if

    associate(diags => marbl_diags%diags(:))
      ! vertical integrals
      diags(marbl_interior_diag_ind%Jint_Ctot)%field_2d = c0
      diags(marbl_interior_diag_ind%Jint_100m_Ctot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(dic_ind,k) + dtracer(doc_ind,k) +                     &
                sum(dtracer(zooplankton(:)%C_ind,k)) +                        &
                sum(dtracer(autotrophs(:)%C_ind,k))
        do auto_ind = 1, autotroph_cnt
          n = autotrophs(auto_ind)%CaCO3_ind
          if (n > 0) then
            work1 = work1 + dtracer(n,k)
          endif
        end do
        diags(marbl_interior_diag_ind%Jint_Ctot)%field_2d =                            &
                 diags(marbl_interior_diag_ind%Jint_Ctot)%field_2d +                   &
                 delta_z(k) * work1 + (POC%sed_loss(k) + P_CaCO3%sed_loss(k))

        if (ztop < 100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Ctot)%field_2d =                     &
                              diags(marbl_interior_diag_ind%Jint_100m_Ctot)%field_2d + &
                              min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Ctot)%field_2d =                     &
                              diags(marbl_interior_diag_ind%Jint_100m_Ctot)%field_2d + &
                              (POC%sed_loss(k) + P_CaCO3%sed_loss(k))
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_carbon_fluxes

  !***********************************************************************

  subroutine store_diagnostics_nitrogen_fluxes(marbl_domain, zw, POC, denitrif,&
                                    sed_denitrif, autotroph_secondary_species, &
                                    dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type

    use marbl_parms     , only : Q

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zw  ! km
    type(column_sinking_particle_type), intent(in) :: POC
    real(r8), dimension(:), intent(in) :: denitrif ! km
    real(r8), dimension(:), intent(in) :: sed_denitrif ! km
    type(autotroph_secondary_species_type), dimension(:,:), intent(in) ::     &
                                              autotroph_secondary_species
    real(r8), dimension(:,:), intent(in) :: dtracer ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: k, n
    real(r8), dimension(km) :: delta_z
    real(r8) :: ztop, work1

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if

    associate(diags => marbl_diags%diags(:))
      ! vertical integrals
      diags(marbl_interior_diag_ind%Jint_Ntot)%field_2d = c0
      diags(marbl_interior_diag_ind%Jint_100m_Ntot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(no3_ind,k) + dtracer(nh4_ind,k) +                     &
                dtracer(don_ind,k) + dtracer(donr_ind,k) +                    &
                Q * sum(dtracer(zooplankton(:)%C_ind,k)) +                    &
                Q * sum(dtracer(autotrophs(:)%C_ind,k))
        ! add back column and sediment denitrification
        work1 = work1 + denitrif(k) + sed_denitrif(k)
        ! subtract out N fixation
        do n = 1, autotroph_cnt
          if (autotrophs(n)%Nfixer) then
            work1 = work1 - autotroph_secondary_species(n,k)%Nfix
          end if
        end do
        diags(marbl_interior_diag_ind%Jint_Ntot)%field_2d =                            &
                                   diags(marbl_interior_diag_ind%Jint_Ntot)%field_2d + &
                                   delta_z(k) * work1 + POC%sed_loss(k) * Q

        if (ztop < 100.0e2_r8) then
            diags(marbl_interior_diag_ind%Jint_100m_Ntot)%field_2d =                   &
                              diags(marbl_interior_diag_ind%Jint_100m_Ntot)%field_2d + &
                           min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
            diags(marbl_interior_diag_ind%Jint_100m_Ntot)%field_2d =                   &
                              diags(marbl_interior_diag_ind%Jint_100m_Ntot)%field_2d + &
                              POC%sed_loss(k) * Q
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_nitrogen_fluxes

  !***********************************************************************

  subroutine store_diagnostics_phosphorus_fluxes(marbl_domain, zw, POC, &
       dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_share_mod, only : autotroph_type, zooplankton_type

    use marbl_parms     , only : Q
    use marbl_parms     , only : Qp_zoo_pom

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zw  ! km
    type(column_sinking_particle_type), intent(in) :: POC
    real(r8), dimension(:,:), intent(in) :: dtracer ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags


    integer(int_kind) :: k, n
    real(r8), dimension(km) :: delta_z
    real(r8) :: ztop, work1

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if

    associate(diags => marbl_diags%diags(:))
      ! vertical integrals
      diags(marbl_interior_diag_ind%Jint_Ptot)%field_2d = c0
      diags(marbl_interior_diag_ind%Jint_100m_Ptot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(po4_ind,k) + dtracer(dop_ind,k) + dtracer(dopr_ind,k)
        do n = 1, zooplankton_cnt
          work1 = work1 + Qp_zoo_pom * dtracer(zooplankton(n)%C_ind,k)
        end do
        do n = 1, autotroph_cnt
          work1 = work1 + autotrophs(n)%Qp * dtracer(autotrophs(n)%C_ind,k)
        end do
        diags(marbl_interior_diag_ind%Jint_Ptot)%field_2d =                            &
                            diags(marbl_interior_diag_ind%Jint_Ptot)%field_2d +        &
                            delta_z(k) * work1 + POC%sed_loss(k) * Qp_zoo_pom

        if (ztop < 100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Ptot)%field_2d =                     &
                              diags(marbl_interior_diag_ind%Jint_100m_Ptot)%field_2d + &
                              min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Ptot)%field_2d =                     &
                              diags(marbl_interior_diag_ind%Jint_100m_Ptot)%field_2d + &
                              POC%sed_loss(k) * Qp_zoo_pom
        end if
        ztop = zw(k)
    end do
    end associate

  end subroutine store_diagnostics_phosphorus_fluxes

  !***********************************************************************

  subroutine store_diagnostics_silicon_fluxes(marbl_domain, zw, P_SiO2, &
                                              dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_parms    , only : Q
    use marbl_parms    , only : Qp_zoo_pom

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zw  ! km
    type(column_sinking_particle_type), intent(in) :: P_SiO2
    real(r8), dimension(:,:), intent(in) :: dtracer ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: k, n
    real(r8), dimension(km) :: delta_z
    real(r8) :: ztop, work1

    if (partial_bottom_cells) then
       delta_z = marbl_domain%dzt
    else
       delta_z = marbl_domain%dz
    end if

    associate(diags => marbl_diags%diags(:))
      ! vertical integrals
      diags(marbl_interior_diag_ind%Jint_Sitot)%field_2d = c0
      diags(marbl_interior_diag_ind%Jint_100m_Sitot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(sio3_ind,k)
        do n = 1, autotroph_cnt
          if (autotrophs(n)%Si_ind > 0) then
            work1 = work1 + dtracer(autotrophs(n)%Si_ind,k)
          end if
        end do
        diags(marbl_interior_diag_ind%Jint_Sitot)%field_2d =                           &
                                  diags(marbl_interior_diag_ind%Jint_Sitot)%field_2d + &
                                  delta_z(k) * work1 + P_SiO2%sed_loss(k)

        if (ztop < 100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Sitot)%field_2d =                    &
                             diags(marbl_interior_diag_ind%Jint_100m_Sitot)%field_2d + &
                             min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_interior_diag_ind%Jint_100m_Sitot)%field_2d =                    &
                             diags(marbl_interior_diag_ind%Jint_100m_Sitot)%field_2d + &
                             P_SiO2%sed_loss(k)
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_silicon_fluxes

end module ecosys_diagnostics_mod
