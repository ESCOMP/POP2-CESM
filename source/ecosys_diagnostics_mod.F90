! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! Will be part of MARBL library
module ecosys_diagnostics_mod

  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : log_kind
  use marbl_kinds_mod, only : char_len

  use marbl_share_mod, only : autotrophs
  use marbl_share_mod, only : autotroph_cnt
  use marbl_share_mod, only : zooplankton
  use marbl_share_mod, only : zooplankton_cnt
  use marbl_share_mod, only : ecosys_tracer_cnt
  use marbl_share_mod, only : ecosys_ciso_tracer_cnt

  use marbl_interface_types, only : carbonate_type
  use marbl_interface_types, only : zooplankton_secondary_species_type
  use marbl_interface_types, only : autotroph_secondary_species_type
  use marbl_interface_types, only : photosynthetically_available_radiation_type
  use marbl_interface_types, only : dissolved_organic_matter_type
  use marbl_interface_types, only : marbl_diagnostics_type
  use marbl_interface_types, only : marbl_column_domain_type
  use marbl_interface_types, only : marbl_gcm_state_type

  use marbl_parms, only : c0
  use marbl_parms, only : c1
  use marbl_parms, only : po4_ind
  use marbl_parms, only : no3_ind
  use marbl_parms, only : sio3_ind
  use marbl_parms, only : nh4_ind
  use marbl_parms, only : fe_ind
  use marbl_parms, only : dic_ind
  use marbl_parms, only : doc_ind
  use marbl_parms, only : don_ind
  use marbl_parms, only : dop_ind
  use marbl_parms, only : dopr_ind
  use marbl_parms, only : donr_ind
  use marbl_parms, only : docr_ind

  Implicit None
  Public

  !-----------------------------------------------------------------------
  !  Largest possible size for each class of diagnostics
  !-----------------------------------------------------------------------

  ! FIXME - the following should be counted and not be parameters
  integer, public, parameter :: max_interior_diags = 112 + (40*autotroph_cnt) + (8*zooplankton_cnt)
  integer, public, parameter :: max_forcing_diags  = 40
  integer, public, parameter :: max_restore_diags  = ecosys_tracer_cnt

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
    integer(int_kind) :: Jint_Fetot
    integer(int_kind) :: Jint_100m_Fetot

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
    integer(int_kind) :: DOCr_remin
    integer(int_kind) :: DON_prod
    integer(int_kind) :: DON_remin
    integer(int_kind) :: DONr_remin
    integer(int_kind) :: DOP_prod
    integer(int_kind) :: DOP_remin
    integer(int_kind) :: DOPr_remin
    integer(int_kind) :: Fe_scavenge
    integer(int_kind) :: Fe_scavenge_rate

    ! Particulate 3D diags
    integer(int_kind) :: POC_FLUX_IN
    integer(int_kind) :: POC_PROD
    integer(int_kind) :: POC_REMIN
    integer(int_kind) :: POC_REMIN_DIC
    integer(int_kind) :: PON_REMIN_NH4
    integer(int_kind) :: POP_REMIN_PO4
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

     !  ciso ids for nonstandard 3d fields
     integer (int_kind) :: CISO_PO13C_FLUX_IN                                 ! po13c flux into cell
     integer (int_kind) :: CISO_PO14C_FLUX_IN                                 ! po14c flux into cell
     integer (int_kind) :: CISO_PO13C_PROD                                    ! po13c production
     integer (int_kind) :: CISO_PO14C_PROD                                    ! po14c production
     integer (int_kind) :: CISO_PO13C_REMIN                                   ! po13c remineralization
     integer (int_kind) :: CISO_PO14C_REMIN                                   ! po14c remineralization
     integer (int_kind) :: CISO_Ca13CO3_PROD                                  ! ca13co3 production
     integer (int_kind) :: CISO_Ca14CO3_PROD                                  ! ca14co3 production
     integer (int_kind) :: CISO_Ca13CO3_REMIN                                 ! ca13co3 remineralization
     integer (int_kind) :: CISO_Ca14CO3_REMIN                                 ! ca14co3 remineralization
     integer (int_kind) :: CISO_Ca13CO3_FLUX_IN                               ! ca13co3 flux into cell
     integer (int_kind) :: CISO_Ca14CO3_FLUX_IN                               ! ca14co3 flux into cell
     integer (int_kind) :: CISO_photo13C_TOT                                  ! total 13C fixation
     integer (int_kind) :: CISO_photo14C_TOT                                  ! total 14C fixation
     integer (int_kind) :: CISO_photo13C_TOT_zint                             ! total 13C fixation vertical integral
     integer (int_kind) :: CISO_photo14C_TOT_zint                             ! total 14C fixation vertical integral

     ! ciso ids for  MORE nonstandard 3d fields
     integer (int_kind), dimension(autotroph_cnt) :: CISO_eps_autotroph       ! epsilon for each autotroph
     integer (int_kind), dimension(autotroph_cnt) :: CISO_mui_to_co2star      ! mui_to_co2star for each autotroph
     integer (int_kind), dimension(autotroph_cnt) :: CISO_Ca13CO3_form        ! Ca13CO3 formation
     integer (int_kind), dimension(autotroph_cnt) :: CISO_Ca14CO3_form        ! Ca14CO3 formation
     integer (int_kind), dimension(autotroph_cnt) :: CISO_Ca13CO3_form_zint   ! Ca13CO3 formation vertical integral 0-100 m
     integer (int_kind), dimension(autotroph_cnt) :: CISO_Ca14CO3_form_zint   ! Ca14CO3 formation vertical integral 0-100 m
     integer (int_kind), dimension(autotroph_cnt) :: CISO_photo13C            ! 13C fixation
     integer (int_kind), dimension(autotroph_cnt) :: CISO_photo14C            ! 14C fixation
     integer (int_kind), dimension(autotroph_cnt) :: CISO_photo13C_zint       ! 13C fixation vertical integral
     integer (int_kind), dimension(autotroph_cnt) :: CISO_photo14C_zint       ! 14C fixation vertical integral
     integer (int_kind), dimension(autotroph_cnt) :: CISO_d13C                ! if for d13C of autotroph carbon
     integer (int_kind), dimension(autotroph_cnt) :: CISO_d14C                ! if for d14C of autotroph carbon
     integer (int_kind), dimension(autotroph_cnt) :: CISO_autotrophCaCO3_d14C ! if for d14C of autotrophCaCO3
     integer (int_kind), dimension(autotroph_cnt) :: CISO_autotrophCaCO3_d13C ! if for d13C of autotrophCaCO3

     integer (int_kind) :: CISO_eps_aq_g                                      ! eps_aq_g
     integer (int_kind) :: CISO_eps_dic_g                                     ! eps_dic_g
     integer (int_kind) :: CISO_DO13C_prod                                    ! do13c production
     integer (int_kind) :: CISO_DO14C_prod                                    ! do14c production
     integer (int_kind) :: CISO_DO13C_remin                                   ! do13c remineralization
     integer (int_kind) :: CISO_DO14C_remin                                   ! do14c remineralization
     integer (int_kind) :: CISO_Jint_13Ctot                                   ! vertically integrated source sink term, 13Ctot
     integer (int_kind) :: CISO_Jint_14Ctot                                   ! vertically integrated source sink term, 14Ctot
     integer (int_kind) :: CISO_Jint_100m_13Ctot                              ! vertically integrated source sink term, 0-100m, 13Ctot !FIXME - this is not being computed now
     integer (int_kind) :: CISO_Jint_100m_14Ctot                              ! vertically integrated source sink term, 0-100m, 14Ctot !FIXME - this is not being computed now
     integer (int_kind) :: CISO_zooC_d13C                                     ! if for d13C of zooC
     integer (int_kind) :: CISO_zooC_d14C                                     ! if for d14C of zooC
     integer (int_kind) :: CISO_DOC_d13C                                      ! if for d13C of DOC
     integer (int_kind) :: CISO_DOC_d14C                                      ! if for d14C of DOC
     integer (int_kind) :: CISO_DIC_d13C                                      ! if for d13C of DIC
     integer (int_kind) :: CISO_DIC_d14C                                      ! if for d14C of DIC
     integer (int_kind) :: calcToSed_13C                                      ! calcite flux sedimentary burial
     integer (int_kind) :: calcToSed_14C                                      ! calcite flux sedimentary burial
     integer (int_kind) :: pocToSed_13C                                       ! poc burial flux to sediments
     integer (int_kind) :: pocToSed_14C                                       ! poc burial flux to sediments

  end type marbl_interior_diagnostics_indexing_type
  type(marbl_interior_diagnostics_indexing_type), public :: marbl_interior_diag_ind

  !***********************************************************************

  type marbl_forcing_diagnostics_indexing_type
     integer(int_kind) :: ECOSYS_IFRAC
     integer(int_kind) :: ECOSYS_XKW
     integer(int_kind) :: ECOSYS_ATM_PRESS
     integer(int_kind) :: PV_O2
     integer(int_kind) :: SCHMIDT_O2
     integer(int_kind) :: O2SAT
     integer(int_kind) :: O2_GAS_FLUX
     integer(int_kind) :: CO2STAR
     integer(int_kind) :: DCO2STAR
     integer(int_kind) :: pCO2SURF
     integer(int_kind) :: DpCO2
     integer(int_kind) :: PV_CO2
     integer(int_kind) :: SCHMIDT_CO2
     integer(int_kind) :: DIC_GAS_FLUX
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
     integer(int_kind) :: DON_RIV_FLUX
     integer(int_kind) :: DONr_RIV_FLUX
     integer(int_kind) :: DOP_RIV_FLUX
     integer(int_kind) :: DOPr_RIV_FLUX
     integer(int_kind) :: DSI_RIV_FLUX
     integer(int_kind) :: DFE_RIV_FLUX
     integer(int_kind) :: DIC_RIV_FLUX
     integer(int_kind) :: ALK_RIV_FLUX
     integer(int_kind) :: DOC_RIV_FLUX
     integer(int_kind) :: DOCr_RIV_FLUX
  end type marbl_forcing_diagnostics_indexing_type
  type(marbl_forcing_diagnostics_indexing_type), public :: marbl_forcing_diag_ind

  !***********************************************************************

contains

  !***********************************************************************

  subroutine marbl_diagnostics_init( &
       marbl_interior_diags,                &
       marbl_restore_diags,                 &
       marbl_forcing_diags,                 &
       num_elements_interior,               &
       num_elements_forcing,                &
       num_levels,                          &
       tracer_d_module, ciso_on)

    use marbl_interface_types , only : marbl_tracer_metadata_type

    ! intent(in)s
    integer                          , intent(in) :: num_elements_interior
    integer                          , intent(in) :: num_elements_forcing
    integer                          , intent(in) :: num_levels
    ! descriptors for each tracer
    type (marbl_tracer_metadata_type), intent(in) :: tracer_d_module(:)
    logical (log_kind)               , intent(in) :: ciso_on 

    ! intent(inout)s
    type(marbl_diagnostics_type), intent(inout) :: marbl_interior_diags
    type(marbl_diagnostics_type), intent(inout) :: marbl_restore_diags
    type(marbl_diagnostics_type), intent(inout) :: marbl_forcing_diags

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
    call marbl_forcing_diags%construct(max_forcing_diags, num_elements_forcing, num_levels)

    associate(                          &
         ind => marbl_forcing_diag_ind, &
         diags => marbl_forcing_diags   &
         )

    lname    = 'Ice Fraction for ecosys fluxes'
    sname    = 'ECOSYS_IFRAC'
    units    = 'fraction'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ECOSYS_IFRAC)

    lname    = 'XKW for ecosys fluxes'
    sname    = 'ECOSYS_XKW'
    units    = 'cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ECOSYS_XKW)

    lname    = 'Atmospheric Pressure for ecosys fluxes'
    sname    = 'ECOSYS_ATM_PRESS'
    units    = 'atmospheres'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ECOSYS_ATM_PRESS)

    lname    = 'PV_O2'
    sname    = 'PV_O2'
    units    = 'cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PV_O2)

    lname    = 'O2 Schmidt Number'
    sname    = 'SCHMIDT_O2'
    units    = 'none'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SCHMIDT_O2)

    lname    = 'O2 Saturation'
    sname    = 'O2SAT'
    units    = 'mmol/m^3'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2SAT)

    lname    = 'Dissolved Oxygen Surface Flux'
    sname    = 'STF_O2'
    units    = 'mmol/m^3 cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2_GAS_FLUX)

    lname    = 'CO2 Star'
    sname    = 'CO2STAR'
    units    = 'mmol/m^3'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CO2STAR)

    lname    = 'D CO2 Star'
    sname    = 'DCO2STAR'
    units    = 'mmol/m^3'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DCO2STAR)

    lname    = 'surface pCO2'
    sname    = 'pCO2SURF'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pCO2SURF)

    lname    = 'D pCO2'
    sname    = 'DpCO2'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DpCO2)

    lname    = 'CO2 Piston Velocity'
    sname    = 'PV_CO2'
    units    = 'cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PV_CO2)

    lname    = 'CO2 Schmidt Number'
    sname    = 'SCHMIDT_CO2'
    units    = 'none'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SCHMIDT_CO2)

    lname    = 'DIC Surface Gas Flux'
    sname    = 'FG_CO2'
    units    = 'mmol/m^3 cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DIC_GAS_FLUX)

    lname    = 'Surface pH'
    sname    = 'PH'
    units    = 'none'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PH)

    lname    = 'Atmospheric CO2'
    sname    = 'ATM_CO2'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ATM_CO2)

    lname    = 'CO2 Star, Alternative CO2'
    sname    = 'CO2STAR_ALT_CO2'
    units    = 'mmol/m^3'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CO2STAR_ALT_CO2)

    lname    = 'D CO2 Star, Alternative CO2'
    sname    = 'DCO2STAR_ALT_CO2'
    units    = 'mmol/m^3'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DCO2STAR_ALT_CO2)

    lname    = 'surface pCO2, Alternative CO2'
    sname    = 'pCO2SURF_ALT_CO2'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pCO2SURF_ALT_CO2)

    lname    = 'D pCO2, Alternative CO2'
    sname    = 'DpCO2_ALT_CO2'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DpCO2_ALT_CO2)

    lname    = 'DIC Surface Gas Flux, Alternative CO2'
    sname    = 'FG_ALT_CO2'
    units    = 'mmol/m^3 cm/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DIC_GAS_FLUX_ALT_CO2)

    lname    = 'Surface pH, Alternative CO2'
    sname    = 'PH_ALT_CO2'
    units    = 'none'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PH_ALT_CO2)

    lname    = 'Atmospheric Alternative CO2'
    sname    = 'ATM_ALT_CO2'
    units    = 'ppmv'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ATM_ALT_CO2)

    lname    = 'Atmospheric Iron Flux'
    sname    = 'IRON_FLUX'
    units    = 'mmol/m^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%IRON_FLUX)

    lname    = 'Dust Flux'
    sname    = 'DUST_FLUX'
    units    = 'g/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DUST_FLUX)

    lname    = 'Flux of NOx from Atmosphere'
    sname    = 'NOx_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%NOx_FLUX)

    lname    = 'Flux of NHy from Atmosphere'
    sname    = 'NHy_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%NHy_FLUX)

    lname    = 'Flux of DIN from rivers'
    sname    = 'DIN_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DIN_RIV_FLUX)

    lname    = 'Flux of DIP from rivers'
    sname    = 'DIP_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DIP_RIV_FLUX)

    lname    = 'Flux of DON from rivers'
    sname    = 'DON_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DON_RIV_FLUX)

    lname    = 'Flux of DONr from rivers'
    sname    = 'DONr_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DONr_RIV_FLUX)

    lname    = 'Flux of DOP from rivers'
    sname    = 'DOP_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOP_RIV_FLUX)

    lname    = 'Flux of DOPr from rivers'
    sname    = 'DOPr_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOPr_RIV_FLUX)

    lname    = 'Flux of DSI from rivers'
    sname    = 'DSI_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DSI_RIV_FLUX)

    lname    = 'Flux of DFE from rivers'
    sname    = 'DFE_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DFE_RIV_FLUX)

    lname    = 'Flux of DIC from rivers'
    sname    = 'DIC_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DIC_RIV_FLUX)

    lname    = 'Flux of ALK from rivers'
    sname    = 'ALK_RIV_FLUX'
    units    = 'alk/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ALK_RIV_FLUX)

    lname    = 'Flux of DOC from rivers'
    sname    = 'DOC_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOC_RIV_FLUX)

    lname    = 'Flux of DOCr from rivers'
    sname    = 'DOCr_RIV_FLUX'
    units    = 'nmol/cm^2/s'
    vgrid    = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOCr_RIV_FLUX)

    end associate

    !-----------------------------------------------------------------
    ! Interior diagnostics
    !-----------------------------------------------------------------

    ! Allocate memory for interior diagnostics
    call marbl_interior_diags%construct(max_interior_diags, num_elements_interior, num_levels)

    associate(                           &
         ind => marbl_interior_diag_ind, &
         diags => marbl_interior_diags   &
         )

    ! General 2D diags
    lname = 'Calcite Saturation Depth'
    sname = 'zsatcalc'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zsatcalc)

    lname = 'Aragonite Saturation Depth'
    sname = 'zsatarag'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zsatarag)

    lname = 'Vertical Minimum of O2'
    sname = 'O2_ZMIN'
    units = 'mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2_ZMIN)

    lname = 'Depth of Vertical Minimum of O2'
    sname = 'O2_ZMIN_DEPTH'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2_ZMIN_DEPTH)

    lname = 'Total C Fixation Vertical Integral'
    sname = 'photoC_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_TOT_zint)

    lname = 'Total C Fixation from NO3 Vertical Integral'
    sname = 'photoC_NO3_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_NO3_TOT_zint)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot'
    sname = 'Jint_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot, 0-100m'
    sname = 'Jint_100m_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_100m_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot'
    sname = 'Jint_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot, 0-100m'
    sname = 'Jint_100m_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_100m_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot'
    sname = 'Jint_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot, 0-100m'
    sname = 'Jint_100m_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_100m_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot'
    sname = 'Jint_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_Sitot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot, 0-100m'
    sname = 'Jint_100m_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_100m_Sitot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Fetot'
    sname = 'Jint_Fetot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_Fetot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Fetot, 0-100m'
    sname = 'Jint_100m_Fetot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Jint_100m_Fetot)

    ! Particulate 2D diags
    lname = 'CaCO3 Flux to Sediments'
    sname = 'calcToSed'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%calcToSed)

    lname = 'POC Flux to Sediments'
    sname = 'pocToSed'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pocToSed)

    lname = 'nitrogen burial Flux to Sediments'
    sname = 'ponToSed'
    units = 'nmolN/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ponToSed)

    lname = 'nitrogen loss in Sediments'
    sname = 'SedDenitrif'
    units = 'nmolN/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SedDenitrif)

    lname = 'non-oxic,non-dentr remin in Sediments'
    sname = 'OtherRemin'
    units = 'nmolC/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%OtherRemin)

    lname = 'phosphorus Flux to Sediments'
    sname = 'popToSed'
    units = 'nmolP/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%popToSed)

    lname = 'biogenic Si Flux to Sediments'
    sname = 'bsiToSed'
    units = 'nmolSi/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%bsiToSed)

    lname = 'dust Flux to Sediments'
    sname = 'dustToSed'
    units = 'g/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%dustToSed)

    lname = 'pFe Flux to Sediments'
    sname = 'pfeToSed'
    units = 'nmolFe/cm^2/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pfeToSed)

    ! Autotroph 2D diags
    do n=1,autotroph_cnt
       lname = trim(autotrophs(n)%lname) // ' C Fixation Vertical Integral'
       sname = 'photoC_' // trim(autotrophs(n)%sname) // '_zint'
       units = 'mmol/m^3 cm/s'
       vgrid = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_zint(n))

       lname = trim(autotrophs(n)%lname) // ' C Fixation from NO3 Vertical Integral'
       sname = 'photoC_NO3_' // trim(autotrophs(n)%sname) // '_zint'
       units = 'mmol/m^3 cm/s'
       vgrid = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_NO3_zint(n))

       if (autotrophs(n)%CaCO3_ind.gt.0) then
          lname = trim(autotrophs(n)%lname) //                                  &
               ' CaCO3 Formation Vertical Integral'
          sname = trim(autotrophs(n)%sname) // '_CaCO3_form_zint'
          units = 'mmol/m^3 cm/s'
          vgrid = 'none'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CaCO3_form_zint(n))
       else
          ind%CaCO3_form_zint(n) = -1
       end if
    end do

    lname = 'Total CaCO3 Formation Vertical Integral'
    sname = 'CaCO3_form_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%tot_CaCO3_form_zint)

    ! General 3D diags
    lname = 'Carbonate Ion Concentration'
    sname = 'CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CO3)

    lname = 'Bicarbonate Ion Concentration'
    sname = 'HCO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%HCO3)

    lname = 'Carbonic Acid Concentration'
    sname = 'H2CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%H2CO3)

    lname = 'pH'
    sname = 'pH_3D'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ph_3D)

    lname = 'Carbonate Ion Concentration, Alternative CO2'
    sname = 'CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CO3_ALT_CO2)

    lname = 'Bicarbonate Ion Concentration, Alternative CO2'
    sname = 'HCO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%HCO3_ALT_CO2)

    lname = 'Carbonic Acid Concentration, Alternative CO2'
    sname = 'H2CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%H2CO3_ALT_CO2)

    lname = 'pH, Alternative CO2'
    sname = 'pH_3D_ALT_CO2'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%ph_3D_ALT_CO2)

    lname = 'CO3 concentration at calcite saturation'
    sname = 'co3_sat_calc'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%co3_sat_calc)

    lname = 'CO3 concentration at aragonite saturation'
    sname = 'co3_sat_arag'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%co3_sat_arag)

    lname = 'Nitrification'
    sname = 'NITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%NITRIF)

    lname = 'Denitrification'
    sname = 'DENITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DENITRIF)

    lname = 'O2 Production'
    sname = 'O2_PRODUCTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2_PRODUCTION)

    lname = 'O2 Consumption'
    sname = 'O2_CONSUMPTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%O2_CONSUMPTION)

    lname = 'Apparent O2 Utilization'
    sname = 'AOU'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%AOU)

    lname = 'PAR Average over Model Cell'
    sname = 'PAR_avg'
    units = 'W/m^2'
    vgrid = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PAR_avg)

    lname = 'Total Autotroph Grazing'
    sname = 'graze_auto_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_graze_TOT)

    lname = 'Total C Fixation'
    sname = 'photoC_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_TOT)

    lname = 'Total C Fixation from NO3'
    sname = 'photoC_NO3_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_NO3_TOT)

    lname = 'DOC Production'
    sname = 'DOC_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOC_prod)

    lname = 'DOC Remineralization'
    sname = 'DOC_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOC_remin)

    lname = 'DOCr Remineralization'
    sname = 'DOCr_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOCr_remin)

    lname = 'DON Production'
    sname = 'DON_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DON_prod)

    lname = 'DON Remineralization'
    sname = 'DON_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DON_remin)

    lname = 'DONr Remineralization'
    sname = 'DONr_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DONr_remin)

    lname = 'DOP Production'
    sname = 'DOP_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOP_prod)

    lname = 'DOP Remineralization'
    sname = 'DOP_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOP_remin)

    lname = 'DOPr Remineralization'
    sname = 'DOPr_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOPr_remin)

    lname = 'Iron Scavenging'
    sname = 'Fe_scavenge'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Fe_scavenge)

    lname = 'Iron Scavenging Rate'
    sname = 'Fe_scavenge_rate'
    units = '1/y'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Fe_scavenge_rate)

    ! Particulate 3D diags
    lname = 'POC Flux into Cell'
    sname = 'POC_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%POC_FLUX_IN)

    lname = 'POC Production'
    sname = 'POC_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%POC_PROD)

    lname = 'POC Remineralization'
    sname = 'POC_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%POC_REMIN)

    lname = 'POC Remineralization routed to DIC'
    sname = 'POC_REMIN_DIC'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%POC_REMIN_DIC)

    lname = 'PON Remineralization routed to NH4'
    sname = 'PON_REMIN_NH4'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PON_REMIN_NH4)

    lname = 'POP Remineralization routed to PO4'
    sname = 'POP_REMIN_PO4'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%POP_REMIN_PO4)

    lname = 'CaCO3 Flux into Cell'
    sname = 'CaCO3_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CaCO3_FLUX_IN)

    lname = 'CaCO3 Production'
    sname = 'CaCO3_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CaCO3_PROD)

    lname = 'CaCO3 Remineralization'
    sname = 'CaCO3_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CaCO3_REMIN)

    lname = 'SiO2 Flux into Cell'
    sname = 'SiO2_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SiO2_FLUX_IN)

    lname = 'SiO2 Production'
    sname = 'SiO2_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SiO2_PROD)

    lname = 'SiO2 Remineralization'
    sname = 'SiO2_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SiO2_REMIN)

    lname = 'Dust Flux into Cell'
    sname = 'dust_FLUX_IN'
    units = 'ng/s/m^2'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%dust_FLUX_IN)

    lname = 'Dust Remineralization'
    sname = 'dust_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%dust_REMIN)

    lname = 'P_iron Flux into Cell'
    sname = 'P_iron_FLUX_IN'
    units = 'mmol/m^3 cm/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%P_iron_FLUX_IN)

    lname = 'P_iron Production'
    sname = 'P_iron_PROD'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%P_iron_PROD)

    lname = 'P_iron Remineralization'
    sname = 'P_iron_REMIN'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%P_iron_REMIN)

    ! Autotroph 3D diags
    do n=1,autotroph_cnt

       lname = trim(autotrophs(n)%lname) // ' N Limitation'
       sname = trim(autotrophs(n)%sname) // '_N_lim'
       units = 'none'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%N_lim(n))

       lname = trim(autotrophs(n)%lname) // ' P Limitation'
       sname = trim(autotrophs(n)%sname) // '_P_lim'
       units = 'none'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%P_lim(n))

       lname = trim(autotrophs(n)%lname) // ' Fe Limitation'
       sname = trim(autotrophs(n)%sname) // '_Fe_lim'
       units = 'none'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Fe_lim(n))

       if (autotrophs(n)%kSiO3 > c0) then
          lname = trim(autotrophs(n)%lname) // ' SiO3 Limitation'
          sname = trim(autotrophs(n)%sname) // '_SiO3_lim'
          units = 'none'
          vgrid = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%SiO3_lim(n))
       else
          ind%SiO3_lim(n) = -1
       end if

       lname = trim(autotrophs(n)%lname) // ' Light Limitation'
       sname = trim(autotrophs(n)%sname) // '_light_lim'
       units = 'none'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%light_lim(n))

       lname = trim(autotrophs(n)%lname) // ' C Fixation'
       sname = 'photoC_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC(n))

       lname = trim(autotrophs(n)%lname) // ' C Fixation from NO3'
       sname = 'photoC_NO3_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoC_NO3(n))

       lname = trim(autotrophs(n)%lname) // ' Fe Uptake'
       sname = 'photoFe_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoFe(n))

       lname = trim(autotrophs(n)%lname) // ' NO3 Uptake'
       sname = 'photoNO3_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoNO3(n))

       lname = trim(autotrophs(n)%lname) // ' NH4 Uptake'
       sname = 'photoNH4_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%photoNH4(n))

       lname = trim(autotrophs(n)%lname) // ' DOP Uptake'
       sname = 'DOP_' // trim(autotrophs(n)%sname) // '_uptake'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%DOP_uptake(n))

       lname = trim(autotrophs(n)%lname) // ' PO4 Uptake'
       sname = 'PO4_' // trim(autotrophs(n)%sname) // '_uptake'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%PO4_uptake(n))

       lname = trim(autotrophs(n)%lname) // ' Grazing'
       sname = 'graze_' // trim(autotrophs(n)%sname)
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_graze(n))

       lname = trim(autotrophs(n)%lname) // ' Grazing to POC'
       sname = 'graze_' // trim(autotrophs(n)%sname) // '_poc'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_graze_poc(n))

       lname = trim(autotrophs(n)%lname) // ' Grazing to DOC'
       sname = 'graze_' // trim(autotrophs(n)%sname) // '_doc'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_graze_doc(n))

       lname = trim(autotrophs(n)%lname) // ' Grazing to ZOO'
       sname = 'graze_' // trim(autotrophs(n)%sname) // '_zoo'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_graze_zoo(n))

       lname = trim(autotrophs(n)%lname) // ' Loss'
       sname = trim(autotrophs(n)%sname) // '_loss'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_loss(n))

       lname = trim(autotrophs(n)%lname) // ' Loss to POC'
       sname = trim(autotrophs(n)%sname) // '_loss_poc'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_loss_poc(n))

       lname = trim(autotrophs(n)%lname) // ' Loss to DOC'
       sname = trim(autotrophs(n)%sname) // '_loss_doc'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_loss_doc(n))

       lname = trim(autotrophs(n)%lname) // ' Aggregate'
       sname = trim(autotrophs(n)%sname) // '_agg'
       units = 'mmol/m^3/s'
       vgrid = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%auto_agg(n))

       if (autotrophs(n)%Si_ind.gt.0) then
          lname = trim(autotrophs(n)%lname) // ' Si Uptake'
         !sname = trim(autotrophs(n)%sname) // '_bSi_form'
          sname = trim(autotrophs(n)%sname) // 'bSi_form' ! FIXME: eventually add _
          units = 'mmol/m^3/s'
          vgrid = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%bSi_form(n))
       else
          ind%bSi_form(n) = -1
       end if

       if (autotrophs(n)%CaCO3_ind.gt.0) then
          lname = trim(autotrophs(n)%lname) // ' CaCO3 Formation'
          sname = trim(autotrophs(n)%sname) // '_CaCO3_form'
          units = 'mmol/m^3/s'
          vgrid = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CaCO3_form(n))
       else
          ind%CaCO3_form(n) = -1
       end if

       if (autotrophs(n)%Nfixer) then
          lname = trim(autotrophs(n)%lname) // ' N Fixation'
          sname = trim(autotrophs(n)%sname) // '_Nfix'
          units = 'mmol/m^3/s'
          vgrid = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%Nfix(n))
       else
          ind%Nfix(n) = -1
       end if
    end do

    lname    = 'Total Si Uptake'
    sname    = 'bSi_form'
    units    = 'mmol/m^3/s'
    vgrid    = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%tot_bSi_form)

    lname    = 'Total CaCO3 Formation'
    sname    = 'CaCO3_form'
    units    = 'mmol/m^3/s'
    vgrid    = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%tot_CaCO3_form)

    lname    = 'Total N Fixation'
    sname    = 'Nfix'
    units    = 'mmol/m^3/s'
    vgrid    = 'layer_avg'
    truncate = .true.
    call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%tot_Nfix)

    ! Zooplankton 3D diags
    do n        =1,zooplankton_cnt
       lname    = trim(zooplankton(n)%lname) // ' Loss'
       sname    = trim(zooplankton(n)%sname) // '_loss'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_loss(n))

       lname    = trim(zooplankton(n)%lname) // ' Loss to POC'
       sname    = trim(zooplankton(n)%sname) // '_loss_poc'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_loss_poc(n))

       lname    = trim(zooplankton(n)%lname) // ' Loss to DOC'
       sname    = trim(zooplankton(n)%sname) // '_loss_doc'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_loss_doc(n))

       lname    = trim(zooplankton(n)%lname) // ' grazing loss'
       sname    = 'graze_' // trim(zooplankton(n)%sname)
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_graze(n))

       lname    = trim(zooplankton(n)%lname) // ' grazing loss to POC'
       sname    = 'graze_' // trim(zooplankton(n)%sname) // '_poc'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_graze_poc(n))

       lname    = trim(zooplankton(n)%lname) // ' grazing loss to DOC'
       sname    = 'graze_' // trim(zooplankton(n)%sname) // '_doc'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_graze_doc(n))

       lname    = trim(zooplankton(n)%lname) // ' grazing loss to ZOO'
       sname    = 'graze_' // trim(zooplankton(n)%sname) // '_zoo'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%zoo_graze_zoo(n))

       lname    = trim(zooplankton(n)%lname) // ' grazing gain'
       sname    = 'x_graze_' // trim(zooplankton(n)%sname)
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%x_graze_zoo(n))

    end do

    if (ciso_on) then

       !  nonstandard 3D fields

       lname    = 'PO13C Flux into Cell'
       sname    = 'CISO_PO13C_FLUX_IN'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO13C_FLUX_IN)

       lname    = 'PO13C Production'
       sname    = 'CISO_PO13C_PROD'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO13C_PROD)

       lname    = 'PO13C Remineralization'
       sname    = 'CISO_PO13C_REMIN'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO13C_REMIN)
       
       lname    = 'DO13C Production'
       sname    = 'CISO_DO13C_prod'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DO13C_prod)

       lname    = 'DO13C Remineralization'
       sname    = 'CISO_DO13C_remin'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DO13C_remin)

       lname    = 'Ca13CO3 flux into cell'
       sname    = 'CISO_Ca13CO3_FLUX_IN'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca13CO3_FLUX_IN)

       lname    = 'Ca13CO3 Production'
       sname    = 'CISO_Ca13CO3_PROD'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca13CO3_PROD)

       lname    = 'Ca13CO3 Remineralization'
       sname    = 'CISO_Ca13CO3_REMIN'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca13CO3_REMIN)

       lname    = 'Total 13C Fixation'
       sname    = 'CISO_photo13C_TOT'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo13C_TOT)

       lname    = 'd13C of DIC'
       sname    = 'CISO_DIC_d13C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DIC_d13C)

       lname    = 'd13C of DOC'
       sname    = 'CISO_DOC_d13C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DOC_d13C)

       lname    = 'd13C of zooC'
       sname    = 'CISO_zooC_d13C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_zooC_d13C)

       lname    = 'PO14C Flux into Cell'
       sname    = 'CISO_PO14C_FLUX_IN'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO14C_FLUX_IN)

       lname    = 'PO14C Production'
       sname    = 'CISO_PO14C_PROD'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO14C_PROD)

       lname    = 'PO14C Remineralization'
       sname    = 'CISO_PO14C_REMIN'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_PO14C_REMIN)

       lname    = 'DO14C Production'
       sname    = 'CISO_DO14C_prod'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DO14C_prod)

       lname    = 'DO14C Remineralization'
       sname    = 'CISO_DO14C_remin'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DO14C_remin)

       lname    = 'Ca14CO3 flux into cell'
       sname    = 'CISO_Ca14CO3_FLUX_IN'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca14CO3_FLUX_IN)

       lname    = 'Ca14CO3 Production'
       sname    = 'CISO_Ca14CO3_PROD'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca14CO3_PROD)

       lname    = 'Ca14CO3 Remineralization'
       sname    = 'CISO_Ca14CO3_REMIN'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca14CO3_REMIN)

       lname    = 'Total 14C Fixation'
       sname    = 'CISO_photo14C_TOT'
       units    = 'mmol/m^3/s'
       vgrid    = 'layer_avg'
       truncate = .true.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo14C_TOT)

       lname    = 'd14C of DIC'
       sname    = 'CISO_DIC_d14C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DIC_d14C)

       lname    = 'd14C of DOC'
       sname    = 'CISO_DOC_d14C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_DOC_d14C)

       lname    = 'd14C of zooC'
       sname    = 'CISO_zooC_d14C'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_zooC_d14C)

       !  Nonstandard 2D fields

       lname    = 'Total 13C Fixation Vertical Integral'
       sname    = 'CISO_photo13C_TOT_zint'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo13C_TOT_zint)

       lname    = 'Total 14C Fixation Vertical Integral'
       sname    = 'CISO_photo14C_TOT_zint'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo14C_TOT_zint)

       lname    = '13Ctot Source Sink Term Vertical Integral'
       sname    = 'CISO_Jint_13Ctot'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Jint_13Ctot)

       lname    = '14Ctot Source Sink Term Vertical Integral'
       sname    = 'CISO_Jint_14Ctot'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Jint_14Ctot)

       lname    = '13Ctot Source Sink Term Vertical Integral, 0-100m'
       sname    = 'CISO_Jint_100m_13Ctot'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Jint_100m_13Ctot)

       lname    = '14Ctot Source Sink Term Vertical Integral, 0-100m'
       sname    = 'CISO_Jint_100m_14Ctot'
       units    = 'mmol/m^3 cm/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Jint_100m_14Ctot)

       !  Nonstandard autotroph 2D and 3D fields for each autotroph

       do n = 1, autotroph_cnt

          !FIXME - the comments seem to be needed below - need to fix this
!!$          if (autotrophs(n)%Ca13CO3_ind > 0) then
             lname    = trim(autotrophs(n)%lname) // ' Ca13CO3 Formation'
             sname    = 'CISO_' // trim(autotrophs(n)%sname) // '_Ca13CO3_form'
             units    = 'mmol/m^3/s'
             vgrid    = 'layer_avg'
             truncate = .true.
             call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca13CO3_form(n))

             lname    = trim(autotrophs(n)%lname) // ' Ca13CO3 Formation Vertical Integral'
             sname    = trim(sname) // '_zint'
             units    = 'mmol/m^3 cm/s' 
             vgrid    = 'none'
             truncate = .false.
             call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca13CO3_form_zint(n))
!!$          else
!!$             ind%CISO_Ca13CO3_form(n) = -1
!!$             ind%CISO_Ca13CO3_form_zint(n) = -1
!!$          end if

!!$          if (autotrophs(n)%Ca14CO3_ind > 0) then
             lname    = trim(autotrophs(n)%lname) // ' Ca14CO3 Formation'
             sname    = 'CISO_' // trim(autotrophs(n)%sname) // '_Ca14CO3_form'
             units    = 'mmol/m^3/s'
             vgrid    = 'layer_avg'
             truncate = .true.
             call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca14CO3_form(n))

             lname    = trim(autotrophs(n)%lname) // ' Ca14CO3 Formation Vertical Integral'
             sname    = trim(sname) // '_zint'
             units    = 'mmol/m^3 cm/s' 
             vgrid    = 'none'
             truncate = .false.
             call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_Ca14CO3_form_zint(n))
!!$          else
!!$             ind%CISO_Ca14CO3_form(n) = -1
!!$             ind%CISO_Ca14CO3_form_zint(n) = -1
!!$          endif

          lname    = trim(autotrophs(n)%lname) // ' d13C of CaCO3'
          sname    = 'CISO_autotrophCaCO3_d13C_' // trim(autotrophs(n)%sname)
          units    = 'mmol/m^3/s'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_autotrophCaCO3_d13C(n))

          lname    = trim(autotrophs(n)%lname) // ' d14C of CaCO3'
          sname    = 'CISO_autotrophCaCO3_d14C_' // trim(autotrophs(n)%sname)
          units    = 'mmol/m^3/s'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_autotrophCaCO3_d14C(n))

          lname    = trim(autotrophs(n)%lname) // ' 13C Fixation'
          sname    = 'CISO_photo13C_' // trim(autotrophs(n)%sname)
          units    = 'mmol/m^3/s'
          vgrid    = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo13C(n))

          lname    = trim(autotrophs(n)%lname) // ' 14C Fixation'
          sname    = 'CISO_photo14C_' // trim(autotrophs(n)%sname)
          units    = 'mmol/m^3/s'
          vgrid    = 'layer_avg'
          truncate = .true.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo14C(n))

          lname    = trim(autotrophs(n)%lname) // ' 13C Fixation Vertical Integral'
          sname    = 'CISO_photo13C_' // trim(autotrophs(n)%sname) // '_zint'
          units    = 'mmol/m^3 cm/s'
          vgrid    = 'none'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo13C_zint(n))

          lname    = trim(autotrophs(n)%lname) // ' 14C Fixation Vertical Integral'
          sname    = 'CISO_photo14C_' // trim(autotrophs(n)%sname) // '_zint'
          units    = 'mmol/m^3 cm/s'
          vgrid    = 'none'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_photo14C_zint(n))

          lname    = trim(autotrophs(n)%lname) // ' discrimination factor (eps)'
          sname    = 'CISO_eps_autotroph_' // trim(autotrophs(n)%sname)
          units    = 'permil'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_eps_autotroph(n))

          lname    = trim(autotrophs(n)%lname) // ' d13C'
          sname    = 'CISO_d13C_' // trim(autotrophs(n)%sname)
          units    = 'permil'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_d13C(n))

          lname    = trim(autotrophs(n)%lname) // ' d14C'
          sname    = 'CISO_d14C_' // trim(autotrophs(n)%sname)
          units    = 'permil'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_d14C(n))

          lname    = trim(autotrophs(n)%lname) // ' instanteous growth rate over [CO2*]'
          sname    = 'CISO_mui_to_co2star_' // trim(autotrophs(n)%sname)
          units    = 'm^3/mmol C/s'
          vgrid    = 'layer_avg'
          truncate = .false.
          call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_mui_to_co2star(n))

       end do

       !  More nonstandard 3D fields

       lname    = 'Equilibrium fractionation (CO2_gaseous <-> CO2_aq)'
       sname    = 'CISO_eps_aq_g'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_eps_aq_g)

       lname    = 'Equilibrium fractionation between total DIC and gaseous CO2'
       sname    = 'CISO_eps_dic_g'
       units    = 'permil'
       vgrid    = 'layer_avg'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%CISO_eps_dic_g)

       !  Vars to sum up burial in sediments (2D)

       lname    = 'Ca13CO3 Flux to Sediments'
       sname    = 'calcToSed_13C'
       units    = 'nmolC/cm^2/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%calcToSed_13C)

       lname    = 'PO13C Flux to Sediments'
       sname    = 'pocToSed_13C'
       units    = 'nmolC/cm^2/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pocToSed_13C)

       lname    = 'Ca14CO3 Flux to Sediments'
       sname    = 'calcToSed_14C'
       units    = 'nmolC/cm^2/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%calcToSed_14C)

       lname    = 'PO14C Flux to Sediments'
       sname    = 'pocToSed_14C'
       units    = 'nmolC/cm^2/s'
       vgrid    = 'none'
       truncate = .false.
       call diags%add_diagnostic(lname, sname, units, vgrid, truncate, ind%pocToSed_14C)

    end if  ! end of if ciso_on

    end associate

    !-----------------------------------------------------------------
    ! Restoring diagnostics
    !-----------------------------------------------------------------

    ! Allocate memory for restore diagnostics
    call marbl_restore_diags%construct(max_restore_diags,                     &
         num_elements_interior, num_levels)

    associate(                        &
         diags => marbl_restore_diags &
         )

    do n = 1,ecosys_tracer_cnt
       lname = trim(tracer_d_module(n)%long_name) // " Restoring"
       sname = trim(tracer_d_module(n)%short_name) // "_RESTORE"
       units = 'mmol/m^3'
       vgrid = 'layer_avg'
       
       ! Note that tmp_id is a temp variable because restoring diagnostics
       ! have same indexing as the ecosys tracers
       call diags%add_diagnostic(lname, sname, units, vgrid, .false., tmp_id)
    end do

    end associate

    call marbl_interior_diags%set_to_zero()
    call marbl_restore_diags%set_to_zero()
    call marbl_forcing_diags%set_to_zero()

  end subroutine marbl_diagnostics_init

  !***********************************************************************

  subroutine store_diagnostics_carbonate(marbl_domain, carbonate, marbl_diags)

    type(marbl_column_domain_type)      , intent(in)    :: marbl_domain
    type(carbonate_type)         , intent(in)    :: carbonate(:)
    type(marbl_diagnostics_type) , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: k
    real(r8) :: zsat_calcite, zsat_aragonite
    !-----------------------------------------------------------------------
    
    associate(                                               &
         km                => marbl_domain%km,               &
         diags             => marbl_diags%diags(:),          &
         ind               => marbl_interior_diag_ind,       &
         CO3               => carbonate(:)%CO3,              &
         CO3_sat_calcite   => carbonate(:)%CO3_sat_calcite,  &
         CO3_sat_aragonite => carbonate(:)%CO3_sat_aragonite &
         )

    ! Find depth where CO3 = CO3_sat_calcite or CO3_sat_argonite
    diags(ind%zsatcalc)%field_2d(1) =  marbl_compute_saturation_depth(marbl_domain, CO3, CO3_sat_calcite)
    diags(ind%zsatarag)%field_2d(1) =  marbl_compute_saturation_depth(marbl_domain, CO3, CO3_sat_aragonite)

    do k = 1, km
       diags(ind%CO3)%field_3d(k, 1)           = carbonate(k)%CO3
       diags(ind%HCO3)%field_3d(k, 1)          = carbonate(k)%HCO3
       diags(ind%H2CO3)%field_3d(k, 1)         = carbonate(k)%H2CO3
       diags(ind%pH_3D)%field_3d(k, 1)         = carbonate(k)%pH
       diags(ind%CO3_ALT_CO2)%field_3d(k, 1)   = carbonate(k)%CO3_ALT_CO2
       diags(ind%HCO3_ALT_CO2)%field_3d(k, 1)  = carbonate(k)%HCO3_ALT_CO2
       diags(ind%H2CO3_ALT_CO2)%field_3d(k, 1) = carbonate(k)%H2CO3_ALT_CO2
       diags(ind%pH_3D_ALT_CO2)%field_3d(k, 1) = carbonate(k)%pH_ALT_CO2
       diags(ind%co3_sat_calc)%field_3d(k, 1)  = carbonate(k)%CO3_sat_calcite
       diags(ind%co3_sat_arag)%field_3d(k, 1)  = carbonate(k)%CO3_sat_aragonite
    end do

    end associate

  end subroutine store_diagnostics_carbonate

  !***********************************************************************

  function marbl_compute_saturation_depth(marbl_domain, CO3, sat_val)

    type(marbl_column_domain_type) , intent(in) :: marbl_domain
    real(r8)                , intent(in) :: CO3(:)
    real(r8)                , intent(in) :: sat_val(:)

    real(r8) :: marbl_compute_saturation_depth

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real(r8) :: anomaly(marbl_domain%km) ! CO3 concentration above saturation at level
    integer  :: k
    !-----------------------------------------------------------------------

    associate(                    &
         kmt => marbl_domain%kmt, &
         zt  => marbl_domain%zt,  &
         zw  => marbl_domain%zw   &
         )

    anomaly(:) = CO3(:) - sat_val(:)

    if (all(anomaly(1:kmt).gt.c0)) then
       marbl_compute_saturation_depth = zw(kmt)
    elseif (anomaly(1).le.c0) then
       marbl_compute_saturation_depth = c0
    else
       do k=2,kmt
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

    real(r8)                     , intent(in)    :: nitrif(:)
    real(r8)                     , intent(in)    :: denitrif(:)
    type(marbl_diagnostics_type) , intent(inout) :: marbl_diags

    associate(                            &
         diags => marbl_diags%diags(:),   &
         ind   => marbl_interior_diag_ind &
         )

    diags(ind%NITRIF)%field_3d(:, 1)   = nitrif
    diags(ind%DENITRIF)%field_3d(:, 1) = denitrif

    end associate

  end subroutine store_diagnostics_nitrification

  !***********************************************************************

  subroutine store_diagnostics_autotrophs(marbl_domain, &
       autotroph_secondary_species, marbl_diags)

    type(marbl_column_domain_type)         , intent(in)    :: marbl_domain
    type(autotroph_secondary_species_type) , intent(in)    :: autotroph_secondary_species(:,:) ! autotroph_cnt, km
    type(marbl_diagnostics_type)           , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt, &
         delta_z => marbl_domain%delta_z(:)  &
         )

    diags(ind%tot_CaCO3_form_zint)%field_2d(1) = c0
    diags(ind%tot_bSi_form)%field_3d(:, 1) = c0
    diags(ind%tot_Nfix)%field_3d(:, 1) = c0
    diags(ind%tot_CaCO3_form)%field_3d(:, 1) = c0
    do n = 1, autotroph_cnt
       diags(ind%N_lim(n))%field_3d(:, 1)       = autotroph_secondary_species(n,:)%VNtot
       diags(ind%Fe_lim(n))%field_3d(:, 1)      = autotroph_secondary_species(n,:)%VFe
       diags(ind%P_lim(n))%field_3d(:, 1)       = autotroph_secondary_species(n,:)%VPtot

       if (ind%SiO3_lim(n).ne.-1) then
          diags(ind%SiO3_lim(n))%field_3d(:, 1) = autotroph_secondary_species(n,:)%VSiO3
       end if

       diags(ind%light_lim(n))%field_3d(:, 1)   = autotroph_secondary_species(n,:)%light_lim
       diags(ind%photoNO3(n))%field_3d(:, 1)    = autotroph_secondary_species(n,:)%NO3_V
       diags(ind%photoNH4(n))%field_3d(:, 1)    = autotroph_secondary_species(n,:)%NH4_V
       diags(ind%PO4_uptake(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%PO4_V
       diags(ind%DOP_uptake(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%DOP_V
       diags(ind%photoFE(n))%field_3d(:, 1)     = autotroph_secondary_species(n,:)%photoFe

       if (ind%bSi_form(n).ne.-1) then
          diags(ind%bSi_form(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%photoSi
          diags(ind%tot_bSi_form)%field_3d(:, 1) = diags(ind%tot_bSi_form)%field_3d(:, 1) + &
               diags(ind%bSi_form(n))%field_3d(:, 1)
       endif

       if (ind%CaCO3_form(n).ne.-1) then
          diags(ind%CaCO3_form(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%CaCO3_PROD
          diags(ind%tot_CaCO3_form)%field_3d(:, 1) = diags(ind%tot_CaCO3_form)%field_3d(:, 1) + &
               diags(ind%CaCO3_form(n))%field_3d(:, 1)
       end if

       if (ind%Nfix(n).ne.-1) then
          diags(ind%Nfix(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%Nfix
          diags(ind%tot_Nfix)%field_3d(:, 1) = diags(ind%tot_Nfix)%field_3d(:, 1) + &
               diags(ind%Nfix(n))%field_3d(:, 1)
       end if

       diags(ind%auto_graze(n))%field_3d(:, 1)     = autotroph_secondary_species(n,:)%auto_graze
       diags(ind%auto_graze_poc(n))%field_3d(:, 1) = autotroph_secondary_species(n,:)%auto_graze_poc
       diags(ind%auto_graze_doc(n))%field_3d(:, 1) = autotroph_secondary_species(n,:)%auto_graze_doc
       diags(ind%auto_graze_zoo(n))%field_3d(:, 1) = autotroph_secondary_species(n,:)%auto_graze_zoo
       diags(ind%auto_loss(n))%field_3d(:, 1)      = autotroph_secondary_species(n,:)%auto_loss
       diags(ind%auto_loss_poc(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%auto_loss_poc
       diags(ind%auto_loss_doc(n))%field_3d(:, 1)  = autotroph_secondary_species(n,:)%auto_loss_doc
       diags(ind%auto_agg(n))%field_3d(:, 1)       = autotroph_secondary_species(n,:)%auto_agg
       diags(ind%photoC(n))%field_3d(:, 1)         = autotroph_secondary_species(n,:)%photoC

       diags(ind%photoC_NO3(n))%field_3d(:, 1) = c0
       where (autotroph_secondary_species(n,:)%VNtot > c0)
          diags(ind%photoC_NO3(n))%field_3d(:, 1) = autotroph_secondary_species(n,:)%photoC * &
               (autotroph_secondary_species(n,:)%VNO3 / autotroph_secondary_species(n,:)%VNtot)
       end where

       ! vertical integrals
       if (ind%CaCO3_form_zint(n).ne.-1) then
          call compute_vertical_integrals(autotroph_secondary_species(n,:)%CaCO3_PROD, &
               delta_z, kmt, full_depth_integral=diags(ind%CaCO3_form_zint(n))%field_2d(1))

          diags(ind%tot_CaCO3_form_zint)%field_2d(1) = diags(ind%tot_CaCO3_form_zint)%field_2d(1) + &
               diags(ind%CaCO3_form_zint(n))%field_2d(1)
       end if

       diags(ind%photoC_zint(n))%field_2d(1) = c0

       call compute_vertical_integrals(autotroph_secondary_species(n,:)%photoC, &
            delta_z, kmt, full_depth_integral=diags(ind%photoC_zint(n))%field_2d(1))

       call compute_vertical_integrals(diags(ind%photoC_NO3(n))%field_3d(:, 1), &
            delta_z, kmt, full_depth_integral=diags(ind%photoC_NO3_zint(n))%field_2d(1))
    end do ! do n

    end associate

  end subroutine store_diagnostics_autotrophs

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_autotroph_sums(marbl_domain, &
       autotroph_secondary_species, marbl_diags)

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    type(autotroph_secondary_species_type), dimension(:,:), intent(in) :: autotroph_secondary_species ! autotroph_cnt, km
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: n

    associate(                              &
         ind   => marbl_interior_diag_ind,  &
         diags   => marbl_diags%diags(:),   &
         delta_z => marbl_domain%delta_z(:) &
         )

    diags(ind%auto_graze_TOT)%field_3d(:, 1) = sum(autotroph_secondary_species%auto_graze, dim=1)
    diags(ind%photoC_TOT)%field_3d(:, 1) = sum(autotroph_secondary_species%photoC, dim=1)

    diags(ind%photoC_NO3_TOT)%field_3d(:, 1) = c0
    do n = 1, autotroph_cnt
       where (autotroph_secondary_species(n,:)%VNtot > c0)
          diags(ind%photoC_NO3_TOT)%field_3d(:, 1) = diags(ind%photoC_NO3_TOT)%field_3d(:, 1) + &
               (autotroph_secondary_species(n,:)%VNO3 /  &
               autotroph_secondary_species(n,:)%VNtot) * &
               autotroph_secondary_species(n,:)%photoC
       end where
    end do

    diags(ind%photoC_TOT_zint)%field_2d(1) = sum(delta_z * sum(autotroph_secondary_species%photoC, dim=1))
    diags(ind%photoC_NO3_TOT_zint)%field_2d(1) = sum(delta_z * diags(ind%photoC_NO3_TOT)%field_3d(:, 1))

    end associate

  end subroutine store_diagnostics_autotroph_sums

  !***********************************************************************

  subroutine store_diagnostics_particulates(marbl_domain, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, PON_remin, PON_sed_loss, POP_remin, &
       POP_sed_loss, sed_denitrif, other_remin, marbl_diags)

    !-----------------------------------------------------------------------
    ! - Set tavg variables.
    ! - Accumulte losses of BGC tracers to sediments
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : column_sinking_particle_type
    use marbl_parms     , only : POCremin_refract
    use marbl_parms     , only : PONremin_refract
    use marbl_parms     , only : POPremin_refract

    type(marbl_column_domain_type)     , intent(in)    :: marbl_domain
    type(column_sinking_particle_type) , intent(in)    :: POC
    type(column_sinking_particle_type) , intent(in)    :: P_CaCO3
    type(column_sinking_particle_type) , intent(in)    :: P_SiO2
    type(column_sinking_particle_type) , intent(in)    :: dust
    type(column_sinking_particle_type) , intent(in)    :: P_iron
    real(r8), dimension(:)             , intent(in)    :: PON_remin    ! km
    real(r8), dimension(:)             , intent(in)    :: PON_sed_loss ! km
    real(r8), dimension(:)             , intent(in)    :: POP_remin    ! km
    real(r8), dimension(:)             , intent(in)    :: POP_sed_loss ! km
    real(r8), dimension(:)             , intent(in)    :: sed_denitrif ! km
    real(r8), dimension(:)             , intent(in)    :: other_remin  ! km
    type(marbl_diagnostics_type)       , intent(inout) :: marbl_diags

    associate(                               &
         ind     => marbl_interior_diag_ind, &
         diags   => marbl_diags%diags(:),    &
         delta_z => marbl_domain%delta_z(:)  &
         )

    diags(ind%POC_FLUX_IN)%field_3d(:, 1)    = POC%sflux_in + POC%hflux_in
    diags(ind%POC_PROD)%field_3d(:, 1)       = POC%prod
    diags(ind%POC_REMIN)%field_3d(:, 1)      = POC%remin

    diags(ind%POC_REMIN_DIC)%field_3d(:, 1)  = POC%remin * (c1 - POCremin_refract)
    diags(ind%PON_REMIN_NH4)%field_3d(:, 1)  = PON_remin * (c1 - PONremin_refract)
    diags(ind%POP_REMIN_PO4)%field_3d(:, 1)  = POP_remin * (c1 - POPremin_refract)

    diags(ind%CaCO3_FLUX_IN)%field_3d(:, 1)  = P_CaCO3%sflux_in + P_CaCO3%hflux_in
    diags(ind%CaCO3_PROD)%field_3d(:, 1)     = P_CaCO3%prod
    diags(ind%CaCO3_REMIN)%field_3d(:, 1)    = P_CaCO3%remin

    diags(ind%SiO2_FLUX_IN)%field_3d(:, 1)   = P_SiO2%sflux_in + P_SiO2%hflux_in
    diags(ind%SiO2_PROD)%field_3d(:, 1)      = P_SiO2%prod
    diags(ind%SiO2_REMIN)%field_3d(:, 1)     = P_SiO2%remin

    diags(ind%dust_FLUX_IN)%field_3d(:, 1)   = dust%sflux_in + dust%hflux_in
    diags(ind%dust_REMIN)%field_3d(:, 1)     = P_SiO2%remin

    diags(ind%P_iron_FLUX_IN)%field_3d(:, 1) = P_iron%sflux_in + P_iron%hflux_in
    diags(ind%P_iron_PROD)%field_3d(:, 1)    = P_iron%prod
    diags(ind%P_iron_REMIN)%field_3d(:, 1)   = P_iron%remin

    diags(ind%calcToSed)%field_2d(1)   = sum(P_CaCO3%sed_loss)
    diags(ind%bsiToSed)%field_2d(1)    = sum(P_SiO2%sed_loss)
    diags(ind%pocToSed)%field_2d(1)    = sum(POC%sed_loss)
    diags(ind%SedDenitrif)%field_2d(1) = sum(sed_denitrif * delta_z)
    diags(ind%OtherRemin)%field_2d(1)  = sum(other_remin * delta_z)
    diags(ind%ponToSed)%field_2d(1)    = sum(PON_sed_loss)
    diags(ind%popToSed)%field_2d(1)    = sum(POP_sed_loss)
    diags(ind%dustToSed)%field_2d(1)   = sum(dust%sed_loss)
    diags(ind%pfeToSed)%field_2d(1)    = sum(P_iron%sed_loss)

    end associate

  end subroutine store_diagnostics_particulates

   !***********************************************************************

  subroutine store_diagnostics_oxygen(marbl_domain, marbl_gcm_state, &
       column_o2, o2_production, o2_consumption, marbl_diags)

    use marbl_oxygen, only : o2sat_scalar

    type(marbl_column_domain_type) , intent(in)    :: marbl_domain
    type(marbl_gcm_state_type)     , intent(in)    :: marbl_gcm_state
    real(r8)                       , intent(in)    :: column_o2(:)
    real(r8)                       , intent(in)    :: o2_production(:)
    real(r8)                       , intent(in)    :: o2_consumption(:)
    type(marbl_diagnostics_type)   , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: k, min_ind
    !-----------------------------------------------------------------------

    ! Find min_o2 and depth_min_o2
    associate(                                          &
         temperature => marbl_gcm_state%temperature(:), &
         salinity    => marbl_gcm_state%salinity(:),    &
         kmt         => marbl_domain%kmt,               &
         zt          => marbl_domain%zt,                &
         diags       => marbl_diags%diags(:),           &
         ind         => marbl_interior_diag_ind         &
         )

    min_ind = minloc(column_o2(1:kmt), dim=1)

    diags(ind%O2_ZMIN)%field_2d(1)       = column_o2(min_ind)
    diags(ind%O2_ZMIN_DEPTH)%field_2d(1) = zt(min_ind)

    diags(ind%O2_PRODUCTION)%field_3d(:, 1)  = o2_production
    diags(ind%O2_CONSUMPTION)%field_3d(:, 1) = o2_consumption

    do k=1,kmt
       diags(ind%AOU)%field_3d(k, 1) = O2SAT_scalar(temperature(k), salinity(k)) - column_o2(k)
    end do

    end associate

  end subroutine store_diagnostics_oxygen

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_photosynthetically_available_radiation( marbl_domain, &
       PAR_col_frac, PAR_avg, marbl_diags)

    type(marbl_column_domain_type) , intent(in)    :: marbl_domain
    real(r8)                       , intent(in)    :: PAR_col_frac(:)
    real(r8)                       , intent(in)    :: PAR_avg(:,:)
    type(marbl_diagnostics_type)   , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: k
    !-----------------------------------------------------------------------

    associate(                                      &
         km           => marbl_domain%km,           &
         PAR_nsubcols => marbl_domain%PAR_nsubcols, &
         diags        => marbl_diags%diags(:),      &
         ind          => marbl_interior_diag_ind    &
         )

    do k=1,km
       diags(ind%PAR_avg)%field_3d(k, 1) = sum(PAR_col_frac(:)*PAR_avg(k,:))
    end do

    end associate

  end subroutine store_diagnostics_photosynthetically_available_radiation

  !***********************************************************************

  subroutine store_diagnostics_misc(marbl_diags)
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags
  end subroutine store_diagnostics_misc

  !***********************************************************************

  subroutine store_diagnostics_zooplankton(zooplankton_secondary_species, marbl_diags)

    type(zooplankton_secondary_species_type) , intent(in)    :: zooplankton_secondary_species(:,:)
    type(marbl_diagnostics_type)             , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    !-----------------------------------------------------------------------

    associate(                            &
         diags => marbl_diags%diags(:),   &
         ind   => marbl_interior_diag_ind &
         )

    do n = 1, zooplankton_cnt
       diags(ind%zoo_loss(n))%field_3d(:, 1)      = zooplankton_secondary_species(n,:)%zoo_loss
       diags(ind%zoo_loss_poc(n))%field_3d(:, 1)  = zooplankton_secondary_species(n,:)%zoo_loss_poc
       diags(ind%zoo_loss_doc(n))%field_3d(:, 1)  = zooplankton_secondary_species(n,:)%zoo_loss_doc
       diags(ind%zoo_graze(n))%field_3d(:, 1)     = zooplankton_secondary_species(n,:)%zoo_graze
       diags(ind%zoo_graze_poc(n))%field_3d(:, 1) = zooplankton_secondary_species(n,:)%zoo_graze_poc
       diags(ind%zoo_graze_doc(n))%field_3d(:, 1) = zooplankton_secondary_species(n,:)%zoo_graze_doc
       diags(ind%zoo_graze_zoo(n))%field_3d(:, 1) = zooplankton_secondary_species(n,:)%zoo_graze_zoo
       diags(ind%x_graze_zoo(n))%field_3d(:, 1)   = zooplankton_secondary_species(n,:)%x_graze_zoo
    end do

    end associate

  end subroutine store_diagnostics_zooplankton

  !***********************************************************************

  subroutine store_diagnostics_dissolved_organic_matter(marbl_domain, &
       dissolved_organic_matter, fe_scavenge, fe_scavenge_rate, marbl_diags)

    type(marbl_column_domain_type)      , intent(in)    :: marbl_domain
    type(dissolved_organic_matter_type) , intent(in)    :: dissolved_organic_matter(:) ! (km)
    real(r8)                            , intent(in)    :: fe_scavenge(:)              ! (km)
    real(r8)                            , intent(in)    :: fe_scavenge_rate(:)         ! (km)
    type(marbl_diagnostics_type)        , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: k
    !-----------------------------------------------------------------------

    associate(                            &
         km    => marbl_domain%km,        &
         diags => marbl_diags%diags(:),   &
         ind   => marbl_interior_diag_ind &
         )

    do k = 1, km
       diags(ind%DOC_prod)%field_3d(k, 1)         = dissolved_organic_matter(k)%DOC_prod
       diags(ind%DOC_remin)%field_3d(k, 1)        = dissolved_organic_matter(k)%DOC_remin
       diags(ind%DOCr_remin)%field_3d(k, 1)       = dissolved_organic_matter(k)%DOCr_remin
       diags(ind%DON_prod)%field_3d(k, 1)         = dissolved_organic_matter(k)%DON_prod
       diags(ind%DON_remin)%field_3d(k, 1)        = dissolved_organic_matter(k)%DON_remin
       diags(ind%DONr_remin)%field_3d(k, 1)       = dissolved_organic_matter(k)%DONr_remin
       diags(ind%DOP_prod)%field_3d(k, 1)         = dissolved_organic_matter(k)%DOP_prod
       diags(ind%DOP_remin)%field_3d(k, 1)        = dissolved_organic_matter(k)%DOP_remin
       diags(ind%DOPr_remin)%field_3d(k, 1)       = dissolved_organic_matter(k)%DOPr_remin

       diags(ind%Fe_scavenge)%field_3d(k, 1)      = Fe_scavenge(k)
       diags(ind%Fe_scavenge_rate)%field_3d(k, 1) = Fe_scavenge_rate(k)
    end do

    end associate

  end subroutine store_diagnostics_dissolved_organic_matter

  !***********************************************************************

  subroutine store_diagnostics_carbon_fluxes(marbl_domain, &
       POC, P_CaCO3, dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type

    type(marbl_column_domain_type)     , intent(in)    :: marbl_domain
    type(column_sinking_particle_type) , intent(in)    :: POC
    type(column_sinking_particle_type) , intent(in)    :: P_CaCO3
    real(r8)                           , intent(in)    :: dtracer(:,:) ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type)       , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n, auto_ind
    real(r8), dimension(marbl_domain%km) :: work
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt,        &
         delta_z => marbl_domain%delta_z(:)  &
         )

    ! vertical integrals
    work = dtracer(dic_ind,:) + dtracer(doc_ind,:) +                          &
         dtracer(docr_ind,:) + sum(dtracer(zooplankton(:)%C_ind,:), dim=1) +  &
         sum(dtracer(autotrophs(:)%C_ind,:),dim=1)
    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%CaCO3_ind
       if (n.gt.0) then
          work = work + dtracer(n,:)
       end if
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%Jint_Ctot)%field_2d(1),                &
         near_surface_integral=diags(ind%Jint_100m_Ctot)%field_2d(1),         &
         integrated_terms = POC%sed_loss + P_CaCO3%sed_loss)

    end associate

  end subroutine store_diagnostics_carbon_fluxes

  !***********************************************************************

  subroutine store_diagnostics_nitrogen_fluxes(marbl_domain, &
       PON_sed_loss, denitrif, sed_denitrif, autotroph_secondary_species, dtracer, &
       marbl_diags)

    use marbl_parms, only : Q

    type(marbl_column_domain_type)         , intent(in)    :: marbl_domain
    real(r8)                               , intent(in)    :: PON_sed_loss(:) ! km
    real(r8)                               , intent(in)    :: denitrif(:)     ! km
    real(r8)                               , intent(in)    :: sed_denitrif(:) ! km
    type(autotroph_secondary_species_type) , intent(in)    :: autotroph_secondary_species(:,:)
    real(r8)                               , intent(in)    :: dtracer(:,:)      ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type)           , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    real(r8), dimension(marbl_domain%km) :: work
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt,        &
         delta_z => marbl_domain%delta_z(:)  &
         )

    ! vertical integrals
    work = dtracer(no3_ind,:) + dtracer(nh4_ind,:) +                          &
           dtracer(don_ind,:) + dtracer(donr_ind,:) +                         &
           Q * sum(dtracer(zooplankton(:)%C_ind,:), dim=1) +                  &
           Q * sum(dtracer(autotrophs(:)%C_ind,:), dim=1) +                   &
           denitrif(:) + sed_denitrif(:)
    ! subtract out N fixation
    do n = 1, autotroph_cnt
       if (autotrophs(n)%Nfixer) then
          work = work - autotroph_secondary_species(n,:)%Nfix
       end if
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%Jint_Ntot)%field_2d(1),                &
         near_surface_integral=diags(ind%Jint_100m_Ntot)%field_2d(1),         &
         integrated_terms = PON_sed_loss)

    end associate

  end subroutine store_diagnostics_nitrogen_fluxes

  !***********************************************************************

  subroutine store_diagnostics_phosphorus_fluxes(marbl_domain, &
       POP_sed_loss, dtracer, marbl_diags)

    use marbl_parms,  only : Qp_zoo_pom

    type(marbl_column_domain_type) , intent(in)    :: marbl_domain
    real(r8)                       , intent(in)    :: POP_sed_loss(:) ! km
    real(r8)                       , intent(in)    :: dtracer(:,:)    ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type)   , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    real(r8), dimension(marbl_domain%km) :: work
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt,        &
         delta_z => marbl_domain%delta_z(:)  &
         )

    ! vertical integrals
    work = dtracer(po4_ind,:) + dtracer(dop_ind,:) + dtracer(dopr_ind,:)
    do n = 1, zooplankton_cnt
       work = work + Qp_zoo_pom * dtracer(zooplankton(n)%C_ind,:)
    end do
    do n = 1, autotroph_cnt
       work = work + autotrophs(n)%Qp * dtracer(autotrophs(n)%C_ind,:)
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%Jint_Ptot)%field_2d(1),                &
         near_surface_integral=diags(ind%Jint_100m_Ptot)%field_2d(1),         &
         integrated_terms = POP_sed_loss)

    end associate

  end subroutine store_diagnostics_phosphorus_fluxes

  !***********************************************************************

  subroutine store_diagnostics_silicon_fluxes(marbl_domain, &
       P_SiO2, dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type

    type(marbl_column_domain_type)     , intent(in)    :: marbl_domain
    type(column_sinking_particle_type) , intent(in)    :: P_SiO2
    real(r8)                           , intent(in)    :: dtracer(:,:) ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type)       , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    real(r8), dimension(marbl_domain%km) :: work
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt,        &
         zw      => marbl_domain%zw(:),      &
         delta_z => marbl_domain%delta_z(:)  &
         )

    ! vertical integrals
    work = dtracer(sio3_ind,:)
    do n = 1, autotroph_cnt
       if (autotrophs(n)%Si_ind > 0) then
          work = work + dtracer(autotrophs(n)%Si_ind,:)
       end if
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%Jint_Sitot)%field_2d(1),               &
         near_surface_integral=diags(ind%Jint_100m_Sitot)%field_2d(1),        &
         integrated_terms = P_SiO2%sed_loss)

    end associate

  end subroutine store_diagnostics_silicon_fluxes

  !***********************************************************************

  subroutine store_diagnostics_iron_fluxes(marbl_domain, &
       P_iron, dust, fesedflux, dtracer, marbl_diags)

    use marbl_share_mod , only : column_sinking_particle_type
    use marbl_parms     , only : Qfe_zoo
    use marbl_parms     , only : dust_to_Fe

    type(marbl_column_domain_type)            , intent(in)    :: marbl_domain
    type(column_sinking_particle_type) , intent(in)    :: P_iron
    type(column_sinking_particle_type) , intent(in)    :: dust
    real(r8), dimension(:)             , intent(in)    :: fesedflux  ! km
    real(r8), dimension(:,:)           , intent(in)    :: dtracer ! ecosys_tracer_cnt, km
    type(marbl_diagnostics_type)       , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: n
    real(r8), dimension(marbl_domain%km) :: work
    !-----------------------------------------------------------------------

    associate(                               &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind, &
         kmt     => marbl_domain%kmt,        &
         zw      => marbl_domain%zw(:),      &
         delta_z => marbl_domain%delta_z(:)  &
         )

    ! vertical integrals
    work = dtracer(Fe_ind, :) + sum(dtracer(autotrophs(:)%Fe_ind, :),dim=1) + &
           Qfe_zoo * sum(dtracer(zooplankton(:)%C_ind, :),dim=1) -            &
           dust%remin(:) * dust_to_Fe
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%Jint_Fetot)%field_2d(1),               &
         near_surface_integral=diags(ind%Jint_100m_Fetot)%field_2d(1),        &
         integrated_terms = P_iron%sed_loss - fesedflux)

    end associate

  end subroutine store_diagnostics_iron_fluxes

  !***********************************************************************

  subroutine store_diagnostics_sflux(marbl_forcing_input, &
       marbl_forcing_output, marbl_diags)

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.

    use constants, only : mpercm

    use marbl_interface_types , only : marbl_forcing_input_type
    use marbl_interface_types , only : marbl_forcing_output_type
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : ndep_shr_stream_scale_factor
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
    use marbl_share_mod       , only : iron_flux     
    use marbl_share_mod       , only : din_riv_flux     
    use marbl_share_mod       , only : dip_riv_flux     
    use marbl_share_mod       , only : don_riv_flux     
    use marbl_share_mod       , only : dop_riv_flux     
    use marbl_share_mod       , only : dsi_riv_flux     
    use marbl_share_mod       , only : dfe_riv_flux     
    use marbl_share_mod       , only : dic_riv_flux     
    use marbl_share_mod       , only : alk_riv_flux     
    use marbl_parms           , only : ind_nox_flux
    use marbl_parms           , only : ind_nhy_flux
    use marbl_parms           , only : ind_din_riv_flux
    use marbl_parms           , only : ind_dsi_riv_flux
    use marbl_parms           , only : ind_dfe_riv_flux
    use marbl_parms           , only : ind_dic_riv_flux
    use marbl_parms           , only : ind_alk_riv_flux
    use marbl_parms           , only : po4_ind
    use marbl_parms           , only : no3_ind
    use marbl_parms           , only : sio3_ind
    use marbl_parms           , only : nh4_ind
    use marbl_parms           , only : fe_ind
    use marbl_parms           , only : o2_ind
    use marbl_parms           , only : dic_ind
    use marbl_parms           , only : dic_alt_co2_ind
    use marbl_parms           , only : alk_ind
    use marbl_parms           , only : doc_ind
    use marbl_parms           , only : don_ind
    use marbl_parms           , only : dop_ind
    use marbl_parms           , only : dopr_ind
    use marbl_parms           , only : donr_ind
    use marbl_parms           , only : docr_ind

    type(marbl_forcing_input_type)    , intent(in)    :: marbl_forcing_input
    type(marbl_forcing_output_type)   , intent(inout) :: marbl_forcing_output
    type(marbl_diagnostics_type)      , intent(inout) :: marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_diagnostics:store_diagnostics_sflux'
    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

    associate(                                                                &
         diags           => marbl_diags%diags(:),                             &
         ind             => marbl_forcing_diag_ind,                           &
         xkw             => marbl_forcing_input%xkw,                          &
         xco2            => marbl_forcing_input%xco2,                         &
         xco2_alt_co2    => marbl_forcing_input%xco2_alt_co2,                 &
         ifrac           => marbl_forcing_input%ifrac,                        &
         ap_used         => marbl_forcing_input%atm_press,                    &
         dust_flux_in    => marbl_forcing_input%dust_flux,                    &
         marbl_stf       => marbl_forcing_input%marbl_stf,                    &
         ph_prev         => marbl_forcing_output%ph_prev,                     &
         ph_prev_alt_co2 => marbl_forcing_output%ph_prev_alt_co2,             &
         iron_flux_in    => marbl_forcing_output%iron_flux,                   &
         flux_co2        => marbl_forcing_output%flux_co2,                    &
         flux_alt_co2    => marbl_forcing_output%flux_alt_co2,                &
         co2star         => marbl_forcing_output%co2star,                     &
         dco2star        => marbl_forcing_output%dco2star,                    &
         pco2surf        => marbl_forcing_output%pco2surf,                    &
         dpco2           => marbl_forcing_output%dpco2,                       &
         co2star_alt     => marbl_forcing_output%co2star_alt,                 &
         dco2star_alt    => marbl_forcing_output%dco2star_alt,                &
         pco2surf_alt    => marbl_forcing_output%pco2surf_alt,                &
         dpco2_alt       => marbl_forcing_output%dpco2_alt,                   &
         pv_co2          => marbl_forcing_output%pv_co2,                      &
         pv_o2           => marbl_forcing_output%pv_o2,                       &
         schmidt_co2     => marbl_forcing_output%schmidt_co2,                 &
         schmidt_o2      => marbl_forcing_output%schmidt_o2,                  &
         o2sat           => marbl_forcing_output%o2sat,                       &
         stf_module      => marbl_forcing_output%stf_module                   &
         )

    if (lflux_gas_o2 .or. lflux_gas_co2) then

       diags(ind%ECOSYS_IFRAC)%field_2d(:)     = ifrac(:)
       diags(ind%ECOSYS_XKW)%field_2d(:)       = xkw(:)
       diags(ind%ECOSYS_ATM_PRESS)%field_2d(:) = AP_USED(:)

    endif  ! lflux_gas_o2 .or. lflux_gas_co2

    if (lflux_gas_o2) then

       diags(ind%PV_O2)%field_2d(:)      = PV_O2(:)
       diags(ind%SCHMIDT_O2)%field_2d(:) = SCHMIDT_O2(:)
       diags(ind%O2SAT)%field_2d(:)      = O2SAT(:)
       
    endif  ! lflux_gas_o2

    !-----------------------------------------------------------------------
    !  compute CO2 flux, computing disequilibrium one row at a time
    !-----------------------------------------------------------------------
    
    if (lflux_gas_co2) then
       
       diags(ind%CO2STAR)%field_2d(:)              = CO2STAR(:)
       diags(ind%DCO2STAR)%field_2d(:)             = DCO2STAR(:)
       diags(ind%pCO2SURF)%field_2d(:)             = pCO2SURF(:)
       diags(ind%DpCO2)%field_2d(:)                = DpCO2(:)
       
       diags(ind%CO2STAR_ALT_CO2)%field_2d(:)      = CO2STAR_ALT(:)
       diags(ind%DCO2STAR_ALT_CO2)%field_2d(:)     = DCO2STAR_ALT(:)
       diags(ind%pCO2SURF_ALT_CO2)%field_2d(:)     = pCO2SURF_ALT(:)
       diags(ind%DpCO2_ALT_CO2)%field_2d(:)        = DpCO2_ALT(:)
       
       diags(ind%PV_CO2)%field_2d(:)               = PV_CO2(:)
       diags(ind%SCHMIDT_CO2)%field_2d(:)          = SCHMIDT_CO2(:)
       diags(ind%DIC_GAS_FLUX)%field_2d(:)         = FLUX_CO2(:)
       diags(ind%PH)%field_2d(:)                   = PH_PREV(:)
       diags(ind%ATM_CO2)%field_2d(:)              = XCO2(:)
      
       diags(ind%DIC_GAS_FLUX_ALT_CO2)%field_2d(:) = FLUX_ALT_CO2(:)
       diags(ind%PH_ALT_CO2)%field_2d(:)           = PH_PREV_ALT_CO2(:)
       diags(ind%ATM_ALT_CO2)%field_2d(:)          = XCO2_ALT_CO2(:)
       
    endif  !  lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    ! multiply IRON flux by mpercm (.01) to convert from model units (cm/s)(mmol/m^3) to mmol/s/m^2

    if (iron_flux%has_data) then
       diags(ind%IRON_FLUX)%field_2d(:) = IRON_FLUX_IN(:) * mpercm
    endif

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------

    if (nox_flux_monthly%has_data) then
       diags(ind%NOx_FLUX)%field_2d(:) = MARBL_STF(:, ind_nox_flux)
    endif

    if (nhy_flux_monthly%has_data) then
       diags(ind%NHy_FLUX)%field_2d(:) = MARBL_STF(:, ind_nhy_flux)
    endif

    if (trim(ndep_data_type) == 'shr_stream') then
       diags(ind%NOx_FLUX)%field_2d(:) = &
            ndep_shr_stream_scale_factor * MARBL_STF(:, ind_nox_flux)
       diags(ind%NHy_FLUX)%field_2d(:) = &
            ndep_shr_stream_scale_factor * MARBL_STF(:, ind_nhy_flux)
    endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------

    if (din_riv_flux%has_data) then
       diags(ind%DIN_RIV_FLUX)%field_2d(:) = MARBL_STF(:, ind_din_riv_flux)
    endif
    if (dsi_riv_flux%has_data) then
       diags(ind%DSI_RIV_FLUX)%field_2d(:) = MARBL_STF(:, ind_dsi_riv_flux)
    endif
    if (dfe_riv_flux%has_data) then
       diags(ind%DFE_RIV_FLUX)%field_2d(:) = MARBL_STF(:, ind_dfe_riv_flux)
    endif
    if (dic_riv_flux%has_data) then
       diags(ind%DIC_RIV_FLUX)%field_2d(:) = MARBL_STF(:, ind_dic_riv_flux)
    endif
    if (alk_riv_flux%has_data) then
       diags(ind%ALK_RIV_FLUX)%field_2d(:) = MARBL_STF(:, ind_alk_riv_flux)
    endif
    diags(ind%O2_GAS_FLUX)%field_2d(:)   = STF_MODULE(:, o2_ind)
    diags(ind%DIP_RIV_FLUX)%field_2d(:)  = STF_MODULE(:, po4_ind)
    diags(ind%DON_RIV_FLUX)%field_2d(:)  = STF_MODULE(:, don_ind)
    diags(ind%DONr_RIV_FLUX)%field_2d(:) = STF_MODULE(:, donr_ind)
    diags(ind%DOP_RIV_FLUX)%field_2d(:)  = STF_MODULE(:, dop_ind)
    diags(ind%DOPr_RIV_FLUX)%field_2d(:) = STF_MODULE(:, dopr_ind)
    diags(ind%DOC_RIV_FLUX)%field_2d(:)  = STF_MODULE(:, doc_ind)
    diags(ind%DOCr_RIV_FLUX)%field_2d(:) = STF_MODULE(:, docr_ind)

    ! multiply DUST flux by mpercm (.01) to convert from model units (cm/s)(mmol/m^3) to mmol/s/m^2
    diags(ind%DUST_FLUX)%field_2d(:) = DUST_FLUX_IN(:)*mpercm

    end associate

  end subroutine store_diagnostics_sflux

  !*****************************************************************************

  subroutine store_diagnostics_ciso_interior(&
       marbl_domain,        &
       autotroph_d13C,      &
       autotroph_d14C,      &
       autotrophCaCO3_d13C, &
       autotrophCaCO3_d14C, &
       DIC_d13C,            &
       DIC_d14C,            &
       DOC_d13C,            &
       DOC_d14C,            &
       zooC_d13C,           &
       zooC_d14C,           &
       photo13C,            &
       photo14C,            &
       eps_autotroph,       &
       mui_to_co2star,      &
       Ca13CO3_prod,        &
       Ca14CO3_prod,        &
       DO13C_prod,          &
       DO14C_prod,          &
       DO13C_remin,         &
       DO14C_remin,         &
       eps_aq_g,            &
       eps_dic_g,           &
       PO13C,               &
       PO14C,               &
       P_Ca13CO3,           &
       P_Ca14CO3,           &
       dtracer,             &
       marbl_diags)

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Update marbl_interior_ciso_diags data type 
    !---------------------------------------------------------------------

    use marbl_interface_types , only : marbl_diagnostics_type
    use marbl_interface_types , only : marbl_column_domain_type
    use marbl_share_mod       , only : column_sinking_particle_type
    use marbl_parms           , only : di13c_ind
    use marbl_parms           , only : do13c_ind
    use marbl_parms           , only : zoo13C_ind
    use marbl_parms           , only : di14c_ind
    use marbl_parms           , only : do14c_ind
    use marbl_parms           , only : zoo14C_ind

    implicit none

    type(marbl_column_domain_type), intent(in)    :: marbl_domain

    real (r8), intent(in),  dimension(autotroph_cnt, marbl_domain%km) :: &
         autotroph_d13C      , & ! d13C of autotroph C
         autotroph_d14C      , & ! d14C of autotroph C
         autotrophCaCO3_d13C , & ! d13C of autotrophCaCO3
         autotrophCaCO3_d14C , & ! d14C of autotrophCaCO3
         photo13C            , & ! Carbon autotroph 13C-fixation (mmol C/m^3/sec)
         photo14C            , & ! Carbon autotroph 14C-fixation (mmol C/m^3/sec)
         eps_autotroph       , & ! Permil fractionation (or discrimination factor) for Carbon autotroph types sp, diat, diaz
         mui_to_co2star      , & ! Carbon autotroph instanteous growth rate over [CO2*] (m^3 /mmol C /s)
         Ca13CO3_prod        , & ! prod. of 13C CaCO3 by small phyto (mmol CaCO3/m^3/sec)
         Ca14CO3_prod            ! prod. of 13C CaCO3 by small phyto (mmol CaCO3/m^3/sec)

    real (r8), intent(in),  dimension(marbl_domain%km) :: &
         DIC_d13C    , & ! d13C of DIC
         DOC_d13C    , & ! d13C of DOC
         zooC_d13C   , & ! d13C of zooC
         DIC_d14C    , & ! d14C of DIC
         DOC_d14C    , & ! d14C of DOC
         zooC_d14C   , & ! d14C of zooC
         DO13C_prod  , & ! production of 13C DOC (mmol C/m^3/sec)
         DO13C_remin , & ! remineralization of 13C DOC (mmol C/m^3/sec)
         DO14C_prod  , & ! production of 13C DOC (mmol C/m^3/sec)
         DO14C_remin , & ! remineralization of 13C DOC (mmol C/m^3/sec)
         eps_aq_g    , & ! equilibrium fractionation (CO2_gaseous <-> CO2_aq)
         eps_dic_g       ! equilibrium fractionation between total DIC and gaseous CO2

    real (r8), intent(in), dimension(ecosys_ciso_tracer_cnt, marbl_domain%km) :: &
         dtracer   ! computed source/sink terms

    type(column_sinking_particle_type), intent(in) :: &
         PO13C,        &  ! base units = nmol 13C
         PO14C,        &  ! base units = nmol 14C
         P_Ca13CO3,    &  ! base units = nmol 13C CaCO3
         P_Ca14CO3        ! base units = nmol 14C CaCO3

    type(marbl_diagnostics_type), intent(inout) :: &
         marbl_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: k, n, auto_ind
    real (r8)          :: work(marbl_domain%km)
    !-----------------------------------------------------------------------

    associate( &
         km      => marbl_domain%km,         &
         kmt     => marbl_domain%kmt,        &
         zw      => marbl_domain%zw(:),      &
         delta_z => marbl_domain%delta_z(:), &
         diags   => marbl_diags%diags(:),    &
         ind     => marbl_interior_diag_ind  &
         )

    diags(ind%calcToSed_13C)%field_2d(1) = sum(P_Ca13CO3%sed_loss)
    diags(ind%calcToSed_14C)%field_2d(1) = sum(P_Ca14CO3%sed_loss)

    diags(ind%pocToSed_13C)%field_2d(1)  = sum(PO13C%sed_loss)
    diags(ind%pocToSed_14C)%field_2d(1)  = sum(PO14C%sed_loss)

    diags(ind%CISO_photo13C_TOT)%field_3d(:, 1) = sum(photo13C, dim=1)
    diags(ind%CISO_photo14C_TOT)%field_3d(:, 1) = sum(photo14C, dim=1)

    diags(ind%CISO_photo13C_TOT_zint)%field_2d(1) = sum(delta_z * sum(photo13C, dim=1))
    diags(ind%CISO_photo14C_TOT_zint)%field_2d(1) = sum(delta_z * sum(photo14C, dim=1))

    ! Vertical integrals - CISO_Jint_13Ctot and Jint_100m_13Ctot

    diags(ind%CISO_Jint_13Ctot)%field_3d(:, 1) = c0
    work(:) = dtracer(di13c_ind,:) + dtracer(do13c_ind,:) + dtracer(zoo13C_ind,:) &
         + sum(dtracer(autotrophs(:)%C13_ind,:), dim=1)
    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%Ca13CO3_ind
       if (n > 0) then
          work = work + dtracer(n,:)
       end if
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%CISO_Jint_13Ctot)%field_2d(1),         &
         near_surface_integral=diags(ind%CISO_Jint_100m_13Ctot)%field_2d(1),  &
         integrated_terms = PO13C%sed_loss + P_Ca13CO3%sed_loss)

    ! Vertical integral - CISO_Jint_14Ctot and Jint_100m_14Ctot

    diags(ind%CISO_Jint_14Ctot)%field_3d(:, 1) = c0
    work(:) = dtracer(di14c_ind,:) + dtracer(do14c_ind,:) + dtracer(zoo14C_ind,:) &
         + sum(dtracer(autotrophs(:)%C14_ind,:), dim=1)
    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%Ca14CO3_ind
       if (n > 0) then
          work = work + dtracer(n,:)
       end if
    end do
    call compute_vertical_integrals(work, delta_z, kmt,                       &
         full_depth_integral=diags(ind%CISO_Jint_14Ctot)%field_2d(1),         &
         near_surface_integral=diags(ind%CISO_Jint_100m_14Ctot)%field_2d(1),  &
         integrated_terms = PO14C%sed_loss + P_Ca14CO3%sed_loss)

    ! Other vertical integrals

    do n = 1,autotroph_cnt
       diags(ind%CISO_photo13C_zint(n))%field_2d(1) = c0
       call compute_vertical_integrals(photo13C(n,:), delta_z, kmt,           &
            full_depth_integral=diags(ind%CISO_photo13C_zint(n))%field_2d(1))
       
       diags(ind%CISO_photo14C_zint(n))%field_2d(1) = c0
       call compute_vertical_integrals(photo14C(n,:), delta_z, kmt,           &
            full_depth_integral=diags(ind%CISO_photo14C_zint(n))%field_2d(1))
       
       diags(ind%CISO_Ca13CO3_form_zint(n))%field_2d(1) = c0
       call compute_vertical_integrals(Ca13CO3_prod(n,:), delta_z, kmt,       &
            full_depth_integral=diags(ind%CISO_Ca13CO3_form_zint(n))%field_2d(1))
       
       diags(ind%CISO_Ca14CO3_form_zint(n))%field_2d(1) = c0
       call compute_vertical_integrals(Ca14CO3_prod(n,:), delta_z, kmt,       &
            full_depth_integral=diags(ind%CISO_Ca14CO3_form_zint(n))%field_2d(1))
    end do

    do k = 1,km
       do n = 1, autotroph_cnt
          diags(ind%CISO_d13C(n))%field_3d(k, 1)                = autotroph_d13C(n,k)
          diags(ind%CISO_d14C(n))%field_3d(k, 1)                = autotroph_d14C(n,k)

          diags(ind%CISO_autotrophCaCO3_d13C(n))%field_3d(k, 1) = autotrophCaCO3_d13C(n,k)
          diags(ind%CISO_autotrophCaCO3_d14C(n))%field_3d(k, 1) = autotrophCaCO3_d14C(n,k)

          diags(ind%CISO_photo13C(n))%field_3d(k, 1)            = photo13C(n,k)
          diags(ind%CISO_photo14C(n))%field_3d(k, 1)            = photo14C(n,k)

          diags(ind%CISO_eps_autotroph(n))%field_3d(k, 1)       = eps_autotroph(n,k)

          diags(ind%CISO_mui_to_co2star(n))%field_3d(k, 1)      = mui_to_co2star(n,k)

          if (autotrophs(n)%imp_calcifier) then
             diags(ind%CISO_Ca13CO3_form(n))%field_3d(k, 1)     = Ca13CO3_prod(n,k)
             diags(ind%CISO_Ca14CO3_form(n))%field_3d(k, 1)     = Ca14CO3_prod(n,k)
          end if
       end do  ! end loop over autotrophs
    end do  ! end loop over k
    
    do k = 1,km
       diags(ind%CISO_DIC_d13C)%field_3d(k, 1)        = DIC_d13C(k)
       diags(ind%CISO_DIC_d14C)%field_3d(k, 1)        = DIC_d14C(k)

       diags(ind%CISO_DOC_d13C)%field_3d(k, 1)        = DOC_d13C(k)
       diags(ind%CISO_DOC_d14C)%field_3d(k, 1)        = DOC_d14C(k)

       diags(ind%CISO_DO13C_prod)%field_3d(k, 1)      = DO13C_prod(k)   
       diags(ind%CISO_DO14C_prod)%field_3d(k, 1)      = DO14C_prod(k)      

       diags(ind%CISO_DO13C_remin)%field_3d(k, 1)     = DO13C_remin(k)     
       diags(ind%CISO_DO14C_remin)%field_3d(k, 1)     = DO14C_remin(k)     

       diags(ind%CISO_zooC_d13C)%field_3d(k, 1)       = zooC_d13C(k)
       diags(ind%CISO_zooC_d14C)%field_3d(k, 1)       = zooC_d14C(k)

       diags(ind%CISO_eps_aq_g)%field_3d(k, 1)        = eps_aq_g(k)        
       diags(ind%CISO_eps_dic_g)%field_3d(k, 1)       = eps_dic_g(k)       

       diags(ind%CISO_Ca13CO3_flux_in)%field_3d(k, 1) = P_Ca13CO3%sflux_in(k) + P_Ca13CO3%hflux_in(k)
       diags(ind%CISO_Ca14CO3_flux_in)%field_3d(k, 1) = P_Ca14CO3%sflux_in(k) + P_Ca14CO3%hflux_in(k)

       diags(ind%CISO_Ca13CO3_prod)%field_3d(k, 1)    = P_Ca13CO3%prod(k)
       diags(ind%CISO_Ca14CO3_prod)%field_3d(k, 1)    = P_Ca14CO3%prod(k)

       diags(ind%CISO_Ca13CO3_remin)%field_3d(k, 1)   = P_Ca13CO3%remin(k)
       diags(ind%CISO_Ca14CO3_remin)%field_3d(k, 1)   = P_Ca14CO3%remin(k)

       diags(ind%CISO_PO13C_flux_in)%field_3d(k, 1)   = PO13C%sflux_in(k) + PO13C%hflux_in(k)
       diags(ind%CISO_PO14C_flux_in)%field_3d(k, 1)   = PO14C%sflux_in(k) + PO14C%hflux_in(k)

       diags(ind%CISO_PO13C_prod)%field_3d(k, 1)      = PO13C%prod(k)
       diags(ind%CISO_PO14C_prod)%field_3d(k, 1)      = PO14C%prod(k)

       diags(ind%CISO_PO13C_remin)%field_3d(k, 1)     = PO13C%remin(k)
       diags(ind%CISO_PO14C_remin)%field_3d(k, 1)     = PO14C%remin(k)
    end do

    end associate

  end subroutine store_diagnostics_ciso_interior

  !*****************************************************************************

  subroutine compute_vertical_integrals(integrand, delta_z, kmt, &
       full_depth_integral, near_surface_integral, integrated_terms)

    real(kind=r8) , intent(in)             :: integrand(:)
    real(kind=r8) , intent(in)             :: delta_z(:)
    integer       , intent(in)             :: kmt
    ! For some vertical integral diagnostics, we need to add additional terms
    ! that have already been integrated, so they are separated from the
    ! integrand
    real(kind=r8) , intent(in)  , optional :: integrated_terms(:)
    real(kind=r8) , intent(out) , optional :: full_depth_integral
    real(kind=r8) , intent(out) , optional :: near_surface_integral

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: k
    real(kind=r8)      :: integrated_terms_used(size(integrand))
    real(kind=r8)      :: zw
    real(kind=r8)      :: ztop
    real(kind=r8)      :: shallow_depth = 100.0e2_r8
    !-----------------------------------------------------------------------

    if (present(integrated_terms)) then
       integrated_terms_used = integrated_terms
    else
       integrated_terms_used = c0
    end if

    if (present(full_depth_integral)) then
       full_depth_integral = sum(delta_z(1:kmt)*integrand(1:kmt) + integrated_terms_used(1:kmt))
    end if

    if (present(near_surface_integral)) then
       ! initialize integral to zero
       near_surface_integral = c0
       ztop = c0
       zw = c0
       do k=1,kmt
          zw = zw + delta_z(k)
          near_surface_integral = near_surface_integral +                     &
                              min(shallow_depth-ztop,delta_z(k))*integrand(k)
          if (zw.le.shallow_depth) then
             near_surface_integral = near_surface_integral + integrated_terms_used(k)
          else
             exit
          end if
          ztop = zw
       end do
    end if

  end subroutine compute_vertical_integrals

end module ecosys_diagnostics_mod
