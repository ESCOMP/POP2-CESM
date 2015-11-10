module ecosys_diagnostics_mod

  use ecosys_constants, only : ecosys_tracer_cnt

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

  use marbl_interface_types, only : carbonate_type
  use marbl_interface_types, only : zooplankton_secondary_species_type
  use marbl_interface_types, only : autotroph_secondary_species_type
  use marbl_interface_types, only : photosynthetically_available_radiation_type
  use marbl_interface_types, only : dissolved_organic_matter_type
  use marbl_interface_types, only : marbl_diagnostics_type
  use marbl_interface_types, only : marbl_column_domain_type
  use marbl_interface_types, only : diag_cnt

  use grid, only : partial_bottom_cells
  use domain_size, only : km

  Implicit None
  Public

  ! FIXME: MNL - These are copied from ecosys_mod and need to be put somewhere
  ! accessible from both modules
  integer (int_kind), private, parameter :: &
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

!-----------------------------------------------------------------------
!  indices for diagnostic values written to tavg files
!-----------------------------------------------------------------------

  type, public :: marbl_diagnostics_indexing_type
    ! Total number of diagnostics (used for looping)
    integer(int_kind) :: count

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
    integer(int_kind) :: Fe_scavenge
    integer(int_kind) :: Fe_scavenge_rate
  end type marbl_diagnostics_indexing_type

  type(marbl_diagnostics_indexing_type) :: marbl_diag_ind

  integer(int_kind), parameter ::   auto_diag_cnt_2d = 3
  integer (int_kind), parameter ::  &
      photoC_zint_diag_ind     =  1,  &
      photoC_NO3_zint_diag_ind =  2,  &
      CaCO3_form_zint_diag_ind =  3
  integer(int_kind), parameter ::   auto_diag_cnt_3d = 23
  integer (int_kind), parameter ::  &
      N_lim_diag_ind           =  1,  &
      P_lim_diag_ind           =  2,  &
      Fe_lim_diag_ind          =  3,  &
      SiO3_lim_diag_ind        =  4,  &
      light_lim_diag_ind       =  5,  &
      photoC_diag_ind          =  6,  &
      photoC_NO3_diag_ind      =  7,  &
      photoFe_diag_ind         =  8,  &
      photoNO3_diag_ind        =  9,  &
      photoNH4_diag_ind        = 10,  &
      DOP_uptake_diag_ind      = 11,  &
      PO4_uptake_diag_ind      = 12,  &
      auto_graze_diag_ind      = 13,  &
      auto_graze_poc_diag_ind  = 14,  &
      auto_graze_doc_diag_ind  = 15,  &
      auto_graze_zoo_diag_ind  = 16,  &
      auto_loss_diag_ind       = 17,  &
      auto_loss_poc_diag_ind   = 18,  &
      auto_loss_doc_diag_ind   = 19,  &
      auto_agg_diag_ind        = 20,  &
      bSi_form_diag_ind        = 21,  &
      CaCO3_form_diag_ind      = 22,  &
      Nfix_diag_ind            = 23

  integer(int_kind), parameter ::    zoo_diag_cnt_2d =  0
  integer(int_kind), parameter ::    zoo_diag_cnt_3d =  8
  integer (int_kind), parameter ::   &
      zoo_loss_diag_ind        =  1, &
      zoo_loss_poc_diag_ind    =  2, &
      zoo_loss_doc_diag_ind    =  3, &
      zoo_graze_diag_ind       =  4, &
      zoo_graze_poc_diag_ind   =  5, &
      zoo_graze_doc_diag_ind   =  6, &
      zoo_graze_zoo_diag_ind   =  7, &
      x_graze_zoo_diag_ind     =  8

  integer(int_kind), parameter ::   part_diag_cnt_2d =   9
  integer (int_kind), parameter ::   &
      calcToSed_diag_ind       = 1, &
      pocToSed_diag_ind        = 2, &
      ponToSed_diag_ind        = 3, &
      SedDenitrif_diag_ind     = 4, &
      OtherRemin_diag_ind      = 5, &
      popToSed_diag_ind        = 6, &
      bsiToSed_diag_ind        = 7, &
      dustToSed_diag_ind       = 8, &
      pfeToSed_diag_ind        = 9

  integer(int_kind), parameter ::   part_diag_cnt_3d =  14
  integer (int_kind), parameter ::   &
      POC_FLUX_IN_diag_ind     =  1, &
      POC_PROD_diag_ind        =  2, &
      POC_REMIN_diag_ind       =  3, &
      CaCO3_FLUX_IN_diag_ind   =  4, &
      CaCO3_PROD_diag_ind      =  5, &
      CaCO3_REMIN_diag_ind     =  6, &
      SiO2_FLUX_IN_diag_ind    =  7, &
      SiO2_PROD_diag_ind       =  8, &
      SiO2_REMIN_diag_ind      =  9, &
      dust_FLUX_IN_diag_ind    = 10, &
      dust_REMIN_diag_ind      = 11, &
      P_iron_FLUX_IN_diag_ind  = 12, &
      P_iron_PROD_diag_ind     = 13, &
      P_iron_REMIN_diag_ind    = 14

  integer(int_kind), parameter ::   forcing_diag_cnt =  38
  integer (int_kind), parameter ::   &
      ECOSYS_IFRAC_diag_ind         =  1, &
      ECOSYS_XKW_diag_ind           =  2, &
      ECOSYS_ATM_PRESS_diag_ind     =  3, &
      PV_O2_diag_ind                =  4, &
      SCHMIDT_O2_diag_ind           =  5, &
      O2SAT_diag_ind                =  6, &
      O2_GAS_FLUX_diag_ind          =  7, &
      CO2STAR_diag_ind              =  8, &
      DCO2STAR_diag_ind             =  9, &
      pCO2SURF_diag_ind             = 10, &
      DpCO2_diag_ind                = 11, &
      PV_CO2_diag_ind               = 12, &
      SCHMIDT_CO2_diag_ind          = 13, &
      DIC_GAS_FLUX_diag_ind         = 14, &
      PH_diag_ind                   = 15, &
      ATM_CO2_diag_ind              = 16, &
      CO2STAR_ALT_CO2_diag_ind      = 17, &
      DCO2STAR_ALT_CO2_diag_ind     = 18, &
      pCO2SURF_ALT_CO2_diag_ind     = 19, &
      DpCO2_ALT_CO2_diag_ind        = 20, &
      DIC_GAS_FLUX_ALT_CO2_diag_ind = 21, &
      PH_ALT_CO2_diag_ind           = 22, &
      ATM_ALT_CO2_diag_ind          = 23, &
      IRON_FLUX_diag_ind            = 24, &
      DUST_FLUX_diag_ind            = 25, &
      NOx_FLUX_diag_ind             = 26, &
      NHy_FLUX_diag_ind             = 27, &
      DIN_RIV_FLUX_diag_ind         = 28, &
      DIP_RIV_FLUX_diag_ind         = 29, &
      DoN_RIV_FLUX_diag_ind         = 30, &
      DoNr_RIV_FLUX_diag_ind        = 31, &
      DOP_RIV_FLUX_diag_ind         = 32, &
      DOPr_RIV_FLUX_diag_ind        = 33, &
      DSI_RIV_FLUX_diag_ind         = 34, &
      DFE_RIV_FLUX_diag_ind         = 35, &
      DIC_RIV_FLUX_diag_ind         = 36, &
      ALK_RIV_FLUX_diag_ind         = 37, &
      DOC_RIV_FLUX_diag_ind         = 38

contains

  subroutine ecosys_diagnostics_init(marbl_diags)

    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer :: n, auto_ind, zoo_ind
    character(len=char_len) :: lname, sname, units, vgrid
    logical :: truncate

    ! Allocate memory in marbl_diagnostics_type
    call marbl_diags%construct(auto_diag_cnt_2d, auto_diag_cnt_3d,            &
                               zoo_diag_cnt_2d, zoo_diag_cnt_3d,              &
                               part_diag_cnt_2d, part_diag_cnt_3d,            &
                               ecosys_tracer_cnt, autotroph_cnt,              &
                               zooplankton_cnt)

    ! Setup meta-data
    associate(                                                                &
              auto_diags_2d => marbl_diags%auto_diags_2d(:,:),                &
              auto_diags_3d => marbl_diags%auto_diags_3d(:,:),                &
              zoo_diags_2d  => marbl_diags%zoo_diags_2d(:,:),                 &
              zoo_diags_3d  => marbl_diags%zoo_diags_3d(:,:),                 &
              part_diags_2d => marbl_diags%part_diags_2d(:),                  &
              part_diags_3d => marbl_diags%part_diags_3d(:),                  &
              restore_diags => marbl_diags%restore_diags(:)                   &
             )

    lname = 'Calcite Saturation Depth'
    sname = 'zsatcalc'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%zsatcalc)

    lname = 'Aragonite Saturation Depth'
    sname = 'zsatarag'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%zsatarag)

    lname = 'Vertical Minimum of O2'
    sname = 'O2_ZMIN'
    units = 'mmol/m^3'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%O2_ZMIN)

    lname = 'Depth of Vertical Minimum of O2'
    sname = 'O2_ZMIN_DEPTH'
    units = 'cm'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%O2_ZMIN_DEPTH)

    lname = 'Total C Fixation Vertical Integral'
    sname = 'photoC_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%photoC_TOT_zint)

    lname = 'Total C Fixation from NO3 Vertical Integral'
    sname = 'photoC_NO3_TOT_zint'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%photoC_NO3_TOT_zint)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot'
    sname = 'Jint_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ctot, 0-100m'
    sname = 'Jint_100m_Ctot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_100m_Ctot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot'
    sname = 'Jint_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ntot, 0-100m'
    sname = 'Jint_100m_Ntot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_100m_Ntot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot'
    sname = 'Jint_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Ptot, 0-100m'
    sname = 'Jint_100m_Ptot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_100m_Ptot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot'
    sname = 'Jint_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_Sitot)

    lname = 'Vertical Integral of Conservative Subterms of Source Sink Term for Sitot, 0-100m'
    sname = 'Jint_100m_Sitot'
    units = 'mmol/m^3 cm/s'
    vgrid = 'none'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Jint_100m_Sitot)

    ! 3D Diags
    lname = 'Carbonate Ion Concentration'
    sname = 'CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%CO3)

    lname = 'Bicarbonate Ion Concentration'
    sname = 'HCO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%HCO3)

    lname = 'Carbonic Acid Concentration'
    sname = 'H2CO3'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%H2CO3)

    lname = 'pH'
    sname = 'pH_3D'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%ph_3D)

    lname = 'Carbonate Ion Concentration, Alternative CO2'
    sname = 'CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%CO3_ALT_CO2)

    lname = 'Bicarbonate Ion Concentration, Alternative CO2'
    sname = 'HCO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%HCO3_ALT_CO2)

    lname = 'Carbonic Acid Concentration, Alternative CO2'
    sname = 'H2CO3_ALT_CO2'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%H2CO3_ALT_CO2)

    lname = 'pH, Alternative CO2'
    sname = 'pH_3D_ALT_CO2'
    units = 'none'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%ph_3D_ALT_CO2)

    lname = 'CO3 concentration at calcite saturation'
    sname = 'co3_sat_calc'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%co3_sat_calc)

    lname = 'CO3 concentration at aragonite saturation'
    sname = 'co3_sat_arag'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%co3_sat_arag)

    lname = 'Nitrification'
    sname = 'NITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%NITRIF)

    lname = 'Denitrification'
    sname = 'DENITRIF'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DENITRIF)

    lname = 'O2 Production'
    sname = 'O2_PRODUCTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%O2_PRODUCTION)

    lname = 'O2 Consumption'
    sname = 'O2_CONSUMPTION'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%O2_CONSUMPTION)

    lname = 'Apparent O2 Utilization'
    sname = 'AOU'
    units = 'mmol/m^3'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%AOU)

    lname = 'PAR Average over Model Cell'
    sname = 'PAR_avg'
    !units = 'W/m^2'
    units = 'w/m^2' ! FIXME: watt / m^2 should be W/m^2?
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%PAR_avg)

    lname = 'Total Autotroph Grazing'
    sname = 'graze_auto_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%auto_graze_TOT)

    lname = 'Total C Fixation'
    sname = 'photoC_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%photoC_TOT)

    lname = 'Total C Fixation from NO3'
    sname = 'photoC_TOT'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%photoC_NO3_TOT)

    lname = 'DOC Production'
    sname = 'DOC_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOC_prod)

    lname = 'DOC Remineralization'
    sname = 'DOC_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOC_remin)

    lname = 'DON Production'
    sname = 'DON_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DON_prod)

    lname = 'DON Remineralization'
    sname = 'DON_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DON_remin)

    lname = 'DOP Production'
    sname = 'DOP_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOP_prod)

    lname = 'DOP Remineralization'
    sname = 'DOP_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOP_remin)

    lname = 'DOFe Production'
    sname = 'DOFe_prod'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOFe_prod)

    lname = 'DOFe Remineralization'
    sname = 'DOFe_remin'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .true.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%DOFe_remin)

    lname = 'Iron Scavenging'
    sname = 'Fe_scavenge'
    units = 'mmol/m^3/s'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Fe_scavenge)

    lname = 'Iron Scavenging Rate'
    sname = 'Fe_scavenge_rate'
    units = '1/y'
    vgrid = 'layer_avg'
    truncate = .false.
    call marbl_diags%add_diagnostic(lname, sname, units, vgrid, truncate,     &
                                    marbl_diag_ind%Fe_scavenge_rate)

    marbl_diag_ind%count = diag_cnt

    do n=1,part_diag_cnt_2d
      select case (n)
        case (calcToSed_diag_ind)
          lname = 'CaCO3 Flux to Sediments'
          sname = 'calcToSed'
          units = 'nmolC/cm^2/s'
        case (pocToSed_diag_ind)
          lname = 'POC Flux to Sediments'
          sname = 'pocToSed'
          units = 'nmolC/cm^2/s'
        case (ponToSed_diag_ind)
          lname = 'nitrogen burial Flux to Sediments'
          sname = 'ponToSed'
          units = 'nmolN/cm^2/s'
        case (SedDenitrif_diag_ind)
          lname = 'nitrogen loss in Sediments'
          sname = 'SedDenitrif'
          units = 'nmolN/cm^2/s'
        case (OtherRemin_diag_ind)
          lname = 'non-oxic,non-dentr remin in Sediments'
          sname = 'OtherRemin'
          units = 'nmolC/cm^2/s'
        case (popToSed_diag_ind)
          ! FIXME: should be phosphorus? (missing second h?)
!          lname = 'phosphorus Flux to Sediments'
          lname = 'phosporus Flux to Sediments'
          sname = 'popToSed'
          units = 'nmolP/cm^2/s'
        case (bsiToSed_diag_ind)
          lname = 'biogenic Si Flux to Sediments'
          sname = 'bsiToSed'
          units = 'nmolSi/cm^2/s'
        case (dustToSed_diag_ind)
          lname = 'dust Flux to Sediments'
          sname = 'dustToSed'
          units = 'g/cm^2/s'
        case (pfeToSed_diag_ind)
          lname = 'pFe Flux to Sediments'
          sname = 'pfeToSed'
          units = 'nmolFe/cm^2/s'
        case DEFAULT
          print*, "ERROR in ecosys_diagnostics_init():"
          print*, n, " is not a valid index for marbl_diags%part_diags_2d"
          sname = 'ERRORERRORERROR'
      end select
      part_diags_2d(n)%long_name = trim(lname)
      part_diags_2d(n)%short_name = trim(sname)
      part_diags_2d(n)%units = trim(units)
      part_diags_2d(n)%compute_now = .true.
    end do

    do auto_ind=1,autotroph_cnt
      do n=1,auto_diag_cnt_2d
        select case (n)
          case (photoC_zint_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) //                       &
                    ' C Fixation Vertical Integral'
            sname = 'photoC_' // trim(autotrophs(auto_ind)%sname) // '_zint'
            units = 'mmol/m^3 cm/s'
          case (photoC_NO3_zint_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) //                       &
                    ' C Fixation from NO3 Vertical Integral'
            sname = 'photoC_NO3_' // trim(autotrophs(auto_ind)%sname) // '_zint'
            units = 'mmol/m^3 cm/s'
          case (CACO3_form_zint_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) //                       &
                    ' CaCO3 Formation Vertical Integral'
            sname = trim(autotrophs(auto_ind)%sname) // '_CaCO3_form_zint'
            units = 'mmol/m^3 cm/s'
          case DEFAULT
            print*, "ERROR in ecosys_diagnostics_init():"
            print*, n, " is not a valid index for marbl_diags%diags_2d"
            sname = 'ERRORERRORERROR'
        end select
        auto_diags_2d(n, auto_ind)%long_name = trim(lname)
        auto_diags_2d(n, auto_ind)%short_name = trim(sname)
        auto_diags_2d(n, auto_ind)%units = trim(units)
        auto_diags_2d(n, auto_ind)%compute_now = .true.
      end do
    end do

    do n=1,part_diag_cnt_3d
      ! Default assumption: layer average vertical coordinates, not truncated
      vgrid = 'layer_avg'
      truncate = .false.
      select case (n)
        case (POC_FLUX_IN_diag_ind)
          lname = 'POC Flux into Cell'
          sname = 'POC_FLUX_IN'
          units = 'mmol/m^3 cm/s'
        case (POC_PROD_diag_ind)
          lname = 'POC Production'
          sname = 'POC_PROD'
          units = 'mmol/m^3/s'
        case (POC_REMIN_diag_ind)
          lname = 'POC Remineralization'
          sname = 'POC_REMIN'
          units = 'mmol/m^3/s'
        case (CaCO3_FLUX_IN_diag_ind)
!          lname = 'CaCO3 Flux into Cell'
          lname = 'CaCO3 flux into cell' ! FIXME lower case? consistent with rest?
          sname = 'CaCO3_FLUX_IN'
          units = 'mmol/m^3 cm/s'
        case (CaCO3_PROD_diag_ind)
          lname = 'CaCO3 Production'
          sname = 'CaCO3_PROD'
          units = 'mmol/m^3/s'
        case (CaCO3_REMIN_diag_ind)
          lname = 'CaCO3 Remineralization'
          sname = 'CaCO3_REMIN'
          units = 'mmol/m^3/s'
        case (SiO2_FLUX_IN_diag_ind)
          lname = 'SiO2 Flux into Cell'
          sname = 'SiO2_FLUX_IN'
          units = 'mmol/m^3 cm/s'
        case (SiO2_PROD_diag_ind)
          lname = 'SiO2 Production'
          sname = 'SiO2_PROD'
          units = 'mmol/m^3/s'
        case (SiO2_REMIN_diag_ind)
          lname = 'SiO2 Remineralization'
          sname = 'SiO2_REMIN'
          units = 'mmol/m^3/s'
        case (dust_FLUX_IN_diag_ind)
          lname = 'Dust Flux into Cell'
          sname = 'dust_FLUX_IN'
          units = 'ng/s/m^2'
        case (dust_REMIN_diag_ind)
          lname = 'Dust Remineralization'
          sname = 'dust_REMIN'
          units = 'mmol/m^3/s'
        case (P_iron_FLUX_IN_diag_ind)
          lname = 'P_iron Flux into Cell'
          sname = 'P_iron_FLUX_IN'
          units = 'mmol/m^3 cm/s'
        case (P_iron_PROD_diag_ind)
          lname = 'P_iron Production'
          sname = 'P_iron_PROD'
          units = 'mmol/m^3/s'
        case (P_iron_REMIN_diag_ind)
          lname = 'P_iron Remineralization'
          sname = 'P_iron_REMIN'
          units = 'mmol/m^3/s'
        case DEFAULT
          print*, "ERROR in ecosys_diagnostics_init():"
          print*, n, " is not a valid index for marbl_diags%diags_2d"
          sname = 'ERRORERRORERROR'
      end select
      part_diags_3d(n)%long_name = trim(lname)
      part_diags_3d(n)%short_name = trim(sname)
      part_diags_3d(n)%units = trim(units)
      part_diags_3d(n)%vertical_grid = trim(vgrid)
      part_diags_3d(n)%ltruncated_vertical_extent = truncate
      part_diags_3d(n)%compute_now = .true.
    end do

    do zoo_ind=1,zooplankton_cnt
      do n=1,zoo_diag_cnt_3d
        ! Default assumption: layer average vertical coordinates, not truncated
        vgrid = 'layer_avg'
        truncate = .false.
        select case (n)
          case (zoo_loss_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' Loss'
            sname = trim(zooplankton(zoo_ind)%sname) // '_loss'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_loss_poc_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' Loss to POC'
            sname = trim(zooplankton(zoo_ind)%sname) // '_loss_poc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_loss_doc_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' Loss to DOC'
            sname = trim(zooplankton(zoo_ind)%sname) // '_loss_doc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_graze_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' grazing loss'
            sname = 'graze_' // trim(zooplankton(zoo_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_graze_poc_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' grazing loss to POC'
            sname = 'graze_' // trim(zooplankton(zoo_ind)%sname) // '_poc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_graze_doc_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' grazing loss to DOC'
            sname = 'graze_' // trim(zooplankton(zoo_ind)%sname) // '_doc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (zoo_graze_zoo_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' grazing loss to ZOO'
            sname = 'graze_' // trim(zooplankton(zoo_ind)%sname) // '_zoo'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (x_graze_zoo_diag_ind)
            lname = trim(zooplankton(zoo_ind)%lname) // ' grazing gain'
            sname = 'x_graze_' // trim(zooplankton(zoo_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case DEFAULT
            print*, "ERROR in ecosys_diagnostics_init():"
            print*, n, " is not a valid index for marbl_diags%diags_2d"
            sname = 'ERRORERRORERROR'
        end select
        zoo_diags_3d(n,zoo_ind)%long_name = trim(lname)
        zoo_diags_3d(n,zoo_ind)%short_name = trim(sname)
        zoo_diags_3d(n,zoo_ind)%units = trim(units)
        zoo_diags_3d(n,zoo_ind)%vertical_grid = trim(vgrid)
        zoo_diags_3d(n,zoo_ind)%ltruncated_vertical_extent = truncate
        zoo_diags_3d(n,zoo_ind)%compute_now = .true.
      end do
    end do

    do auto_ind=1,autotroph_cnt
      do n=1,auto_diag_cnt_3d
        ! Default assumption: layer average vertical coordinates, not truncated
        vgrid = 'layer_avg'
        truncate = .false.
        select case (n)
          case (N_lim_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' N Limitation'
            sname = trim(autotrophs(auto_ind)%sname) // '_N_lim'
            units = 'none'
            truncate = .true.
          case (P_lim_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' P Limitation'
            sname = trim(autotrophs(auto_ind)%sname) // '_P_lim'
            units = 'none'
            truncate = .true.
          case (Fe_lim_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Fe Limitation'
            sname = trim(autotrophs(auto_ind)%sname) // '_Fe_lim'
            units = 'none'
            truncate = .true.
          case (SiO3_lim_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' SiO3 Limitation'
            sname = trim(autotrophs(auto_ind)%sname) // '_SiO3_lim'
            units = 'none'
            truncate = .true.
          case (light_lim_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Light Limitation'
            sname = trim(autotrophs(auto_ind)%sname) // '_light_lim'
            units = 'none'
            truncate = .true.
          case (photoC_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' C Fixation'
            sname = 'photoC_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (photoC_NO3_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' C Fixation from NO3'
            sname = 'photoC_NO3_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (photoFe_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Fe Uptake'
            sname = 'photoFe_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (photoNO3_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' NO3 Uptake'
            sname = 'photoNO3_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (photoNH4_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' NH4 Uptake'
            sname = 'photoNH4_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (DOP_uptake_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' DOP Uptake'
            sname = 'DOP_' // trim(autotrophs(auto_ind)%sname) // '_uptake'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (PO4_uptake_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' PO4 Uptake'
            sname = 'PO4_' // trim(autotrophs(auto_ind)%sname) // '_uptake'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_graze_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Grazing'
            sname = 'graze_' // trim(autotrophs(auto_ind)%sname)
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_graze_poc_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Grazing to POC'
            sname = 'graze_' // trim(autotrophs(auto_ind)%sname) // '_poc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_graze_doc_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Grazing to DOC'
            sname = 'graze_' // trim(autotrophs(auto_ind)%sname) // '_doc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_graze_zoo_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Grazing to ZOO'
            sname = 'graze_' // trim(autotrophs(auto_ind)%sname) // '_zoo'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_loss_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Loss'
            sname = trim(autotrophs(auto_ind)%sname) // '_loss'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_loss_poc_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Loss to POC'
            sname = trim(autotrophs(auto_ind)%sname) // '_loss_poc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_loss_doc_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Loss to DOC'
            sname = trim(autotrophs(auto_ind)%sname) // '_loss_doc'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (auto_agg_diag_ind)
            lname = trim(autotrophs(auto_ind)%lname) // ' Aggregate'
            sname = trim(autotrophs(auto_ind)%sname) // '_agg'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (bSi_form_diag_ind)
            lname=trim(autotrophs(auto_ind)%lname) // ' Si Uptake' ! FIXME: formation?
            sname = trim(autotrophs(auto_ind)%sname) // 'bSi_form' ! FIXME: _?
            units = 'mmol/m^3/s'
            truncate = .true.
          case (CaCO3_form_diag_ind)
            lname=trim(autotrophs(auto_ind)%lname) // ' CaCO3 Formation'
            sname = trim(autotrophs(auto_ind)%sname) // '_CaCO3_form'
            units = 'mmol/m^3/s'
            truncate = .true.
          case (Nfix_diag_ind)
            lname=trim(autotrophs(auto_ind)%lname) // ' N Fixation'
            sname = trim(autotrophs(auto_ind)%sname) // '_Nfix'
            units = 'mmol/m^3/s'
            truncate = .true.
          case DEFAULT
            print*, "ERROR in ecosys_diagnostics_init():"
            print*, n, " is not a valid index for marbl_diags%diags_2d"
            sname = 'ERRORERRORERROR'
        end select
        auto_diags_3d(n,auto_ind)%long_name = trim(lname)
        auto_diags_3d(n,auto_ind)%short_name = trim(sname)
        auto_diags_3d(n,auto_ind)%units = trim(units)
        auto_diags_3d(n,auto_ind)%vertical_grid = trim(vgrid)
        auto_diags_3d(n,auto_ind)%ltruncated_vertical_extent = truncate
        auto_diags_3d(n,auto_ind)%compute_now = .true.
      end do
    end do

    end associate

  end subroutine ecosys_diagnostics_init

  subroutine store_diagnostics_carbonate(marbl_domain, carbonate, zt, zw,     &
                                         marbl_diags)

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
      diags(marbl_diag_ind%zsatcalc)%field_2d =                               &
        compute_saturation_depth(marbl_domain, zt, zw, CO3, CO3_sat_calcite)
      diags(marbl_diag_ind%zsatarag)%field_2d =                               &
        compute_saturation_depth(marbl_domain, zt, zw, CO3, CO3_sat_aragonite)

      diags(marbl_diag_ind%CO3)%field_3d(:) = carbonate%CO3
      diags(marbl_diag_ind%HCO3)%field_3d(:) = carbonate%HCO3
      diags(marbl_diag_ind%H2CO3)%field_3d(:) = carbonate%H2CO3
      diags(marbl_diag_ind%pH_3D)%field_3d(:) = carbonate%pH
      diags(marbl_diag_ind%CO3_ALT_CO2)%field_3d(:) = carbonate%CO3_ALT_CO2
      diags(marbl_diag_ind%HCO3_ALT_CO2)%field_3d(:) = carbonate%HCO3_ALT_CO2
      diags(marbl_diag_ind%H2CO3_ALT_CO2)%field_3d(:) = carbonate%H2CO3_ALT_CO2
      diags(marbl_diag_ind%pH_3D_ALT_CO2)%field_3d(:) = carbonate%pH_ALT_CO2
      diags(marbl_diag_ind%co3_sat_calc)%field_3d(:) = carbonate%CO3_sat_calcite
      diags(marbl_diag_ind%co3_sat_arag)%field_3d(:) = carbonate%CO3_sat_aragonite
    end associate

  end subroutine store_diagnostics_carbonate

  function compute_saturation_depth(marbl_domain, zt, zw, CO3, sat_val)

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zt, zw, CO3, sat_val
    real(r8) :: compute_saturation_depth

    real(r8) :: anomaly(km) ! CO3 concentration above saturation at level
    integer :: k

    associate(&
              column_kmt => marbl_domain%kmt                                  &
             )
      anomaly   = CO3 - sat_val

      if (all(anomaly(1:column_kmt).gt.c0)) then
        compute_saturation_depth = zw(column_kmt)
      elseif (anomaly(1).le.c0) then
        compute_saturation_depth = c0
      else
        do k=2,column_kmt
          if (anomaly(k).le.c0) exit
        end do
        ! saturation depth is location of root of anomaly
        compute_saturation_depth = linear_root(zt(k-1:k), anomaly(k-1:k))
      end if
    end associate

  end function compute_saturation_depth

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

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_nitrification(nitrif, denitrif, marbl_diags)

    real(r8), dimension(:), intent(in) :: nitrif
    real(r8), dimension(:), intent(in) :: denitrif
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_diag_ind%NITRIF)%field_3d(:) = nitrif
      diags(marbl_diag_ind%DENITRIF)%field_3d(:) = denitrif
    end associate

  end subroutine store_diagnostics_nitrification

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_autotrophs(marbl_domain,                       &
                                          autotroph_secondary_species,        &
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

    associate(&
              auto_diags_2d => marbl_diags%auto_diags_2d(:,:)%field,          &
              auto_diags_3d => marbl_diags%auto_diags_3d(:,:)                 &
             )
      do n = 1, autotroph_cnt
         auto_diags_3d(N_lim_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%VNtot
         auto_diags_3d(Fe_lim_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%VFe
         auto_diags_3d(P_lim_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%VPtot

         if (autotrophs(n)%kSiO3 > c0) then
            auto_diags_3d(SiO3_lim_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%VSiO3
         end if

         auto_diags_3d(light_lim_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%light_lim
         auto_diags_3d(photoNO3_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%NO3_V
         auto_diags_3d(photoNH4_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%NH4_V
         auto_diags_3d(PO4_uptake_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%PO4_V
         auto_diags_3d(DOP_uptake_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%DOP_V
         auto_diags_3d(photoFE_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%photoFe

         if (autotrophs(n)%Si_ind > 0) then
           auto_diags_3d(bSi_form_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%photoSi
         endif

         auto_diags_3d(CaCO3_form_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%CaCO3_PROD
         auto_diags_3d(Nfix_diag_ind, n)%field(:) = autotroph_secondary_species(n,:)%Nfix
         auto_diags_3d(auto_graze_diag_ind, n)%field(:)      = autotroph_secondary_species(n,:)%auto_graze
         auto_diags_3d(auto_graze_poc_diag_ind, n)%field(:)  = autotroph_secondary_species(n,:)%auto_graze_poc
         auto_diags_3d(auto_graze_doc_diag_ind, n)%field(:)  = autotroph_secondary_species(n,:)%auto_graze_doc
         auto_diags_3d(auto_graze_zoo_diag_ind, n)%field(:)  = autotroph_secondary_species(n,:)%auto_graze_zoo
         auto_diags_3d(auto_loss_diag_ind, n)%field(:)       = autotroph_secondary_species(n,:)%auto_loss
         auto_diags_3d(auto_loss_poc_diag_ind, n)%field(:)   = autotroph_secondary_species(n,:)%auto_loss_poc
         auto_diags_3d(auto_loss_doc_diag_ind, n)%field(:)   = autotroph_secondary_species(n,:)%auto_loss_doc
         auto_diags_3d(auto_agg_diag_ind, n)%field(:)        = autotroph_secondary_species(n,:)%auto_agg
         auto_diags_3d(photoC_diag_ind, n)%field(:)          = autotroph_secondary_species(n,:)%photoC

         auto_diags_3d(photoC_NO3_diag_ind, n)%field(:) = c0
         where (autotroph_secondary_species(n,:)%VNtot > c0)
           auto_diags_3d(photoC_NO3_diag_ind, n)%field(:) =                   &
 autotroph_secondary_species(n,:)%photoC *                                    &
 (autotroph_secondary_species(n,:)%VNO3 / autotroph_secondary_species(n,:)%VNtot)
         end where

         ! vertical integrals
         auto_diags_2d(CaCO3_form_zint_diag_ind, n) = c0
         auto_diags_2d(photoC_zint_diag_ind, n) = c0
         auto_diags_2d(photoC_NO3_zint_diag_ind, n) = c0
         do k = 1,marbl_domain%kmt
           auto_diags_2d(CaCO3_form_zint_diag_ind, n) =                       &
                     auto_diags_2d(CaCO3_form_zint_diag_ind, n) +             &
                     delta_z(k) * autotroph_secondary_species(n,k)%CaCO3_PROD
           auto_diags_2d(photoC_zint_diag_ind, n) =                           &
                         auto_diags_2d(photoC_zint_diag_ind, n) +             &
                         delta_z(k) * autotroph_secondary_species(n,k)%photoC
           auto_diags_2d(photoC_NO3_zint_diag_ind, n) =                       &
                   auto_diags_2d(photoC_NO3_zint_diag_ind, n) +               &
                   delta_z(k) * auto_diags_3d(photoC_NO3_diag_ind,n)%field(k)
         end do ! do k
      end do ! do n
    end associate

  end subroutine store_diagnostics_autotrophs

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_autotroph_sums(marbl_domain,                   &
                                              autotroph_secondary_species,    &
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
      diags(marbl_diag_ind%auto_graze_TOT)%field_3d(:) =                      &
                           sum(autotroph_secondary_species%auto_graze, dim=1)
      diags(marbl_diag_ind%photoC_TOT)%field_3d(:) =                          &
                           sum(autotroph_secondary_species%photoC, dim=1)

      diags(marbl_diag_ind%photoC_NO3_TOT)%field_3d(:) = c0
      do n = 1, autotroph_cnt
        where (autotroph_secondary_species(n,:)%VNtot > c0)
          diags(marbl_diag_ind%photoC_NO3_TOT)%field_3d(:) =                  &
                           diags(marbl_diag_ind%photoC_NO3_TOT)%field_3d(:) + &
                           (autotroph_secondary_species(n,:)%VNO3 /           &
                           autotroph_secondary_species(n,:)%VNtot) *          &
                           autotroph_secondary_species(n,:)%photoC
        end where
      end do

      diags(marbl_diag_ind%photoC_TOT_zint)%field_2d = sum(                   &
                    delta_z * sum(autotroph_secondary_species%photoC, dim=1))
      diags(marbl_diag_ind%photoC_NO3_TOT_zint)%field_2d = sum(               &
                    delta_z * diags(marbl_diag_ind%photoC_NO3_TOT)%field_3d)
    end associate

  end subroutine store_diagnostics_autotroph_sums

  !-----------------------------------------------------------------------

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
    associate(&
              part_diags_2d => marbl_diags%part_diags_2d(:)%field,            &
              part_diags_3d => marbl_diags%part_diags_3d(:)                   &
              )
      part_diags_3d(POC_FLUX_IN_diag_ind)%field(:) = POC%sflux_in + POC%hflux_in
      part_diags_3d(POC_PROD_diag_ind)%field(:) = POC%prod
      part_diags_3d(POC_REMIN_diag_ind)%field(:) = POC%remin

      part_diags_3d(CaCO3_FLUX_IN_diag_ind)%field(:) = P_CaCO3%sflux_in + P_CaCO3%hflux_in
      part_diags_3d(CaCO3_PROD_diag_ind)%field(:) = P_CaCO3%prod
      part_diags_3d(CaCO3_REMIN_diag_ind)%field(:) = P_CaCO3%remin

      part_diags_3d(SiO2_FLUX_IN_diag_ind)%field(:) = P_SiO2%sflux_in + P_SiO2%hflux_in
      part_diags_3d(SiO2_PROD_diag_ind)%field(:) = P_SiO2%prod
      part_diags_3d(SiO2_REMIN_diag_ind)%field(:) = P_SiO2%remin

      part_diags_3d(dust_FLUX_IN_diag_ind)%field(:) = dust%sflux_in + dust%hflux_in
      part_diags_3d(dust_REMIN_diag_ind)%field(:) = P_SiO2%remin

      part_diags_3d(P_iron_FLUX_IN_diag_ind)%field(:) = P_iron%sflux_in + P_iron%hflux_in
      part_diags_3d(P_iron_PROD_diag_ind)%field(:) = P_iron%prod
      part_diags_3d(P_iron_REMIN_diag_ind)%field(:) = P_iron%remin

      part_diags_2d(calcToSed_diag_ind) = sum(P_CaCO3%sed_loss)
      part_diags_2d(bsiToSed_diag_ind) = sum(P_SiO2%sed_loss)
      part_diags_2d(pocToSed_diag_ind) = sum(POC%sed_loss)
      part_diags_2d(SedDenitrif_diag_ind) = sum(sed_denitrif * delta_z)
      part_diags_2d(OtherRemin_diag_ind) = sum(other_remin * delta_z)
      part_diags_2d(ponToSed_diag_ind) = sum(POC%sed_loss * Q)
      part_diags_2d(popToSed_diag_ind) = sum(POC%sed_loss * Qp_zoo_pom)
      part_diags_2d(dustToSed_diag_ind) = sum(dust%sed_loss)
      part_diags_2d(pfeToSed_diag_ind) = sum(P_iron%sed_loss)
    end associate

  end subroutine store_diagnostics_particulates

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_oxygen(marbl_domain, column_zt,                &
                                      column_o2, o2_production, o2_consumption,&
                                      marbl_diags)

    use marbl_oxygen, only : o2sat_scalar

    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: column_zt
    real(r8), dimension(:), intent(in) :: column_o2
    real(r8), dimension(:), intent(in) :: o2_production
    real(r8), dimension(:), intent(in) :: o2_consumption
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: k, min_ind

    ! Find min_o2 and depth_min_o2
    associate(diags => marbl_diags%diags(:))
      min_ind = minloc(column_o2(1:marbl_domain%kmt), dim=1)
      diags(marbl_diag_ind%O2_ZMIN)%field_2d = column_o2(min_ind)
      diags(marbl_diag_ind%O2_ZMIN_DEPTH)%field_2d = column_zt(min_ind)

      diags(marbl_diag_ind%O2_PRODUCTION)%field_3d(:) = o2_production
      diags(marbl_diag_ind%O2_CONSUMPTION)%field_3d(:) = o2_consumption

      diags(marbl_diag_ind%AOU)%field_3d(:) = -column_o2
      do k=1,marbl_domain%kmt
        if (marbl_domain%land_mask) then
          diags(marbl_diag_ind%AOU)%field_3d(k) =                             &
        O2SAT_scalar(marbl_domain%temperature(k), marbl_domain%salinity(k)) - &
        column_o2(k)
        end if
      end do
    end associate

  end subroutine store_diagnostics_oxygen

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_photosynthetically_available_radiation( &
       PAR, marbl_diags)

    type(photosynthetically_available_radiation_type), dimension(:), intent(in) :: PAR
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_diag_ind%PAR_avg)%field_3d(:) = PAR%avg
    end associate

  end subroutine store_diagnostics_photosynthetically_available_radiation

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_misc(marbl_diags)
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags
  end subroutine store_diagnostics_misc

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_zooplankton(zooplankton_secondary_species,     &
                                           marbl_diags)

    type(zooplankton_secondary_species_type), dimension(:,:), intent(in) ::   &
                                              zooplankton_secondary_species
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    integer(int_kind) :: n

    associate(zoo_diags_3d => marbl_diags%zoo_diags_3d(:,:))
      do n = 1, zooplankton_cnt
        zoo_diags_3d(zoo_loss_diag_ind, n)%field(:)       = zooplankton_secondary_species(n,:)%zoo_loss
        zoo_diags_3d(zoo_loss_poc_diag_ind, n)%field(:)   = zooplankton_secondary_species(n,:)%zoo_loss_poc
        zoo_diags_3d(zoo_loss_doc_diag_ind, n)%field(:)   = zooplankton_secondary_species(n,:)%zoo_loss_doc
        zoo_diags_3d(zoo_graze_diag_ind, n)%field(:)      = zooplankton_secondary_species(n,:)%zoo_graze
        zoo_diags_3d(zoo_graze_poc_diag_ind, n)%field(:)  = zooplankton_secondary_species(n,:)%zoo_graze_poc
        zoo_diags_3d(zoo_graze_doc_diag_ind, n)%field(:)  = zooplankton_secondary_species(n,:)%zoo_graze_doc
        zoo_diags_3d(zoo_graze_zoo_diag_ind, n)%field(:)  = zooplankton_secondary_species(n,:)%zoo_graze_zoo
        zoo_diags_3d(x_graze_zoo_diag_ind, n)%field(:)    = zooplankton_secondary_species(n,:)%x_graze_zoo
      end do
    end associate

  end subroutine store_diagnostics_zooplankton

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_dissolved_organic_matter(&
       dissolved_organic_matter, fe_scavenge, fe_scavenge_rate, marbl_diags)

    type(dissolved_organic_matter_type), dimension(:), intent(in) ::          &
                                            dissolved_organic_matter
    real(r8), dimension(:), intent(in) :: fe_scavenge, fe_scavenge_rate
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags => marbl_diags%diags(:))
      diags(marbl_diag_ind%DOC_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DOC_prod
      diags(marbl_diag_ind%DOC_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOC_remin
      diags(marbl_diag_ind%DON_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DON_prod
      diags(marbl_diag_ind%DON_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DON_remin
      diags(marbl_diag_ind%DOP_prod)%field_3d(:) =                            &
                                          dissolved_organic_matter%DOP_prod
      diags(marbl_diag_ind%DOP_remin)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOP_remin
      diags(marbl_diag_ind%DOFe_prod)%field_3d(:) =                           &
                                          dissolved_organic_matter%DOFe_prod
      diags(marbl_diag_ind%DOFe_remin)%field_3d(:) =                          &
                                          dissolved_organic_matter%DOFe_remin

      diags(marbl_diag_ind%Fe_scavenge)%field_3d(:) = Fe_scavenge
      diags(marbl_diag_ind%Fe_scavenge_rate)%field_3d(:) = Fe_scavenge_rate
    end associate

  end subroutine store_diagnostics_dissolved_organic_matter

  !-----------------------------------------------------------------------

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
      diags(marbl_diag_ind%Jint_Ctot)%field_2d = c0
      diags(marbl_diag_ind%Jint_100m_Ctot)%field_2d = c0

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
        diags(marbl_diag_ind%Jint_Ctot)%field_2d =                            &
                 diags(marbl_diag_ind%Jint_Ctot)%field_2d +                   &
                 delta_z(k) * work1 + (POC%sed_loss(k) + P_CaCO3%sed_loss(k))

        if (ztop < 100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Ctot)%field_2d =                     &
                              diags(marbl_diag_ind%Jint_100m_Ctot)%field_2d + &
                              min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Ctot)%field_2d =                     &
                              diags(marbl_diag_ind%Jint_100m_Ctot)%field_2d + &
                              (POC%sed_loss(k) + P_CaCO3%sed_loss(k))
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_carbon_fluxes

  !-----------------------------------------------------------------------

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
      diags(marbl_diag_ind%Jint_Ntot)%field_2d = c0
      diags(marbl_diag_ind%Jint_100m_Ntot)%field_2d = c0

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
        diags(marbl_diag_ind%Jint_Ntot)%field_2d =                            &
                                   diags(marbl_diag_ind%Jint_Ntot)%field_2d + &
                                   delta_z(k) * work1 + POC%sed_loss(k) * Q

        if (ztop < 100.0e2_r8) then
            diags(marbl_diag_ind%Jint_100m_Ntot)%field_2d =                   &
                              diags(marbl_diag_ind%Jint_100m_Ntot)%field_2d + &
                           min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
            diags(marbl_diag_ind%Jint_100m_Ntot)%field_2d =                   &
                              diags(marbl_diag_ind%Jint_100m_Ntot)%field_2d + &
                              POC%sed_loss(k) * Q
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_nitrogen_fluxes

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_phosphorus_fluxes(marbl_domain, zw, POC,       &
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
      diags(marbl_diag_ind%Jint_Ptot)%field_2d = c0
      diags(marbl_diag_ind%Jint_100m_Ptot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(po4_ind,k) + dtracer(dop_ind,k) + dtracer(dopr_ind,k)
        do n = 1, zooplankton_cnt
          work1 = work1 + Qp_zoo_pom * dtracer(zooplankton(n)%C_ind,k)
        end do
        do n = 1, autotroph_cnt
          work1 = work1 + autotrophs(n)%Qp * dtracer(autotrophs(n)%C_ind,k)
        end do
        diags(marbl_diag_ind%Jint_Ptot)%field_2d =                            &
                            diags(marbl_diag_ind%Jint_Ptot)%field_2d +        &
                            delta_z(k) * work1 + POC%sed_loss(k) * Qp_zoo_pom

        if (ztop < 100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Ptot)%field_2d =                     &
                              diags(marbl_diag_ind%Jint_100m_Ptot)%field_2d + &
                              min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Ptot)%field_2d =                     &
                              diags(marbl_diag_ind%Jint_100m_Ptot)%field_2d + &
                              POC%sed_loss(k) * Qp_zoo_pom
        end if
        ztop = zw(k)
    end do
    end associate

  end subroutine store_diagnostics_phosphorus_fluxes

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_silicon_fluxes(marbl_domain, zw, P_SiO2,       &
                                              dtracer, marbl_diags)

    use marbl_share_mod, only : column_sinking_particle_type

    use marbl_parms     , only : Q
    use marbl_parms     , only : Qp_zoo_pom

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
      diags(marbl_diag_ind%Jint_Sitot)%field_2d = c0
      diags(marbl_diag_ind%Jint_100m_Sitot)%field_2d = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(sio3_ind,k)
        do n = 1, autotroph_cnt
          if (autotrophs(n)%Si_ind > 0) then
            work1 = work1 + dtracer(autotrophs(n)%Si_ind,k)
          end if
        end do
        diags(marbl_diag_ind%Jint_Sitot)%field_2d =                           &
                                  diags(marbl_diag_ind%Jint_Sitot)%field_2d + &
                                  delta_z(k) * work1 + P_SiO2%sed_loss(k)

        if (ztop < 100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Sitot)%field_2d =                    &
                             diags(marbl_diag_ind%Jint_100m_Sitot)%field_2d + &
                             min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags(marbl_diag_ind%Jint_100m_Sitot)%field_2d =                    &
                             diags(marbl_diag_ind%Jint_100m_Sitot)%field_2d + &
                             P_SiO2%sed_loss(k)
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_silicon_fluxes

  !-----------------------------------------------------------------------

end module ecosys_diagnostics_mod
