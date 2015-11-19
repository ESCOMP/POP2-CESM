module ecosys_constants

  use kinds_mod, only : int_kind

  implicit none

  public

  integer(int_kind), parameter :: ecosys_tracer_cnt = ECOSYS_NT

!-----------------------------------------------------------------------
!  indices for diagnostic values written to tavg files
!-----------------------------------------------------------------------

  integer(int_kind), parameter :: ecosys_diag_cnt = 45
  integer (int_kind), parameter ::  &
      CO3_diag_ind                 =  1, &
      HCO3_diag_ind                =  2, &
      H2CO3_diag_ind               =  3, &
      pH_3D_diag_ind               =  4, &
      CO3_ALT_CO2_diag_ind         =  5, &
      HCO3_ALT_CO2_diag_ind        =  6, &
      H2CO3_ALT_CO2_diag_ind       =  7, &
      pH_3D_ALT_CO2_diag_ind       =  8, &
      co3_sat_calc_diag_ind        =  9, &
      zsatcalc_diag_ind            = 10, &
      co3_sat_arag_diag_ind        = 11, &
      zsatarag_diag_ind            = 12, &
      NITRIF_diag_ind              = 13, &
      DENITRIF_diag_ind            = 14, &
      O2_ZMIN_diag_ind             = 15, &
      O2_ZMIN_DEPTH_diag_ind       = 16, &
      O2_PRODUCTION_diag_ind       = 17, &
      O2_CONSUMPTION_diag_ind      = 18, &
      AOU_diag_ind                 = 19, &
      PAR_avg_diag_ind             = 20, &
      auto_graze_TOT_diag_ind      = 21, &
      photoC_TOT_diag_ind          = 22, &
      photoC_TOT_zint_diag_ind     = 23, &
      photoC_NO3_TOT_diag_ind      = 24, &
      photoC_NO3_TOT_zint_diag_ind = 25, &
      DOC_prod_diag_ind            = 26, &
      DOC_remin_diag_ind           = 27, &
      DON_prod_diag_ind            = 28, &
      DON_remin_diag_ind           = 29, &
      DOP_prod_diag_ind            = 30, &
      DOP_remin_diag_ind           = 31, &
      DOFe_prod_diag_ind           = 32, &
      DOFe_remin_diag_ind          = 33, &
      Fe_scavenge_diag_ind         = 34, &
      Fe_scavenge_rate_diag_ind    = 35, &
      Jint_Ctot_diag_ind           = 36, &
      Jint_100m_Ctot_diag_ind      = 37, &
      Jint_Ntot_diag_ind           = 38, &
      Jint_100m_Ntot_diag_ind      = 39, &
      Jint_Ptot_diag_ind           = 40, &
      Jint_100m_Ptot_diag_ind      = 41, &
      Jint_Sitot_diag_ind          = 42, &
      Jint_100m_Sitot_diag_ind     = 43, &
      Jint_Fetot_diag_ind          = 44, &
      Jint_100m_Fetot_diag_ind     = 45

  integer(int_kind), parameter ::   auto_diag_cnt = 26
  integer (int_kind), parameter ::  &
      N_lim_diag_ind           =  1,  &
      P_lim_diag_ind           =  2,  &
      Fe_lim_diag_ind          =  3,  &
      SiO3_lim_diag_ind        =  4,  &
      light_lim_diag_ind       =  5,  &
      photoC_diag_ind          =  6,  &
      photoC_zint_diag_ind     =  7,  &
      photoC_NO3_diag_ind      =  8,  &
      photoC_NO3_zint_diag_ind =  9,  &
      photoFe_diag_ind         = 10,  &
      photoNO3_diag_ind        = 11,  &
      photoNH4_diag_ind        = 12,  &
      DOP_uptake_diag_ind      = 13,  &
      PO4_uptake_diag_ind      = 14,  &
      auto_graze_diag_ind      = 15,  &
      auto_graze_poc_diag_ind  = 16,  &
      auto_graze_doc_diag_ind  = 17,  &
      auto_graze_zoo_diag_ind  = 18,  &
      auto_loss_diag_ind       = 19,  &
      auto_loss_poc_diag_ind   = 20,  &
      auto_loss_doc_diag_ind   = 21,  &
      auto_agg_diag_ind        = 22,  &
      bSi_form_diag_ind        = 23,  &
      CaCO3_form_diag_ind      = 24,  &
      CaCO3_form_zint_diag_ind = 25,  &
      Nfix_diag_ind            = 26

  integer(int_kind), parameter ::    zoo_diag_cnt =  8
  integer (int_kind), parameter ::   &
      zoo_loss_diag_ind        =  1, &
      zoo_loss_poc_diag_ind    =  2, &
      zoo_loss_doc_diag_ind    =  3, &
      zoo_graze_diag_ind       =  4, &
      zoo_graze_poc_diag_ind   =  5, &
      zoo_graze_doc_diag_ind   =  6, &
      zoo_graze_zoo_diag_ind   =  7, &
      x_graze_zoo_diag_ind     =  8

  integer(int_kind), parameter ::   part_diag_cnt =  23
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
      P_iron_REMIN_diag_ind    = 14, &
      calcToSed_diag_ind       = 15, &
      bsiToSed_diag_ind        = 16, &
      pocToSed_diag_ind        = 17, &
      SedDenitrif_diag_ind     = 18, &
      OtherRemin_diag_ind      = 19, &
      ponToSed_diag_ind        = 20, &
      popToSed_diag_ind        = 21, &
      dustToSed_diag_ind       = 22, &
      pfeToSed_diag_ind        = 23

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

end module ecosys_constants
