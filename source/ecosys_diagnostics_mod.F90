module ecosys_diagnostics_mod

  use ecosys_constants, only : ecosys_tracer_cnt

  use marbl_kinds_mod, only : c0
  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : log_kind

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
  use marbl_interface_types, only : marbl_gcm_state_type

  use grid, only : partial_bottom_cells
  use domain_size, only : km

  use marbl_parms, only : po4_ind 
  use marbl_parms, only : no3_ind         
  use marbl_parms, only : sio3_ind        
  use marbl_parms, only : nh4_ind         
  use marbl_parms, only : fe_ind          
  use marbl_parms, only : o2_ind          
  use marbl_parms, only : dic_ind         
  use marbl_parms, only : dic_alt_co2_ind 
  use marbl_parms, only : alk_ind         
  use marbl_parms, only : doc_ind         
  use marbl_parms, only : don_ind         
  use marbl_parms, only : dofe_ind        
  use marbl_parms, only : dop_ind         
  use marbl_parms, only : dopr_ind        
  use marbl_parms, only : donr_ind        

  Implicit None
  Public

!-----------------------------------------------------------------------
!  indices for diagnostic values written to tavg files
!-----------------------------------------------------------------------

  integer(int_kind), parameter :: ecosys_diag_cnt_2d =  6
  integer (int_kind), parameter ::  &
      zsatcalc_diag_ind            = 1, &
      zsatarag_diag_ind            = 2, &
      O2_ZMIN_diag_ind             = 3, &
      O2_ZMIN_DEPTH_diag_ind       = 4, &
      photoC_TOT_zint_diag_ind     = 5, &
      photoC_NO3_TOT_zint_diag_ind = 6

  integer(int_kind), parameter :: ecosys_diag_cnt_3d = 37
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
      co3_sat_arag_diag_ind        = 10, &
      NITRIF_diag_ind              = 11, &
      DENITRIF_diag_ind            = 12, &
      O2_PRODUCTION_diag_ind       = 13, &
      O2_CONSUMPTION_diag_ind      = 14, &
      AOU_diag_ind                 = 15, &
      PAR_avg_diag_ind             = 16, &
      auto_graze_TOT_diag_ind      = 17, &
      photoC_TOT_diag_ind          = 18, &
      photoC_NO3_TOT_diag_ind      = 19, &
      DOC_prod_diag_ind            = 20, &
      DOC_remin_diag_ind           = 21, &
      DON_prod_diag_ind            = 22, &
      DON_remin_diag_ind           = 23, &
      DOP_prod_diag_ind            = 24, &
      DOP_remin_diag_ind           = 25, &
      DOFe_prod_diag_ind           = 26, &
      DOFe_remin_diag_ind          = 27, &
      Fe_scavenge_diag_ind         = 28, &
      Fe_scavenge_rate_diag_ind    = 29, &
      Jint_Ctot_diag_ind           = 30, &
      Jint_100m_Ctot_diag_ind      = 31, &
      Jint_Ntot_diag_ind           = 32, &
      Jint_100m_Ntot_diag_ind      = 33, &
      Jint_Ptot_diag_ind           = 34, &
      Jint_100m_Ptot_diag_ind      = 35, &
      Jint_Sitot_diag_ind          = 36, &
      Jint_100m_Sitot_diag_ind     = 37

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

contains

  subroutine store_diagnostics_carbonate(carbonate, zsat_calcite,             &
                                         zsat_aragonite, marbl_diags)

    type(carbonate_type), dimension(:), intent(in) :: carbonate
    real(r8), dimension(:), intent(in) :: zsat_calcite
    real(r8), dimension(:), intent(in) :: zsat_aragonite
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d                                &
              )
      diags_3d(:, CO3_diag_ind) = carbonate%CO3
      diags_3d(:, HCO3_diag_ind) = carbonate%HCO3
      diags_3d(:, H2CO3_diag_ind) = carbonate%H2CO3
      diags_3d(:, pH_3D_diag_ind) = carbonate%pH
      diags_3d(:, CO3_ALT_CO2_diag_ind) = carbonate%CO3_ALT_CO2
      diags_3d(:, HCO3_ALT_CO2_diag_ind) = carbonate%HCO3_ALT_CO2
      diags_3d(:, H2CO3_ALT_CO2_diag_ind) = carbonate%H2CO3_ALT_CO2
      diags_3d(:, pH_3D_ALT_CO2_diag_ind) = carbonate%pH_ALT_CO2
      diags_3d(:, co3_sat_calc_diag_ind) = carbonate%CO3_sat_calcite
      diags_3d(:, co3_sat_arag_diag_ind) = carbonate%CO3_sat_aragonite
      diags_2d(:, zsatcalc_diag_ind) = zsat_calcite
      diags_2d(:, zsatarag_diag_ind) = zsat_aragonite
    end associate

  end subroutine store_diagnostics_carbonate

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_nitrification(nitrif, denitrif, marbl_diags)

    real(r8), dimension(:), intent(in) :: nitrif
    real(r8), dimension(:), intent(in) :: denitrif
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    associate(diags_3d => marbl_diags%diags_3d)
      diags_3d(:, NITRIF_diag_ind) = nitrif
      diags_3d(:, DENITRIF_diag_ind) = denitrif
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

    associate(auto_diags => marbl_diags%auto_diags)
      do n = 1, autotroph_cnt
         auto_diags(:, N_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VNtot
         auto_diags(:, Fe_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VFe
         auto_diags(:, P_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VPtot

         if (autotrophs(n)%kSiO3 > c0) then
            auto_diags(:, SiO3_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VSiO3
         end if

         auto_diags(:, light_lim_diag_ind, n) = autotroph_secondary_species(n,:)%light_lim
         auto_diags(:, photoNO3_diag_ind, n) = autotroph_secondary_species(n,:)%NO3_V
         auto_diags(:, photoNH4_diag_ind, n) = autotroph_secondary_species(n,:)%NH4_V
         auto_diags(:, PO4_uptake_diag_ind, n) = autotroph_secondary_species(n,:)%PO4_V
         auto_diags(:, DOP_uptake_diag_ind, n) = autotroph_secondary_species(n,:)%DOP_V
         auto_diags(:, photoFE_diag_ind, n) = autotroph_secondary_species(n,:)%photoFe

         if (autotrophs(n)%Si_ind > 0) then
           auto_diags(:, bSi_form_diag_ind, n) = autotroph_secondary_species(n,:)%photoSi
         endif

         auto_diags(:, CaCO3_form_diag_ind, n) = autotroph_secondary_species(n,:)%CaCO3_PROD
         auto_diags(:, Nfix_diag_ind, n) = autotroph_secondary_species(n,:)%Nfix
         auto_diags(:, auto_graze_diag_ind, n)      = autotroph_secondary_species(n,:)%auto_graze
         auto_diags(:, auto_graze_poc_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_poc
         auto_diags(:, auto_graze_doc_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_doc
         auto_diags(:, auto_graze_zoo_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_zoo
         auto_diags(:, auto_loss_diag_ind, n)       = autotroph_secondary_species(n,:)%auto_loss
         auto_diags(:, auto_loss_poc_diag_ind, n)   = autotroph_secondary_species(n,:)%auto_loss_poc
         auto_diags(:, auto_loss_doc_diag_ind, n)   = autotroph_secondary_species(n,:)%auto_loss_doc
         auto_diags(:, auto_agg_diag_ind, n)        = autotroph_secondary_species(n,:)%auto_agg
         auto_diags(:, photoC_diag_ind, n)          = autotroph_secondary_species(n,:)%photoC

         auto_diags(:, CaCO3_form_zint_diag_ind, n) = c0
         auto_diags(:, photoC_zint_diag_ind, n) = c0
         auto_diags(:, photoC_NO3_zint_diag_ind, n) = c0
         auto_diags(:, photoC_NO3_diag_ind, n) = c0
         where (autotroph_secondary_species(n,:)%VNtot > c0)
           auto_diags(:, photoC_NO3_diag_ind, n) =                            &
 autotroph_secondary_species(n,:)%photoC *                                    &
 (autotroph_secondary_species(n,:)%VNO3 / autotroph_secondary_species(n,:)%VNtot)
         end where

         ! vertical integrals
         do k = 1,marbl_domain%kmt
           auto_diags(k, CaCO3_form_zint_diag_ind, n) = delta_z(k) *          &
                                  autotroph_secondary_species(n,k)%CaCO3_PROD
           auto_diags(k, photoC_zint_diag_ind, n) = delta_z(k) *              &
                                      autotroph_secondary_species(n,k)%photoC
           auto_diags(k, photoC_NO3_zint_diag_ind, n) = delta_z(k)*           &
                                        auto_diags(k, photoC_NO3_diag_ind, n)
         end do
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

    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d                                &
              )
      diags_3d(:, auto_graze_TOT_diag_ind) =                                  &
                           sum(autotroph_secondary_species%auto_graze, dim=1)
      diags_3d(:, photoC_TOT_diag_ind) =                                      &
                           sum(autotroph_secondary_species%photoC, dim=1)
      diags_2d(:, photoC_TOT_zint_diag_ind) =                                 &
                 delta_z * sum(autotroph_secondary_species%photoC, dim=1)

      diags_3d(:, photoC_NO3_TOT_diag_ind) = c0
      do n = 1, autotroph_cnt
        where (autotroph_secondary_species(n,:)%VNtot > c0)
          diags_3d(:, photoC_NO3_TOT_diag_ind) =                          &
                                   diags_3d(:, photoC_NO3_TOT_diag_ind) + &
                                 (autotroph_secondary_species(n,:)%VNO3 / &
                                autotroph_secondary_species(n,:)%VNtot) * &
                                autotroph_secondary_species(n,:)%photoC
        end where
      end do
      diags_2d(:, photoC_NO3_TOT_zint_diag_ind) =                             &
                              delta_z * diags_3d(:, photoC_NO3_TOT_diag_ind)
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
    associate(part_diags => marbl_diags%part_diags)
      part_diags(:, POC_FLUX_IN_diag_ind) = POC%sflux_in + POC%hflux_in
      part_diags(:, POC_PROD_diag_ind) = POC%prod
      part_diags(:, POC_REMIN_diag_ind) = POC%remin

      part_diags(:, CaCO3_FLUX_IN_diag_ind) = P_CaCO3%sflux_in + P_CaCO3%hflux_in
      part_diags(:, CaCO3_PROD_diag_ind) = P_CaCO3%prod
      part_diags(:, CaCO3_REMIN_diag_ind) = P_CaCO3%remin

      part_diags(:, SiO2_FLUX_IN_diag_ind) = P_SiO2%sflux_in + P_SiO2%hflux_in
      part_diags(:, SiO2_PROD_diag_ind) = P_SiO2%prod
      part_diags(:, SiO2_REMIN_diag_ind) = P_SiO2%remin

      part_diags(:, dust_FLUX_IN_diag_ind) = dust%sflux_in + dust%hflux_in
      part_diags(:, dust_REMIN_diag_ind) = P_SiO2%remin

      part_diags(:, P_iron_FLUX_IN_diag_ind) = P_iron%sflux_in + P_iron%hflux_in
      part_diags(:, P_iron_PROD_diag_ind) = P_iron%prod
      part_diags(:, P_iron_REMIN_diag_ind) = P_iron%remin

      part_diags(:, calcToSed_diag_ind) = P_CaCO3%sed_loss
      part_diags(:, bsiToSed_diag_ind) = P_SiO2%sed_loss
      part_diags(:, pocToSed_diag_ind) = POC%sed_loss
      part_diags(:, SedDenitrif_diag_ind) = sed_denitrif * delta_z
      part_diags(:, OtherRemin_diag_ind) = other_remin * delta_z
      part_diags(:, ponToSed_diag_ind) = (POC%sed_loss * Q)
      part_diags(:, popToSed_diag_ind) = (POC%sed_loss * Qp_zoo_pom)
      part_diags(:, dustToSed_diag_ind) = dust%sed_loss
      part_diags(:, pfeToSed_diag_ind) = P_iron%sed_loss
    end associate

  end subroutine store_diagnostics_particulates

  !-----------------------------------------------------------------------

  subroutine store_diagnostics_oxygen(marbl_domain, marbl_gcm_state, column_zt,                &
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
    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d                                &
              )
      min_ind = minloc(column_o2(1:marbl_domain%kmt), dim=1)
      if ((min_ind.gt.0).and.(min_ind.le.marbl_domain%kmt)) then
        diags_2d(:, O2_ZMIN_diag_ind) = column_o2(min_ind)
        diags_2d(:, O2_ZMIN_DEPTH_diag_ind) = column_zt(min_ind)
      else
        diags_2d(:, O2_ZMIN_diag_ind) = column_o2(1)
        diags_2d(:, O2_ZMIN_DEPTH_diag_ind) = column_zt(1)
      end if

      diags_3d(:, O2_PRODUCTION_diag_ind) = o2_production
      diags_3d(:, O2_CONSUMPTION_diag_ind) = o2_consumption

      diags_3d(:, AOU_diag_ind) = -column_o2
      do k=1,marbl_domain%kmt
        if (marbl_domain%land_mask) then
          diags_3d(k, AOU_diag_ind) =                                         &
        O2SAT_scalar(marbl_gcm_state%temperature(k), marbl_gcm_state%salinity(k)) - &
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

    associate(diags_3d => marbl_diags%diags_3d)
      diags_3d(:, PAR_avg_diag_ind) = PAR%avg
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

    associate(zoo_diags => marbl_diags%zoo_diags)
      do n = 1, zooplankton_cnt
        zoo_diags(:, zoo_loss_diag_ind, n)       = zooplankton_secondary_species(n,:)%zoo_loss
        zoo_diags(:, zoo_loss_poc_diag_ind, n)   = zooplankton_secondary_species(n,:)%zoo_loss_poc
        zoo_diags(:, zoo_loss_doc_diag_ind, n)   = zooplankton_secondary_species(n,:)%zoo_loss_doc
        zoo_diags(:, zoo_graze_diag_ind, n)      = zooplankton_secondary_species(n,:)%zoo_graze
        zoo_diags(:, zoo_graze_poc_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_poc
        zoo_diags(:, zoo_graze_doc_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_doc
        zoo_diags(:, zoo_graze_zoo_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_zoo
        zoo_diags(:, x_graze_zoo_diag_ind, n)    = zooplankton_secondary_species(n,:)%x_graze_zoo
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

    associate(diags_3d => marbl_diags%diags_3d)
      diags_3d(:, DOC_prod_diag_ind)         = dissolved_organic_matter%DOC_prod
      diags_3d(:, DOC_remin_diag_ind)        = dissolved_organic_matter%DOC_remin
      diags_3d(:, DON_prod_diag_ind)         = dissolved_organic_matter%DON_prod
      diags_3d(:, DON_remin_diag_ind)        = dissolved_organic_matter%DON_remin
      diags_3d(:, DOP_prod_diag_ind)         = dissolved_organic_matter%DOP_prod
      diags_3d(:, DOP_remin_diag_ind)        = dissolved_organic_matter%DOP_remin
      diags_3d(:, DOFe_prod_diag_ind)        = dissolved_organic_matter%DOFe_prod
      diags_3d(:, DOFe_remin_diag_ind)       = dissolved_organic_matter%DOFe_remin
      diags_3d(:, Fe_scavenge_diag_ind)      = Fe_scavenge
      diags_3d(:, Fe_scavenge_rate_diag_ind) = Fe_scavenge_rate
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

    associate(diags_3d => marbl_diags%diags_3d)
      ! vertical integrals
      diags_3d(:, Jint_Ctot_diag_ind) = c0
      diags_3d(:, Jint_100m_Ctot_diag_ind) = c0

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
        diags_3d(k, Jint_Ctot_diag_ind) = delta_z(k) * work1 +                &
                                      (POC%sed_loss(k) + P_CaCO3%sed_loss(k))

        if (ztop < 100.0e2_r8) then
          diags_3d(k, Jint_100m_Ctot_diag_ind) =                              &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_3d(k, Jint_100m_Ctot_diag_ind) =                              &
                                       diags_3d(k, Jint_100m_Ctot_diag_ind) + &
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

    associate(diags_3d => marbl_diags%diags_3d)
      ! vertical integrals
      diags_3d(:, Jint_Ntot_diag_ind) = c0
      diags_3d(:, Jint_100m_Ntot_diag_ind) = c0

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
        diags_3d(k, Jint_Ntot_diag_ind) = delta_z(k) * work1 +                &
                                          POC%sed_loss(k) * Q

        if (ztop < 100.0e2_r8) then
            diags_3d(k, Jint_100m_Ntot_diag_ind) =                            &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
            diags_3d(k, Jint_100m_Ntot_diag_ind) =                            &
                   diags_3d(k, Jint_100m_Ntot_diag_ind) + POC%sed_loss(k) * Q
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

    associate(diags_3d => marbl_diags%diags_3d)
      ! vertical integrals
      diags_3d(:, Jint_Ptot_diag_ind) = c0
      diags_3d(:, Jint_100m_Ptot_diag_ind) = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(po4_ind,k) + dtracer(dop_ind,k) + dtracer(dopr_ind,k)
        do n = 1, zooplankton_cnt
          work1 = work1 + Qp_zoo_pom * dtracer(zooplankton(n)%C_ind,k)
        end do
        do n = 1, autotroph_cnt
          work1 = work1 + autotrophs(n)%Qp * dtracer(autotrophs(n)%C_ind,k)
        end do
        diags_3d(k, Jint_Ptot_diag_ind) = delta_z(k) * work1 +                &
                                                 POC%sed_loss(k) * Qp_zoo_pom

        if (ztop < 100.0e2_r8) then
          diags_3d(k, Jint_100m_Ptot_diag_ind) =                              &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_3d(k, Jint_100m_Ptot_diag_ind) =                              &
                                       diags_3d(k, Jint_100m_Ptot_diag_ind) + &
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

    associate(diags_3d => marbl_diags%diags_3d)
      ! vertical integrals
      diags_3d(:, Jint_Sitot_diag_ind) = c0
      diags_3d(:, Jint_100m_Sitot_diag_ind) = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(sio3_ind,k)
        do n = 1, autotroph_cnt
          if (autotrophs(n)%Si_ind > 0) then
            work1 = work1 + dtracer(autotrophs(n)%Si_ind,k)
          end if
        end do
        diags_3d(k, Jint_Sitot_diag_ind) = delta_z(k) * work1 +               &
                                                           P_SiO2%sed_loss(k)

        if (ztop < 100.0e2_r8) then
          diags_3d(k, Jint_100m_Sitot_diag_ind) =                             &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_3d(k, Jint_100m_Sitot_diag_ind) =                             &
                                      diags_3d(k, Jint_100m_Sitot_diag_ind) + &
                                                           P_SiO2%sed_loss(k)
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_silicon_fluxes

  !-----------------------------------------------------------------------

end module ecosys_diagnostics_mod
