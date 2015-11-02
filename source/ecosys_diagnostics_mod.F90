module ecosys_diagnostics_mod

  use ecosys_constants, only : ecosys_tracer_cnt

  use marbl_kinds_mod, only : c0
  use marbl_kinds_mod, only : c1
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

  integer(int_kind), parameter :: ecosys_diag_cnt_2d = 14
  integer (int_kind), parameter ::  &
      zsatcalc_diag_ind            =  1, &
      zsatarag_diag_ind            =  2, &
      O2_ZMIN_diag_ind             =  3, &
      O2_ZMIN_DEPTH_diag_ind       =  4, &
      photoC_TOT_zint_diag_ind     =  5, &
      photoC_NO3_TOT_zint_diag_ind =  6, &
      Jint_Ctot_diag_ind           =  7, &
      Jint_100m_Ctot_diag_ind      =  8, &
      Jint_Ntot_diag_ind           =  9, &
      Jint_100m_Ntot_diag_ind      = 10, &
      Jint_Ptot_diag_ind           = 11, &
      Jint_100m_Ptot_diag_ind      = 12, &
      Jint_Sitot_diag_ind          = 13, &
      Jint_100m_Sitot_diag_ind     = 14

  integer(int_kind), parameter :: ecosys_diag_cnt_3d = 29
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
      Fe_scavenge_rate_diag_ind    = 29

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

  subroutine store_diagnostics_carbonate(marbl_domain, carbonate, zt, zw,     &
                                         marbl_diags)

    type(carbonate_type), dimension(:), intent(in) :: carbonate
    type(marbl_column_domain_type), intent(in) :: marbl_domain
    real(r8), dimension(:), intent(in) :: zt, zw
    type(marbl_diagnostics_type), intent(inout) :: marbl_diags

    real(r8) :: zsat_calcite, zsat_aragonite

    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d,                               &
              CO3 => carbonate%CO3,                                           &
              CO3_sat_calcite => carbonate%CO3_sat_calcite,                   &
              CO3_sat_aragonite => carbonate%CO3_sat_aragonite                &
              )
      ! Find depth where CO3 = CO3_sat_calcite or CO3_sat_argonite
      diags_2d(zsatcalc_diag_ind) = compute_saturation_depth(marbl_domain,    &
                                                           zt, zw, CO3,       &
                                                           CO3_sat_calcite)
      diags_2d(zsatarag_diag_ind) = compute_saturation_depth(marbl_domain,    &
                                                           zt, zw, CO3,       &
                                                           CO3_sat_aragonite)

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

    real(kind=r8), dimension(2), intent(in) :: x,y
    real(kind=r8) :: linear_root

    real(kind=r8) :: m_inv

    if (y(1)*y(2).gt.c0) then
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

    associate(&
              auto_diags_2d => marbl_diags%auto_diags_2d,&
              auto_diags_3d => marbl_diags%auto_diags_3d &
             )
      do n = 1, autotroph_cnt
         auto_diags_3d(:, N_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VNtot
         auto_diags_3d(:, Fe_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VFe
         auto_diags_3d(:, P_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VPtot

         if (autotrophs(n)%kSiO3 > c0) then
            auto_diags_3d(:, SiO3_lim_diag_ind, n) = autotroph_secondary_species(n,:)%VSiO3
         end if

         auto_diags_3d(:, light_lim_diag_ind, n) = autotroph_secondary_species(n,:)%light_lim
         auto_diags_3d(:, photoNO3_diag_ind, n) = autotroph_secondary_species(n,:)%NO3_V
         auto_diags_3d(:, photoNH4_diag_ind, n) = autotroph_secondary_species(n,:)%NH4_V
         auto_diags_3d(:, PO4_uptake_diag_ind, n) = autotroph_secondary_species(n,:)%PO4_V
         auto_diags_3d(:, DOP_uptake_diag_ind, n) = autotroph_secondary_species(n,:)%DOP_V
         auto_diags_3d(:, photoFE_diag_ind, n) = autotroph_secondary_species(n,:)%photoFe

         if (autotrophs(n)%Si_ind > 0) then
           auto_diags_3d(:, bSi_form_diag_ind, n) = autotroph_secondary_species(n,:)%photoSi
         endif

         auto_diags_3d(:, CaCO3_form_diag_ind, n) = autotroph_secondary_species(n,:)%CaCO3_PROD
         auto_diags_3d(:, Nfix_diag_ind, n) = autotroph_secondary_species(n,:)%Nfix
         auto_diags_3d(:, auto_graze_diag_ind, n)      = autotroph_secondary_species(n,:)%auto_graze
         auto_diags_3d(:, auto_graze_poc_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_poc
         auto_diags_3d(:, auto_graze_doc_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_doc
         auto_diags_3d(:, auto_graze_zoo_diag_ind, n)  = autotroph_secondary_species(n,:)%auto_graze_zoo
         auto_diags_3d(:, auto_loss_diag_ind, n)       = autotroph_secondary_species(n,:)%auto_loss
         auto_diags_3d(:, auto_loss_poc_diag_ind, n)   = autotroph_secondary_species(n,:)%auto_loss_poc
         auto_diags_3d(:, auto_loss_doc_diag_ind, n)   = autotroph_secondary_species(n,:)%auto_loss_doc
         auto_diags_3d(:, auto_agg_diag_ind, n)        = autotroph_secondary_species(n,:)%auto_agg
         auto_diags_3d(:, photoC_diag_ind, n)          = autotroph_secondary_species(n,:)%photoC

         auto_diags_3d(:, photoC_NO3_diag_ind, n) = c0
         where (autotroph_secondary_species(n,:)%VNtot > c0)
           auto_diags_3d(:, photoC_NO3_diag_ind, n) =                         &
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
                        auto_diags_2d(photoC_NO3_zint_diag_ind, n) +          &
                        delta_z(k) * auto_diags_3d(k, photoC_NO3_diag_ind, n)
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

    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d                                &
              )
      diags_3d(:, auto_graze_TOT_diag_ind) =                                  &
                           sum(autotroph_secondary_species%auto_graze, dim=1)
      diags_3d(:, photoC_TOT_diag_ind) =                                      &
                           sum(autotroph_secondary_species%photoC, dim=1)

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

      diags_2d(photoC_TOT_zint_diag_ind) = sum(                               &
                 delta_z * sum(autotroph_secondary_species%photoC, dim=1))
      diags_2d(photoC_NO3_TOT_zint_diag_ind) = sum(                           &
                              delta_z * diags_3d(:, photoC_NO3_TOT_diag_ind))
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
              part_diags_2d => marbl_diags%part_diags_2d,&
              part_diags_3d => marbl_diags%part_diags_3d &
              )
      part_diags_3d(:, POC_FLUX_IN_diag_ind) = POC%sflux_in + POC%hflux_in
      part_diags_3d(:, POC_PROD_diag_ind) = POC%prod
      part_diags_3d(:, POC_REMIN_diag_ind) = POC%remin

      part_diags_3d(:, CaCO3_FLUX_IN_diag_ind) = P_CaCO3%sflux_in + P_CaCO3%hflux_in
      part_diags_3d(:, CaCO3_PROD_diag_ind) = P_CaCO3%prod
      part_diags_3d(:, CaCO3_REMIN_diag_ind) = P_CaCO3%remin

      part_diags_3d(:, SiO2_FLUX_IN_diag_ind) = P_SiO2%sflux_in + P_SiO2%hflux_in
      part_diags_3d(:, SiO2_PROD_diag_ind) = P_SiO2%prod
      part_diags_3d(:, SiO2_REMIN_diag_ind) = P_SiO2%remin

      part_diags_3d(:, dust_FLUX_IN_diag_ind) = dust%sflux_in + dust%hflux_in
      part_diags_3d(:, dust_REMIN_diag_ind) = P_SiO2%remin

      part_diags_3d(:, P_iron_FLUX_IN_diag_ind) = P_iron%sflux_in + P_iron%hflux_in
      part_diags_3d(:, P_iron_PROD_diag_ind) = P_iron%prod
      part_diags_3d(:, P_iron_REMIN_diag_ind) = P_iron%remin

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
    associate(                                                                &
              diags_2d => marbl_diags%diags_2d,                               &
              diags_3d => marbl_diags%diags_3d                                &
              )
      min_ind = minloc(column_o2(1:marbl_domain%kmt), dim=1)
      diags_2d(O2_ZMIN_diag_ind) = column_o2(min_ind)
      diags_2d(O2_ZMIN_DEPTH_diag_ind) = column_zt(min_ind)

      diags_3d(:, O2_PRODUCTION_diag_ind) = o2_production
      diags_3d(:, O2_CONSUMPTION_diag_ind) = o2_consumption

      diags_3d(:, AOU_diag_ind) = -column_o2
      do k=1,marbl_domain%kmt
        if (marbl_domain%land_mask) then
          diags_3d(k, AOU_diag_ind) =                                         &
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

    associate(zoo_diags_3d => marbl_diags%zoo_diags_3d)
      do n = 1, zooplankton_cnt
        zoo_diags_3d(:, zoo_loss_diag_ind, n)       = zooplankton_secondary_species(n,:)%zoo_loss
        zoo_diags_3d(:, zoo_loss_poc_diag_ind, n)   = zooplankton_secondary_species(n,:)%zoo_loss_poc
        zoo_diags_3d(:, zoo_loss_doc_diag_ind, n)   = zooplankton_secondary_species(n,:)%zoo_loss_doc
        zoo_diags_3d(:, zoo_graze_diag_ind, n)      = zooplankton_secondary_species(n,:)%zoo_graze
        zoo_diags_3d(:, zoo_graze_poc_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_poc
        zoo_diags_3d(:, zoo_graze_doc_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_doc
        zoo_diags_3d(:, zoo_graze_zoo_diag_ind, n)  = zooplankton_secondary_species(n,:)%zoo_graze_zoo
        zoo_diags_3d(:, x_graze_zoo_diag_ind, n)    = zooplankton_secondary_species(n,:)%x_graze_zoo
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

    associate(diags_2d => marbl_diags%diags_2d)
      ! vertical integrals
      diags_2d(Jint_Ctot_diag_ind) = c0
      diags_2d(Jint_100m_Ctot_diag_ind) = c0

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
        diags_2d(Jint_Ctot_diag_ind) = diags_2d(Jint_Ctot_diag_ind) +         &
                 delta_z(k) * work1 + (POC%sed_loss(k) + P_CaCO3%sed_loss(k))

        if (ztop < 100.0e2_r8) then
          diags_2d(Jint_100m_Ctot_diag_ind) =                                 &
                                   diags_2d(Jint_100m_Ctot_diag_ind) +        &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_2d(Jint_100m_Ctot_diag_ind) =                                 &
                                      diags_2d(Jint_100m_Ctot_diag_ind) +     &
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

    associate(diags_2d => marbl_diags%diags_2d)
      ! vertical integrals
      diags_2d(Jint_Ntot_diag_ind) = c0
      diags_2d(Jint_100m_Ntot_diag_ind) = c0

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
        diags_2d(Jint_Ntot_diag_ind) = diags_2d(Jint_Ntot_diag_ind) +         &
                                     delta_z(k) * work1 + POC%sed_loss(k) * Q

        if (ztop < 100.0e2_r8) then
            diags_2d(Jint_100m_Ntot_diag_ind) =                               &
                                   diags_2d(Jint_100m_Ntot_diag_ind) +        &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
            diags_2d(Jint_100m_Ntot_diag_ind) =                               &
                      diags_2d(Jint_100m_Ntot_diag_ind) + POC%sed_loss(k) * Q
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

    associate(diags_2d => marbl_diags%diags_2d)
      ! vertical integrals
      diags_2d(Jint_Ptot_diag_ind) = c0
      diags_2d(Jint_100m_Ptot_diag_ind) = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(po4_ind,k) + dtracer(dop_ind,k) + dtracer(dopr_ind,k)
        do n = 1, zooplankton_cnt
          work1 = work1 + Qp_zoo_pom * dtracer(zooplankton(n)%C_ind,k)
        end do
        do n = 1, autotroph_cnt
          work1 = work1 + autotrophs(n)%Qp * dtracer(autotrophs(n)%C_ind,k)
        end do
        diags_2d(Jint_Ptot_diag_ind) = diags_2d(Jint_Ptot_diag_ind) +         &
                            delta_z(k) * work1 + POC%sed_loss(k) * Qp_zoo_pom

        if (ztop < 100.0e2_r8) then
          diags_2d(Jint_100m_Ptot_diag_ind) =                                 &
                                   diags_2d(Jint_100m_Ptot_diag_ind) +        &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_2d(Jint_100m_Ptot_diag_ind) =                                 &
             diags_2d(Jint_100m_Ptot_diag_ind) + POC%sed_loss(k) * Qp_zoo_pom
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

    associate(diags_2d => marbl_diags%diags_2d)
      ! vertical integrals
      diags_2d(Jint_Sitot_diag_ind) = c0
      diags_2d(Jint_100m_Sitot_diag_ind) = c0

      ztop = c0
      do k = 1,marbl_domain%kmt
        work1 = dtracer(sio3_ind,k)
        do n = 1, autotroph_cnt
          if (autotrophs(n)%Si_ind > 0) then
            work1 = work1 + dtracer(autotrophs(n)%Si_ind,k)
          end if
        end do
        diags_2d(Jint_Sitot_diag_ind) = diags_2d(Jint_Sitot_diag_ind) +       &
                                      delta_z(k) * work1 + P_SiO2%sed_loss(k)

        if (ztop < 100.0e2_r8) then
          diags_2d(Jint_100m_Sitot_diag_ind) =                                &
                                   diags_2d(Jint_100m_Sitot_diag_ind) +       &
                                   min(100.0e2_r8 - ztop, delta_z(k)) * work1
        end if
        if (zw(k).le.100.0e2_r8) then
          diags_2d(Jint_100m_Sitot_diag_ind) =                                &
                      diags_2d(Jint_100m_Sitot_diag_ind) + P_SiO2%sed_loss(k)
        end if
        ztop = zw(k)
      end do
    end associate

  end subroutine store_diagnostics_silicon_fluxes

  !-----------------------------------------------------------------------

end module ecosys_diagnostics_mod
