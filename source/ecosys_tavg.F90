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

  use shr_sys_mod, only : shr_sys_abort

  use ecosys_constants, only : ecosys_tracer_cnt
  use ecosys_diagnostics_mod, only : marbl_diag_ind

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

  integer (int_kind), allocatable, dimension(:) :: tavg_ecosys

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

  integer (int_kind) :: tavg_POC_ACCUM      ! tavg id for poc accumulation

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

    character(char_len) :: err_msg, gloc, coords, sname
    integer(int_kind) :: n, ndims
    logical :: define_field

!-----------------------------------------------------------------------
!   Define tavg fields for surface flux terms
!-----------------------------------------------------------------------

    allocate(tavg_ecosys(marbl_diag_ind%count))
    associate(                                                                &
              diags         => marbl_diags%diags(:),                          &
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
!   Define 2D tavg fields
!-----------------------------------------------------------------------

    do n=1,marbl_diag_ind%count
      if (trim(diags(n)%vertical_grid).eq.'none') then
        ndims = 2
        gloc = '2110'
        coords = 'TLONG TLAT time'
      else
        ndims = 3
        if (trim(diags(n)%vertical_grid).eq.'layer_avg') then
          if (diags(n)%ltruncated_vertical_extent) then
            gloc = '3114'
            coords = 'TLONG TLAT z_t_150m time'
          else
            gloc = '3111'
            coords = 'TLONG TLAT z_t time'
          end if
        elseif (trim(diags(n)%vertical_grid).eq.'layer_iface') then
            gloc = '3113'
            coords = 'TLONG TLAT z_w_bot time'
        else
          write(err_msg,*) "'", trim(diags(n)%vertical_grid),                 &
                           "' is not a valid vertical grid"
          call shr_sys_abort(err_msg)
        end if
      end if
      call define_tavg_field(tavg_ecosys(n),trim(diags(n)%short_name),        &
                             ndims, long_name=trim(diags(n)%long_name),       &
                             units=trim(diags(n)%units), grid_loc=gloc,       &
                             coordinates=coords)
    end do

    call define_tavg_field(tavg_tot_CaCO3_form_zint, 'CaCO3_form_zint', 2, &
                           long_name='Total CaCO3 Formation Vertical Integral', &
                           units='mmol/m^3 cm/s', grid_loc='2110', &
                           coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!  Define 3D tavg fields for autotrophs
!-----------------------------------------------------------------------

!       if (autotrophs(auto_ind)%Si_ind > 0) then
!          sname = trim(autotrophs(auto_ind)%sname) // 'bSi_form'
!          call define_tavg_field(tavg_auto_3d(bSi_form_diag_ind,auto_ind), sname, 3, &
!                                 long_name=trim(autotrophs(auto_ind)%lname) // ' Si Uptake', &
!                                 units='mmol/m^3/s', grid_loc='3114', &
!                                 coordinates='TLONG TLAT z_t_150m time')
!       endif

!       if (autotrophs(auto_ind)%CaCO3_ind > 0) then
!          sname = trim(autotrophs(auto_ind)%sname) // '_CaCO3_form'
!          call define_tavg_field(tavg_auto_3d(CaCO3_form_diag_ind,auto_ind), sname, 3, &
!                                 long_name=trim(autotrophs(auto_ind)%lname) // ' CaCO3 Formation', &
!                                 units='mmol/m^3/s', grid_loc='3114', &
!                                 coordinates='TLONG TLAT z_t_150m time')
!       endif

!       if (autotrophs(auto_ind)%Nfixer) then
!          call define_tavg_field(tavg_auto_3d(Nfix_diag_ind,auto_ind), &
!                                 trim(autotrophs(auto_ind)%sname) // '_Nfix', 3, &
!                                 long_name=trim(autotrophs(auto_ind)%lname) // ' N Fixation', &
!                                 units='mmol/m^3/s', grid_loc='3114',   &
!                                 coordinates='TLONG TLAT z_t_150m time')
!       endif
!    end do

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

    call define_tavg_field(tavg_POC_ACCUM,'POC_ACCUM',3,                &
                           long_name='POC Accumulation',                &
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
              diags         => marbl_diagnostics%diags(:),                    &
              restore_diags => marbl_diagnostics%restore_diags(:)             &
             )

    ! Accumulate general diagnostics
    do n=1,marbl_diag_ind%count
      if (trim(diags(n)%vertical_grid).eq.'none') then
        call accumulate_tavg_field(diags(n)%field_2d, tavg_ecosys(n),         &
                                   bid, i, c)
      else
        call accumulate_tavg_field(diags(n)%field_3d(:), tavg_ecosys(n),      &
                                   bid, i, c)
      end if
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


