! Stays in POP (ecosys driver)
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

  use kinds_mod             , only : r8
  use kinds_mod             , only : int_kind
  use kinds_mod             , only : char_len

  use blocks                , only : nx_block
  use blocks                , only : ny_block

  use domain                , only : nblocks_clinic
  use domain_size           , only : km

  use tavg                  , only : define_tavg_field
  use tavg                  , only : accumulate_tavg_field

  use shr_sys_mod           , only : shr_sys_abort

  use ecosys_constants      , only : ecosys_tracer_cnt

  use marbl_share_mod       , only : autotrophs
  use marbl_share_mod       , only : zooplankton
  use marbl_share_mod       , only : autotroph_cnt
  use marbl_share_mod       , only : zooplankton_cnt

  use marbl_interface_types , only : marbl_diagnostics_type
  use marbl_interface_types , only : max_interior_diags
  use marbl_interface_types , only : max_forcing_diags

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_tavg_init
  public :: ecosys_tavg_accumulate
  public :: ecosys_tavg_accumulate_flux

  !-----------------------------------------------------------------------
  !  define tavg id for interior diagnostics, diagnostics related to
  !  restoring, surface forcing diagnostics, and duplicate surface forcing
  !  variables
  !-----------------------------------------------------------------------

  integer (int_kind), dimension(max_forcing_diags ) :: tavg_ids_forcing
  integer (int_kind), dimension(max_interior_diags) :: tavg_ids_interior
  integer (int_kind), dimension(ecosys_tracer_cnt ) :: tavg_ids_restore
  integer (int_kind) ::      &
      tavg_ECOSYS_IFRAC_2,   &! ice fraction duplicate
      tavg_ECOSYS_XKW_2,     &! xkw duplicate
      tavg_O2_GAS_FLUX_2,    &! O2 flux duplicate
      tavg_DpCO2_2,          &! delta pco2 duplicate
      tavg_DIC_GAS_FLUX_2     ! dic flux duplicate

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_tavg_init(marbl_interior_diags, marbl_restore_diags, marbl_forcing_diags)

    ! !DESCRIPTION:
    !  call define_tavg_field for all tavg fields
    !
    type(marbl_diagnostics_type)    , intent(in)  :: marbl_interior_diags
    type(marbl_diagnostics_type)    , intent(in)  :: marbl_restore_diags
    type(marbl_diagnostics_type)    , intent(in)  :: marbl_forcing_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'

    !-----------------------------------------------------------------------
    !  Define tavg fields for MARBL diagnostics
    !-----------------------------------------------------------------------

    call ecosys_tavg_define_from_diag(marbl_diags=marbl_interior_diags,       &
                                      tavg_ids=tavg_ids_interior)
    call ecosys_tavg_define_from_diag(marbl_diags=marbl_restore_diags,        &
                                      tavg_ids=tavg_ids_restore)
    call ecosys_tavg_define_from_diag(marbl_diags=marbl_forcing_diags,        &
                                      tavg_ids=tavg_ids_forcing)

    call define_tavg_field(tavg_ECOSYS_IFRAC_2,'ECOSYS_IFRAC_2',2,      &
                           long_name='Ice Fraction for ecosys fluxes',  &
                           units='fraction', grid_loc='2110',           &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_ECOSYS_XKW_2,'ECOSYS_XKW_2',2,          &
                           long_name='XKW for ecosys fluxes',           &
                           units='cm/s', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_O2_GAS_FLUX_2,'STF_O2_2',2,             &
                           long_name='Dissolved Oxygen Surface Flux',   &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_DpCO2_2,'DpCO2_2',2,                    &
                           long_name='D pCO2',                          &
                           units='ppmv', grid_loc='2110',               &
                           coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_DIC_GAS_FLUX_2,'FG_CO2_2',2,            &
                           long_name='DIC Surface Gas Flux',            &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

  end subroutine ecosys_tavg_init

  !***********************************************************************

  subroutine ecosys_tavg_accumulate(i, c, bid, marbl_interior_diags, marbl_restore_diags)

    integer, intent(in) :: i, c, bid ! column indices and block index
    type(marbl_diagnostics_type), intent(in) :: marbl_interior_diags
    type(marbl_diagnostics_type), intent(in) :: marbl_restore_diags

    integer :: n

    call ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_interior_diags, tavg_ids_interior)
    call ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_restore_diags , tavg_ids_restore)

  end subroutine ecosys_tavg_accumulate

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_diags, tavg_ids)

    ! Accumulate diagnostics

    integer                      , intent(in) :: i, c, bid ! column indices and block index
    type(marbl_diagnostics_type) , intent(in) :: marbl_diags
    integer, dimension(:)        , intent(in) :: tavg_ids

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer :: n
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))
      do n=1,marbl_diags%diag_cnt
        if (trim(diags(n)%vertical_grid).eq.'none') then
          call accumulate_tavg_field(diags(n)%field_2d, tavg_ids(n), bid, i, c)
        else
          call accumulate_tavg_field(diags(n)%field_3d(:), tavg_ids(n), bid, i, c)
        end if
      end do
    end associate

  end subroutine ecosys_tavg_accumulate_from_diag

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_flux(flux_diags, marbl_forcing_diags)

    ! Compute diagnostics for surface fluxes

    use ecosys_diagnostics_mod, only : ind => marbl_forcing_diag_ind

    implicit none
    real (r8)                    , intent(in)    :: flux_diags (:, :, :, :)
    type(marbl_diagnostics_type) , intent(inout) :: marbl_forcing_diags(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer :: i, iblock
    integer :: nblocks_clinic
    !-----------------------------------------------------------------------

    nblocks_clinic = size(flux_diags,4)

    !$OMP PARALLEL DO PRIVATE(iblock)
    do iblock=1,nblocks_clinic
       do i = 1,marbl_forcing_diags(iblock)%diag_cnt
          call accumulate_tavg_field(FLUX_DIAGS(:,:,i,iblock), tavg_ids_forcing(i), iblock, i)
       end do
       call accumulate_tavg_field(FLUX_diags(:,:,ind%ECOSYS_IFRAC,iblock),    &
            tavg_ECOSYS_IFRAC_2, iblock, i)
       call accumulate_tavg_field(FLUX_diags(:,:,ind%ECOSYS_XKW,iblock),      &
            tavg_ECOSYS_XKW_2, iblock, i)
       call accumulate_tavg_field(FLUX_diags(:,:,ind%O2_GAS_FLUX,iblock),     &
            tavg_O2_GAS_FLUX_2, iblock, i)
       call accumulate_tavg_field(FLUX_diags(:,:,ind%DpCO2,iblock),           &
            tavg_DpCO2_2, iblock, i)
       call accumulate_tavg_field(FLUX_diags(:,:,ind%DIC_GAS_FLUX,iblock),    &
            tavg_DIC_GAS_FLUX_2, iblock, i)
    end do
    !$OMP END PARALLEL DO
    
  end subroutine ecosys_tavg_accumulate_flux

  !***********************************************************************

  subroutine ecosys_tavg_define_from_diag(marbl_diags, tavg_ids)

    type(marbl_diagnostics_type)    , intent(in)    :: marbl_diags
    integer(int_kind), dimension(:) , intent(inout) :: tavg_ids

    character(char_len) :: err_msg, gloc, coords
    integer :: n, ndims

    associate(diags => marbl_diags%diags(:))

      do n=1,marbl_diags%diag_cnt
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
               write(err_msg,*) "'", trim(diags(n)%vertical_grid), &
                    "' is not a valid vertical grid"
               call shr_sys_abort(err_msg)
            end if
         end if
         call define_tavg_field(tavg_ids(n),             &
              trim(diags(n)%short_name),                 &
              ndims,                                     &
              long_name=trim(diags(n)%long_name),        &
              units=trim(diags(n)%units), grid_loc=gloc, &
              coordinates=coords)
      end do
    end associate
    
  end subroutine ecosys_tavg_define_from_diag

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


