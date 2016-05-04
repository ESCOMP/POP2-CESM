! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
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
  use tavg                  , only : define_tavg_field
  use tavg                  , only : accumulate_tavg_field
  use shr_sys_mod           , only : shr_sys_abort
  use marbl_interface       , only : marbl_interface_class
  use marbl_diagnostics_mod , only : marbl_diagnostics_type

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

  integer (int_kind), allocatable :: tavg_ids_interior_forcing(:)
  integer (int_kind), allocatable :: tavg_ids_interior_restore(:)
  integer (int_kind), allocatable :: tavg_ids_surface_forcing(:)

  integer (int_kind) :: tavg_ECOSYS_IFRAC_2 ! ice fraction duplicate
  integer (int_kind) :: tavg_ECOSYS_XKW_2   ! xkw duplicate
  integer (int_kind) :: tavg_O2_GAS_FLUX_2  ! O2 flux duplicate
  integer (int_kind) :: tavg_DpCO2_2        ! delta pco2 duplicate
  integer (int_kind) :: tavg_DIC_GAS_FLUX_2 ! dic flux duplicate

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_tavg_init(marbl_instance)

    ! !DESCRIPTION:
    !  call define_tavg_field for all tavg fields

    implicit none

    type(marbl_interface_class)  , intent(in) :: marbl_instance

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  allocate memory for module variables
    !-----------------------------------------------------------------------

    associate(&
         cnt_interior_forcing => marbl_instance%interior_forcing_diags%diag_cnt, &
         cnt_interior_restore => marbl_instance%interior_restore_diags%diag_cnt, &
         cnt_surface_forcing  => marbl_instance%surface_forcing_diags%diag_cnt   &
         )

    allocate(tavg_ids_interior_forcing(cnt_interior_forcing))
    allocate(tavg_ids_interior_restore(cnt_interior_restore))
    allocate(tavg_ids_surface_forcing(cnt_surface_forcing))

    end associate

    !-----------------------------------------------------------------------
    !  Define tavg fields for MARBL diagnostics
    !-----------------------------------------------------------------------

    associate(&
         interior_forcing_diags => marbl_instance%interior_forcing_diags, &
         interior_restore_diags => marbl_instance%interior_restore_diags, &
         surface_forcing_diags => marbl_instance%surface_forcing_diags    &
         )

    call ecosys_tavg_define_from_diag(marbl_diags=interior_forcing_diags, &
         tavg_ids=tavg_ids_interior_forcing)

    call ecosys_tavg_define_from_diag(marbl_diags=interior_restore_diags,  &
         tavg_ids=tavg_ids_interior_restore)
    
    call ecosys_tavg_define_from_diag(marbl_diags=surface_forcing_diags,  &
         tavg_ids=tavg_ids_surface_forcing)

    end associate

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

  subroutine ecosys_tavg_accumulate(i, c, bid, &
       marbl_interior_forcing_diags,           &
       marbl_interior_restore_diags,           &
       marbl_surface_forcing_diags)

    implicit none

    integer , dimension(:), intent(in) :: i, c ! column indices
    integer ,               intent(in) :: bid ! block index

    type(marbl_diagnostics_type), optional, intent(in) :: marbl_interior_forcing_diags
    type(marbl_diagnostics_type), optional, intent(in) :: marbl_interior_restore_diags
    type(marbl_diagnostics_type), optional, intent(in) :: marbl_surface_forcing_diags
    !-----------------------------------------------------------------------

    if (present(marbl_interior_forcing_diags)) then
      call ecosys_tavg_accumulate_from_diag(i, c, bid, &
           marbl_diags = marbl_interior_forcing_diags,  &
           tavg_ids = tavg_ids_interior_forcing, &
           num_elements =marbl_interior_forcing_diags%num_elements)
    end if

    if (present(marbl_interior_restore_diags)) then
      call ecosys_tavg_accumulate_from_diag(i, c, bid, &
           marbl_diags = marbl_interior_restore_diags,   &
           tavg_ids = tavg_ids_interior_restore, &
           num_elements = marbl_interior_restore_diags%num_elements)
    end if

    if (present(marbl_surface_forcing_diags)) then
      call ecosys_tavg_accumulate_from_diag(i, c, bid, &
           marbl_diags = marbl_surface_forcing_diags,   &
           tavg_ids = tavg_ids_surface_forcing, &
           num_elements = marbl_surface_forcing_diags%num_elements)
    end if

  end subroutine ecosys_tavg_accumulate

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_diags, tavg_ids, num_elements)

    ! Accumulate diagnostics

    implicit none

    integer, dimension(:)        , intent(in) :: i, c ! column indices
    integer                      , intent(in) :: bid ! block index
    type(marbl_diagnostics_type) , intent(in) :: marbl_diags
    integer, dimension(:)        , intent(in) :: tavg_ids
    integer                      , intent(in) :: num_elements

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: n, ne
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))

    do n=1,marbl_diags%diag_cnt
       do ne = 1,num_elements
          if (trim(diags(n)%vertical_grid).eq.'none') then
             call accumulate_tavg_field(diags(n)%field_2d(ne)  , tavg_ids(n), bid, i(ne), c(ne))
          else
             call accumulate_tavg_field(diags(n)%field_3d(:,ne), tavg_ids(n), bid, i(ne), c(ne))
          end if
       end do
    end do

    end associate

  end subroutine ecosys_tavg_accumulate_from_diag

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_flux(flux_diags, marbl_instances)

    ! Compute diagnostics for surface fluxes

    use marbl_diagnostics_mod, only : ind => marbl_surface_forcing_diag_ind

    implicit none

    real (r8)                  , intent(in) :: flux_diags (:, :, :, :)
    type(marbl_interface_class), intent(in) :: marbl_instances(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: i, iblock
    integer :: nblocks_clinic
    !-----------------------------------------------------------------------

    nblocks_clinic = size(flux_diags,4)

    !$OMP PARALLEL DO PRIVATE(iblock,i)
    do iblock=1,nblocks_clinic

       associate (diag_cnt => marbl_instances(iblock)%surface_forcing_diags%diag_cnt)

       do i = 1,diag_cnt
          call accumulate_tavg_field(FLUX_DIAGS(:,:,i,iblock), &
               tavg_ids_surface_forcing(i), iblock, i)
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

       end associate

    end do
    !$OMP END PARALLEL DO
    
  end subroutine ecosys_tavg_accumulate_flux

  !***********************************************************************

  subroutine ecosys_tavg_define_from_diag(marbl_diags, tavg_ids)

    implicit none

    type(marbl_diagnostics_type)    , intent(in)    :: marbl_diags
    integer(int_kind), dimension(:) , intent(inout) :: tavg_ids

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(char_len) :: err_msg, gloc, coords
    integer :: n, ndims
    !-----------------------------------------------------------------------

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
         call define_tavg_field(tavg_ids(n),      &
              trim(diags(n)%short_name),          &
              ndims,                              &
              long_name=trim(diags(n)%long_name), &
              units=trim(diags(n)%units),         &
              grid_loc=gloc,                      &
              coordinates=coords)
      end do
    end associate
    
  end subroutine ecosys_tavg_define_from_diag

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


