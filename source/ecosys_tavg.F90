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

  use ecosys_diagnostics_operators_mod, only : max_marbl_streams

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_tavg_init
  public :: ecosys_tavg_accumulate_interior
  public :: ecosys_tavg_accumulate_surface
  public :: ecosys_tavg_accumulate_scalar_rmeans

  !-----------------------------------------------------------------------
  !  define tavg id for interior diagnostics, diagnostics related to
  !  restoring, surface forcing diagnostics, and duplicate surface forcing
  !  variables
  !-----------------------------------------------------------------------

  integer (int_kind), allocatable :: tavg_ids_interior_forcing(:,:)
  integer (int_kind), allocatable :: tavg_ids_surface_forcing(:,:)

  integer (int_kind) :: tavg_ECOSYS_IFRAC_2 ! ice fraction duplicate
  integer (int_kind) :: tavg_ECOSYS_XKW_2   ! xkw duplicate
  integer (int_kind) :: tavg_O2_GAS_FLUX_2  ! O2 flux duplicate
  integer (int_kind) :: tavg_DpCO2_2        ! delta pco2 duplicate
  integer (int_kind) :: tavg_DIC_GAS_FLUX_2 ! dic flux duplicate

  integer (int_kind), allocatable :: tavg_ids_scalar_rmean_interior(:)
  integer (int_kind), allocatable :: tavg_ids_scalar_rmean_surface(:)

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_tavg_init(marbl_instance)

    ! !DESCRIPTION:
    !  call define_tavg_field for all tavg fields

    use ecosys_tracers_and_saved_state_mod, only : marbl_tracer_cnt
    use ecosys_diagnostics_operators_mod,   only : ecosys_diagnostics_operators_init

    implicit none

    type(marbl_interface_class)  , intent(in) :: marbl_instance

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'
    character(*), parameter :: marbl_diag_file = 'marbl_diagnostics_operators'
    character(char_len) :: sname, lname
    integer (int_kind) :: n, num_interior_diags, num_surface_diags
    integer (int_kind) :: rmean_var_cnt
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  1. allocate memory for module variables
    !  2. Read marbl_diag_file
    !  3. define tavg fields for MARBL diagnostics
    !-----------------------------------------------------------------------

    associate(&
         interior_forcing => marbl_instance%interior_forcing_diags, &
         surface_forcing => marbl_instance%surface_forcing_diags    &
         )

      num_interior_diags = size(interior_forcing%diags)
      num_surface_diags  = size(surface_forcing%diags)
      allocate(tavg_ids_interior_forcing(num_interior_diags, max_marbl_streams))
      allocate(tavg_ids_surface_forcing(num_surface_diags, max_marbl_streams))

      call ecosys_diagnostics_operators_init(marbl_diag_file, num_surface_diags + num_interior_diags)

      call ecosys_tavg_define_from_diag(marbl_diags=interior_forcing, &
           tavg_ids=tavg_ids_interior_forcing)

      call ecosys_tavg_define_from_diag(marbl_diags=surface_forcing,  &
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

    rmean_var_cnt = size(marbl_instance%glo_scalar_rmean_interior)
    allocate(tavg_ids_scalar_rmean_interior(rmean_var_cnt))
    do n = 1, rmean_var_cnt
      call define_tavg_field(tavg_ids_scalar_rmean_interior(n), &
                             marbl_instance%glo_scalar_rmean_interior(n)%sname, 0)
    end do

    rmean_var_cnt = size(marbl_instance%glo_scalar_rmean_surface)
    allocate(tavg_ids_scalar_rmean_surface(rmean_var_cnt))
    do n = 1, rmean_var_cnt
      call define_tavg_field(tavg_ids_scalar_rmean_surface(n), &
                             marbl_instance%glo_scalar_rmean_surface(n)%sname, 0)
    end do

  end subroutine ecosys_tavg_init

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_interior(i, c, bid, &
       marbl_interior_forcing_diags)

    implicit none

    integer , dimension(:), intent(in) :: i, c ! column indices
    integer ,               intent(in) :: bid ! block index

    type(marbl_diagnostics_type), intent(in) :: marbl_interior_forcing_diags
    !-----------------------------------------------------------------------

    call ecosys_tavg_accumulate_from_diag(i, c, bid, &
         marbl_diags = marbl_interior_forcing_diags,  &
         tavg_ids = tavg_ids_interior_forcing, &
         num_elements =marbl_interior_forcing_diags%num_elements)

  end subroutine ecosys_tavg_accumulate_interior

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_diags, tavg_ids, num_elements)

    ! Accumulate diagnostics

    implicit none

    integer, dimension(:)        , intent(in) :: i, c ! column indices
    integer                      , intent(in) :: bid ! block index
    type(marbl_diagnostics_type) , intent(in) :: marbl_diags
    integer, dimension(:,:)      , intent(in) :: tavg_ids
    integer                      , intent(in) :: num_elements

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: m, n, ne
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))

      do n=1,size(diags)
         do ne = 1,num_elements
           ! FIXME: once we are reading variable names from marbl_diags_sname,
           !        loop from 1 to marbl_diags_stream_cnt(i), where
           !        marbl_diags_sname(i) == diags(n)%short_name
           do m=1, 1
             if (allocated(diags(n)%field_2d)) then
               call accumulate_tavg_field(diags(n)%field_2d(ne)  , tavg_ids(n,m), bid, i(ne), c(ne))
             else
               call accumulate_tavg_field(diags(n)%field_3d(:,ne), tavg_ids(n,m), bid, i(ne), c(ne))
             end if
           end do
         end do
      end do

    end associate

  end subroutine ecosys_tavg_accumulate_from_diag

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_surface(surface_forcing_diags, STF, marbl_instances)

    ! Accumulate diagnostics for surface fluxes

    use ecosys_tracers_and_saved_state_mod, only : marbl_tracer_cnt, o2_ind
    use marbl_diagnostics_mod,              only : ind => marbl_surface_forcing_diag_ind

    implicit none

    real (r8)                  , intent(in) :: surface_forcing_diags (:, :, :, :)
    real (r8)                  , intent(in) :: STF(:,:,:,:)
    type(marbl_interface_class), intent(in) :: marbl_instances(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: i, m, n, iblock
    integer :: nblocks_clinic
    !-----------------------------------------------------------------------

    nblocks_clinic = size(surface_forcing_diags,4)

    !$OMP PARALLEL DO PRIVATE(iblock,i)
    do iblock=1,nblocks_clinic

       associate (diags => marbl_instances(iblock)%surface_forcing_diags%diags)

       do i = 1,size(diags)
         ! FIXME: once we are reading variable names from marbl_diags_sname,
         !        loop from 1 to marbl_diags_stream_cnt(i), where
         !        marbl_diags_sname(i) == diags(n)%short_name
         do m=1, 1
            call accumulate_tavg_field(surface_forcing_diags(:,:,i,iblock),     &
                 tavg_ids_surface_forcing(i,m), iblock, 1)
         end do
       end do

       call accumulate_tavg_field(surface_forcing_diags(:,:,ind%ECOSYS_IFRAC,iblock), &
            tavg_ECOSYS_IFRAC_2, iblock, 1)

       call accumulate_tavg_field(surface_forcing_diags(:,:,ind%ECOSYS_XKW,iblock),   &
            tavg_ECOSYS_XKW_2, iblock, 1)

       call accumulate_tavg_field(STF(:,:,o2_ind,iblock), tavg_O2_GAS_FLUX_2, iblock, 1)

       call accumulate_tavg_field(surface_forcing_diags(:,:,ind%DpCO2,iblock),        &
            tavg_DpCO2_2, iblock, 1)

       call accumulate_tavg_field(surface_forcing_diags(:,:,ind%DIC_GAS_FLUX,iblock), &
            tavg_DIC_GAS_FLUX_2, iblock, 1)

       end associate

    end do
    !$OMP END PARALLEL DO

  end subroutine ecosys_tavg_accumulate_surface

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_scalar_rmeans(marbl_instance, field_source)

    ! Accumulate diagnostics for scalar running means

    implicit none

    type(marbl_interface_class), intent(in) :: marbl_instance
    character (*),               intent(in) :: field_source   ! 'interior' or 'surface'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: n
    !-----------------------------------------------------------------------

    if (trim(field_source) == 'interior') then
      do n = 1, size(marbl_instance%glo_scalar_rmean_interior)
        call accumulate_tavg_field(marbl_instance%glo_scalar_rmean_interior(n)%rmean, &
                                   tavg_ids_scalar_rmean_interior(n))
      end do
    else
      do n = 1, size(marbl_instance%glo_scalar_rmean_surface)
        call accumulate_tavg_field(marbl_instance%glo_scalar_rmean_surface(n)%rmean, &
                                   tavg_ids_scalar_rmean_surface(n))
      end do
    end if

  end subroutine ecosys_tavg_accumulate_scalar_rmeans

  !***********************************************************************

  subroutine ecosys_tavg_define_from_diag(marbl_diags, tavg_ids)

    use tavg, only : tavg_method_avg

    implicit none

    type(marbl_diagnostics_type),      intent(in)    :: marbl_diags
    integer(int_kind), dimension(:,:), intent(inout) :: tavg_ids

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(char_len) :: err_msg, gloc, coords
    integer :: m, n, ndims
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))

      do n=1,size(diags)
         ! FIXME: Continue to next n if diags(n)%short_name not in marbl_diags_sname
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

         ! FIXME: once we are reading variable names from marbl_diags_sname,
         !        loop from 1 to marbl_diags_stream_cnt(i), where
         !        marbl_diags_sname(i) == diags(n)%short_name
         do m=1, 1
           call define_tavg_field(tavg_ids(n,m),    &
                trim(diags(n)%short_name),          &
                ndims,                              &
                tavg_method = tavg_method_avg,      &
                long_name=trim(diags(n)%long_name), &
                units=trim(diags(n)%units),         &
                grid_loc=gloc,                      &
                coordinates=coords,                 &
                transpose_field=(ndims .eq. 3))
         end do
      end do
    end associate

  end subroutine ecosys_tavg_define_from_diag

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


