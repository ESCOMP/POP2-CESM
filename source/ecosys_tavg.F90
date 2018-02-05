
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
  use marbl_interface_public_types , only : marbl_diagnostics_type

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

    use ecosys_diagnostics_operators_mod,   only : ecosys_diagnostics_operators_init
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_surface_stream_cnt
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_interior_stream_cnt

    implicit none

    type(marbl_interface_class)  , intent(in) :: marbl_instance

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_tavg:ecosys_tavg_init'
    character(*), parameter :: marbl_diag_file = 'marbl_diagnostics_operators'
    character(char_len) :: sname, lname
    integer (int_kind) :: n
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

      call ecosys_diagnostics_operators_init(marbl_diag_file, surface_forcing, interior_forcing)

      allocate(tavg_ids_interior_forcing(size(interior_forcing%diags), max_marbl_streams))
      allocate(tavg_ids_surface_forcing(size(surface_forcing%diags), max_marbl_streams))

      call ecosys_tavg_define_from_diag(marbl_diags=interior_forcing, &
           stream_cnt=marbl_diags_interior_stream_cnt, &
           tavg_ids=tavg_ids_interior_forcing)

      call ecosys_tavg_define_from_diag(marbl_diags=surface_forcing,  &
           stream_cnt=marbl_diags_surface_stream_cnt, &
           tavg_ids=tavg_ids_surface_forcing)

    end associate

    call define_tavg_field(tavg_O2_GAS_FLUX_2,'STF_O2_2',2,             &
                           long_name='Dissolved Oxygen Surface Flux',   &
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

  subroutine ecosys_tavg_accumulate_surface(marbl_col_to_pop_i, marbl_col_to_pop_j, &
             STF, marbl_instance, bid)

    use ecosys_diagnostics_operators_mod,   only : marbl_diags_surface_stream_cnt
    use ecosys_tracers_and_saved_state_mod, only : o2_ind

    implicit none

    integer,                     intent(in) :: marbl_col_to_pop_i(:)
    integer,                     intent(in) :: marbl_col_to_pop_j(:)
    real (r8)                  , intent(in) :: STF(:,:,:)
    type(marbl_interface_class), intent(in) :: marbl_instance
    integer,                     intent(in) :: bid

    !-----------------------------------------------------------------------

    ! Accumulate surface_forcing_diags
    call ecosys_tavg_accumulate_from_diag(marbl_col_to_pop_i(:), &
         marbl_col_to_pop_j(:), bid, &
         marbl_diags = marbl_instance%surface_forcing_diags, &
         marbl_stream_cnt = marbl_diags_surface_stream_cnt, &
         tavg_ids = tavg_ids_surface_forcing, &
         num_elements = marbl_instance%surface_forcing_diags%num_elements)

    call accumulate_tavg_field(STF(:,:,o2_ind), tavg_O2_GAS_FLUX_2, bid, 1)

  end subroutine ecosys_tavg_accumulate_surface

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_interior(i, c, marbl_instance, bid)

    use ecosys_diagnostics_operators_mod, only : marbl_diags_interior_stream_cnt

    implicit none

    integer,                     intent(in) :: i, c ! column indices
    type(marbl_interface_class), intent(in) :: marbl_instance
    integer,                     intent(in) :: bid ! block index

    !-----------------------------------------------------------------------

    ! Accumulate diagnostics from marbl_interior_forcing_diags
    call ecosys_tavg_accumulate_from_diag((/i/), (/c/), bid, &
         marbl_diags = marbl_instance%interior_forcing_diags, &
         marbl_stream_cnt = marbl_diags_interior_stream_cnt, &
         tavg_ids = tavg_ids_interior_forcing, &
         num_elements = marbl_instance%interior_forcing_diags%num_elements)

  end subroutine ecosys_tavg_accumulate_interior

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_diags, marbl_stream_cnt, tavg_ids, num_elements)

    ! Accumulate diagnostics

    implicit none

    integer, dimension(:)        , intent(in) :: i, c ! column indices
    integer                      , intent(in) :: bid ! block index
    type(marbl_diagnostics_type) , intent(in) :: marbl_diags
    integer, dimension(:)        , intent(in) :: marbl_stream_cnt
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
         do m=1,marbl_stream_cnt(n)
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

  subroutine ecosys_tavg_define_from_diag(marbl_diags, stream_cnt, tavg_ids)

    use tavg, only : tavg_method_avg

    implicit none

    type(marbl_diagnostics_type),      intent(in)    :: marbl_diags
    integer(int_kind), dimension(:),   intent(in)    :: stream_cnt
    integer(int_kind), dimension(:,:), intent(inout) :: tavg_ids

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(char_len) :: err_msg, gloc, coords, short_name
    integer :: m, n, ndims
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))

      do n=1,size(diags)
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

         do m=1,stream_cnt(n)
           if (m .eq. 1) then
             write(short_name, "(A)") trim(diags(n)%short_name)
           else
             write(short_name, "(A,'_',I0)") trim(diags(n)%short_name), m
           end if
           call define_tavg_field(tavg_ids(n,m),    &
                short_name,                         &
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

