
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
  use marbl_logging         , only : marbl_log_type
  use marbl_interface_public_types , only : marbl_diagnostics_type

  use ecosys_diagnostics_operators_mod, only : max_marbl_diags_stream_cnt

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: ecosys_tavg_init
  public :: ecosys_tavg_accumulate_interior
  public :: ecosys_tavg_accumulate_surface
  public :: ecosys_tavg_accumulate_scalar_rmeans
  public :: ecosys_tavg_set_compute_now

  !-----------------------------------------------------------------------
  !  define tavg id for interior tendency diagnostics, diagnostics related
  !  to restoring, surface flux diagnostics, and duplicate surface flux
  !  variables
  !-----------------------------------------------------------------------

  integer (int_kind), allocatable :: tavg_ids_interior_tendency(:,:)
  integer (int_kind), allocatable :: tavg_ids_surface_flux(:,:)
  integer (int_kind) :: tavg_O2_GAS_FLUX_2  ! O2 flux duplicate

  integer (int_kind), allocatable :: tavg_ids_scalar_rmean_interior(:)
  integer (int_kind), allocatable :: tavg_ids_scalar_rmean_surface(:)

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_tavg_init(marbl_instance, status_log)

    ! !DESCRIPTION:
    !  call define_tavg_field for all tavg fields

    use ecosys_diagnostics_operators_mod,   only : ecosys_diagnostics_operators_init
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_stream_cnt_surface
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_stream_cnt_interior

    type(marbl_interface_class), intent(in)    :: marbl_instance
    type(marbl_log_type),        intent(inout) :: status_log

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
         interior_tendency => marbl_instance%interior_tendency_diags, &
         surface_flux      => marbl_instance%surface_flux_diags       &
         )

      call ecosys_diagnostics_operators_init(marbl_diag_file, surface_flux, interior_tendency)

      allocate(tavg_ids_interior_tendency(size(interior_tendency%diags), max_marbl_diags_stream_cnt))
      allocate(tavg_ids_surface_flux(size(surface_flux%diags), max_marbl_diags_stream_cnt))

      call ecosys_tavg_define_from_diag(marbl_diags=interior_tendency, &
           stream_cnt=marbl_diags_stream_cnt_interior, &
           tavg_ids=tavg_ids_interior_tendency, &
           status_log=status_log)

      call ecosys_tavg_define_from_diag(marbl_diags=surface_flux,  &
           stream_cnt=marbl_diags_stream_cnt_surface, &
           tavg_ids=tavg_ids_surface_flux, &
           status_log=status_log)

    end associate

    call define_tavg_field(tavg_O2_GAS_FLUX_2,'STF_O2_2',2,             &
                           long_name='Dissolved Oxygen Surface Flux',   &
                           units='mmol/m^3 cm/s', grid_loc='2110',      &
                           coordinates='TLONG TLAT time')

    rmean_var_cnt = size(marbl_instance%glo_scalar_rmean_interior_tendency)
    allocate(tavg_ids_scalar_rmean_interior(rmean_var_cnt))
    do n = 1, rmean_var_cnt
      call define_tavg_field(tavg_ids_scalar_rmean_interior(n), &
                             marbl_instance%glo_scalar_rmean_interior_tendency(n)%sname, 0)
    end do

    rmean_var_cnt = size(marbl_instance%glo_scalar_rmean_surface_flux)
    allocate(tavg_ids_scalar_rmean_surface(rmean_var_cnt))
    do n = 1, rmean_var_cnt
      call define_tavg_field(tavg_ids_scalar_rmean_surface(n), &
                             marbl_instance%glo_scalar_rmean_surface_flux(n)%sname, 0)
    end do

  end subroutine ecosys_tavg_init

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_surface(marbl_col_to_pop_i, marbl_col_to_pop_j, &
             STF, marbl_instance, bid)

    use ecosys_diagnostics_operators_mod,   only : marbl_diags_stream_cnt_surface
    use ecosys_tracers_and_saved_state_mod, only : o2_ind

    integer,                     intent(in) :: marbl_col_to_pop_i(:)
    integer,                     intent(in) :: marbl_col_to_pop_j(:)
    real (r8)                  , intent(in) :: STF(:,:,:)
    type(marbl_interface_class), intent(in) :: marbl_instance
    integer,                     intent(in) :: bid

    !-----------------------------------------------------------------------

    ! Accumulate surface_flux_diags
    call ecosys_tavg_accumulate_from_diag(marbl_col_to_pop_i(:), &
         marbl_col_to_pop_j(:), bid, &
         marbl_diags = marbl_instance%surface_flux_diags, &
         marbl_diags_stream_cnt = marbl_diags_stream_cnt_surface, &
         tavg_ids = tavg_ids_surface_flux, &
         num_elements = marbl_instance%surface_flux_diags%num_elements)

    call accumulate_tavg_field(STF(:,:,o2_ind), tavg_O2_GAS_FLUX_2, bid, 1)

  end subroutine ecosys_tavg_accumulate_surface

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_interior(i, c, marbl_instance, bid)

    use ecosys_diagnostics_operators_mod, only : marbl_diags_stream_cnt_interior

    integer,                     intent(in) :: i, c ! column indices
    type(marbl_interface_class), intent(in) :: marbl_instance
    integer,                     intent(in) :: bid ! block index

    !-----------------------------------------------------------------------

    ! Accumulate diagnostics from marbl_interior_tendency_diags
    call ecosys_tavg_accumulate_from_diag((/i/), (/c/), bid, &
         marbl_diags = marbl_instance%interior_tendency_diags, &
         marbl_diags_stream_cnt = marbl_diags_stream_cnt_interior, &
         tavg_ids = tavg_ids_interior_tendency, &
         num_elements = marbl_instance%interior_tendency_diags%num_elements)

  end subroutine ecosys_tavg_accumulate_interior

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_from_diag(i, c, bid, marbl_diags, marbl_diags_stream_cnt, tavg_ids, num_elements)

    ! Accumulate diagnostics

    integer, dimension(:)        , intent(in) :: i, c ! column indices
    integer                      , intent(in) :: bid ! block index
    type(marbl_diagnostics_type) , intent(in) :: marbl_diags
    integer, dimension(:)        , intent(in) :: marbl_diags_stream_cnt
    integer, dimension(:,:)      , intent(in) :: tavg_ids
    integer                      , intent(in) :: num_elements

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: m, n
    !-----------------------------------------------------------------------

    associate(diags => marbl_diags%diags(:))

    do n=1,size(diags)
      if (allocated(diags(n)%field_2d)) then
         do m=1,marbl_diags_stream_cnt(n)
           call accumulate_tavg_field(diags(n)%field_2d(:), tavg_ids(n,m), bid, i(:), c(:))
         end do
       else
         do m=1,marbl_diags_stream_cnt(n)
           call accumulate_tavg_field(diags(n)%field_3d(:,:), tavg_ids(n,m), bid, i(:), c(:))
         end do
       end if
    end do

    end associate

  end subroutine ecosys_tavg_accumulate_from_diag

  !***********************************************************************

  subroutine ecosys_tavg_accumulate_scalar_rmeans(marbl_instance, field_source)

    ! Accumulate diagnostics for scalar running means

    type(marbl_interface_class), intent(in) :: marbl_instance
    character (*),               intent(in) :: field_source   ! 'interior_tendency' or 'surface_flux'

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: n
    !-----------------------------------------------------------------------

    if (trim(field_source) == 'interior_tendency') then
      do n = 1, size(marbl_instance%glo_scalar_rmean_interior_tendency)
        call accumulate_tavg_field(marbl_instance%glo_scalar_rmean_interior_tendency(n)%rmean, &
                                   tavg_ids_scalar_rmean_interior(n))
      end do
    else
      do n = 1, size(marbl_instance%glo_scalar_rmean_surface_flux)
        call accumulate_tavg_field(marbl_instance%glo_scalar_rmean_surface_flux(n)%rmean, &
                                   tavg_ids_scalar_rmean_surface(n))
      end do
    end if

  end subroutine ecosys_tavg_accumulate_scalar_rmeans

  !***********************************************************************

  subroutine ecosys_tavg_define_from_diag(marbl_diags, stream_cnt, tavg_ids, status_log)

    use tavg, only : tavg_method_avg
    use pop_constants, only : cmperm
    use domain_size, only : km
    use grid, only : zw

    type(marbl_diagnostics_type),      intent(in)    :: marbl_diags
    integer(int_kind), dimension(:),   intent(in)    :: stream_cnt
    integer(int_kind), dimension(:,:), intent(inout) :: tavg_ids
    type(marbl_log_type),              intent(inout) :: status_log

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(char_len) :: err_msg, gloc, coords, short_name, vert_grid
    integer :: m, n, ndims

    real (r8) :: ref_depth_cm
    integer :: ref_k
    !-----------------------------------------------------------------------


    do n=1,size(marbl_diags%diags)
      vert_grid = marbl_diags%get_str(n, 'vertical_grid', status_log)
      if (trim(vert_grid).eq.'none') then
        ndims = 2
        gloc = '2110'
        coords = 'TLONG TLAT time'
        ! find layer containing ref_depth, i.e., zw(k-1) .le. ref_depth .lt. zw(k)
        ref_depth_cm = cmperm * marbl_diags%diags(n)%ref_depth
        if (ref_depth_cm .lt. zw(km)) then
          do ref_k = 1, km
            if (ref_depth_cm .lt. zw(ref_k)) exit
          end do
        else
          ref_k = km
        end if
      else
        ndims = 3
        if (trim(vert_grid) .eq. 'layer_avg') then
            if (marbl_diags%diags(n)%ltruncated_vertical_extent) then
              gloc = '3114'
              coords = 'TLONG TLAT z_t_150m time'
            else
              gloc = '3111'
              coords = 'TLONG TLAT z_t time'
            end if
        elseif (trim(vert_grid).eq.'layer_iface') then
            gloc = '3112'
            coords = 'TLONG TLAT z_w_top time'
        else
            write(err_msg,*) "'", trim(vert_grid), &
                "' is not a valid vertical grid"
            call shr_sys_abort(err_msg)
        end if
      end if

      do m=1,stream_cnt(n)
        if (m .eq. 1) then
          write(short_name, "(A)") trim(marbl_diags%get_str(n, 'short_name', status_log))
        else
          write(short_name, "(A,'_',I0)") trim(marbl_diags%get_str(n, 'short_name', status_log)), m
        end if
        if (ndims .eq. 2) then
          call define_tavg_field(tavg_ids(n,m),                          &
              short_name,                                                &
              ndims,                                                     &
              tavg_method = tavg_method_avg,                             &
              long_name=marbl_diags%get_str(n, 'long_name', status_log), &
              units=marbl_diags%get_str(n, 'units', status_log),         &
              grid_loc=gloc,                                             &
              mask_k=ref_k,                                              &
              coordinates=coords,                                        &
              transpose_field=(ndims .eq. 3))
        else
          call define_tavg_field(tavg_ids(n,m),                          &
              short_name,                                                &
              ndims,                                                     &
              tavg_method = tavg_method_avg,                             &
              long_name=marbl_diags%get_str(n, 'long_name', status_log), &
              units=marbl_diags%get_str(n, 'units', status_log),         &
              grid_loc=gloc,                                             &
              coordinates=coords,                                        &
              transpose_field=(ndims .eq. 3))
        end if
      end do
    end do

  end subroutine ecosys_tavg_define_from_diag

  !***********************************************************************

  subroutine ecosys_tavg_set_compute_now(marbl_diags, field_source, status_log)

    use tavg, only : set_in_tavg_contents
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_stream_cnt_surface
    use ecosys_diagnostics_operators_mod,   only : marbl_diags_stream_cnt_interior

    type(marbl_diagnostics_type), intent(inout) :: marbl_diags
    character(len=*),             intent(in)    :: field_source
    type(marbl_log_type),         intent(inout) :: status_log

    character(len=char_len) :: err_msg
    logical :: in_tavg_contents
    integer :: n, m

    select case (trim(field_source))
      case ('surface_flux')
        do n=1,size(marbl_diags%diags)
          in_tavg_contents = .false.
          do m=1,marbl_diags_stream_cnt_surface(n)
            in_tavg_contents = in_tavg_contents .or. set_in_tavg_contents(tavg_ids_surface_flux(n,m))
          end do
          call marbl_diags%set(n, 'compute_now', in_tavg_contents, status_log)
        end do
      case ('interior_tendency')
        do n=1,size(marbl_diags%diags)
          in_tavg_contents = .false.
          do m=1,marbl_diags_stream_cnt_interior(n)
            in_tavg_contents = in_tavg_contents .or. set_in_tavg_contents(tavg_ids_interior_tendency(n,m))
          end do
          call marbl_diags%set(n, 'compute_now', in_tavg_contents, status_log)
        end do
      case DEFAULT
        write(err_msg, "(3A)") "'", trim(field_source), "' is not a valid field source"
        call shr_sys_abort(err_msg)
    end select

    do n=1,size(marbl_diags%diags)
      write(err_msg, *) "compute_now for '", trim(marbl_diags%get_str(n, 'short_name', status_log)), "': ", marbl_diags%diags(n)%compute_now
      call status_log%log_noerror(err_msg, 'tmp_subname')
    end do

  end subroutine ecosys_tavg_set_compute_now

end module ecosys_tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
