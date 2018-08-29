module ecosys_running_mean_saved_state_mod

  !-----------------------------------------------------------------------
  ! running-mean related saved state variables
  !-----------------------------------------------------------------------

  use kinds_mod, only : int_kind, r8

  implicit none
  private

  !-----------------------------------------------------------------------
  ! module public variables
  !-----------------------------------------------------------------------

  integer (int_kind), allocatable, target, dimension(:) :: glo_avg_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:) :: glo_avg_rmean_ind_surface
  integer (int_kind), allocatable, target, dimension(:) :: glo_scalar_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:) :: glo_scalar_rmean_ind_surface

  public :: ecosys_running_mean_state_init
  public :: ecosys_running_mean_saved_state_get_var_vals
  public :: ecosys_running_mean_saved_state_update

  !***********************************************************************

contains

  !***********************************************************************

  subroutine ecosys_running_mean_state_init(marbl_instance, ecosys_restart_filename)

    use marbl_interface,  only : marbl_interface_class
    use running_mean_mod, only : running_mean_get_var

    type(marbl_interface_class), intent(in) :: marbl_instance
    character(len=*),            intent(in) :: ecosys_restart_filename

    integer :: n, iblock
    real(r8) :: rmean_val

    ! FIXME : move the setup of running means of global averages into MARBL

    call init_rmean_var(marbl_instance%glo_avg_rmean_interior_tendency, ecosys_restart_filename, &
         glo_avg_rmean_ind_interior)

    call init_rmean_var(marbl_instance%glo_avg_rmean_surface_flux, ecosys_restart_filename, &
         glo_avg_rmean_ind_surface)

    call init_rmean_var(marbl_instance%glo_scalar_rmean_interior_tendency, ecosys_restart_filename, &
         glo_scalar_rmean_ind_interior)

    call init_rmean_var(marbl_instance%glo_scalar_rmean_surface_flux, ecosys_restart_filename, &
         glo_scalar_rmean_ind_surface)

  end subroutine ecosys_running_mean_state_init

  !***********************************************************************

  subroutine ecosys_running_mean_saved_state_get_var_vals(field_source, lscalar, array_out)

    use running_mean_mod , only : running_mean_get_var

    character(len=*), intent(in)  :: field_source
    logical,          intent(in)  :: lscalar
    real(r8),         intent(out) :: array_out(:)

    integer, pointer :: glo_rmean_ind(:)
    integer :: n

    call get_glo_rmean_ind_pointer(field_source, lscalar, glo_rmean_ind)

    do n=1,size(glo_rmean_ind)
      call running_mean_get_var(glo_rmean_ind(n), vals_0d = array_out(n))
    end do

  end subroutine ecosys_running_mean_saved_state_get_var_vals

  !***********************************************************************

  subroutine ecosys_running_mean_saved_state_update(field_source, lscalar, array_out)

    use running_mean_mod , only : running_mean_update_var

    character(len=*), intent(in)  :: field_source
    logical,          intent(in)  :: lscalar
    real(r8),         intent(out) :: array_out(:)

    integer, pointer :: glo_rmean_ind(:)
    integer :: n

    call get_glo_rmean_ind_pointer(field_source, lscalar, glo_rmean_ind)

    do n=1,size(glo_rmean_ind)
      call running_mean_update_var(glo_rmean_ind(n), vals_0d=array_out(n))
    end do

  end subroutine ecosys_running_mean_saved_state_update

  !***********************************************************************

  subroutine init_rmean_var(marbl_running_mean_var, ecosys_restart_filename, rmean_ind)

    use marbl_interface_public_types, only : marbl_running_mean_0d_type
    use running_mean_mod, only : running_mean_define_var
    use running_mean_mod, only : running_mean_init_var

    type(marbl_running_mean_0d_type), intent(in)  :: marbl_running_mean_var(:)
    character(len=*),                 intent(in)  :: ecosys_restart_filename
    integer (int_kind), allocatable,  intent(out) :: rmean_ind(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: rmean_var_cnt
    integer (int_kind) :: n

    !-----------------------------------------------------------------------

    rmean_var_cnt = size(marbl_running_mean_var, 1)
    allocate(rmean_ind(rmean_var_cnt))

    do n = 1, rmean_var_cnt
       call running_mean_define_var(name=trim(marbl_running_mean_var(n)%sname), rank=0, &
          timescale=marbl_running_mean_var(n)%timescale, index=rmean_ind(n))

       if (marbl_running_mean_var(n)%linit_by_val) then
          call running_mean_init_var(rmean_ind(n), vals_0d=marbl_running_mean_var(n)%init_val)
       else
          call running_mean_init_var(rmean_ind(n), filename=ecosys_restart_filename)
       end if
    end do

  end subroutine init_rmean_var

  !***********************************************************************

  subroutine get_glo_rmean_ind_pointer(field_source, lscalar, ptr_out)

    character(len=*), intent(in)  :: field_source
    logical,          intent(in)  :: lscalar
    integer, pointer, intent(out) :: ptr_out(:)

    if (trim(field_source) .eq. 'interior_tendency') then
      if (lscalar) then
        ptr_out => glo_scalar_rmean_ind_interior
      else
        ptr_out => glo_avg_rmean_ind_interior
      end if
    else
      if (lscalar) then
        ptr_out => glo_scalar_rmean_ind_surface
      else
        ptr_out => glo_avg_rmean_ind_surface
      end if
    end if

  end subroutine get_glo_rmean_ind_pointer

  !***********************************************************************

end module ecosys_running_mean_saved_state_mod
