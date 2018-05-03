module ecosys_running_mean_saved_state_mod

  use kinds_mod, only : int_kind, r8

  implicit none
  private

  ! FIXME : move the following variables into MARBL

  integer (int_kind), allocatable, target, dimension(:), public :: glo_avg_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:), public :: glo_avg_rmean_ind_surface
  integer (int_kind), allocatable, target, dimension(:), public :: glo_scalar_rmean_ind_interior
  integer (int_kind), allocatable, target, dimension(:), public :: glo_scalar_rmean_ind_surface

  public :: ecosys_running_mean_state_init
  public :: ecosys_running_mean_saved_state_get_var_val
  public :: ecosys_running_mean_saved_state_update

contains

  subroutine ecosys_running_mean_state_init(ecosys_restart_filename, marbl_instances)

    use marbl_interface,  only : marbl_interface_class
    use running_mean_mod, only : running_mean_get_var

    character(len=*),            intent(in)  :: ecosys_restart_filename
    type(marbl_interface_class), intent(inout) :: marbl_instances(:)

    integer :: n, iblock
    real(r8) :: rmean_val

    ! FIXME : move the setup of running means of global averages into MARBL

    call init_rmean_var(marbl_instances(1)%glo_avg_rmean_interior, ecosys_restart_filename, &
         glo_avg_rmean_ind_interior)

    call init_rmean_var(marbl_instances(1)%glo_avg_rmean_surface, ecosys_restart_filename, &
         glo_avg_rmean_ind_surface)

    call init_rmean_var(marbl_instances(1)%glo_scalar_rmean_interior, ecosys_restart_filename, &
         glo_scalar_rmean_ind_interior)

    call init_rmean_var(marbl_instances(1)%glo_scalar_rmean_surface, ecosys_restart_filename, &
         glo_scalar_rmean_ind_surface)

    ! copy values from POP's running mean to MARBL interface

    do n = 1, size(glo_avg_rmean_ind_interior(:))
       call running_mean_get_var(glo_avg_rmean_ind_interior(n), vals_0d=rmean_val)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_avg_rmean_interior(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_avg_rmean_ind_surface(:))
       call running_mean_get_var(glo_avg_rmean_ind_surface(n), vals_0d=rmean_val)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_avg_rmean_surface(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_scalar_rmean_ind_interior(:))
       call running_mean_get_var(glo_scalar_rmean_ind_interior(n), vals_0d=rmean_val)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_scalar_rmean_interior(n)%rmean = rmean_val
       end do
    end do

    do n = 1, size(glo_scalar_rmean_ind_surface(:))
       call running_mean_get_var(glo_scalar_rmean_ind_surface(n), vals_0d=rmean_val)
       do iblock = 1, size(marbl_instances)
          marbl_instances(iblock)%glo_scalar_rmean_surface(n)%rmean = rmean_val
       end do
    end do

  end subroutine ecosys_running_mean_state_init

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

  subroutine ecosys_running_mean_saved_state_get_var_val(index, block, k, vals_0d, &
       vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
       vals_3d_1klev_1block, vals_3d_1klev_blocks)

    use running_mean_mod , only : running_mean_get_var

    integer (int_kind), intent(in)                       :: index
    integer (int_kind), intent(in), optional             :: block
    integer (int_kind), intent(in), optional             :: k

 ! !OUTPUT PARAMETERS:

    real (r8), intent(out), optional                     :: vals_0d
    real (r8), intent(out), optional                     :: vals_1d_1klev
    real (r8), intent(out), dimension(:), optional       :: vals_1d_klevs
    real (r8), intent(out), dimension(:,:), optional     :: vals_2d_1block
    real (r8), intent(out), dimension(:,:,:), optional   :: vals_2d_blocks
    real (r8), intent(out), dimension(:,:), optional     :: vals_3d_1klev_1block
    real (r8), intent(out), dimension(:,:,:), optional   :: vals_3d_1klev_blocks

    call running_mean_get_var(index, block, k, vals_0d, vals_1d_1klev, &
         vals_1d_klevs, vals_2d_1block, vals_2d_blocks, vals_3d_1klev_1block, &
         vals_3d_1klev_blocks)

  end subroutine ecosys_running_mean_saved_state_get_var_val

  !***********************************************************************

  subroutine ecosys_running_mean_saved_state_update(index, block, k, vals_0d, &
       vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
       vals_3d_1klev_1block, vals_3d_1klev_blocks)

    use running_mean_mod , only : running_mean_update_var

    integer (int_kind), intent(in)                      :: index
    integer (int_kind), intent(in), optional            :: block
    integer (int_kind), intent(in), optional            :: k
    real (r8), intent(in), optional                     :: vals_0d
    real (r8), intent(in), optional                     :: vals_1d_1klev
    real (r8), intent(in), dimension(:), optional       :: vals_1d_klevs
    real (r8), intent(in), dimension(:,:), optional     :: vals_2d_1block
    real (r8), intent(in), dimension(:,:,:), optional   :: vals_2d_blocks
    real (r8), intent(in), dimension(:,:), optional     :: vals_3d_1klev_1block
    real (r8), intent(in), dimension(:,:,:), optional   :: vals_3d_1klev_blocks

    call running_mean_update_var(index, block, k, vals_0d, &
         vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
         vals_3d_1klev_1block, vals_3d_1klev_blocks)

  end subroutine ecosys_running_mean_saved_state_update

  !***********************************************************************

end module ecosys_running_mean_saved_state_mod