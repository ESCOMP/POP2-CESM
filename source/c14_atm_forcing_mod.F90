module c14_atm_forcing_mod

  ! !DESCRIPTION
  !  This module provides an interface for providing atmospheric boundary
  !  conditions for c14.

  use kinds_mod,              only: char_len, r8, int_kind
  use blocks,                 only: nx_block, ny_block
  use io_tools,               only: document
  use exit_mod,               only: sigAbort, exit_POP
  use forcing_timeseries_mod, only: forcing_timeseries_dataset

  implicit none
  private

  public :: c14_atm_forcing_init
  public :: c14_atm_forcing_update_data
  public :: c14_atm_forcing_get_data

  !-----------------------------------------------------------------------
  !  module private data
  !-----------------------------------------------------------------------

  character(char_len) :: &
    c14_atm_forcing_opt             ! option for D14C varying or constant forcing

  real (r8) :: &
    c14_atm_forcing_const           ! atmospheric 14CO2 constant [permil]

  real (r8), dimension(3) :: &
    c14_atm_forcing_lat_band_vals   ! atmospheric 14CO2 constant [permil]

  character (char_len) :: &
    c14_atm_forcing_filename        ! filename for varying atm D14C

  integer (int_kind) :: &
    c14_atm_forcing_model_year,   & ! arbitrary model year
    c14_atm_forcing_data_year       ! year in data that corresponds to c14_atm_forcing_model_year

  type (forcing_timeseries_dataset) :: &
    c14_atm_forcing_dataset         ! data structure for atm 14C timeseries

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine c14_atm_forcing_init(caller, c14_atm_forcing_opt_in, &
      c14_atm_forcing_const_in, c14_atm_forcing_lat_band_vals_in, &
      c14_atm_forcing_filename_in, c14_atm_forcing_model_year_in, c14_atm_forcing_data_year_in)

    use kinds_mod,              only: log_kind
    use forcing_timeseries_mod, only: forcing_timeseries_init_dataset
    use forcing_timeseries_mod, only: forcing_timeseries_dataset_var_size

    character(*), intent(in) :: &
      caller,                          & ! name of module calling c14_atm_forcing_init
      c14_atm_forcing_opt_in             ! option for D14C varying or constant forcing

    real (r8), intent(in) :: &
      c14_atm_forcing_const_in           ! atmospheric 14CO2 constant [permil]

    real (r8), dimension(3), intent(in) :: &
      c14_atm_forcing_lat_band_vals_in   ! atmospheric 14CO2 constant [permil]

    character (*), intent(in) :: &
      c14_atm_forcing_filename_in        ! filename for varying atm D14C

    integer (int_kind), intent(in) :: &
      c14_atm_forcing_model_year_in,   & ! arbitrary model year
      c14_atm_forcing_data_year_in       ! year in data that corresponds to c14_atm_forcing_model_year_in

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(len=*), parameter  :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_init'

    logical(log_kind)            :: first_call = .true.
    logical(log_kind)            :: val_mismatch
    integer(int_kind)            :: lat_ind

    !-----------------------------------------------------------------------

    call document(subname, 'first_call', first_call)
    call document(subname, 'caller', caller)

    !-----------------------------------------------------------------------
    !  if this is not the first call, verify that the
    !  arguments match previously passed in values
    !-----------------------------------------------------------------------

    if (.not. first_call) then
      val_mismatch = .false.
      if (trim(c14_atm_forcing_opt) /= trim(c14_atm_forcing_opt_in)) then
        call document(subname, 'c14_atm_forcing_opt', c14_atm_forcing_opt)
        call document(subname, 'c14_atm_forcing_opt_in', c14_atm_forcing_opt_in)
        val_mismatch = .true.
      else
        select case (trim(c14_atm_forcing_opt))
        case ('const')
          if (c14_atm_forcing_const /= c14_atm_forcing_const_in) then
            call document(subname, 'c14_atm_forcing_const', c14_atm_forcing_const)
            call document(subname, 'c14_atm_forcing_const_in', c14_atm_forcing_const_in)
            val_mismatch = .true.
          end if
        case ('lat_bands')
          do lat_ind = 1, 3
            if (c14_atm_forcing_lat_band_vals(lat_ind) /= c14_atm_forcing_lat_band_vals_in(lat_ind)) then
              call document(subname, 'lat_ind', lat_ind)
              call document(subname, 'c14_atm_forcing_lat_band_vals(lat_ind)', c14_atm_forcing_lat_band_vals(lat_ind))
              call document(subname, 'c14_atm_forcing_lat_band_vals_in(lat_ind)', c14_atm_forcing_lat_band_vals_in(lat_ind))
              val_mismatch = .true.
            end if
          end do
        case ('file')
          if (trim(c14_atm_forcing_filename) /= trim(c14_atm_forcing_filename_in)) then
            call document(subname, 'c14_atm_forcing_filename', c14_atm_forcing_filename)
            call document(subname, 'c14_atm_forcing_filename_in', c14_atm_forcing_filename_in)
            val_mismatch = .true.
          end if
          if (c14_atm_forcing_model_year /= c14_atm_forcing_model_year_in) then
            call document(subname, 'c14_atm_forcing_model_year', c14_atm_forcing_model_year)
            call document(subname, 'c14_atm_forcing_model_year_in', c14_atm_forcing_model_year_in)
            val_mismatch = .true.
          end if
          if (c14_atm_forcing_data_year /= c14_atm_forcing_data_year_in) then
            call document(subname, 'c14_atm_forcing_data_year', c14_atm_forcing_data_year)
            call document(subname, 'c14_atm_forcing_data_year_in', c14_atm_forcing_data_year_in)
            val_mismatch = .true.
          end if
        end select
      end if
      if (val_mismatch) then
        call exit_POP(sigAbort, 'c14_atm_forcing mismatch')
      end if
      return
    end if

    c14_atm_forcing_opt           = c14_atm_forcing_opt_in
    c14_atm_forcing_const         = c14_atm_forcing_const_in
    c14_atm_forcing_lat_band_vals = c14_atm_forcing_lat_band_vals_in
    c14_atm_forcing_filename      = c14_atm_forcing_filename_in
    c14_atm_forcing_model_year    = c14_atm_forcing_model_year_in
    c14_atm_forcing_data_year     = c14_atm_forcing_data_year_in

    call document(subname, 'c14_atm_forcing_opt', c14_atm_forcing_opt)

    select case (trim(c14_atm_forcing_opt))

    case ('const')

      call document(subname, 'c14_atm_forcing_const', c14_atm_forcing_const)

    case ('lat_bands')

      do lat_ind = 1, 3
        call document(subname, 'lat_ind', lat_ind)
        call document(subname, 'c14_atm_forcing_lat_band_vals(lat_ind)', c14_atm_forcing_lat_band_vals(lat_ind))
      end do

    case ('file')

      !-----------------------------------------------------------------------
      !  read in D14C data from file for option file
      !  ensure that read in data has dimlen=3 in 1st dimension
      !-----------------------------------------------------------------------

      call forcing_timeseries_init_dataset(c14_atm_forcing_filename, &
        varnames      = (/ 'Delta14co2_in_air' /), &
        model_year    = c14_atm_forcing_model_year, &
        data_year     = c14_atm_forcing_data_year, &
        taxmode_start = 'endpoint', &
        taxmode_end   = 'extrapolate', &
        dataset       = c14_atm_forcing_dataset)

      if (forcing_timeseries_dataset_var_size(c14_atm_forcing_dataset, varind=1, dim=1) /= 3) then
        call document(subname, 'size(c14_atm_forcing_data, dim=1)', &
          forcing_timeseries_dataset_var_size(c14_atm_forcing_dataset, varind=1, dim=1))
        call exit_POP(sigAbort, 'size(c14_atm_forcing_data, dim=1) /= 3')
      end if

      call document(subname, 'c14_atm_forcing_model_year', c14_atm_forcing_model_year)
      call document(subname, 'c14_atm_forcing_data_year', c14_atm_forcing_data_year)

    case default

      call exit_POP(sigAbort, 'unknown c14_atm_forcing_opt')

    end select

    first_call = .false.

    !-----------------------------------------------------------------------

  end subroutine c14_atm_forcing_init

  !*****************************************************************************

  subroutine c14_atm_forcing_update_data

    use forcing_timeseries_mod, only: forcing_timeseries_dataset_update_data

    if (trim(c14_atm_forcing_opt) == 'file') &
      call forcing_timeseries_dataset_update_data(c14_atm_forcing_dataset)

  end subroutine c14_atm_forcing_update_data

  !*****************************************************************************

  subroutine c14_atm_forcing_get_data(iblock, D14C)

    use forcing_timeseries_mod, only: forcing_timeseries_dataset_get_var

    integer (int_kind), intent(in) :: iblock

    real (r8), dimension(nx_block, ny_block), intent(out) :: D14C

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_get_data'

    real (r8) :: D14C_curr(3)

    !-----------------------------------------------------------------------

    select case (trim(c14_atm_forcing_opt))

    case ('const')

      D14C = c14_atm_forcing_const

    case ('lat_bands')

      call c14_atm_forcing_comp_D14C_bands(iblock, c14_atm_forcing_lat_band_vals, D14C)

    case ('file')

      !-----------------------------------------------------------------------
      !  Generate hemisphere values for current time step.
      !-----------------------------------------------------------------------

      call forcing_timeseries_dataset_get_var(c14_atm_forcing_dataset, varind=1, data_2d=D14C_curr)

      call c14_atm_forcing_comp_D14C_bands(iblock, D14C_curr, D14C)

    end select

  end subroutine c14_atm_forcing_get_data

  !*****************************************************************************

  subroutine c14_atm_forcing_comp_D14C_bands(iblock, D14C_curr, D14C)

! !DESCRIPTION:
!  Set atmospheric D14C to values based on latitude band
!  Global field of D14C is determined by:
!   -Northern Hemisphere value is used for 30N - 90 N
!   -Equator value is used for 30 S - 30 N
!   -Southern Hemisphere value is used for 30 S - 90 S

! !USES:

    use grid,                   only: TLATD

    integer (int_kind), intent(in) :: iblock

    real (r8), intent(in) :: D14C_curr(3)

    real (r8), dimension(nx_block, ny_block), intent(out) :: D14C

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_comp_D14C_bands'

    integer (int_kind) :: i, j

    !-----------------------------------------------------------------------

    do j = 1, ny_block
      do i = 1, nx_block
        if (TLATD(i,j,iblock) > 30.0_r8) then
          D14C(i,j) = D14C_curr(1)
        else if (TLATD(i,j,iblock) > -30.0_r8) then
          D14C(i,j) = D14C_curr(2)
        else
          D14C(i,j) = D14C_curr(3)
        end if
      end do
    end do

    !-----------------------------------------------------------------------

  end subroutine c14_atm_forcing_comp_D14C_bands

  !*****************************************************************************

end module c14_atm_forcing_mod
