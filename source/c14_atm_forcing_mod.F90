module c14_atm_forcing_mod

  ! !DESCRIPTION
  !  This module provides an interface for providing atmospheric boundary
  !  conditions for c14.

  use kinds_mod, only: char_len, r8, int_kind, log_kind
  use blocks, only: nx_block, ny_block
  use io_tools, only: document
  use exit_mod, only: sigAbort, exit_POP

  implicit none
  private

  public :: c14_atm_forcing_init
  public :: c14_atm_forcing_get_data

  !-----------------------------------------------------------------------
  !  module private data
  !-----------------------------------------------------------------------

  character(char_len) :: &
    c14_atm_forcing_opt             ! option for D14C varying or constant forcing

  real (r8) :: &
    c14_atm_forcing_const           !  atmospheric 14CO2 constant [permil]

  character (char_len), dimension(3) :: &
    c14_atm_forcing_filename        ! filenames for varying atm D14C (one each for NH, SH, EQ)

  integer (int_kind) :: &
    c14_atm_forcing_model_year,   & !  arbitrary model year
    c14_atm_forcing_data_year,    & !  year in data that corresponds to c14_atm_forcing_model_year
    c14_atm_forcing_data_nbval      !  number of values in c14_atm_forcing_filename

  real (r8), dimension(:,:) , allocatable :: &
    c14_atm_forcing_data,         & !  atmospheric DC14 values in datafile (sh, eq, nh, in permil)
    c14_atm_forcing_data_yr         !  date of atmospheric DC14 values in datafile (sh, eq, nh)

  integer (int_kind), dimension(:), allocatable :: &
    c14_atm_forcing_data_ind        ! data index for D14C data

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine c14_atm_forcing_init(caller, c14_atm_forcing_opt_in, &
      c14_atm_forcing_const_in, c14_atm_forcing_filename_in, &
      c14_atm_forcing_model_year_in, c14_atm_forcing_data_year_in)

    character(*), intent(in) :: &
      caller,                          & ! name of module calling c14_atm_forcing_init
      c14_atm_forcing_opt_in             ! option for D14C varying or constant forcing

    real (r8), intent(in) :: &
      c14_atm_forcing_const_in           !  atmospheric 14CO2 constant [permil]

    character (*), dimension(3), intent(in) :: &
      c14_atm_forcing_filename_in        ! filenames for varying atm D14C (one each for NH, SH, EQ)

    integer (int_kind), intent(in) :: &
      c14_atm_forcing_model_year_in,   & !  arbitrary model year
      c14_atm_forcing_data_year_in       !  year in data that corresponds to c14_atm_forcing_model_year_in

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(len=*), parameter  :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_init'

    logical(log_kind)            :: first_call = .true.
    logical(log_kind)            :: val_mismatch

    integer(int_kind)            :: il

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
      end if
      select case (c14_atm_forcing_opt)
      case ('const')
        if (c14_atm_forcing_const /= c14_atm_forcing_const_in) then
          call document(subname, 'c14_atm_forcing_const', c14_atm_forcing_const)
          call document(subname, 'c14_atm_forcing_const_in', c14_atm_forcing_const_in)
          val_mismatch = .true.
        end if
      case ('file')
        do il = 1, 3
          if (trim(c14_atm_forcing_filename(il)) /= trim(c14_atm_forcing_filename_in(il))) then
            call document(subname, 'il', il)
            call document(subname, 'c14_atm_forcing_filename(il)', c14_atm_forcing_filename(il))
            call document(subname, 'c14_atm_forcing_filename_in(il)', c14_atm_forcing_filename_in(il))
            val_mismatch = .true.
          end if
        end do
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
      case default
        call exit_POP(sigAbort, 'unknown c14_atm_forcing_opt ' /&
          &/ trim(c14_atm_forcing_opt) /&
          &/ ' in ' /&
          &/ trim(subname))
      end select
      if (val_mismatch) then
        call exit_POP(sigAbort, 'forcing mismatch in ' /&
          &/ trim(subname))
      end if
      return
    end if

    c14_atm_forcing_opt         = c14_atm_forcing_opt_in
    c14_atm_forcing_const       = c14_atm_forcing_const_in
    c14_atm_forcing_filename(:) = c14_atm_forcing_filename_in(:)
    c14_atm_forcing_model_year  = c14_atm_forcing_model_year_in
    c14_atm_forcing_data_year   = c14_atm_forcing_data_year_in

    call document(subname, 'c14_atm_forcing_opt', c14_atm_forcing_opt)

    select case (c14_atm_forcing_opt)

    case ('const')

      call document(subname, 'c14_atm_forcing_const', c14_atm_forcing_const)

    case ('file')

      !-----------------------------------------------------------------------
      !  read in D14C data from files for option file
      !-----------------------------------------------------------------------

      call c14_atm_forcing_read_file
      call document(subname, 'c14_atm_forcing_model_year', c14_atm_forcing_model_year)
      call document(subname, 'c14_atm_forcing_data_year', c14_atm_forcing_data_year)

    case default

      call exit_POP(sigAbort, 'unknown c14_atm_forcing_opt')

    end select

    first_call = .false.

    !-----------------------------------------------------------------------

  end subroutine c14_atm_forcing_init

  !*****************************************************************************

  subroutine c14_atm_forcing_read_file

! !DESCRIPTION:
!  Read atmospheric D14C data from file
!
!  Have the master_task do the following :
!     1) get length of data
!     2) allocate memory for data
!     3) read in data, checking for consistent lengths
!  Then, outside master_task conditional
!     1) broadcast length of data
!     2) have non-mastertasks allocate memory for data
!     3) broadcast data

    use io_types, only: stdout, get_unit
    use communicate, only: my_task, master_task
    use broadcast, only: broadcast_array, broadcast_scalar
    use domain, only: nblocks_clinic

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_read_file'

    integer (int_kind) :: &
      iunit,       & ! i/o unit number
      stat,        & ! i/o status code
      irec,        & ! counter for looping
      skiplines,   & ! number of comment lines at beginning of ascii file
      il             ! looping index

    character (char_len) :: &
      sglchr         ! variable to read characters from file into

    integer (int_kind) :: &
      c14_atm_forcing_data_nbval_tmp

    logical (log_kind) :: &
      nbval_mismatch

    !-----------------------------------------------------------------------

    call get_unit(iunit)

    !-----------------------------------------------------------------------
    ! ensure that three datafiles have same number of entries
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
      do il = 1, 3
        call document(subname, 'c14_atm_forcing_filename', c14_atm_forcing_filename(il))
      end do
      nbval_mismatch = .false.
      do il = 1, 3
        open(iunit, file=c14_atm_forcing_filename(il), status='old', iostat=stat)
        if (stat /= 0) then
          write(stdout,*) 'open failed for ', trim(c14_atm_forcing_filename(il))
          go to 99
        end if
        read(iunit, FMT=*, iostat=stat) skiplines, c14_atm_forcing_data_nbval_tmp
        if (stat /= 0) then
          write(stdout,*) '1st line read failed for ', trim(c14_atm_forcing_filename(il))
          go to 99
        end if
        close(iunit)
        if (il == 1) then
          c14_atm_forcing_data_nbval = c14_atm_forcing_data_nbval_tmp
        else
          if (c14_atm_forcing_data_nbval /= c14_atm_forcing_data_nbval_tmp) nbval_mismatch = .true.
        end if
      end do
    end if

    call broadcast_scalar(nbval_mismatch, master_task)
    if (nbval_mismatch) then
      call document(subname, 'D14C data files must all have the same number of values')
      call exit_POP(sigAbort, 'stopping in ' /&
                             &/ subname)
    end if

    call broadcast_scalar(c14_atm_forcing_data_nbval, master_task)
    allocate(c14_atm_forcing_data_yr(c14_atm_forcing_data_nbval,3))
    allocate(c14_atm_forcing_data(c14_atm_forcing_data_nbval,3))

    !-----------------------------------------------------------------------
    !  READ in C14 data from files - three files, for SH, EQ, NH
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
      do il = 1, 3
        open(iunit, file=c14_atm_forcing_filename(il), status='old', iostat=stat)
        if (stat /= 0) then
          write(stdout,*) 'open failed for ', trim(c14_atm_forcing_filename(il))
          go to 99
        end if
        read(iunit, FMT=*, iostat=stat) skiplines, c14_atm_forcing_data_nbval_tmp
        if (stat /= 0) then
          write(stdout,*) '1st line read failed for ', trim(c14_atm_forcing_filename(il))
          go to 99
        end if
        do irec = 1, skiplines
          read(iunit, FMT=*, iostat=stat) sglchr
          if (stat /= 0) then
            write(stdout, fmt=*) 'skipline read failed for ', trim(c14_atm_forcing_filename(il))
            go to 99
          end if
        end do
        do irec = 1, c14_atm_forcing_data_nbval
          read(iunit, FMT=*, iostat=stat) c14_atm_forcing_data_yr(irec,il), c14_atm_forcing_data(irec,il)
          if (stat /= 0) then
            write(stdout, fmt=*) 'data read failed for ', trim(c14_atm_forcing_filename(il))
            go to 99
          end if
        end do
        close(iunit)
      end do
    end if

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                           &/ subname)

    !---------------------------------------------------------------------
    !  Broadcast the variables to other tasks beside master_task
    !---------------------------------------------------------------------

    call broadcast_array(c14_atm_forcing_data, master_task)
    call broadcast_array(c14_atm_forcing_data_yr, master_task)

    !---------------------------------------------------------------------
    !  set initial index into forcing_data arrays
    !---------------------------------------------------------------------

    allocate(c14_atm_forcing_data_ind(nblocks_clinic))
    c14_atm_forcing_data_ind(:) = -1

    !---------------------------------------------------------------------

  end subroutine c14_atm_forcing_read_file

  !*****************************************************************************

  subroutine c14_atm_forcing_get_data(iblock, D14C)

    integer (int_kind), intent(in) :: iblock

    real (r8), dimension(nx_block, ny_block), intent(out) :: D14C

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_get_data'

    !-----------------------------------------------------------------------

    select case (c14_atm_forcing_opt)

    case ('const')

      D14C = c14_atm_forcing_const

    case ('file')

      call c14_atm_forcing_comp_varying_D14C(iblock, D14C)

    case default

      call document(subname, 'c14_atm_forcing_opt', c14_atm_forcing_opt)
      call exit_POP(sigAbort, 'unknown c14_atm_forcing_opt')

    end select

  end subroutine c14_atm_forcing_get_data

  !*****************************************************************************

  subroutine c14_atm_forcing_comp_varying_D14C(iblock, D14C)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of CO2 when temporarily
!  varying data is read from files
!  1. Linearly interpolate hemispheric values to current time step
!  2. Make global field of D14C, determined by:
!   -Northern Hemisphere value is used for 20N - 90 N
!   -Southern Hemisphere value is used for 20 S - 90 S
!   -Equator value is used for 20 S- 20 N

! !USES:

    use grid, only: TLATD
    use constants, only: c0, c1
    use time_management, only: iyear, iday_of_year, frac_day, days_in_year

    integer (int_kind), intent(in) :: iblock

    real (r8), dimension(nx_block, ny_block), intent(out) :: D14C

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'c14_atm_forcing_mod:c14_atm_forcing_comp_varying_D14C'

    integer (int_kind) :: &
      data_ind,     & ! c14_atm_forcing_data_ind for this block
      i, j, il        ! loop indices

    real (r8) :: &
      model_date,   & ! date of current model timestep
      mapped_date,  & ! model_date mapped to data timeline
      weight,       & ! weighting for temporal interpolation
      D14c_curr_sh, & ! current atmospheric D14C value for SH (interpolated from data to model date)
      D14c_curr_nh, & ! current atmospheric D14C value for NH (interpolated from data to model date)
      D14c_curr_eq    ! current atmospheric D14C value for EQ (interpolated from data to model date)


    !-----------------------------------------------------------------------
    !  Generate mapped_date and check to see if it is too large.
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    mapped_date = model_date - c14_atm_forcing_model_year + c14_atm_forcing_data_year
    do il=1,3
      if (mapped_date >= c14_atm_forcing_data_yr(c14_atm_forcing_data_nbval,il)) then
        call exit_POP(sigAbort, 'model date maps to date after end of D14C data in files.')
      end if
    end do

    !-----------------------------------------------------------------------
    !  Set atmospheric D14C concentrations to zero before D14C record begins
    !-----------------------------------------------------------------------

    if (mapped_date < c14_atm_forcing_data_yr(1,1)) then
      D14C(:,:) = c0
      c14_atm_forcing_data_ind(iblock) = 1
      call document(subname, 'Model date less than start of D14C data --> D14C=0')
      return
    end if

    !-----------------------------------------------------------------------
    !  On first time step, perform linear search to find data_ind.
    !-----------------------------------------------------------------------

    if (c14_atm_forcing_data_ind(iblock) == -1) then
      do data_ind = c14_atm_forcing_data_nbval-1,1,-1
        if (mapped_date >= c14_atm_forcing_data_yr(data_ind,1)) exit
      end do
      c14_atm_forcing_data_ind(iblock) = data_ind
    end if

    !-----------------------------------------------------------------------
    !  See if c14_atm_forcing_data_ind need to be updated,
    !  but do not set it to c14_atm_forcing_data_nbval.
    !-----------------------------------------------------------------------

    data_ind = c14_atm_forcing_data_ind(iblock)
    if (data_ind < c14_atm_forcing_data_nbval-1) then
      if (mapped_date >= c14_atm_forcing_data_yr(data_ind+1,1)) then
        data_ind = data_ind + 1
        c14_atm_forcing_data_ind(iblock) = data_ind
      end if
    end if

    !-----------------------------------------------------------------------
    !  Generate hemisphere values for current time step.
    !-----------------------------------------------------------------------

    weight = (mapped_date - c14_atm_forcing_data_yr(data_ind,1)) &
      / (c14_atm_forcing_data_yr(data_ind+1,1) - c14_atm_forcing_data_yr(data_ind,1))

    d14c_curr_sh = weight * c14_atm_forcing_data(data_ind+1,1) \
      + (c1-weight) * c14_atm_forcing_data(data_ind,1)
    d14c_curr_eq = weight * c14_atm_forcing_data(data_ind+1,2) \
      + (c1-weight) * c14_atm_forcing_data(data_ind,2)
    d14c_curr_nh = weight * c14_atm_forcing_data(data_ind+1,3) \
      + (c1-weight) * c14_atm_forcing_data(data_ind,3)

    !-----------------------------------------------------------------------
    !  Merge hemisphere values for D14C
    !      -Northern Hemisphere value is used for >20N - 90 N
    !      -Southern Hemisphere value is used for >20 S - 90 S
    !      -Equatorial value is used for 20 S to 20 N
    !-----------------------------------------------------------------------

    do j = 1, ny_block
      do i = 1, nx_block
        if (TLATD(i,j,iblock) < -20.0_r8) then
          D14C(i,j) = d14c_curr_sh
        else if (TLATD(i,j,iblock) > 20.0_r8) then
          D14C(i,j) = d14c_curr_nh
        else
          D14C(i,j) = d14c_curr_eq
        end if
      end do
    end do

    !-----------------------------------------------------------------------

  end subroutine c14_atm_forcing_comp_varying_D14C

  !*****************************************************************************

end module c14_atm_forcing_mod
