module forcing_timeseries_mod

  ! !DESCRIPTION
  !  This module provides an interface for reading and interpolating 1d and 2d timeseries

  use kinds_mod, only: int_kind, char_len, r8
  use io_tools,  only: document
  use exit_mod,  only: sigAbort, exit_POP

  implicit none
  private

  public :: forcing_timeseries_dataset

  public :: forcing_timeseries_init_dataset
  public :: forcing_timeseries_dataset_var_size
  public :: forcing_timeseries_dataset_update_data
  public :: forcing_timeseries_dataset_get_var

  !*****************************************************************************

  integer (int_kind), parameter :: forcing_timeseries_taxmode_endpoint    = 1
  integer (int_kind), parameter :: forcing_timeseries_taxmode_extrapolate = 2

  real (r8), parameter :: max_yr_extension = 4.0_r8 ! maximum number of years that a timeseries will be extrapolated

  type :: forcing_timeseries_var
    character (char_len)   :: name
    integer (int_kind)     :: rank
    real (r8), allocatable :: data_1d(:)        ! dimension is time
    real (r8)              :: data_1d_curr      ! data_1d interpolated to the current model time
    real (r8), allocatable :: data_2d(:,:)      ! dimensions are (*,time)
    real (r8), allocatable :: data_2d_curr(:)   ! data_2d interpolated to the current model time
  end type forcing_timeseries_var

  ! a dataset is a collection of variables from one file, all on the same time dimension

  type :: forcing_timeseries_dataset
    character (char_len)                       :: filename       ! name of file containing dataset
    real (r8), allocatable                     :: time_yr(:)     ! time values [years]
    real (r8)                                  :: model_year     ! arbitrary model year
    real (r8)                                  :: data_year      ! year in time_yr corresponding to model_year
    integer (int_kind)                         :: taxmode_start  ! taxmode at the start of timeseries
    integer (int_kind)                         :: taxmode_end    ! taxmode at the end of timeseries
    type (forcing_timeseries_var), allocatable :: vars(:)
    integer (int_kind)                         :: data_ind       ! current time index
    real (r8)                                  :: last_update_model_date ! model_date from most recent update
  end type forcing_timeseries_dataset

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine forcing_timeseries_init_dataset(filename, varnames, &
    model_year, data_year, taxmode_start, taxmode_end, dataset)

! !DESCRIPTION:
!  Initialize dataset variable, populating it with data from varnames in filename

    character (*), intent(in) :: &
      filename           ! name of file containing timeseries

    character (*), intent(in) :: &
      varnames(:)         ! names of variable being read in

    integer (int_kind), intent(in) :: &
      model_year,       & ! arbitrary model year
      data_year           ! year in time_yr corresponding to model_year

    character (*), intent(in) :: &
      taxmode_start,    & ! taxmode at the start of timeseries
      taxmode_end         ! taxmode at the end of timeseries

    type (forcing_timeseries_dataset), intent(out) :: &
      dataset             ! forcing_timeseries_dataset being initialized

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_init_dataset'

    integer (int_kind) :: &
      nvars,            & ! number of vars in dataset
      varind              ! var loop index

    !-----------------------------------------------------------------------

    call document(subname, 'filename', filename)
    dataset%filename = trim(filename)

    call document(subname, 'model_year', model_year)
    dataset%model_year = model_year

    call document(subname, 'data_year', data_year)
    dataset%data_year = data_year

    call document(subname, 'taxmode_start', taxmode_start)
    select case (trim(taxmode_start))
    case ('endpoint')
      dataset%taxmode_start = forcing_timeseries_taxmode_endpoint
    case ('extrapolate')
      dataset%taxmode_start = forcing_timeseries_taxmode_extrapolate
    case default
      call exit_POP(sigAbort, 'unknown taxmode_start')
    end select

    call document(subname, 'taxmode_end', taxmode_end)
    select case (trim(taxmode_end))
    case ('endpoint')
      dataset%taxmode_end = forcing_timeseries_taxmode_endpoint
    case ('extrapolate')
      dataset%taxmode_end = forcing_timeseries_taxmode_extrapolate
    case default
      call exit_POP(sigAbort, 'unknown taxmode_end')
    end select

    nvars = size(varnames)
    allocate(dataset%vars(nvars))

    do varind = 1, nvars
      dataset%vars(varind)%name = trim(varnames(varind))
    end do

    call forcing_timeseries_read_dataset(dataset)

    dataset%data_ind = -1
    dataset%last_update_model_date = -1

  end subroutine forcing_timeseries_init_dataset

  !*****************************************************************************

  subroutine forcing_timeseries_read_dataset(dataset)

! !DESCRIPTION:
!  Read timeseries data, and corresponding years, from a netCDF file
!  Supported ranks are 1, 2

    use netcdf
    use broadcast, only: broadcast_array, broadcast_scalar
    use communicate, only: my_task, master_task

    type (forcing_timeseries_dataset), intent(inout) :: &
      dataset             ! forcing_timeseries_dataset being read in

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_read_dataset'

    character (len=char_len) :: &
      timename,         & ! name of 'time' dimension
      units,            & ! netCDF units attribute
      valid_unit_pref,  & ! accepted value for prefix in time units
      ref_date            ! reference date in time units

    integer (int_kind) :: &
      nvars,            & ! number of vars in dataset
      stat,             & ! status of netCDF call
      ncid,             & ! netCDF file id
      varind,           & ! var loop index
      varid,            & ! netCDF variable id
      last_dimid,       & ! dimid of last dimension of var
      dim,              & ! dimension loop index
      time_dimlen,      & ! dimlen of last dimension of var
      time_ndims,       & ! number of dimensions for time
      timeind,          & ! time index
      pos                 ! character position in a string

    integer (int_kind) :: &
      var_dimids(2),    & ! netCDF dimension ids for varname
      var_dimlens(2)      ! netCDF dimension lengths for varname

    real (r8) :: &
      ref_year            ! year in ref_date

    !-----------------------------------------------------------------------

    nvars = size(dataset%vars)

    !-----------------------------------------------------------------------
    !  perform netCDF I/O on master_task
    !  jump out of master_task conditional if an error is encountered
    !-----------------------------------------------------------------------

    if (my_task == master_task) then

      stat = nf90_open(dataset%filename, 0, ncid)
      if (stat /= 0) then
        call document(subname, 'nf_open: ', nf90_strerror(stat))
        go to 99
      end if

      do varind = 1, nvars
        associate (var => dataset%vars(varind))
          call document(subname, 'reading timeseries ' /&
            &/ trim(var%name) /&
            &/ ' from ' /&
            &/ trim(dataset%filename))

          stat = nf90_inq_varid(ncid, var%name, varid)
          if (stat /= 0) then
            call document(subname, 'nf_inq_varid', nf90_strerror(stat))
            go to 99
          end if

          stat = nf90_inquire_variable(ncid, varid, ndims=var%rank)
          if (stat /= 0) then
            call document(subname, 'nf_inq_varndims', nf90_strerror(stat))
            go to 99
          end if
          if ((var%rank < 1) .or. (var%rank > 2)) then
            call document(subname, 'non-supported rank', var%rank)
            stat = 1
            go to 99
          end if

          stat = nf90_inquire_variable(ncid, varid, dimids=var_dimids)
          if (stat /= 0) then
            call document(subname, 'nf_inq_vardimid', nf90_strerror(stat))
            go to 99
          end if

          ! if varind==1, store the last ('time') dimension
          ! if varind>1, confirm that the last ('time') dimension agrees with last_dimid

          if (varind == 1) then
            last_dimid = var_dimids(var%rank)
          else
            if (var_dimids(var%rank) .ne. last_dimid) then
              call document(subname, 'last dimid does not match last dimid from previous vars')
              stat = 1
              go to 99
            end if
          end if

          ! get dimension lengths

          do dim = 1, var%rank
            call document(subname, 'dim', dim)
            stat = nf90_inquire_dimension(ncid, var_dimids(dim), len=var_dimlens(dim))
            if (stat /= 0) then
              call document(subname, 'nf_inq_dimlen(dim)', nf90_strerror(stat))
              go to 99
            end if
            call document(subname, 'dimlen', var_dimlens(dim))
          end do

          ! allocate space for data, and read in data

          select case (var%rank)
          case (1)
            allocate(var%data_1d(var_dimlens(1)))
            stat = nf90_get_var(ncid, varid, var%data_1d)
          case (2)
            allocate(var%data_2d(var_dimlens(1), var_dimlens(2)))
            stat = nf90_get_var(ncid, varid, var%data_2d)
          end select

          if (stat /= 0) then
            call document(subname, 'nf_get_var', nf90_strerror(stat))
            go to 99
          end if
        end associate
      end do ! do varind = 1, nvars

      ! read in 'time'
      ! this is often named 'time', but we don't assume it to be
      ! confirm that 'time' is a 1d var, on the same dimension as the last dimension of varname

      timename = ' '
      stat = nf90_inquire_dimension(ncid, last_dimid, name=timename)
      if (stat /= 0) then
        call document(subname, 'nf_inq_dimname(last_dimid)', nf90_strerror(stat))
        go to 99
      end if

      call document(subname, 'timename', timename)

      stat = nf90_inq_varid(ncid, timename, varid)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varid for time: ', nf90_strerror(stat))
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, ndims=time_ndims)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varndims for time: ', nf90_strerror(stat))
        go to 99
      end if
      if (time_ndims /= 1) then
        call document(subname, 'non-supported time_ndims', time_ndims)
        stat = 1
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, dimids=var_dimids)
      if (stat /= 0) then
        call document(subname, 'nf_inq_vardimid for time: ', nf90_strerror(stat))
        go to 99
      end if

      if (var_dimids(1) /= last_dimid) then
        call document(subname, 'time dimid /= last_dimid')
        stat = 1
        go to 99
      end if

      stat = nf90_inquire_dimension(ncid, var_dimids(1), len=time_dimlen)
      if (stat /= 0) then
        call document(subname, 'nf_inq_dimlen(time_dimids(1))', nf90_strerror(stat))
        go to 99
      end if

      allocate(dataset%time_yr(time_dimlen))
      stat = nf90_get_var(ncid, varid, dataset%time_yr)
      if (stat /= 0) then
        call document(subname, 'nf_get_var for time: ', nf90_strerror(stat))
        go to 99
      end if

      ! ensure that 'time' is increasing

      do timeind = 1, time_dimlen-1
        if (dataset%time_yr(timeind+1) <= dataset%time_yr(timeind)) then
          call document(subname, 'timeind', timeind)
          call document(subname, 'time_yr(timeind)', dataset%time_yr(timeind))
          call document(subname, 'time_yr(timeind+1)', dataset%time_yr(timeind+1))
          call document(subname, 'time_yr values violate assumption that time is strictly increasing')
          stat = 1
          go to 99
        end if
      end do

      ! get units of 'time' and convert into absolute years, by adding ref_year

      units = ' '
      stat = nf90_get_att(ncid, varid, 'units', units)
      if (stat /= 0) then
        call document(subname, 'nf_get_att for time@units: ', nf90_strerror(stat))
        go to 99
      end if

      ! find character position of reference year in units string
      ! supported units are 'years since ...' and 'common_years since ...'

      units = adjustl(units) ! remove leading space(s)

      pos = 0
      valid_unit_pref = 'years since'
      if (index(units, trim(valid_unit_pref)) == 1) then
        pos = len_trim(valid_unit_pref) + 1
      end if
      valid_unit_pref = 'common_years since'
      if (index(units, trim(valid_unit_pref)) == 1) then
        pos = len_trim(valid_unit_pref) + 1
      end if
      if (pos == 0) then
        call document(subname, 'non-supported time@units', units)
        stat = 1
        go to 99
      end if

      ref_date = adjustl(units(pos:len_trim(units))) ! remove prefix and space(s) immediately after it
      call ref_date_to_ref_year(ref_date, ref_year)

      dataset%time_yr(:) = dataset%time_yr(:) + ref_year
    end if       ! my_task == master_task

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                           &/ subname)

    !---------------------------------------------------------------------
    !  Broadcast the variables to other tasks beside master_task
    !  allocate space for _curr value
    !---------------------------------------------------------------------

    call broadcast_scalar(time_dimlen, master_task)
    if (my_task /= master_task) allocate(dataset%time_yr(time_dimlen))
    call broadcast_array(dataset%time_yr, master_task)

    do varind = 1, nvars
      associate (var => dataset%vars(varind))
        call broadcast_scalar(var%rank, master_task)
        select case (var%rank)
        case (1)
          if (my_task == master_task) var_dimlens(1) = size(var%data_1d)
          call broadcast_scalar(var_dimlens(1), master_task)
          if (my_task /= master_task) allocate(var%data_1d(var_dimlens(1)))
          call broadcast_array(var%data_1d, master_task)
        case (2)
          if (my_task == master_task) var_dimlens(1:2) = shape(var%data_2d)
          call broadcast_array(var_dimlens(1:2), master_task)
          if (my_task /= master_task) allocate(var%data_2d(var_dimlens(1), var_dimlens(2)))
          call broadcast_array(var%data_2d, master_task)
          allocate(var%data_2d_curr(var_dimlens(1)))
        end select
      end associate
    end do

  end subroutine forcing_timeseries_read_dataset

  !*****************************************************************************

  subroutine ref_date_to_ref_year(ref_date, ref_year)

    use constants, only: c1

    character (*), intent(in) :: ref_date

    real (r8), intent(out) :: ref_year

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:ref_date_to_ref_year'

    real (r8), dimension(12), parameter :: day0 = & ! starting day of each month
      (/   0.0_r8,  31.0_r8,  59.0_r8,  90.0_r8, 120.0_r8, 151.0_r8, &
         181.0_r8, 212.0_r8, 243.0_r8, 273.0_r8, 304.0_r8, 334.0_r8 /)

    integer (int_kind) :: &
      ref_date_len, & ! length of input string
      pos, dpos,    & ! character positions in string
      mon,          & ! month number
      day             ! day number

    !-----------------------------------------------------------------------

    ref_date_len = len(ref_date)

    ! search for a separator between the year and month
    pos = 1
    dpos = scan(ref_date(pos:ref_date_len), '-')

    ! read the year
    if (dpos == 0) then
      read(ref_date(pos:ref_date_len),*) ref_year
    else
      read(ref_date(pos:pos+dpos-2),*) ref_year
    end if

    ! no separator found, so there is nothing after the year
    if (dpos == 0) return

    ! search for a separator between the month and day
    pos = pos+dpos
    dpos = scan(ref_date(pos:ref_date_len), '-')

    ! read the month and add proper fraction to ref_year
    if (dpos == 0) then
      read(ref_date(pos:ref_date_len),*) mon
    else
      read(ref_date(pos:pos+dpos-2),*) mon
    end if
    ref_year = ref_year + day0(mon) / 365.0_r8

    ! no separator found, so there is nothing after the month
    if (dpos == 0) return

    ! search for a separator between the day and anything else
    pos = pos+dpos
    dpos = scan(ref_date(pos:ref_date_len), '- ')

    ! read the day and add proper fraction to ref_year
    if (dpos == 0) then
      read(ref_date(pos:ref_date_len),*) day
    else
      read(ref_date(pos:pos+dpos-2),*) day
    end if
    ref_year = ref_year + (day-c1) / 365.0_r8

    ! ignore anything potentially after the day

  end subroutine ref_date_to_ref_year

  !*****************************************************************************

  function forcing_timeseries_dataset_var_size(dataset, varind, dim)

! !DESCRIPTION:
!  return the size of a variable in a dataset

    type (forcing_timeseries_dataset), intent(in) :: &
      dataset             ! forcing_timeseries_dataset that inquiry is about

    integer (int_kind), intent(in) :: &
      varind,           & ! index of variable that inquiry is about
      dim                 ! dimension of variable that inquiry is about

    integer (int_kind) :: &
      forcing_timeseries_dataset_var_size ! result

    !---------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_dataset_var_size'

    !-----------------------------------------------------------------------

    associate (var => dataset%vars(varind))

      ! ensure consistency between arguments and var corresponding to varind

      if ((dim < 1) .or. (dim > var%rank)) then
        call document(subname, 'varname', var%name)
        call document(subname, 'var%rank', var%rank)
        call document(subname, 'dim', dim)
        call exit_POP(sigAbort, 'requested dim out of range for var')
      end if

      select case (var%rank)
      case (1)
        forcing_timeseries_dataset_var_size = size(var%data_1d, dim=dim)
      case (2)
        forcing_timeseries_dataset_var_size = size(var%data_2d, dim=dim)
      end select

    end associate

    !---------------------------------------------------------------------

  end function forcing_timeseries_dataset_var_size

  !*****************************************************************************

  subroutine forcing_timeseries_dataset_update_data(dataset)

! !DESCRIPTION:
!  update the vars in a dataset, interpolating in time, storing the result in curr var components

    use time_management, only: iyear, iday_of_year, frac_day, days_in_year
    use constants,       only: c0

    type (forcing_timeseries_dataset), intent(inout) :: &
      dataset             ! forcing_timeseries_dataset containing data

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_dataset_update_data'

    integer (int_kind) :: &
      nbval,        & ! number of values in timeseries
      varind,       & ! variable index
      data_ind,     & ! data_ind for this block
      w0_ind,       & ! time index into timeseries of data values corresponding to weight=0
      w1_ind          ! time index into timeseries of data values corresponding to weight=1

    real (r8) :: &
      model_date,   & ! date of current model timestep
      mapped_date,  & ! model_date mapped to data timeline
      weight          ! weighting for temporal interpolation

    !-----------------------------------------------------------------------
    !  Generate mapped_date
    !    return immediately if update_data has already been called for the current model_date
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    if (model_date .eq. dataset%last_update_model_date) return
    dataset%last_update_model_date = model_date

    mapped_date = model_date - dataset%model_year + dataset%data_year

    nbval = size(dataset%time_yr)

    if (mapped_date < dataset%time_yr(1)) then

      !-----------------------------------------------------------------------
      !  Set interpolation weights and indices if mapped_date is before first date in file
      !-----------------------------------------------------------------------

      if ((dataset%taxmode_start == forcing_timeseries_taxmode_extrapolate) &
         .and. (mapped_date < dataset%time_yr(1) - max_yr_extension)) then
        call document(subname, 'filename', dataset%filename)
        call document(subname, 'model_date', model_date)
        call document(subname, 'mapped_date', mapped_date)
        call document(subname, 'time_yr(1)', dataset%time_yr(1))
        call exit_POP(sigAbort, 'model date maps too far before first date in file')
      end if

      select case (dataset%taxmode_start)
      case (forcing_timeseries_taxmode_endpoint)
        weight = c0
      case (forcing_timeseries_taxmode_extrapolate)
        weight = (mapped_date - dataset%time_yr(1)) / (dataset%time_yr(2) - dataset%time_yr(1))
      end select

      w0_ind = 1
      w1_ind = 2

      dataset%data_ind = 1

    else if (mapped_date < dataset%time_yr(nbval)) then

      !-----------------------------------------------------------------------
      !  Set interpolation weights and indices if mapped_date is between first and last dates in file
      !-----------------------------------------------------------------------

      if (dataset%data_ind == -1) then

        !-----------------------------------------------------------------------
        !  If data_ind is not set, perform linear search to find data_ind.
        !-----------------------------------------------------------------------

        do data_ind = nbval-1,1,-1
          if (mapped_date >= dataset%time_yr(data_ind)) exit
        end do
        dataset%data_ind = data_ind

      else

        !-----------------------------------------------------------------------
        !  See if data_ind need to be updated, but do not set it to nbval.
        !-----------------------------------------------------------------------

        data_ind = dataset%data_ind
        if (data_ind < nbval-1) then
          if (mapped_date >= dataset%time_yr(data_ind+1)) then
            data_ind = data_ind + 1
            dataset%data_ind = data_ind
          end if
        end if

      end if

      weight = (mapped_date - dataset%time_yr(data_ind)) / (dataset%time_yr(data_ind+1) - dataset%time_yr(data_ind))

      w0_ind = data_ind
      w1_ind = data_ind+1

    else

      !-----------------------------------------------------------------------
      !  Set interpolation weights and indices if mapped_date is after last date in file
      !-----------------------------------------------------------------------

      if ((dataset%taxmode_end == forcing_timeseries_taxmode_extrapolate) &
         .and. (mapped_date > dataset%time_yr(nbval) + max_yr_extension)) then
        call document(subname, 'filename', dataset%filename)
        call document(subname, 'model_date', model_date)
        call document(subname, 'mapped_date', mapped_date)
        call document(subname, 'nbval', nbval)
        call document(subname, 'time_yr(nbval)', dataset%time_yr(nbval))
        call exit_POP(sigAbort, 'model date maps too far beyond last date in file')
      end if

      select case (dataset%taxmode_end)
      case (forcing_timeseries_taxmode_endpoint)
        weight = c0
      case (forcing_timeseries_taxmode_extrapolate)
        weight = (dataset%time_yr(nbval) - mapped_date) / (dataset%time_yr(nbval) - dataset%time_yr(nbval-1))
      end select

      w0_ind = nbval
      w1_ind = nbval-1

      dataset%data_ind = nbval-1

    end if

    !-----------------------------------------------------------------------
    !  Perform linear interpolation
    !-----------------------------------------------------------------------

    do varind = 1, size(dataset%vars)
      associate (var => dataset%vars(varind))
        select case (var%rank)
        case (1)
          var%data_1d_curr    = var%data_1d(w0_ind)   + weight * (var%data_1d(w1_ind)   - var%data_1d(w0_ind))
        case (2)
          var%data_2d_curr(:) = var%data_2d(:,w0_ind) + weight * (var%data_2d(:,w1_ind) - var%data_2d(:,w0_ind))
        end select
      end associate
    end do

    !---------------------------------------------------------------------

  end subroutine forcing_timeseries_dataset_update_data

  !*****************************************************************************

  subroutine forcing_timeseries_dataset_get_var(dataset, varind, data_1d, data_2d)

! !DESCRIPTION:
!  return the current value of a variable in a dataset

    type (forcing_timeseries_dataset), intent(inout) :: &
      dataset             ! forcing_timeseries_dataset containing data

    integer (int_kind), intent(in) :: &
      varind              ! index of variable that inquiry is about

    real (r8), intent(out), optional :: &
      data_1d,          & ! where values from 1d timeseries are stored
      data_2d(:)          ! where values from 2d timeseries are stored

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_dataset_get_var'

    !-----------------------------------------------------------------------

    associate (var => dataset%vars(varind))

      ! ensure consistency between arguments and var corresponding to varind

      select case (var%rank)
      case (1)
        if (.not. present(data_1d)) then
          call document(subname, 'data_1d argument not passed in for rank=1 var', var%name)
          call exit_POP(sigAbort, 'data rank mismatch')
        end if
        if (present(data_2d)) then
          call document(subname, 'data_2d argument passed in for rank=1 var', var%name)
          call exit_POP(sigAbort, 'data rank mismatch')
        end if
      case (2)
        if (present(data_1d)) then
          call document(subname, 'data_1d argument passed in for rank=2 var', var%name)
          call exit_POP(sigAbort, 'data rank mismatch')
        end if
        if (.not. present(data_2d)) then
          call document(subname, 'data_2d argument not passed in for rank=2 var', var%name)
          call exit_POP(sigAbort, 'data rank mismatch')
        end if
        if (size(data_2d) /= size(var%data_2d_curr)) then
          call document(subname, 'var', var%name)
          call document(subname, 'size(var%data_2d_curr)', size(var%data_2d_curr))
          call document(subname, 'size(data_2d)', size(data_2d))
          call exit_POP(sigAbort, 'dimension size mismatch')
        end if
      end select

      if (present(data_1d)) data_1d    = var%data_1d_curr

      if (present(data_2d)) data_2d(:) = var%data_2d_curr(:)

    end associate

    !---------------------------------------------------------------------

  end subroutine forcing_timeseries_dataset_get_var

  !*****************************************************************************

end module forcing_timeseries_mod
