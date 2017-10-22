module forcing_timeseries_mod

  ! !DESCRIPTION
  !  This module provides an interface for reading 1d and 2d timeseries

  use kinds_mod, only: r8, int_kind

  implicit none
  private

  public :: forcing_timeseries_read_file

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine forcing_timeseries_read_file(filename, varname, &
    data_nbval, data_yr, data_1d, data_2d)

! !DESCRIPTION:
!  Read timeseries data, and corresponding year, from a netCDF file
!  Supported ranks are 1, 2

    use netcdf
    use kinds_mod, only: char_len
    use broadcast, only: broadcast_array, broadcast_scalar
    use communicate, only: my_task, master_task
    use exit_mod, only: sigAbort, exit_POP
    use io_tools, only: document

    character (*), intent(in) :: &
      filename,        & ! name of file containing timeseries
      varname            ! name of variable being read in

    integer (int_kind), intent(out) :: &
      data_nbval         ! lenght of timeseries

    real (r8), dimension(:), allocatable, intent(out) :: &
      data_yr            ! time values corresponding to varname, in years

    real (r8), dimension(:), allocatable, intent(out), optional :: &
      data_1d            ! where read in data is to be stored

    real (r8), dimension(:,:), allocatable, intent(out), optional :: &
      data_2d            ! where read in data is to be stored

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'forcing_timeseries_mod:forcing_timeseries_read_file'

    character (len=char_len) :: &
      timename,        & ! name of 'time' dimension
      units,           & ! netCDF units attribute
      valid_unit_pref, & ! accepted value for prefix in time units
      ref_date           ! reference date in time units

    integer (int_kind) :: &
      data_rank,       & ! rank of data argument
      stat,            & ! status of netCDF call
      ncid,            & ! netCDF file id
      varid,           & ! netCDF variable id
      ndims,           & ! number of dimensions for varid
      dim,             & ! dimension loop index
      dimlen,          & ! netCDF dimension length
      pos                ! character position in a string

    integer (int_kind), dimension(2) :: &
      data_dimids,     & ! netCDF dimension ids for varname
      time_dimids,     & ! netCDF dimension ids for time
      data_dimlens       ! netCDF dimension lengths for varname

    real (r8) :: &
      ref_year           ! year in ref_date

    !-----------------------------------------------------------------------

    call document(subname, 'reading timeseries ' /&
      &/ trim(varname) /&
      &/ ' from ' /&
      &/ trim(filename))

    data_rank = 0

    if (present(data_1d)) then
      data_rank = 1
    end if

    if (present(data_2d)) then
      if (data_rank /= 0) then
        call exit_POP(sigAbort, 'cannot pass more than 1 data argument to ' /&
          &/ subname)
      end if
      data_rank = 2
    end if

    if (data_rank == 0) then
      call exit_POP(sigAbort, 'must pass one data argument to ' /&
        &/ subname)
    end if

    call document(subname, 'data_rank', data_rank)

    !-----------------------------------------------------------------------
    !  perform netCDF I/O on master_task
    !  jump out of master_task conditional if an error is encountered
    !-----------------------------------------------------------------------

    if (my_task == master_task) then

      stat = nf90_open(filename, 0, ncid)
      if (stat /= 0) then
        call document(subname, 'nf_open: ', nf90_strerror(stat))
        go to 99
      end if

      stat = nf90_inq_varid(ncid, varname, varid)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varid', nf90_strerror(stat))
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, ndims=ndims)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varndims', nf90_strerror(stat))
        go to 99
      end if
      if (ndims /= data_rank) then
        call document(subname, 'ndims', ndims)
        call document(subname, 'ndims /= data_rank')
        stat = 1
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, dimids=data_dimids)
      if (stat /= 0) then
        call document(subname, 'nf_inq_vardimid', nf90_strerror(stat))
        go to 99
      end if

      ! get dimension lengths

      do dim = 1, data_rank
        call document(subname, 'dim', dim)
        stat = nf90_inquire_dimension(ncid, data_dimids(dim), len=data_dimlens(dim))
        if (stat /= 0) then
          call document(subname, 'nf_inq_dimlen(dim)', nf90_strerror(stat))
          go to 99
        end if
        call document(subname, 'dimlen', data_dimlens(dim))
      end do

      ! allocate space for data_yr and data, and read in data

      select case (data_rank)
      case (1)
        allocate(data_yr(data_dimlens(1)))
        allocate(data_1d(data_dimlens(1)))
        stat = nf90_get_var(ncid, varid, data_1d)
      case (2)
        allocate(data_yr(data_dimlens(2)))
        allocate(data_2d(data_dimlens(1), data_dimlens(2)))
        stat = nf90_get_var(ncid, varid, data_2d)
      end select

      if (stat /= 0) then
        call document(subname, 'nf_get_var', nf90_strerror(stat))
        go to 99
      end if

      ! get name of last dimension of varname
      ! this is often 'time', but we don't assume it to be
      ! confirm that 'time' is a 1d var, on the same dimension as the last dimension of varname

      timename = ' '
      stat = nf90_inquire_dimension(ncid, data_dimids(data_rank), name=timename)
      if (stat /= 0) then
        call document(subname, 'nf_inq_dimname(data_rank)', nf90_strerror(stat))
        go to 99
      end if

      stat = nf90_inq_varid(ncid, timename, varid)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varid for time: ', nf90_strerror(stat))
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, ndims=ndims)
      if (stat /= 0) then
        call document(subname, 'nf_inq_varndims for time: ', nf90_strerror(stat))
        go to 99
      end if
      if (ndims /= 1) then
        call document(subname, 'ndims /= 1 for time')
        stat = 1
        go to 99
      end if

      stat = nf90_inquire_variable(ncid, varid, dimids=time_dimids)
      if (stat /= 0) then
        call document(subname, 'nf_inq_vardimid for time: ', nf90_strerror(stat))
        go to 99
      end if

      if (time_dimids(1) /= data_dimids(data_rank)) then
        call document(subname, 'time dimid /= data_dimid(data_rank)')
        stat = 1
        go to 99
      end if

      ! read 'time' variable

      stat = nf90_get_var(ncid, varid, data_yr)
      if (stat /= 0) then
        call document(subname, 'nf_get_var for time: ', nf90_strerror(stat))
        go to 99
      end if

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

      data_yr(:) = data_yr(:) + ref_year
    end if

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                           &/ subname)

    !---------------------------------------------------------------------
    !  Broadcast the variables to other tasks beside master_task
    !---------------------------------------------------------------------

    call broadcast_array(data_dimlens(1:data_rank), master_task)

    if (my_task /= master_task) then
      select case (data_rank)
      case (1)
        allocate(data_yr(data_dimlens(1)))
        allocate(data_1d(data_dimlens(1)))
      case (2)
        allocate(data_yr(data_dimlens(2)))
        allocate(data_2d(data_dimlens(1), data_dimlens(2)))
      end select
    end if

    call broadcast_array(data_yr, master_task)

    select case (data_rank)
    case (1)
      data_nbval = data_dimlens(1)
      call broadcast_array(data_1d, master_task)
    case (2)
      data_nbval = data_dimlens(2)
      call broadcast_array(data_2d, master_task)
    end select

  end subroutine forcing_timeseries_read_file

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
    dpos = scan(ref_date(pos:ref_date_len), '-')

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

end module forcing_timeseries_mod
