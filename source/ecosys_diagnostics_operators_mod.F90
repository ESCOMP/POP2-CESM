!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_diagnostics_operators_mod

  !BOP
  ! !MODULE: ecosys_tavg

  ! !DESCRIPTION:
  !  This module is used to read the marbl_diagnostics_operators file
  !  and provide information to ecosys_tavg
  !
  !  Written by: Michael Levy, NCAR, Dec 2017


  ! !REVISION HISTORY:
  !  SVN:$Id$

  ! !USES:
  use kinds_mod, only : int_kind
  use kinds_mod, only : char_len
  use tavg,      only : tavg_method_unknown
  use exit_mod,  only : sigAbort, exit_POP

  implicit none
  private

  !-----------------------------------------------------------------------
  !  Variables used to help define / accumulate tavg_fields in ecosys_tavg
  !-----------------------------------------------------------------------

  integer(int_kind),       parameter   :: max_marbl_streams = 2
  integer(int_kind),       allocatable :: marbl_diags_surface_stream_cnt(:), marbl_diags_interior_stream_cnt(:)
  integer(int_kind),       allocatable :: marbl_diags_surface_operators(:,:), marbl_diags_interior_operators(:,:)

  public :: ecosys_diagnostics_operators_init
  public :: max_marbl_streams
  public :: marbl_diags_surface_stream_cnt
  public :: marbl_diags_interior_stream_cnt
  public :: marbl_diags_surface_operators
  public :: marbl_diags_interior_operators

contains

  !-----------------------------------------------------------------------

  subroutine ecosys_diagnostics_operators_init(file_in, surface_diags, interior_diags)
    ! Initializes module variables based on contents of marbl_diagnostics_operators file
    ! (Indexing in marbl_diags_surface_* and marbl_diags_interior_* should match the
    ! indexing in surface_diags and interior_diags, respectively)
    use io_types,                     only : nml_in
    use communicate,                  only : my_task, master_task
    use broadcast,                    only : broadcast_scalar, broadcast_array
    use marbl_interface_public_types, only : marbl_diagnostics_type


    character(len=*),             intent(in) :: file_in
    type(marbl_diagnostics_type), intent(in) :: surface_diags
    type(marbl_diagnostics_type), intent(in) :: interior_diags

    ! Local variables
    character(len=char_len) :: single_line, single_diag, single_op, err_msg
    integer :: io_err, diag_cnt, diag_cnt_in
    integer :: surface_diag_ind, interior_diag_ind

    ! Set up surface diag variables
    diag_cnt_in = size(surface_diags%diags)
    allocate(marbl_diags_surface_stream_cnt(diag_cnt_in))
    allocate(marbl_diags_surface_operators(diag_cnt_in, max_marbl_streams))

    ! Set up interior diag variables
    diag_cnt_in = size(interior_diags%diags)
    allocate(marbl_diags_interior_stream_cnt(diag_cnt_in))
    allocate(marbl_diags_interior_operators(diag_cnt_in, max_marbl_streams))

    ! Intialize module variables / arrays
    marbl_diags_surface_stream_cnt = 0
    marbl_diags_surface_operators = tavg_method_unknown
    marbl_diags_interior_stream_cnt = 0
    marbl_diags_interior_operators = tavg_method_unknown

    if (my_task .eq. master_task) then
      open(unit=nml_in, file=file_in, iostat=io_err)
      if (io_err .ne. 0) then
        call exit_POP(sigAbort, "Error opening marbl_diagnostics_operators file")
      end if
    end if
    ! Set io_err on non-master tasks as well
    io_err = 0

    ! Initialize local variables used for reading file / error-checking
    diag_cnt = 0
    single_line = ''

    ! Read file
    do while (io_err .eq. 0)
      ! (1) Broadcast line just read in on master_task to all tasks
      !     (Initial entry to this loop will broadcast blank line which will be ignored)
      call broadcast_scalar(single_line, master_task)

      ! (2) Remove leading spaces from line just read, and treat lines beginning with '#' as empty
      single_line = adjustl(single_line)
      if (single_line(1:1) .eq. '#') then
        single_line = ''
      end if

      ! (3) process non-empty lines
      if (len_trim(single_line) .gt. 0) then
        ! (a) get the diagnostic name and (all) operators
        !     - get_varname_and_operators() will abort if the line is formatted incorrectly
        !       (at this point, that just means "there is no ':' separating diag_name from operators")
        call get_varname_and_operators(single_line, single_diag, single_op)

        ! (b) If single_diag is not a valid diagnostic name, abort
        surface_diag_ind = get_diag_ind(single_diag, surface_diags)
        if (surface_diag_ind .eq. -1) then
          interior_diag_ind = get_diag_ind(single_diag, interior_diags)
          if (interior_diag_ind .eq. -1) then
            write(err_msg, "(3A)") "Can not find ", trim(single_diag), " in list of diagnostics from MARBL"
            call exit_POP(sigAbort, err_msg)
          end if
        else
          interior_diag_ind = -1
        end if

        ! (c) Increase diag_cnt, make sure this does not exceed max allowable number
        !     - In POP, this max will be the total number of diagnostics returned from MARBL,
        !       so exceeding this cap will indicate a formatting error (diagnostics appearing
        !       multiple times or an unknown diagnostic being included in the file)
        diag_cnt = diag_cnt+1
        if (diag_cnt .gt. diag_cnt_in) then
          write(err_msg,"(A,I0,A)") "ERROR: read in line number ", diag_cnt, " but no memory to store it in array!"
          call exit_POP(sigAbort, err_msg)
        end if

        ! (d) Save the diagnostic name as well as all operators associated with it
        !     - parse_op() will abort if the operator count exceeds max_marbl_streams
        if (surface_diag_ind .ne. 0) then
          call parse_op(single_op, marbl_diags_surface_operators(surface_diag_ind,:), &
                        marbl_diags_surface_stream_cnt(surface_diag_ind))
        else
          call parse_op(single_op, marbl_diags_interior_operators(interior_diag_ind,:), &
                        marbl_diags_interior_stream_cnt(interior_diag_ind))
        end if
      end if

      ! (4) Read next line on master, iostat value out (that's how loop ends)
      if (my_task .eq. master_task) then
        read(nml_in, "(A)", iostat=io_err) single_line
      end if
      call broadcast_scalar(io_err, master_task)
    end do
    ! Abort if iostat did not return "End of File" status code

    ! Close the file on master task
    if (my_task .eq. master_task) close(nml_in)

  end subroutine ecosys_diagnostics_operators_init

  !-----------------------------------------------------------------------

  function get_diag_ind(diag_name, marbl_diags)
    ! Return index of diags such that diags%short_name == diag_name
    ! (Return -1 if no match is found)
    use marbl_interface_public_types, only : marbl_diagnostics_type

    character(len=*), intent(in) :: diag_name
    type(marbl_diagnostics_type), intent(in) :: marbl_diags
    integer :: get_diag_ind

    integer :: n

    get_diag_ind = -1
    do n=1,size(marbl_diags%diags)
      if (trim(diag_name) .eq. trim(marbl_diags%diags(n)%short_name)) then
        get_diag_ind = n
        return
      end if
    end do

  end function get_diag_ind

  !-----------------------------------------------------------------------

  subroutine parse_op(op_in, ops_out, op_cnt)
    ! Split "operator1, operator2, ..., operatorN" into the N elements of a string array

    character(len=*),  intent(in)  :: op_in
    integer(int_kind), intent(out) :: ops_out(:)
    integer(int_kind), intent(out) :: op_cnt

    character(len=char_len) :: err_msg
    integer :: n, last_comma

    op_cnt = 1
    last_comma = 0
    do n=1,len_trim(op_in)
      if (op_in(n:n) .eq. ',') then
        ops_out(op_cnt) = get_tavg_method(adjustl(op_in(last_comma+1:n-1)))
        last_comma = n
        op_cnt = op_cnt + 1
        if (op_cnt .gt. max_marbl_streams) then
          write(err_msg,"(3A,I0,A,I0)") "ERROR: '", trim(op_in), "' contains at least ", op_cnt, &
                                        " operators but memory only provided for ", size(ops_out)
          call exit_POP(sigAbort, err_msg)
        end if
      end if
    end do
    ops_out(op_cnt) = get_tavg_method(adjustl(op_in(last_comma+1:len_trim(op_in))))

  end subroutine parse_op

  !-----------------------------------------------------------------------

  subroutine get_varname_and_operators(line_in, var_out, op_out)
    ! Split a "DIAGNOSTIC_NAME : operator1, operator2, ..., operatorN" into
    ! strings containing "DIAGNOSTIC_NAME" and "operator1, operator2, ..., operatorN"

    character(len=*), intent(in)  :: line_in
    character(len=*), intent(out) :: var_out, op_out

    character(len=char_len) :: err_msg
    integer :: n

    do n=1,len_trim(line_in)
      if (line_in(n:n) .eq. ':') then
        var_out = adjustl(line_in(1:n-1))
        op_out = adjustl(line_in(n+1:len(line_in)))
        ! This return is only called if a ':' is found
        return
      end if
    end do
    ! Exiting the loop implies no ':' found
    write(err_msg,"(3A)") "ERROR: could not parse the line '", trim(line_in), "'"
    call exit_POP(sigAbort, err_msg)

  end subroutine get_varname_and_operators

  !-----------------------------------------------------------------------

  function get_tavg_method(op_in)
    ! Convert a string (such as 'average' to tavg_method integer)

    use tavg, only : tavg_method_avg
    use tavg, only : tavg_method_min
    use tavg, only : tavg_method_max
    use tavg, only : tavg_method_constant

    character(len=*), intent(in) :: op_in
    integer(int_kind) :: get_tavg_method

    select case(trim(op_in))
      case ('average')
        get_tavg_method = tavg_method_avg
      case ('minimum')
        get_tavg_method = tavg_method_min
      case ('maximum')
        get_tavg_method = tavg_method_max
      case DEFAULT
        get_tavg_method = tavg_method_unknown
    end select
  end function get_tavg_method

  !-----------------------------------------------------------------------

end module ecosys_diagnostics_operators_mod