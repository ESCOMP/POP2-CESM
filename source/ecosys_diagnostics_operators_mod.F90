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
    character(len=char_len) :: single_line, diag_name, diag_ops, err_msg
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
        !     - get_diag_name_and_operators() will abort if the line is formatted incorrectly
        !       (at this point, that just means "there is no ':' separating diag_name from operators")
        call get_diag_name_and_operators(single_line, diag_name, diag_ops)

        ! (b) If diag_name is not a valid diagnostic name, abort
        surface_diag_ind = get_diag_ind(diag_name, surface_diags)
        if (surface_diag_ind .eq. 0) then
          interior_diag_ind = get_diag_ind(diag_name, interior_diags)
          if (interior_diag_ind .eq. 0) then
            write(err_msg, "(3A)") "Can not find ", trim(diag_name), " in list of diagnostics from MARBL"
            call exit_POP(sigAbort, err_msg)
          end if
        else
          interior_diag_ind = 0
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
        !     - parse_diag_ops() will abort if the operator count exceeds max_marbl_streams
        if (surface_diag_ind .ne. 0) then
          call parse_diag_ops(diag_ops, marbl_diags_surface_operators(surface_diag_ind,:), &
                              marbl_diags_surface_stream_cnt(surface_diag_ind))
        else
          call parse_diag_ops(diag_ops, marbl_diags_interior_operators(interior_diag_ind,:), &
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
    ! (Return 0 if no match is found)
    use marbl_interface_public_types, only : marbl_diagnostics_type

    character(len=*), intent(in) :: diag_name
    type(marbl_diagnostics_type), intent(in) :: marbl_diags
    integer :: get_diag_ind

    integer :: n

    get_diag_ind = 0
    do n=1,size(marbl_diags%diags)
      if (trim(diag_name) .eq. trim(marbl_diags%diags(n)%short_name)) then
        get_diag_ind = n
        return
      end if
    end do

  end function get_diag_ind

  !-----------------------------------------------------------------------

  subroutine parse_diag_ops(diag_ops, diag_ops_array, num_ops)
    ! Split "operator1, operator2, ..., operatorN" into the N elements of a string array

    character(len=*),  intent(in)  :: diag_ops
    integer(int_kind), intent(out) :: diag_ops_array(:)
    integer(int_kind), intent(out) :: num_ops

    character(len=char_len) :: err_msg
    integer :: str_pos, comma_pos, str_len

    num_ops = 0
    str_pos = 1
    str_len = len(diag_ops)
    do while (len_trim(diag_ops(str_pos:str_len)) .gt. 0)
      ! Update operator count
      num_ops = num_ops + 1
      if (num_ops .gt. max_marbl_streams) then
        write(err_msg,"(3A,I0,A,I0)") "ERROR: '", trim(diag_ops), "' contains at least ", num_ops, &
                                      " operators but memory only provided for ", size(diag_ops_array)
        call exit_POP(sigAbort, err_msg)
      end if

      ! Find first comma after str_pos
      comma_pos = index(diag_ops(str_pos:str_len), ',')
      if (comma_pos .eq. 0) comma_pos = str_len - str_pos + 2 ! imaginary comma after last character
      diag_ops_array(num_ops) = get_tavg_method(adjustl(diag_ops(str_pos:str_pos+comma_pos-2)))
      str_pos = str_pos + comma_pos
    end do

  end subroutine parse_diag_ops

  !-----------------------------------------------------------------------

  subroutine get_diag_name_and_operators(line_in, diag_name, diag_ops)
    ! Split a "DIAGNOSTIC_NAME : operator1, operator2, ..., operatorN" into
    ! strings containing "DIAGNOSTIC_NAME" and "operator1, operator2, ..., operatorN"

    character(len=*), intent(in)  :: line_in
    character(len=*), intent(out) :: diag_name, diag_ops

    character(len=char_len) :: err_msg
    integer :: n

    n = index(line_in, ':')
    if (n .ne. 0) then
      diag_name = adjustl(line_in(1:n-1))
      diag_ops  = adjustl(line_in(n+1:len(line_in)))
      ! Make sure neither string is empty
      if (len_trim(diag_name) .eq. 0) then
        write(err_msg, "(3A)") "ERROR: the line '", trim(line_in), "' does not contain a diagnostic"
        call exit_POP(sigAbort, err_msg)
      end if
      if (len_trim(diag_ops) .eq. 0) then
        write(err_msg, "(3A)") "ERROR: the line '", trim(line_in), "' does not contain an operator"
        call exit_POP(sigAbort, err_msg)
      end if
    else
      write(err_msg,"(3A)") "ERROR: could not parse the line '", trim(line_in), "'"
      call exit_POP(sigAbort, err_msg)
    end if

  end subroutine get_diag_name_and_operators

  !-----------------------------------------------------------------------

  function get_tavg_method(diag_ops)
    ! Convert a string (such as 'average' to tavg_method integer)

    use tavg, only : tavg_method_avg
    use tavg, only : tavg_method_min
    use tavg, only : tavg_method_max
    use tavg, only : tavg_method_constant

    character(len=*), intent(in) :: diag_ops
    integer(int_kind) :: get_tavg_method

    select case(trim(diag_ops))
      case ('average')
        get_tavg_method = tavg_method_avg
      case ('minimum')
        get_tavg_method = tavg_method_min
      case ('maximum')
        get_tavg_method = tavg_method_max
      case ('instantaneous')
        get_tavg_method = tavg_method_constant
      case DEFAULT
        get_tavg_method = tavg_method_unknown
    end select
  end function get_tavg_method

  !-----------------------------------------------------------------------

end module ecosys_diagnostics_operators_mod