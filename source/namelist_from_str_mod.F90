module namelist_from_str_mod

  ! This module contains functions to convert from one long string to either
  ! an array of strings each containing a single namelist or a single variable.
  ! There is also a function to return a requested namelist (from the former),
  ! which is used for the read() calls. This is needed due to what I believe is
  ! a bug in gfortran, where reading a string that does not contain the requested
  ! namelist returns a successful status code so we can't just loop over all
  ! elements of nl_buffer until the namelist has been read.

  implicit none
  private

  ! Need to know what carriage return is on the system; use #define if we
  ! come across a machine that doesn't use achar(10)

  character,       parameter :: cr = achar(10)

  public :: namelist_split_by_line
  public :: namelist_split_by_nl
  public :: namelist_find

  !***********************************************************************

contains

  !***********************************************************************

  subroutine namelist_split_by_line(str_in, array_out)
  ! This routine takes a string (str_in) containing the entire contents of a
  ! a namelist file and returns an array of strings (array_out) where each
  ! element contains a line from the file.

    ! FIXME: This routine depends on the namelist file conforming
    !        to very specific formatting - a more general / robust
    !        solution would be preferred

    ! FIXME: Strip comments out of str_in (without accidentally removing
    !        strings that happen to contain exclamation points)

    character(len=*), intent(in) :: str_in
    ! array_out is intent(inout) because we initialized to '' previously
    ! (and also to save memory)
    character(len=*), dimension(:), intent(inout) :: array_out

    character(len=len(str_in)) :: str_tmp
    integer :: old_pos, line_cnt, i, j

    ! each line needs to be stored in different element of array_out
    old_pos = 1
    line_cnt = 1
    do i=1,len_trim(str_in)
      if ((str_in(i:i) .eq. cr) .or. (i.eq.len_trim(str_in))) then
        ! FIXME: add error checking in case 
        !        (i+1-old_pos) > nl_buffer_size
        array_out(line_cnt) = str_in(old_pos:i-1)
        line_cnt = line_cnt+1
        old_pos = i+1
      end if
    end do

  end subroutine namelist_split_by_line

  !***********************************************************************

  subroutine namelist_split_by_nl(str_in, array_out)
  ! This routine takes a string (str_in) containing the entire contents of a
  ! a namelist file and returns an array of strings (array_out) where each
  ! element contains a single namelist. It also removes all carriage returns
  ! from the elements of array_out

    ! FIXME: This routine depends on the namelist file conforming
    !        to very specific formatting - a more general / robust
    !        solution would be preferred

    ! FIXME: Strip comments out of str_in (without accidentally removing
    !        strings that happen to contain exclamation points)

    character(len=*), intent(in) :: str_in
    ! array_out is intent(inout) because we initialized to '' previously
    ! (and also to save memory)
    character(len=*), dimension(:), intent(inout) :: array_out

    character(len=len(str_in)) :: str_tmp
    integer :: old_pos, nl_cnt, i, j

    ! each namelist needs to be stored in different element of array_out
    old_pos = 1
    nl_cnt = 1
    do i=1,len_trim(str_in)-1
      if (str_in(i:i+1) .eq. '/' // cr) then
        ! FIXME: add error checking in case 
        !        (i+1-old_pos) > nl_buffer_size
        array_out(nl_cnt) = str_in(old_pos:i)
        nl_cnt = nl_cnt+1
        old_pos = i+2
      end if
    end do

    ! We need to strip carriage returns from the namelist, replacing them with
    ! empty space
    do j= 1,nl_cnt
      str_tmp = array_out(j)
      do i=1,len_trim(str_tmp)
        if (str_tmp(i:i).eq.cr) then
          str_tmp(i:i) = ' '
        end if
      end do
      ! Remove whitespace from beginning of string (if any)
      array_out(j) = trim(adjustl(str_tmp))
      ! FIXME: add error checking in case first character is not '&'
    end do
  end subroutine namelist_split_by_nl

  !*****************************************************************************

  function namelist_find(nl_buffer, nl_name)

    character(len=*), intent(in) :: nl_buffer(:)
    character(len=*), intent(in) :: nl_name
    character(len=len(nl_buffer)) :: namelist_find

    character(len=len(nl_buffer)) :: single_namelist
    integer :: j, n

    ! Will return empty string if namelist not found
    namelist_find = ''

    ! Look for correct namelist in array
    do j = 1, size(nl_buffer)
      single_namelist = nl_buffer(j)
      n = len_trim(nl_name)
      if (single_namelist(2:n+1).eq.trim(nl_name)) then
         namelist_find = single_namelist
         exit
      end if
    end do

!    FIXME: call exit_pop if namelist is not found
!    if (trim(namelist_find).eq.'') then
!      write(log_message, "(2A)") trim(nl_name), ' is not included in nl_buffer'
!      call marbl_status_log%log_error(log_message, subname)
!      return
!    end if

  end function namelist_find

  !*****************************************************************************

end module namelist_from_str_mod
