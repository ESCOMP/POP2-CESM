!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module utils_mod

  !BOP
  ! !MODULE: utils_mod

  ! !DESCRIPTION:
  !  This module is meant to be a home to a variety of utilities that
  !  are software-focused rather than science-focused, such as a tool
  !  to split a string into an array of strings whenever a specified
  !  separator is encountered (similar to str.split() in python)
  !
  !  Written by: Michael Levy, NCAR, Dec 2017

  ! !REVISION HISTORY:
  !  SVN:$Id$

  use kinds_mod, only : char_len
  use exit_mod,  only : sigAbort, exit_POP

  implicit none
  public

contains

  !-----------------------------------------------------------------------

  subroutine utils_split(str_in, separator, array_out)
    ! Return an array of substrings of str_in, where separator is the delimiter
    ! between the substrings and both leading / trailing white space is removed.
    ! E.g. after
    !
    !   call utils_split('foo, bar, baz ', ',', array_out)
    !
    ! array_out = (/'foo', 'bar', 'baz')
    use shr_string_mod, only : shr_string_countChar

    character(len=*),               intent(in)  :: str_in
    character,                      intent(in)  :: separator
    character(len=*), dimension(:), intent(out) :: array_out

    ! local variables
    character(len=*), parameter :: subname = 'utils_mod:utils_split'
    character(len=char_len) :: log_message
    integer :: substr_num, substr_start, substr_ind

    ! abort if array_out is not long enough to contain all substrings
    if (size(array_out) .lt. shr_string_countChar(str_in, separator)+1) then
      write(log_message, "(A,I0,A,I0)") "There are ", shr_string_countChar(str_in, separator)+1, &
            " substrings, but array_out is size ", size(array_out)
      call exit_POP(sigAbort, log_message)
    end if

    ! initialize substring counter and position tracker
    array_out = ''
    substr_num = 1
    substr_start = 1

    ! loop through str_in character by character
    do substr_ind = 1, len_trim(str_in)
      ! When separator is encountered, copy previous string to array_out(substr_num)
      if (str_in(substr_ind:substr_ind) .eq. separator) then
        array_out(substr_num) = adjustl(str_in(substr_start:substr_ind-1))
        ! next substring will go in next index of array_out
        substr_num = substr_num + 1
        ! next substring starts with next character
        substr_start = substr_ind + 1
      ! Also make sure to copy the last substring over
      elseif (substr_ind .eq. len_trim(str_in)) then
        array_out(substr_num) = adjustl(str_in(substr_start:substr_ind-1))
      end if
    end do

  end subroutine utils_split

end module utils_mod