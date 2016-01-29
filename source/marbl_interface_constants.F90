module marbl_interface_constants
  ! constants that need to be shared between the marbl public api and marbl private data.

  ! This module should hav NO dependencies on other modules other than marbl_kinds?

  implicit none

  private

  
  ! FIXME(bja, 2015-01) nl_buffer_size shouldn't be a hard coded
  ! constant, but runtime configurable?! Just not sure what the best
  ! approach is at the moment....
  integer, public, parameter :: marbl_nl_buffer_size = 262144

  ! string sizes
  integer, public, parameter :: marbl_str_length = 1024

  ! status codes
  integer, public, parameter :: marbl_status_ok = 0
  integer, public, parameter :: marbl_status_could_not_read_namelist = 1
  
end module marbl_interface_constants
