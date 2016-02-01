module marbl_kinds_mod

  implicit none

  public

  integer, parameter, public ::                 &
       char_len       = 384                    ,&
       char_len_long  = 512                    ,&
       log_kind       = kind(.true.)           ,&
       int_kind       = kind(1)                ,&
       i4             = selected_int_kind(6)   ,&
       i8             = selected_int_kind(13)  ,&
       r8             = selected_real_kind(13) ,&  
       r4             = selected_real_kind(6)  

  ! FIXME(bja, 2015-02) these don't belong in the kinds module, but
  ! not sure where they should go. Only used internally to marbl? Used
  ! as part of the API interface? Maybe they do belong in the kinds
  ! module because they are a specific kind of c0 and c1....
  real(kind=r8), parameter, public :: &
      c0     =    0.0_r8       , &
      c1     =    1.0_r8
  
end module marbl_kinds_mod
