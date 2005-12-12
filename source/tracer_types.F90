  module tracer_types
 
!-----------------------------------------------------------------------
!   This module contains f90 derived type information for passive tracers
!
!   CVS:$Id$
!   CVS:$Name$
!
!-----------------------------------------------------------------------
 
      use kinds_mod
      
      implicit none
      save
 
!-----------------------------------------------------------------------
!     define constants f90 types used for passive tracer tavg registration
!-----------------------------------------------------------------------

      integer (int_kind), parameter ::  &
         tavg_passive_interior_type = 1      &
      ,  tavg_passive_stf_type = 2


      type tavg_passive_nonstd
         character(char_len) :: sname
         integer (int_kind) :: grid, type, ndims
         real(r4), dimension(:,:,:), pointer :: data
      end type tavg_passive_nonstd
 
 
  end module tracer_types

