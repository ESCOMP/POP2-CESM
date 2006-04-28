!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  module tracer_types

!BOP
! !MODULE: tracer_types
!
! !DESCRIPTION:
!  This module contains f90 derived type information for passive tracers
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:
 
  use kinds_mod
      
  implicit none
  save
 
!EOP
!BOC
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

