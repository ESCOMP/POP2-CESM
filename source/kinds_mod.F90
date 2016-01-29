!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module kinds_mod

!BOP
! !MODULE: kinds_mod
!
! !DESCRIPTION:
!  This module defines default numerical data types for all common data
!  types like integer, character, logical, real4 and real8.
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:
!  uses no other modules

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   ! Note: we want char_len to be 256 and char_len_long to be 512 but this is
   !       a quick work-around because some of variables containing the name of
   !       restart files, history files, and other output files are set to
   !       char_len. A later update should change these to char_len_long and
   !       then reset char_len to its original 256 value. (MNL; 10 Aug 2015)
   integer, parameter, public ::               &
      char_len       = 384                    ,&
      char_len_long  = 512                    ,&
      log_kind       = kind(.true.)           ,&
      int_kind       = kind(1)                ,&
      i4             = selected_int_kind(6)   ,&
      i8             = selected_int_kind(13)  ,&
      r4             = selected_real_kind(6)  ,&
      r8             = selected_real_kind(13)

   integer, parameter, public ::               &
#ifdef TAVG_R8
      rtavg          = r8 ! nonstandard r8 for debugging purposes only
#else
      rtavg          = r4 ! standard, single-precision
#endif


!EOP
!BOC
!EOC
!***********************************************************************

 end module kinds_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
