!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tracer_share

!BOP
! !MODULE: tracer_share

! !DESCRIPTION:
!  This module includes subroutines and functions that are used by 
!  several modules, to avoid code duplication.
!  Subroutines and Functions in this module can be called by 
!  individual tracer modules.

! !REVISION HISTORY:
!  SVN:$Id:$

! !USES:

   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use constants, only: c0
 
   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public ::                                 &
      SCHMIDT_CO2

!EOP
!BOC


!EOC
!***********************************************************************

 contains

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_CO2
! !INTERFACE:

 function SCHMIDT_CO2(SST, LAND_MASK)

! !DESCRIPTION:
!  Compute Schmidt number of CO2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
!  pp. 7373-7382, May 15, 1992
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(in) :: SST

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: SCHMIDT_CO2

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a = 2073.1_r8, &
      b = 125.62_r8, &
      c = 3.6276_r8, &
      d = 0.043219_r8

   where (LAND_MASK)
      SCHMIDT_CO2 = a + SST * (-b + SST * (c + SST * (-d)))
   elsewhere
      SCHMIDT_CO2 = c0
   end where

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_CO2

!*****************************************************************************
 

 end module tracer_share

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
