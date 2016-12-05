!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module qflux_mod

!BOP
! !MODULE: qflux_mod
!
! !DESCRIPTION:
!  This module supports time-averaged qflux computations
!
! !REVISION HISTORY:
! SVN:$Id$
!

! !USES:

   use kinds_mod
   use blocks
   use constants
   use exit_mod
   use time_management
   use tavg

   implicit none
   public
   save

!EOP
!BOC

   integer (int_kind) :: tavg_QFLUX       ! tavg id for QFLUX 


   real (r8) ::  &
      q0_mean_budget  ! <Q_0> Robert Filter budget diagnostic term

   
!EOC
!***********************************************************************

      contains

!***********************************************************************

   subroutine init_qflux

   character (char_len) :: string
 
   
!-----------------------------------------------------------------------
!
!  time-averaged field definition
!
!-----------------------------------------------------------------------
 
   string = 'Internal Ocean Heat Flux Due to Ice Formation; ' /&
             &/ 'heat of fusion > 0 or ice-melting potential < 0 '
 
   call define_tavg_field(tavg_QFLUX,'QFLUX',2,               &
                          long_name=trim(string),             &
                          tavg_method=tavg_method_qflux,      &
                          units='Watts/meter^2',              &
                          coordinates  ='TLONG TLAT time',    &
                          grid_loc ='2111')

   string = 'Internal Ocean Heat Flux Due to Ice Formation for Robert Filter diagnostics; ' /&
             &/ 'heat of fusion > 0 or ice-melting potential < 0 '
 

   end subroutine init_qflux
 
 
!***********************************************************************

 end module qflux_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
