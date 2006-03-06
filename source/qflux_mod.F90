!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


 module qflux_mod

!***********************************************************************

!-----------------------------------------------------------------------
!
!     This module computes running time-averages of qflux
!    
!
!     CVS:$Id$
!
!-----------------------------------------------------------------------

   use kinds_mod
   use constants
   use exit_mod
   use time_management
   use tavg
   use shr_sys_mod

   implicit none
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_qflux
 
   integer (int_kind) :: &
      tavg_QFLUX          ! tavg id for QFLUX 


!***********************************************************************

      contains

!***********************************************************************

   subroutine init_qflux

   character (char_len) :: string
 
   
   tavg_sum_qflux = c0

   string = 'Internal Ocean Heat Flux Due to Ice Formation; ' /&
             &/ 'heat of fusion > 0 or ice-melting potential < 0 '
 
   call define_tavg_field(tavg_QFLUX,'QFLUX',2,               &
                          long_name=trim(string),             &
                          tavg_method=tavg_method_qflux,      &
                          missing_value=undefined_nf_r4,      &
                          units='Watts/meter^2',              &
                          coordinates  ='TLONG TLAT time',    &
                          grid_loc ='2111')

   end subroutine init_qflux
 
 
!***********************************************************************

 end module qflux_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
