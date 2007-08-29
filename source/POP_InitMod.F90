!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_InitMod

!BOP
! !MODULE: POP_InitMod
! !DESCRIPTION:
!  This module contains the POP initialization method and initializes
!  everything needed by a POP simulation.  Primarily it is a driver
!  that calls individual initialization routines for each POP module.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use kinds_mod, only: int_kind
   use initial
   use domain, only: distrb_clinic
   use timers, only: get_timer
   use time_management, only: init_time_flag
#if coupled
   use forcing_coupled, only: pop_init_coupled, pop_send_to_coupler
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_Initialize, POP_Initialize1, POP_Initialize2,  &
             POP_Initialize_coupling
             

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------
   integer (int_kind), public :: &
      fstop_now,                 &! flag id for stop_now flag
      timer_total,               &! timer for entire run phase
      nscan


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_Initialize 
! !INTERFACE:

 subroutine POP_Initialize(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines, with coupling
!  calls inbetween initialization calls. When invoking CCSM with "coupling
!  at the top," the routines called in this subroutine will be accessed 
!  elsewhere, and this routine will not be called.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  call pop initialization routines in two stages, with coupling inbetween
!
!-----------------------------------------------------------------------

   call POP_Initialize1(errorCode)

!-----------------------------------------------------------------------
!
!  exchange initial information with coupler
!
!-----------------------------------------------------------------------
   call POP_Initialize_coupling

!-----------------------------------------------------------------------
!
!  complete the initialiation of the model
!
!-----------------------------------------------------------------------

   call POP_Initialize2(errorCode)


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize 
!***********************************************************************
!BOP
! !IROUTINE: POP_Initialize1
! !INTERFACE:

 subroutine POP_Initialize1(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!  initialize return flag
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!  call pop initialization routines
!
!-----------------------------------------------------------------------

   call pop_init_phase1

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize1

!***********************************************************************

!BOP
! !IROUTINE: POP_Initialize2
! !INTERFACE:

 subroutine POP_Initialize2(errorCode)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
      errorCode              ! Returns an error code if any init fails

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  initialize return flag
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
!  complete pop initialization process
!
!-----------------------------------------------------------------------

   call pop_init_phase2

!-----------------------------------------------------------------------
!
!  initialize variables used by the pop driver
!
!-----------------------------------------------------------------------
   nscan = 0

!-----------------------------------------------------------------------
!
!  initialize driver-level flags and timers
!
!-----------------------------------------------------------------------
   fstop_now  = init_time_flag('stop_now')

   call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize2

!***********************************************************************

!BOP
! !IROUTINE: POP_Initialize_coupling
! !INTERFACE:

 subroutine POP_Initialize_coupling

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a POP run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module
! !USES


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


#if coupled
!-----------------------------------------------------------------------
!
!  exchange initial information with coupler
!
!-----------------------------------------------------------------------

   call pop_init_coupled

!-----------------------------------------------------------------------
!
!  send initial state information to the coupler
!
!-----------------------------------------------------------------------
   call pop_send_to_coupler

#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Initialize_coupling

 end module POP_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
