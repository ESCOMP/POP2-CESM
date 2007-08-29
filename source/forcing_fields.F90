!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module forcing_fields

!BOP
! !MODULE: forcing_fields

! !DESCRIPTION:
!  Contains the forcing fields necessary for supporting high-level coupling
!  These fields originally resided in module forcing

! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use kinds_mod
   use blocks,      only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic,nt


! !PUBLIC DATA MEMBERS:

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      public, target :: &
      SMF,  &!  surface momentum fluxes (wind stress)
      SMFT   !  surface momentum fluxes on T points if avail

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      public, target :: &
      STF,  &!  surface tracer fluxes
      TFW    ! tracer content in freshwater flux


   logical (log_kind), public :: &
      lsmft_avail   ! true if SMFT is an available field

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      public, target ::  &
      IFRAC,             &! ice fraction; not initialized in this routine
      U10_SQR,           &! 10m wind speed squared; not initialized in this routine
      ATM_PRESS,         &! atmospheric pressure forcing
      FW,FW_OLD           ! freshwater flux at T points (cm/s)
                          ! FW_OLD is at time n-1


!***********************************************************************

 end module forcing_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
