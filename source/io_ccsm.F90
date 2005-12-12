!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_ccsm

!BOP
! !MODULE: io_ccsm
!
! !DESCRIPTION:
!  This module provides a kludge interface for writing time_bound 
!  to ccsm netCDF output files
!
! !REVISION HISTORY:
!  CVS:$Id$
!  CVS:$Name$

! !USES:

   use kinds_mod
   use blocks
   use communicate
   use broadcast
   use exit_mod
   use domain
   use constants
   use io_netcdf
   use io_binary
   use io_types

   implicit none
   public  ! to get io_types without having to explicitly use io_types
           ! module directly
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: data_set_ccsm

!EOP
!BOC
!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: data_set_ccsm
! !INTERFACE:

 subroutine data_set_ccsm (data_file, operation, time_dim, d2_dim,  &
                           time_bound_id, lower_time_bound, upper_time_bound)

! !DESCRIPTION:
!  This routine is kludge interface to define_time_bound_netcdf and 
!  write_time_bound_netcdf
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent (in)   :: operation

! !INPUT/OUTPUT PARAMETERS:

   integer (i4),intent (inout)        :: time_bound_id
   type (datafile),   intent (inout)  :: data_file
   type (io_dim),     intent (inout)  :: time_dim, d2_dim
   real (r8), intent(in), optional    ::  &
      lower_time_bound,                   &
      upper_time_bound

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: num_writes  ! place-holder until more than one
                                     ! time level written to a single file


!-----------------------------------------------------------------------
!
!  select operation to perform
!
!-----------------------------------------------------------------------


   if (data_file%data_format=='bin') then
         call exit_POP(sigAbort, &
             '(data_set_ccsm) ERROR: cannot call this routine with bin format')
   endif

   select case (trim(operation))

!-----------------------------------------------------------------------
!
!  define an io field
!
!-----------------------------------------------------------------------

   case ('define')

      call define_time_bound_netcdf(data_file, time_dim, d2_dim,time_bound_id)

!-----------------------------------------------------------------------
!
!  write an io field
!
!-----------------------------------------------------------------------

   case ('write')

      if (.not. present(lower_time_bound) .or. .not. present(upper_time_bound)) then
         call exit_POP(sigAbort, &
            '(data_set_ccsm) ERROR: must specify lower and upper time bounds')
      endif

      num_writes =  1  ! for now, only support one time value per output file
      call write_time_bound_netcdf(data_file, time_bound_id, num_writes, &
                                   lower_time_bound=lower_time_bound,    &
                                   upper_time_bound=upper_time_bound)

!-----------------------------------------------------------------------
!
!  unknown operation
!
!-----------------------------------------------------------------------

   case default

      if (my_task == master_task) &
         write(stdout,*) 'data_set_ccsm operation: ',trim(operation)
      call exit_POP(sigAbort,'data_set_ccsm: Unknown operation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine data_set_ccsm

!***********************************************************************


 end module io_ccsm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
