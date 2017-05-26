!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module var_consistency_mod

!BOP
! !MODULE: var_consistency_mod
!
! !DESCRIPTION:
!  This module contains tools to verify that variables are consistent across all POP tasks.

! !REVISION HISTORY:
!     SVN:$Id$

! !USES:

   use kinds_mod,     only: log_kind, int_kind, r8
   use communicate,   only: my_task, master_task
   use broadcast,     only: broadcast_scalar, broadcast_array
   use io_tools,      only: stdout
   use POP_CommMod,   only: POP_Barrier
   use exit_mod,      only: exit_POP, sigAbort

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:
   public :: var_consistency_check

!EOP
!BOC

!-----------------------------------------------------------------------
!  generic interface definitions
!-----------------------------------------------------------------------

   interface var_consistency_check
      module procedure var_consistency_check_0D_int, &
                       var_consistency_check_1D_int, &
                       var_consistency_check_0D_log, &
                       var_consistency_check_1D_log, &
                       var_consistency_check_0D_r8,  &
                       var_consistency_check_1D_r8
   end interface

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_0D_int
! !INTERFACE:

   subroutine var_consistency_check_0D_int(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*),      intent(in) :: message
   integer (int_kind), intent(in) :: var

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   integer (int_kind) :: var_master ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master = var
   call broadcast_scalar(var_master, master_task)

   mismatch = .false.

   if (var /= var_master) then
     mismatch = .true.
     write(stdout, '(a,i0,1x,a,1x,2(1x,i0))') 'var_consistency_check mismatch: my_task=', &
       my_task, trim(message), var, var_master
   end if

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_0D_int

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_1D_int
! !INTERFACE:

   subroutine var_consistency_check_1D_int(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*),      intent(in) :: message
   integer (int_kind), intent(in) :: var(:)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   integer (int_kind) :: n
   integer (int_kind) :: var_master(size(var)) ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master(:) = var(:)
   call broadcast_array(var_master(:), master_task)

   mismatch = .false.

   do n = 1, size(var)
     if (var(n) /= var_master(n)) then
       mismatch = .true.
       write(stdout, '(a,i0,1x,a,a,i0,1x,2(1x,i0))') 'var_consistency_check mismatch: my_task=', &
         my_task, trim(message), ', n=', n, var(n), var_master(n)
     end if
   end do

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_1D_int

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_0D_log
! !INTERFACE:

   subroutine var_consistency_check_0D_log(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*),      intent(in) :: message
   logical (log_kind), intent(in) :: var

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   logical (log_kind) :: var_master ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master = var
   call broadcast_scalar(var_master, master_task)

   mismatch = .false.

   if (var .neqv. var_master) then
     mismatch = .true.
     write(stdout, '(a,i0,1x,a,1x,2(1x,l1))') 'var_consistency_check mismatch: my_task=', &
       my_task, trim(message), var, var_master
   end if

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_0D_log

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_1D_log
! !INTERFACE:

   subroutine var_consistency_check_1D_log(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*),      intent(in) :: message
   logical (log_kind), intent(in) :: var(:)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   integer (int_kind) :: n
   logical (log_kind) :: var_master(size(var)) ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master(:) = var(:)
   call broadcast_array(var_master(:), master_task)

   mismatch = .false.

   do n = 1, size(var)
     if (var(n) .neqv. var_master(n)) then
       mismatch = .true.
       write(stdout, '(a,i0,1x,a,a,i0,1x,2(1x,l1))') 'var_consistency_check mismatch: my_task=', &
         my_task, trim(message), ', n=', n, var(n), var_master(n)
     end if
   end do

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_1D_log

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_0D_r8
! !INTERFACE:

   subroutine var_consistency_check_0D_r8(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: message
   real (r8),     intent(in) :: var

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   real (r8)          :: var_master ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master = var
   call broadcast_scalar(var_master, master_task)

   mismatch = .false.

   if (var /= var_master) then
     mismatch = .true.
     write(stdout, '(a,i0,1x,a,1x,2(1x,1pe23.16))') 'var_consistency_check mismatch: my_task=', &
       my_task, trim(message), var, var_master
   end if

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_0D_r8

!***********************************************************************
!BOP
! !IROUTINE: var_consistency_check_1D_r8
! !INTERFACE:

   subroutine var_consistency_check_1D_r8(message, var)

! !DESCRIPTION:
!  Verify that var is consistent across all tasks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: message
   real (r8),     intent(in) :: var(:)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: mismatch
   integer (int_kind) :: n
   real (r8)          :: var_master(size(var)) ! var, broadcast from master_task

!-----------------------------------------------------------------------

   ! broadcast var from master_task to all tasks, for comparison
   if (my_task == master_task) var_master(:) = var(:)
   call broadcast_array(var_master(:), master_task)

   mismatch = .false.

   do n = 1, size(var)
     if (var(n) /= var_master(n)) then
       mismatch = .true.
       write(stdout, '(a,i0,1x,a,a,i0,1x,2(1x,1pe23.16))') 'var_consistency_check mismatch: my_task=', &
         my_task, trim(message), ', n=', n, var(n), var_master(n)
     end if
   end do

   call POP_Barrier()

   if (mismatch) call exit_POP(sigAbort, 'ERROR: var_consistency_check mismatch')

!-----------------------------------------------------------------------

   end subroutine var_consistency_check_1D_r8

!***********************************************************************

 end module var_consistency_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
