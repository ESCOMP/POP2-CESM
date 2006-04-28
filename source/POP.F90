!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !ROUTINE: POP
! !INTERFACE:

#ifdef SINGLE_EXEC
 subroutine ccsm_ocn()
#else
 program POP
#endif

! !DESCRIPTION:
!  This is the main driver for the Parallel Ocean Program (POP).
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_InitMod
   use POP_FinalMod
   use kinds_mod, only: int_kind, r8
   use communicate, only: my_task, master_task
   use exit_mod
   use domain, only: distrb_clinic
   use timers, only: timer_print_all, get_timer, timer_start, timer_stop
   use time_management, only: init_time_flag, check_time_flag, sigAbort,    &
       nsteps_run, stdout, sigExit, exit_pop, set_time_flag
   use step_mod, only: step
   use initial, only: initialize_pop
   use diagnostics, only: check_KE
   use output, only: output_driver
   use solvers, only: solv_sum_iters
   use forcing_coupled, only: lcoupled

   implicit none

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      timer_total,       &! timer number for total time
      timer_step,        &! timer number for step
      ierr,              &! error flag
      fstop_now,         &! flag id for stop_now flag
      nscan

   integer (POP_i4) :: &
      errorCode         ! error code

!-----------------------------------------------------------------------
!
!  initialize the model run
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_Initialize(errorCode)

   fstop_now = init_time_flag('stop_now')
   nscan = 0

!-----------------------------------------------------------------------
!
!  start up the main timer
!
!-----------------------------------------------------------------------

   call get_timer(timer_total,'TOTAL',1,distrb_clinic%nprocs)
   call timer_start(timer_total)

   call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!
!  advance the model in time
!
!-----------------------------------------------------------------------

   advance: do while (.not. check_time_flag(fstop_now))

      call timer_start(timer_step)
      call step
      call timer_stop(timer_step)

      if (lcoupled .and. check_time_flag(fstop_now)) exit advance

      nscan = nscan + solv_sum_iters

      !***
      !*** exit if energy is blowing
      !***

      if (check_KE(100.0_r8)) then
         call set_time_flag(fstop_now,.true.)
         call output_driver
         call exit_POP(sigAbort,'ERROR: k.e. > 100 ')
      endif

!-----------------------------------------------------------------------
!
!     write restart dumps and archiving
!
!-----------------------------------------------------------------------

      call output_driver

   enddo advance

!-----------------------------------------------------------------------
!
!  write an end restart if we are through the stepping loop 
!  without an error
!
!-----------------------------------------------------------------------

   nscan = nscan/nsteps_run
   if (my_task == master_task) & 
      write(stdout,*) ' average # scans =', nscan

!-----------------------------------------------------------------------
!
!  print timing information and clean up various environments if 
!  they have been used
!
!-----------------------------------------------------------------------

   call timer_stop(timer_total)

   call POP_Final(errorCode)

!-----------------------------------------------------------------------
!EOC

#ifdef SINGLE_EXEC
 end subroutine ccsm_ocn
#else
 end program POP
#endif

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
