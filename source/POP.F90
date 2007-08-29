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

#ifdef SINGLE_EXEC
   use MPH_module, only : MPH_get_argument
#endif
   use POP_KindsMod
   use POP_InitMod, only: POP_Initialize, fstop_now, nscan, timer_total
   use POP_FinalMod
   use kinds_mod, only: int_kind, r8
   use communicate, only: my_task, master_task
   use exit_mod
   use timers, only: timer_print_all, get_timer, timer_start, timer_stop
   use time_management, only: init_time_flag, check_time_flag, sigAbort, &
       nsteps_run, stdout, sigExit, exit_pop, set_time_flag
   use step_mod, only: step
   use diagnostics, only: check_KE
   use output, only: output_driver
   use solvers, only: solv_sum_iters
   use registry
#if coupled
   use forcing_coupled, only: pop_coupling
#endif

   implicit none

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode         ! error code

#ifdef SINGLE_EXEC
   integer (int_kind) :: &
      nThreads

   call MPH_get_argument("THREADS", nThreads, "ocn")
#ifdef _OPENMP
   call OMP_SET_NUM_THREADS(nThreads)
#endif
#endif

!-----------------------------------------------------------------------
!
!  initialize the model run 
!
!-----------------------------------------------------------------------

   call POP_Initialize(errorCode)

!-----------------------------------------------------------------------
!
!  start up the main timer
!
!-----------------------------------------------------------------------

   call timer_start(timer_total)

!-----------------------------------------------------------------------
!
!  advance the model in time
!
!-----------------------------------------------------------------------

   advance: do while (.not. check_time_flag(fstop_now))

      call pop_coupling
      if ( registry_match('lcoupled') .and. check_time_flag(fstop_now)) exit advance

      call step

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
