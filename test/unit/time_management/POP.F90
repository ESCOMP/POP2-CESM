!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_time_managementTest

!----------------------------------------------------------------------
!
!  this program tests the POP time_management.F90 module over
!  a standard range of options.
!
!----------------------------------------------------------------------

   use kinds_mod
   use io
   use io_tools
   use time_management

   implicit none


!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      stop_now,          &! flag id for stop_now flag
      thats_enough,      &! total number of elapsed days
      iyear_end,         &
      imonth_end,        &
      iday_end

   integer (i4) :: iostat,ierr,nml_error

   logical IsOpen

   namelist /driver_nml/ iyear_end, imonth_end, iday_end


!----------------------------------------------------------------------
!
!  Read driver namelist
!
!----------------------------------------------------------------------
   

   if (my_task == master_task) then

     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     write(6,*) ' nml_error = ', nml_error

      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif

      do while (nml_error > 0)
         read(nml_in, nml=driver_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)

   endif ! master_task

!----------------------------------------------------------------------
!
!  initialize time management
!
!----------------------------------------------------------------------
   call register_string('init_ts')

   call init_time1
   call init_time2

   call init_time_flag('stop_now', stop_now) 
   call ymd2eday (iyear_end, imonth_end, iday_end, thats_enough)

!-----------------------------------------------------------------------
!
!  test the advancement of the model in time
!
!-----------------------------------------------------------------------

   stepping_loop: do while (.not. check_time_flag(stop_now) )

      call time_manager(.true., .true.)
 
      call report

     !if (elapsed_days_this_run > thats_enough)  exit stepping_loop

   enddo stepping_loop


!----------------------------------------------------------------------

 end program POP_time_managementTest

 subroutine report

   use time_management

     if (eoy ) then
        write(6,*) '=================================================='
     endif

     if (eoy .or. eom .or. eod) write(6,*) ' '
      
    !write(6,1100) iyear, imonth, iday, seconds_this_day, eoy, eom, eod

 1100 format (1x, i4.4,'-',i2.2,'-',i2.2, 2x, f8.2, 3L3)

 end subroutine report

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
