 module io_tools
!===============
 
!-----------------------------------------------------------------------
!   This module contains routines to facilitate io
!
!   CVS:$Id$
!   CVS:$Name$
!
!-----------------------------------------------------------------------
 
   use kinds_mod
   use grid
   use io    
   use shr_sys_mod
      
   implicit none
   save
 
!-----------------------------------------------------------------------
!     interfaces
!-----------------------------------------------------------------------
 
   interface document
       module procedure document_char, &
                        document_int,  &
                        document_log,  &
                        document_dbl,  &
                        document_real 
   end interface
 
 contains
 
   subroutine document_char (sub_name,message)
!  ========================
 
   character (*) :: sub_name, message
 
   if (my_task == master_task) then
      write(stdout,1000)  sub_name, message
      call shr_sys_flush (stdout)
   endif
 
1000  format(/,5x,'(',a,')  ', a ) 
 
   end subroutine document_char
 
 
   subroutine document_int (sub_name,message,ival)
!  =======================
 
   character (*)       :: sub_name, message
   integer (int_kind)  :: ival
 
   if (my_task == master_task) then
      write(stdout,1000)  sub_name, message, ival
      call shr_sys_flush (stdout)
   endif
 
1000  format(/,5x,'(',a,')  ', a,1x,i10)
 
   end subroutine document_int
 
 
   subroutine document_log (sub_name,message,lval)
!  =======================
 
   character (*)      :: sub_name, message
   logical (log_kind) :: lval
 
   if (my_task == master_task) then
      write(stdout,1000)  sub_name, message, lval
      call shr_sys_flush (stdout)
   endif
 
1000  format(/,5x,'(',a,')  ', a,1x,L3)
 
   end subroutine document_log
 
   subroutine document_dbl (sub_name,message,dval)
!  =======================
 
   character (*) :: sub_name, message
   real (r8)     :: dval
 
   if (my_task == master_task) then
      write(stdout,1000)  sub_name, message, dval
      call shr_sys_flush (stdout)
   endif
 
1000  format(/,5x,'(',a,')  ', a, 1x, 1pe15.5)
 
   end subroutine document_dbl
 
 
   subroutine document_real (sub_name,message,rval)
!  ========================
 
   character (*) :: sub_name, message
   real(r4)    :: rval
 
   if (my_task == master_task) then
      write(stdout,1000)  sub_name, message, rval
      call shr_sys_flush (stdout)
   endif
 
1000  format(/,5x,'(',a,')  ', a, 1x, 1pe15.5)
 
   end subroutine document_real
 
 
!===================
 end module io_tools
!===================

