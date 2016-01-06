!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module software_eng_mod

!BOP
! !MODULE: software_eng_mod
!
! !DESCRIPTION:
!  Contains a public logical variable that can be read in via the se_nml
!  namelist to help test answer-changing tags on branches. Use as follows:
!    tag_001 -- existing tag before answer-changing tag; use to create
!               baselines_001 (same results regardless of value of
!               lchange_ans)
!    tag_002 -- branch tag where lchange_ans = .false. => bit-for-bit with
!               tag_001 and lchange_ans = .true. => new modifications; make
!               sure lchange_ans = .false. matches baselines_001 then create
!               baselines_002 using lchange_ans = .true.
!    tag_003 -- branch tag where lchange_ans = .false. fork has been entirely
!               removed; compare to baselines_002 (same results regardless of
!               value of lchange_ans)
!  When the branch is merged back to master, tag_002 will be merged first to
!  ensure branch is still bit-for-bit when lchange_ans is .false. and then new
!  baselines will be constructed with lchange_ans set to .true. After those
!  two steps, tag_003 can be merged to the master and compared to the new
!  baselines.
!
! !REVISION HISTORY:
! $Id$

! !USES:

  use kinds_mod, only : log_kind

  use io_types, only : stdout
  use io_types, only : nml_in
  use io_types, only : nml_filename

  use communicate, only : my_task
  use communicate, only : master_task

  use broadcast, only : broadcast_scalar

  use exit_mod, only : exit_POP
  use exit_mod, only : sigAbort

  implicit none
  private

  public :: init_software_eng

  logical(kind=log_kind), public :: lchange_ans

  Contains

    subroutine init_software_eng()

      integer :: nml_error
      namelist /se_nml/lchange_ans

      ! Default value for namelist variable
      lchange_ans = .false.

      ! Read namelist on master_task
      if (my_task.eq.master_task) then
        open(nml_in, file=nml_filename, status='old', iostat=nml_error)
        if (nml_error /= 0) then
           nml_error = -1
        else
           nml_error =  1
        endif
        do while (nml_error > 0)
          read(nml_in, nml=se_nml, iostat=nml_error)
        end do
        if (nml_error.eq.0) close(nml_in)
      end if

      call broadcast_scalar(nml_error, master_task)
      if (nml_error.ne.0) then
        call exit_POP(sigAbort, 'ERROR reading se_nml')
      end if

      ! Write namelist to log
      if (my_task.eq.master_task) then
        write(stdout,*) ' se_nml namelist settings:'
        write(stdout,se_nml)
      end if

      ! Broadcast namelist variables
      call broadcast_scalar(lchange_ans, master_task)

    end subroutine init_software_eng

end module software_eng_mod
