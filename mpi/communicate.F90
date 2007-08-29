!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP

 module communicate

! !MODULE: communicate
! !DESCRIPTION:
!  This module contains the necessary routines and variables for
!  communicating between processors.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use kinds_mod
   use cpl_interface_mod, only : cpl_interface_init,cpl_interface_finalize
   use cpl_fields_mod, only : cpl_fields_ocnname

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public  :: init_communicate,          &
              exit_message_environment,  &
              abort_message_environment, &
              get_num_procs,             &
              create_communicator

! !PUBLIC DATA MEMBERS:

   integer (int_kind), public :: &
      MPI_COMM_OCN,             &! MPI communicator for ocn comms
      mpi_dbl,                  &! MPI type for dbl_kind
      my_task,                  &! MPI task number for this task
      master_task,              &! task number of master task
      cpl_task                   ! task number of coupler master task

   integer (int_kind), parameter, public :: &
      mpitag_bndy_2d        = 1,    &! MPI tags for various
      mpitag_bndy_3d        = 2000, &!
      mpitag_gs             = 1000   ! communication patterns

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_communicate
! !INTERFACE:

 subroutine init_communicate

! !DESCRIPTION:
!  This routine sets up MPI environment and defines ocean
!  communicator.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   include 'mpif.h'   ! MPI Fortran include file

   integer (int_kind) :: ierr  ! MPI error flag

!-----------------------------------------------------------------------
!
!  initiate mpi environment and create communicator for internal
!  ocean communications
!
!-----------------------------------------------------------------------

   call create_ocn_communicator

   master_task = 0
   call MPI_COMM_RANK  (MPI_COMM_OCN, my_task, ierr)

!-----------------------------------------------------------------------
!
!  On some 64-bit machines where real_kind and dbl_kind are
!  identical, the MPI implementation uses MPI_REAL for both.
!  In these cases, set MPI_DBL to MPI_REAL.
!
!-----------------------------------------------------------------------

   MPI_DBL = MPI_DOUBLE_PRECISION

!-----------------------------------------------------------------------
!EOC

 end subroutine init_communicate

!***********************************************************************
!BOP
! !IROUTINE: get_num_procs
! !INTERFACE:

 function get_num_procs()

! !DESCRIPTION:
!  This function returns the number of processor assigned to
!  MPI_COMM_OCN
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind) :: get_num_procs

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr

!-----------------------------------------------------------------------

   call MPI_COMM_SIZE(MPI_COMM_OCN, get_num_procs, ierr)

!-----------------------------------------------------------------------
!EOC

 end function get_num_procs

!***********************************************************************
!BOP
! !IROUTINE: exit_message_environment
! !INTERFACE:

 subroutine exit_message_environment(ierr)

! !DESCRIPTION:
!  This routine exits the message environment properly when model
!  stops.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'   ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr   ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

#ifdef coupled
   call cpl_interface_finalize(cpl_fields_ocnname)
#else
   call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine exit_message_environment

!***********************************************************************
!BOP
! !IROUTINE: abort_message_environment
! !INTERFACE:

 subroutine abort_message_environment(ierr)

! !DESCRIPTION:
!  This routine aborts the message environment when model stops.
!  It will attempt to abort the entire MPI COMM WORLD.
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'   ! MPI Fortran include file

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: ierr   ! MPI error flag

!EOP
!BOC
!-----------------------------------------------------------------------

#ifdef coupled
   call MPI_BARRIER(MPI_COMM_OCN,ierr)
   ierr = 13
   call MPI_ABORT(0,ierr)
   call cpl_interface_finalize(cpl_fields_ocnname)
#else
   call MPI_BARRIER(MPI_COMM_OCN, ierr)
   call MPI_ABORT(MPI_COMM_WORLD, ierr)
   call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine abort_message_environment

!***********************************************************************
!BOP
! !IROUTINE: create_ocn_communicator
! !INTERFACE:

 subroutine create_ocn_communicator

! !DESCRIPTION:
!  This routine queries all the tasks in MPI_COMM_WORLD to see
!  which belong to the ocean.  In standalone mode, this should
!  be all tasks, but in coupled mode POP needs to determine
!  which tasks are assigned to the ocean component.
!
!  this routine should be called after mpi_init, but before
!  setting up any internal mpi setups (since these will require
!  the internal communicators returned by this routine)
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     ierr                    ! error flag for MPI comms

!-----------------------------------------------------------------------
!

#ifdef coupled
   call cpl_interface_init(cpl_fields_ocnname, MPI_COMM_OCN)
#else
   call MPI_INIT(ierr)
   call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_OCN, ierr)
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine create_ocn_communicator

!***********************************************************************
!BOP
! !IROUTINE: create_communicator
! !INTERFACE:

 subroutine create_communicator(new_comm, num_procs)

! !DESCRIPTION:
!  This routine creates a separate communicator for a subset of
!  processors under default ocean communicator.
!
!  this routine should be called from init_domain1 when the
!  domain configuration (e.g. nprocs_btrop) has been determined
!
! !REVISION HISTORY:
!  same as module

! !INCLUDES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      num_procs         ! num of procs in new distribution

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      new_comm          ! new communicator for this distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     MPI_GROUP_OCN,         &! group of processors assigned to ocn
     MPI_GROUP_NEW           ! group of processors assigned to new dist

   integer (int_kind) :: &
     ierr                    ! error flag for MPI comms

   integer (int_kind), dimension(3) :: &
     range                   ! range of tasks assigned to new dist
                             !  (assumed 0,num_procs-1)

!-----------------------------------------------------------------------
!
!  determine group of processes assigned to distribution
!
!-----------------------------------------------------------------------

   call MPI_COMM_GROUP (MPI_COMM_OCN, MPI_GROUP_OCN, ierr)

   range(1) = 0
   range(2) = num_procs-1
   range(3) = 1

!-----------------------------------------------------------------------
!
!  create subroup and communicator for new distribution
!  note: MPI_COMM_CREATE must be called by all procs in MPI_COMM_OCN
!
!-----------------------------------------------------------------------

   call MPI_GROUP_RANGE_INCL(MPI_GROUP_OCN, 1, range, &
                             MPI_GROUP_NEW, ierr)

   call MPI_COMM_CREATE (MPI_COMM_OCN, MPI_GROUP_NEW,  &
                         new_comm, ierr)

!-----------------------------------------------------------------------
!EOC

 end subroutine create_communicator

!***********************************************************************

 end module communicate

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
