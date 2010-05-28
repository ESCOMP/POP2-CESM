!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module io_pio

!BOP
! !MODULE: io_pio 
! !DESCRIPTION:
!  Interfaces for pio initialization
!
! !REVISION HISTORY:
! SVN:$ID: 
!

! !USES:

  use kinds_mod
  use broadcast
  use communicate
  use exit_mod
  use POP_IOUnitsMod
  use io_types
  use pio

  implicit none
  private
  save

!PUBLIC MEMBER FUNCTIONS:

  public io_pio_init
  public ::  io_pio_initdecomp  

!PUBLIC DATA MEMBERS

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

  integer, private      :: io_pio_stride, io_pio_num_iotasks, &
	                   io_pio_root, io_pio_type

  type(iosystem_desc_t), public :: io_pio_subsystem

  integer (i4), parameter :: nmax = 10

  type ptr_ioDesc_int_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_int_type
  type (ptr_ioDesc_int_type), dimension(:), allocatable :: ptr_ioDesc_i

  type ptr_ioDesc_real_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_real_type
  type (ptr_ioDesc_real_type), dimension(:), allocatable :: ptr_ioDesc_r

  type ptr_ioDesc_double_type
     type (IO_desc_t), pointer :: ioDesc(:)
  end type ptr_ioDesc_double_type
  type (ptr_ioDesc_double_type), dimension(:), allocatable :: ptr_ioDesc_d

  integer(i4), parameter :: iunset = -999
  integer(i4), dimension(nmax) :: nsize3d_i = iunset	
  integer(i4), dimension(nmax) :: nsize3d_r = iunset	
  integer(i4), dimension(nmax) :: nsize3d_d = iunset	

  integer(i4), dimension(nmax) :: ksize3d_i = iunset	
  integer(i4), dimension(nmax) :: ksize3d_r = iunset	
  integer(i4), dimension(nmax) :: ksize3d_d = iunset	

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: io_pio_init - initialize io for input or output
! !INTERFACE: 
   subroutine io_pio_init(mode, filename, File, clobber, cdf64)
!
! !DESCRIPTION:
!    Read the pio_inparm namelist and initialize the io subsystem
!
! !REVISION HISTORY:
!    2009-Feb-17 - J. Edwards - initial version
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   character(len=*)     , intent(in)    :: mode
   character(len=*)     , intent(in)    :: filename
   type(file_desc_t)    , intent(inout) :: File
   logical,optional     , intent(in)    :: clobber
   logical,optional     , intent(in)    :: cdf64
!
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error          ! namelist read error flag

   character(len=16)  :: io_pio_type_name

   logical :: exists
   logical :: lclobber
   logical :: lcdf64
   integer :: status
   integer :: nmode
   logical, save :: first_call = .true.
   character(*),parameter :: subName = '(io_pio_init) '

   !  Input namelist
    
   namelist /io_pio_nml/    &
        io_pio_num_iotasks, &
        io_pio_stride,      &
        io_pio_type_name 

!-----------------------------------------------------------------------
!
!  read and define namelist inputs
!
!-----------------------------------------------------------------------

   if (first_call) then

      ! io_pio_root must be set to master_task since since currrently non-standard 
      ! variables are only written from master_task

      io_pio_root =  master_task 

      io_pio_num_iotasks   = -1  ! set based on io_stride value when initialized < 0
      io_pio_stride    = -1  ! set based on num_iotasks value when initialized < 0
      io_pio_type_name = 'netcdf'
      
      if (my_task == master_task) then
         
#ifdef CCSMCOUPLED
         call get_unit(nml_in)
#endif
         open (nml_in, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nml_in, nml=io_pio_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nml_in)

         call io_pio_set_params(io_pio_type_name)

         write(stdout,*) 'POP2 PIO parameter settings...'
         write(stdout,*) '  io_pio_stride    = ',io_pio_stride
         write(stdout,*) '  io_pio_num_iotasks   = ',io_pio_num_iotasks
         write(stdout,*) '  io pio_type_name = ',io_pio_type_name
      endif

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call exit_POP(sigAbort,' error reading pio_nml')
      endif
      call broadcast_scalar(io_pio_num_iotasks, master_task)
      call broadcast_scalar(io_pio_root,        master_task)
      call broadcast_scalar(io_pio_stride,      master_task)
      call broadcast_scalar(io_pio_type,        master_task)
      
      call pio_init(my_task, MPI_COMM_OCN, io_pio_num_iotasks, &
           0, io_pio_stride, PIO_REARR_BOX, io_pio_subsystem, base=io_pio_root)

      first_call = .false.
   end if

   if (trim(mode) == 'write') then
      lclobber = .false.
      if (present(clobber)) lclobber=clobber
   
      lcdf64 = .false.
      if (present(cdf64)) lcdf64=cdf64

      if (File%fh<0) then
         ! filename not open
         inquire(file=trim(filename),exist=exists)
         if (exists) then
            if (lclobber) then
               nmode = pio_clobber
               if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
               status = pio_createfile(io_pio_subsystem, File, io_pio_type, trim(filename), nmode)
               if (my_task == master_task) then
                  write(stdout,*) subname,' create file ',trim(filename)
               end if
            else
               status = pio_openfile(io_pio_subsystem, File, io_pio_type, trim(filename), pio_write)
               if (my_task == master_task) then
                  write(stdout,*) subname,' open file ',trim(filename)
               end if
            endif
         else
            nmode = pio_noclobber
            if (lcdf64) nmode = ior(nmode,PIO_64BIT_OFFSET)
            status = pio_createfile(io_pio_subsystem, File, io_pio_type, trim(filename), nmode)
            if (my_task == master_task) then
               write(stdout,*) subname,' create file ',trim(filename)
            end if
         endif
      else
         ! filename is already open, just return
      endif
   end if

   if (trim(mode) == 'read') then
      inquire(file=trim(filename),exist=exists)
      if (exists) then
         status = pio_openfile(io_pio_subsystem, File, io_pio_type, trim(filename), pio_nowrite)
      else
         if(my_task==master_task) then
            write(stdout,*) 'io_pio_ropen ERROR: file invalid ',trim(filename)
         end if
         call exit_POP(sigAbort, 'aborting in io_pio_ropen with invalid file')
      endif
   end if
      
  end subroutine io_pio_init

!===============================================================================
!BOP
!
! !IROUTINE: io_pio_set_params - set pio parameters
!
! !INTERFACE: 
    subroutine io_pio_set_params(io_pio_type_name)
!
! !DESCRIPTION:
!    Set the pio parameters for the subsystem
!
! !USES:
!
   use shr_string_mod,only: shr_string_toUpper
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   character(len=*), intent(in) :: io_pio_type_name
!
!EOP
!
   character(len=16) :: tmpname
   integer (kind=int_kind) :: npes

   tmpname = shr_string_toupper(io_pio_type_name)

   if (trim(tmpname) == 'NETCDF') then
      io_pio_type = iotype_netcdf
   else if (trim(tmpname) == 'PNETCDF') then
      io_pio_type = iotype_pnetcdf
   else
      if (my_task == master_task) then
         write(stdout,*)' Bad io_type argument - using iotype_netcdf'
      end if
      io_pio_type = iotype_netcdf
   end if

   npes = get_num_procs()
   if      (io_pio_stride>0 .and. io_pio_num_iotasks<0) then
      io_pio_num_iotasks = npes/io_pio_stride
   else if (io_pio_num_iotasks>0 .and. io_pio_stride<0) then
      io_pio_stride = npes/io_pio_num_iotasks
   else if (io_pio_num_iotasks<0 .and. io_pio_stride<0) then
      io_pio_stride = max(min(npes,4),npes/8)
      io_pio_num_iotasks = npes/io_pio_stride
   end if

   if (io_pio_root<0) then
      io_pio_root = 1
   endif
   io_pio_root = min(io_pio_root,npes-1)
   
    if(io_pio_root + (io_pio_stride)*(io_pio_num_iotasks-1) >= npes .or. &
       io_pio_stride<=0 .or. io_pio_num_iotasks<=0 .or. io_pio_root < 0 .or. &
       io_pio_root > npes-1) then
       if (my_task == master_task) then
          write(stdout,*)&
               'io_pio_stride or io_pio_num_iotasks out of bounds, resetting to defaults ',&
               io_pio_stride, io_pio_num_iotasks, io_pio_root
       end if
       io_pio_stride = max(1,npes/4)
       io_pio_num_iotasks = npes/io_pio_stride
       io_pio_root = min(1,npes-1)
    end if
    if (my_task == master_task) then
       write(stdout,*)'Using io_type=',tmpname,' stride=',io_pio_stride,&
            ' iotasks=',io_pio_num_iotasks,' root=',io_pio_root
    end if

   end subroutine io_pio_set_params

!================================================================================

   subroutine io_pio_initdecomp (basetype, ndim3, kdim3, iodesc)

      use blocks, only : block, nx_block, ny_block, get_block
      use domain, only : nblocks_clinic, blocks_clinic
      use POP_DomainSizeMod, only : POP_nxGlobal, POP_nyGlobal  

      integer (i4)          , intent(in) :: basetype
      integer(kind=int_kind), intent(in) :: ndim3
      integer(kind=int_kind), intent(in) :: kdim3
      type(io_desc_t)       , pointer    :: iodesc

      integer (kind=int_kind) :: &
          iblk,ib,ie,jb,je,lon,lat,i,j,n,k,index

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof3d(:)

      logical, save :: first_time = .true.

      logical :: set_iodesc

      if (first_time) then
         allocate(ptr_ioDesc_i(nmax))
         allocate(ptr_ioDesc_r(nmax))
         allocate(ptr_ioDesc_d(nmax))
         do i = 1,nmax
            allocate(ptr_ioDesc_i(i)%ioDesc(1))
            allocate(ptr_ioDesc_r(i)%ioDesc(1))
            allocate(ptr_ioDesc_d(i)%ioDesc(1))
         end do
         first_time = .false.
      end if
      
      if (basetype == PIO_INT) then
         do i = 1,nmax
            if (nsize3d_i(i) == ndim3 .and. ksize3d_i(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_i(i) == iunset .and. ksize3d_i(i) == iunset) then
               index = i
               nsize3d_i(index) = ndim3 
	       ksize3d_i(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      else if (basetype == PIO_REAL) then
         do i = 1,nmax
            if (nsize3d_r(i) == ndim3 .and. ksize3d_r(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_r(i) == iunset .and. ksize3d_r(i) == iunset) then
               index = i
               nsize3d_r(index) = ndim3 
	       ksize3d_r(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      else if (basetype == PIO_DOUBLE) then
         do i = 1,nmax  
            if (nsize3d_d(i) == ndim3 .and. ksize3d_d(i) == kdim3) then
               index = i
               set_ioDesc = .false.
               exit
            else if (nsize3d_d(i) == iunset .and. ksize3d_d(i) == iunset) then
               index = i
               nsize3d_d(index) = ndim3 
	       ksize3d_d(index) = kdim3
               set_ioDesc = .true.
               exit
            end if
         end do
      end if

      if (set_ioDesc) then

	 if ((ndim3 == 0 .and. kdim3 /= 0) .or. (ndim3 /=0 .and. kdim3 == 0)) then
            call exit_POP(sigAbort,' io_pio_initdecomp: ndim3 and kdim3 must both be zero or nonzero')
         end if

	 if (ndim3 > kdim3) then
            call exit_POP(sigAbort,' io_pio_initdecomp: ndim3 must be less than or equal to kdim3')
         end if

         if (ndim3 == 0) then
            allocate(dof3d(nx_block*ny_block*nblocks_clinic))
            n=0
            do iblk = 1, nblocks_clinic
               this_block = get_block(blocks_clinic(iblk),iblk)         
               ib = this_block%ib
               ie = this_block%ie
               jb = this_block%jb
               je = this_block%je
               
               do j=1,ny_block
               do i=1,nx_block  
                  n = n+1
                  if (j < jb .or. j>je) then
                     dof3d(n)=0
                  else if (i < ib .or. i > ie) then
                     dof3d(n) = 0
                  else
                     lon = this_block%i_glob(i)
                     lat = this_block%j_glob(j)
                     dof3d(n) = ((lat-1)*POP_nxGlobal + lon)
                  endif
               enddo !i
               enddo !j
            end do
         else
            allocate(dof3d(nx_block*ny_block*nblocks_clinic*kdim3))
            n=0
            do iblk = 1, nblocks_clinic
               this_block = get_block(blocks_clinic(iblk),iblk)         
               ib = this_block%ib
               ie = this_block%ie
               jb = this_block%jb
               je = this_block%je
               
               do k=1,kdim3
               do j=1,ny_block
               do i=1,nx_block  
                  n = n+1
                  if (j < jb .or. j>je) then
                     dof3d(n)=0
                  else if (i < ib .or. i > ie) then
                     dof3d(n) = 0
                  else
                     if (k > ndim3) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*POP_nxGlobal + lon) + (k-1)*POP_nxGlobal*POP_nyGlobal 
                     end if
                  endif
               enddo !i
               enddo !j
               enddo !kdim3
            end do
         end if

         if (basetype == PIO_INT) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_i(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_i(index)%ioDesc(1))
            end if
         else if (basetype == PIO_REAL) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_r(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_r(index)%ioDesc(1))
            end if
         else if (basetype == PIO_DOUBLE) then
            if (ndim3 == 0) then
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal/), &
                    dof3d, ptr_ioDesc_d(index)%ioDesc(1))
            else
               call pio_initdecomp(io_pio_subsystem, basetype, (/POP_nxGlobal,POP_nyGlobal,ndim3/), &
                    dof3d, ptr_ioDesc_d(index)%ioDesc(1))
            end if
         end if

         deallocate(dof3d)
      end if

      if (basetype == PIO_INT) then
         iodesc => ptr_ioDesc_i(index)%ioDesc(1)
      elseif (basetype == PIO_REAL) then
         iodesc => ptr_ioDesc_r(index)%ioDesc(1)
      elseif (basetype == PIO_DOUBLE) then
         iodesc => ptr_ioDesc_d(index)%ioDesc(1)
      end if

    end subroutine io_pio_initdecomp

!================================================================================

end module io_pio      
