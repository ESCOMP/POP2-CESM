!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_DistributionMod

!BOP
! !MODULE: POP_DistributionMod
!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_DistributionMod.F90 67 2007-03-21 14:39:24Z pwjones $
!  2006-09-06: Phil Jones
!              new naming conventions, more general functionality
!                and some new query capabilities

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_BlocksMod

   implicit none
   private
   save

! !PUBLIC TYPES:

   type, public :: POP_distrb  ! distribution data type
      private

      integer (POP_i4) ::   &
         numProcs          ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (POP_i4), dimension(:), pointer :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_DistributionCreate,      &
             POP_DistributionDestroy,     &
             POP_DistributionGet,         &
             POP_DistributionGetBlockLoc, &
             POP_DistributionGetBlockID

! !DEFINED PARAMETERS:

   ! supported methods for block distribution
   integer (POP_i4), parameter, public ::  &
      POP_distributionMethodNone      = 0, &
      POP_distributionMethodCartesian = 1, &
      POP_distributionMethodRake      = 2

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionCreate
! !INTERFACE:

 function POP_DistributionCreate(distrbMethod, numProcs, workPerBlock, &
                                 errorCode) result(newDistrb)

! !DESCRIPTION:
!  This routine creates a new data structure that describes the 
!  distribution of blocks across processors.  Based on the input
!  choice of method for distribution, an appropriate creation
!  routine is called. Currently only two distributions are supported:
!  2-d Cartesian distribution (cartesian) and a load-balanced
!  distribution using a rake algorithm based on the amount of 
!  work per block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      numProcs,            &! number of processors in this distribution
      distrbMethod          ! method for distributing blocks

   integer (POP_i4), dimension(:), intent(in) :: &
      workPerBlock          ! estimated amount of work per block

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

   type (POP_distrb) :: &
      newDistrb             ! resulting structure describing
                            !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   select case (distrbMethod)

   case(POP_distributionMethodCartesian)

      newDistrb = POP_DistributionCreateCartesian(numProcs, &
                                              workPerBlock, errorCode) 

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
               'POP_DistributionCreate: error in cartesian create')
         return
      endif

   case(POP_distributionMethodRake)

      newDistrb = POP_DistributionCreateRake(numProcs, &
                                             workPerBlock, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
               'POP_DistributionCreate: error in rake create')
         return
      endif

   case default

      call POP_ErrorSet(errorCode, &
               'POP_DistributionCreate: unknown distribution method')
      return

   end select

!-----------------------------------------------------------------------
!EOC

 end function POP_DistributionCreate

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionDestroy
! !INTERFACE:

 subroutine POP_DistributionDestroy(distribution, errorCode)

! !DESCRIPTION:
!  This routine destroys a defined distribution by deallocating
!  all memory associated with the distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   type (POP_distrb), intent(inout) :: &
      distribution          ! distribution to destroy

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: istat  ! status flag for deallocate

!----------------------------------------------------------------------
!
!  reset scalars
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   distribution%numProcs       = 0
   distribution%communicator   = 0
   distribution%numLocalBlocks = 0

!----------------------------------------------------------------------
!
!  deallocate arrays
!
!----------------------------------------------------------------------

   deallocate(distribution%blockLocation, stat=istat)

   if (istat > 0 ) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionDestroy: error deallocating blockLocation')
      return
   endif

   deallocate(distribution%blockLocalID , stat=istat)

   if (istat > 0 ) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionDestroy: error deallocating blockLocalID')
      return
   endif

   deallocate(distribution%blockGlobalID, stat=istat)

   if (istat > 0 ) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionDestroy: error deallocating blockGlobalID')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_DistributionDestroy

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionGet
! !INTERFACE:

 subroutine POP_DistributionGet(distribution, errorCode,            &
                            numProcs, communicator, numLocalBlocks, &
                            blockLocation, blockLocalID, blockGlobalID)


! !DESCRIPTION:
!  This routine extracts information from a distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

! !OUTPUT PARAMETERS:

      integer (POP_i4), intent(out) ::   &
         errorCode           ! returned error code

      integer (POP_i4), intent(out), optional ::   &
         numProcs          ,&! number of processors in this dist
         communicator      ,&! communicator to use in this dist
         numLocalBlocks      ! number of blocks distributed to this
                             !   local processor

      integer (POP_i4), dimension(:), pointer, optional :: &
         blockLocation     ,&! processor location for all blocks
         blockLocalID      ,&! local  block id for all blocks
         blockGlobalID       ! global block id for each local block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  depending on which optional arguments are present, extract the
!  requested info
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (present(numProcs))       numProcs       = distribution%numProcs
   if (present(communicator))   communicator   = distribution%communicator
   if (present(numLocalBlocks)) numLocalBlocks = distribution%numLocalBlocks

   if (present(blockLocation)) then
      if (associated(distribution%blockLocation)) then
         blockLocation => distribution%blockLocation
      else
         call POP_ErrorSet(errorCode, &
            'POP_DistributionGet: blockLocation not allocated')
         return
      endif
   endif

   if (present(blockLocalID)) then
      if (associated(distribution%blockLocalID)) then
         blockLocalID = distribution%blockLocalID
      else
         call POP_ErrorSet(errorCode, &
            'POP_DistributionGet: blockLocalID not allocated')
         return
      endif
   endif

   if (present(blockGlobalID)) then
      if (associated(distribution%blockGlobalID)) then
         blockGlobalID = distribution%blockGlobalID
      else
         call POP_ErrorSet(errorCode, &
            'POP_DistributionGet: blockGlobalID not allocated')
         return
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_DistributionGet

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionGetBlockLoc
! !INTERFACE:

 subroutine POP_DistributionGetBlockLoc(distribution, blockID, &
                                        processor, localID, errorCode)


! !DESCRIPTION:
!  Given a distribution of blocks and a global block ID, return
!  the processor and local index for the block.  A zero for both
!  is returned in the case that the block has been eliminated from
!  the distribution (i.e. has no active points).
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (POP_i4), intent(in) :: &
      blockID                ! global block id for which location requested

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) ::  &
      errorCode              ! returned error code

   integer (POP_i4), intent(out) ::  &
      processor,            &! processor on which block resides
      localID                ! local index for this block on this proc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid blockID
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (blockID < 0 .or. blockID > POP_numBlocks) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionGetBlockLoc: invalid block id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the location from the distribution data structure
!
!-----------------------------------------------------------------------

   processor = distribution%blockLocation(blockID)
   localID   = distribution%blockLocalID (blockID)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_DistributionGetBlockLoc

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionGetBlockID
! !INTERFACE:

 subroutine POP_DistributionGetBlockID(distribution, localID, &
                                       blockID, errorCode)


! !DESCRIPTION:
!  Given a distribution of blocks and a local block index, return
!  the global block id for the block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (POP_distrb), intent(in) :: &
      distribution           ! input distribution for which information
                             !  is requested

   integer (POP_i4), intent(in) ::  &
      localID                ! local index for this block on this proc

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      blockID                ! global block id for this local block

   integer (POP_i4), intent(out) ::  &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for valid localID
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (localID < 0 .or. localID > distribution%numLocalBlocks) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionGetBlockID: invalid local id')
      return
   endif

!-----------------------------------------------------------------------
!
!  extract the global ID from the distribution data structure
!
!-----------------------------------------------------------------------

   blockID   = distribution%blockGlobalID (localID)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_DistributionGetBlockID

!***********************************************************************
!BOP
! !IROUTINE: POP_DistributionCreateCartesian
! !INTERFACE:

 function POP_DistributionCreateCartesian(numProcs, workPerBlock, &
                                          errorCode) result(newDistrb)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      numProcs            ! number of processors in this distribution

   integer (POP_i4), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

   type (POP_distrb) :: &
      newDistrb           ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: &
      i, j,                  &! dummy loop indices
      istat,                 &! status flag for allocation
      iblock, jblock,        &!
      is, ie, js, je,        &! start, end block indices for each proc
      processor,             &! processor position in cartesian decomp
      globalID,              &! global block ID
      localID,               &! block location on this processor
      numProcsX,             &! num of procs in x for global domain
      numProcsY,             &! num of procs in y for global domain
      numBlocksXPerProc,     &! num of blocks per processor in x
      numBlocksYPerProc       ! num of blocks per processor in y

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   call POP_CommCreateCommunicator(newDistrb%communicator, numProcs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   newDistrb%numProcs = numProcs

   call POP_DistributionProcDecomp(numProcs, &
                                      numProcsX, numProcsY, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
        'POP_DistributionCreateCartesian: error in proc decomposition')
      return
   endif

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (newDistrb%blockLocation(POP_numBlocks), &
             newDistrb%blockLocalID (POP_numBlocks), stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateCartesian: error allocating blockLoc')
      return
   endif

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   numBlocksXPerProc = (POP_numBlocksX-1)/numProcsX + 1
   numBlocksYPerProc = (POP_numBlocksY-1)/numProcsY + 1

   do j=1,numProcsY
   do i=1,numProcsX
      processor = (j-1)*numProcsX + i    ! number the processors 
                                         ! left to right, bot to top

      is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
      ie =  i   *numBlocksXPerProc       ! ending   block in i
      if (ie > POP_numBlocksX) ie = POP_numBlocksX
      js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
      je =  j   *numBlocksYPerProc       ! ending   block in j
      if (je > POP_numBlocksY) je = POP_numBlocksY

      localID        = 0  ! initialize counter for local index
      do jblock = js,je
      do iblock = is,ie
         globalID = (jblock - 1)*POP_numBlocksX + iblock
         if (workPerBlock(globalID) /= 0) then
            localID = localID + 1
            newDistrb%blockLocation(globalID) = processor
            newDistrb%blockLocalID (globalID) = localID
         else  ! no work - eliminate block from distribution
            newDistrb%blockLocation(globalID) = 0
            newDistrb%blockLocalID (globalID) = 0
         endif
      end do
      end do

      ! if this is the local processor, set number of local blocks
      if (POP_myTask == processor - 1) then
         newDistrb%numLocalBlocks = localID
      endif

   end do
   end do

!----------------------------------------------------------------------
!
!  now store the local info
!
!----------------------------------------------------------------------

   if (newDistrb%numLocalBlocks > 0) then
      allocate (newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
                stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateCartesian: error allocating blockLoc')
         return
      endif

      do j=1,numProcsY
      do i=1,numProcsX
         processor = (j-1)*numProcsX + i

         if (processor == POP_myTask + 1) then
            is = (i-1)*numBlocksXPerProc + 1   ! starting block in i
            ie =  i   *numBlocksXPerProc       ! ending   block in i
            if (ie > POP_numBlocksX) ie = POP_numBlocksX
            js = (j-1)*numBlocksYPerProc + 1   ! starting block in j
            je =  j   *numBlocksYPerProc       ! ending   block in j
            if (je > POP_numBlocksY) je = POP_numBlocksY

            localID        = 0  ! initialize counter for local index
            do jblock = js,je
            do iblock = is,ie
               globalID = (jblock - 1)*POP_numBlocksX + iblock
               if (workPerBlock(globalID) /= 0) then
                  localID = localID + 1
                  newDistrb%blockGlobalID (localID) = globalID
               endif
            end do
            end do

         endif

      end do
      end do

   endif

!----------------------------------------------------------------------
!EOC

 end function POP_DistributionCreateCartesian

!**********************************************************************
!BOP
! !IROUTINE: POP_DistributionCreateRake
! !INTERFACE:

 function POP_DistributionCreateRake(numProcs, workPerBlock, &
                                     errorCode) result(newDistrb)


! !DESCRIPTION:
!  This  function distributes blocks across processors in a
!  load-balanced manner based on the amount of work per block.
!  A rake algorithm is used in which the blocks are first distributed
!  in a Cartesian distribution and then a rake is applied in each
!  Cartesian direction.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      numProcs                ! number of processors in this distribution

   integer (POP_i4), dimension(:), intent(in) :: &
      workPerBlock        ! amount of work per block

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode           ! returned error code

   type (POP_distrb) :: &
      newDistrb           ! resulting structure describing
                          ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) ::    &
      i,j,n              ,&! dummy loop indices
      pid                ,&! dummy for processor id
      istat              ,&! status flag for allocates
      localBlock         ,&! local block position on processor
      numOcnBlocks       ,&! number of ocean blocks
      maxWork            ,&! max amount of work in any block
      numProcsX          ,&! num of procs in x for global domain
      numProcsY            ! num of procs in y for global domain

   integer (POP_i4), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      workTmp            ,&! work per row or column for rake algrthm
      procTmp              ! temp processor id for rake algrthm

   type (POP_distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   dist = POP_DistributionCreateCartesian(numProcs, &
                                    workPerBlock, errorCode)
   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error creating initial dist')
      return
   endif

!----------------------------------------------------------------------
!
!  if the number of blocks is close to the number of processors,
!  only do a 1-d rake on the entire distribution
!
!----------------------------------------------------------------------

   numOcnBlocks = count(workPerBlock /= 0)

   if (numOcnBlocks <= 2*numProcs) then

      allocate(priority(POP_numBlocks), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error allocating priority array')
         return
      endif

      !*** initialize priority array

      do j=1,POP_numBlocksY
      do i=1,POP_numBlocksX
         n=(j-1)*POP_numBlocksX + i
         if (workPerBlock(n) > 0) then
            priority(n) = maxWork + n - workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(POP_numBlocks), procTmp(POP_numBlocks), stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error allocating tmps')
         return
      endif

      workTmp(:) = 0
      do i=1,numProcs
         procTmp(i) = i
         do n=1,POP_numBlocks
            if (dist%blockLocation(n) == i) then
               workTmp(i) = workTmp(i) + workPerBlock(n)
            endif
         end do
      end do

      call POP_DistributionRake (workTmp, procTmp, workPerBlock, &
                                 priority, dist, errorCode)

      if (errorCode /= POP_Success) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error in linear rake')
         return
      endif

      deallocate(workTmp, procTmp, stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error deallocating tmps')
         return
      endif

!----------------------------------------------------------------------
!
!  otherwise re-distribute blocks using a rake in each direction
!
!----------------------------------------------------------------------

   else

      maxWork = maxval(workPerBlock)

      call POP_DistributionProcDecomp(dist%numProcs, &
                                   numProcsX, numProcsY, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
           'POP_DistributionCreateRake: error in proc decomposition')
         return
      endif

!----------------------------------------------------------------------
!
!     load-balance using a rake algorithm in the x-direction first
!
!----------------------------------------------------------------------

      allocate(priority(POP_numBlocks), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error allocating priority array')
         return
      endif
   
      !*** set highest priority such that eastern-most blocks
      !*** and blocks with the least amount of work are
      !*** moved first

      do j=1,POP_numBlocksY
      do i=1,POP_numBlocksX
         n=(j-1)*POP_numBlocksX + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(POP_numBlocksX + i) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(numProcsX), procTmp(numProcsX), stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error allocating tmps')
         return
      endif

      do j=1,numProcsY

         workTmp(:) = 0
         do i=1,numProcsX
            pid = (j-1)*numProcsX + i
            procTmp(i) = pid
            do n=1,POP_numBlocks
               if (dist%blockLocation(n) == pid) then
                  workTmp(i) = workTmp(i) + workPerBlock(n)
               endif
            end do
         end do

         call POP_DistributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist, errorCode)

         if (errorCode /= POP_Success) then 
            call POP_ErrorSet(errorCode, &
               'POP_DistributionCreateRake: error in first rake')
            return
         endif

      end do
   
      deallocate(workTmp, procTmp, stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error deallocating tmps')
         return
      endif

!----------------------------------------------------------------------
!
!     use a rake algorithm in the y-direction now
!
!----------------------------------------------------------------------

      !*** set highest priority for northern-most blocks

      do j=1,POP_numBlocksY
      do i=1,POP_numBlocksX
         n=(j-1)*POP_numBlocksX + i
         if (workPerBlock(n) > 0) then
            priority(n) = (maxWork + 1)*(POP_numBlocksY + j) - &
                          workPerBlock(n)
         else
            priority(n) = 0
         endif
      end do
      end do

      allocate(workTmp(numProcsY), procTmp(numProcsY), stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error allocating tmps 2')
         return
      endif


      do i=1,numProcsX

         workTmp(:) = 0
         do j=1,numProcsY
            pid = (j-1)*numProcsX + i
            procTmp(j) = pid
            do n=1,POP_numBlocks
               if (dist%blockLocation(n) == pid) then
                  workTmp(j) = workTmp(j) + workPerBlock(n)
               endif
            end do
         end do

         call POP_DistributionRake (workTmp, procTmp, workPerBlock, &
                                    priority, dist, errorCode)

         if (errorCode /= POP_Success) then 
            call POP_ErrorSet(errorCode, &
               'POP_DistributionCreateRake: error in second rake')
            return
         endif


      end do

      deallocate(workTmp, procTmp, priority, stat=istat)

      if (istat > 0) then 
         call POP_ErrorSet(errorCode, &
            'POP_DistributionCreateRake: error deallocating arrays')
         return
      endif
   endif  ! 1d or 2d rake

!----------------------------------------------------------------------
!
!  create new distribution with info extracted from the temporary
!  distribution
!
!----------------------------------------------------------------------

   newDistrb%numProcs     = numProcs
   newDistrb%communicator = dist%communicator

   allocate(newDistrb%blockLocation(POP_numBlocks), &
            newDistrb%blockLocalID(POP_numBlocks), stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error allocating blockLocation')
      return
   endif

   allocate(procTmp(numProcs), stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error allocating last procTmp')
      return
   endif

   procTmp = 0
   localBlock = 0
   do n=1,POP_numBlocks
      pid = dist%blockLocation(n)  ! processor id
      newDistrb%blockLocation(n) = pid

      if (pid > 0) then
         procTmp(pid) = procTmp(pid) + 1
         newDistrb%blockLocalID (n) = procTmp(pid)
      else
         newDistrb%blockLocalID (n) = 0
      endif
   end do

   newDistrb%numLocalBlocks = procTmp(POP_myTask+1)

   if (minval(procTmp) < 1) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: processors left with no blocks')
      return
   endif

   deallocate(procTmp, stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error allocating last procTmp')
      return
   endif

   allocate(newDistrb%blockGlobalID(newDistrb%numLocalBlocks), &
            stat=istat)

   if (istat > 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error allocating blockGlobalID')
      return
   endif

   localBlock = 0
   do n=1,POP_numBlocks
      if (newDistrb%blockLocation(n) == POP_myTask+1) then
         localBlock = localBlock + 1
         newDistrb%blockGlobalID(localBlock) = n
      endif
   end do

!----------------------------------------------------------------------

   call POP_DistributionDestroy(dist, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionCreateRake: error destroying temp distrib')
      return
   endif

!----------------------------------------------------------------------
!EOC

 end function POP_DistributionCreateRake

!**********************************************************************
!BOP
! !IROUTINE: POP_DistributionProcDecomp
! !INTERFACE:

 subroutine POP_DistributionProcDecomp(numProcs, &
                                   numProcsX, numProcsY, errorCode)

! !DESCRIPTION:
!  This subroutine attempts to find an optimal (nearly square)
!  2d processor decomposition for a given number of processors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      numProcs                       ! total number or processors

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      numProcsX, numProcsY,         &! number of procs in each dimension
      errorCode                      ! returned error code

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: &
      iguess, jguess               ! guesses for numProcsX,Y

   real (POP_r8) :: &
      squareRootN                  ! square root of numProcs

!----------------------------------------------------------------------
!
!  start with an initial guess that is closest to square decomp
!
!----------------------------------------------------------------------

   squareRootN = sqrt(real(numProcs))
   numProcsX = 0
   numProcsY = 0

   iguess = nint(squareRootN)

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   procLoop: do
      jguess = numProcs/iguess

      if (iguess*jguess == numProcs) then ! valid decomp

         !***
         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         !***

         if (mod(POP_numBlocksX,iguess) == 0 .and. &
             mod(POP_numBlocksY,jguess) == 0) then
            numProcsX = iguess
            numProcsY = jguess
            exit procLoop

         !***
         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         !***

         else if (mod(POP_numBlocksX,jguess) == 0 .and. &
                  mod(POP_numBlocksY,iguess) == 0) then
            numProcsX = jguess
            numProcsY = iguess
            exit procLoop

         !***
         !*** A valid decomposition, but keep searching for
         !***  a better one
         !***

         else
            if (numProcsX == 0) then
               numProcsX = iguess
               numProcsY = jguess
            endif
            iguess = iguess - 1
            if (iguess == 0) then
               exit procLoop
            else
               cycle procLoop
            endif
         endif

      else ! invalid decomp - keep trying

         iguess = iguess - 1
         if (iguess == 0) then
            exit procLoop
         else
            cycle procLoop
         endif
      endif
   end do procLoop

   if (numProcsX == 0) then
      call POP_ErrorSet(errorCode, &
         'POP_DistributionProcDecomp: Unable to find 2d config')
      return
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine POP_DistributionProcDecomp

!**********************************************************************
!BOP
! !IROUTINE: POP_DistributionRake
! !INTERFACE:

 subroutine POP_DistributionRake (procWork, procID, blockWork, &
                                  priority, distribution, errorCode)

! !DESCRIPTION:
!  This subroutine performs a rake algorithm to distribute the work
!  along a vector of processors.  In the rake algorithm, a work
!  threshold is first set.  Then, moving from left to right, work
!  above that threshold is raked to the next processor in line.
!  The process continues until the end of the vector is reached
!  and then the threshold is reduced by one for a second rake pass.
!  In this implementation, a priority for moving blocks is defined
!  such that the rake algorithm chooses the highest priority
!  block to be moved to the next processor.  This can be used
!  for example to always choose the eastern-most block or to
!  ensure a block does not stray too far from its neighbors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in), dimension(:) :: &
      blockWork          ,&! amount of work per block
      procID               ! global processor number

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout), dimension(:) :: &
      procWork           ,&! amount of work per processor
      priority             ! priority for moving a given block

   type (POP_distrb), intent(inout) :: &
      distribution         ! distribution to change

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (POP_i4) :: &
      i, n,                  &! dummy loop indices
      np1,                   &! n+1 corrected for cyclical wrap
      iproc, inext,          &! processor ids for current and next 
      numProcs, numBlocks,   &! number of blocks, processors
      lastPriority,          &! priority for most recent block
      minPriority,           &! minimum priority
      lastLoc,               &! location for most recent block
      meanWork, maxWork,     &! mean,max work per processor
      diffWork, residual,    &! work differences and residual work
      numTransfers            ! counter for number of block transfers

!----------------------------------------------------------------------
!
!  initialization
!
!----------------------------------------------------------------------

   errorCode = POP_Success

   numProcs  = size(procWork)
   numBlocks = size(blockWork)

   !*** compute mean,max work per processor

   meanWork = sum(procWork)/numProcs + 1
   maxWork  = maxval(procWork)
   residual = mod(meanWork,numProcs)

   minPriority = 1000000
   do n=1,numProcs
      iproc = procID(n)
      do i=1,numBlocks
         if (distribution%blockLocation(i) == iproc) then
            minPriority = min(minPriority,priority(i))
         endif
      end do
   end do

!----------------------------------------------------------------------
!
!  do two sets of transfers
!
!----------------------------------------------------------------------

   transferLoop: do

!----------------------------------------------------------------------
!
!     do rake across the processors
!
!----------------------------------------------------------------------

      numTransfers = 0
      do n=1,numProcs
         if (n < numProcs) then
            np1   = n+1
         else
            np1   = 1
         endif
         iproc = procID(n)
         inext = procID(np1)

         if (procWork(n) > meanWork) then !*** pass work to next

            diffWork = procWork(n) - meanWork

            rake1: do while (diffWork > 1)

               !*** attempt to find a block with the required
               !*** amount of work and with the highest priority
               !*** for moving (eg boundary blocks first)

               lastPriority = 0
               lastLoc = 0

               do i=1,numBlocks
                  if (distribution%blockLocation(i) == iproc) then
                     if (priority(i) > lastPriority ) then
                        lastPriority = priority(i)
                        lastLoc = i
                     endif
                  endif
               end do
               if (lastLoc == 0) exit rake1 ! could not shift work

               numTransfers = numTransfers + 1
               distribution%blockLocation(lastLoc) = inext
               if (np1 == 1) priority(lastLoc) = minPriority
               diffWork = diffWork - blockWork(lastLoc)

               procWork(n  ) = procWork(n  )-blockWork(lastLoc)
               procWork(np1) = procWork(np1)+blockWork(lastLoc)
            end do rake1
         endif

      end do

!----------------------------------------------------------------------
!
!     increment meanWork by one and repeat
!
!----------------------------------------------------------------------

      meanWork = meanWork + 1
      if (numTransfers == 0 .or. meanWork > maxWork) exit transferLoop

   end do transferLoop

!----------------------------------------------------------------------
!EOC

end subroutine POP_DistributionRake

!***********************************************************************

end module POP_DistributionMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
