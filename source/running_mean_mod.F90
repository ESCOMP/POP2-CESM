!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module running_mean_mod

!BOP
! !MODULE: running_mean_mod
! !DESCRIPTION:
!  Provide an interface to define, initialize, accumulate, and extract
!  running means of user-specified variables.
!
!  Running mean variables are declared with a string name, rank, and
!  timescale. The definition subroutine returns an integer index which is
!  subsequently used by the user to refer to the running mean variable in
!  calls to the init, update, and get subroutines. This index is retrievable
!  with the get_index subroutine.
!
!  Running means variables can be initialized by reading values from a file,
!  or by having values passed to the init subroutine. netcdf is the only
!  file format supported for file reads.
!
!  If a running mean variable is not initialized, its running mean is set
!  to the value provided in the first call to the update subroutine. It is
!  a fatal error to call the get subroutine for an uninitialized running
!  mean variable.
!
!  Subroutines provide support for 0D (scalar), 1D (km), 2D (nx_global,ny_global),
!  and 3D (nx_global,ny_global,km) variables. Subroutines are only provided
!  for real variables of kind r8.
!
!  Subroutines provide support to be called from within threaded regions
!  (i.e. works on a single block) as well as from non-threaded region
!  (i.e. works on all blocks with a single call).
!
!  If running_mean_test_mode is set to .true., then a test variable of
!  each rank will be defined, initialized, and updated. Two 2D test
!  variables are implemented, one that update within surface forcing
!  subroutines. All of the other test variables are updated after the time
!  manager flags are updated. The tavg variables RUNNING_MEAN_TEST_VAR_0D,
!  RUNNING_MEAN_TEST_VAR_2D, and RUNNING_MEAN_TEST_VAR_2D_SFLUX can be added
!  to tavg_contents.
!
! !REVISION HISTORY:
!  SVN:$Id:$
!

! !USES:

   use kinds_mod, only: log_kind, int_kind, char_len, r8
   use io_tools, only: document
   use exit_mod, only: exit_POP, sigAbort
   use domain_size, only: km, nx_global, ny_global
   use blocks, only: nx_block, ny_block
   use domain, only: nblocks_clinic
   use constants, only: c0, p5, c1, c2, c3, c4

   use io_types, only: datafile, io_dim, io_field_desc
   use io_types, only: construct_file, construct_io_dim, construct_io_field
   use io_types, only: destroy_file, destroy_io_field
   use io, only: data_set

   use time_management, only: lpre_time_manager, avg_ts_next, back_to_back_next, avg_ts, back_to_back, dtt

   use tavg, only: define_tavg_field, accumulate_tavg_field

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: running_mean_init
   public :: running_mean_test_update_var
   public :: running_mean_test_update_sflux_var
   public :: running_mean_define_var
   public :: running_mean_get_var_index
   public :: running_mean_var_exists_in_file
   public :: running_mean_init_var
   public :: running_mean_update_var
   public :: running_mean_get_var
   public :: running_mean_write_restart

!EOP
!BOC
!-----------------------------------------------------------------------
!  module private types and data
!-----------------------------------------------------------------------

   logical (log_kind) :: running_mean_test_mode

   integer (int_kind) :: test_index_0d
   integer (int_kind) :: test_index_1d
   integer (int_kind) :: test_index_2d
   integer (int_kind) :: test_index_2d_sflux
   integer (int_kind) :: test_index_3d
   integer (int_kind) :: tavg_test_0d
   integer (int_kind) :: tavg_test_2d
   integer (int_kind) :: tavg_test_2d_sflux

   integer (int_kind), parameter :: running_mean_cnt_max = 100
   integer (int_kind) :: running_mean_cnt = 0

   type :: running_mean_type
      character (char_len)                        :: name
      integer (int_kind)                          :: rank
      real (r8)                                   :: timescale ! [seconds]
      character (char_len)                        :: file_varname ! name of variable in restart file
      logical (log_kind)                          :: linit_0d ! have vals been initialized?
      logical (log_kind), dimension(:), pointer   :: linit_1d ! k dimension
      logical (log_kind), dimension(:), pointer   :: linit_2d ! block dimension
      logical (log_kind), dimension(:,:), pointer :: linit_3d ! k and block dimensions
      real (r8)                                   :: vals_0d
      real (r8), dimension(:), pointer            :: vals_1d
      real (r8), dimension(:,:,:), pointer        :: vals_2d ! includes block index
      real (r8), dimension(:,:,:,:), pointer      :: vals_3d ! includes block index
   end type

   type (running_mean_type), dimension(running_mean_cnt_max), target :: &
      running_mean_array   ! array of running mean variables

!-----------------------------------------------------------------------
!  generic interface definitions
!-----------------------------------------------------------------------

   interface running_mean_init_var
      module procedure running_mean_init_var_filename, &
                       running_mean_init_var_vals
   end interface

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: running_mean_init
! !INTERFACE:

   subroutine running_mean_init(init_ts_file_fmt, read_restart_filename)

! !DESCRIPTION:
!  Initialize running mean module, only necessary for test mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: block ! loop index
   integer (int_kind) :: k     ! loop index

   real (r8) :: tau
   real (r8), dimension(nx_block,ny_block) :: FIELD

!-----------------------------------------------------------------------

   running_mean_test_mode = .false.

   if (running_mean_test_mode) then
      tau = 365.0_r8*86400.0_r8 ! 1 year
      call running_mean_define_var('RUNNING_MEAN_TEST_VAR_0D', 0, tau, test_index_0d)
      call running_mean_define_var('RUNNING_MEAN_TEST_VAR_1D', 1, c2*tau, test_index_1d)
      call running_mean_define_var('RUNNING_MEAN_TEST_VAR_2D', 2, c3*tau, test_index_2d)
      call running_mean_define_var('RUNNING_MEAN_TEST_VAR_2D_SFLUX', 2, c3*tau, test_index_2d_sflux)
      call running_mean_define_var('RUNNING_MEAN_TEST_VAR_3D', 3, c4*tau, test_index_3d)

      call define_tavg_field(tavg_test_0d, 'RUNNING_MEAN_TEST_VAR_0D', 0, &
                             long_name='RUNNING_MEAN_TEST_VAR_0D', units='1', &
                             grid_loc='0000')

      call define_tavg_field(tavg_test_2d, 'RUNNING_MEAN_TEST_VAR_2D', 2, &
                             long_name='RUNNING_MEAN_TEST_VAR_2D', units='1', &
                             grid_loc='2110', coordinates='TLONG TLAT time')

      call define_tavg_field(tavg_test_2d_sflux, 'RUNNING_MEAN_TEST_VAR_2D_SFLUX', 2, &
                             long_name='RUNNING_MEAN_TEST_VAR_2D_SFLUX', units='1', &
                             grid_loc='2110', coordinates='TLONG TLAT time')

      ! initialize test vars from ts restart file, if present
      if (read_restart_filename /= 'undefined') then

         call document('running_mean_init', 'read_restart_filename', read_restart_filename)
         if (init_ts_file_fmt /= 'nc') then
           call document('running_mean_init', 'init_ts_file_fmt', init_ts_file_fmt)
           call exit_POP(sigAbort, 'unsupported init_ts_file_fmt')
         endif
         call running_mean_init_var(test_index_0d, read_restart_filename)
         call running_mean_init_var(test_index_1d, read_restart_filename)
         call running_mean_init_var(test_index_2d, read_restart_filename)
         call running_mean_init_var(test_index_2d_sflux, read_restart_filename)
         call running_mean_init_var(test_index_3d, read_restart_filename)

      else ! otherwise initialize running means to zero

         call running_mean_init_var(test_index_0d, vals_0d=c0)
         do k=1,km
            call running_mean_init_var(test_index_1d, k=k, vals_1d_1klev=c0)
         end do
         FIELD = c0
         do block=1,nblocks_clinic
            call running_mean_init_var(test_index_2d, block=block, vals_2d_1block=FIELD)
            call running_mean_init_var(test_index_2d_sflux, block=block, vals_2d_1block=FIELD)
         end do
         do block=1,nblocks_clinic
            do k=1,km
               call running_mean_init_var(test_index_3d, k=k, block=block, vals_3d_1klev_1block=FIELD)
            end do
         end do
      endif

   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_init

!***********************************************************************
!BOP
! !IROUTINE: running_mean_test_update_sflux_var
! !INTERFACE:

   subroutine running_mean_test_update_sflux_var

! !DESCRIPTION:
!  Call running_mean_update_var for sflux test vars, only necessary for test mode.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,nblocks_clinic) :: FIELD

!-----------------------------------------------------------------------

   if (running_mean_test_mode) then
      ! update fields, corresponding tavg call is in running_mean_test_update_var
      FIELD = c1
      call running_mean_update_var(test_index_2d_sflux, vals_2d_blocks=FIELD)
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_test_update_sflux_var

!***********************************************************************
!BOP
! !IROUTINE: running_mean_test_update_var
! !INTERFACE:

   subroutine running_mean_test_update_var(k, block)

! !DESCRIPTION:
!  Call running_mean_update_var for test vars, only necessary for test mode.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k
   integer (int_kind), intent(in) :: block

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8) :: val_0d
   real (r8), dimension(nx_block,ny_block) :: FIELD

!-----------------------------------------------------------------------

   if (running_mean_test_mode) then
      ! update fields
      FIELD = c1
      if (k == 1) then
         if (block == 1) call running_mean_update_var(test_index_0d, vals_0d=c1)
         call running_mean_update_var(test_index_2d, block=block, vals_2d_1block=FIELD)
      endif
      if (block == 1) call running_mean_update_var(test_index_1d, k=k, vals_1d_1klev=c1)
      call running_mean_update_var(test_index_3d, k=k, block=block, vals_3d_1klev_1block=FIELD)

      ! accumulate corresponding 0d and 2d tavg vars
      if (k == 1) then
         if (block == 1) then
            call running_mean_get_var(test_index_0d, vals_0d=val_0d)
            call accumulate_tavg_field(val_0d, tavg_test_0d)
         endif

         call running_mean_get_var(test_index_2d, block=block, vals_2d_1block=FIELD)
         call accumulate_tavg_field(FIELD, tavg_test_2d, block, k)

         call running_mean_get_var(test_index_2d_sflux, block=block, vals_2d_1block=FIELD)
         call accumulate_tavg_field(FIELD, tavg_test_2d_sflux, block, k)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_test_update_var

!***********************************************************************
!BOP
! !IROUTINE: running_mean_define_var
! !INTERFACE:

   subroutine running_mean_define_var(name, rank, timescale, index)

! !DESCRIPTION:
!  Define a running mean. It is a fatal error to attempt to define a
!  previously defined name.
!
!  This should only be called once per task.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)       :: name      ! variable name
   integer (int_kind), intent(in)  :: rank      ! rank of variable
   real (r8), intent(in)           :: timescale ! running mean timescale

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: index     ! returned index

!EOP
!BOC
!-----------------------------------------------------------------------

   call document('running_mean_define_var', 'name', name)

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   if (rank < 0 .or. rank > 3) then
      call document('running_mean_define_var', 'rank', rank)
      call exit_POP(sigAbort, 'unsupported rank')
   endif

   ! avoid spaces in name, these cause problems for creating restart file
   ! variable names based on name
   if (scan(name, ' ') > 0) then
      call exit_POP(sigAbort, 'spaces are not allowed in running_mean names')
   endif

   ! check to see if name is already registered
   call running_mean_get_var_index(name, index, exit_on_err=.false.)
   if (index > 0) then
      call exit_POP(sigAbort, 'running mean name already defined')
   endif

   running_mean_cnt = running_mean_cnt + 1
   if (running_mean_cnt > running_mean_cnt_max) then
      call document('running_mean_define_var', 'running_mean_cnt_max', &
                    running_mean_cnt_max)
      call exit_POP(sigAbort, 'too many running mean variables defined')
   endif

!-----------------------------------------------------------------------
!  setup new running mean
!-----------------------------------------------------------------------

   running_mean_array(running_mean_cnt)%name = name
   running_mean_array(running_mean_cnt)%rank = rank
   running_mean_array(running_mean_cnt)%timescale = timescale

   running_mean_array(running_mean_cnt)%file_varname = 'running_mean_' /&
      &/ name

   nullify(running_mean_array(running_mean_cnt)%linit_1d, &
           running_mean_array(running_mean_cnt)%linit_2d, &
           running_mean_array(running_mean_cnt)%linit_3d, &
           running_mean_array(running_mean_cnt)%vals_1d, &
           running_mean_array(running_mean_cnt)%vals_2d, &
           running_mean_array(running_mean_cnt)%vals_3d)

   select case (rank)
   case(0)
      running_mean_array(running_mean_cnt)%linit_0d = .false.
      running_mean_array(running_mean_cnt)%vals_0d = c0
   case(1)
      allocate(running_mean_array(running_mean_cnt)%linit_1d(km))
      running_mean_array(running_mean_cnt)%linit_1d(:) = .false.
      allocate(running_mean_array(running_mean_cnt)%vals_1d(km))
      running_mean_array(running_mean_cnt)%vals_1d(:) = c0
   case(2)
      allocate(running_mean_array(running_mean_cnt)%linit_2d(nblocks_clinic))
      running_mean_array(running_mean_cnt)%linit_2d(:) = .false.
      allocate(running_mean_array(running_mean_cnt)%vals_2d(nx_block,ny_block,nblocks_clinic))
      running_mean_array(running_mean_cnt)%vals_2d(:,:,:) = c0
   case(3)
      allocate(running_mean_array(running_mean_cnt)%linit_3d(km,nblocks_clinic))
      running_mean_array(running_mean_cnt)%linit_3d(:,:) = .false.
      allocate(running_mean_array(running_mean_cnt)%vals_3d(nx_block,ny_block,km,nblocks_clinic))
      running_mean_array(running_mean_cnt)%vals_3d(:,:,:,:) = c0
   end select

   index = running_mean_cnt

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_define_var

!***********************************************************************
!BOP
! !IROUTINE: running_mean_get_var_index
! !INTERFACE:

   subroutine running_mean_get_var_index(name, index, exit_on_err)

! !DESCRIPTION:
!  Search for a running mean variable by name and return its index.
!  If the name is not found then exit_POP is called, unless
!     the optional argument exit_on_err is set to .false., in which
!     case index is set to 0.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)                :: name         ! name of variable to be looked up
   logical (log_kind), intent(in), optional :: exit_on_err  ! Is exit_POP called if name not found?

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out)          :: index

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   logical (log_kind) :: loc_exit_on_err   ! local copy of exit_on_err
   integer (int_kind) :: i                 ! loop index

!-----------------------------------------------------------------------

   if (.not. present(exit_on_err)) then
      loc_exit_on_err = .true.
   else
      loc_exit_on_err = exit_on_err
   endif

   index = 0
   do i=1,running_mean_cnt
      if (running_mean_array(i)%name == name) then
         index = i
         exit
      endif
   end do

   if (index == 0 .and. loc_exit_on_err) then
      call document('running_mean_get_var_index', 'name', name)
      call exit_POP(sigAbort, 'name not found')
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_get_var_index

!***********************************************************************
!BOP
! !IROUTINE: running_mean_var_exists_in_file
! !INTERFACE:

 function running_mean_var_exists_in_file(index, filename)

! !DESCRIPTION:
!  Determine if a running mean variable exists in a file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: index
   character (*), intent(in)      :: filename

! !OUTPUT PARAMETERS:

   logical (log_kind) :: running_mean_var_exists_in_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile) :: file ! io file descriptor

!-----------------------------------------------------------------------

   file = construct_file('nc', full_name=trim(filename))

   call data_set(file, 'open_read')

   call data_set(file, 'field_exists', fieldname=running_mean_array(index)%file_varname, &
                 field_exists=running_mean_var_exists_in_file)

   call data_set(file, 'close')

   call destroy_file(file)

!-----------------------------------------------------------------------
!EOC

 end function running_mean_var_exists_in_file

!***********************************************************************
!BOP
! !IROUTINE: running_mean_init_var_filename
! !INTERFACE:

   subroutine running_mean_init_var_filename(index, filename)

! !DESCRIPTION:
!  Initialize a running mean variable from a file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: index
   character (*), intent(in)      :: filename

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type (datafile)      :: file         ! io file descriptor
   type (io_dim)        :: i_dim        ! dimension descriptor
   type (io_dim)        :: j_dim        ! dimension descriptor
   type (io_dim)        :: k_dim        ! dimension descriptor
   type (io_field_desc) :: vals_desc    ! io field descriptor

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   if (index < 1 .or. index > running_mean_cnt) then
      call document('running_mean_init_var_filename', 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
   endif

!-----------------------------------------------------------------------
!  initialize values from file
!-----------------------------------------------------------------------

   call document('running_mean_init_var_filename', 'name', running_mean_array(index)%name)
   call document('running_mean_init_var_filename', 'file_varname', running_mean_array(index)%file_varname)

   file = construct_file('nc', full_name=trim(filename))

   call data_set(file, 'open_read')

   select case (running_mean_array(index)%rank)
   case(0)
      running_mean_array(index)%linit_0d = .true.
      vals_desc = &
         construct_io_field(running_mean_array(index)%file_varname, &
                            d0d_array=running_mean_array(index)%vals_0d)
   case(1)
      running_mean_array(index)%linit_1d(:) = .true.
      k_dim = construct_io_dim('k', km)
      vals_desc = &
         construct_io_field(running_mean_array(index)%file_varname, dim1=k_dim, &
                            d1d_array=running_mean_array(index)%vals_1d)
   case(2)
      running_mean_array(index)%linit_2d(:) = .true.
      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)
      vals_desc = &
         construct_io_field(running_mean_array(index)%file_varname, dim1=i_dim, dim2=j_dim, &
                            d2d_array=running_mean_array(index)%vals_2d)
   case(3)
      running_mean_array(index)%linit_3d(:,:) = .true.
      i_dim = construct_io_dim('i', nx_global)
      j_dim = construct_io_dim('j', ny_global)
      k_dim = construct_io_dim('k', km)
      vals_desc = &
         construct_io_field(running_mean_array(index)%file_varname, dim1=i_dim, dim2=j_dim, dim3=k_dim, &
                            d3d_array=running_mean_array(index)%vals_3d)
   end select

   call data_set(file, 'define', vals_desc)
   call data_set(file, 'read', vals_desc)
   call destroy_io_field(vals_desc)
   call data_set(file, 'close')
   call destroy_file(file)

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_init_var_filename

!***********************************************************************
!BOP
! !IROUTINE: running_mean_init_var_vals
! !INTERFACE:

   subroutine running_mean_init_var_vals(index, block, k, vals_0d, &
      vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
      vals_3d_1klev_1block, vals_3d_1klev_blocks)

! !DESCRIPTION:
!  Initialize a running mean variable from specified values.
!  A scalar value can be provided to initialize a variable of any rank.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in)                      :: index
   integer (int_kind), intent(in), optional            :: block
   integer (int_kind), intent(in), optional            :: k
   real (r8), intent(in), optional                     :: vals_0d
   real (r8), intent(in), optional                     :: vals_1d_1klev
   real (r8), intent(in), dimension(:), optional       :: vals_1d_klevs
   real (r8), intent(in), dimension(:,:), optional     :: vals_2d_1block
   real (r8), intent(in), dimension(:,:,:), optional   :: vals_2d_blocks
   real (r8), intent(in), dimension(:,:), optional     :: vals_3d_1klev_1block
   real (r8), intent(in), dimension(:,:,:), optional   :: vals_3d_1klev_blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock ! loop index

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   if (index < 1 .or. index > running_mean_cnt) then
      call document('running_mean_init_var_vals', 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
   endif

   call document('running_mean_init_var_vals', 'name', running_mean_array(index)%name)

   select case (running_mean_array(index)%rank)
   case(0)
      if (.not. present(vals_0d)) then
         call exit_POP(sigAbort, 'vals_0d must be supplied for 0d vars')
      endif
   case(1)
      if (.not. present(vals_0d) .and. .not. present(vals_1d_1klev) .and. .not. present(vals_1d_klevs)) then
         call exit_POP(sigAbort, 'vals_0d, vals_1d_1klev, or vals_1d_klevs must be supplied for 1d vars')
      endif
      if (present(vals_1d_1klev)) then
         if (.not. present(k)) then
            call exit_POP(sigAbort, 'k must be supplied if vals_1d_1klev is supplied')
         else
            call document('running_mean_init_var_vals', 'k', k)
         endif
      endif
   case(2)
      if (.not. present(vals_0d) .and. .not. present(vals_2d_1block) .and. .not. present(vals_2d_blocks)) then
         call exit_POP(sigAbort, 'vals_0d, vals_2d_1block, or vals_2d_blocks must be supplied for 2d vars')
      endif
      if (present(vals_2d_1block)) then
         if (.not. present(block)) then
            call exit_POP(sigAbort, 'block must be supplied if vals_2d_1block is supplied')
         else
            call document('running_mean_init_var_vals', 'block', block)
         endif
      endif
   case(3)
      if (.not. present(vals_0d) .and. .not. present(vals_3d_1klev_1block) .and. .not. present(vals_3d_1klev_blocks)) then
         call exit_POP(sigAbort, 'vals_0d, vals_3d_1klev_1block, or vals_3d_1klev_blocks must be supplied for 3d vars')
      endif
      if (present(vals_3d_1klev_1block)) then
         if (.not. present(k)) then
            call exit_POP(sigAbort, 'k must be supplied if vals_3d_1klev_1block is supplied')
         else
            call document('running_mean_init_var_vals', 'k', k)
         endif
         if (.not. present(block)) then
            call exit_POP(sigAbort, 'block must be supplied if vals_3d_1klev_1block is supplied')
         else
            call document('running_mean_init_var_vals', 'block', block)
         endif
      endif
      if (present(vals_3d_1klev_blocks)) then
         if (.not. present(k)) then
            call exit_POP(sigAbort, 'k must be supplied if vals_3d_1klev_blocks is supplied')
         else
            call document('running_mean_init_var_vals', 'k', k)
         endif
      endif
   end select

!-----------------------------------------------------------------------
!  initialize values from subroutine arguments
!-----------------------------------------------------------------------

   select case (running_mean_array(index)%rank)
   case(0)
      running_mean_array(index)%linit_0d = .true.
      running_mean_array(index)%vals_0d = vals_0d
   case(1)
      if (present(vals_0d)) then
         running_mean_array(index)%linit_1d(:) = .true.
         running_mean_array(index)%vals_1d(:) = vals_0d
      else if (present(vals_1d_1klev)) then
         running_mean_array(index)%linit_1d(k) = .true.
         running_mean_array(index)%vals_1d(k) = vals_1d_1klev
      else
         running_mean_array(index)%linit_1d(:) = .true.
         running_mean_array(index)%vals_1d(:) = vals_1d_klevs(:)
      endif
   case(2)
      if (present(vals_0d)) then
         running_mean_array(index)%linit_2d(:) = .true.
         running_mean_array(index)%vals_2d(:,:,:) = vals_0d
      else if (present(vals_2d_1block)) then
         running_mean_array(index)%linit_2d(block) = .true.
         running_mean_array(index)%vals_2d(:,:,block) = vals_2d_1block(:,:)
      else
         do iblock=1,nblocks_clinic
            running_mean_array(index)%linit_2d(iblock) = .true.
            running_mean_array(index)%vals_2d(:,:,iblock) = vals_2d_blocks(:,:,iblock)
         end do
      endif
   case(3)
      if (present(vals_0d)) then
         running_mean_array(index)%linit_3d(:,:) = .true.
         running_mean_array(index)%vals_3d(:,:,:,:) = vals_0d
      else if (present(vals_3d_1klev_1block)) then
         running_mean_array(index)%linit_3d(k,block) = .true.
         running_mean_array(index)%vals_3d(:,:,k,block) = vals_3d_1klev_1block(:,:)
      else
         do iblock=1,nblocks_clinic
            running_mean_array(index)%linit_3d(k,iblock) = .true.
            running_mean_array(index)%vals_3d(:,:,k,iblock) = vals_3d_1klev_blocks(:,:,iblock)
         end do
      endif
   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_init_var_vals

!***********************************************************************
!BOP
! !IROUTINE: running_mean_update_var
! !INTERFACE:

   subroutine running_mean_update_var(index, block, k, vals_0d, &
      vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
      vals_3d_1klev_1block, vals_3d_1klev_blocks)

! !DESCRIPTION:
!  Update a running mean variable using specified values.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in)                      :: index
   integer (int_kind), intent(in), optional            :: block
   integer (int_kind), intent(in), optional            :: k
   real (r8), intent(in), optional                     :: vals_0d
   real (r8), intent(in), optional                     :: vals_1d_1klev
   real (r8), intent(in), dimension(:), optional       :: vals_1d_klevs
   real (r8), intent(in), dimension(:,:), optional     :: vals_2d_1block
   real (r8), intent(in), dimension(:,:,:), optional   :: vals_2d_blocks
   real (r8), intent(in), dimension(:,:), optional     :: vals_3d_1klev_1block
   real (r8), intent(in), dimension(:,:,:), optional   :: vals_3d_1klev_blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock ! loop index
   integer (int_kind) :: k_loc  ! loop index
   real (r8)          :: weight ! weight applied to current running mean

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   if (index < 1 .or. index > running_mean_cnt) then
      call document('running_mean_update_var', 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
   endif

   select case (running_mean_array(index)%rank)
   case(0)
      if (.not. present(vals_0d)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_0d must be supplied for 0d vars')
      endif
   case(1)
      if (.not. present(vals_1d_1klev) .and. .not. present(vals_1d_klevs)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_1d_1klev or vals_1d_klevs must be supplied for 1d vars')
      endif
      if (present(vals_1d_1klev) .and. .not. present(k)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'k must be supplied if vals_1d_1klev is supplied')
      endif
   case(2)
      if (.not. present(vals_2d_1block) .and. .not. present(vals_2d_blocks)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_2d_1block or vals_2d_blocks must be supplied for 2d vars')
      endif
      if (present(vals_2d_1block) .and. .not. present(block)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'block must be supplied if vals_2d_1block is supplied')
      endif
   case(3)
      if (.not. present(vals_3d_1klev_1block) .and. .not. present(vals_3d_1klev_blocks)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_3d_1klev_1block or vals_3d_1klev_blocks must be supplied for 3d vars')
      endif
      if (present(vals_3d_1klev_1block) .and. .not. present(k)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'k must be supplied if vals_3d_1klev_1block is supplied')
      endif
      if (present(vals_3d_1klev_1block) .and. .not. present(block)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'block must be supplied if vals_3d_1klev_1block is supplied')
      endif
      if (present(vals_3d_1klev_blocks) .and. .not. present(k)) then
         call document('running_mean_update_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'k must be supplied if vals_3d_1klev_blocks is supplied')
      endif
   end select

!-----------------------------------------------------------------------
!  update values from subroutine arguments
!-----------------------------------------------------------------------

   if (lpre_time_manager) then
      if (avg_ts_next .or. back_to_back_next) then
         weight = exp(-p5*dtt/running_mean_array(index)%timescale)
      else
         weight = exp(-dtt/running_mean_array(index)%timescale)
      endif
   else
      if (avg_ts .or. back_to_back) then
         weight = exp(-p5*dtt/running_mean_array(index)%timescale)
      else
         weight = exp(-dtt/running_mean_array(index)%timescale)
      endif
   endif

   select case (running_mean_array(index)%rank)
   case(0)
      if (running_mean_array(index)%linit_0d) then
         running_mean_array(index)%vals_0d = running_mean_array(index)%vals_0d &
            + (c1 - weight) * (vals_0d - running_mean_array(index)%vals_0d)
      else
         running_mean_array(index)%vals_0d = vals_0d
         running_mean_array(index)%linit_0d = .true.
      endif
   case(1)
      if (present(vals_1d_1klev)) then
         if (running_mean_array(index)%linit_1d(k)) then
            running_mean_array(index)%vals_1d(k) = running_mean_array(index)%vals_1d(k) &
               + (c1 - weight) * (vals_1d_1klev - running_mean_array(index)%vals_1d(k))
         else
            running_mean_array(index)%vals_1d(k) = vals_1d_1klev
            running_mean_array(index)%linit_1d(k) = .true.
         endif
      else
         do k_loc=1,km
            if (running_mean_array(index)%linit_1d(k_loc)) then
               running_mean_array(index)%vals_1d(k_loc) = running_mean_array(index)%vals_1d(k_loc) &
                  + (c1 - weight) * (vals_1d_klevs(k_loc) - running_mean_array(index)%vals_1d(k_loc))
            else
               running_mean_array(index)%vals_1d(k_loc) = vals_1d_klevs(k_loc)
               running_mean_array(index)%linit_1d(k_loc) = .true.
            endif
         end do
      endif
   case(2)
      if (present(vals_2d_1block)) then
         if (running_mean_array(index)%linit_2d(block)) then
            running_mean_array(index)%vals_2d(:,:,block) = running_mean_array(index)%vals_2d(:,:,block) &
               + (c1 - weight) * (vals_2d_1block(:,:) - running_mean_array(index)%vals_2d(:,:,block))
         else
            running_mean_array(index)%vals_2d(:,:,block) = vals_2d_1block(:,:)
            running_mean_array(index)%linit_2d(block) = .true.
         endif
      else
         do iblock=1,nblocks_clinic
            if (running_mean_array(index)%linit_2d(iblock)) then
               running_mean_array(index)%vals_2d(:,:,iblock) = running_mean_array(index)%vals_2d(:,:,iblock) &
                  + (c1 - weight) * (vals_2d_blocks(:,:,iblock) - running_mean_array(index)%vals_2d(:,:,iblock))
            else
               running_mean_array(index)%vals_2d(:,:,iblock) = vals_2d_blocks(:,:,iblock)
               running_mean_array(index)%linit_2d(iblock) = .true.
            endif
         end do
      endif
   case(3)
      if (present(vals_3d_1klev_1block)) then
         if (running_mean_array(index)%linit_3d(k,block)) then
            running_mean_array(index)%vals_3d(:,:,k,block) = running_mean_array(index)%vals_3d(:,:,k,block) &
               + (c1 - weight) * (vals_3d_1klev_1block(:,:) - running_mean_array(index)%vals_3d(:,:,k,block))
         else
            running_mean_array(index)%vals_3d(:,:,k,block) = vals_3d_1klev_1block(:,:)
            running_mean_array(index)%linit_3d(k,block) = .true.
         endif
      else
         do iblock=1,nblocks_clinic
            if (running_mean_array(index)%linit_3d(k,iblock)) then
               running_mean_array(index)%vals_3d(:,:,k,iblock) = running_mean_array(index)%vals_3d(:,:,k,iblock) &
                  + (c1 - weight) * (vals_3d_1klev_blocks(:,:,iblock) - running_mean_array(index)%vals_3d(:,:,k,iblock))
            else
               running_mean_array(index)%vals_3d(:,:,k,iblock) = vals_3d_1klev_blocks(:,:,iblock)
               running_mean_array(index)%linit_3d(k,iblock) = .true.
            endif
         end do
      endif
   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_update_var

!***********************************************************************
!BOP
! !IROUTINE: running_mean_get_var
! !INTERFACE:

   subroutine running_mean_get_var(index, block, k, vals_0d, &
      vals_1d_1klev, vals_1d_klevs, vals_2d_1block, vals_2d_blocks, &
      vals_3d_1klev_1block, vals_3d_1klev_blocks)

! !DESCRIPTION:
!  Update a running mean variable using specified values.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in)                       :: index
   integer (int_kind), intent(in), optional             :: block
   integer (int_kind), intent(in), optional             :: k

! !OUTPUT PARAMETERS:

   real (r8), intent(out), optional                     :: vals_0d
   real (r8), intent(out), optional                     :: vals_1d_1klev
   real (r8), intent(out), dimension(:), optional       :: vals_1d_klevs
   real (r8), intent(out), dimension(:,:), optional     :: vals_2d_1block
   real (r8), intent(out), dimension(:,:,:), optional   :: vals_2d_blocks
   real (r8), intent(out), dimension(:,:), optional     :: vals_3d_1klev_1block
   real (r8), intent(out), dimension(:,:,:), optional   :: vals_3d_1klev_blocks

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock ! loop index
   integer (int_kind) :: k_loc  ! loop index

!-----------------------------------------------------------------------
!  error checking
!-----------------------------------------------------------------------

   if (index < 1 .or. index > running_mean_cnt) then
      call document('running_mean_get_var', 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
   endif

   select case (running_mean_array(index)%rank)
   case(0)
      if (.not. present(vals_0d)) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_0d must be supplied for 0d vars')
      endif
      if (.not. running_mean_array(index)%linit_0d) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
      endif
   case(1)
      if (.not. present(vals_1d_1klev) .and. .not. present(vals_1d_klevs)) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_1d_1klev or vals_1d_klevs must be supplied for 1d vars')
      endif
      if (present(vals_1d_1klev)) then
         if (.not. present(k)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call exit_POP(sigAbort, 'k must be supplied if vals_1d_1klev is supplied')
         endif
         if (.not. running_mean_array(index)%linit_1d(k)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call document('running_mean_get_var', 'k', k)
            call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
         endif
      else
         do k_loc=1,km
            if (.not. running_mean_array(index)%linit_1d(k_loc)) then
               call document('running_mean_get_var', 'name', running_mean_array(index)%name)
               call document('running_mean_get_var', 'k_loc', k_loc)
               call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
            endif
         end do
      endif
   case(2)
      if (.not. present(vals_2d_1block) .and. .not. present(vals_2d_blocks)) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_2d_1block or vals_2d_blocks must be supplied for 2d vars')
      endif
      if (present(vals_2d_1block)) then
         if (.not. present(block)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call exit_POP(sigAbort, 'block must be supplied if vals_2d_1block is supplied')
         endif
         if (.not. running_mean_array(index)%linit_2d(block)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call document('running_mean_get_var', 'block', block)
            call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
         endif
      else
         do iblock=1,nblocks_clinic
            if (.not. running_mean_array(index)%linit_2d(iblock)) then
               call document('running_mean_get_var', 'name', running_mean_array(index)%name)
               call document('running_mean_get_var', 'iblock', iblock)
               call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
            endif
         end do
      endif
   case(3)
      if (.not. present(k)) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'k must be supplied for 3d vars')
      endif
      if (.not. present(vals_3d_1klev_1block) .and. .not. present(vals_3d_1klev_blocks)) then
         call document('running_mean_get_var', 'name', running_mean_array(index)%name)
         call exit_POP(sigAbort, 'vals_3d_1klev_1block or vals_3d_1klev_blocks must be supplied for 3d vars')
      endif
      if (present(vals_3d_1klev_1block)) then
         if (.not. present(block)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call exit_POP(sigAbort, 'block must be supplied if vals_3d_1klev_1block is supplied')
         endif
         if (.not. running_mean_array(index)%linit_3d(k,block)) then
            call document('running_mean_get_var', 'name', running_mean_array(index)%name)
            call document('running_mean_get_var', 'k', k)
            call document('running_mean_get_var', 'block', block)
            call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
         endif
      else
         do iblock=1,nblocks_clinic
            if (.not. running_mean_array(index)%linit_3d(k,iblock)) then
               call document('running_mean_get_var', 'name', running_mean_array(index)%name)
               call document('running_mean_get_var', 'k', k)
               call document('running_mean_get_var', 'iblock', iblock)
               call exit_POP(sigAbort, 'running_mean_get_var must not be called on an uninitialized running mean')
            endif
         end do
      endif
   end select

!-----------------------------------------------------------------------
!  get values
!-----------------------------------------------------------------------

   select case (running_mean_array(index)%rank)
   case(0)
      vals_0d = running_mean_array(index)%vals_0d
   case(1)
      if (present(vals_1d_1klev)) then
         vals_1d_1klev = running_mean_array(index)%vals_1d(k)
      else
         vals_1d_klevs(:) = running_mean_array(index)%vals_1d(:)
      endif
   case(2)
      if (present(vals_2d_1block)) then
         vals_2d_1block(:,:) = running_mean_array(index)%vals_2d(:,:,block)
      else
         do iblock=1,nblocks_clinic
            vals_2d_blocks(:,:,iblock) = running_mean_array(index)%vals_2d(:,:,iblock)
         end do
      endif
   case(3)
      if (present(vals_3d_1klev_1block)) then
         vals_3d_1klev_1block(:,:) = running_mean_array(index)%vals_3d(:,:,k,block)
      else
         do iblock=1,nblocks_clinic
            vals_3d_1klev_blocks(:,:,iblock) = running_mean_array(index)%vals_3d(:,:,k,iblock)
         end do
      endif
   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_get_var

!***********************************************************************
!BOP
! !IROUTINE: running_mean_write_restart
! !INTERFACE:

   subroutine running_mean_write_restart(restart_file, action)

! !DESCRIPTION:
!  Write necessary data to the restart file
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: index
   type (io_dim)        :: i_dim        ! dimension descriptor
   type (io_dim)        :: j_dim        ! dimension descriptor
   type (io_dim)        :: k_dim        ! dimension descriptor
   type (io_field_desc), dimension(:), allocatable, save &
                        :: vals_desc    ! io field descriptor

!-----------------------------------------------------------------------

   if (trim(action) == 'define') then
      allocate(vals_desc(running_mean_cnt))

      do index=1,running_mean_cnt
         select case (running_mean_array(index)%rank)
         case(0)
            vals_desc(index) = &
               construct_io_field(running_mean_array(index)%file_varname, &
                                  d0d_array=running_mean_array(index)%vals_0d)
         case(1)
            k_dim = construct_io_dim('k', km)
            vals_desc(index) = &
               construct_io_field(running_mean_array(index)%file_varname, dim1=k_dim, &
                                  d1d_array=running_mean_array(index)%vals_1d)
         case(2)
            i_dim = construct_io_dim('i', nx_global)
            j_dim = construct_io_dim('j', ny_global)
            vals_desc(index) = &
               construct_io_field(running_mean_array(index)%file_varname, dim1=i_dim, dim2=j_dim, &
                                  d2d_array=running_mean_array(index)%vals_2d)
         case(3)
            i_dim = construct_io_dim('i', nx_global)
            j_dim = construct_io_dim('j', ny_global)
            k_dim = construct_io_dim('k', km)
            vals_desc(index) = &
               construct_io_field(running_mean_array(index)%file_varname, dim1=i_dim, dim2=j_dim, dim3=k_dim, &
                                  d3d_array=running_mean_array(index)%vals_3d)
         end select

         call data_set(restart_file, 'define', vals_desc(index))
      end do
   endif

   if (trim(action) == 'write') then
      do index=1,running_mean_cnt
         call data_set(restart_file, 'write', vals_desc(index))
         call destroy_io_field(vals_desc(index))
      end do

      deallocate(vals_desc)
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine running_mean_write_restart

end module running_mean_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
