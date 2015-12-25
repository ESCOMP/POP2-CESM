!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module IRF_mod

!BOP
! !MODULE: IRF_mod
!
! Impulse response functions (IRFs) are diagnostics tracers in a
! GCM. The tracers are set to specified patterns (e.g. point impulses
! localized in space) and the GCM computes tendencies (e.g. advective and
! mixing fluxes) based on these patterns. The tracers are reset to their
! respective patterns at every model timestep. The time-mean of the computed
! tendencies are the output of interest. This diagnostic output provides a
! representation of the GCM's circulation operators (advection and mixing).
! Two target applications of IRFs in CESM-POP are to provide inputs for
! offline tracer transport models and to provide a representation of model
! circulation for the preconditioner in the Newton-Krylov based solver
! for fast tracer spinup.
!
! IRF_NT is set as a key-value pair in CESM's OCN_TRACER_MODULE_OPT in env_build.xml
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use blocks, only: nx_block, ny_block
   use domain_size, only: km
   use kinds_mod

   implicit none

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: IRF_tracer_cnt, &
             IRF_init,       &
             IRF_reset

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: IRF_tracer_cnt = IRF_NT

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: IRF_init
! !INTERFACE:

 subroutine IRF_init(tracer_d_module, TRACER_MODULE)

! !DESCRIPTION:
!  Initialize IRF tracer module. This involves setting metadata, reading
!  the module's namelist and setting initial conditions.
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast
   use constants, only: c0, blank_fmt, delim_fmt, ndelim_fmt
   use communicate, only: my_task, master_task
   use domain, only: nblocks_clinic
   use domain_size, only: max_blocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use io_tools, only: document
   use io_types, only: stdout, nml_in, nml_filename
   use passive_tracer_tools, only: read_field
   use prognostic, only: oldtime, curtime, tracer_field
   use grid, only: KMT

   use netcdf

! !INPUT/OUTPUT PARAMETERS:

!   type (tracer_field), dimension(IRF_tracer_cnt), intent(inout) :: &
   type (tracer_field), dimension(:), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

!   real(r8), dimension(nx_block,ny_block,km,IRF_tracer_cnt,3,max_blocks_clinic), &
   real(r8), dimension(:,:,:,:,:,:), &
      intent(inout) :: TRACER_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'IRF_mod:IRF_init'

   character(char_len) :: &
      irf_tracer_file              ! filename for tracer impulse patterns

   integer(int_kind) :: &
      irf_tracer_file_ind_start, & ! starting index for tracers from irf_tracer_file
      stat,                      & ! status of netCDF call
      ncid,                      & ! netCDF file id
      varid,                     & ! netCDF variable id
      irf_tracer_file_nchar,     & ! character length of variable names in irf_tracer_file
      char_pos,                  & ! character position in a string
      n,                         & ! index for looping over tracers
      iblock,                    & ! index for looping over blocks
      k,                         & ! index for looping over depth levels
      nml_error                    ! namelist i/o error flag

   integer (int_kind), dimension(2) :: &
      varnames_dimids              ! dimids for tracer names variable

   namelist /irf_nml/ &
      irf_tracer_file, irf_tracer_file_ind_start

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   irf_tracer_file = 'unknown'
   irf_tracer_file_ind_start = 1

!-----------------------------------------------------------------------
!  read namelist input
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=irf_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'irf_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' IRF:'
      write(stdout,blank_fmt)
      write(stdout,*) ' IRF namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,irf_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(irf_tracer_file, master_task)
   call broadcast_scalar(irf_tracer_file_ind_start, master_task)

!-----------------------------------------------------------------------
!  read in tracer names on master_task
!  jump out of master_task conditional if an error is encountered
!  broadcast values that were read in
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      stat = nf90_open(irf_tracer_file, 0, ncid)
      if (stat /= 0) then
         write(stdout,*) 'error from nf90_open: ', nf90_strerror(stat)
         go to 99
      endif

      stat = nf90_inq_varid(ncid, 'var_names', varid)
      if (stat /= 0) then
         write(stdout,*) 'error from nf90_inq_varid for var_names: ', &
                         nf90_strerror(stat)
         go to 99
      endif

      stat = nf90_inquire_variable(ncid, varid, dimids=varnames_dimids)
      if (stat /= 0) then
         write(stdout,*) 'error from nf90_inquire_variable for dimids for var_names: ', &
                         nf90_strerror(stat)
         go to 99
      endif

      stat = nf90_inquire_dimension(ncid, varnames_dimids(1), len=irf_tracer_file_nchar)
      if (stat /= 0) then
         write(stdout,*) 'error from nf90_inquire_dimension for len of 1st dim of var_names: ', &
                         nf90_strerror(stat)
         go to 99
      endif

      do n = 1, IRF_tracer_cnt
         tracer_d_module(n)%short_name = ''
         stat = nf90_get_var(ncid, varid, tracer_d_module(n)%short_name, &
                             start=(/ 1, n+irf_tracer_file_ind_start-1 /), &
                             count=(/ irf_tracer_file_nchar, 1 /))
         if (stat /= 0) then
            write(stdout,*) 'error from nf90_get_var for varname ', n, ': ', &
                            nf90_strerror(stat)
            go to 99
         endif

         ! replace trailing null characters with spaces
         char_pos = index(tracer_d_module(n)%short_name, char(0))
         if ( char_pos > 0 ) then
            tracer_d_module(n)%short_name(char_pos:char_len) = ' '
         endif
      end do

      stat = nf90_close(ncid)
      if (stat /= 0) then
         write(stdout,*) 'error from nf90_close: ', nf90_strerror(stat)
         go to 99
      endif

   endif

99 call broadcast_scalar(stat, master_task)
   if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                          &/ subname)

   do n = 1, IRF_tracer_cnt
      call broadcast_scalar(tracer_d_module(n)%short_name, master_task)
   end do

!-----------------------------------------------------------------------
!  initialize other tracer_d values
!-----------------------------------------------------------------------

   do n = 1, IRF_tracer_cnt
      tracer_d_module(n)%long_name  = tracer_d_module(n)%short_name
      tracer_d_module(n)%units      = '1'
      tracer_d_module(n)%tend_units = '1/s'
      tracer_d_module(n)%flux_units = '1/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  initialize tracers
!
!  Because the IRF tracers are reset to their initial state at the end
!  of every timestep, they do not really evolve in time. So it is never
!  necessary to read them from the model restart file. So the code is
!  hardwired to always initialize the tracers from the irf_tracer_file.
!
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do n = 1, IRF_tracer_cnt
      call read_field('nc', irf_tracer_file, &
                      tracer_d_module(n)%short_name, &
                      TRACER_MODULE(:,:,:,n,oldtime,:))
   end do

   !$OMP PARALLEL DO PRIVATE(iblock,n,k)
   do iblock=1,nblocks_clinic
      do n = 1,IRF_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
            TRACER_MODULE(:,:,k,n,curtime,iblock) = &
               TRACER_MODULE(:,:,k,n,oldtime,iblock)
         end do
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine IRF_init

!***********************************************************************
!BOP
! !IROUTINE: IRF_reset(TRACER_MODULE_OLD, TRACER_MODULE_NEW)
! !INTERFACE:

 subroutine IRF_reset(TRACER_MODULE_OLD, TRACER_MODULE_NEW)

! !DESCRIPTION:
!  reset IRF newtime values to oldtime values
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

!   real(r8), dimension(nx_block,ny_block,km,IRF_tracer_cnt), intent(in) :: &
   real(r8), dimension(:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD  ! IRF tracers at oldtime

! !INPUT/OUTPUT PARAMETERS:

!   real(r8), dimension(nx_block,ny_block,km,IRF_tracer_cnt), intent(inout) :: &
   real(r8), dimension(:,:,:,:), intent(inout) :: &
      TRACER_MODULE_NEW  ! IRF tracers at new

!EOP
!BOC
!-----------------------------------------------------------------------

   TRACER_MODULE_NEW = TRACER_MODULE_OLD

!-----------------------------------------------------------------------
!EOC

 end subroutine IRF_reset

!***********************************************************************

end module IRF_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
