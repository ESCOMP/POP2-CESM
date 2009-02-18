!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module tavg

!BOP
! !MODULE: tavg
! !DESCRIPTION:
!  This module contains data types and routines for computing running 
!  time-averages of selected fields and writing this data to files.
!
! !REVISION HISTORY:
!  SVN:$Id$
!  

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use prognostic
   use grid
   use time_management
   use global_reductions
   use broadcast
   use io
   use io_types
   use exit_mod
 
   !*** ccsm
   use gather_scatter
   use operators, only: zcurl
   use io_ccsm
   use io_tools
   use diag_bsf, only: pcg_diag_bsf_solver, init_diag_bsf
   use diags_on_lat_aux_grid
   use registry
   use timers
#ifdef CCSMCOUPLED
   use shr_sys_mod
#endif

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_tavg,              &
             define_tavg_field,      &
             tavg_increment_sum_qflux,&
             accumulate_tavg_field,  &
             tavg_requested,         &
             set_in_tavg_contents,   &
             write_tavg,             &
             read_tavg,              &
             tavg_set_flag
 
   !*** ccsm
   public :: tavg_id,                &
             tavg_global_sum_2D


! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      ltavg_on      = .false., & ! tavg file output wanted
      ltavg_restart = .false.    ! run started from restart

   integer (i4), parameter, public ::  &
      tavg_method_unknown = 0,         &
      tavg_method_avg     = 1,         &
      tavg_method_min     = 2,         &
      tavg_method_max     = 3,         &
      tavg_method_qflux   = 4

   real (r8),public ::  &
      tavg_sum           ! accumulated time (in seconds)

   !*** ccsm
   real (r8),public ::  &
      tavg_sum_qflux


!-----------------------------------------------------------------------
!
!  tavg field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

! !PUBLIC TYPES:
   !*** ccsm
   type,public :: tavg_field_desc_ccsm
      character(char_len)     :: short_name     ! short name for field
      character(char_len)     :: long_name      ! long descriptive name
      character(char_len)     :: units          ! units
      character(char_len)     :: coordinates    ! coordinates
      character(char_len)     :: nftype         ! indicates data type 
      character(4)            :: grid_loc       ! location in grid
      real (rtavg)            :: fill_value     ! _FillValue
      real (rtavg)            :: missing_value  ! value on land pts
      real (rtavg)            :: scale_factor   ! r4 scale factor
      real (r4), dimension(2) :: valid_range    ! min/max
      integer (i4)            :: ndims          ! num dims (2 or 3)
      integer (i4)            :: buf_loc        ! location in buffer
      integer (i4)            :: method         ! method for averaging
      integer (i4)            :: field_loc      ! grid location and field
      integer (i4)            :: field_type     ! type for io, ghost cells
   end type

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  tavg field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

   type :: tavg_field_desc
      character(char_len)     :: short_name     ! short name for field
      character(char_len)     :: long_name      ! long descriptive name
      character(char_len)     :: units          ! units
      character(4)            :: grid_loc       ! location in grid
      real (r4)               :: missing_value  ! value on land pts
      real (r4), dimension(2) :: valid_range    ! min/max
      integer (i4)            :: ndims          ! num dims (2 or 3)
      integer (i4)            :: buf_loc        ! location in buffer
      integer (i4)            :: method         ! method for averaging
      integer (i4)            :: field_loc      ! grid location and field
      integer (i4)            :: field_type     !  type for io, ghost cells
   end type

   integer (int_kind), parameter :: &
      max_avail_tavg_fields = 600    ! limit on available fields - can
                                     !   be pushed as high as necessary

   integer (int_kind) ::                &
      num_avail_tavg_fields      = 0,   &! current number of defined fields
      num_requested_tavg_fields,        &! number of fields requested
      tavg_flag                          ! time flag for writing tavg files

   !*** ccsm
   type (tavg_field_desc_ccsm), dimension(max_avail_tavg_fields) :: &
      avail_tavg_fields

   type (io_dim) ::   &
      i_dim, j_dim,   &! dimension descriptors for horiz dims
      k_dim,          &! dimension descriptor for vert levels (z_t or z_w grid)
      time_dim         ! dimension descriptor for (unlimited) time dim
 

!-----------------------------------------------------------------------
!
!  buffers for holding running tavg variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_bufsize_2d,   &    ! size of buffer for 2d fields
      tavg_bufsize_3d         ! size of buffer for 3d fields

   real (rtavg), dimension(:,:,:,:), allocatable :: &
      TAVG_BUF_2D         ! buffer for holding accumulated sums

   real (rtavg), dimension(:,:,:,:,:), allocatable :: &
      TAVG_BUF_3D         ! buffer for holding accumulated sums

   integer (i4), dimension(:), allocatable :: &
      TAVG_BUF_2D_METHOD,  &! method for each requested 2d field
      TAVG_BUF_3D_METHOD    ! method for each requested 3d field

   real (rtavg), dimension (:,:,:), allocatable ::  &
         TAVG_TEMP          ! work array in write_restart
!-----------------------------------------------------------------------
!
!  variables for writing data
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      tavg_freq_iopt,  &! frequency option for writing tavg
      tavg_freq,       &! frequency of tavg output
      tavg_start_iopt, &! start after option
      tavg_start        ! start tavg after tavg_start

   character (char_len) ::  &
      tavg_infile,          & ! filename for restart input
      tavg_outfile,         & ! root filename for tavg output
      tavg_outfile_orig       ! root filename for tavg output (original)

   character (char_len) ::    &
      tavg_fmt_in,            & ! format (nc or bin) for reading
      tavg_fmt_out              ! format (nc or bin) for writing

   type (datafile) :: tavg_file_desc    ! IO file descriptor
 

!-----------------------------------------------------------------------
!
!  scalars
!
!-----------------------------------------------------------------------

   real (r8) ::        &
      dtavg             ! current time step

   character (10) :: &
      beg_date       ! date on which the current accumulated sum
                     ! was started (not the tavg_start date)

!-----------------------------------------------------------------------
!
!  coupled code -- ccsm-related
!
!-----------------------------------------------------------------------
   logical (log_kind) :: &
      lccsm,             &
      ldiag_bsf,         &
      ltavg_fmt_out_nc, &! true if netCDF output format
      implied_time_dim

   logical (log_kind), dimension (:,:,:,:), allocatable ::  &
      MASK_22
 
   integer (int_kind), parameter ::      &
      max_avail_tavg_nstd_fields =  20,  &! limit on available fields
      max_avail_tavg_labels      =  10,  &
      max_num_ccsm_coordinates   =  10,  &
      max_num_ccsm_time_invar    =  50,  &
      max_num_ccsm_scalars       = 100
 
   integer (int_kind) ::     &
      num_avail_tavg_nstd_fields   = 0, &! current number of defined nonstandard fields
      num_avail_tavg_labels        = 0, &! current number of ccsm labels
      num_ccsm_coordinates         = 0, &
      num_ccsm_time_invar          = 0, &
      num_ccsm_scalars             = 0, &
      zt_150m_levs

   real (r4), dimension(:), allocatable, target :: &
      ZT_150m_R4   ! single precision array
 
   integer (int_kind) ::  &
      tavg_BSF,           &
      tavg_MOC,           &
      tavg_N_HEAT,        &
      tavg_N_SALT

   integer (int_kind), dimension(:,:), allocatable    ::  &
      KMU_G            ! k index of deepest grid cell on global U grid

   integer (i4) ::  &
      time_bound_id,&
      moc_id,       &
      n_heat_id,    &
      n_salt_id

   integer (i4), dimension(max_avail_tavg_labels) ::  &
      avail_tavg_labels_id

   type (tavg_field_desc_ccsm) ::  &
      avail_tavg_nstd_fields (max_avail_tavg_nstd_fields),  &
      avail_tavg_labels      (max_avail_tavg_labels)

   type (io_field_desc) ::                          &
      ccsm_coordinates (max_num_ccsm_coordinates),  &
      ccsm_time_invar  (max_num_ccsm_time_invar),   &
      ccsm_scalars     (max_num_ccsm_scalars),      &
      time_coordinate  (1)

   type (io_dim) ::       &
      z_dim,              &! dimension descriptor for vert (z_t or z_w grid)
      zt_dim,             &! dimension descriptor for vert (z_t grid)
      zt_150m_dim,        &! dimension descriptor for near-surf vert (z_t grid)
      zw_dim,             &! dimension descriptor for vert (z_w grid)
      tr_dim,             &! dimension descriptor 
      nchar_dim,          &! dimension descriptor for character arrays
      d2_dim,             &! dimension descriptor  
      lat_aux_grid_dim,   &! dimension descriptor  
      moc_z_dim,          &! dimension descriptor  
      moc_comp_dim,       &! dimension descriptor  
      transport_comp_dim, &! dimension descriptor  
      transport_reg_dim    ! dimension descriptor  
 
   type(io_dim)                                          ::  &
      time_bound_dims   (2),                                 &
      io_dims_nstd_ccsm (5,max_avail_tavg_nstd_fields),      &
      io_dims_labels    (2,max_avail_tavg_labels)

   integer, dimension (max_avail_tavg_nstd_fields) ::  &
      ndims_nstd_ccsm
 
   real (r8) ::        &
     lower_time_bound, &! lower time bound for time_bounds variable
     upper_time_bound

   integer (int_kind) ::  &
      tavg_debug  = 0      ! debug level [0,1]  1 ==> messages
 
!-----------------------------------------------------------------------
!
!     variables for local spatial averaging of some time-averaged fields
!
!-----------------------------------------------------------------------

   logical (log_kind) ::   &
      ltavg_nino_diags    ! true if requesting nino diagnostics and if
                          ! other requirements are met

   integer (int_kind), parameter ::  &
      n_reg_0D = 4                    ! number of regions

   real (r8), dimension(:), allocatable :: &
      SAVG_0D,                             &! local- and time-mean value 
      SAVG_0D_AREA                          ! area of the region

   real (r8), dimension(:,:,:,:), allocatable ::  &
      SAVG_0D_MASK                          ! mask for the region, i.e. 0 and 1 mean
                                            ! outside and inside the region, respectively

   character (char_len), dimension(:), allocatable ::  &
      SAVG_0D_NAME                          ! name of the region

   integer (int_kind) :: &
      timer_write_std,   &
      timer_write_nstd,  &
      timer_tavg_ccsm_diags_bsf, &
      timer_tavg_ccsm_diags_moc, &
      timer_tavg_ccsm_diags_trans

!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_tavg
! !INTERFACE:

 subroutine init_tavg

! !DESCRIPTION:
!  This routine initializes tavg options and reads in contents file to
!  determine which fields for which the user wants time-averaged data.
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
   save

   integer (POP_i4) ::     &
      errorCode

   integer (i4) ::         &
      n,                   &! dummy index
      i,ip1,j,k,           &! dummy indices
      iblock,              &! local block index
      loc,                 &! location of field in buffer
      nu,                  &! unit for contents input file
      cindex,              &! character index for manipulating strings
      nml_error,           &! namelist i/o error flag
      contents_error        ! error flag for contents file read

   integer (int_kind) ::   &
      max_days            ! maximum number of days per month in a year

   character (char_len) :: &
      tavg_freq_opt,       &! choice for frequency of tavg output
      tavg_start_opt,      &! choice for starting averaging
      tavg_contents,       &! filename for choosing fields for output
      char_temp             ! temporary for manipulating fields

   character (33), parameter :: &
      freq_fmt = "('tavg diagnostics every ',i6,a8)"

   character (44), parameter :: &
      start_fmt = "('tavg sums accumulated starting at ',a5,i8)"

   type (block) ::        &
      this_block          ! block information for current block

   namelist /tavg_nml/ tavg_freq_opt, tavg_freq, tavg_infile,       &
                       tavg_outfile, tavg_contents, tavg_start_opt, &
                       tavg_start, tavg_fmt_in, tavg_fmt_out,       &
                       ltavg_nino_diags


!-----------------------------------------------------------------------
!
!  determine if this is a ccsm coupled run
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   lccsm = registry_match('lcoupled')

!-----------------------------------------------------------------------
!
!  read tavg file output frequency and filenames from namelist
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a12)') ' Tavg:'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout)
   endif

   tavg_freq_iopt = freq_opt_never
   tavg_freq      = 100000
   tavg_start_iopt = start_opt_nstep
   tavg_start      = 0
   tavg_infile    = 'unknown_tavg_infile'
   tavg_outfile   = 't'
   tavg_contents  = 'unknown_tavg_contents'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=tavg_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading tavg_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a28)') ' tavg_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,tavg_nml)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout)
   endif

   if (my_task == master_task) then
      select case (tavg_freq_opt)
      case ('never')
         tavg_freq_iopt = freq_opt_never
         write(stdout,'(a20)') 'tavg diagnostics off'
      case ('nyear')
         tavg_freq_iopt = freq_opt_nyear
         write(stdout,freq_fmt) tavg_freq,' years  '
      case ('nmonth')
         tavg_freq_iopt = freq_opt_nmonth
         write(stdout,freq_fmt) tavg_freq,' months '
      case ('nday')
         tavg_freq_iopt = freq_opt_nday
         write(stdout,freq_fmt) tavg_freq,' days   '
      case ('nhour')
         tavg_freq_iopt = freq_opt_nhour
         write(stdout,freq_fmt) tavg_freq,' hours  '
      case ('nsecond')
         tavg_freq_iopt = freq_opt_nsecond
         write(stdout,freq_fmt) tavg_freq,' seconds'
      case ('nstep')
         tavg_freq_iopt = freq_opt_nstep
         write(stdout,freq_fmt) tavg_freq,' steps  '
      case default
         tavg_freq_iopt = -1000
      end select

      if (tavg_freq_iopt /= freq_opt_never) then
         select case (tavg_start_opt)
         case ('nstep')
            tavg_start_iopt = start_opt_nstep
            write(stdout,start_fmt) 'step ', tavg_start
         case ('nday')
            tavg_start_iopt = start_opt_nday
            write(stdout,start_fmt) 'day  ', tavg_start
         case ('nyear')
            tavg_start_iopt = start_opt_nyear
            write(stdout,start_fmt) 'year ', tavg_start
         case ('date')
            tavg_start_iopt = start_opt_date
            write(stdout,start_fmt) '     ', tavg_start
         case default
            tavg_start_iopt = -1000
         end select
      endif

   endif

   call POP_IOUnitsFlush(POP_stdout)

   call broadcast_scalar(tavg_freq_iopt, master_task)

   if (tavg_freq_iopt == -1000) then
      call exit_POP(sigAbort,'unknown option for tavg file frequency')
   else if (tavg_freq_iopt /= freq_opt_never) then
      call broadcast_scalar(tavg_freq,         master_task)
      call broadcast_scalar(tavg_start_iopt,   master_task)
      call broadcast_scalar(tavg_start,        master_task)
      call broadcast_scalar(tavg_infile,       master_task)
      call broadcast_scalar(tavg_outfile,      master_task)
      call broadcast_scalar(tavg_contents,     master_task)
      call broadcast_scalar(tavg_fmt_in,       master_task)
      call broadcast_scalar(tavg_fmt_out,      master_task)
      call broadcast_scalar(ltavg_nino_diags,  master_task)

      if (tavg_start_iopt == -1000) then
         call exit_POP(sigAbort,'unknown option for tavg start option')
      endif

   endif

   if (trim(tavg_fmt_out) == 'nc') then
      ltavg_fmt_out_nc = .true.
   else
      ltavg_fmt_out_nc = .false.
   endif

   tavg_outfile_orig = char_blank
   tavg_outfile_orig = trim(tavg_outfile)

!-----------------------------------------------------------------------
!
!  define time-averaged BSF field; it may or may not be requested in 
!  the contents file
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_BSF,'BSF',2,                                 &
                          long_name='Diagnostic barotropic streamfunction', &
                          units='Sv', grid_loc='2220',                      &
                          coordinates='ULONG ULAT time')


!-----------------------------------------------------------------------
!
!  initialize time flag for writing tavg files
!
!-----------------------------------------------------------------------

   tavg_flag = init_time_flag('tavg',default=.false.,    &
                              freq_opt = tavg_freq_iopt, &
                              freq     = tavg_freq)

!-----------------------------------------------------------------------
!
!  read contents file to determine which fields to dump
!
!-----------------------------------------------------------------------

   if (tavg_freq_iopt /= freq_opt_never) then

      tavg_bufsize_2d = 0
      tavg_bufsize_3d = 0

      !*** count the number of lines in tavg_contents file
      call tavg_count_contents_ccsm(num_requested_tavg_fields,tavg_contents)

      call get_unit(nu)
      if (my_task == master_task) then
         open(nu, file=tavg_contents, status='old')
         write(stdout,'(a38)') 'tavg diagnostics requested for fields:'
         call POP_IOUnitsFlush(POP_stdout)
      endif

      call broadcast_scalar(num_requested_tavg_fields, master_task)

      contents_error = 0

      do n=1,num_requested_tavg_fields

         if (my_task == master_task) then
            read(nu,'(a)',iostat=contents_error) char_temp
            char_temp = adjustl(char_temp)
            cindex = index(char_temp,' ')
            char_temp(cindex:) = ' '
            write(stdout,*) '  ',trim(char_temp)
            call POP_IOUnitsFlush(POP_stdout)
         endif

         call broadcast_scalar(contents_error, master_task)
         if (contents_error /= 0) then
            call exit_POP(sigAbort,'error reading tavg contents')
         endif

         call broadcast_scalar(char_temp, master_task)

         !*** activate requested tavg fields
         !    determines tavg_bufsize_2d, tavg_bufsize_3d
         call request_tavg_field(trim(char_temp))  
      end do

      call release_unit(nu)

      !*** allocate and initialize running tavg buffers

      allocate(                                                            &
         TAVG_BUF_2D(nx_block,ny_block,   nblocks_clinic,tavg_bufsize_2d), &
         TAVG_BUF_3D(nx_block,ny_block,km,nblocks_clinic,tavg_bufsize_3d), &
         TAVG_BUF_2D_METHOD(tavg_bufsize_2d),                              &
         TAVG_BUF_3D_METHOD(tavg_bufsize_3d))

      allocate(TAVG_TEMP(nx_block,ny_block,nblocks_clinic))

      tavg_sum = c0

      call time_stamp('now','ymd',date_string=beg_date)
      if (tavg_freq_iopt == freq_opt_nstep) &
         write(beg_date,'(i10)') nsteps_total

      do n = 1,num_avail_tavg_fields  ! check all available fields
         loc = abs(avail_tavg_fields(n)%buf_loc)
         if (loc /= 0) then  ! field is actually requested and in buffer
            if (avail_tavg_fields(n)%ndims == 2) then
               TAVG_BUF_2D_METHOD(loc) = avail_tavg_fields(n)%method
            else if (avail_tavg_fields(n)%ndims == 3) then
               TAVG_BUF_3D_METHOD(loc) = avail_tavg_fields(n)%method
            endif
         endif
      end do
 
      !*** initialize buffers based on requested method
      call tavg_reset_field_all

   endif !tavg_freq_iopt



   !*** define dimensions for tavg output files
   i_dim     = construct_io_dim('i',nx_global)
   j_dim     = construct_io_dim('j',ny_global)
   k_dim     = construct_io_dim('k',km)
   time_dim  = construct_io_dim('time',0) ! used only to set %active=.false.
 
!-----------------------------------------------------------------------
!
!  ccsm-specific initializations
!
!-----------------------------------------------------------------------
 
   if (lccsm) then
 
     !*** initialze barotropic stream function

     call init_diag_bsf(set_in_tavg_contents(tavg_BSF), errorCode)

     if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
           'init_tavg: error in init_diag_bsf')
        call exit_POP(sigAbort,'init_tavg: error in init_diag_bsf')
        return
     endif

     ldiag_bsf     = registry_match('ldiag_bsf')

     !*** initialze moc and heat/salt transport diagnostics
     call init_lat_aux_grid
     call init_moc_ts_transport_arrays

     !*** how many levels have their midpoint shallower than 150m
     zt_150m_levs = count(zt < 150.0e2_r8)
     if (.not. allocated(ZT_150m_R4)) &
       allocate(ZT_150m_R4(zt_150m_levs))

     !*** define dimensions for tavg output files
     i_dim      = construct_io_dim('nlon',nx_global)
     j_dim      = construct_io_dim('nlat',ny_global)
     zt_dim     = construct_io_dim('z_t',km)
     zt_150m_dim= construct_io_dim('z_t_150m',zt_150m_levs)
     zw_dim     = construct_io_dim('z_w',km)
     time_dim   = construct_io_dim('time',0)    ! "0" ==> unlimited dimension
     tr_dim     = construct_io_dim('tracers',nt)
     nchar_dim  = construct_io_dim('nchar',char_len)
     d2_dim     = construct_io_dim('d2',2)
     if (moc .or. n_heat_trans .or. n_salt_trans) then
       lat_aux_grid_dim  = construct_io_dim('lat_aux_grid',n_lat_aux_grid+1)
       transport_reg_dim = construct_io_dim('transport_reg',n_transport_reg)
       if (moc) then
         moc_z_dim    = construct_io_dim('moc_z',km+1)
         moc_comp_dim = construct_io_dim('moc_comp',n_moc_comp)
       endif
       if (n_heat_trans .or. n_salt_trans) &
         transport_comp_dim = construct_io_dim('transport_comp',n_transport_comp)
     endif

!-----------------------------------------------------------------------
!
!  set up masking arrays. allocate global arrays, use them, then deallocate
!
!-----------------------------------------------------------------------

     allocate (MASK_22(nx_block,ny_block,km,nblocks_clinic))

 
     !--------------------------------------------------------
     !    create MASK_22, layer by layer
     !--------------------------------------------------------
 
 
     do iblock = 1,nblocks_clinic 
     this_block = get_block(blocks_clinic(iblock),iblock)  
       do k=1,km
         MASK_22(:,:,k,iblock) = .false.
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
            if (k > KMT(i  ,j,   iblock)   .and.  &
                k > KMT(i+1,j,   iblock)   .and.  &
                k > KMT(i  ,j+1, iblock)   .and.  &
                k > KMT(i+1,j+1, iblock))   then
                  MASK_22(i,j,k,iblock) = .true.
            endif
         enddo ! i
         enddo ! j
       enddo ! k
     enddo ! iblock
 
   endif ! lccsm
 
  !*** check if local spatial averaging is possible based on tavg options
  !*** and namelist settings

   max_days = maxval(days_in_month) 

   if (set_in_tavg_contents(tavg_id('TEMP'))  ) then
    if ( (tavg_freq_iopt == freq_opt_nmonth  &
                .and. tavg_freq == 1)  .or.  &
         (tavg_freq_iopt == freq_opt_nday    &
                .and. tavg_freq <= max_days) .or.    &
         (tavg_freq_iopt == freq_opt_nhour           &
                .and. tavg_freq <= max_days*24) .or. &
         (tavg_freq_iopt == freq_opt_nsecond         &
                 .and. tavg_freq <= max_days*24*seconds_in_hour) .or. &
         (tavg_freq_iopt == freq_opt_nstep           &
                 .and. tavg_freq <= max_days*nsteps_per_day) ) then
         !*** ok to have ltavg_nino_diags enabled
    else
         !*** not ok to have ltavg_nino_diags enabled; disable it
         if (ltavg_nino_diags)  &
           call exit_POP(sigAbort,'init_tavg: nino diagnostics cannot be computed')
    endif

   else
         !*** not ok to have ltavg_nino_diags enabled; disable it
         if (ltavg_nino_diags)  call exit_POP(sigAbort,  &
           'init_tavg: nino diagnostics cannot be computed -- TEMP is not in contents file')
   endif ! tavg_id

   call tavg_init_local_spatial_avg
 

!-----------------------------------------------------------------------
!
!  finally, read restart file if necessary
!  must do this after i_dim, j_dim, etc are defined
!-----------------------------------------------------------------------

   !*** make sure tavg flag is set correctly
   call tavg_set_flag(flagonly=.true.)

   if (ltavg_on .and. ltavg_restart) then
      !*** do not read restart if last restart was at a tavg dump
      !*** interval (should start new tavg sums in this case)

      if (.not. time_to_do(tavg_freq_iopt, tavg_freq)) then
         call read_tavg
      endif
   endif

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------
   if (moc) then
     if (.not. set_in_tavg_contents(tavg_id('WVEL')) .or.  &
         .not. set_in_tavg_contents(tavg_id('VVEL')) )     &
        call exit_POP (SigAbort,                        &
       '(init_tavg) for moc diagnostics, WVEL and VVEL must be requested in tavg_contents file')
   endif

   if ( n_heat_trans .or. n_salt_trans) then
     if (.not. set_in_tavg_contents(tavg_id('ADVT'))   .or. &
         .not. set_in_tavg_contents(tavg_id('ADVS'))   .or. &
         .not. set_in_tavg_contents(tavg_id('VNT'))    .or. &
         .not. set_in_tavg_contents(tavg_id('VNS'))    .or. &
         .not. set_in_tavg_contents(tavg_id('HDIFT'))  .or. &
         .not. set_in_tavg_contents(tavg_id('HDIFS')) ) then
          call exit_POP (SigAbort, &
         '(init_tavg) for diag_gm_bolus, ADVT ADVS VNT VNS HDIFT HDIFS must all  be requested in tavg_contents file')
     endif
   endif

   if (n_heat_trans .or. n_salt_trans) then
   if (registry_match('diag_gm_bolus')) then
     if (.not. set_in_tavg_contents(tavg_id('ADVT_ISOP'))   .or. &
         .not. set_in_tavg_contents(tavg_id('ADVS_ISOP'))   .or. &
         .not. set_in_tavg_contents(tavg_id('VNT_ISOP'))    .or. &
         .not. set_in_tavg_contents(tavg_id('VNS_ISOP'))  ) then
          call exit_POP (SigAbort, &
         '(init_tavg) for diag_gm_bolus, ADVT_ISOP ADVS_ISOP VNT_ISOP VNS_ISOP must all  be requested in tavg_contents file')
     endif
   endif
   endif


!-----------------------------------------------------------------------
!
!  initialize timers
!
!-----------------------------------------------------------------------

   call get_timer(timer_write_std,'TAVG_WRITE_STD', nblocks_clinic, distrb_clinic%nprocs)
   call get_timer(timer_write_nstd,'TAVG_WRITE_NONSTD', nblocks_clinic, distrb_clinic%nprocs)
   if (ldiag_bsf)  &
   call get_timer(timer_tavg_ccsm_diags_bsf,'TAVG_CCSM_DIAGS_BSF', nblocks_clinic, distrb_clinic%nprocs)
   if (moc)  &
   call get_timer(timer_tavg_ccsm_diags_moc,'TAVG_CCSM_DIAGS_MOC', nblocks_clinic, distrb_clinic%nprocs)
   if (n_heat_trans .or. n_salt_trans)  &
   call get_timer(timer_tavg_ccsm_diags_trans,'TAVG_CCSM_DIAGS_TRANS', nblocks_clinic, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_tavg

!***********************************************************************
!BOP
! !IROUTINE: tavg_set_flag
! !INTERFACE:

 subroutine tavg_set_flag(flagonly)

! !DESCRIPTION:
!  This routine checks the time avg option and tavg start condition
!  to see whether tavg sums should be accumulated.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in), optional :: &
      flagonly        ! if true, only sets ltavg_on without advancing
                      !  time interval

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n ! loop index

   logical (log_kind) :: update_time

!-----------------------------------------------------------------------
!
!  if tavg requested and tavg not already turned on, check to see
!  if it is time to start time averaging
!
!-----------------------------------------------------------------------

   if (.not. ltavg_on .and. tavg_freq_iopt /= freq_opt_never) then

      ltavg_on = time_to_start(tavg_start_iopt, tavg_start)
      call time_stamp('now','ymd',date_string=beg_date)
      if (tavg_freq_iopt == freq_opt_nstep) &
         write(beg_date,'(i10)') nsteps_total

      !*** if it is time to start, make sure requested fields
      !*** get triggered by the requested function
 
      if (ltavg_on) then
         do n=1,num_avail_tavg_fields
            if (avail_tavg_fields(n)%buf_loc < 0) &
                avail_tavg_fields(n)%buf_loc =    &
                abs(avail_tavg_fields(n)%buf_loc)
         end do
         lower_time_bound = tday00
      endif
   endif

!-----------------------------------------------------------------------
!
!  setup time step and total integrated time for time average
!  adjust for averaging timesteps: if this is an averaging timestep,
!  the past values only contribute for 1/4 of a step and the
!  values for the step just before an averaging timestep contribute
!  for 1 1/4 steps.
!
!-----------------------------------------------------------------------

   update_time = .true.
   if (present(flagonly)) then
      if (flagonly) update_time = .false.
   endif

   if (ltavg_on .and. update_time) then
      if (avg_ts .or. back_to_back) then
         dtavg = p5*dtt
      else
         dtavg = dtt
      endif

      tavg_sum = tavg_sum + dtavg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_set_flag

!***********************************************************************
!BOP
! !IROUTINE: write_tavg
! !INTERFACE:

 subroutine write_tavg(restart_type)

! !DESCRIPTION:
!  This routine writes requested tavg fields to a file.  The fields are
!  normalized by the time interval before writing.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      restart_type           ! tells tavg whether to write restart

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
                       
   integer (int_kind) ::  &
      ndims

   integer (i4) ::  &
      nu,           &! i/o unit for output file
      iblock,       &! dummy block index
      nfield,       &! dummy field index
      nindex,       &! dummy field index
      loc,          &! buffer location for field
      io_phase,     &!'define' or 'write'
      k,n,          &! indices
      nn,           &! index
      i,j,          &! indices
      nstd_field_id

   character (char_len) ::  &
      string,               &! dummy character string
      file_suffix,          &! suffix to append to tavg file name
      hist_string,          &! string containing file history
      tavg_filename,        &! filename for tavg data
      tavg_pointer_file      ! filename for pointer file containing
                             !   location/name of last restart file

   character (8) ::  &
      date_created   ! string with (real) date this file created

   character (10) ::  &
      time_created   ! string with (real) date this file created

   logical (log_kind) ::  &
      ltavg_write,        &! time to write a file
      lreg_tavg_dump       ! time to reset time averages and write reg tavg file

   type (io_field_desc), dimension(:), allocatable ::  &
      tavg_fields

   type (io_field_desc)  ::  &
      qflux_field

   type (block) ::        &
      this_block          ! block information for current block

   num_avail_tavg_nstd_fields = 0

!-----------------------------------------------------------------------
!
!  is it time to write a file - if yes, create a file suffix
!
!-----------------------------------------------------------------------

   ltavg_write    = .false.
   lreg_tavg_dump = .false.

   if (ltavg_on) then
      ltavg_write = check_time_flag(tavg_flag)

      !*** regular tavg dump
      if (ltavg_write) then
         lreg_tavg_dump = .true.
         tavg_outfile = tavg_outfile_orig
         if (lccsm) then
           call tavg_create_suffix_ccsm(file_suffix)
         else
           call tavg_create_suffix(file_suffix)
         endif
      endif

      !*** tavg restart
      if (trim(restart_type) /= 'none') then
         if (.not. ltavg_write) then
            ltavg_write = .true.

            !*** modify tavg_outfile to conform to ccsm requirements
            !*** for ccsm, always want date string suffix
            if (lccsm) then
              call tavg_create_outfile_ccsm(tavg_outfile_orig,tavg_outfile)
              call tavg_create_suffix_ccsm(file_suffix,date_string='ymds')
              string = trim(file_suffix)
            else
              string = trim(runid)
            endif

            select case (trim(restart_type))
            case('even')
               file_suffix = trim(string)/&
                                          &/'.even'
            case('odd')
               file_suffix = trim(string)/&
                                          &/'.odd'
            case('end')
               file_suffix = trim(string)/&
                                          &/'.end'
            case('restart')
               ! do nothing
            case default
              if (lccsm) then
               file_suffix = trim(string)/&
                                          &/'.restart'
              else
               call tavg_create_suffix(file_suffix)
               file_suffix = trim(file_suffix)/&
                                               &/'.restart'
              endif
            end select
         
         endif ! .not. ltavg_write
      endif ! restart_type
   endif ! ltavg_on

!-----------------------------------------------------------------------
!
!  do the rest only if it is time to do a tavg dump (regular or restart)
!
!-----------------------------------------------------------------------

   if (ltavg_write) then

   if (lccsm .and. ltavg_fmt_out_nc .and. lreg_tavg_dump) then
     time_dim%active = .true.
   else
     time_dim%active = .false.
   endif


!-----------------------------------------------------------------------
!
!     compute global averages of tavg fields
!     do this before normalization
!
!-----------------------------------------------------------------------

      call tavg_global

!-----------------------------------------------------------------------
!
!     normalize time averages
!
!-----------------------------------------------------------------------

      call tavg_norm_field_all ('normalize')
 
!-----------------------------------------------------------------------
!
!  compute ccsm diagnostcs from tavg quantities
!  do this after normalizing TAVG_BUF quantities
!
!-----------------------------------------------------------------------

      if (lreg_tavg_dump) then

        !*** barotropic stream function
        call tavg_bsf_ccsm
     
        !*** MOC diagnostics
        call tavg_moc_ccsm
 
        !*** northward heat/salt transport diagnostics
        call tavg_transport_ccsm
      endif

      !*** compute local means
      call tavg_local_spatial_avg

!-----------------------------------------------------------------------
!
!     create data file descriptor
!
!-----------------------------------------------------------------------

      if (ltavg_fmt_out_nc) then
         tavg_file_desc = construct_file(tavg_fmt_out,                    &
                                      root_name  = trim(tavg_outfile),    &
                                      file_suffix= trim(file_suffix),     &
                                      title      = trim(runid),           &
                                      conventions='CF-1.0; ' /&
                   &/'http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-current.htm',&
                                      record_length = rec_type_real,      &
                                      recl_words=nx_global*ny_global)
      else
         call date_and_time(date=date_created, time=time_created)
         hist_string = char_blank
         write(hist_string,'(a23,a8,1x,a10)') & 
         'POP TAVG file created: ',date_created,time_created
         tavg_file_desc = construct_file(tavg_fmt_out,                    &
                                      root_name  = trim(tavg_outfile),    &
                                      file_suffix= trim(file_suffix),     &
                                      title      ='POP TAVG file',        &
                                      conventions='POP TAVG conventions', &
                                      history    = trim(hist_string),     &
                                      record_length = rec_type_real,      &
                                      recl_words=nx_global*ny_global)
      endif

!-----------------------------------------------------------------------
!
!     add scalar fields to file as file attributes
!
!-----------------------------------------------------------------------

      call add_attrib_file(tavg_file_desc, 'tavg_sum'    , tavg_sum)
      call add_attrib_file(tavg_file_desc, 'nsteps_total', nsteps_total)

      if (ltavg_fmt_out_nc .and. lreg_tavg_dump) then
        call tavg_add_attrib_file_ccsm (tavg_file_desc) 
      else
        call add_attrib_file(tavg_file_desc, 'tday'        , tday)
        call add_attrib_file(tavg_file_desc, 'iyear'       , iyear)
        call add_attrib_file(tavg_file_desc, 'imonth'      , imonth)
        call add_attrib_file(tavg_file_desc, 'iday'        , iday)
        call add_attrib_file(tavg_file_desc, 'beg_date'    , beg_date)
      endif

!-----------------------------------------------------------------------
!
!     open output file
!
!-----------------------------------------------------------------------

      call data_set (tavg_file_desc, 'open')

!-----------------------------------------------------------------------
!
!     writing fields to file requires two phases;
!     in this first phase, we define all the fields to be written
!
!-----------------------------------------------------------------------

      if (ltavg_fmt_out_nc .and. lreg_tavg_dump) then
        call tavg_define_time_ccsm       (tavg_file_desc)
        call tavg_define_coord_vars_ccsm (tavg_file_desc)
        call tavg_define_time_invar_ccsm (tavg_file_desc)
        call tavg_define_scalars_ccsm    (tavg_file_desc)
        call tavg_define_labels_ccsm     (tavg_file_desc)
      endif

      allocate(tavg_fields(num_avail_tavg_fields))

      !-----------------------------------------------------------------------
      ! Apply topography masking
      !-----------------------------------------------------------------------

      if (ltavg_fmt_out_nc .and. lreg_tavg_dump) then
        do nfield = 1,num_avail_tavg_fields  ! check all available fields
           loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer
           if (loc > 0) then  ! field is actually requested and in buffer
              if (avail_tavg_fields(nfield)%ndims == 2) then
              !*** mask 2D data fields
                   call tavg_mask(TAVG_BUF_2D(:,:,:,loc),avail_tavg_fields(nfield),1)
              else if (avail_tavg_fields(nfield)%ndims == 3) then
                   !*** mask 3D data fields
                   do k=1,km
                     TAVG_TEMP(:,:,:)=TAVG_BUF_3D(:,:,k,:,loc)
                     call tavg_mask(TAVG_TEMP(:,:,:),avail_tavg_fields(nfield),k)
                     TAVG_BUF_3D(:,:,k,:,loc)=TAVG_TEMP(:,:,:)
                   enddo ! k
              endif ! ndims
           endif ! loc
        enddo ! nfield
      endif ! lccsm .and. ltavg_fmt_out_nc .and. lreg_tavg_dump


      do nfield = 1,num_avail_tavg_fields  ! check all available fields

         loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer

         if (loc > 0) then  ! field is actually requested and in buffer

            !*** construct io_field descriptors for each field

            if (avail_tavg_fields(nfield)%ndims == 2) then

               tavg_fields(nfield) = construct_io_field(               &
                              avail_tavg_fields(nfield)%short_name,    &
                              dim1=i_dim, dim2=j_dim,time_dim=time_dim,&
                    long_name=avail_tavg_fields(nfield)%long_name,     &
                    units    =avail_tavg_fields(nfield)%units    ,     &
                    grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                   field_loc =avail_tavg_fields(nfield)%field_loc,     &
                  field_type =avail_tavg_fields(nfield)%field_type,    &
                 coordinates =avail_tavg_fields(nfield)%coordinates,   &
                  valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                   d2d_array =TAVG_BUF_2D(:,:,:,loc) ) !nonstandard, debugging only
#else
                   r2d_array =TAVG_BUF_2D(:,:,:,loc) )
#endif

            else if (avail_tavg_fields(nfield)%ndims == 3) then

               if (ltavg_fmt_out_nc .and. lreg_tavg_dump) then
                 select case (trim(avail_tavg_fields(nfield)%grid_loc(4:4)))
                   case('1')
                     z_dim = zt_dim
                   case('2')
                     z_dim = zw_dim
                 end select
               else
                 z_dim = k_dim
               endif ! lccsm

               if (avail_tavg_fields(nfield)%grid_loc(4:4) == '3') &
                 z_dim = zt_150m_dim

               tavg_fields(nfield) = construct_io_field(               &
                              avail_tavg_fields(nfield)%short_name,    &
                              dim1=i_dim, dim2=j_dim, dim3=z_dim,      &
                              time_dim=time_dim,                       &
                    long_name=avail_tavg_fields(nfield)%long_name,     &
                    units    =avail_tavg_fields(nfield)%units    ,     &
                    grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                   field_loc =avail_tavg_fields(nfield)%field_loc,     &
                  field_type =avail_tavg_fields(nfield)%field_type,    &
                 coordinates =avail_tavg_fields(nfield)%coordinates,   &
                  valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                   d3d_array =TAVG_BUF_3D(:,:,:,:,loc) ) !nonstandard, debugging only
#else
                   r3d_array =TAVG_BUF_3D(:,:,:,:,loc) ) 
#endif
            endif ! 2D/3D test

            if (lccsm .and. lreg_tavg_dump) &
               call tavg_add_attrib_io_field_ccsm (tavg_fields(nfield),nfield)

            call data_set (tavg_file_desc, 'define', tavg_fields(nfield))
         endif
      end do


!-----------------------------------------------------------------------
!
!     define nonstandard fields
!
!-----------------------------------------------------------------------


      if (ltavg_fmt_out_nc) then

        do nn = 1, num_avail_tavg_nstd_fields
                     
          call data_set_nstd_ccsm (                                 &
                           tavg_file_desc, 'define',                &
                           nstd_field_id,                           &
                           ndims_nstd_ccsm(nn),                     &
                           io_dims_nstd_ccsm(:,nn),                 &
                short_name=avail_tavg_nstd_fields(nn)%short_name,   &
                 long_name=avail_tavg_nstd_fields(nn)%long_name,    &
                     units=avail_tavg_nstd_fields(nn)%units,        &
               coordinates=avail_tavg_nstd_fields(nn)%coordinates,  &
             missing_value=avail_tavg_nstd_fields(nn)%missing_value,&
                fill_value=avail_tavg_nstd_fields(nn)%fill_value,   &
                    nftype=avail_tavg_nstd_fields(nn)%nftype        )


          if (nn == tavg_MOC) then
              moc_id = nstd_field_id
          elseif (nn == tavg_N_HEAT) then
              n_heat_id = nstd_field_id
          elseif (nn == tavg_N_SALT) then
              n_salt_id = nstd_field_id
          endif

        enddo !nn 

      endif ! lccsm

!-----------------------------------------------------------------------
!
!     write fields to file
!     in this second phase, we actually write the data for all the fields.
!     after writing all fields, the field descriptors are destroyed and the
!     file can be closed
!
!-----------------------------------------------------------------------

      if (ltavg_fmt_out_nc .and. lreg_tavg_dump) then
          call timer_start(timer_write_nstd)
          call tavg_write_vars_ccsm (tavg_file_desc,1,time_coordinate)
          call tavg_write_vars_ccsm (tavg_file_desc,num_ccsm_coordinates,ccsm_coordinates)
          call tavg_write_vars_ccsm (tavg_file_desc,num_ccsm_time_invar, ccsm_time_invar)
          call tavg_write_vars_ccsm (tavg_file_desc,num_ccsm_scalars,    ccsm_scalars)

          !*** ccsm nonstandard fields (uses data_set_nstd_ccsm to write)
          !*** time_bound, labels, transport diags
          call tavg_write_vars_nstd_ccsm (tavg_file_desc) 
          call timer_stop(timer_write_nstd)
      endif
 

      !*** standard 2D and 3D fields
      call timer_start(timer_write_std)
      do nfield = 1,num_avail_tavg_fields  ! check all available fields
         loc = avail_tavg_fields(nfield)%buf_loc ! locate field in buffer
         if (loc > 0) then  ! field is actually requested and in buffer
            call data_set (tavg_file_desc, 'write', tavg_fields(nfield))
            call destroy_io_field(tavg_fields(nfield))
         endif
      end do
      call timer_stop(timer_write_std)


!-----------------------------------------------------------------------
!
!     after writing all fields, the file can be closed
!
!-----------------------------------------------------------------------

      deallocate(tavg_fields)
      call data_set (tavg_file_desc, 'close')

      if (my_task == master_task) then
         write(stdout,blank_fmt)
         write(stdout,*) 'tavg file written: ', trim(tavg_file_desc%full_name)
         call POP_IOUnitsFlush(POP_stdout)
      endif

!-----------------------------------------------------------------------
!
!     if pointer files are used, write tavg filenames to pointer file
!     do this only for tavg restarts - not tavg dumps
!
!-----------------------------------------------------------------------

      if (luse_pointer_files .and. .not. lreg_tavg_dump) then
         call get_unit(nu)
         if (my_task == master_task) then
            tavg_pointer_file = trim(pointer_filename)/&
                                                       &/'.tavg'

            open(nu,file=tavg_pointer_file,form='formatted', &
                    status='unknown')
            write(nu,'(a)') trim(tavg_file_desc%full_name)
            close(nu)
         endif
         call release_unit(nu)
      endif

!-----------------------------------------------------------------------
!
!     reset time averages if this is a regular tavg dump (as opposed
!     to a restart dump)  if this is a restart dump, denormalize
!     in case a normal restart dump is written on the same timestep
!
!-----------------------------------------------------------------------

      if (lreg_tavg_dump) then
         tavg_sum = c0
         call time_stamp('now', 'ymd',date_string=beg_date)
         if (tavg_freq_iopt == freq_opt_nstep) &
            write(beg_date,'(i10)') nsteps_total
 
         call tavg_reset_field_all
      else
         call tavg_norm_field_all ('denormalize')
      endif


!-----------------------------------------------------------------------
!
!     get rid of file descriptor
!
!-----------------------------------------------------------------------

      call destroy_file(tavg_file_desc)
   endif ! ltavg_write

!-----------------------------------------------------------------------
!EOC

 end subroutine write_tavg

!***********************************************************************
!BOP
! !IROUTINE: read_tavg
! !INTERFACE:

 subroutine read_tavg

! !DESCRIPTION:
!  This routine reads a time average restart dump to continue
!  running time averages of requested fields.
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

   integer (i4) ::     &
     nu,               &   ! i/o unit
     iblock,           &   ! dummy block index
     n,                &   ! dummy for indexing character string
     in_fields,        &   ! num of fields in restart file
     nfield,           &   ! dummy field counter
     hdr_error,        &   ! error file for reading restart hdr
     in_nsteps_total,  &   ! nsteps_total according to tavg file
     in_iyear,         &   ! iyear according to tavg file
     in_imonth,        &   ! imonth according to tavg file
     in_iday,          &   ! iday according to tavg file
     loc                   ! buffer location

   real (r8) ::        &
     in_tday               ! tday according to tavg file

   character (char_len) ::  &
     header_filename,   &  ! filename for restart contents
     char_temp,         &  ! for string manipulation
     tavg_pointer_file     ! filename for pointer file containing
                           !   location/name of last restart file

   type (io_field_desc), dimension(:), allocatable :: &
      tavg_fields          ! io field description for each field in file


!-----------------------------------------------------------------------
!
!  if pointer files are used, pointer file and must be read to get 
!  actual filenames
!
!-----------------------------------------------------------------------

   call get_unit(nu)

   if (luse_pointer_files) then

      if (my_task == master_task) then
         tavg_pointer_file = char_blank
         tavg_pointer_file = trim(pointer_filename)/&
                                                   &/'.tavg'
         write(stdout,*) 'Reading pointer file: ', &
                         trim(tavg_pointer_file)
         open(nu, file=trim(tavg_pointer_file), form='formatted', &
                  status='old')
         read(nu,'(a)') tavg_infile
         close(nu)
      endif
      call broadcast_scalar(tavg_infile, master_task)

   endif

   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  define input file
!
!-----------------------------------------------------------------------

   tavg_file_desc = construct_file (tavg_fmt_in,                   &
                                    full_name=trim(tavg_infile),   &
                                    record_length = rec_type_real, &
                                    recl_words=nx_global*ny_global)

!-----------------------------------------------------------------------
!
!  define scalar fields in file as file attributes to be read during
!  open
!
!-----------------------------------------------------------------------

   if (lccsm) then
   else
      call add_attrib_file(tavg_file_desc, 'tday'        , tday)
      call add_attrib_file(tavg_file_desc, 'iyear'       , iyear)
      call add_attrib_file(tavg_file_desc, 'imonth'      , imonth)
      call add_attrib_file(tavg_file_desc, 'iday'        , iday)
      call add_attrib_file(tavg_file_desc, 'beg_date'    , beg_date)
   endif
   call add_attrib_file(tavg_file_desc, 'nsteps_total', nsteps_total)
   call add_attrib_file(tavg_file_desc, 'tavg_sum'    , tavg_sum)

!-----------------------------------------------------------------------
!
!  open input file
!  this will also extract scalar variables which are defined as
!  file attributes
!
!-----------------------------------------------------------------------

   call data_set (tavg_file_desc, 'open_read')

   if (lccsm) then
   else
      call extract_attrib_file(tavg_file_desc, 'beg_date', beg_date)
      !call extract_attrib_file(tavg_file_desc, 'tday', in_tday)
      !call extract_attrib_file(tavg_file_desc, 'iyear', in_iyear)
      !call extract_attrib_file(tavg_file_desc, 'imonth', in_imonth)
      !call extract_attrib_file(tavg_file_desc, 'iday', in_iday)
   endif
   call extract_attrib_file(tavg_file_desc, 'nsteps_total', &
                                          in_nsteps_total)
   call extract_attrib_file(tavg_file_desc, 'tavg_sum', tavg_sum)

   !*** report nsteps total and tavg_sum
   if (my_task == master_task) then
      write(stdout,'(i6,a29,i6,a35)') &
      in_nsteps_total,' nsteps_total in tavg restart', &
      nsteps_total,   ' nsteps_total in current simulation'
      write(stdout,*) ' tavg_sum = ', tavg_sum, ' in tavg restart'
      call POP_IOUnitsFlush(POP_stdout)
   endif

   !*** check nsteps total for validity
   if (in_nsteps_total /= nsteps_total) then
      if (my_task == master_task) then
      write(stdout,'(i6,a29,i6,a35)') &
         in_nsteps_total,' nsteps_total in tavg restart', &
         nsteps_total,   ' nsteps_total in current simulation'
      endif
      call exit_POP(sigAbort,'TAVG:restart file has wrong time step?')
   endif

!-----------------------------------------------------------------------
!
!  define requested fields to read in from file
!  NOTE: This requires that the tavg_contents file is consistent
!  with the tavg restart file.  There are currently no checks on this.
!
!-----------------------------------------------------------------------

   !*** define dimensions

   allocate(tavg_fields(num_avail_tavg_fields))

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
         if (avail_tavg_fields(nfield)%ndims == 2) then

            tavg_fields(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            dim1=i_dim, dim2=j_dim,                  &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                 d2d_array =TAVG_BUF_2D(:,:,:,loc) ) !nonstandard, debugging only
#else
                 r2d_array =TAVG_BUF_2D(:,:,:,loc) )
#endif
         else if (avail_tavg_fields(nfield)%ndims == 3) then

            if (avail_tavg_fields(nfield)%grid_loc(4:4) == '3') then
               z_dim = zt_150m_dim
            else
               z_dim = k_dim
            endif

            tavg_fields(nfield) = construct_io_field(                &
                            avail_tavg_fields(nfield)%short_name,    &
                            dim1=i_dim, dim2=j_dim, dim3=z_dim,      &
                  long_name=avail_tavg_fields(nfield)%long_name,     &
                  units    =avail_tavg_fields(nfield)%units    ,     &
                  grid_loc =avail_tavg_fields(nfield)%grid_loc ,     &
                 field_loc =avail_tavg_fields(nfield)%field_loc,     &
                field_type =avail_tavg_fields(nfield)%field_type,    &
                valid_range=avail_tavg_fields(nfield)%valid_range,   &
#ifdef TAVG_R8
                 d3d_array =TAVG_BUF_3D(:,:,:,:,loc) ) !nonstandard, debugging only
#else
                 r3d_array =TAVG_BUF_3D(:,:,:,:,loc) )
#endif
         endif

         call data_set (tavg_file_desc, 'define', tavg_fields(nfield))
      endif
   end do

!-----------------------------------------------------------------------
!
!  now we actually read each field
!  after reading, get rid of io field descriptors and close file
!
!-----------------------------------------------------------------------

   do nfield = 1,num_avail_tavg_fields
      loc = avail_tavg_fields(nfield)%buf_loc
      if (loc > 0) then
         call data_set (tavg_file_desc, 'read', tavg_fields(nfield))
         call destroy_io_field(tavg_fields(nfield))
      endif
   end do

   deallocate(tavg_fields)
   call data_set (tavg_file_desc, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' file read: ', tavg_infile
   endif

   call destroy_file(tavg_file_desc)
   call release_unit(nu)

!-----------------------------------------------------------------------
!
!  de-normalize sums
!
!-----------------------------------------------------------------------

   call tavg_norm_field_all ('denormalize')

!-----------------------------------------------------------------------

   call tavg_global   ! print global sums of time averages

!-----------------------------------------------------------------------
!EOC

 end subroutine read_tavg

!***********************************************************************
!BOP
! !IROUTINE: tavg_global
! !INTERFACE:

 subroutine tavg_global

! !DESCRIPTION:
!  Calculates and print global integrals of time average fields
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

   integer (i4) ::     &
      k,               &   ! vertical level index
      ifield,          &   ! field identifier
      iblock,          &   ! block index
      nfield,          &   ! dummy field index
      field_loc,       &   ! field location (center,Nface,Eface,NEcorner)
      field_type           ! field type (scalar, vector, angle)

   real (r8) ::        &
      tavg_field_sum,  &   ! sum of tavg field
      tavg_norm            ! normalization for average

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK               ! temp for holding area_weighted field

   real (r8), dimension (:,:), allocatable ::  &
      RMASK              ! topography mask for global sum

   character (char_len) ::  &
      time_string,          &
      date_string

!-----------------------------------------------------------------------
!
!  calculate globally-integrated time average of each chosen 2d field
!
!-----------------------------------------------------------------------

   allocate (RMASK(nx_block,ny_block), &
             WORK (nx_block,ny_block,nblocks_clinic))

   call time_stamp ('now','mdy',date_string=date_string,time_string=time_string)

   if (my_task == master_task) then
     write (stdout,blank_fmt)
     write (stdout,*) 'Global Time Averages: ' // trim(date_string) // ' ' // trim(time_string)
   endif

   fields_loop: do nfield=1,num_avail_tavg_fields
      ifield = avail_tavg_fields(nfield)%buf_loc
      if (ifield > 0) then

         field_loc  = avail_tavg_fields(nfield)%field_loc
         field_type = avail_tavg_fields(nfield)%field_type

         if (avail_tavg_fields(nfield)%method == tavg_method_avg) then
            tavg_norm = tavg_sum
         else if (avail_tavg_fields(nfield)%method == tavg_method_qflux) then
            tavg_norm = tavg_sum_qflux
         else
            tavg_norm = c1
         endif

         if (tavg_norm == c0) then
            if (my_task == master_task) then
              write(stdout,*) 'Cannot compute global integral of ', &
                               trim (avail_tavg_fields(nfield)%short_name)
              call POP_IOUnitsFlush(POP_stdout)
            endif
            cycle fields_loop 
            !*** call exit_POP (SigAbort,'ERROR: tavg_norm = 0 in tavg_global')
         endif

         !*** 2-d fields

         if (avail_tavg_fields(nfield)%ndims == 2) then

            !$OMP PARALLEL DO
            do iblock = 1,nblocks_clinic
               select case(field_loc)
               case(field_loc_center)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    TAREA(:,:,iblock)*RCALCT(:,:,iblock)
               case(field_loc_NEcorner)
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               case default ! make U cell the default for all other cases
                  WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,ifield)* &
                                    UAREA(:,:,iblock)*RCALCU(:,:,iblock)
               end select
            end do
            !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*area_u)
            end select

         !*** 3-d fields

         else
      
            !$OMP PARALLEL DO PRIVATE(k)
            do iblock = 1,nblocks_clinic
               WORK(:,:,iblock) = c0

               select case(field_loc)

               case(field_loc_center)
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMT(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        TAREA(:,:,iblock)*RMASK
                  end do

               case(field_loc_NEcorner)
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMU(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               case default ! make U cell the default for all other cases
                  do k=1,km
                     RMASK = merge(c1, c0, k <= KMU(:,:,iblock)) 
                     WORK(:,:,iblock) = WORK(:,:,iblock) + dz(k)* &
                                        TAVG_BUF_3D(:,:,k,iblock,ifield)* &
                                        UAREA(:,:,iblock)*RMASK
                  end do

               end select
            end do
            !$OMP END PARALLEL DO

            tavg_field_sum = global_sum(WORK, distrb_clinic, field_loc)

            select case(field_loc)
            case(field_loc_center)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_t)
            case(field_loc_NEcorner)
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            case default ! make U cell the default for all other cases
               tavg_field_sum = tavg_field_sum/(tavg_norm*volume_u)
            end select

         endif

         if (my_task == master_task) then
            write (stdout,*) trim(avail_tavg_fields(nfield)%short_name), &
                             ': ', tavg_field_sum
         endif
      endif
   end do fields_loop

   deallocate (RMASK, WORK)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_global

!***********************************************************************
!BOP
! !IROUTINE: tavg_global_sum_2D
! !INTERFACE:

 function tavg_global_sum_2D (id)

! !DESCRIPTION:
!  Calculates the global sum of a requested 2D time-averaged field
!  Presently, this is a special-purpose routine used only by the
!  budget_diagnostics module.  It could be extended and generalized,
!  and used with a revised version of subroutine tavg_global, if time
!  and interest allow.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) ::  &
      id                               ! identifier of time-averaged field 

! !OUTPUT PARAMETERS:

   real (r8) :: &
      tavg_global_sum_2D     ! result of this function: the global sum of
                             !   a requested 2D time-averaged field


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) ::     &
      loc,             &   ! field buffer location
      iblock,          &   ! block index
      nfield,          &   ! dummy field index
      field_loc,       &   ! field location (center,Nface,Eface,NEcorner)
      field_type           ! field type (scalar, vector, angle)

   real (r8) ::            &
      tavg_norm                ! normalization for average

   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) ::  &
      WORK               ! temp for holding area_weighted field

   real (r8), dimension (nx_block, ny_block) ::  &
      RMASK              ! topography mask for global sum


!-----------------------------------------------------------------------
!
!  does this field exist?
!
!-----------------------------------------------------------------------

   if (.not. tavg_requested(id)) then
      call document ('tavg_global_sum_2D', 'id  ', id)
      call document ('tavg_global_sum_2D', 'tavg_requested(id) ', tavg_requested(id))
      call exit_POP(sigAbort,'(tavg_global_sum_2D) ERROR: invalid field request')
   endif

!-----------------------------------------------------------------------
!
!  is this a 2D field?
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%ndims /= 2) then
      call exit_POP(sigAbort,'(tavg_global_sum_2D) ERROR: invalid dimension')
   endif

!-----------------------------------------------------------------------
!
!  identify the requested field
!
!-----------------------------------------------------------------------

   field_loc  = avail_tavg_fields(id)%field_loc
   field_type = avail_tavg_fields(id)%field_type

!-----------------------------------------------------------------------
!
!  identify the buffer location in TAVG_BUF_2D
!
!-----------------------------------------------------------------------

   loc = avail_tavg_fields(id)%buf_loc
   if (loc <= 0) &
     call exit_POP(sigAbort, &
                    'tavg_global_sum_2D: invalid loc')

!-----------------------------------------------------------------------
!
!  compute the global average of the 2D field
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%method == tavg_method_avg) then
      tavg_norm = tavg_sum
   else if (avail_tavg_fields(id)%method == tavg_method_qflux) then
      tavg_norm = tavg_sum_qflux
   else
      tavg_norm = c1
   endif

   if (tavg_norm == c0) &
     call exit_POP (SigAbort,'ERROR: tavg_norm = 0 in tavg_global_sum_2D')

   !$OMP PARALLEL DO
   do iblock = 1,nblocks_clinic
      select case(field_loc)
        case(field_loc_center)
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               TAREA(:,:,iblock)*RCALCT(:,:,iblock)
        case(field_loc_NEcorner)
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               UAREA(:,:,iblock)*RCALCU(:,:,iblock)
        case default ! make U cell the default for all other cases
           WORK(:,:,iblock)  = TAVG_BUF_2D(:,:,iblock,loc)* &
                               UAREA(:,:,iblock)*RCALCU(:,:,iblock)
      end select
   end do
   !$OMP END PARALLEL DO

   tavg_global_sum_2D = global_sum(WORK, distrb_clinic, field_loc)

   select case(field_loc)
     case(field_loc_center)
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_t)
     case(field_loc_NEcorner)
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_u)
     case default ! make U cell the default for all other cases
        tavg_global_sum_2D = tavg_global_sum_2D/(tavg_norm*area_u)
   end select



!-----------------------------------------------------------------------
!EOC

 end function tavg_global_sum_2D

!***********************************************************************
!BOP
! !IROUTINE: tavg_increment_sum_qflux
! !INTERFACE:

 subroutine tavg_increment_sum_qflux(const)

! !DESCRIPTION:
!  Increment the scalar tavg_sum_qflux with the weight for this timestep.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), intent(in) ::  &
      const
!EOP
!BOC
!-----------------------------------------------------------------------

   tavg_sum_qflux = tavg_sum_qflux + const ! const = tlast_ice

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_increment_sum_qflux

!***********************************************************************
!BOP
! !IROUTINE: accumulate_tavg_field
! !INTERFACE:

 subroutine accumulate_tavg_field(ARRAY,field_id,block,k,const)

! !DESCRIPTION:
!  This routine updates a tavg field.  If the time average of the
!  field is requested, it accumulates a time sum of a field by 
!  multiplying by the time step and accumulating the sum into the 
!  tavg buffer array.  If the min or max of a field is requested, it
!  checks the current value and replaces the min, max if the current
!  value is less than or greater than the stored value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block,           &! local block address (in baroclinic distribution)
      k,               &! vertical level
      field_id          ! index into available fields for tavg field info

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      ARRAY             ! array of data for this block to add to 
                        !  accumulated sum in tavg buffer
   real (r8), optional, intent(in) ::  &
      const
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc,            &! location of field in tavg buffer
      ndims               ! rank of field (2=2d,3=3d)

!-----------------------------------------------------------------------
!
!  get buffer location and field info from avail_tavg_field array
!
!-----------------------------------------------------------------------

   bufloc = avail_tavg_fields(field_id)%buf_loc

   if (bufloc <= 0) &
     call exit_POP(sigAbort, &
                    'tavg: attempt to accumulate bad tavg field')

   ndims = avail_tavg_fields(field_id)%ndims

   if ((ndims == 3) .and. &
       (avail_tavg_fields(field_id)%grid_loc(4:4) == '3') .and. &
       (k > zt_150m_levs)) return

!-----------------------------------------------------------------------
!
!  update the field into the tavg buffer
!
!-----------------------------------------------------------------------

   select case (avail_tavg_fields(field_id)%method)

   case (tavg_method_avg)  ! accumulate running time sum for time avg
      if (ndims == 2) then
         TAVG_BUF_2D(:,:,block,bufloc) = &
         TAVG_BUF_2D(:,:,block,bufloc) + dtavg*ARRAY
      else
         TAVG_BUF_3D(:,:,k,block,bufloc) = &
         TAVG_BUF_3D(:,:,k,block,bufloc) + dtavg*ARRAY
      endif
   case (tavg_method_qflux)  
         TAVG_BUF_2D(:,:,block,bufloc) =  &
         TAVG_BUF_2D(:,:,block,bufloc) + const*max (c0,ARRAY)
   case (tavg_method_min)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY < TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY < TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif

   case (tavg_method_max)  ! replace with current minimum value
      if (ndims == 2) then
         where (ARRAY > TAVG_BUF_2D(:,:,block,bufloc))
            TAVG_BUF_2D(:,:,block,bufloc) = ARRAY
         end where
      else
         where (ARRAY > TAVG_BUF_3D(:,:,k,block,bufloc))
            TAVG_BUF_3D(:,:,k,block,bufloc) = ARRAY
         end where
      endif

   case default
   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine accumulate_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: tavg_reset_field_all
! !INTERFACE:

 subroutine tavg_reset_field_all
 
! !DESCRIPTION:
!  Resets all TAVG_BUF_2D and TAVG_BUF_3D arrays 
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
 
   integer (int_kind) ::  &
      iblock,             &
      nfield


  !$OMP PARALLEL DO PRIVATE(nfield)
  do iblock=1,nblocks_clinic
     do nfield=1,tavg_bufsize_2d
       select case (TAVG_BUF_2D_METHOD(nfield))
        case (tavg_method_avg)
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
        case (tavg_method_qflux)
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
           tavg_sum_qflux = c0
        case (tavg_method_min)
           TAVG_BUF_2D(:,:,  iblock,nfield) = bignum
        case (tavg_method_max)
           TAVG_BUF_2D(:,:,  iblock,nfield) = -bignum
        case default
           TAVG_BUF_2D(:,:,  iblock,nfield) = c0
        end select
     end do
 
     do nfield=1,tavg_bufsize_3d
        select case (TAVG_BUF_3D_METHOD(nfield))
        case (tavg_method_avg)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
        case (tavg_method_min)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = bignum
        case (tavg_method_max)
           TAVG_BUF_3D(:,:,:,iblock,nfield) = -bignum
        case default
           TAVG_BUF_3D(:,:,:,iblock,nfield) = c0
        end select
     end do
  end do
  !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_reset_field_all



!***********************************************************************
!BOP
! !IROUTINE: tavg_norm_field_all
! !INTERFACE:

 subroutine tavg_norm_field_all (norm_flag)
 
! !DESCRIPTION:
!  Normalizes or de-normalizes all TAVG_BUF_2D and TAVG_BUF_3D arrays 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 
   character (*), intent(in) ::  &
      norm_flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
   real (r8) ::  &
      factor,    &
      factorq

   integer (int_kind) ::  &
      iblock,             &
      nfield

   
   select case (trim(norm_flag))
 
      case ('normalize')
        if (tavg_sum /= 0) then
          factor  = c1/tavg_sum
        else
          call exit_POP(sigAbort, &
           'ERROR in tavg_norm_field_all; attempt to divide by zero tavg_sum')
        endif
        if (tavg_sum_qflux /= 0) then
          factorq = c1/tavg_sum_qflux
        elseif (tavg_freq_iopt == freq_opt_nhour   .or.  &
                tavg_freq_iopt == freq_opt_nsecond .or.  &
                tavg_freq_iopt == freq_opt_nstep         ) then
          ! do nothing; these frequencies are not compatable with time-averaged qflux
          factorq = c1
        else
          call exit_POP(sigAbort, &
           'ERROR in tavg_norm_field_all; attempt to divide by zero tavg_sum_qflux')
        endif
      case ('denormalize')
          factor  = tavg_sum
          factorq = tavg_sum_qflux
      case default
          call exit_POP(sigAbort,'ERROR in tavg_norm_field_all; unknown option')
   end select
 
 
   !$OMP PARALLEL DO PRIVATE(nfield)
   do iblock=1,nblocks_clinic
      do nfield=1,tavg_bufsize_2d
         select case (TAVG_BUF_2D_METHOD(nfield))
           case (tavg_method_avg)
             TAVG_BUF_2D(:,:,  iblock,nfield) = &
             TAVG_BUF_2D(:,:,  iblock,nfield)*factor
           case (tavg_method_qflux) 
             TAVG_BUF_2D(:,:,  iblock,nfield) = &
             TAVG_BUF_2D(:,:,  iblock,nfield)*factorq
         end select
      end do
      do nfield=1,tavg_bufsize_3d
         if (TAVG_BUF_3D_METHOD(nfield) == tavg_method_avg) then
            TAVG_BUF_3D(:,:,:,iblock,nfield) = &
            TAVG_BUF_3D(:,:,:,iblock,nfield)*factor
         endif
      end do


   end do
   !$OMP END PARALLEL DO
 

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_norm_field_all


!***********************************************************************
!BOP
! !IROUTINE: define_tavg_field
! !INTERFACE:

 subroutine define_tavg_field(id, short_name, ndims, tavg_method,       &
                                  long_name, units,                     &
                                  grid_loc, valid_range,                &
                                  field_loc, field_type,coordinates,    &
                                  scale_factor,                         &
                                  nftype,                               &
                                  nstd_fields,                          &
                                  num_nstd_fields, max_nstd_fields      )

! !DESCRIPTION:
!  Initializes description of an available field and returns location
!  in the available fields array for use in later tavg calls.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      id                ! location in avail_fields array for use in
                        ! later tavg routines 

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      short_name                ! short name for field

   integer (i4), intent(in) :: &
      ndims                     ! number of dims of the field

   integer (i4), intent(in), optional :: &
      field_loc,              &! location in grid 
      field_type,             &! type of field (scalar, vector, angle)
      tavg_method              ! id for method of averaging
                               ! default is tavg_method_avg

   character(*), intent(in), optional :: &
      long_name,              &! long descriptive name for field
      units,                  &! physical units for field
      coordinates              ! CF coordinates

   character(4), intent(in), optional :: &
      grid_loc                 ! location in grid (in 4-digit code)

   real (rtavg), intent(in), optional :: &
      scale_factor             ! scale factor

   real (r4), dimension(2), intent(in), optional :: &
      valid_range              ! min/max

   character(*), intent(in), optional :: &
      nftype                    ! type string 

   integer (int_kind), intent(in), optional ::  &
      max_nstd_fields

! !INPUT/OUTPUT PARAMETERS:

   type (tavg_field_desc_ccsm), dimension(:), intent(inout), optional ::  &
      nstd_fields

   integer (int_kind), intent(inout), optional ::  &
      num_nstd_fields

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (tavg_field_desc_ccsm) ::  &
      tavg_field

   integer (int_kind) ::  &
      num_fields

   logical (log_kind) ::  &
      error,              &
      nonstandard_fields

!-----------------------------------------------------------------------
!
!  increment the number of defined fields and make sure it does not
!  exceed the maximum
!  return the id as the current number
!
!-----------------------------------------------------------------------

   error = .false.
   nonstandard_fields = .false.

   if (present (nstd_fields) ) then
      nonstandard_fields = .true.
      num_nstd_fields = num_nstd_fields + 1
      num_fields = num_nstd_fields
      if (num_nstd_fields > max_nstd_fields) error = .true.
   else
      num_avail_tavg_fields = num_avail_tavg_fields + 1
      num_fields = num_avail_tavg_fields
      if (num_avail_tavg_fields > max_avail_tavg_fields) error = .true.
   endif

   if (error) call exit_POP(sigAbort, &
                    'tavg: defined tavg fields > max allowed')
 
   id = num_fields

   if (my_task == master_task .and. nsteps_run <= 1) then
     write(stdout,*) 'define_tavg_field: id = ', id, ' ', short_name
     call POP_IOUnitsFlush(POP_stdout)
   endif

!-----------------------------------------------------------------------
!
!  now fill the local field descriptor
!
!-----------------------------------------------------------------------

   tavg_field%ndims      = ndims
   tavg_field%short_name = short_name
   tavg_field%buf_loc    = 0  ! will be reset later

   if (present(long_name)) then
      tavg_field%long_name = long_name
   else
      tavg_field%long_name = char_blank
   endif

   if (present(units)) then
      tavg_field%units = units
   else
      tavg_field%units = char_blank
   endif

   if (present(grid_loc)) then
      tavg_field%grid_loc = grid_loc
   else
      tavg_field%grid_loc = '    '
   endif

   if (present(tavg_method)) then
      tavg_field%method = tavg_method
   else
      tavg_field%method = tavg_method_avg
   endif


   if (present(scale_factor)) then
      tavg_field%scale_factor = scale_factor
      if (scale_factor /= 0.0_rtavg) then
        tavg_field%fill_value    = undefined_nf/scale_factor
        tavg_field%missing_value = undefined_nf/scale_factor
      endif
   else
      tavg_field%scale_factor  = undefined_nf
      tavg_field%fill_value    = undefined_nf
      tavg_field%missing_value = undefined_nf
   endif


   if (present(valid_range)) then
      tavg_field%valid_range = valid_range
   else
      tavg_field%valid_range = undefined
   endif

   if (present(coordinates)) then
      tavg_field%coordinates = coordinates
   else
      tavg_field%coordinates = char_blank
   endif

   !*** set field location, field type used by i/o, ghost cell update
   !*** and other communication routines.  because ghost cells for tavg
   !*** fields are not typically used, the default is field_xxx_noupdate

   if (present(field_loc)) then
      tavg_field%field_loc = field_loc
   else
      tavg_field%field_loc = field_loc_noupdate
      if (present(grid_loc)) then
         !*** try to decode field location from grid_loc
         if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '1') then
            tavg_field%field_loc = field_loc_center
         else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '2') then
            tavg_field%field_loc = field_loc_NEcorner
         else if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '2') then
            tavg_field%field_loc = field_loc_Nface
         else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '1') then
            tavg_field%field_loc = field_loc_Eface
         endif
      endif
   endif

   if (present(field_type)) then
      tavg_field%field_type = field_type
   else
      tavg_field%field_type = field_type_noupdate
   endif

   if (present(nftype)) then
      tavg_field%nftype = trim(nftype)
   else
      tavg_field%nftype = 'float'
   endif


   if (nonstandard_fields) then
     nstd_fields(id) = tavg_field
   else
     avail_tavg_fields(id) = tavg_field
   endif

   if (my_task == master_task .and. tavg_debug > 0) then
     call document ('define_tavg_field',  trim(tavg_field%short_name))
     call document ('define_tavg_field', 'buffer id number ', id)
     call document ('define_tavg_field',  trim(tavg_field%long_name))
     call document ('define_tavg_field',  trim(tavg_field%units))
     call document ('define_tavg_field',  trim(tavg_field%grid_loc))
     call document ('define_tavg_field', 'fill_value',   tavg_field%fill_value)
     call document ('define_tavg_field', 'missing_value',tavg_field%missing_value)
     call document ('define_tavg_field', 'scale_factor', tavg_field%scale_factor)
     call document ('define_tavg_field',  trim(tavg_field%coordinates))
     call document ('define_tavg_field', 'field_type',   tavg_field%field_type)
     call document ('define_tavg_field',  trim(tavg_field%nftype))
     call document ('define_tavg_field', 'end subroutine ')
   endif
!-----------------------------------------------------------------------
!EOC

 end subroutine define_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: request_tavg_field
! !INTERFACE:

 subroutine request_tavg_field(short_name)

! !DESCRIPTION:
!  This field marks an available field as requested and computes
!  the location in the tavg buffer array.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      short_name                ! the short name of the field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind)  ::  &
      n,                   &! loop index
      id                    ! location of field in avail_fields array

!-----------------------------------------------------------------------
!
!  search for field with same name
!
!-----------------------------------------------------------------------

   id = 0

   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == short_name) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   if (id == 0) then
      write(stdout,*) 'Requested ', trim(short_name)
#ifdef CCSMCOUPLED
      call shr_sys_flush(stdout)
#endif
      call exit_POP(sigAbort,'subroutine request_tavg_field: requested field unknown')
   endif

!-----------------------------------------------------------------------
!
!  set the position in the buffer and advance the buffer position
!  for the next field
!
!-----------------------------------------------------------------------

   if (avail_tavg_fields(id)%ndims == 2) then
      tavg_bufsize_2d = tavg_bufsize_2d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_2d
   else if (avail_tavg_fields(id)%ndims == 3) then
      tavg_bufsize_3d = tavg_bufsize_3d + 1
      avail_tavg_fields(id)%buf_loc = tavg_bufsize_3d
   endif

!-----------------------------------------------------------------------
!
!  if tavg is on, but not started yet, set the buf_loc to -buf_loc
!  to show that it is requested, but will not return true for
!  requested_tavg_field
!
!-----------------------------------------------------------------------

   if (.not. ltavg_on) then
      avail_tavg_fields(id)%buf_loc = -avail_tavg_fields(id)%buf_loc
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine request_tavg_field

!***********************************************************************
!BOP
! !IROUTINE: tavg_requested
! !INTERFACE:

 function tavg_requested(id)

! !DESCRIPTION:
!  This function determines whether an available (defined) tavg field
!  has been requested by a user (through the input contents file) and 
!  returns true if it has.  Note that if tavg has been turned off, 
!  the function will always return false.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      tavg_requested     ! result of checking whether the field has
                         !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if zero, the field has not been
!  requested
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
      call exit_POP(sigAbort,'tavg_requested: invalid tavg id')
   endif

   if (avail_tavg_fields(id)%buf_loc > 0) then
      tavg_requested = .true.
   else
      tavg_requested = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function tavg_requested

!***********************************************************************
!BOP
! !IROUTINE: set_in_tavg_contents
! !INTERFACE:

 function set_in_tavg_contents(id)

! !DESCRIPTION:
!  This function determines whether a tavg field has been set in
!  the input contents file and returns true if it has.  This function is
!  different from tavg_requested in that ltavg_on status is irrelevent.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id                   ! id returned by the define function which
                           !   gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) ::      &
      set_in_tavg_contents    ! result of checking whether the field has
                               !   been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check the buffer location - if not zero, then the field is in the
!  tavg contents file
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_tavg_fields) then
     set_in_tavg_contents = .false.
   elseif (avail_tavg_fields(id)%buf_loc /= 0) then
     set_in_tavg_contents = .true.
   else
     set_in_tavg_contents = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function set_in_tavg_contents

!***********************************************************************
!BOP
! !IROUTINE: tavg_id
! !INTERFACE:

 function tavg_id(short_name,quiet)

! !DESCRIPTION:
!  This function determines whether a tavg field has been defined
!  by some module.  If so, then the id of that field is returned.
!  This function does not cause a model exit if the field has not
!  been defined; error-checking must be done by the calling routine.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)  :: &
      short_name                  ! the short name of the tavg field

   logical (log_kind), intent(in), optional :: &
      quiet                        ! do not print error message

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      tavg_id               ! id of the tavg field, if it exists

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  &
     id,                  &
     n

   logical (log_kind) ::  &
     msg

   id = 0

   srch_loop: do n=1,num_avail_tavg_fields
      if (trim(avail_tavg_fields(n)%short_name) == trim(short_name)) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   msg = .true.
   if (present(quiet)) then
     if (quiet) msg = .false.
   endif
     
   if (id == 0 .and. msg) then
      if (my_task == master_task)  & 
          write(stdout,*) 'Field ', trim(short_name), ' has not been defined.'
   endif

   tavg_id = id

!-----------------------------------------------------------------------
!EOC

 end function tavg_id

!***********************************************************************
!BOP
! !IROUTINE: tavg_create_suffix
! !INTERFACE:

 subroutine tavg_create_suffix(file_suffix)

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix           ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variable
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      cindx1, cindx2,    &! indices into character strings
      len_date            ! length of date string

   character (char_len) :: &
      char_temp            ! temp character space (for removing spaces)

   character (10) :: &
      cstep_beg,     &! beginning step  of this particular average
      cstep_end,     &! ending    step  of this particular average
      cdate           ! character string with yyyymmdd and optional 
                      ! separator (eg yyyy-mm-dd)

   character (4) :: &
      cyear_beg,    &! beginning year  of this particular average
      cyear_end      ! end       year  of this particular average

   character (2) :: &
      cmonth_beg,   &! beginning month of this particular average
      cmonth_end,   &! end       month of this particular average
      cday_beg,     &! beginning day   of this particular average
      cday_end       ! end       day   of this particular average

!-----------------------------------------------------------------------
!
!  start suffix with runid
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   cindx2 = len_trim(runid) + 1
   file_suffix(1:cindx2) = trim(runid)/&
                                       &/'.'
   cindx1 = cindx2 + 1
   
!-----------------------------------------------------------------------
!
!  extract beginning year, month, day or time step from beg_date
!  and determine end date
!
!-----------------------------------------------------------------------

   cdate = adjustl(beg_date)

   !***
   !*** use step numbers if tavg freq option is nstep
   !***

   cstep_beg  = trim(cdate) ! in case beg_date actually step number
   write(cstep_end,'(i10)') nsteps_total - 1
   cdate  = adjustl(cstep_end)
   cstep_end = trim(cdate)

   call time_stamp('last','ymd',date_string=cdate)  ! last date

   if (date_separator == ' ') then  ! no date separator
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(5:6)
      cday_beg   = beg_date(7:8)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(5:6)
      cday_end   = cdate(7:8)
   else
      cyear_beg  = beg_date(1:4)
      cmonth_beg = beg_date(6:7)
      cday_beg   = beg_date(9:10)
      cyear_end  = cdate(1:4)
      cmonth_end = cdate(6:7)
      cday_end   = cdate(9:10)
   endif

!-----------------------------------------------------------------------
!
!  create time portion of suffix based on frequency option
!  note that concatenation operator split across lines to avoid
!   problems with some cpp preprocessors
!
!-----------------------------------------------------------------------

   select case (tavg_freq_iopt)
   case (freq_opt_nyear)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 3
         file_suffix(cindx1:cindx2) = cyear_end
      else
         cindx2 = cindx1 + 8
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                                &/'-'/&
                                                      &/cyear_end
      endif

   case (freq_opt_nmonth)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 5
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end
      else
         cindx2 = cindx1 + 12
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end
      endif

   case (freq_opt_nday)
      if (tavg_freq == 1) then
         cindx2 = cindx1 + 7
         file_suffix(cindx1:cindx2) = cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      else
         cindx2 = cindx1 + 16
         file_suffix(cindx1:cindx2) = cyear_beg/&
                                    &/cmonth_beg/&
                                    &/cday_beg/&
                                    &/'-'/&
                                    &/cyear_end/&
                                    &/cmonth_end/&
                                    &/cday_end
      endif

   case (freq_opt_nstep)
      cindx2 = cindx1 + len_trim(cstep_beg) + len_trim(cstep_end)
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   case default  ! use nstep for other options
      cindx2 = len_trim(cstep_beg) + len_trim(cstep_end) + 1
      file_suffix(cindx1:cindx2) = trim(cstep_beg)/&
                                 &/'-'/&
                                 &/trim(cstep_end)

   end select
 
!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_suffix

!***********************************************************************
!BOP
! !IROUTINE: tavg_create_suffix_ccsm
! !INTERFACE:

 subroutine tavg_create_suffix_ccsm(file_suffix,date_string) 

! !DESCRIPTION:
!  Creates a suffix to append to output filename based on frequency 
!  option and averaging interval.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in), optional :: &
      date_string

! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix           ! suffix to append to root filename

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      char_temp,           &! temp character space
      ccsm_date_string


!-----------------------------------------------------------------------
!
!  clear character strings
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   char_temp   = char_blank


!-----------------------------------------------------------------------
!
!  for a ccsm tavg file, append a date/time string to the root name
!
!-----------------------------------------------------------------------

      select case (tavg_freq_iopt)
      case (freq_opt_nyear)
        char_temp = 'y'

      case (freq_opt_nmonth)
        char_temp = 'ym'

      case (freq_opt_nday)
        char_temp = 'ymd'

      case (freq_opt_nhour)
        char_temp = 'ymds'

      case (freq_opt_nsecond)
        char_temp = 'ymds'

      case (freq_opt_nstep)
        char_temp = 'ymds'
 
      case default
        char_temp = 'ymds'
      end select

!-----------------------------------------------------------------------
!     override tavg_freq_iopt if this is a tavg restart file
!-----------------------------------------------------------------------

      if (present (date_string) ) then
         char_temp = trim(date_string)
      endif

      call ccsm_date_stamp (ccsm_date_string, char_temp)
 
      file_suffix = trim(ccsm_date_string)

 
!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_suffix_ccsm

!***********************************************************************
!BOP
! !IROUTINE: tavg_create_outfile_ccsm
! !INTERFACE:

 subroutine tavg_create_outfile_ccsm(tavg_outfile_orig,tavg_outfile)

! !DESCRIPTION:
!    Modifies tavg_output to conform to the ccsm3 requirements.
!    CCSM requires restart files for tavg history files to 
!    be of the form pathname/casename.rh.datestring, and
!    pathname must originally contain the string hist. This routine
!    changes the string hist to rest.
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   character (char_len), intent(in) :: &
      tavg_outfile_orig           

! !INPUT/OUTPUT PARAMETERS:

   character (char_len), intent(inout) :: &
      tavg_outfile           

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   integer (i4) :: &
      iii           ! local index

   character (char_len) :: &
      tavg_outfile_temp     ! temp tavg output filename


!-----------------------------------------------------------------------
!    store original filename in tavg_outfile_temp
!-----------------------------------------------------------------------
     tavg_outfile_temp = trim(tavg_outfile_orig)

     iii = len_trim(tavg_outfile_temp)

!-----------------------------------------------------------------------
!    error checking
!-----------------------------------------------------------------------
     if (iii .lt. 1) &
       call exit_POP (sigAbort, &
       '(tavg_create_outfile_ccsm) error forming tavg_outfile_temp')

!-----------------------------------------------------------------------
!    replace the final character ('r') with the string 'rh'
!-----------------------------------------------------------------------
     tavg_outfile_temp(iii:iii+1) = 'rh'

!-----------------------------------------------------------------------
!    replace the string 'hist' in the pathname with 'rest'
!-----------------------------------------------------------------------
     hist_search: do iii=1,len_trim(tavg_outfile_temp)
       if (tavg_outfile_temp(iii:iii+3) == 'hist') then
           tavg_outfile_temp(iii:iii+3) = 'rest'
           exit hist_search
       endif
     enddo hist_search

     tavg_outfile = trim(tavg_outfile_temp)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_create_outfile_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_count_contents_ccsm
! !INTERFACE:

 subroutine tavg_count_contents_ccsm (num_fields,tavg_contents)
 
! !DESCRIPTION:
!  This routine counts the number of lines in the tavg_contents file
!
! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:
 
  character (char_len),intent(in)   :: tavg_contents
 
! !OUTPUT PARAMETERS:
 
  integer   (int_kind), intent(out) :: num_fields
 
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  integer   (int_kind)  :: nu
  integer   (int_kind)  :: ios
  integer   (int_kind)  :: mt
  integer   (int_kind)  :: grid_code_int
  character (char_len)  :: file_line
 

!---------------------------------------------------------------------
! get unit and open contents file
!---------------------------------------------------------------------
  call get_unit(nu)
  if (my_task == master_task) then
     open(nu, file=tavg_contents, status='old')
  endif

  num_fields = 0

  count_loop: do
 
!---------------------------------------------------------------------
!  read line from file, checking for end-of-file
!---------------------------------------------------------------------
   if (my_task == master_task) then
      read (nu,'(A)',iostat=ios) file_line
   endif
   call broadcast_scalar(ios, master_task)
   if (ios < 0) exit count_loop
   call broadcast_scalar(file_line, master_task)
 
!---------------------------------------------------------------------
!  increment 
!---------------------------------------------------------------------
    num_fields = num_fields + 1

  enddo count_loop
 
  if (my_task == master_task) then
     write(stdout,*) 'There are ',num_fields,  &
                    ' tavg fields requested via tavg_contents.'
     call POP_IOUnitsFlush(POP_stdout)
  endif

!---------------------------------------------------------------------
! close file and release unit
!---------------------------------------------------------------------
 close(nu)
 call release_unit(nu)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_count_contents_ccsm

!***********************************************************************
!BOP
! !IROUTINE: tavg_add_attrib_file_ccsm
! !INTERFACE:

 subroutine tavg_add_attrib_file_ccsm (tavg_file_desc)
 
! !DESCRIPTION:
!  This routine adds attributes to the ccsm tavg history file
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor
 
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
 
 character (char_len) ::   &
    start_time,    & 
    current_date,  &
    current_time,  &
    cell_methods,  &
    calendar


 call add_attrib_file(tavg_file_desc, 'contents', 'Diagnostic and Prognostic Variables')
 call add_attrib_file(tavg_file_desc, 'source', 'CCSM POP2, the CCSM Ocean Component')
 call add_attrib_file(tavg_file_desc, 'revision', &
   '$Id$')

 if (allow_leapyear) then
    write(calendar,'(a,i5,a,i5,a)') &
       ' Leap years allowed. Normal years have',days_in_norm_year, &
       ' days. Leap years have ' ,days_in_leap_year, ' days.'
 else
    write(calendar,'(a,i5,a)')'All years have exactly', days_in_norm_year, ' days.'
 endif
 call add_attrib_file(tavg_file_desc, 'calendar', trim(calendar))

 
 call date_and_time (date=current_date,time=current_time)
 start_time = char_blank
 write(start_time,1000) current_date(1:4), current_date(5:6),  &
                        current_date(7:8), current_time(1:2),  &
                        current_time(3:4), current_time(5:8)
 call add_attrib_file(tavg_file_desc, 'start_time', trim(start_time))


 cell_methods = char_blank
 cell_methods = 'cell_methods = time: mean ==> the variable values ' /&
     &/  'are averaged over the time interval between the previous '      /&
     &/  'time coordinate and the current one.          '                 /&
     &/  'cell_methods  absent  ==> the variable values '                 /&
     &/  'are at the time given by the current time coordinate. '
 call add_attrib_file(tavg_file_desc, 'cell_methods', trim(cell_methods))

1000  format('This dataset was created on ',a,'-',a,'-',a,' at ',a,':',a,':',a)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_add_attrib_file_ccsm 
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_add_attrib_io_field_ccsm
! !INTERFACE:

 subroutine tavg_add_attrib_io_field_ccsm (tavg_field,nfield)
 
! !DESCRIPTION:
!  This routine adds ccsm-required attributes to ccsm tavg fields
!  cell_methods, scale_factor, and _FillValue
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
 integer (i4), intent(in) ::  &
    nfield

! !INPUT/OUTPUT PARAMETERS:
 type (io_field_desc), intent(inout) ::  &
    tavg_field    ! IO file descriptor
 
!EOP
!BOC

    call add_attrib_io_field(tavg_field, 'cell_methods', 'time: mean')

    if (avail_tavg_fields(nfield)%scale_factor /= undefined_nf) then
      call add_attrib_io_field(tavg_field, 'scale_factor',avail_tavg_fields(nfield)%scale_factor)
    endif

    call add_attrib_io_field(tavg_field,'_FillValue',avail_tavg_fields(nfield)%fill_value )
    call add_attrib_io_field(tavg_field,'missing_value',avail_tavg_fields(nfield)%missing_value )


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_add_attrib_io_field_ccsm 

!***********************************************************************
!BOP
! !IROUTINE: tavg_define_coord_vars_ccsm
! !INTERFACE:

 subroutine tavg_define_coord_vars_ccsm(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF coordinates that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real (r4), dimension(km) ::  &
      ZT_R4,      &! single precision array
      ZW_R4        ! single precision array

   real (r4), dimension(0:km) ::  &
      MOC_Z_R4

   real (r4), dimension(1000) ::  &
      LAT_AUX_GRID_R4

   integer (int_kind) ::  &
      ii, n, ndims   

   save  



   ii=0

   !*** z_t
   ii=ii+1
   ZT_R4 = zt 
   ccsm_coordinates(ii) = construct_io_field('z_t',zt_dim,                    &
                         long_name='depth from surface to midpoint of layer', &
                         units    ='centimeters',                             &
                         r1d_array =ZT_R4)

   call add_attrib_io_field(ccsm_coordinates(ii), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_min', ZT_R4(1))
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_max', ZT_R4(km))

   !*** z_t
   ii=ii+1
   ZT_150m_R4 = zt(1:zt_150m_levs)
   ccsm_coordinates(ii) = construct_io_field('z_t_150m',zt_150m_dim,          &
                         long_name='depth from surface to midpoint of layer', &
                         units    ='centimeters',                             &
                         r1d_array =ZT_150m_R4)

   call add_attrib_io_field(ccsm_coordinates(ii), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_min', ZT_150m_R4(1))
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_max', ZT_150m_R4(zt_150m_levs))

   !*** z_w
   ii=ii+1
   ZW_R4(1) = c0 
   ZW_R4(2:km) = zw(1:km-1)
   ccsm_coordinates(ii) = construct_io_field('z_w',zw_dim,                    &
                         long_name='depth from surface to top of layer',      &
                         units    ='centimeters',                             &
                         r1d_array =ZW_R4)

   call add_attrib_io_field(ccsm_coordinates(ii), 'positive', 'down')
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_min', ZW_R4(1 ))
   call add_attrib_io_field(ccsm_coordinates(ii), 'valid_max', ZW_R4(km))


   if (moc .or. n_heat_trans  .or. n_salt_trans) then
     !*** lat_aux_grid
     ii=ii+1

     if (n_lat_aux_grid+1 > 1000) &
        call exit_POP(sigAbort,'Must increase dimension of LAT_AUX_GRID_R4')

     LAT_AUX_GRID_R4 = 0
     LAT_AUX_GRID_R4(1:n_lat_aux_grid+1) = lat_aux_edge(1:n_lat_aux_grid+1)
     ccsm_coordinates(ii) = construct_io_field('lat_aux_grid',lat_aux_grid_dim, &
                           long_name='latitude grid for transport diagnostics', &
                           units    ='degrees_north',                           &
                           r1d_array =LAT_AUX_GRID_R4(1:n_lat_aux_grid+1))
     call add_attrib_io_field(ccsm_coordinates(ii), 'valid_min', LAT_AUX_GRID_R4(1 ))
     call add_attrib_io_field(ccsm_coordinates(ii), 'valid_max', LAT_AUX_GRID_R4(n_lat_aux_grid+1))
   endif

   if (moc) then
     !*** moc_z
     ii=ii+1

     MOC_Z_R4(0) = c0 
     MOC_Z_R4(1:km) = zw(1:km)
     ccsm_coordinates(ii) = construct_io_field('moc_z',moc_z_dim,               &
                           long_name='depth from surface to top of layer',      &
                           units    ='centimeters',                             &
                           r1d_array =MOC_Z_R4)

     call add_attrib_io_field(ccsm_coordinates(ii), 'positive', 'down')
     call add_attrib_io_field(ccsm_coordinates(ii), 'valid_min', MOC_Z_R4(0 ))
     call add_attrib_io_field(ccsm_coordinates(ii), 'valid_max', MOC_Z_R4(km))

   endif ! moc


   ! after all coordinates are defined, define the total number of coordinates
   num_ccsm_coordinates = ii
 
   if (num_ccsm_coordinates > max_num_ccsm_coordinates) &
      call exit_POP(sigAbort,'ERROR reset max_num_ccsm_coordinates')

   do n=1,num_ccsm_coordinates
    call data_set (tavg_file_desc, 'define', ccsm_coordinates(n))
   enddo


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_coord_vars_ccsm
 
!***********************************************************************
!BOP
! !IROUTINE: tavg_define_time_invar_ccsm
! !INTERFACE:

 subroutine tavg_define_time_invar_ccsm(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF time-invariant variables that are 
!  used in the ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ii, n

   real (r4), dimension(km)   ::  &
      DZ_R4   

   real (r4), dimension(0:km-1) ::  &
      DZW_R4

   real (r4), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      TLON_DEG, TLAT_DEG, ULON_DEG, ULAT_DEG

   real (r4)          :: missing_value
   integer (int_kind) :: missing_value_i

   save

   missing_value   = undefined_nf_r4
   missing_value_i = undefined_nf_int

   ii=0

   !*** dz
   ii=ii+1
   DZ_R4 = dz 
   ccsm_time_invar(ii) = construct_io_field('dz',zt_dim,                    &
                                long_name='thickness of layer k',           &
                                units    ='centimeters',                    &
                                r1d_array =DZ_R4)
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** dzw
   ii=ii+1
   DZW_R4(0:) = dzw(0:km-1)
   ccsm_time_invar(ii) = construct_io_field('dzw',zw_dim,                   &
                                long_name='midpoint of k to midpoint of k+1',      &
                                units    ='centimeters',                           &
                                r1d_array =DZW_R4)
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** ULONG
   ii=ii+1

    ULON_DEG = ULON*radian
    ccsm_time_invar(ii) = construct_io_field(  &
        'ULONG', dim1=i_dim, dim2=j_dim,              &
         long_name='array of u-grid longitudes',      &
         units    ='degrees_east',                    &
         r2d_array =ULON_DEG(:,:,:) )

   !*** ULAT
   ii=ii+1

    ULAT_DEG = ULAT*radian
    ccsm_time_invar(ii) = construct_io_field(  &
        'ULAT', dim1=i_dim, dim2=j_dim,               &
         long_name='array of u-grid latitudes',       &
         units    ='degrees_north',                   &
         r2d_array =ULAT_DEG(:,:,:) )

   !*** TLONG
   ii=ii+1

    TLON_DEG = TLON*radian
    ccsm_time_invar(ii) = construct_io_field(  &
        'TLONG', dim1=i_dim, dim2=j_dim,              &
         long_name='array of t-grid longitudes',      &
         units    ='degrees_east',                    &
         r2d_array =TLON_DEG(:,:,:) )

   !*** TLAT
   ii=ii+1

    TLAT_DEG = TLAT*radian
    ccsm_time_invar(ii) = construct_io_field(  &
        'TLAT', dim1=i_dim, dim2=j_dim,               &
         long_name='array of t-grid latitudes',       &
         units    ='degrees_north',                   &
         r2d_array =TLAT_DEG(:,:,:) )

   !*** KMT
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'KMT', dim1=i_dim, dim2=j_dim,                        &
         long_name='k Index of Deepest Grid Cell on T Grid',  &
         coordinates = "TLONG TLAT",                          &
         i2d_array =KMT(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value_i )

   !*** KMU
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'KMU', dim1=i_dim, dim2=j_dim,                        &
         long_name='k Index of Deepest Grid Cell on U Grid',  &
         coordinates = "ULONG ULAT",                          &
         i2d_array =KMU(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value_i )


   !*** REGION_MASK
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'REGION_MASK', dim1=i_dim, dim2=j_dim,                &
         long_name='basin index number (signed integers)',    &
         coordinates = "TLONG TLAT",                          &
         i2d_array =REGION_MASK(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value_i )

   !*** UAREA
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'UAREA', dim1=i_dim, dim2=j_dim,                      &
         long_name='area of U cells',                         &
         units    ='centimeter^2',                            &
         coordinates = "ULONG ULAT",                          &
         d2d_array =UAREA(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** TAREA
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'TAREA', dim1=i_dim, dim2=j_dim,                      &
         long_name='area of T cells',                         &
         units    ='centimeter^2',                            &
         coordinates = "TLONG TLAT",                          &
         d2d_array =TAREA(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HU
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HU', dim1=i_dim, dim2=j_dim,                         &
         long_name='ocean depth at U points',                 &
         units='centimeter',                                  &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HU(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HT
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HT', dim1=i_dim, dim2=j_dim,                         &
         long_name='ocean depth at T points',                 &
         units='centimeter',                                  &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HT(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** DXU
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'DXU', dim1=i_dim, dim2=j_dim,                        &
         long_name='x-spacing centered at U points',          &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =DXU(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** DYU
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'DYU', dim1=i_dim, dim2=j_dim,                        &
         long_name='y-spacing centered at U points',          &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =DYU(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** DXT
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'DXT', dim1=i_dim, dim2=j_dim,                        &
         long_name='x-spacing centered at T points',          &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =DXT(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** DYT
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'DYT', dim1=i_dim, dim2=j_dim,                        &
         long_name='y-spacing centered at T points',          &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =DYT(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HTN
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HTN', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on North sides of T cell',    &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HTN(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HTE
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HTE', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on East sides of T cell',     &
         units='centimeters',                                 &
         coordinates = "TLONG TLAT",                          &
         d2d_array =HTE(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HUS
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HUS', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on South sides of U cell',    &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HUS(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** HUW
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'HUW', dim1=i_dim, dim2=j_dim,                        &
         long_name='cell widths on West sides of U cell',     &
         units='centimeters',                                 &
         coordinates = "ULONG ULAT",                          &
         d2d_array =HUW(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** ANGLE
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(          &
        'ANGLE', dim1=i_dim, dim2=j_dim,                      &
         long_name='angle grid makes with latitude line',     &
         units='radians',                                     &
         coordinates = "ULONG ULAT",                          &
         d2d_array =ANGLE(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )

   !*** ANGLET
   ii=ii+1

    ccsm_time_invar(ii) = construct_io_field(               &
        'ANGLET', dim1=i_dim, dim2=j_dim,                          &
         long_name='angle grid makes with latitude line on T grid',&
         units='radians',                                          &
         coordinates = "TLONG TLAT",                               &
         d2d_array =ANGLET(:,:,:) )
   call add_attrib_io_field(ccsm_time_invar(ii),'missing_value',missing_value )



   ! after all time-invariant arrays are defined, define the total number
   num_ccsm_time_invar = ii

   if (num_ccsm_time_invar  > max_num_ccsm_time_invar) &
      call exit_POP(sigAbort,'ERROR reset max_num_ccsm_time_invar')

   do n=1,num_ccsm_time_invar
    call data_set (tavg_file_desc, 'define', ccsm_time_invar(n))
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_time_invar_ccsm

!***********************************************************************
!BOP
! !IROUTINE: tavg_define_scalars_ccsm
! !INTERFACE:

 subroutine tavg_define_scalars_ccsm(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF scalars that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

   real (r8) ::  &
      d0d_days_in_norm_year,  &
      d0d_days_in_leap_year,  &
      d0d_nsurface_t,         &
      d0d_nsurface_u

   integer (int_kind) ::  &
      ii,n

   save  


   ii=0

   !*** days_in_norm_year
   d0d_days_in_norm_year = days_in_norm_year
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('days_in_norm_year',      &
                         long_name='Calendar Length',              &
                         units    ='days',                         &
                         d0d_array =d0d_days_in_norm_year)

   !*** grav
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('grav',                   &
                         long_name='Acceleration Due to Gravity',  &
                         units    ='centimeter/s^2',               &
                         d0d_array =grav)

   !*** omega
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('omega',                   &
                         long_name='Earths Angular Velocity',       &
                         units    ='1/second',                      &
                         d0d_array =omega)

   !*** radius
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('radius',   &
                         long_name='Earths Radius',  &
                         units    ='centimeters',    &
                         d0d_array =radius)


   !*** cp_sw
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('cp_sw',                 &
                         long_name='Specific Heat of Sea Water',  &
                         units    ='erg/g/K',                     &
                         d0d_array =cp_sw)

   !*** sound
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('sound',                 &
                         long_name='Speed of Sound',              &
                         units    ='centimeter/s',                &
                         d0d_array =sound)

   !*** vonkar
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('vonkar',                &
                         long_name='von Karman Constant',         &
                         d0d_array =vonkar)

   !*** cp_air
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('cp_air',                &
                         long_name='Heat Capacity of Air',        &
                         units    ='joule/kg/degK',               &
                         d0d_array =cp_air)

   !*** rho_air
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('rho_air',               &
                         long_name='Ambient Air Density',         &
                         units    ='kg/m^3',                      &
                         d0d_array =rho_air)

   !*** rho_sw
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('rho_sw',                &
                         long_name='Density of Sea Water',        &
                         units    ='gram/centimeter^3',           &
                         d0d_array =rho_sw)

   !*** rho_fw
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('rho_fw',                &
                         long_name='Density of Fresh Water',      &
                         units    ='gram/centimeter^3',           &
                         d0d_array =rho_fw)

   !*** stefan_boltzmann
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('stefan_boltzmann',      &
                         long_name='Stefan-Boltzmann Constant',    &
                         units    ='watt/m^2/degK^4',             &
                         d0d_array =stefan_boltzmann)

   !*** latent_heat_vapor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('latent_heat_vapor',     &
                         long_name='Latent Heat of Vaporization', &
                         units    ='erg/g',                       &
                         d0d_array =latent_heat_vapor)

   !*** latent_heat_fusion
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('latent_heat_fusion',    &
                         long_name='Latent Heat of Fusion',       &
                         units    ='erg/g',                       &
                         d0d_array =latent_heat_fusion)

   !*** ocn_ref_salinity
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('ocn_ref_salinity',      &
                         long_name='Ocean Reference Salinity',    &
                         units    ='g/kg',                        &
                         d0d_array =ocn_ref_salinity)

   !*** sea_ice_salinity
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('sea_ice_salinity',      &
                         long_name='Salinity of Sea Ice',         &
                         units    ='g/kg',                        &
                         d0d_array =sea_ice_salinity)

   !*** T0_Kelvin
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('T0_Kelvin',             &
                         long_name='Zero Point for Celsius',      &
                         units    ='degK',                        &
                         d0d_array =T0_Kelvin)

   !*** salt_to_ppt
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('salt_to_ppt',           &
                         long_name='Convert Salt in gram/gram to g/kg',&
                         d0d_array =salt_to_ppt)

   !*** ppt_to_salt
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('ppt_to_salt',           &
                         long_name='Convert Salt in g/kg to gram/gram',&
                         d0d_array =ppt_to_salt)
   !*** mass_to_Sv
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('mass_to_Sv',           &
                         long_name='Convert Mass Flux to Sverdrups',&
                         d0d_array =mass_to_Sv)

   !*** heat_to_PW
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('heat_to_PW',           &
                         long_name='Convert Heat Flux to Petawatts', &
                         d0d_array =heat_to_PW)

   !*** salt_to_Svppt
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('salt_to_Svppt',           &
                         long_name='Convert Salt Flux to Sverdrups*g/kg', & 
                         d0d_array =salt_to_Svppt)
   !*** salt_to_mmday
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('salt_to_mmday',           &
                         long_name='Convert Salt to Water (millimeters/day)', &
                         d0d_array =salt_to_mmday)

   !*** momentum_factor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('momentum_factor',           &
                         long_name='Convert Windstress to Velocity Flux', &
                         d0d_array =momentum_factor)

   !*** hflux_factor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('hflux_factor',           &
        long_name='Convert Heat and Solar Flux to Temperature Flux',&
                         d0d_array =hflux_factor)

   !*** fwflux_factor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('fwflux_factor',           &
  long_name='Convert Net Fresh Water Flux to Salt Flux (in model units)', &
                         d0d_array =fwflux_factor)
 
   !*** salinity_factor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('salinity_factor',           &
  long_name=' ', &
                         d0d_array =salinity_factor)
 
   !*** sflux_factor
   ii=ii+1
   ccsm_scalars(ii) = construct_io_field('sflux_factor',           &
  long_name='Convert Salt Flux to Salt Flux (in model units)', &
                         d0d_array =sflux_factor)
 
   !*** nsurface_t
   ii=ii+1
   d0d_nsurface_t = nsurface_t
   ccsm_scalars(ii) = construct_io_field('nsurface_t',           &
  long_name='Number of Ocean T Points at Surface', &
                         d0d_array =d0d_nsurface_t)
 
   !*** nsurface_u
   ii=ii+1
   d0d_nsurface_u = nsurface_u
   ccsm_scalars(ii) = construct_io_field('nsurface_u',           &
  long_name='Number of Ocean U Points at Surface', &
                         d0d_array =d0d_nsurface_u)
 
 
   ! after all scalars are defined, define the total number of scalars
   num_ccsm_scalars = ii

   if (num_ccsm_scalars > max_num_ccsm_scalars) &
      call exit_POP(sigAbort,'ERROR reset num_ccsm_scalars')


   do n=1,num_ccsm_scalars
    call data_set (tavg_file_desc, 'define', ccsm_scalars(n))
   enddo


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_scalars_ccsm

!***********************************************************************
!BOP
! !IROUTINE: tavg_define_labels_ccsm
! !INTERFACE:

 subroutine tavg_define_labels_ccsm(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF labels that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      id, ii, n, ndims
 
   integer (i4) ::  &
      nstd_field_id
         

   save  

!-----------------------------------------------------------------------
!
! Note that avail_tavg_labels is type tavg_field_desc_ccsm, not io_field_desc
!
!-----------------------------------------------------------------------
   ii = 0
   num_avail_tavg_labels = 0

   if (moc) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=moc_comp_dim
     call define_tavg_field(id, 'moc_components',ndims,            &
                            long_name='MOC component names',       &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (n_heat_trans .or. n_salt_trans) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=transport_comp_dim
     call define_tavg_field(id, 'transport_components',ndims,      &
                            long_name='T,S transport components',  &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (moc .or. n_heat_trans .or. n_salt_trans) then
     ii = ii + 1
     ndims = 2
     io_dims_labels(1,ii)=nchar_dim
     io_dims_labels(2,ii)=transport_reg_dim
     call define_tavg_field(id, 'transport_regions',ndims,         &
                            long_name='regions for all transport diagnostics',  &
                            nftype='char',                         &
                            nstd_fields=avail_tavg_labels,         &
                            num_nstd_fields=num_avail_tavg_labels, &
                            max_nstd_fields=max_avail_tavg_labels  )
   endif

   if (ii /= num_avail_tavg_labels) then
      write (stdout,*) 'ii = ', ii
      write (stdout,*) 'num_avail_tavg_labels = ', num_avail_tavg_labels
      call exit_POP(sigAbort,  &
         '(tavg_define_labels_ccsm) mismatch in number of labels')
   endif

   do n=1,num_avail_tavg_labels
    call data_set_nstd_ccsm (                          &
                     tavg_file_desc, 'define',         &
                     nstd_field_id,                    &
                     avail_tavg_labels(n)%ndims,       &
                     io_dims_labels(:,n),              &
            short_name=avail_tavg_labels(n)%short_name,&
             long_name=avail_tavg_labels(n)%long_name, &
                 units=avail_tavg_labels(n)%units,     &
                nftype=avail_tavg_labels(n)%nftype     )

    avail_tavg_labels_id(n) = nstd_field_id

   enddo


!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_labels_ccsm
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_define_time_ccsm
! !INTERFACE:

 subroutine tavg_define_time_ccsm(tavg_file_desc)

! !DESCRIPTION:
!  This routine defines the netCDF time variables that are used in the
!  ccsm netCDF output tavg files
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:
   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(1) ::  &
      TIME1D  

   character(char_len) ::  &
      nftype              

   integer (int_kind) ::  &
      ndims

   save  

   !*** time
   TIME1D(1)=tday00
   time_coordinate(1) = construct_io_field('time',time_dim,               &
                     long_name ='time',                                   &
                     units     ='days since 0000-01-01 00:00:00',         &
                     d1d_array = TIME1D                                   )

   call add_attrib_io_field(time_coordinate(1), 'bounds', 'time_bound')
   call add_attrib_io_field(time_coordinate(1), 'calendar', 'noleap')
   call data_set (tavg_file_desc, 'define', time_coordinate(1))


   !*** time_bound
   !*** NB: call data_set_nstd_ccsm
   time_bound_dims(1) = d2_dim
   time_bound_dims(2) = time_dim
   ndims = 2
   nftype = 'double'

   call data_set_nstd_ccsm (tavg_file_desc, 'define',               &
                            time_bound_id,                          & 
                            ndims, time_bound_dims,                 &
                short_name='time_bound',                            &
                 long_name='boundaries for time-averaging interval',&
                     units='days since 0000-01-01 00:00:00',        &
                    nftype=nftype                                   )

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_define_time_ccsm
 

!***********************************************************************
!BOP
! !IROUTINE: tavg_write_vars_ccsm
! !INTERFACE:

 subroutine tavg_write_vars_ccsm (tavg_file_desc, nvars, ccsm_vars) 

! !DESCRIPTION:
!  This routine writes the ccsm variables (coordinates, scalars, and 
!  time-independent variables) to the ccsm version of the netCDF 
!  tavg output files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) ::  &
      nvars

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent(inout) ::  &
      tavg_file_desc    ! IO file descriptor

   type (io_field_desc), dimension(nvars), intent(inout)  :: &
      ccsm_vars        

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      n, ndims             ! local index

   do n=1,nvars
      call data_set (tavg_file_desc, 'write', ccsm_vars(n))
      call destroy_io_field(ccsm_vars(n))
   enddo

!-----------------------------------------------------------------------
!EOC
 end subroutine tavg_write_vars_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_write_vars_nstd_ccsm
! !INTERFACE:

 subroutine tavg_write_vars_nstd_ccsm(tavg_file_desc) 

! !DESCRIPTION:
!  This routine writes the nonstandard ccsm variables to the
!  ccsm netCDF output tavg file.  These variables exist in
!  a variable of type (tavg_field_desc_ccsm), but do not have the
!  corresponding data bundled together via the function
!  construct_io_field. Therefore, this routine associates the
!  variable's data with the varible's definition, and sends
!  the information in separate pieces to be written.
!  Calls data_set_nstd_ccsm
!  
! 
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent(inout) :: tavg_file_desc    ! IO file descriptor

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   character (char_len) ::  &
      nftype

   integer (int_kind), parameter ::  &
      max_writes = 1

   integer (int_kind) ::  &
      n, ndims, nnn,      &! local index
      num_writes,         &
      ii, indx

   integer (i4) ::  &
      id

   real (r8), dimension(2,max_writes) ::  &
      data_2d_r8

   character (char_len), dimension(30) :: &
      data_1d_ch

   type(io_dim), dimension(2) ::  & ! local array
       io_dims


   num_writes=1  ! place-holder until this model supports multiple times/file


   !*** time bound 
   upper_time_bound = tday00
   ndims=2
   data_2d_r8(1,num_writes) = lower_time_bound
   data_2d_r8(2,num_writes) = upper_time_bound
   nftype = 'double'
   implied_time_dim = .false.
 
   call data_set_nstd_ccsm (tavg_file_desc, 'write',  &
                            time_bound_id,            &
                            ndims, time_bound_dims,   &
                            nftype,                   &
           implied_time_dim=implied_time_dim,         &
                data_2d_r8 = data_2d_r8)

   lower_time_bound = tday00 


   !*** moc-related variables
   if (moc) then
 
      !*** determine index of moc_components

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'moc_components')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Eulerian Mean'
        if (n_moc_comp >= 2) &
          data_1d_ch(2) = 'Eddy-Induced (bolus)'
        if (n_moc_comp >= 3) & 
          data_1d_ch(3) = 'Submeso'
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        implied_time_dim = .false.
 

      call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                               id, ndims, io_dims, nftype, &
              implied_time_dim=implied_time_dim,           &
                   data_1d_ch = data_1d_ch                 )
       endif

      !*** determine index of transport_regions

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'transport_regions')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Global Ocean - Marginal Seas'
        if (n_transport_reg > 1 .and. nreg2_transport >=  1) then
          data_1d_ch(2) = trim(transport_region_info(1)%name)
          do nnn=2,nreg2_transport
             data_1d_ch(2)=trim(data_1d_ch(2)) /&
                                                 &/ ' + ' /&
                                                 &/trim(transport_region_info(nnn)%name)
          enddo
        endif
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        implied_time_dim = .false.
 

         call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                                  id, ndims, io_dims, nftype, &
                 implied_time_dim=implied_time_dim,           &
                      data_1d_ch = data_1d_ch                 )
       endif


       implied_time_dim = .true.
       call write_nstd_netcdf(                          &
          tavg_file_desc,moc_id,1,5,io_dims_nstd_ccsm(:,1),&
          'float',                                         &  ! <-- generalize later
          implied_time_dim=implied_time_dim,               &
          indata_4d_r4=TAVG_MOC_G)
   endif

   !*** tranpsort variables
   if (n_heat_trans .or. n_salt_trans) then
 
      !*** determine index of transport_components

      indx=0
      do ii=1,num_avail_tavg_labels
      if (trim(avail_tavg_labels(ii)%short_name) == 'transport_components')   &
          indx = ii
      enddo

     
      if (indx /= 0) then
        ndims=2
        data_1d_ch(1) = 'Total'
        data_1d_ch(2) = 'Eulerian-Mean Advection'
        if (registry_match('init_gm')) then 
          data_1d_ch(3) = 'Eddy-Induced Advection (bolus) + Diffusion'
        else
          data_1d_ch(3) = 'Diffusion'
        endif
        if ( n_transport_comp >= 4 ) & 
          data_1d_ch(4) = 'Eddy-Induced (bolus) Advection'
        if ( n_transport_comp >= 5 ) & 
          data_1d_ch(5) = 'Submeso Advection'
        io_dims(:) = io_dims_labels(:,indx)
        id = avail_tavg_labels_id(indx)
        nftype = 'char'
        implied_time_dim = .false.
 

      call data_set_nstd_ccsm (tavg_file_desc, 'write',    &
                               id, ndims, io_dims, nftype, &
              implied_time_dim=implied_time_dim,           &
                   data_1d_ch = data_1d_ch                 )
       endif


       implied_time_dim = .true.
       call write_nstd_netcdf(                          &
          tavg_file_desc,moc_id,1,5,io_dims_nstd_ccsm(:,1),&
          'float',                                         &  ! <-- generalize later
          implied_time_dim=implied_time_dim,               &
          indata_4d_r4=TAVG_MOC_G)
   endif

   if (n_heat_trans) then
       !*** N_HEAT
       implied_time_dim = .true.
       call write_nstd_netcdf(                          &
          tavg_file_desc,n_heat_id,1,4,io_dims_nstd_ccsm(:,2),& !<-- generalize later
          'float',                                         &  ! <-- generalize later
          implied_time_dim=implied_time_dim,               &
          indata_3d_r4=TAVG_N_HEAT_TRANS_G)
   endif

   if (n_salt_trans) then
       !*** N_SALT
       implied_time_dim = .true.
       call write_nstd_netcdf(                          &
          tavg_file_desc,n_salt_id,1,4,io_dims_nstd_ccsm(:,3),& !<-- generalize later
          'float',                                         &  ! <-- generalize later
          implied_time_dim=implied_time_dim,               &
          indata_3d_r4=TAVG_N_SALT_TRANS_G)
   endif


!-----------------------------------------------------------------------
!EOC
 
 end subroutine tavg_write_vars_nstd_ccsm


!***********************************************************************
!BOP
! !IROUTINE:  tavg_mask
! !INTERFACE:

 subroutine tavg_mask(ARRAY,tavg_field,k)

! !DESCRIPTION:
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (rtavg), dimension (:,:,:), intent(inout) ::  &
      ARRAY

   type (tavg_field_desc_ccsm), intent(in) ::  &
      tavg_field

   integer (int_kind), intent(in) ::  &
      k

   type (block) ::        &
      this_block          ! block information for current block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) :: iblock, i, j


   select case (trim(tavg_field%grid_loc(2:3)))

     case('11')
          
       do iblock=1,nblocks_clinic
        ARRAY(:,:,iblock) = merge (tavg_field%fill_value,ARRAY(:,:,iblock),  &
                              k > KMT(:,:,iblock))
       enddo ! iblock
 
     case('22')

       do iblock=1,nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
             if (MASK_22(i,j,k,iblock)) ARRAY(i,j,iblock) = tavg_field%fill_value
         enddo ! i
         enddo ! j
       enddo ! iblock
 
     case default
         !*** if field_desc is not on standard tracer or 
         !    velocity grid, never mind

   end select


!EOC
 
 end subroutine tavg_mask


!***********************************************************************
!BOP
! !IROUTINE: tavg_bsf_ccsm
! !INTERFACE:

 subroutine tavg_bsf_ccsm

! !DESCRIPTION:
!
!  Compute barotropic stream function diagnostically from other tavg
!  quantities
!
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode

   real (r8), dimension (nx_block, ny_block,max_blocks_clinic) ::  &
      WORK1,WORK2,WORK3,  &
      PSI_T, PSI_U

   integer (int_kind) ::   &
      tavg_id_SU,          &
      tavg_id_SV,          &
      tavg_id_BSF,         &
      tavg_loc_SU,         &
      tavg_loc_SV,         &
      tavg_loc_BSF

   integer (int_kind) ::   &
      iblock,              &
      i,ii,                &
      j,jj

   type (block) ::        &
      this_block          ! block information for current block

!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (.not. ldiag_bsf) return
 
   !*** start bsf timer
   call timer_start(timer_tavg_ccsm_diags_bsf)

   !*** return if attempting to compute every nstep timesteps
   if (tavg_freq_iopt == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: BSF diagnostic is not computed if tavg_freq_iopt == freq_opt_nstep'
        call POP_IOUnitsFlush(POP_stdout)
     endif
     return
   endif


   !*** zero out location identifiers
   tavg_loc_SU = 0; tavg_loc_SV  = 0; tavg_loc_BSF = 0

   !*** determine the tavg ids for tavg fields required by this module
   tavg_id_SU  = tavg_id('SU')
   tavg_id_SV  = tavg_id('SV')
   tavg_id_BSF = tavg_id('BSF')

   if (tavg_id_SU == 0 .or. tavg_id_SV == 0 .or. tavg_id_BSF == 0) then
     !*** do not abort; write warning message and proceed
     call document ('write_tavg', 'WARNING: cannot compute BSF diagnostic')
   else
 
     tavg_loc_SU  = avail_tavg_fields(tavg_id_SU)%buf_loc
     tavg_loc_SV  = avail_tavg_fields(tavg_id_SV)%buf_loc
     tavg_loc_BSF = avail_tavg_fields(tavg_id_BSF)%buf_loc

     !$OMP PARALLEL DO PRIVATE (this_block,iblock)
     do iblock=1,nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)

       WORK2(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_SU)
       WORK3(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_SV)
       call zcurl (1,WORK1(:,:,iblock),WORK2(:,:,iblock),  &
                    WORK3(:,:,iblock),this_block)
     end do
     !$OMP END PARALLEL DO
       
     WORK2 = c0
     call pcg_diag_bsf_solver (WORK2,WORK1, errorCode)
     if (errorCode /= POP_Success) then
        call POP_ErrorSet(errorCode, &
           'tavg_bsf_ccsm: error in pcg_diag_bsf')
        call exit_POP(sigAbort,'tavg_bsf_ccsm: error in pcg_diag_bsf')
        return
     endif
 
     do iblock=1,nblocks_clinic
     TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF) =  &
         merge(c0,WORK2(:,:,iblock), .not. CALCT(:,:,iblock))
     enddo ! iblock

     !*** convert to Sv
     TAVG_BUF_2D(:,:,:,tavg_loc_BSF) =  &
       TAVG_BUF_2D(:,:,:,tavg_loc_BSF)*mass_to_Sv

     do iblock=1,nblocks_clinic
     PSI_T(:,:,iblock)= TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF)
     enddo

     !$OMP PARALLEL DO PRIVATE (iblock, i,j,ii,jj)
     do iblock=1,nblocks_clinic

        PSI_T(:,:,iblock) = TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF)
        do j=2,ny_block-1
        do i=2,nx_block-1
          if (KMT(i,j,iblock) == 0) then
            ii_loop: do ii=i-1,i+1
            jj_loop: do jj=j-1,j+1
              if (KMT(ii,jj,iblock) /= 0) then
                PSI_T(i,j,iblock) = PSI_T(ii,jj,iblock)
                exit ii_loop
              endif
            enddo jj_loop
            enddo ii_loop
           endif
        enddo ! i
        enddo ! j
 
        call tgrid_to_ugrid (PSI_U(:,:,iblock), PSI_T(:,:,iblock),iblock)
 
        TAVG_BUF_2D(:,:,iblock,tavg_loc_BSF) = PSI_U(:,:,iblock)


     end do
     !$OMP END PARALLEL DO

   endif ! test on BSF diagnostic
 
   !*** stop bsf timer
   call timer_stop(timer_tavg_ccsm_diags_bsf)

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_bsf_ccsm

 
!***********************************************************************
!BOP
! !IROUTINE: tavg_moc_ccsm
! !INTERFACE:

 subroutine tavg_moc_ccsm

! !DESCRIPTION:
!
!  Compute meridional overturning circulation diagnostically from 
!  other tavg quantities
!
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      tavg_id_WVEL,         &
      tavg_id_VVEL,         &
      tavg_id_WISOP,        &
      tavg_id_VISOP,        &
      tavg_loc_WVEL,        &
      tavg_loc_VVEL,        &
      tavg_loc_WISOP,       &
      tavg_loc_VISOP,       &
      tavg_id_WSUBM,        &
      tavg_id_VSUBM,        &
      tavg_loc_WSUBM,       &
      tavg_loc_VSUBM 

   logical (log_kind) ::  &
      ldiag_gm_bolus,     &  ! local logical for diag_gm_bolus
      lsubmeso               ! local logical for submesoscale_mixing

   if (.not. moc) return
 
   !*** start timer
   call timer_start(timer_tavg_ccsm_diags_moc)

   ldiag_gm_bolus = .false.
   if ( registry_match('diag_gm_bolus') )  ldiag_gm_bolus = .true.

   lsubmeso = .false.
   if ( registry_match('init_submeso') )   lsubmeso = .true.

   tavg_loc_WVEL      = 0 ; tavg_loc_VVEL      = 0
   tavg_loc_WISOP     = 0 ; tavg_loc_VISOP     = 0
   tavg_loc_WSUBM     = 0 ; tavg_loc_VSUBM     = 0
 
   tavg_id_WVEL   = tavg_id('WVEL')
   tavg_id_VVEL   = tavg_id('VVEL')

   !*** error checking
   if (tavg_id_WVEL  == 0 .or. tavg_id_VVEL  == 0 ) then
     call document ('tavg_moc_ccsm', &
         'Error in moc diagnostics computations: none of the following should be zero')
     call document ('tavg_moc_ccsm', 'tavg_id_WVEL', tavg_id_WVEL)
     call document ('tavg_moc_ccsm', 'tavg_id_VVEL', tavg_id_VVEL)
     call exit_POP (SigAbort, 'Fatal error')
   endif
 
   tavg_loc_WVEL = avail_tavg_fields(tavg_id_WVEL)%buf_loc
   tavg_loc_VVEL = avail_tavg_fields(tavg_id_VVEL)%buf_loc

   !*** error checking
   if (tavg_loc_WVEL  == 0 .or. tavg_loc_VVEL  == 0 ) then
     call document ('tavg_moc_ccsm', &
         'Error in moc diagnostics computations: none of the following should be zero')
     call document ('tavg_moc_ccsm', 'tavg_loc_WVEL', tavg_loc_WVEL)
     call document ('tavg_moc_ccsm', 'tavg_loc_VVEL', tavg_loc_VVEL)
     call exit_POP (SigAbort, 'Fatal error')
   endif

   if (ldiag_gm_bolus) then
 
     tavg_id_WISOP  = tavg_id('WISOP')
     tavg_id_VISOP  = tavg_id('VISOP')
 
     !*** error checking
     if (tavg_id_WISOP  == 0 .or. tavg_id_VISOP  == 0 ) then
       call document ('tavg_moc_ccsm', &
           'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc_ccsm', 'tavg_id_WISOP',  tavg_id_WISOP)
       call document ('tavg_moc_ccsm', 'tavg_id_VISOP',  tavg_id_VISOP)
       call exit_POP (SigAbort, 'Fatal error')
     endif

     tavg_loc_WISOP = avail_tavg_fields(tavg_id_WISOP)%buf_loc
     tavg_loc_VISOP = avail_tavg_fields(tavg_id_VISOP)%buf_loc

     !*** error checking
     if (tavg_loc_WISOP  == 0 .or. tavg_loc_VISOP  == 0 ) then
       call document ('tavg_moc_ccsm', &
           'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc_ccsm', 'tavg_loc_WISOP',  tavg_loc_WISOP)
       call document ('tavg_moc_ccsm', 'tavg_loc_VISOP',  tavg_loc_VISOP)
       call exit_POP (SigAbort, 'Fatal error')
     endif

   endif ! ldiag_gm_bolus

   if (lsubmeso) then

     tavg_id_WSUBM  = tavg_id('WSUBM')
     tavg_id_VSUBM  = tavg_id('VSUBM')

     !*** error checking
     if (tavg_id_WSUBM  == 0 .or. tavg_id_VSUBM  == 0 ) then
       call document ('tavg_moc_ccsm', &
           'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc_ccsm', 'tavg_id_WSUBM',  tavg_id_WSUBM)
       call document ('tavg_moc_ccsm', 'tavg_id_VSUBM',  tavg_id_VSUBM)
       call exit_POP (SigAbort, 'Fatal error')
     endif

     tavg_loc_WSUBM = avail_tavg_fields(tavg_id_WSUBM)%buf_loc
     tavg_loc_VSUBM = avail_tavg_fields(tavg_id_VSUBM)%buf_loc

     !*** error checking
     if (tavg_loc_WSUBM  == 0 .or. tavg_loc_VSUBM  == 0 ) then
       call document ('tavg_moc_ccsm', &
           'Error in moc diagnostics computations: none of the following should be zero')
       call document ('tavg_moc_ccsm', 'tavg_loc_WSUBM',  tavg_loc_WSUBM)
       call document ('tavg_moc_ccsm', 'tavg_loc_VSUBM',  tavg_loc_VSUBM)
       call exit_POP (SigAbort, 'Fatal error')
     endif

   endif ! lsubmeso

   !*** define MOC tavg variable and dimensions
   !*** note that this call fills avail_tavg_nstd_fields

   call define_tavg_field(                               &
        tavg_MOC, 'MOC', 5,                              &
        long_name='Meridional Overturning Circulation',  &
        units='Sverdrups',                               &
        coordinates='lat_aux_grid moc_z moc_components transport_region time',&
        nftype='float'   ,                               &
        nstd_fields=avail_tavg_nstd_fields,              &
        num_nstd_fields=num_avail_tavg_nstd_fields,      &
        max_nstd_fields=max_avail_tavg_nstd_fields       )
 
   !*** return if attempting to compute every nstep timesteps
   if (tavg_freq_iopt == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: MOC diagnostic is not computed if tavg_freq_iopt == freq_opt_nstep'
#ifdef CCSMCOUPLED
        call shr_sys_flush(stdout)
#else
        call POP_IOUnitsFlush(POP_stdout)
#endif
     endif
     return
   endif

   io_dims_nstd_ccsm(1,tavg_MOC) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_MOC) = moc_z_dim
   io_dims_nstd_ccsm(3,tavg_MOC) = moc_comp_dim
   io_dims_nstd_ccsm(4,tavg_MOC) = transport_reg_dim
   io_dims_nstd_ccsm(5,tavg_MOC) = time_dim
   ndims_nstd_ccsm  (  tavg_MOC) = 5

   if ( ldiag_gm_bolus  .and.  lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                  W_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_WISOP),  &
                  V_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_VISOP),  &
                 W_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_WSUBM),  &
                 V_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_VSUBM))
  elseif ( ldiag_gm_bolus  .and.  .not.lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                  W_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_WISOP),  &
                  V_I = TAVG_BUF_3D(:,:,:,:,tavg_loc_VISOP))
  elseif ( .not.ldiag_gm_bolus  .and.  lsubmeso ) then
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ),  &
                 W_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_WSUBM),  &
                 V_SM = TAVG_BUF_3D(:,:,:,:,tavg_loc_VSUBM))
  else
     call compute_moc ( TAVG_BUF_3D(:,:,:,:,tavg_loc_WVEL ),  &
                        TAVG_BUF_3D(:,:,:,:,tavg_loc_VVEL ))
  endif

   ! stop timer
   call timer_stop(timer_tavg_ccsm_diags_moc)
 
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_moc_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_transport_ccsm
! !INTERFACE:

 subroutine tavg_transport_ccsm

! !DESCRIPTION:
!
!  Compute northward transport of heat/salt diagnostically from 
!  other tavg quantities
!
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::    &
      tavg_id_ADVT,         &
      tavg_id_ADVS ,        &
      tavg_id_VNT,          &
      tavg_id_VNS,          &
      tavg_id_HDIFT,        &
      tavg_id_HDIFS,        &
      tavg_id_ADVT_ISOP,    &
      tavg_id_ADVS_ISOP,    &
      tavg_id_VNT_ISOP,     &
      tavg_id_VNS_ISOP,     &
      tavg_loc_ADVT,        &
      tavg_loc_ADVS ,       &
      tavg_loc_VNT,         &
      tavg_loc_VNS,         &
      tavg_loc_HDIFT,       &
      tavg_loc_HDIFS,       &
      tavg_loc_ADVT_ISOP,   &
      tavg_loc_ADVS_ISOP,   &
      tavg_loc_VNT_ISOP,    &
      tavg_loc_VNS_ISOP,    &
      tavg_id_ADVT_SUBM,    &
      tavg_id_ADVS_SUBM,    &
      tavg_id_VNT_SUBM,     &
      tavg_id_VNS_SUBM,     &
      tavg_loc_ADVT_SUBM,   &
      tavg_loc_ADVS_SUBM,   &
      tavg_loc_VNT_SUBM,    &
      tavg_loc_VNS_SUBM    

   integer (int_kind) ::    &
      indx1,                &! index
      indx2,                &! index
      indx3,                &! index
      indx4,                &! index
      indx5,                &! index
      indx6,                &
      indx7,                &
      n

   logical (log_kind) ::  &
      ldiag_gm_bolus,     &  ! local logical for diag_gm_bolus
      lsubmeso               ! local logical for submesoscale_mixing

   if (.not. (n_heat_trans .or. n_salt_trans)) return

   !*** start timer
   call timer_start(timer_tavg_ccsm_diags_trans)

   ldiag_gm_bolus = .false.
   if ( registry_match('diag_gm_bolus') )  ldiag_gm_bolus = .true.

   lsubmeso = .false.
   if ( registry_match('init_submeso') )   lsubmeso = .true.

   tavg_loc_ADVT      = 0 ; tavg_loc_ADVS      = 0
   tavg_loc_VNT       = 0 ; tavg_loc_VNS       = 0
   tavg_loc_HDIFT     = 0 ; tavg_loc_HDIFS     = 0
   tavg_loc_ADVT_ISOP = 0 ; tavg_loc_ADVS_ISOP = 0
   tavg_loc_VNT_ISOP  = 0 ; tavg_loc_VNS_ISOP  = 0
   tavg_loc_ADVT_SUBM = 0 ; tavg_loc_ADVS_SUBM = 0
   tavg_loc_VNT_SUBM  = 0 ; tavg_loc_VNS_SUBM  = 0
 

  !*** error checking 
   tavg_loc_ADVT   = avail_tavg_fields(tavg_id('ADVT' ))%buf_loc
   tavg_loc_ADVS   = avail_tavg_fields(tavg_id('ADVS' ))%buf_loc
   tavg_loc_VNT    = avail_tavg_fields(tavg_id('VNT'  ))%buf_loc
   tavg_loc_VNS    = avail_tavg_fields(tavg_id('VNS'  ))%buf_loc
   tavg_loc_HDIFT  = avail_tavg_fields(tavg_id('HDIFT'))%buf_loc
   tavg_loc_HDIFS  = avail_tavg_fields(tavg_id('HDIFS'))%buf_loc

   !*** error checking

   if (tavg_loc_ADVT  == 0 .or. tavg_loc_ADVS  == 0 .or.  &
       tavg_loc_VNT   == 0 .or. tavg_loc_VNS   == 0 .or.  &
       tavg_loc_HDIFT == 0 .or. tavg_loc_HDIFS == 0) then
       call document ('tavg_transport_ccsm', &
         'Error in heat/salt transport diags: none of the following should be zero')
       call document ('tavg_transport_ccsm', 'tavg_loc_ADVT',  tavg_loc_ADVT)
       call document ('tavg_transport_ccsm', 'tavg_loc_ADVS',  tavg_loc_ADVS)
       call document ('tavg_transport_ccsm', 'tavg_loc_VNT ',  tavg_loc_VNT )
       call document ('tavg_transport_ccsm', 'tavg_loc_VNS ',  tavg_loc_VNS )
       call document ('tavg_transport_ccsm', 'tavg_loc_HDIFT', tavg_loc_HDIFT)
       call document ('tavg_transport_ccsm', 'tavg_loc_HDIFS', tavg_loc_HDIFS)
       call exit_POP (SigAbort, 'Fatal error')
   endif

   if (ldiag_gm_bolus) then
       !*** note: error checking for tavg_id is done in POP_checks
       tavg_loc_ADVT_ISOP = avail_tavg_fields(tavg_id('ADVT_ISOP'))%buf_loc
       tavg_loc_ADVS_ISOP = avail_tavg_fields(tavg_id('ADVS_ISOP'))%buf_loc
       tavg_loc_VNT_ISOP  = avail_tavg_fields(tavg_id('VNT_ISOP') )%buf_loc
       tavg_loc_VNS_ISOP  = avail_tavg_fields(tavg_id('VNS_ISOP') )%buf_loc

       !*** error checking 
       if (tavg_loc_ADVT_ISOP  == 0 .or. tavg_loc_ADVS_ISOP  == 0 .or.  &
           tavg_loc_VNT_ISOP   == 0 .or. tavg_loc_VNS_ISOP   == 0 ) then
          call document ('tavg_transport_ccsm', &
            'Error in heat/salt transport diags: none of the following should be zero')
          call document ('tavg_transport_ccsm', 'tavg_loc_ADVT_ISOP',  tavg_loc_ADVT)
          call document ('tavg_transport_ccsm', 'tavg_loc_ADVS_ISOP',  tavg_loc_ADVS)
          call document ('tavg_transport_ccsm', 'tavg_loc_VNT_ISOP ',  tavg_loc_VNT )
          call document ('tavg_transport_ccsm', 'tavg_loc_VNS_ISOP ',  tavg_loc_VNS )
            call exit_POP (SigAbort, 'Fatal error')
       endif ! tavg_id testing
   endif ! ldiag_gm_bolus

   if (lsubmeso) then
       tavg_loc_ADVT_SUBM = avail_tavg_fields(tavg_id('ADVT_SUBM'))%buf_loc
       tavg_loc_ADVS_SUBM = avail_tavg_fields(tavg_id('ADVS_SUBM'))%buf_loc
       tavg_loc_VNT_SUBM  = avail_tavg_fields(tavg_id('VNT_SUBM' ))%buf_loc
       tavg_loc_VNS_SUBM  = avail_tavg_fields(tavg_id('VNS_SUBM' ))%buf_loc

       !*** error checking
       if (tavg_loc_ADVT_SUBM  == 0 .or. tavg_loc_ADVS_SUBM  == 0 .or.  &
           tavg_loc_VNT_SUBM   == 0 .or. tavg_loc_VNS_SUBM   == 0 ) then
          call document ('tavg_transport_ccsm', &
            'Error in heat/salt transport diags: none of the following should be zero')
          call document ('tavg_transport_ccsm', 'tavg_loc_ADVT_SUBM',  tavg_loc_ADVT_SUBM)
          call document ('tavg_transport_ccsm', 'tavg_loc_ADVS_SUBM',  tavg_loc_ADVS_SUBM)
          call document ('tavg_transport_ccsm', 'tavg_loc_VNT_SUBM ',  tavg_loc_VNT_SUBM )
          call document ('tavg_transport_ccsm', 'tavg_loc_VNS_SUBM ',  tavg_loc_VNS_SUBM )
            call exit_POP (SigAbort, 'Fatal error')
       endif ! tavg_id testing
   endif ! lsubmeso

  !*** define heat transport diagnostics fields and dimensions
  !*** note that this call fills avail_tavg_nstd_fields
   call define_tavg_field     (                     &
        tavg_N_HEAT, 'N_HEAT', 4,                   &
        long_name='Northward Heat Transport',       &
        units='Pwatt',                              &
        coordinates='lat_aux_grid transport_components transport_regions time',&
        nftype='float',                             &
        nstd_fields=avail_tavg_nstd_fields,         &
        num_nstd_fields=num_avail_tavg_nstd_fields, &
        max_nstd_fields=max_avail_tavg_nstd_fields  )

   io_dims_nstd_ccsm(1,tavg_N_HEAT) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_N_HEAT) = transport_comp_dim
   io_dims_nstd_ccsm(3,tavg_N_HEAT) = transport_reg_dim
   io_dims_nstd_ccsm(4,tavg_N_HEAT) = time_dim
   ndims_nstd_ccsm  (  tavg_N_HEAT) = 4
 
  !*** define salt transport diagnostics fields and dimensions
  !*** note that this call fills avail_tavg_nstd_fields
   call define_tavg_field     (                     &
        tavg_N_SALT, 'N_SALT', 4,                   &
        long_name='Northward Salt Transport',       &
        units='gram centimeter^3/kg/s',             &
        coordinates='lat_aux_grid transport_components transport_regions time',&
        nftype='float',                             &
        nstd_fields=avail_tavg_nstd_fields,         &
        num_nstd_fields=num_avail_tavg_nstd_fields, &
        max_nstd_fields=max_avail_tavg_nstd_fields  )

   io_dims_nstd_ccsm(1,tavg_N_SALT) = lat_aux_grid_dim
   io_dims_nstd_ccsm(2,tavg_N_SALT) = transport_comp_dim
   io_dims_nstd_ccsm(3,tavg_N_SALT) = transport_reg_dim
   io_dims_nstd_ccsm(4,tavg_N_SALT) = time_dim
   ndims_nstd_ccsm  (  tavg_N_SALT)  = 4

   !*** return if attempting to compute every nstep timesteps
   if (tavg_freq_iopt == freq_opt_nstep) then
     if (nsteps_run <= 1 .and. my_task == master_task) then
        write(stdout,*)  &
        'WARNING: transport diagnostics are not computed if tavg_freq_iopt == freq_opt_nstep'
#ifdef CCSMCOUPLED
        call shr_sys_flush(stdout)
#else
        call POP_IOUnitsFlush(POP_stdout)
#endif
     endif
     return
   endif

   do n=1,2

     select case (n)
       case(1)
         indx1 = tavg_loc_ADVT
         indx2 = tavg_loc_HDIFT
         indx3 = tavg_loc_VNT
         if (ldiag_gm_bolus) then
           indx4 = tavg_loc_ADVT_ISOP
           indx5 = tavg_loc_VNT_ISOP
         endif
         if ( lsubmeso ) then
           indx6 = tavg_loc_ADVT_SUBM
           indx7 = tavg_loc_VNT_SUBM
         endif
       case(2)
         indx1 = tavg_loc_ADVS
         indx2 = tavg_loc_HDIFS
         indx3 = tavg_loc_VNS
         if (ldiag_gm_bolus) then
           indx4 = tavg_loc_ADVS_ISOP
           indx5 = tavg_loc_VNS_ISOP
         endif
         if ( lsubmeso ) then
           indx6 = tavg_loc_ADVS_SUBM
           indx7 = tavg_loc_VNS_SUBM
         endif
     end select
 
     if ( ldiag_gm_bolus  .and.  lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
           ADV_I = TAVG_BUF_2D(:,:,:,  indx4),  &
           FN_I  = TAVG_BUF_3D(:,:,:,:,indx5),  &
          ADV_SM = TAVG_BUF_2D(:,:,:,  indx6),  &
           FN_SM = TAVG_BUF_3D(:,:,:,:,indx7))
     elseif ( ldiag_gm_bolus  .and.  .not.lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
           ADV_I = TAVG_BUF_2D(:,:,:,  indx4),  &
           FN_I  = TAVG_BUF_3D(:,:,:,:,indx5))
     elseif ( .not.ldiag_gm_bolus  .and.  lsubmeso ) then
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3),  &
          ADV_SM = TAVG_BUF_2D(:,:,:,  indx6),  &
           FN_SM = TAVG_BUF_3D(:,:,:,:,indx7))
     else
       call compute_tracer_transports (n,       &
                   TAVG_BUF_2D(:,:,:,  indx1),  &
                   TAVG_BUF_2D(:,:,:,  indx2),  &
                   TAVG_BUF_3D(:,:,:,:,indx3))
     endif
   enddo
 
  !*** stop timer
   call timer_stop(timer_tavg_ccsm_diags_trans)
!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_transport_ccsm


!***********************************************************************
!BOP
! !IROUTINE: tavg_init_local_spatial_avg
! !INTERFACE:

 subroutine tavg_init_local_spatial_avg

! !DESCRIPTION:
!
!  Initialize geographical masks and arrays for local spatial
!     averaging of some time-averaged fields
!
!       region key:
!        n_reg_0D = 1 --> Nino 1+2 region
!                 = 2 --> Nino  3  region
!                 = 3 --> Nino 3.4 region
!                 = 4 --> Nino  4  region
!
!       TLON_MIN_0D, TLON_MAX_0D, TLAT_MIN_0D, and TLAT_MAX_0D are all
!        in degrees. Also, TLON_MIN_0D and TLON_MAX_0D are in degrees
!        east.
!
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------
   type (block) ::        &
      this_block          ! block information for current block

   real (r8), dimension(n_reg_0D) ::  &
     TLON_MIN_0D = (/ 270.0_r8, 210.0_r8,      &
                      190.0_r8, 160.0_r8 /),   &
     TLON_MAX_0D = (/ 280.0_r8, 270.0_r8,      &
                      240.0_r8, 210.0_r8 /),   &
     TLAT_MIN_0D = (/ -10.0_r8,  -5.0_r8,      &
                       -5.0_r8,  -5.0_r8 /),   &
     TLAT_MAX_0D = (/   0.0_r8,   5.0_r8,      &
                        5.0_r8,   5.0_r8 /)     

   integer (int_kind) ::  &
      iblock,            &! local block number
      n_reg               ! loop index


   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
     TLATD, TLOND        ! lat/lon of T points in degrees


   if (.not. ltavg_nino_diags) return
!-----------------------------------------------------------------------
!
!     allocate and initialize arrays
!
!-----------------------------------------------------------------------

   allocate ( SAVG_0D(n_reg_0D),       &
              SAVG_0D_AREA(n_reg_0D),  &
              SAVG_0D_NAME(n_reg_0D),  &
              SAVG_0D_MASK(nx_block,ny_block,nblocks_clinic,n_reg_0D) )

   SAVG_0D_NAME = (/'NINO_1_PLUS_2 ',  &
                    'NINO_3        ',  &
                    'NINO_3_POINT_4',  &
                    'NINO_4        '/)

   SAVG_0D      = c0
   SAVG_0D_AREA = c0
   SAVG_0D_MASK = c0

   !*** determine masks for each region
   TLATD = TLAT*radian
   TLOND = TLON*radian

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock) 
      do n_reg=1,n_reg_0D
        where (   CALCT(:,:,iblock)              .and.   &
                ( TLOND(:,:,iblock) >= TLON_MIN_0D(n_reg) )  .and.   &
                ( TLOND(:,:,iblock) <= TLON_MAX_0D(n_reg) )  .and.   &
                ( TLATD(:,:,iblock) >= TLAT_MIN_0D(n_reg) )  .and.   &
                ( TLATD(:,:,iblock) <= TLAT_MAX_0D(n_reg) ) )
          SAVG_0D_MASK(:,:,iblock,n_reg) = c1
        endwhere
      enddo ! n_reg
   enddo ! iblock

   !*** compute areas for each region to be used later for normalization
   do n_reg=1,n_reg_0D
     SAVG_0D_AREA(n_reg) = global_sum(TAREA(:,:,:),distrb_clinic,  &
                              field_loc_center,SAVG_0D_MASK(:,:,:,n_reg) )
     if ( SAVG_0D_AREA(n_reg) == c0 ) then
       call exit_POP(sigAbort,  &
          'ERROR in tavg_init_local_spatial_avg: SAVG_0D_AREA is zero.')
     endif
   enddo

!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_init_local_spatial_avg

!***********************************************************************
!BOP
! !IROUTINE: tavg_local_spatial_avg
! !INTERFACE:

 subroutine tavg_local_spatial_avg

! !DESCRIPTION:
!
!  Compute local spatial averages from time-averaged TEMP
!
!
! !REVISION HISTORY:
!  same as module


!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   type (block) ::        &
      this_block          ! block information for current block

   character (char_len) ::  &
      time_string,          &
      date_string

   integer (int_kind) ::  &
      iblock,             &! local block number
      n_reg,              &! loop index
      tavg_id_TEMP,       &! index for TEMP
      tavg_loc_TEMP        ! location for TEMP in TAVG_BUF_3D

   real (r8), dimension (:,:,:), allocatable ::  &
      WORK


   if (.not. ltavg_nino_diags) return

   allocate (WORK(nx_block,ny_block,nblocks_clinic))
   WORK = c0

   tavg_id_TEMP  = tavg_id('TEMP')

   if (tavg_id_TEMP == 0) then
     !*** do not abort; write warning message and proceed
     call document ('tavg_local_spatial_avg', &
       'WARNING: cannot compute spatial averages from time-averaged TEMP')
     return
   else
     tavg_loc_TEMP = avail_tavg_fields(tavg_id_TEMP)%buf_loc
   endif

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock) 
      WORK(:,:,iblock)=  &
          TAVG_BUF_3D(:,:,1,iblock,tavg_loc_TEMP)*TAREA(:,:,iblock)
   enddo ! iblock

   do n_reg=1,n_reg_0D
     SAVG_0D(n_reg) = global_sum ( WORK(:,:,:),distrb_clinic,  &
                         field_loc_center,SAVG_0D_MASK(:,:,:,n_reg))  &
                         / SAVG_0D_AREA(n_reg)
   enddo ! n_reg


   if ( my_task == master_task ) then

     call time_stamp('now','ymd',time_string=time_string)
     call time_stamp('now','ymd',date_string=date_string)

     write (stdout,*) ' '
     write (stdout,*)  &
     ' Local Time- and Space-Averages for Nino Regions: '  &
      // trim(time_string) // ' ' // trim(date_string)

     do n_reg=1,n_reg_0D
       write (stdout,*) trim(SAVG_0D_NAME(n_reg)),': ',  &
                        SAVG_0D(n_reg)
     enddo

   endif
 
   deallocate (WORK)

!-----------------------------------------------------------------------
!EOC

  end subroutine tavg_local_spatial_avg
 
 end module tavg

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
