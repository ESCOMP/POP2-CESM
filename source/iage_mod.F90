!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module iage_mod

!BOP
! !MODULE: iage_mod
!
! !DESCRIPTION:
!  
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:


   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size, only: max_blocks_clinic, km, nx_global, ny_global
   use domain, only: nblocks_clinic, distrb_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: km, max_blocks_clinic, tracer_d
   use kinds_mod
!maltrud added hflux_factor, c10
   use constants, only : c0, c1, c10, hflux_factor, char_blank, delim_fmt
   use tracer_types, only : tavg_passive_interior_type, tavg_passive_stf_type
   use msg_mod, only : msg_write
   use io, only : data_set
   use io_types, only: stdout, datafile, io_field_desc, io_dim,   &
       nml_in, nml_filename, construct_file, construct_io_dim,             &
       construct_io_field, rec_type_dbl, destroy_file, destroy_io_field,  &
       luse_pointer_files, pointer_filename, get_unit,      &
       release_unit, destroy_file, add_attrib_file, destroy_io_field,       &
       extract_attrib_file
   use tavg, only: define_tavg_field, tavg_requested, accumulate_tavg_field, &
       tavg_field_desc_ccsm
   use timers, only : get_timer
   use shr_sys_mod

  implicit none
  private


! !PUBLIC MEMBER FUNCTIONS:

  public ::                    &
       iage_tracer_cnt,        &
       iage_ind_begin,         &
       iage_ind_end,           &
       iage_tracer_names,      &
       iage_name2ind,          &
       iage_ind2name,          &
       iage_init,              &
       iage_tavg,              &
       iage_set_interior,      &
       iage_reset

!EOP
!BOC

!maltrud
     character(char_len) :: mmessage

  !--------------------------------------------------------------------
  !   module variables required by forcing_passive_tracer
  !--------------------------------------------------------------------

  integer(int_kind), parameter :: &
       iage_tracer_cnt = 1,       &
       num_iage_tavg_fields = 8

   character (char_len), dimension(iage_tracer_cnt), parameter ::  &
      iage_tracer_names =  (/ 'IAGE' /)
      
!--------------------------------------------------------------
!     index bounds of passive tracer module variables in TRACER
!--------------------------------------------------------------

   integer (int_kind) ::  &
     iage_ind_begin, iage_ind_end

  !--------------------------------------------------------------------
  !   relative tracer indices
  !--------------------------------------------------------------------

  integer(int_kind), parameter :: &
       iage_ind = 1     ! iage index

  !--------------------------------------------------------------------
  !   derived type & parameter for tracer index lookup
  !--------------------------------------------------------------------

  type ind_name_pair
     integer(int_kind) :: ind
     character(char_len) :: name
  end type ind_name_pair

  type(ind_name_pair), dimension(iage_tracer_cnt), parameter :: &
       ind_name_table = (/ ind_name_pair(iage_ind, 'IAGE') /)

  !--------------------------------------------------------------------
  !   derived type for tracer initialization
  !--------------------------------------------------------------------

  type tracer_read
     character(char_len) :: mod_varname, filename, file_varname, file_fmt
     real(r8) :: scale_factor, default_val
  end type tracer_read

  !--------------------------------------------------------------------
  !   define tavg ids for prognostic variables
  !--------------------------------------------------------------------

   integer (int_kind), dimension(num_iage_tavg_fields) :: &
      tavg_IAGE    ! tavg ids for IAGE

   type (io_field_desc), dimension(iage_tracer_cnt) :: &
      TRACER_CUR, TRACER_OLD    ! tracers at current, old times

   logical(log_kind), dimension(:,:,:), allocatable, save :: &
      LAND_MASK

!EOC
!*****************************************************************************

contains

  !*****************************************************************************

  subroutine iage_init(TRACER_MODULE)

    use broadcast, only : broadcast_scalar, broadcast_array
    use constants, only : c0, c2, c1000, p5, field_loc_center, blank_fmt, &
        field_type_scalar
    use prognostic, only : nx_global, ny_global, curtime, oldtime
    use grid, only : KMT, zt, zw, topo_smooth, fill_points
    use forcing_tools, only : find_forcing_times
    use msg_mod, only : msg_set_state, msg_set_iunit
    use time_management, only : freq_opt_nyear, freq_opt_nmonth, init_time_flag

    !------------------------------------------------------------------
    !   arguments
    !------------------------------------------------------------------

    real(r8), dimension(nx_block, ny_block, km, iage_tracer_cnt,   &
         3, max_blocks_clinic), &
         intent(inout) :: TRACER_MODULE

    !------------------------------------------------------------------
    !   local variables
    !------------------------------------------------------------------

    character(*), parameter :: subname = 'iage_mod:iage_init'

    character(char_len) :: &
         init_iage_option,           & ! option for initialization of iage
         init_iage_init_file,        & ! filename for option 'file'
         init_iage_init_file_fmt       ! file format for option 'file'

    logical(log_kind) :: &
         default,             & ! arg to init_time_flag
         lnml_found             ! Was iage_nml found ?

    integer(int_kind) :: &
         n,                   & ! index for looping over tracers
         k,                   & ! index for looping over depth levels
         l,                   & ! index for looping over time levels
         ind,                 & ! tracer index for tracer name from namelist
         iblock,              & ! index for looping over blocks
         nml_error              ! namelist i/o error flag

    type(tracer_read), dimension(iage_tracer_cnt) :: &
         tracer_init_int,      &! namelist variable for initializing tracers
         tracer_init_ext        ! namelist variable for initializing tracers

    namelist /iage_nml/ &
         init_iage_option, init_iage_init_file, tracer_init_ext, &
         init_iage_init_file_fmt

   type (datafile) ::    &
      in_file             ! data file type for init ts file

   type (io_field_desc) :: &
      io_tracer           ! io field descriptors for input Tracer

   type (io_dim) :: &
      i_dim, j_dim, k_dim, month_dim ! dimension descriptors

   real (r8), dimension(:,:,:,:), allocatable :: &
      TEMP_DATA           ! temp array for reading Tracer data

   character (char_len) ::  &
      restart_filename,     &! modified file name for restart file
      restart_pointer_file   ! file name for restart pointer file

   integer (int_kind) :: &
      nu,                &! i/o unit for pointer file reads
      n_absolute,        &! absolute tracer index
      cindx,cindx2        ! indices into character strings

   character (char_len) :: &
      init_ts_option,      &! option for initializing t,s
      init_ts_suboption,   &! suboption for initializing t,s
      init_ts_file,        &! filename for input T,S file
      init_ts_file_fmt      ! format (bin or nc) for input file

   namelist /init_ts_nml/ init_ts_option, init_ts_file, init_ts_file_fmt,  &
                          init_ts_suboption

   type (datafile) :: &
      restart_file    ! io file descriptor

   character (char_len) ::  &
      short_name, long_name  ! tracer name temporaries


    !------------------------------------------------------------------
    !   initialize nf_wrap & msg_mod
    !------------------------------------------------------------------

    call msg_set_state(my_task == master_task)
    call msg_set_iunit(stdout)

    !------------------------------------------------------------------
    !   initialize tracer_d values
    !------------------------------------------------------------------

    do n=1,iage_tracer_cnt
      n_absolute = n + iage_ind_begin - 1
      tracer_d(n_absolute)%short_name = trim(iage_tracer_names(n))
      tracer_d(n_absolute)%long_name  = 'Ideal Age'
      tracer_d(n_absolute)%units      = 'years'
      tracer_d(n_absolute)%tend_units = 'years/s'
      tracer_d(n_absolute)%flux_units = 'cm years/s'
    enddo

    !------------------------------------------------------------------
    !   default namelist settings
    !------------------------------------------------------------------

    init_iage_option = 'unknown'
    init_iage_init_file = 'unknown'
    init_iage_init_file_fmt = 'bin'


    do n = 1,iage_tracer_cnt
       tracer_init_ext(n)%mod_varname  = 'unknown'
       tracer_init_ext(n)%filename     = 'unknown'
       tracer_init_ext(n)%file_varname = 'unknown'
       tracer_init_ext(n)%scale_factor = c1
       tracer_init_ext(n)%default_val  = c0
       tracer_init_ext(n)%file_fmt     = 'bin'
    end do


   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then  
         nml_error = -1
      else
         nml_error =  1      
      endif
      do while (nml_error > 0)
         read(nml_in, nml=iage_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call msg_write(subname, 'iage_nml not found')
       mmessage = 'ERROR : stopping in '/&
                                         &/ subname
       call exit_POP(sigAbort,mmessage)
    end if

    !------------------------------------------------------------------
    !   broadcast all namelist variables
    !------------------------------------------------------------------

    call broadcast_scalar(init_iage_option , master_task)
    call broadcast_scalar(init_iage_init_file, master_task)
    call broadcast_scalar(init_iage_init_file_fmt, master_task)

    do n = 1,iage_tracer_cnt
       call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
       call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
       call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
       call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
       call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
       call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
    end do



    !------------------------------------------------------------------
    !   read in init_ts namelist in case we need it.
    !------------------------------------------------------------------

   init_ts_option  = 'unknown'
   init_ts_suboption  = 'rest'
   init_ts_file    = 'unknown'
   init_ts_file_fmt= 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then  
         nml_error = -1
      else
         nml_error =  1      
      endif
      do while (nml_error > 0)
         read(nml_in, nml=init_ts_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
      if (trim(init_ts_option)    == 'startup' .and.  &
          trim(init_ts_suboption) == 'spunup') then
               init_ts_option = 'startup_spunup'
      endif
   endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call msg_write(subname, 'init_ts_nml not found')
       mmessage = 'ERROR : stopping in '/&
                                         &/ subname
       call exit_POP(sigAbort,mmessage)
    end if

    !------------------------------------------------------------------
    !   broadcast all namelist variables
    !------------------------------------------------------------------

   call broadcast_scalar(init_ts_option    , master_task)
   call broadcast_scalar(init_ts_suboption , master_task)
   call broadcast_scalar(init_ts_file      , master_task)
   call broadcast_scalar(init_ts_file_fmt  , master_task)

    !------------------------------------------------------------------
    !   initialize tracers
    !------------------------------------------------------------------


    select case (init_iage_option)

    case ('startup', 'zero', 'startup_spunup')
       TRACER_MODULE = c0
       if (my_task == master_task) then
           write(stdout,delim_fmt)
           write(stdout,*) ' Initial 3-d Ideal Age set to all zeros' 
           write(stdout,delim_fmt)
           call shr_sys_flush(stdout)
       endif
       
!###################### debug/temporary #######################
!!!! presently, the restart.F90 module is set up
!!!  such that it reads iage (nt = 3) from the standard
!!!  restart file.  This will change when the LANL
!!!  ecosystem model is incorporated into ccsm pop2.

!!!!case ('restart','continue', 'branch', 'hybrid')
    case ('hide_restart','hide_continue', 'hide_branch', 'hide_hybrid')
!################## end debug/temporary #######################
       restart_filename = char_blank
       if (init_iage_init_file == 'same_as_TS') then
          if (init_ts_option /= 'restart' .and. init_ts_option /= 'branch') then
             call exit_POP(sigAbort,  &
                'init_ts_option and init_iage_option are inconsistent')
          end if
          if (luse_pointer_files) then
             restart_filename = char_blank
             restart_pointer_file = char_blank
             call get_unit(nu)
             if (my_task == master_task) then
                restart_pointer_file = pointer_filename
                cindx = len_trim(pointer_filename) + 1
                cindx2= cindx + 7
                restart_pointer_file(cindx:cindx2) = '.restart'
                write(stdout,*) 'Reading pointer file: ', &
                                 trim(restart_pointer_file)
                call shr_sys_flush(stdout)
                open(nu, file=trim(restart_pointer_file), form='formatted', &
                         status='old')
                read(nu,'(a)') restart_filename
                close(nu)
             endif
             call release_unit(nu)

             call broadcast_scalar(restart_filename, master_task)

!--------------------------------------------------------------
!
!  otherwise use input filename
!
!--------------------------------------------------------------

          else
             cindx2 = len_trim(init_ts_file)
             restart_filename(1:cindx2) = trim(init_ts_file)
          endif

          init_iage_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

          cindx2 = len_trim(init_iage_init_file)
          restart_filename(1:cindx2) = trim(init_iage_init_file)

      end if

      restart_file = construct_file(init_iage_init_file_fmt,          &
                        full_name=trim(restart_filename),          &
                        record_length = rec_type_dbl,              &
                        recl_words=nx_global*ny_global)

      call data_set(restart_file, 'open_read')

   i_dim = construct_io_dim('i', nx_global)
   j_dim = construct_io_dim('j', ny_global)
   k_dim = construct_io_dim('k', km)

   do n=1,iage_tracer_cnt
      n_absolute = n + iage_ind_begin - 1
      short_name = char_blank
      short_name = trim(tracer_d(n_absolute)%short_name)/&
                                                &/'_CUR'
      long_name = char_blank
      long_name = trim(tracer_d(n_absolute)%long_name)/&
                                              &/'at current time'

      TRACER_CUR(n) = construct_io_field(trim(short_name),            &
                   dim1=i_dim, dim2=j_dim, dim3=k_dim,                &
                   long_name=trim(long_name),                         &
                   units    =trim(tracer_d(n_absolute)%units),        &
                   grid_loc ='3111',                                  &
                   field_loc = field_loc_center,                      &
                   field_type = field_type_scalar,                    &
                   d3d_array = TRACER_MODULE(:,:,:,n,curtime,:))
      call data_set (restart_file, 'define', TRACER_CUR(n))
   end do

   do n=1,iage_tracer_cnt
      n_absolute = n + iage_ind_begin - 1
      short_name = char_blank
      short_name = trim(tracer_d(n_absolute)%short_name)/&
                                                &/'_OLD'
      long_name = char_blank
      long_name = trim(tracer_d(n_absolute)%long_name)/&
                                              &/'at old time'
      TRACER_OLD(n) = construct_io_field(trim(short_name),            &
                      dim1=i_dim, dim2=j_dim, dim3=k_dim,             &
                      long_name=trim(long_name),                      &
                      units    =trim(tracer_d(n_absolute)%units),     &
                      grid_loc ='3111',                               &
                      field_loc = field_loc_center,                   &
                      field_type = field_type_scalar,                 &
                      d3d_array = TRACER_MODULE(:,:,:,n,oldtime,:))

      call data_set (restart_file, 'define', TRACER_OLD(n))
   end do

!--------------------------------------------------------------
!
!  now we actually read each field
!  after reading, get rid of io field descriptors and close file
!
!--------------------------------------------------------------


   do n=1,iage_tracer_cnt
      call data_set (restart_file, 'read', TRACER_CUR(n))
   end do
   do n=1,iage_tracer_cnt
      call data_set (restart_file, 'read', TRACER_OLD(n))
   end do

   do n=1,iage_tracer_cnt
      call destroy_io_field (TRACER_CUR(n))
   end do
   do n=1,iage_tracer_cnt
      call destroy_io_field (TRACER_OLD(n))
   end do

   call data_set (restart_file, 'close')

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' file read: ', trim(restart_filename)
   endif

    case ('file')
       call msg_write(subname, 'iage being read from separate file')

       !---------------------------------------------------------------
       !   initialize internal tracer_init array
       !---------------------------------------------------------------

       do n = 1,iage_tracer_cnt
          tracer_init_int(n)%mod_varname  = iage_ind2name(n)
          tracer_init_int(n)%filename     = init_iage_init_file
          tracer_init_int(n)%file_varname = iage_ind2name(n)
          tracer_init_int(n)%scale_factor = c1
          tracer_init_int(n)%default_val  = c0
          tracer_init_int(n)%file_fmt     = init_iage_init_file_fmt
       end do

       !---------------------------------------------------------------
       !   copy non-default values from external tracer_init array
       !---------------------------------------------------------------

       do n = 1,iage_tracer_cnt
          if (trim(tracer_init_ext(n)%mod_varname) /= 'unknown') then
             ind = iage_name2ind(tracer_init_ext(n)%mod_varname)

             if (trim(tracer_init_ext(n)%filename) /= 'unknown') &
                  tracer_init_int(ind)%filename = &
                  tracer_init_ext(n)%filename

             if (trim(tracer_init_ext(n)%file_varname) /= 'unknown') &
                  tracer_init_int(ind)%file_varname = &
                  tracer_init_ext(n)%file_varname

             if (tracer_init_ext(n)%scale_factor /= c1) &
                  tracer_init_int(ind)%scale_factor = &
                  tracer_init_ext(n)%scale_factor

             if (tracer_init_ext(n)%default_val /= c1) &
                  tracer_init_int(ind)%default_val = &
                  tracer_init_ext(n)%default_val
          end if
       end do

       !---------------------------------------------------------------
       !   process internal tracer_init array
       !---------------------------------------------------------------

       do n = 1,iage_tracer_cnt
          if (trim(tracer_init_int(n)%filename) == 'none' .or. &
               trim(tracer_init_int(n)%filename) == 'unknown') then
             mmessage = 'initializing ' /&
                                         &/trim(tracer_init_int(n)%mod_varname)/&
                                         &/' to default_val'
             call msg_write(subname,mmessage)
             do iblock = 1,nblocks_clinic
                TRACER_MODULE(:,:,:,n,curtime,iblock) =  &
                   tracer_init_int(n)%default_val
             enddo
          else
             mmessage = 'initializing ' /&
                                         &/trim(tracer_init_int(n)%mod_varname) /&
                                         &/ ' with ' /&
                                         &/trim(tracer_init_int(n)%file_varname) /&
                                         &/ ' from ' /&
                                         &/trim(tracer_init_int(n)%filename)
             call msg_write(subname,mmessage)

             allocate(TEMP_DATA(nx_block,ny_block,km,max_blocks_clinic))

             in_file = construct_file(tracer_init_int(n)%file_fmt,          &
                               full_name=trim(tracer_init_int(n)%filename), &
                               record_length = rec_type_dbl,                &
                               recl_words=nx_global*ny_global)
             call data_set(in_file,'open_read')

             i_dim = construct_io_dim('i',nx_global)
             j_dim = construct_io_dim('j',ny_global)
             k_dim = construct_io_dim('k',km)

             io_tracer = &
                 construct_io_field(trim(tracer_init_int(n)%file_varname), &
                 dim1=i_dim, dim2=j_dim, dim3=k_dim,                       &
                 field_loc = field_loc_center,                             &
                 field_type = field_type_scalar,                           &
                 d3d_array=TEMP_DATA)

             call data_set(in_file,'define',io_tracer)

             call data_set(in_file,'read'  ,io_tracer)
             do iblock=1,nblocks_clinic
                TRACER_MODULE(:,:,:,n,curtime,iblock) = &
                  TEMP_DATA(:,:,:,iblock)*tracer_init_int(n)%scale_factor
                where (TRACER_MODULE(:,:,:,n,curtime,iblock) < c0) &
                  TRACER_MODULE(:,:,:,n,curtime,iblock) = c0
             end do

             call destroy_io_field(io_tracer)

             deallocate(TEMP_DATA)

             call data_set(in_file,'close')
             call destroy_file(in_file)

             if (my_task == master_task) then
                write(stdout,blank_fmt)
                write(stdout,'(a12,a)') ' file read: ', &
                   trim(tracer_init_int(n)%filename)
             endif

          end if
          do iblock=1,nblocks_clinic
             TRACER_MODULE(:,:,:,n,oldtime,iblock) = &
                TRACER_MODULE(:,:,:,n,curtime,iblock)
          enddo
       end do
 
       if (topo_smooth) then
          do k=1,km
             call fill_points(k,TRACER_MODULE(:,:,k,1,curtime,:))
          enddo
       endif


    case default
       call msg_write(subname, 'init_iage_option = ', init_iage_option)
!###################### debug/temporary #######################
!!!    uncomment this when restart/hybrid/branch option has been 
!!!    activated

!      call exit_POP('ERROR: stopping in ' // subname)
!################## end debug/temporary #######################

    end select


    !------------------------------------------------------------------
    !   apply land mask to tracers
    !------------------------------------------------------------------

!!!!do iblock=1,nblocks_clinic
!   do n = 1,iage_tracer_cnt
!      do k = 1,km
!         where (k > KMT(:,:,iblock))
!            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
!            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
!         end where
!      end do
!   end do
!!!!enddo



    !------------------------------------------------------------------
    !   allocate and initialize LAND_MASK (true for ocean points)
    !------------------------------------------------------------------

    allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
    LAND_MASK = merge(.true., .false., KMT > 0)

    call iage_init_tavg

  end subroutine iage_init

  !*****************************************************************************

  subroutine iage_set_interior(k, DTRACER_MODULE)

    !------------------------------------------------------------------
    !   set interior source/sink term for ideal age tracer
    !------------------------------------------------------------------

    use time_management, only : seconds_in_year

    !------------------------------------------------------------------
    !   arguments
    !------------------------------------------------------------------

    integer(int_kind), intent(in) :: &
         k                   ! vertical level index

    real(r8), dimension(nx_block,ny_block,iage_tracer_cnt), intent(out) :: &
         DTRACER_MODULE      ! computed source/sink term

    !------------------------------------------------------------------
    !  below the surface, ideal age increases by one timestep.
    !------------------------------------------------------------------

    if (k > 1) DTRACER_MODULE = c1 / seconds_in_year

   end subroutine iage_set_interior

  !*****************************************************************************

  subroutine iage_reset(TRACER_MODULE)

    !------------------------------------------------------------------
    !   reset surface value for ideal age tracer
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    !   arguments
    !------------------------------------------------------------------

    real(r8), dimension(nx_block,ny_block,km,iage_tracer_cnt), intent(inout) :: &
         TRACER_MODULE      ! ideal age tracer

    !------------------------------------------------------------------
    !  reset ideal age to zero at the surface.
    !------------------------------------------------------------------

    TRACER_MODULE(:,:,1,:) = c0

   end subroutine iage_reset

  !*****************************************************************************

  function iage_name2ind(name)

    !------------------------------------------------------------------
    !   arguments
    !------------------------------------------------------------------

    character(char_len), intent(in) :: name

    !------------------------------------------------------------------
    !   result declaration
    !------------------------------------------------------------------

    integer(int_kind) :: iage_name2ind

    !------------------------------------------------------------------
    !   local variables
    !------------------------------------------------------------------

    character(*), parameter :: subname = 'iage_mod:iage_name2ind'
    integer(int_kind) :: i

    do i = 1,iage_tracer_cnt
       if (trim(name) == trim(ind_name_table(i)%name)) then
          iage_name2ind = ind_name_table(i)%ind
          return
       end if
    end do

    iage_name2ind = 0

  end function iage_name2ind

  !*****************************************************************************

  function iage_ind2name(ind)

    !------------------------------------------------------------------
    !   arguments
    !------------------------------------------------------------------

    integer(int_kind), intent(in) :: ind

    !------------------------------------------------------------------
    !   result declaration
    !------------------------------------------------------------------

    character(char_len) :: iage_ind2name

    !------------------------------------------------------------------
    !   local variables
    !------------------------------------------------------------------

    character(*), parameter :: subname = 'iage_mod:iage_ind2name'
    integer(int_kind) :: i

    do i = 1,iage_tracer_cnt
       if (ind == ind_name_table(i)%ind) then
          iage_ind2name = trim(ind_name_table(i)%name)
          return
       end if
    end do

    call msg_write(subname, 'lookup failed for ', ind)
!   call exit_POP('ERROR : stopping in ' // subname)
!      call exit_POP('ERROR : stopping in WHEREVER')

  end function iage_ind2name

  !*****************************************************************************


  subroutine iage_tavg(bid, k, TRACER_MODULE)

  implicit none

  integer(int_kind), intent(in) :: bid, k
  
  real(r8), dimension(nx_block,ny_block,iage_tracer_cnt), intent(in) :: &
         TRACER_MODULE

  integer(int_kind) :: &
     ii


  do ii=1,num_iage_tavg_fields
    if (tavg_requested(tavg_IAGE(ii))) then
        call accumulate_tavg_field(TRACER_MODULE(:,:,iage_ind), &
                                   tavg_IAGE(ii),bid,k)
     endif
   enddo

   end subroutine iage_tavg

  !*****************************************************************************

   subroutine iage_init_tavg 

   type (tavg_field_desc_ccsm), dimension (num_iage_tavg_fields) :: tavg_iage_field

   integer (int_kind) ::  &
      nn,                 &
      ii 

   character (char_len) :: basename
  


   do nn=1,iage_tracer_cnt
      basename = char_blank
      basename = trim (iage_tracer_names(nn))

        tavg_iage_field( 1)%short_name  = trim(basename)
        tavg_iage_field( 1)%long_name   = 'Ideal Age'
        tavg_iage_field( 1)%units       = 'years'
        tavg_iage_field( 1)%ndims       = 3
        tavg_iage_field( 1)%coordinates = 'TLONG TLAT z_t time'
        tavg_iage_field( 1)%grid_loc    = '3111'
        

        tavg_iage_field( 2)%short_name  = trim(basename)/&
                                                         &/'_SQR'
        tavg_iage_field( 2)%long_name   = 'Ideal Age Squared'
        tavg_iage_field( 2)%units       = 'years^2'
        tavg_iage_field( 2)%ndims       = 3
        tavg_iage_field( 2)%coordinates = 'TLONG TLAT z_t time'
        tavg_iage_field( 2)%grid_loc    = '3111'


        tavg_iage_field( 3)%short_name  = 'J_'/&
                                               &/trim(basename)
        tavg_iage_field( 3)%long_name   = 'unknown J_IAGE'
        tavg_iage_field( 3)%units       = 'unknown'
        tavg_iage_field( 3)%ndims       = 9999
        tavg_iage_field( 3)%coordinates = 'unknown'
        tavg_iage_field( 3)%grid_loc    = 'unknown'


        tavg_iage_field( 4)%short_name  = 'Jint_'/&
                                                  &/trim(basename)
        tavg_iage_field( 4)%long_name   = 'unknown Jint_IAGE'
        tavg_iage_field( 4)%units       = 'unknown'
        tavg_iage_field( 4)%ndims       = 9999
        tavg_iage_field( 4)%coordinates = 'unknown'
        tavg_iage_field( 4)%grid_loc    = 'unknown'


        tavg_iage_field( 5)%short_name  = 'STF_'/&
                                                 &/trim(basename)
        tavg_iage_field( 5)%long_name   = 'unknown STF_IAGE'
        tavg_iage_field( 5)%units       = 'unknown'
        tavg_iage_field( 5)%ndims       = 9999
        tavg_iage_field( 5)%coordinates = 'unknown'
        tavg_iage_field( 5)%grid_loc    = 'unknown'


        tavg_iage_field( 6)%short_name  = 'RESID_'/&
                                                   &/trim(basename)
        tavg_iage_field( 6)%long_name   = 'Free-Surface Residual Flux '/&
                                                   &/'(' /&
                                                          &/trim(basename)/&
                                                                     &/')'
        tavg_iage_field( 6)%units       = 'years/s'
        tavg_iage_field( 6)%ndims       = 2
        tavg_iage_field( 6)%coordinates = 'TLONG TLAT time'
        tavg_iage_field( 6)%grid_loc    = '2110'


        tavg_iage_field( 7)%short_name  = 'FvPER_'/&
                                                   &/trim(basename)
        tavg_iage_field( 7)%long_name   = 'Virtual Flux of '/&
                                                             &/trim(basename)/&
                                                                              &/', P-E+R'
        tavg_iage_field( 7)%units       = 'years/s'
        tavg_iage_field( 7)%ndims       = 2
        tavg_iage_field( 7)%coordinates = 'TLONG TLAT time'
        tavg_iage_field( 7)%grid_loc    = '2110'


        tavg_iage_field( 8)%short_name  = 'FvICE_'/&
                                           &/trim(basename)
        tavg_iage_field( 8)%long_name   = 'Virtual Flux of '/&
                                                             &/trim(basename)/&
                                                                        &/', Ice Formation'
        tavg_iage_field( 8)%units       = 'years/s'
        tavg_iage_field( 8)%ndims       = 2
        tavg_iage_field( 8)%coordinates = 'TLONG TLAT time'
        tavg_iage_field( 8)%grid_loc    = '2110'

!---------------------------------------------------------------------------------
!  define tavg fields computed from passive_tracer routines
!---------------------------------------------------------------------------------
        do ii=1,num_iage_tavg_fields
        call define_tavg_field(                                            &
                          tavg_IAGE(ii),                                   &
                                      trim(tavg_iage_field(ii)%short_name),  &
                                           tavg_iage_field(ii)%ndims,        &
                          long_name=  trim(tavg_iage_field(ii)%long_name),   &
                          units=      trim(tavg_iage_field(ii)%units),       &
                          grid_loc=   trim(tavg_iage_field(ii)%grid_loc),    &
                          coordinates=trim(tavg_iage_field(ii)%coordinates)  )
        enddo ! ii

  enddo ! nn


  end subroutine iage_init_tavg 

  
  !***********************************************************************

end module iage_mod
