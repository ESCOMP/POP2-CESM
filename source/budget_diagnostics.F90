!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module budget_diagnostics

!BOP
! !MODULE: budget_diagnostics
!
! !DESCRIPTION:
!  This module contains routines for tracers budget diagnostics.
!
! !REVISION HISTORY:
!  SVN:$Id$

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod
   use domain_size
   use diagnostics
   use domain
   use constants
   use prognostic
   use io
   use global_reductions
   use grid 
   use exit_mod
   use ice
   use qflux_mod
   use tavg
   use time_management
   use forcing_sfwf

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: init_budget_diagnostics, &
             diag_for_tracer_budgets, &
             diag_for_tracer_budgets_robert1,  &
             diag_for_tracer_budgets_robert2,  &
             tracer_budgets

   logical (kind=log_kind),public ::  &
      ldiag_global_tracer_budgets,    &! global budget diagnostics for tracers
      ldiag_robert_begin_budget,      &! controls collection of tracer-adjustment terms (Robert filter)
      ldiag_robert_budget_term_use_now ! controls usage of tracer-adjustment terms in tracer budget (Robert filter)
 

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::              &
      tavg_id_SHF,                    &
      tavg_id_SFWF,                   &
      tavg_id_RESID_T,                &
      tavg_id_RESID_S,                &
      tavg_id_FW,                     &
      tavg_id_TFW_T,                  &
      tavg_id_TFW_S,                  &
      tavg_id_QFLUX

   integer (int_kind) ::  &
      budget_stream        ! number of the stream containing budget TAVG fields

   logical (kind=log_kind)     ::  &
      budget_warning_1st_step       ! warning flag for budget diagnostics
 
   real (r8), dimension(:,:,:,:,:), allocatable :: TR_STORE
   real (r8), dimension(:,:,:),     allocatable :: PSURF_STORE
   real (r8), dimension(:,:,:),     allocatable :: WORK1 !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK2 !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK3 !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK4 !work array

   real (r8), dimension(:), allocatable :: robert_budget_term ! accumulate RF budget term

!EOC

!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: init_budget_diagnostics
! !INTERFACE:
      subroutine init_budget_diagnostics

! !DESCRIPTION:
!-----------------------------------------------------------------------
!
!  This subroutine checks if the necessary TAVG_2D fields are
!  specified in the tavg_contents file if the volume/tracer
!  budget diagnostics are on.
!
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC


   integer (kind=int_kind) ::   &
      budget_error_flag,        &! error flag for the missing TAVG_2D fields, etc. 
      nml_error,                &! namelist i/o error flag
      tavg_flag                  ! flag to access tavg frequencies

   character (char_len)    :: exit_string
 
   namelist /budget_diagnostics_nml/ldiag_global_tracer_budgets  


!-----------------------------------------------------------------------
!
!     set budget warning flag if this is the first step of a new run
!
!-----------------------------------------------------------------------
      
   if (first_step) then
       budget_warning_1st_step = .true.
   else
       budget_warning_1st_step = .false.
   endif

 
!-----------------------------------------------------------------------
!
!     read budget_diagnostics namelist
!
!-----------------------------------------------------------------------

   ldiag_global_tracer_budgets = .false.
 
   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old', iostat=nml_error)
      if (nml_error /= 0) then
        nml_error = -1
      else
        nml_error =  1
      endif
      !*** keep reading until find right namelist
      do while (nml_error > 0)
        read(nml_in, nml=budget_diagnostics_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   end if


   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      exit_string = 'FATAL ERROR: reading budget_diagnostics_nml namelist'
      call document ('init_budget_diagnostics', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

 
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Budget Diagnostics:'
      write(stdout,blank_fmt)
      write(stdout,*) ' budget_diagnostics_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,budget_diagnostics_nml)
      write(stdout,blank_fmt)
   endif

   call broadcast_scalar(ldiag_global_tracer_budgets, master_task)

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a)') 'Budget Diagnostic Options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
 
      if ( ldiag_global_tracer_budgets ) then
           write (stdout,*) 'Global tracer budgets diagnostics will be computed'
      else
          write (stdout,*) 'Global tracer budgets diagnostics will not be computed'
      endif
   endif

!-----------------------------------------------------------------------
!
!     initialize Robert filter budget flags
!
!      The S1 term cannot ever be computed at beginning of a startup or
!      continuation run, because the info is not available. 
!
!      Budgets from the first time interval will be incorrect.
!
!      The only time these switches are toggled true is when using Robert 
!      filtering (see tracer_budgets); otherwise, they are always false.
!
!-----------------------------------------------------------------------
      
   ldiag_robert_begin_budget        = .false.
   ldiag_robert_budget_term_use_now = .false.

!-----------------------------------------------------------------------
!
!     if not computing budgets, there is nothing more to do here
!
!-----------------------------------------------------------------------

   if (.not. ldiag_global_tracer_budgets) return
 
   allocate (WORK1(nx_block,ny_block,nblocks_clinic))
   allocate (WORK2(nx_block,ny_block,nblocks_clinic))
   allocate (WORK3(nx_block,ny_block,nblocks_clinic))
   allocate (WORK4(nx_block,ny_block,nblocks_clinic))
   WORK1 = c0
   WORK2 = c0
   WORK3 = c0
   WORK4 = c0
   
   !*** variables associated with Robert filtering
   if (lrobert_filter) then
     allocate (robert_budget_term(nt)) 
     robert_budget_term = c0
   
     allocate (TR_STORE   (nx_block,ny_block,km,nt,nblocks_clinic))
     TR_STORE = c0

     allocate (PSURF_STORE(nx_block,ny_block,       nblocks_clinic))
     PSURF_STORE = c0

   endif
 
   budget_error_flag = 0

   !*** determine the tavg ids for tavg fields required by this module
     tavg_id_SHF        = tavg_id('SHF')
     tavg_id_SFWF       = tavg_id('SFWF')
     tavg_id_RESID_T    = tavg_id('RESID_T')
     tavg_id_RESID_S    = tavg_id('RESID_S')
     tavg_id_FW         = tavg_id('FW')
     tavg_id_TFW_T      = tavg_id('TFW_T')
     tavg_id_TFW_S      = tavg_id('TFW_S')
     tavg_id_QFLUX      = tavg_id('QFLUX')

   !*** determine in which stream the budget diagnostics fields reside
   !    (see also the test below)

   budget_stream =  tavg_in_which_stream(tavg_id_SHF)
 
   tavg_flag = tavg_streams(budget_stream)%field_flag
   if (check_time_flag_int(tavg_flag,freq_opt=.true.) == freq_opt_never) then
      budget_error_flag = -1000
   else
     
    !*** determine if required fields are activated in the tavg_contents file

    if (.not. set_in_tavg_contents(tavg_id_SHF))     budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_SFWF))    budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_RESID_T)) budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_RESID_S)) budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_FW))      budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_TFW_T))   budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_TFW_S))   budget_error_flag = -2000
    if (.not. set_in_tavg_contents(tavg_id_QFLUX))   budget_error_flag = -2000

   endif

   if     ( budget_error_flag == -1000 ) then
     exit_string = 'FATAL ERROR: you cannot select both '    /&
          &/    'tavg_freq_opt == freq_opt_never and ldiag_global_tracer_budgets = .true. '
     call document ('init_budget_diagnostics', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   elseif ( budget_error_flag == -2000 ) then
     exit_string = 'FATAL ERROR: SHF SFWF RESID_T RESID_S '  /&
           &/    'FW TFW_T TFW_S and QFLUX must be included in the tavg_contents file.'
     call document ('init_budget_diagnostics', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   !*** determine if all required fields are activated in the *the same* tavg_contents file

   if (tavg_in_this_stream(tavg_id_SFWF   ,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_RESID_T,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_RESID_S,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_FW     ,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_TFW_T  ,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_TFW_S  ,budget_stream)  .and.   &
       tavg_in_this_stream(tavg_id_QFLUX  ,budget_stream) ) then
          ! ok -- all fields are in the same stream
   else
       budget_error_flag = -3000
   endif
   
   if ( budget_error_flag == -3000 ) then
     exit_string = 'FATAL ERROR: you must select all budget diagnostics' /&
          &/    ' fields in the same stream: SFWF,RESID_T,RESID_S,FW,TFW_T,TFW_S,QFLUX'
     call document ('init_budget_diagnostics', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   budget_error_flag = 0

   if ( my_task == master_task ) then
     select case (check_time_flag_int(tavg_flag,freq_opt=.true.))

      case ( freq_opt_never ) 
        write(stdout,1100) 'no budget interval specified'
        budget_error_flag = -1000
 
      case ( freq_opt_nyear ) 
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' years'
 
      case ( freq_opt_nmonth ) 
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' months'
 
      case ( freq_opt_nday ) 
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' days'

      case ( freq_opt_nhour ) 
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' hours'

      case ( freq_opt_nsecond )
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' seconds'
 
      case ( freq_opt_nstep )
        write(stdout,1101) check_time_flag_int(tavg_flag,freq=.true.),' steps'
 
     end select
   endif

   call broadcast_scalar (budget_error_flag, master_task)

   if ( budget_error_flag == -1000 ) then
     exit_string = 'FATAL ERROR: no budget interval is specified.'
     call document ('init_budget_diagnostics', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif


   call flushm (stdout)

1100  format ('(init_budget_diagnostics): ',a)
1101  format ('(init_budget_diagnostics): tracer budgets are for every ',i4,a)
 
!EOC
   end subroutine init_budget_diagnostics 

!***********************************************************************

!BOP
! !IROUTINE: diag_for_tracer_budgets
! !INTERFACE:

   subroutine diag_for_tracer_budgets (tracer_mean, modified_volume_t,  &
                                       step_call)

! !DESCRIPTION:
!
!  This subroutine computes the global sum of volume * tracer
!  (returned in tracer_mean) and the global ocean volume at T points
!  at curtime (returned in modified_volume_t).

! !REVISION HISTORY:
!  same as module
!  modified for Robert Filter term 2015-09-23 by Nancy J Norton


! !INPUT PARAMETERS:

   logical (kind=log_kind), optional, intent(in) :: step_call

! !OUTPUT PARAMETERS:

   real (r8), dimension(nt), intent(out) :: tracer_mean
   real (r8),                intent(out) :: modified_volume_t

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   real (r8) ::  &
      sum_tmp,   &
      srf_volume

   integer (int_kind) ::  &
      iblock,             &
      k,                  &
      n,                  &
      tavg_TEMP

   logical (log_kind), save ::    &
      first_global_budget = .true.  ! flag for initializing budget diagnostics

   tavg_id_SHF        = tavg_id('SHF')

   if (present(step_call)) then
   if (step_call) then
     tavg_TEMP = tavg_id('TEMP')
     if ( first_global_budget .and. ldiag_global_tracer_budgets .and. accumulate_tavg_now(tavg_TEMP) ) then
       ! continue with diagnostics computation
     else
       return
     endif
   endif ! step_call
   endif ! present(step_call)

!-----------------------------------------------------------------------
!
!     compute the ocean volume based on oldtime and curtime
!
!----------------------------------------------------------------------- 

   k = 1

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1,nblocks_clinic

      if ( sfc_layer_type == sfc_layer_varthick ) then
        WORK2(:,:,iblock) = TAREA(:,:,iblock) * (dz(k) + PSURF(:,:,curtime,iblock) / grav )
        WORK4(:,:,iblock) = TAREA(:,:,iblock) * (dz(k) + PSURF(:,:,oldtime,iblock) / grav )
      else
        WORK2(:,:,iblock) = TAREA(:,:,iblock) * dz(k)
        WORK4(:,:,iblock) = TAREA(:,:,iblock) * dz(k)
      endif

      WORK2(:,:,iblock) = merge(WORK2(:,:,iblock), c0, CALCT(:,:,iblock))
      WORK4(:,:,iblock) = merge(WORK4(:,:,iblock), c0, CALCT(:,:,iblock))
   enddo ! iblock
   !$OMP END PARALLEL DO

   srf_volume = global_sum(WORK2,distrb_clinic, field_loc_center)

   modified_volume_t = volume_t

   if ( sfc_layer_type == sfc_layer_varthick ) then
     modified_volume_t = volume_t - area_t_k(k) * dz(k) &
                       + srf_volume
   endif

!-----------------------------------------------------------------------
!
!     compute the sum of volume * tracers 
!
!-----------------------------------------------------------------------

   do n=1,nt

     do k=1,km

       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1,nblocks_clinic

          if ( tmix_iopt == tmix_matsuno ) then
            if ( k == 1 ) then
              WORK1(:,:,iblock) = TRACER(:,:,k,n,curtime,iblock) * WORK2(:,:,iblock)
            else
              WORK1(:,:,iblock) = TRACER(:,:,k,n,curtime,iblock)
            endif
          else
            if ( k == 1 ) then
              WORK1(:,:,iblock) = p5 * (TRACER(:,:,k,n,oldtime,iblock) * WORK4(:,:,iblock)  &
                          + TRACER(:,:,k,n,curtime,iblock) * WORK2(:,:,iblock))
            else
              WORK1(:,:,iblock) = p5 * (TRACER(:,:,k,n,oldtime,iblock)  &
                       + TRACER(:,:,k,n,curtime,iblock))
            endif
          endif
     enddo ! iblock
     !$OMP END PARALLEL DO

     if ( k == 1 ) then
       sum_tmp    = global_sum(WORK1, distrb_clinic, field_loc_center)
       tracer_mean(n) = sum_tmp
     else
       WORK3(:,:,:) = merge (WORK1(:,:,:)*TAREA(:,:,:), c0, KMT(:,:,:) >= k)
       sum_tmp = global_sum(WORK3, distrb_clinic, field_loc_center)
       tracer_mean(n) = tracer_mean(n) + sum_tmp * dz(k)
     endif

   enddo ! k

   enddo ! n

   if (present(step_call)) then
   if (step_call) then
     if ( first_global_budget .and. ldiag_global_tracer_budgets .and. accumulate_tavg_now(tavg_TEMP) ) then
        if ( my_task == master_task ) then
          write (stdout,1001) volume_t_initial,tracer_mean_initial(1),   &
                              salt_to_ppt*tracer_mean_initial(2)
        endif
       first_global_budget = .false.
     endif ! first_global_budget 
   endif ! step_call
   endif ! present(step_call)

 1001 format (/, 10x, 'VOLUME AND TRACER BUDGET INITIALIZATION:',  &
              /, 10x, '========================================',  &
              /,  5x, ' volume_t (cm^3)           = ', e18.12,     &
              /,  5x, ' SUM [volume*T] (C   cm^3) = ', e18.12,     &
              /,  5x, ' SUM [volume*S] (ppt cm^3) = ', e18.12 )

   end subroutine diag_for_tracer_budgets

!***********************************************************************

!BOP
! !IROUTINE: diag_for_tracer_budgets_robert1
! !INTERFACE:

   subroutine diag_for_tracer_budgets_robert1 
                                       
! !DESCRIPTION:
!
!  This subroutine stores the cur_time tracer value, before
!  the Robert filter is applied. The stored values are used
!  in subroutine diag_for_tracer_budgets_robert2, when they
!  are used to compute the Robert filter budget term,
!  when the Robert Filter is active.

! !REVISION HISTORY:
!  created 2015-10  by Nancy J Norton
!  NOTE: UNDER DEVELOPMENT -- incomplete

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock, k, n


!-----------------------------------------------------------------------
!
!    store the pre-filtered fields at curtime
!
!----------------------------------------------------------------------- 

   do n=1,nt
    do k=1,km
       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1,nblocks_clinic
         TR_STORE(:,:,k,n,iblock) = TRACER(:,:,k,n,curtime,iblock)
       enddo ! iblock
    enddo ! k
   enddo ! n

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic
     PSURF_STORE(:,:,iblock) = PSURF(:,:,curtime,iblock)
   enddo ! iblock

   end subroutine diag_for_tracer_budgets_robert1

!***********************************************************************

!BOP
! !IROUTINE: diag_for_tracer_budgets_robert2 
! !INTERFACE:

   subroutine diag_for_tracer_budgets_robert2
                                       
! !DESCRIPTION:
!
!  This subroutine accumulates the budget term associated with the
!  Robert Filter. subroutine diag_for_tracer_budgets_robert1 must
!  be called prior to the Robert filtering of the tracers.
!    

! !REVISION HISTORY:
!  created 2015-10  by Nancy J Norton
!  NOTE: UNDER DEVELOPMENT -- incomplete

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   real (r8) :: sum_preRF  ! temporary variable
   real (r8) :: sum_postRF ! temporary variable
   real (r8) :: vol_preRF  ! pre-filtering volume
   real (r8) :: vol_postRF ! post-filtering volume

   real (r8), save :: vol_interior ! volume interior ocean 

   integer (int_kind) ::  iblock, k, n

   logical (log_kind), save :: first = .true.


!-----------------------------------------------------------------------
!
!    compute the Robert filter budget adjustment term now.
!      robert_cur*adjustment = TRACER(curtime,post-filter)-TRACER(curtime,pre-filter)
!
!----------------------------------------------------------------------- 

   if (sfc_layer_type == sfc_layer_varthick) then
   else
    !*** temporary budget diagnostics limitation; 
    !    add support for this feature later
    call exit_POP (sigAbort, &
      'FATAL ERROR: budget diagnostics do not work with this option', out_unit=stdout)
   endif

!----------------------------------------------------------------------- 
!
!    compute pre- and post-filtering volumes
!
!----------------------------------------------------------------------- 
    if (first) then

      WORK1 = c0
      vol_interior = c0

      !*** volume vertical interior
      do k=2,km
        !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock = 1,nblocks_clinic
          WORK1(:,:,iblock) = dz(k)*TAREA(:,:,iblock)
        enddo ! iblock
       !$OMP END PARALLEL DO
        vol_interior  = vol_interior + global_sum(WORK1, distrb_clinic, field_loc_center)
     enddo ! k
     first = .false.
   endif

   !*** compute and add pre- and post-filtering surface volumes
   k = 1
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic
     WORK1(:,:,iblock) = (dz(1) + PSURF_STORE(:,:,iblock)/grav)*TAREA(:,:,iblock)
     WORK2(:,:,iblock) = (dz(1)+PSURF(:,:,curtime,iblock)/grav)*TAREA(:,:,iblock)
   enddo ! iblock
   !$OMP END PARALLEL DO

   vol_preRF   = vol_interior + global_sum(WORK1, distrb_clinic, field_loc_center)
   vol_postRF  = vol_interior + global_sum(WORK2, distrb_clinic, field_loc_center)

!----------------------------------------------------------------------- 
!
!    compute pre- and post-filtering volume*tracer
!
!----------------------------------------------------------------------- 

   do n=1,nt
     sum_preRF  = c0
     sum_postRF = c0
     !*** sum over vertical interior; surface is computed separately
     do k=2,km
       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1,nblocks_clinic
         WORK3(:,:,iblock) = dz(k)*TAREA(:,:,iblock)
         WORK1(:,:,iblock) = WORK3(:,:,iblock)*TR_STORE(:,:,k,n,iblock)
         WORK2(:,:,iblock) = WORK3(:,:,iblock)*TRACER  (:,:,k,n,curtime,iblock)
       enddo ! iblock
       !$OMP END PARALLEL DO

       sum_preRF  = sum_preRF +global_sum(WORK1, distrb_clinic, field_loc_center)
       sum_postRF = sum_postRF+global_sum(WORK2, distrb_clinic, field_loc_center)
     enddo ! k

   !*** add pre- and post-filtering surface contribution to volume*TRACER 
   !    assumption: sfc_layer_type == sfc_layer_varthick 
   k = 1
   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic

     !*** surface volume*tracer -- pre-Robert Filter
     WORK1(:,:,iblock) = (dz(1) + PSURF_STORE(:,:,iblock)/grav)*  &
                           TR_STORE(:,:,k,n,iblock)*TAREA(:,:,iblock)
     !*** surface volume*tracer -- post-Robert Filter
     WORK2(:,:,iblock) = (dz(1) + PSURF(:,:,curtime,iblock)/grav)*  &
                         TRACER(:,:,k,n,curtime,iblock)*TAREA(:,:,iblock)
   enddo ! iblock
   !$OMP END PARALLEL DO

   sum_preRF   = sum_preRF  + global_sum(WORK1, distrb_clinic, field_loc_center)
   sum_postRF  = sum_postRF + global_sum(WORK2, distrb_clinic, field_loc_center)
     
   if (vol_preRF .ne. c0 .and. vol_postRF .ne. c0) then
     robert_budget_term(n) = robert_budget_term(n) - sum_preRF/vol_preRF + sum_postRF/vol_postRF
   else
     !*** punt
     robert_budget_term(n) = c0
   endif

!-----------------------------------------------------------------------
!
!  toggle begin_budget switch for next budget cycle; 
!    storage has been completed for this cycle after last tracer
!
!-----------------------------------------------------------------------
     if (n == nt) ldiag_robert_begin_budget = .false.

   enddo ! n

   end subroutine diag_for_tracer_budgets_robert2

!***********************************************************************

!BOP
! !IROUTINE: tracer_budgets
! !INTERFACE:

   subroutine tracer_budgets

! !DESCRIPTION:
!
! This routine computes and prints global volume and tracer budgets for a given
! time-averaging interval. because most of the analysis is based
! on some time-averaged variables, the budgets are correct
! if time-averaging starts from the first time step of an intended
! budget interval. The following variables must be requested in 
! the tavg_contents file:
!
!      SHF, SFWF, RESID_T, RESID_S, FW, TFW_T, TFW_S
!
! the unit conversions are such that:
!
!        volume budget            ---->  cm / s
!        heat budget              ---->  W / m^2
!        virtual salt flux budget ---->  kg of fw   / m^2 / s
!        salt flux budget         ---->  kg of salt / m^2 / s

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   real (r8) ::                     &
      shf_mean,     sfwf_mean,      &
      resid_t_mean, resid_s_mean,   &
      qflux_t_mean, qflux_s_mean,   &
      tfw_t_mean,   tfw_s_mean,     &
      T_change,     S_change,       &
      fw_mean,      volume_change,  &
      robert_adj,                   &
      tavg_norm,    sum

   character*132 :: explanation

!-----------------------------------------------------------------------
!
!   if not yet time to compute and report budgets, return.
!   if it IS time, determine if the S1 budget term needs to be computed
!     next timestep, then continue to copmuting and reporting the budgets.
!
!-----------------------------------------------------------------------

   if (.not. check_time_flag(tavg_streams(budget_stream)%field_flag)) then
    if (lrobert_filter) ldiag_robert_begin_budget = .false.
    !*** not yet time to compute and report budgets
    return
   else
    !*** if using Robert filtering, it will be time to collect S1 for 
    !    the NEXT budget interval during the NEXT pass through step.
    if (lrobert_filter) ldiag_robert_begin_budget = .true.

    !*** now, continue on to computations for this budget interval
   endif

   tavg_norm = c1 / (tavg_sum(budget_stream) * area_t)

!-----------------------------------------------------------------------
!
!   compute means for SFWF, SHF, RESID_T, RESID_S, FW
!
!-----------------------------------------------------------------------

   if ( sfc_layer_type == sfc_layer_varthick .and.  &
        .not. lfw_as_salt_flx ) then
     sfwf_mean = tavg_global_sum_2D(tavg_id_SFWF)*rho_sw*c10
   else
     sfwf_mean = tavg_global_sum_2D(tavg_id_SFWF)
   endif

   shf_mean     =  tavg_global_sum_2D(tavg_id_SHF)

   resid_t_mean = -tavg_global_sum_2D(tavg_id_RESID_T)

   resid_s_mean = -tavg_global_sum_2D(tavg_id_RESID_S)

   fw_mean      =  tavg_global_sum_2D(tavg_id_FW)

   tfw_t_mean   =  tavg_global_sum_2D(tavg_id_TFW_T)

   tfw_s_mean   =  tavg_global_sum_2D(tavg_id_TFW_S)

   !*** compute mean heat and virtual salt fluxes due to ice formation

   qflux_t_mean = c0
   if ( tavg_sum_qflux(tavg_in_which_stream(tavg_id_QFLUX)) /= c0 )  &
        qflux_t_mean = tavg_global_sum_2D(tavg_id('QFLUX'))

   if ( sfc_layer_type == sfc_layer_varthick .and.  &
        .not. lfw_as_salt_flx ) then
     qflux_s_mean = c0
   else
     qflux_s_mean = - c10000*qflux_t_mean*             &
               (ocn_ref_salinity - sea_ice_salinity)/  &
                ocn_ref_salinity/latent_heat_fusion /rho_fw
   endif

   !*** start computing time change of volume and tracers

   call diag_for_tracer_budgets (tracer_mean_final, volume_t_final)

   if ( my_task == master_task ) then
     write (stdout,999) volume_t_final,        &
                        tracer_mean_final(1),  &
                        salt_to_ppt * tracer_mean_final(2)
   endif

   volume_change = (volume_t_final - volume_t_initial)*tavg_norm

   T_change = (tracer_mean_final(1) - tracer_mean_initial(1))  &
             * tavg_norm/hflux_factor

   S_change = (tracer_mean_final(2) - tracer_mean_initial(2))  &
             * tavg_norm

   if ( sfc_layer_type == sfc_layer_varthick .and. &
        .not. lfw_as_salt_flx ) then
     S_change = S_change * rho_sw * c10
   else
     S_change = S_change/ salinity_factor
   endif

!-----------------------------------------------------------------------
!
!     print out the budget diagnostics
!
!-----------------------------------------------------------------------

   if ( my_task == master_task ) then

     write (stdout,1000) 
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
 
     if (budget_warning_1st_step) then
         write (stdout,10001)
         call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
         budget_warning_1st_step = .false.
     endif
 
     if (tmix_iopt == tmix_matsuno) write (stdout,1001 )

     write (stdout,1002)
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

     sum = - volume_change + fw_mean
     explanation = ' Imbalance = -tendency + FW flux '
     write (stdout,1003) volume_change, fw_mean, sum, explanation
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

     write (stdout,1004)
     sum = - T_change+ shf_mean + qflux_t_mean + tfw_t_mean &
           + resid_t_mean
 
     explanation = ' Imbalance           = -tendency + SHF flux'  /&
       &/' + ice formation + FW heat content + free-srf non-consv. '
 
     if (lrobert_filter .and. ldiag_robert_budget_term_use_now) then
       robert_adj = p5*robert_budget_term(1)/hflux_factor
       sum = sum + robert_adj
       explanation = trim(explanation) /&
                                       &/' + Robert adj'
       write (stdout, 10051) T_change, shf_mean, qflux_t_mean,      &
                             tfw_t_mean,  resid_t_mean, robert_adj, &
                             sum, explanation
     else
     write (stdout, 1005) T_change, shf_mean, qflux_t_mean, &
                          tfw_t_mean,  resid_t_mean, sum,   &
                          explanation
     endif
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

     if ( sfc_layer_type == sfc_layer_varthick .and.  &
          .not. lfw_as_salt_flx ) then
       write (stdout,1006)
       call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
     else
       write (stdout,1007)
       call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)
     endif
     sum = - S_change + sfwf_mean + qflux_s_mean + tfw_s_mean  &
           + resid_s_mean
     explanation = ' Imbalance = -tendency + SFWF flux'   /&
     &/          ' + ice formation + FW salt content'     /&
     &/          ' + free-srf non-consv. '

     if (lrobert_filter .and. ldiag_robert_budget_term_use_now) then
       robert_adj = p5*robert_budget_term(2)*rho_sw*c10
       sum = sum + robert_adj
       explanation = trim(explanation) /&
                                       &/' + Robert adj'
       write (stdout, 10081) S_change, sfwf_mean, qflux_s_mean, &
                             tfw_s_mean,  resid_s_mean, robert_adj,&
                             sum, explanation
     else
     write (stdout, 1008) S_change, sfwf_mean, qflux_s_mean, &
                          tfw_s_mean,  resid_s_mean, sum,    &
                          explanation
     endif
     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

   endif

   volume_t_initial    = volume_t_final
   tracer_mean_initial = tracer_mean_final

!-----------------------------------------------------------------------
!
!     update Robert filter diagnostics flags and reset robert_budget_term
!
!-----------------------------------------------------------------------

   if (lrobert_filter) then
     !*** ok to use robert diagnostics in subsequent diagnostics
     ldiag_robert_budget_term_use_now = .true.
     robert_budget_term = c0
   endif

!-----------------------------------------------------------------------
!
!     format statements
!
!-----------------------------------------------------------------------

  999 format (/, 10x, 'AFTER CALL IN tracer_budgets (final):',  &
              /, 10x, '=====================================',  &
              /,  5x, ' volume_t (cm^3)           = ', e18.12,  &
              /,  5x, ' SUM [volume*T] (C   cm^3) = ', e18.12,  &
              /,  5x, ' SUM [volume*S] (ppt cm^3) = ', e18.12 )

 1000 format ( /, 10x, 'VOLUME AND TRACER BUDGETS:',                     &
               /, 10x, '==========================',                     &
           //,  5x, 'Please note the following about the budgets:     ', &
                             //, 5x, '*', 3x,                            &
                    'These budgets are meaningful when time-averaging ', &
            /,  9x, 'starts at the first time step of an intended     ', & 
            /,  9x, 'budget interval.                                 ', &
                             //, 5x, '*', 3x,                            &
                    'Positive-signed quantities ==> ocean gain        ', &     
            /,  9x, 'Negative-signed quantities ==> ocean loss        ', &     
                             //, 5x, '*', 3x,                            &
                    'The imbalance should be small in comparison      ', &
            /,  9x, 'to the tendency term, except for the volume      ', &
            /,  9x, 'budget when a fixed-volume ocean is used         ', &
            /,  9x, '(e.g., when lfw_as_salt_flx = .false.)           ', &
                             //, 5x, '*', 3x,                            &
          'These budgets are computed from single-precision TAVG terms,',&
            /,  9x, 'which contributes to inaccuracy in the imbalance ', &
            /,  9x, 'term.')
 
10001 format(/, 5x, 'WARNING: Because the very first time step of a   ', &
            /,  4x, '  new integration is an Euler-forward time step, ', &
            /,  4x, '  any budget including this first step will be',    &
            /,  4x, '  slightly off balance.' )

 1001 format(  /,  5x,'WARNING: These budget diagnostics were designed', &
            /,  4x, '  for use with the "avg" or "avgfit" mixing      ', &
            /,  4x, '  options. You have selected the Matsuno mixing  ', &
            /,  4x, '  option, for which these diagnostics are only   ', &
            /,  4x, '  approximations.                                ')
 
 1002 format ( /, 5x, 'VOLUME BUDGET FOR THE GIVEN TAVG INTERVAL',       &
                      ' (cm/s):',                                        &
               /, 5x, '-----------------------------------------',       &
                      '--------' )

 1003 format ( /, 7x, ' tendency  = ', e18.12,                           &
               /, 7x, ' FW flux   = ', e18.12,                           &
               /, 7x, ' Imbalance = ', e18.12,                           &
               /, 7x,   a                    )

 1004 format ( /, 5x, 'HEAT BUDGET FOR THE GIVEN TAVG INTERVAL',         &
                      ' (W/m^2):',                                       &
               /, 5x, '---------------------------------------',         &
                      '---------' )

 1005 format ( /, 7x, ' tendency            = ', e18.12,                 &
               /, 7x, ' SHF flux            = ', e18.12,                 &
               /, 7x, ' ice formation       = ', e18.12,                 &
               /, 7x, ' FW heat content     = ', e18.12,                 &
               /, 7x, ' free-srf non-consv. = ', e18.12,                 &
               /, 7x, ' Imbalance           = ', e18.12,                 &
               /, 7x,   a                              )

10051 format  (/, 7x, ' tendency            = ', e18.12,                 &
               /, 7x, ' SHF flux            = ', e18.12,                 &
               /, 7x, ' ice formation       = ', e18.12,                 &
               /, 7x, ' FW heat content     = ', e18.12,                 &
               /, 7x, ' free-srf non-consv. = ', e18.12,                 &
               /, 7x, ' Robert adj. term    = ', e18.12,                 &
               /, 7x, ' Imbalance           = ', e18.12,                 &
               /, 7x,   a                              )


 1006 format ( /, 5x, 'SALT BUDGET FOR THE GIVEN TAVG INTERVAL',         &
                      ' (kg of salt/m^2/s):',                            &
               /, 5x, '---------------------------------------',         &
                      '--------------------' )
 1007 format ( /, 5x, 'IMPLIED FW BUDGET FOR THE GIVEN TAVG INTERVAL',   &
                      ' (kg of freshwater/m^2/s):',                      &
               /, 5x, '---------------------------------------------',   &
                      '--------' )

 1008 format ( /, 7x, ' tendency            = ', e18.12,                 &
               /, 7x, ' SFWF flux           = ', e18.12,                 &
               /, 7x, ' ice formation       = ', e18.12,                 &
               /, 7x, ' FW salt content     = ', e18.12,                 &
               /, 7x, ' free-srf non-consv. = ', e18.12,                 &
               /, 7x, ' Imbalance           = ', e18.12,                 &
               /, 7x,   a                              )

10081 format ( /, 7x, ' tendency            = ', e18.12,                 &
               /, 7x, ' SFWF flux           = ', e18.12,                 &
               /, 7x, ' ice formation       = ', e18.12,                 &
               /, 7x, ' FW salt content     = ', e18.12,                 &
               /, 7x, ' free-srf non-consv. = ', e18.12,                 &
               /, 7x, ' Robert adj. term    = ', e18.12,                 &
               /, 7x, ' Imbalance           = ', e18.12,                 &
               /, 7x,   a                              )

!EOC
      end subroutine tracer_budgets

!***********************************************************************

      end module budget_diagnostics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
