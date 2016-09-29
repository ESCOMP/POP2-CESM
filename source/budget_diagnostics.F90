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
   use registry
   use tavg
   use time_management
   use forcing_sfwf

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: init_budget_diagnostics,        &
             diag_for_tracer_budgets_rf_Sterms, &
             diag_for_tracer_budgets,        &
             tracer_budgets

   logical (kind=log_kind),public ::  &
      ldiag_global_tracer_budgets      ! global budget diagnostics for tracers

! !PRIVATE MEMBER FUNCTIONS:
   logical (kind=log_kind) ::  &
      lrf_diag_apply_adj_now       ! controls tracer-adjustment term for Robert Filter

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local module variables
!
!-----------------------------------------------------------------------

   character (char_len)    :: exit_string

   integer (int_kind) ::              &
      tavg_id_SHF,                    &
      tavg_id_SFWF,                   &
      tavg_id_RESID_T,                &
      tavg_id_RESID_S,                &
      tavg_id_FW,                     &
      tavg_id_TFW_T,                  &
      tavg_id_TFW_S,                  &
      tavg_id_QFLUX

   integer (int_kind), public ::  &
      budget_stream        ! number of the stream containing budget TAVG fields

   logical (kind=log_kind)     :: budget_warning_1st_step   ! warning flag for budget diagnostics
   logical (log_kind)          :: lrf_print_sequencing
   logical (log_kind)          :: lrf_print_budget_term_by_term

   real (r8), dimension(:,:,:),     allocatable :: WORK1    !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK2    !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK3    !work array
   real (r8), dimension(:,:,:),     allocatable :: WORK4    !work array

   real (r8), dimension(:), allocatable :: budget_rf_S1const       !  S_1 constant for use in this budget
   real (r8), dimension(:), allocatable :: budget_rf_Snconst       !  S_n constant 

   real (r8), dimension(:), allocatable :: rf_adj ! rf budget adjustment terms

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

 
   namelist /budget_diagnostics_nml/ldiag_global_tracer_budgets,  &
             lrf_print_budget_term_by_term


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
!     enable/disable sequencing diagnostic print statements
!
!-----------------------------------------------------------------------
      
   lrf_print_sequencing          = .false.
   lrf_print_budget_term_by_term = .false.

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
!     if not computing budgets, there is nothing more to do here
!
!-----------------------------------------------------------------------

   if (.not. ldiag_global_tracer_budgets) return
   
   allocate (WORK1(nx_block,ny_block,nblocks_clinic))
             WORK1 = c0
   allocate (WORK2(nx_block,ny_block,nblocks_clinic))
             WORK2 = c0
   allocate (WORK3(nx_block,ny_block,nblocks_clinic))
             WORK3 = c0
   allocate (WORK4(nx_block,ny_block,nblocks_clinic))
             WORK4 = c0

!-----------------------------------------------------------------------
!
!     ... except initialize and allocate Robert-Filter related variables
!
!-----------------------------------------------------------------------

   !*** initialize Robert filter budget flags

     if (registry_match('ccsm_startup')) then
       lrf_diag_apply_adj_now = .false.
     else
       lrf_diag_apply_adj_now = .true.
     endif

   !*** variables associated with Robert filtering
   if (lrobert_filter) then

     if (.not. (sfc_layer_type == sfc_layer_varthick)) then
       !*** temporary budget diagnostics limitation; 
       !    add support for this feature later
       exit_string = 'FATAL ERROR: budget diagnostics do not work with this option'
       call exit_POP (sigAbort,trim(exit_string), out_unit=stdout)
     endif

     allocate (rf_adj(nt)) ; rf_adj = c0

     allocate (budget_rf_S1const(nt)) ; budget_rf_S1const      = c0
     allocate (budget_rf_Snconst(nt)) ; budget_rf_Snconst      = c0
   endif


   budget_error_flag = 0

   !*** determine the tavg ids for tavg fields required by this module
     tavg_id_SHF          = tavg_id('SHF')
     tavg_id_SFWF         = tavg_id('SFWF')
     tavg_id_RESID_T      = tavg_id('RESID_T')
     tavg_id_RESID_S      = tavg_id('RESID_S')
     tavg_id_FW           = tavg_id('FW')
     tavg_id_TFW_T        = tavg_id('TFW_T')
     tavg_id_TFW_S        = tavg_id('TFW_S')
     tavg_id_QFLUX        = tavg_id('QFLUX')

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
     call exit_POP (sigAbort, trim(exit_string), out_unit=stdout)
   elseif ( budget_error_flag == -2000 ) then
     exit_string = 'FATAL ERROR: SHF SFWF RESID_T RESID_S '  /&
           &/    'FW TFW_T TFW_S and QFLUX must be included in the tavg_contents file.'
     call document ('init_budget_diagnostics', exit_string)
     call exit_POP (sigAbort, trim(exit_string), out_unit=stdout)
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
     call document ('init_budget_diagnostics', trim(exit_string))
     call exit_POP (sigAbort, trim(exit_string), out_unit=stdout)
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
     call exit_POP (sigAbort, trim(exit_string), out_unit=stdout)
   endif


   call flushm (stdout)

1100  format ('(init_budget_diagnostics): ',a)
1101  format ('(init_budget_diagnostics): tracer budgets are for every ',i4,a)
 
!EOC
   end subroutine init_budget_diagnostics 

!***********************************************************************

!BOP
! !IROUTINE: diag_for_tracer_budgets_rf_Sterms
! !INTERFACE:

   subroutine diag_for_tracer_budgets_rf_Sterms (rf_S_in, rf_volume_total_in, collect_S1, collect_Sn)
                                       
! !DESCRIPTION:
!
!  This subroutine stores information used to compute <S_1> and <S_n> 
!  when the Robert Filter is active.
!  When it is time to compute the budgets, the two terms are differenced
!  and scaled appropriately in subroutine tracer_budgets in order to form
!  the rf_adj term used in tracer_budgets.

! !REVISION HISTORY:
!  created 2016-08  by Nancy J Norton

! !INPUT PARAMETERS:

   real (r8), dimension(:), intent(in) :: rf_S_in
   real (r8),               intent(in) :: rf_volume_total_in

   logical (log_kind), intent(in), optional :: collect_S1, collect_Sn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: k,n
   logical (log_kind) :: lcollect_S1, lcollect_Sn
   real    (r8)       :: sum_tmp

   if (present(collect_S1)) then
    lcollect_S1 = .true.
   else
    lcollect_S1 = .false.
   endif

   if (present(collect_Sn)) then
    lcollect_Sn = .true.
   else
    lcollect_Sn = .false.
   endif

   if (lcollect_S1) then

     if (lrf_print_sequencing)  &
       call document ('diag_for_tracer_budgets_rf_Sterms', 'collect <S1> now')

     !*** define <S_1> for this budget interval and compute <S_1> for next budget interval
     do n=1,2
       !*** do not yet divide by area_t-area_t_marg -- see budget_diagnostics
       budget_rf_S1const(n) = rf_S_in(n)*rf_volume_total_in
     enddo ! n

     return

   else if (lcollect_Sn) then

     if (lrf_print_sequencing)  &
       call document ('diag_for_tracer_budgets_rf_Sterms', 'collect <Sn> now')

     !*** compute <S_n>
     do n=1,2
       !*** do not yet divide by area_t-area_t_marg -- see budget_diagnostics
       budget_rf_Snconst(n) = rf_S_in(n)*rf_volume_total_in
     enddo ! n

   else
     !*** error exit
     call exit_POP (sigAbort, 'either collect_S1 or collect_Sn must be present', out_unit=stdout)
   endif

   end subroutine diag_for_tracer_budgets_rf_Sterms

!***********************************************************************

!BOP
! !IROUTINE: diag_for_tracer_budgets
! !INTERFACE:

   subroutine diag_for_tracer_budgets (tracer_mean, modified_volume_t,  &
                                       MASK_BUDGT, bgtarea_t_k, step_call)

! !DESCRIPTION:
!
!  This subroutine computes the global sum of volume * tracer
!  (returned in tracer_mean) and the global ocean volume at T points
!  at curtime (returned in modified_volume_t).

! !REVISION HISTORY:
!  same as module


! !INPUT PARAMETERS:

   real (kind=r8), dimension(km), intent(in) :: bgtarea_t_k

   real (kind=r8), dimension(nx_block,ny_block,km,nblocks_clinic), intent(in) :: MASK_BUDGT
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
     if (lrobert_filter) then
       if (                           ldiag_global_tracer_budgets .and. accumulate_tavg_now(tavg_TEMP) ) then
         ! continue with diagnostics computation
       else
         return
       endif

       if (lrf_print_sequencing)  &
         call document('diag_for_tracer_budgets', 'collect tracer_mean_initial at nsteps_total', nsteps_total)

     else !not robert filter
       if ( first_global_budget .and. ldiag_global_tracer_budgets .and. accumulate_tavg_now(tavg_TEMP) ) then
         ! continue with diagnostics computation
       else
         return
       endif
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

      if (sfc_layer_type == sfc_layer_varthick ) then
        WORK2(:,:,iblock) = TAREA(:,:,iblock) * (dz(k) + PSURF(:,:,curtime,iblock) / grav )
        WORK4(:,:,iblock) = TAREA(:,:,iblock) * (dz(k) + PSURF(:,:,oldtime,iblock) / grav )
      else
        WORK2(:,:,iblock) = TAREA(:,:,iblock) * dz(k)
        WORK4(:,:,iblock) = TAREA(:,:,iblock) * dz(k)
      endif

      if (lrobert_filter) then
        WORK2(:,:,iblock) =WORK2(:,:,iblock)*MASK_BUDGT(:,:,k,iblock)
        WORK4(:,:,iblock) =WORK4(:,:,iblock)*MASK_BUDGT(:,:,k,iblock)
      else 
        WORK2(:,:,iblock) = merge(WORK2(:,:,iblock), c0, CALCT(:,:,iblock))
        WORK4(:,:,iblock) = merge(WORK4(:,:,iblock), c0, CALCT(:,:,iblock))
      endif

   enddo ! iblock
   !$OMP END PARALLEL DO

   srf_volume = global_sum(WORK2,distrb_clinic, field_loc_center)

   !*** test this with fully coupled version (DEBUG)
   if (lrobert_filter) then
     modified_volume_t = volume_t - volume_t_marg
   else
     modified_volume_t = volume_t
   endif

   if (sfc_layer_type == sfc_layer_varthick ) then
    !*** test this with fully coupled version (DEBUG)
     if (lrobert_filter) then
      modified_volume_t = modified_volume_t - bgtarea_t_k(k)*dz(k) + srf_volume
    else
      modified_volume_t = modified_volume_t -    area_t_k(k)*dz(k) + srf_volume
    endif
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
       sum_tmp    = global_sum(WORK1(:,:,:), distrb_clinic, field_loc_center,MASK_BUDGT(:,:,k,:))
       tracer_mean(n) = sum_tmp
     else
       sum_tmp = global_sum(WORK1(:,:,:)*TAREA(:,:,:), distrb_clinic, field_loc_center,MASK_BUDGT(:,:,k,:))
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

!BOP
! !IROUTINE: tracer_budgets (MASK_BUDGT,bgtarea_t_k))
! !INTERFACE:

   subroutine tracer_budgets (MASK_BUDGT, bgtarea_t_k)

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

! !INPUT PARAMETERS:

   real (kind=r8), dimension(nx_block,ny_block,km,nblocks_clinic), intent(in) :: MASK_BUDGT
   real (kind=r8), dimension(km), intent(in) :: bgtarea_t_k

! !REVISION HISTORY:
!  

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   real (r8) ::                      &
      shf_mean,      sfwf_mean,      &
      resid_t_mean,  resid_s_mean,   &
      qflux_t_mean,  qflux_s_mean,   &
      tfw_t_mean,    tfw_s_mean,     &
      T_change,      S_change,       &
      fw_mean,       volume_change,  &
      tavg_norm,     sum, sumr,      &
      sumrmq0,       sumprnt

   integer (int_kind) ::  k, n

   character*132 :: explanation

!-----------------------------------------------------------------------
!
!   if not yet time to compute and report budgets, return.
!
!-----------------------------------------------------------------------

   if (.not. check_time_flag(tavg_streams(budget_stream)%field_flag)) return


   !*** tavg_sum(budget_stream) = ndays in budget interval * seconds_in_day
   if (lrobert_filter) then
     tavg_norm = c1 / (tavg_sum(budget_stream) * (area_t-area_t_marg))
   else
     tavg_norm = c1 / (tavg_sum(budget_stream) * area_t)
   endif
 
!-----------------------------------------------------------------------
!
!   compute means for SFWF, SHF, RESID_T, RESID_S, FW
!
!-----------------------------------------------------------------------

   if (lrobert_filter) then
     if (sfc_layer_type == sfc_layer_varthick .and.  &
          .not. lfw_as_salt_flx ) then
       sfwf_mean = tavg_global_sum_2D(tavg_id_SFWF,MASK_BUDGT)*rho_sw*c10
     else
       sfwf_mean = tavg_global_sum_2D(tavg_id_SFWF,MASK_BUDGT)
     endif

     shf_mean     =  tavg_global_sum_2D(tavg_id_SHF,MASK_BUDGT)
     resid_t_mean = -tavg_global_sum_2D(tavg_id_RESID_T,MASK_BUDGT)
     resid_s_mean = -tavg_global_sum_2D(tavg_id_RESID_S,MASK_BUDGT)
     fw_mean      =  tavg_global_sum_2D(tavg_id_FW,MASK_BUDGT)
     tfw_t_mean   =  tavg_global_sum_2D(tavg_id_TFW_T,MASK_BUDGT)
     tfw_s_mean   =  tavg_global_sum_2D(tavg_id_TFW_S,MASK_BUDGT)

     !    note that tavg_norm is from this budget interval
   else
     if (sfc_layer_type == sfc_layer_varthick .and.  &
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
   endif

   !*** compute mean heat and virtual salt fluxes due to ice formation

   qflux_t_mean = c0
   if ( tavg_sum_qflux(tavg_in_which_stream(tavg_id_QFLUX)) /= c0 ) then
     if (lrobert_filter) then
       qflux_t_mean = p5*tavg_global_sum_2D(tavg_id('QFLUX'),MASK_BUDGT)
     else
       qflux_t_mean = tavg_global_sum_2D(tavg_id('QFLUX'))
     endif
   endif

   if (sfc_layer_type == sfc_layer_varthick .and.  &
        .not. lfw_as_salt_flx ) then
     qflux_s_mean = c0
   else
     qflux_s_mean = - c10000*qflux_t_mean*             &
               (ocn_ref_salinity - sea_ice_salinity)/  &
                ocn_ref_salinity/latent_heat_fusion /rho_fw
   endif

   !*** start computing time change of volume and tracers

   call diag_for_tracer_budgets (tracer_mean_final, volume_t_final, MASK_BUDGT, bgtarea_t_k)

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

   if (sfc_layer_type == sfc_layer_varthick .and. &
        .not. lfw_as_salt_flx ) then
     S_change = S_change * rho_sw * c10
   else
     S_change = S_change/ salinity_factor
   endif

   if (lrobert_filter) then
      !*** form Robert Filter adjustment term, rf_adj

     if (lrf_conserveVT) then
       do n=1,2
         rf_adj(n) =  c0
       enddo ! n
     else 
       do n=1,2
         rf_adj(n) =  p25*(budget_rf_Snconst(n)-budget_rf_S1const(n))
         !*** now divide by area*T
         rf_adj(n) = rf_adj(n)*robert_curtime*tavg_norm
       enddo ! n

       !*** scale T adjustment term
       rf_adj(1) = rf_adj(1)/hflux_factor

       !*** scale S adjustment term
       !    mimic S_change conversion factor in subroutine tracer_budgets
       if (sfc_layer_type == sfc_layer_varthick .and. .not. lfw_as_salt_flx ) then
         rf_adj(2) = rf_adj(2)*rho_sw*c10
       else
         rf_adj(2) = rf_adj(2)/salinity_factor
       endif
     endif 

   endif ! lrobert_filter

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
 
     if (lrobert_filter .and. lrf_diag_apply_adj_now) then
       sumr    = sum  + rf_adj(1)
       explanation = trim(explanation) /&
                                        &/' + Robert adj'
       write (stdout, 10051) T_change, shf_mean, qflux_t_mean,      &
                             tfw_t_mean,  resid_t_mean, rf_adj(1), &
                             sumr, explanation


     else
       write (stdout, 1005) T_change, shf_mean, qflux_t_mean, &
                            tfw_t_mean,  resid_t_mean, sum,   &
                            explanation
     endif

     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

     if (sfc_layer_type == sfc_layer_varthick .and.  &
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

     if (lrobert_filter .and. lrf_diag_apply_adj_now) then
       sumr = sum + rf_adj(2)

       explanation = trim(explanation) /&
                                        &/' + Robert adj' 
       write (stdout, 10081) S_change, sfwf_mean, qflux_s_mean, &
                             tfw_s_mean,  resid_s_mean, rf_adj(2),&
                             sumr, explanation

     else
       write (stdout, 1008) S_change, sfwf_mean, qflux_s_mean, &
                           tfw_s_mean,  resid_s_mean, sum,    &
                           explanation
     endif

     call POP_IOUnitsFlush(POP_stdout); call POP_IOUnitsFlush(stdout)

   endif

!-----------------------------------------------------------------------
!
!     print Robert-Filter T budget diagnostics term-by-term
!
!-----------------------------------------------------------------------

   if (lrobert_filter .and. lrf_print_budget_term_by_term) then
       write(stdout,*) ' '
       sumprnt =                              - T_change
       call document_dblf ('tracer_budgets', '- T_change                                      ', sumprnt)

       sumprnt =                                shf_mean
       call document_dblf ('tracer_budgets', '  shf_mean                                      ', sumprnt)

       sumprnt =                                qflux_t_mean
       call document_dblf ('tracer_budgets', '  qflux_t_mean                                  ', sumprnt)

       sumprnt =                                rf_adj(1)
       call document_dblf ('tracer_budgets', '  rf_adj(1)                                     ', sumprnt)

       write(stdout,*) ' '
       sumprnt =                              - T_change
       call document_dblf ('tracer_budgets', '- T_change                                      ', sumprnt)

       sumprnt =                              - T_change+ shf_mean
       call document_dblf ('tracer_budgets', '- T_change+ shf_mean                            ', sumprnt)

       sumprnt =                              - T_change+ shf_mean  + qflux_t_mean
       call document_dblf ('tracer_budgets', '- T_change+ shf_mean  + qflux_t_mean            ', sumprnt)

       sumprnt =                              - T_change+ shf_mean  + qflux_t_mean + rf_adj(1)
       call document_dblf ('tracer_budgets', '- T_change+ shf_mean  + qflux_t_mean + rf_adj(1)', sumprnt)

       write(stdout,*) ' '
       sumprnt =                             qflux_t_mean/( - T_change+ shf_mean  + rf_adj(1))
       call document_dblf ('tracer_budgets', 'ratioT: qflux_t_mean/(-T_change+shf_mean+rf_adj)', sumprnt)
       write(stdout,*) ' '

       sumprnt =                              - T_change
       call document_dblf ('tracer_budgets', ' (summary) T signal                             ', sumprnt)

       sumprnt =                              - T_change+ shf_mean  + qflux_t_mean + rf_adj(1)
       call document_dblf ('tracer_budgets', ' (summary) T residual                           ', sumprnt)

!-----------------------------------------------------------------------
!
!     print S budget diagnostics term-by-term
!
!-----------------------------------------------------------------------

       write(stdout,*) ' '
       sumprnt =                              - S_change
       call document_dblf ('tracer_budgets', '- S_change                                      ', sumprnt)

       sumprnt =                               sfwf_mean
       call document_dblf ('tracer_budgets', ' sfwf_mean                                      ', sumprnt)

       sumprnt =                               qflux_s_mean
       call document_dblf ('tracer_budgets', ' qflux_s_mean                                   ', sumprnt)

       sumprnt =                               rf_adj(2)
       call document_dblf ('tracer_budgets', ' rf_adj(2)                                      ', sumprnt)

       write(stdout,*) ' '
       sumprnt =                              - S_change
       call document_dblf ('tracer_budgets', '- S_change                                      ', sumprnt)

       sumprnt =                              - S_change+ sfwf_mean
       call document_dblf ('tracer_budgets', '- S_change+ sfwf_mean                           ', sumprnt)

       sumprnt =                              - S_change+ sfwf_mean + qflux_s_mean
       call document_dblf ('tracer_budgets', '- S_change+ sfwf_mean + qflux_s_mean            ', sumprnt)

       sumprnt =                              - S_change+ sfwf_mean + qflux_s_mean + rf_adj(2)
       call document_dblf ('tracer_budgets', '- S_change+ sfwf_mean + qflux_s_mean + rf_adj(2)', sumprnt)

       write(stdout,*) ' '
       sumprnt=                                      qflux_s_mean/(-S_change+sfwf_mean+rf_adj(2))
       call document_dblf ('tracer_budgets', 'ratioS:qflux_s_mean/(-S_change+sfwf_mean+rf_adj)', sumprnt)
       write(stdout,*) ' '

       sumprnt =                              - S_change
       call document_dblf ('tracer_budgets', ' (summary) S signal                             ', sumprnt)

       sumprnt =                              - S_change+ sfwf_mean + qflux_s_mean + rf_adj(2)
       call document_dblf ('tracer_budgets', ' (summary) S residual                           ', sumprnt)

   endif !lrobert_filter

   if (.not. lrobert_filter) then
     !*** cannot cycle; must recompute *_initial terms
     volume_t_initial    = volume_t_final
     tracer_mean_initial = tracer_mean_final
   endif

!-----------------------------------------------------------------------
!
!     update Robert filter diagnostics flags
!
!-----------------------------------------------------------------------

   if (lrobert_filter) then
     !*** ok to use robert diagnostics in subsequent diagnostics
     lrf_diag_apply_adj_now = .true.
   endif

   if (lrf_print_sequencing)  &
     call document('tracer_budgets', 'End of Budget Interval.  nsteps_total', nsteps_total)

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
