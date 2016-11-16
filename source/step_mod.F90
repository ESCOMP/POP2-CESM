!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module step_mod

!BOP
! !MODULE: step_mod

! !DESCRIPTION:
!  Contains the routine for stepping the model forward one timestep
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use blocks
   use domain_size
   use domain
   use constants
   use exit_mod
   use prognostic
   use timers
   use grid
   use diagnostics
   use state_mod, only: state
   use time_management
   use baroclinic
   use barotropic
   use surface_hgt
   use tavg
   use forcing_fields
   use forcing_sfwf, only: lfw_as_salt_flx
   use forcing
   use damping, only : ldamp_uv, damping_uv
   use forcing_shf
   use ice
   use passive_tracers
   use registry
   use communicate
   use global_reductions
   use io_tools
   use io_types
   use budget_diagnostics
   use overflows
   use overflow_type
   use qflux_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: step, init_step

!----------------------------------------------------------------------
!
!   module variables -- general
!
!----------------------------------------------------------------------

   integer (POP_i4), private :: &
      timer_step,              &! timer number for step
      timer_baroclinic,        &! timer for baroclinic parts of step
      timer_barotropic,        &! timer for barotropic part  of step
      timer_3dupdate		! timer for the 3D update after baroclinic component


   integer (POP_i4) :: ierr

!----------------------------------------------------------------------
!
!   module variables -- Robert Filter
!
!----------------------------------------------------------------------

   real (POP_r8), allocatable, dimension(:,:,:,:,:), private :: STORE_RF
   real (POP_r8), allocatable, dimension(:,:,:,:),   private :: VOLSUM_RF
   real (POP_r8), allocatable, dimension(:,:),       private :: WORK1
   real (POP_r8), allocatable, dimension(:,:),       private :: WORK2
   real (POP_r8), allocatable, dimension(:,:),       private :: rf_sum
   real (POP_r8), allocatable, dimension(:),         private :: bgtarea_t_k
   real (POP_r8), allocatable, dimension(:),         private :: rf_S        ! see time_management
   real (POP_r8), allocatable, dimension(:),         private :: rf_S_avg    ! rf conservation term
   real (POP_r8), allocatable, dimension(:),         private :: rf_Svol
   real (POP_r8), allocatable, dimension(:),         private :: rf_Svol_avg ! see time_management
   real (POP_r8), allocatable, dimension(:),         private :: rf_trvol_initial ! optional diagnostic
   real (POP_r8), allocatable, dimension(:),         private :: rf_trvol_final   ! optional diagnostic
   real (POP_r8),                                    private :: rf_qflux_budget_const
   real (POP_r8),                                    private :: rf_sump
   real (POP_r8),                                    private :: rf_volume_2_km
   real (POP_r8),                                    private :: rf_volume_total_cur
   real (POP_r8),                                    private :: rf_volume_total_new
   real (POP_r8),                                    private :: rf_volsfc_cur
   real (POP_r8),                                    private :: rf_volsfc_new
   real (POP_r8),                                    private :: rf_volsfc_invar

   logical (POP_Logical), allocatable, dimension(:,:,:,:)    :: LMASK_BUDGT
   real    (POP_r8),      allocatable, dimension(:,:,:,:), private ::  MASK_BUDGT

   logical (POP_Logical) :: lrf_budget_collect_tracer_mean_initial
   logical (POP_Logical) :: lrf_budget_collect_S1
   logical (POP_Logical) :: lrf_prnt_sequencing = .false.
   logical (POP_Logical) :: lrf_prnt_conserve   = .false.  !sanity check for development purposes




!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: step
! !INTERFACE:

 subroutine step(errorCode)

! !DESCRIPTION:
!  This routine advances the simulation on timestep.
!  It controls logic for leapfrog and/or Matsuno timesteps and performs
!  modified Robert filtering or time-averaging if selected.  
!  Prognostic variables are updated for 
!  the next timestep near the end of the routine.
!  On Matsuno steps, the time (n) velocity and tracer arrays
!  UBTROP,VBTROP,UVEL,VVEL,TRACER contain the predicted new 
!  velocities from the 1st pass for use in the 2nd pass.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------
 
   integer (POP_i4) :: &
      errorCode

   integer (POP_i4) :: &
      i,j,k,n,           &! loop indices
      tmptime,           &! temp space for time index swapping
      iblock,            &! block counter
      ipass,             &! pass counter
      num_passes          ! number of passes through time step
                          ! (Matsuno requires two)
   integer (POP_i4) :: nn ! loop index, ovf_id
   integer (POP_i4) :: ovf_id

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      ZX,ZY,             &! vertically integrated forcing terms
      DH,DHU              ! time change of surface height minus
                          ! freshwater flux at T, U points

   real (POP_r8), dimension(nx_block,ny_block) :: &
      PSURF_FILT_OLD,    &! time filtered PSURF at oldtime
      PSURF_FILT_CUR,    &! time filtered PSURF at curtime
      WORK_MIN,WORK_MAX   ! work variables for enforcing tracer bounds during time filtering

   logical (POP_Logical), save ::    &
      first_call = .true.          ! flag for initializing timers

   type (block) ::        &
      this_block          ! block information for current block

!-----------------------------------------------------------------------
!
!  start step timer
!
!-----------------------------------------------------------------------

   call timer_start(timer_step)

   errorCode = POP_Success

   lpre_time_manager = .true.

!-----------------------------------------------------------------------
!
!  Gather data for comparison with hydrographic data
!
!-----------------------------------------------------------------------
!
!  if(newday) call data_stations
!  if(newday .and. (mod(iday-1,3).eq.0) ) call data_slices
!
!-----------------------------------------------------------------------
!
!  Gather data for comparison with current meter data
!  THIS SECTION NOT FUNCTIONAL AT THIS TIME
!
!-----------------------------------------------------------------------
!
!  if(newday) call data_cmeters
!

!-----------------------------------------------------------------------
!
!     initialize the global budget arrays
!
!-----------------------------------------------------------------------

   if (.not. lrobert_filter) &
   call diag_for_tracer_budgets (tracer_mean_initial,volume_t_initial,  &
                                 MASK_BUDGT, bgtarea_t_k, step_call = .true.)

!-----------------------------------------------------------------------
!
!  read fields for surface forcing
!
!-----------------------------------------------------------------------

   call set_surface_forcing

   if (lidentical_columns) then

     !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
     
       STF(:,:,1,iblock) = global_SHF_coef * RCALCT(:,:,iblock) * hflux_factor
       STF(:,:,2,iblock) = c0 ! * RCALCT(:,:,iblock) * salinity_factor
     
       SHF_QSW(:,:,iblock) = c0 ! * RCALCT(:,:,iblock) * hflux_factor

       SMF(:,:,1,iblock) = global_taux * RCALCT(:,:,iblock)* momentum_factor
       SMF(:,:,2,iblock) = c0 ! * RCALCT(:,:,iblock)* momentum_factor

       SMFT(:,:,:,iblock) = SMF(:,:,:,iblock)

     end do
     !$OMP END PARALLEL DO

   end if

!-----------------------------------------------------------------------
!
!  update timestep counter, set corresponding model time, set
!  time-dependent logical switches to determine program flow.
!
!-----------------------------------------------------------------------

   call time_manager(registry_match('lcoupled'), liceform, licecesm2)

   lpre_time_manager = .false.

   call passive_tracers_send_time



!-----------------------------------------------------------------------
!
!  compute and initialize some time-average diagnostics
!
!-----------------------------------------------------------------------

   call tavg_set_flag(update_time=.true.)
   call tavg_forcing
   if (nt > 2) call passive_tracers_tavg_sflux(STF)
   call movie_forcing


!-----------------------------------------------------------------------
!
!  set timesteps and time-centering parameters for leapfrog or
!  matsuno steps.
!
!-----------------------------------------------------------------------

   mix_pass = 0
   if (matsuno_ts) then
      num_passes = 2
   else
      num_passes = 1
   endif


   do ipass = 1,num_passes



      if (matsuno_ts) mix_pass = mix_pass + 1

      if (leapfrogts) then  ! leapfrog (and averaging) timestep
         mixtime = oldtime
         beta  = alpha
         do k = 1,km
            c2dtt(k) = c2*dt(k)
         enddo
         c2dtu = c2*dtu
         c2dtp = c2*dtp    ! barotropic timestep = baroclinic timestep
         c2dtq = c2*dtu    ! turbulence timestep = mean flow timestep
      else
         mixtime = curtime
         beta  = theta
         do k = 1,km
            c2dtt(k) = dt(k)
         enddo
         c2dtu = dtu
         c2dtp = dtp       ! barotropic timestep = baroclinic timestep
         c2dtq = dtu       ! turbulence timestep = mean flow timestep
      endif

!-----------------------------------------------------------------------
!
!     on 1st pass of matsuno, set time (n-1) variables equal to
!     time (n) variables.
!
!-----------------------------------------------------------------------


      if (mix_pass == 1) then

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock = 1,nblocks_clinic
            UBTROP(:,:,oldtime,iblock) = UBTROP(:,:,curtime,iblock)
            VBTROP(:,:,oldtime,iblock) = VBTROP(:,:,curtime,iblock)
            UVEL(:,:,:,oldtime,iblock) = UVEL(:,:,:,curtime,iblock)
            VVEL(:,:,:,oldtime,iblock) = VVEL(:,:,:,curtime,iblock)
            RHO (:,:,:,oldtime,iblock) = RHO (:,:,:,curtime,iblock)
            TRACER(:,:,:,:,oldtime,iblock) = &
            TRACER(:,:,:,:,curtime,iblock)
         end do
         !$OMP END PARALLEL DO

      endif


!-----------------------------------------------------------------------
!
!     initialize diagnostic flags and sums
!
!-----------------------------------------------------------------------

      call diag_init_sums

!-----------------------------------------------------------------------
!
!     calculate change in surface height dh/dt from surface pressure
!
!-----------------------------------------------------------------------

      call dhdt(DH,DHU)

      call ovf_reg_avgs_prd 

!-----------------------------------------------------------------------
!
!     Integrate baroclinic equations explicitly to find tracers and
!     baroclinic velocities at new time.  Update ghost cells for 
!     forcing terms leading into the barotropic solver.
!
!-----------------------------------------------------------------------

      if(profile_barrier) call POP_Barrier
      call timer_start(timer_baroclinic)
      call baroclinic_driver(ZX,ZY,DH,DHU, errorCode)
      if(profile_barrier) call POP_Barrier
      call timer_stop(timer_baroclinic)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error in baroclinic driver')
         return
      endif

!-----------------------------------------------------------------------
!
!     compute overflow transports
!
!-----------------------------------------------------------------------

      if( overflows_on ) then
         call ovf_driver
      endif
      if ( overflows_on .and. overflows_interactive ) then
         call ovf_rhs_brtrpc_momentum(ZX,ZY)
      endif

      call POP_HaloUpdate(ZX, POP_haloClinic, POP_gridHorzLocNECorner, &
                              POP_fieldKindVector, errorCode,          &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for ZX')
         return
      endif

      call POP_HaloUpdate(ZY, POP_haloClinic, POP_gridHorzLocNECorner, &
                              POP_fieldKindVector, errorCode,          &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for ZY')
         return
      endif

!-----------------------------------------------------------------------
!
!     Solve barotropic equations implicitly to find surface pressure
!     and barotropic velocities.
!
!-----------------------------------------------------------------------

      if(profile_barrier) call POP_Barrier

      if (.not.l1Ddyn) then

        call timer_start(timer_barotropic)
        call barotropic_driver(ZX,ZY,errorCode)
        if(profile_barrier) call POP_Barrier
        call timer_stop(timer_barotropic)

        if (errorCode /= POP_Success) then
           call POP_ErrorSet(errorCode, &
              'Step: error in barotropic')
           return
        endif

      end if

!-----------------------------------------------------------------------
!
!     update tracers using surface height at new time
!     also peform adjustment-like physics (convection, ice formation)
!
!-----------------------------------------------------------------------

      call timer_start(timer_baroclinic)
      call baroclinic_correct_adjust
      call timer_stop(timer_baroclinic)

      if ( overflows_on .and. overflows_interactive ) then
         call ovf_UV_solution
      endif

      if(profile_barrier) call POP_Barrier
      call timer_start(timer_3dupdate)

      call POP_HaloUpdate(UBTROP(:,:,newtime,:), &
                                  POP_haloClinic,                 &
                                  POP_gridHorzLocNECorner,        &
                                  POP_fieldKindVector, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for UBTROP')
         return
      endif

      call POP_HaloUpdate(VBTROP(:,:,newtime,:), &
                                  POP_haloClinic,                 &
                                  POP_gridHorzLocNECorner,        &
                                  POP_fieldKindVector, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for VBTROP')
         return
      endif

      call POP_HaloUpdate(UVEL(:,:,:,newtime,:), & 
                                POP_haloClinic,                 &
                                POP_gridHorzLocNECorner,        &
                                POP_fieldKindVector, errorCode, &
                                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for UVEL')
         return
      endif

      call POP_HaloUpdate(VVEL(:,:,:,newtime,:), &
                                POP_haloClinic,                 &
                                POP_gridHorzLocNECorner,        &
                                POP_fieldKindVector, errorCode, &
                                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for VVEL')
         return
      endif

      call POP_HaloUpdate(RHO(:,:,:,newtime,:), &
                               POP_haloClinic,                 &
                               POP_gridHorzLocCenter,          &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for RHO')
         return
      endif

      call POP_HaloUpdate(TRACER(:,:,:,:,newtime,:), POP_haloClinic, &
                                  POP_gridHorzLocCenter,          &
                                  POP_fieldKindScalar, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for TRACER')
         return
      endif

      call POP_HaloUpdate(QICE(:,:,:), &
                               POP_haloClinic,                 &
                               POP_gridHorzLocCenter,          &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for QICE')
         return
      endif

      call POP_HaloUpdate(AQICE(:,:,:), &
                               POP_haloClinic,                 &
                               POP_gridHorzLocCenter,          &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for AQICE')
         return
      endif


      if(profile_barrier) call POP_Barrier
      call timer_stop(timer_3dupdate)

!-----------------------------------------------------------------------
!
!     add barotropic to baroclinic velocities at new time
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,k,i,j)
      do iblock = 1,nblocks_clinic

         if (l1Ddyn) then
           UBTROP(:,:,newtime,iblock) = c0     
           VBTROP(:,:,newtime,iblock) = c0     
         endif

!CDIR NOVECTOR
         do k=1,km
            do j=1,ny_block
            do i=1,nx_block
               if (k <= KMU(i,j,iblock)) then
                  UVEL(i,j,k,newtime,iblock) = &
                  UVEL(i,j,k,newtime,iblock) + UBTROP(i,j,newtime,iblock)
                  VVEL(i,j,k,newtime,iblock) = &
                  VVEL(i,j,k,newtime,iblock) + VBTROP(i,j,newtime,iblock)
               endif
            enddo
            enddo
         enddo

!-----------------------------------------------------------------------
!
!        Apply damping to UVEL and VVEL
!
!-----------------------------------------------------------------------

         if (ldamp_uv) then
           call damping_uv(UVEL(:,:,:,newtime,iblock),                        &
                           VVEL(:,:,:,newtime,iblock))
         end if

!-----------------------------------------------------------------------
!
!        on matsuno mixing steps update variables and cycle for 2nd pass
!        note: first step is forward only.
!
!-----------------------------------------------------------------------

         if (mix_pass == 1) then

            UBTROP(:,:,curtime,iblock) = UBTROP(:,:,newtime,iblock)
            VBTROP(:,:,curtime,iblock) = VBTROP(:,:,newtime,iblock)
            UVEL(:,:,:,curtime,iblock) = UVEL(:,:,:,newtime,iblock)
            VVEL(:,:,:,curtime,iblock) = VVEL(:,:,:,newtime,iblock)
            RHO (:,:,:,curtime,iblock) = RHO (:,:,:,newtime,iblock)
            TRACER(:,:,:,:,curtime,iblock) = &
            TRACER(:,:,:,:,newtime,iblock)

         endif
      enddo ! block loop
      !$OMP END PARALLEL DO

   end do ! ipass: cycle for 2nd pass in matsuno step

!-----------------------------------------------------------------------
!
!  extrapolate next guess for pressure from three known time levels
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic
      PGUESS(:,:,iblock) = c3*(PSURF(:,:,newtime,iblock) -   &
                               PSURF(:,:,curtime,iblock)) +  &
                               PSURF(:,:,oldtime,iblock)
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  compute some global diagnostics 
!  before updating prognostic variables
!
!-----------------------------------------------------------------------

   call diag_global_preupdate(DH,DHU)

!-----------------------------------------------------------------------
!
!  update prognostic variables for next timestep:
!     on normal timesteps
!        (n) -> (n-1)
!        (n+1) -> (n) 
!     on averaging timesteps
!        [(n) + (n-1)]/2 -> (n-1)
!        [(n+1) + (n)]/2 -> (n)
!
!-----------------------------------------------------------------------

   if (avg_ts .or. back_to_back) then     ! averaging step

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n, &
      !$OMP                     PSURF_FILT_OLD,PSURF_FILT_CUR, &
      !$OMP                     WORK_MIN,WORK_MAX,WORK1,WORK2)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)  

         !*** avg 2-d fields

         UBTROP(:,:,oldtime,iblock) = p5*(UBTROP(:,:,oldtime,iblock) + & 
                                          UBTROP(:,:,curtime,iblock))
         VBTROP(:,:,oldtime,iblock) = p5*(VBTROP(:,:,oldtime,iblock) + &
                                          VBTROP(:,:,curtime,iblock))
         UBTROP(:,:,curtime,iblock) = p5*(UBTROP(:,:,curtime,iblock) + &
                                          UBTROP(:,:,newtime,iblock))
         VBTROP(:,:,curtime,iblock) = p5*(VBTROP(:,:,curtime,iblock) + &
                                          VBTROP(:,:,newtime,iblock))
         GRADPX(:,:,oldtime,iblock) = p5*(GRADPX(:,:,oldtime,iblock) + &
                                          GRADPX(:,:,curtime,iblock))
         GRADPY(:,:,oldtime,iblock) = p5*(GRADPY(:,:,oldtime,iblock) + &
                                          GRADPY(:,:,curtime,iblock))
         GRADPX(:,:,curtime,iblock) = p5*(GRADPX(:,:,curtime,iblock) + &
                                          GRADPX(:,:,newtime,iblock))
         GRADPY(:,:,curtime,iblock) = p5*(GRADPY(:,:,curtime,iblock) + &
                                          GRADPY(:,:,newtime,iblock))
         FW_OLD(:,:,iblock) = p5*(FW(:,:,iblock) + FW_OLD(:,:,iblock))

         !*** avg 3-d fields

         UVEL(:,:,:,oldtime,iblock) = p5*(UVEL(:,:,:,oldtime,iblock) + &
                                          UVEL(:,:,:,curtime,iblock))
         VVEL(:,:,:,oldtime,iblock) = p5*(VVEL(:,:,:,oldtime,iblock) + &
                                          VVEL(:,:,:,curtime,iblock))
         UVEL(:,:,:,curtime,iblock) = p5*(UVEL(:,:,:,curtime,iblock) + &
                                          UVEL(:,:,:,newtime,iblock))
         VVEL(:,:,:,curtime,iblock) = p5*(VVEL(:,:,:,curtime,iblock) + &
                                          VVEL(:,:,:,newtime,iblock))

         do n=1,nt

            do k=2,km
               TRACER(:,:,k,n,oldtime,iblock) =                &
                          p5*(TRACER(:,:,k,n,oldtime,iblock) + &
                              TRACER(:,:,k,n,curtime,iblock))
               TRACER(:,:,k,n,curtime,iblock) =                &
                          p5*(TRACER(:,:,k,n,curtime,iblock) + &
                              TRACER(:,:,k,n,newtime,iblock))
            end do
         end do

         if (sfc_layer_type == sfc_layer_varthick) then

            PSURF_FILT_OLD = p5*(PSURF(:,:,oldtime,iblock) + &
                                 PSURF(:,:,curtime,iblock))
            PSURF_FILT_CUR = p5*(PSURF(:,:,curtime,iblock) + &
                                 PSURF(:,:,newtime,iblock))

            do n = 1,nt
               WORK_MIN = min(TRACER(:,:,1,n,oldtime,iblock), TRACER(:,:,1,n,curtime,iblock))
               WORK_MAX = max(TRACER(:,:,1,n,oldtime,iblock), TRACER(:,:,1,n,curtime,iblock))

               TRACER(:,:,1,n,oldtime,iblock) =                   &
                   p5*((dz(1) + PSURF(:,:,oldtime,iblock)/grav)*  &
                       TRACER(:,:,1,n,oldtime,iblock) +           &
                       (dz(1) + PSURF(:,:,curtime,iblock)/grav)*  &
                       TRACER(:,:,1,n,curtime,iblock) ) 
               TRACER(:,:,1,n,oldtime,iblock) =                   &
                   TRACER(:,:,1,n,oldtime,iblock)/(dz(1) + PSURF_FILT_OLD/grav)

               where (TRACER(:,:,1,n,oldtime,iblock) < WORK_MIN) &
                   TRACER(:,:,1,n,oldtime,iblock) = WORK_MIN
               where (TRACER(:,:,1,n,oldtime,iblock) > WORK_MAX) &
                   TRACER(:,:,1,n,oldtime,iblock) = WORK_MAX


               WORK_MIN = min(TRACER(:,:,1,n,curtime,iblock), TRACER(:,:,1,n,newtime,iblock))
               WORK_MAX = max(TRACER(:,:,1,n,curtime,iblock), TRACER(:,:,1,n,newtime,iblock))

               TRACER(:,:,1,n,curtime,iblock) =                   &
                   p5*((dz(1) + PSURF(:,:,curtime,iblock)/grav)*  &
                       TRACER(:,:,1,n,curtime,iblock) +           &
                       (dz(1) + PSURF(:,:,newtime,iblock)/grav)*  &
                       TRACER(:,:,1,n,newtime,iblock) ) 
               TRACER(:,:,1,n,curtime,iblock) =                   &
                   TRACER(:,:,1,n,curtime,iblock)/(dz(1) + PSURF_FILT_CUR/grav)

               where (TRACER(:,:,1,n,curtime,iblock) < WORK_MIN) &
                   TRACER(:,:,1,n,curtime,iblock) = WORK_MIN
               where (TRACER(:,:,1,n,curtime,iblock) > WORK_MAX) &
                   TRACER(:,:,1,n,curtime,iblock) = WORK_MAX
            enddo

            PSURF(:,:,oldtime,iblock) = PSURF_FILT_OLD
            PSURF(:,:,curtime,iblock) = PSURF_FILT_CUR

         else

            do n=1,nt

               TRACER(:,:,1,n,oldtime,iblock) =                &
                          p5*(TRACER(:,:,1,n,oldtime,iblock) + &
                              TRACER(:,:,1,n,curtime,iblock))
               TRACER(:,:,1,n,curtime,iblock) =                &
                          p5*(TRACER(:,:,1,n,curtime,iblock) + &
                              TRACER(:,:,1,n,newtime,iblock))
            end do

            PSURF (:,:,oldtime,iblock) =                           &
                                  p5*(PSURF (:,:,oldtime,iblock) + &
                                      PSURF (:,:,curtime,iblock))
            PSURF (:,:,curtime,iblock) =                           &
                                  p5*(PSURF (:,:,curtime,iblock) + &
                                      PSURF (:,:,newtime,iblock))

         endif

         do k = 1,km  ! recalculate densities from averaged tracers
            call state(k,k,TRACER(:,:,k,1,oldtime,iblock), &
                           TRACER(:,:,k,2,oldtime,iblock), &
                           this_block,                     &
                         RHOOUT=RHO(:,:,k,oldtime,iblock))
            call state(k,k,TRACER(:,:,k,1,curtime,iblock), &
                           TRACER(:,:,k,2,curtime,iblock), &
                           this_block,                     &
                         RHOOUT=RHO(:,:,k,curtime,iblock))
         enddo 

         !*** correct after avg
         PGUESS(:,:,iblock) = p5*(PGUESS(:,:,iblock) + & 
                                   PSURF(:,:,newtime,iblock)) 
      end do ! block loop
      !$OMP END PARALLEL DO

   else if (lrobert_filter) then   ! robert filter

      call step_RF (DH)

   else  ! non-averaging step
  
      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic

         if (mix_pass == 2) then ! reset time n variables on 2nd pass matsuno

            UBTROP(:,:,curtime,iblock) = UBTROP(:,:,oldtime,iblock)
            VBTROP(:,:,curtime,iblock) = VBTROP(:,:,oldtime,iblock)
            UVEL(:,:,:,curtime,iblock) = UVEL(:,:,:,oldtime,iblock)
            VVEL(:,:,:,curtime,iblock) = VVEL(:,:,:,oldtime,iblock)
            TRACER(:,:,:,:,curtime,iblock) = &
                                     TRACER(:,:,:,:,oldtime,iblock)
            RHO(:,:,:,curtime,iblock) = RHO(:,:,:,oldtime,iblock)

         endif

         FW_OLD(:,:,iblock) = FW(:,:,iblock)

      end do ! block loop
      !$OMP END PARALLEL DO


      tmptime = oldtime
      oldtime = curtime
      curtime = newtime
      newtime = tmptime

   endif


!-----------------------------------------------------------------------
!
!  end of timestep, all variables updated
!  compute and print some more diagnostics
!
!-----------------------------------------------------------------------

   if (registry_match('lcoupled')) then
   if ( liceform .and. check_time_flag(ice_cpl_flag) ) then
     call tavg_increment_sum_qflux(const=tlast_ice)
     !$OMP PARALLEL DO PRIVATE(iblock)
     do iblock = 1,nblocks_clinic
        call ice_flx_to_coupler(TRACER(:,:,:,:,curtime,iblock),iblock)
        call accumulate_tavg_field(QFLUX(:,:,iblock), tavg_id('QFLUX'),  &
                                   iblock,1,const=tlast_ice)
     end do ! block loop
     !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!    time-averaging for ice formation related quantities
!-----------------------------------------------------------------------
     if (nt > 2) call passive_tracers_tavg_FvICE(cp_over_lhfusion, QICE)
   endif
   endif

   call diag_global_afterupdate
   call diag_print
   call diag_transport

   if (eod .and. ldiag_velocity) call diag_velocity


   if (ldiag_global_tracer_budgets) then

     if (lrobert_filter)  then

       !***  collect <S_1> for this budget interval
       if (nsteps_total == 1 .or. nsteps_run == 1 .or. lrf_budget_collect_S1) then
         call diag_for_tracer_budgets_rf_Sterms (rf_S, rf_volume_total_cur, collect_S1 = .true.)
         lrf_budget_collect_S1 = .false.
       endif

       if (check_time_flag(tavg_streams(budget_stream)%field_flag)) then
         !*** the next time around the step loop, which is
         !    timestep after the end of the budget interval

         !  1) collect S1 
         lrf_budget_collect_S1 = .true.

         !  2) compute tracer_mean_initial
         lrf_budget_collect_tracer_mean_initial = .true.
       else
         lrf_budget_collect_tracer_mean_initial = .false.
       endif
     
     endif ! lrobert_filter

     call tracer_budgets (MASK_BUDGT, bgtarea_t_k)

   endif !ldiag_global_tracer_budgets

   if (lrf_prnt_sequencing) call step_RF_doc ('END')

!-----------------------------------------------------------------------
!
!  stop step timer
!
!-----------------------------------------------------------------------

  call timer_stop(timer_step)

!-----------------------------------------------------------------------
!EOC

   end subroutine step

!***********************************************************************

!BOP
! !IROUTINE: step_RF (DH)
! !INTERFACE:

 subroutine step_RF (DH)

! !DESCRIPTION:
!  This routine applies Robert Filtering
!
! !REVISION HISTORY:
!  added July 2016 njn01

! !INPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: DH

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   character (char_len) :: exit_string

   integer (POP_i4) :: &
      i,j,k,n,         &! loop indices
      iblock,          &! block counter
      tmptime           ! temp space for time index swapping

   real (POP_r8)    :: rf_conservation_factor
   real (POP_r8)    :: rf_conserve_factor_new
   real (POP_r8)    :: rf_conserve_factor_cur

   type (block)     :: this_block  ! block information for current block


!-----------------------------------------------------------------------
!
!     optional diagnostics
!
!-----------------------------------------------------------------------

   if (lrf_prnt_sequencing) call step_RF_doc ('begin')

   if (lrf_prnt_conserve)   call step_RF_diag ('initial')

!-----------------------------------------------------------------------
!
!     begin Robert Filtering
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n,WORK1,WORK2)
   do iblock = 1,nblocks_clinic

      !*** filter UBTROP
      WORK1 = & 
      UBTROP(:,:,oldtime,iblock) + UBTROP(:,:,newtime,iblock) - c2*UBTROP(:,:,curtime,iblock)

      if (lrf_nonzero_newtime) &
      UBTROP(:,:,newtime,iblock) = UBTROP(:,:,newtime,iblock) + robert_newtime*WORK1
      UBTROP(:,:,curtime,iblock) = UBTROP(:,:,curtime,iblock) + robert_curtime*WORK1

     !*** filter VBTROP
      WORK1 = & 
      VBTROP(:,:,oldtime,iblock) + VBTROP(:,:,newtime,iblock) - c2*VBTROP(:,:,curtime,iblock)

      if (lrf_nonzero_newtime) &
      VBTROP(:,:,newtime,iblock) = VBTROP(:,:,newtime,iblock) + robert_newtime*WORK1
      VBTROP(:,:,curtime,iblock) = VBTROP(:,:,curtime,iblock) + robert_curtime*WORK1

      !*** filter GRADPX
      WORK1 = &
      GRADPX(:,:,oldtime,iblock) + GRADPX(:,:,newtime,iblock) - c2*GRADPX(:,:,curtime,iblock)

      if (lrf_nonzero_newtime) &
      GRADPX(:,:,newtime,iblock) = GRADPX(:,:,newtime,iblock) + robert_newtime*WORK1
      GRADPX(:,:,curtime,iblock) = GRADPX(:,:,curtime,iblock) + robert_curtime*WORK1

      !*** filter GRADPY
      WORK1 = &
      GRADPY(:,:,oldtime,iblock) + GRADPY(:,:,newtime,iblock) - c2*GRADPY(:,:,curtime,iblock)

      if (lrf_nonzero_newtime) &
      GRADPY(:,:,newtime,iblock) = GRADPY(:,:,newtime,iblock) + robert_newtime*WORK1
      GRADPY(:,:,curtime,iblock) = GRADPY(:,:,curtime,iblock) + robert_curtime*WORK1

      do k=1,km
        !*** filter UVEL
        WORK1 = &
        UVEL(:,:,k,oldtime,iblock) + UVEL(:,:,k,newtime,iblock) - c2*UVEL(:,:,k,curtime,iblock)

        if (lrf_nonzero_newtime) &
        UVEL(:,:,k,newtime,iblock) = UVEL(:,:,k,newtime,iblock) + robert_newtime*WORK1
        UVEL(:,:,k,curtime,iblock) = UVEL(:,:,k,curtime,iblock) + robert_curtime*WORK1

        !*** filter VVEL
        WORK1 = &
        VVEL(:,:,k,oldtime,iblock) + VVEL(:,:,k,newtime,iblock) - c2*VVEL(:,:,k,curtime,iblock)

        if (lrf_nonzero_newtime) &
        VVEL(:,:,k,newtime,iblock) = VVEL(:,:,k,newtime,iblock) + robert_newtime*WORK1
        VVEL(:,:,k,curtime,iblock) = VVEL(:,:,k,curtime,iblock) + robert_curtime*WORK1

      enddo ! k

      !*** filter TRACERS (vertical interior)
      do n=1,nt
        do k=2,km
          !*** form and store Robert Filter term, S
          STORE_RF(:,:,k,n,iblock) = &
          TRACER(:,:,k,n,oldtime,iblock) + TRACER(:,:,k,n,newtime,iblock)  &
                                      - c2*TRACER(:,:,k,n,curtime,iblock)

          !*** filter TRACER
          if (lrf_nonzero_newtime) &
          TRACER(:,:,k,n,newtime,iblock) =  &
          TRACER(:,:,k,n,newtime,iblock) + robert_newtime*STORE_RF(:,:,k,n,iblock)

          TRACER(:,:,k,n,curtime,iblock) = TRACER(:,:,k,n,curtime,iblock) +  &
                                           robert_curtime*STORE_RF(:,:,k,n,iblock)

        enddo ! k
      enddo ! n
     end do ! block loop 
     !$OMP END PARALLEL DO

     !*** accumulate (horizontally masked) volume*S over vertical interior
     do n=1,nt
       rf_Svol(n) = c0 
       do k=2,km
         rf_Svol(n) = rf_Svol(n)  &
                     + global_sum(STORE_RF(:,:,k,n,:)*TAREA(:,:,:)*dz(k),  &
                       distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
       enddo ! k
     enddo ! n
      
      !*** surface TRACERs and PSURF
      if (sfc_layer_type == sfc_layer_varthick) then 
         !====================================

        !*** filter TRACER*surface thickness at the surface
        !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)
        do iblock = 1,nblocks_clinic
        do n=1,nt
            k=1
            !*** form RF term and store for later conservation adjustment
            STORE_RF(:,:,k,n,iblock) = &
            (dz(1) + PSURF(:,:,oldtime,iblock)/grav)*TRACER(:,:,k,n,oldtime,iblock) +  &
            (dz(1) + PSURF(:,:,newtime,iblock)/grav)*TRACER(:,:,k,n,newtime,iblock) -  &
         c2*(dz(1) + PSURF(:,:,curtime,iblock)/grav)*TRACER(:,:,k,n,curtime,iblock)

            !*** filter TRACER*surface thickness
            if (lrf_nonzero_newtime) &
            TRACER(:,:,k,n,newtime,iblock) = (dz(1) + PSURF(:,:,newtime,iblock)/grav)*  &
            TRACER(:,:,k,n,newtime,iblock) + robert_newtime*STORE_RF(:,:,k,n,iblock)

            TRACER(:,:,k,n,curtime,iblock) = (dz(1) + PSURF(:,:,curtime,iblock)/grav)*  &
            TRACER(:,:,k,n,curtime,iblock) + robert_curtime*STORE_RF(:,:,k,n,iblock)

        enddo ! n
        end do ! block loop 
        !$OMP END PARALLEL DO

        !*** break iblock loop and form RF TRACER conservation adjustment term at surface;
        !    compute surface rf_Svol (thickness is already included in STORE_RF(k=1))
        do n=1,nt
          k=1
          rf_Svol(n) = rf_Svol(n) +  &
            global_sum(STORE_RF(:,:,k,n,:)*TAREA(:,:,:), &
                       distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
        enddo ! n
        
        !*** filter PSURF
        !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)
        do iblock = 1,nblocks_clinic
          k=1
          n=1
          STORE_RF(:,:,k,n,iblock) =  &
            PSURF(:,:,oldtime,iblock) + PSURF(:,:,newtime,iblock) - c2*PSURF(:,:,curtime,iblock)

          if (lrf_nonzero_newtime) &
          PSURF(:,:,newtime,iblock) =  &
          PSURF(:,:,newtime,iblock) + robert_newtime*STORE_RF(:,:,k,n,iblock)

          PSURF(:,:,curtime,iblock) =  &
          PSURF(:,:,curtime,iblock) + robert_curtime*STORE_RF(:,:,k,n,iblock)

        end do ! block loop (iblock)
        !$OMP END PARALLEL DO

        !*** compute RF conservation adjustment term for PSURF
        k=1
        n=1
        rf_sump = global_sum(STORE_RF(:,:,k,n,:)*TAREA(:,:,:),  &
                             distrb_clinic,field_loc_center, MASK_BUDGT(:,:,k,:))
        rf_sump = rf_sump/bgtarea_t_k(k)

        !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n,WORK2)
        do iblock = 1,nblocks_clinic

        !*** apply RF conservation adjustment term to PSURF
        k=1
        WORK2(:,:) = merge (rf_sump, c0, LMASK_BUDGT(:,:,k,iblock))

        if (lrf_nonzero_newtime) &
        PSURF(:,:,newtime,iblock) = PSURF(:,:,newtime,iblock) - robert_newtime*WORK2(:,:)
        PSURF(:,:,curtime,iblock) = PSURF(:,:,curtime,iblock) - robert_curtime*WORK2(:,:)

        !*** solve for surface TRACER = (TRACER*thickness)/thickness
        do n=1,nt

          if (lrf_nonzero_newtime) &
          TRACER(:,:,1,n,newtime,iblock) =  &
          TRACER(:,:,1,n,newtime,iblock)/(dz(1) + PSURF(:,:,newtime,iblock)/grav)

          TRACER(:,:,1,n,curtime,iblock) =  &
          TRACER(:,:,1,n,curtime,iblock)/(dz(1) + PSURF(:,:,curtime,iblock)/grav)

        enddo ! n
        end do ! block loop (iblock)
        !$OMP END PARALLEL DO


     !*** note: apply RF conservation term to TRACERS after end of if-block

      else ! sfc_layer_type .ne. sfc_layer_varthick 
           ! ======================================

        exit_string = 'FATAL ERROR: must use sfc_layer_type sfc_layer_varthick with Robert Filter option'
        call document ('step_RF', exit_string)
        call exit_POP (sigAbort, exit_string, out_unit=stdout)


      endif ! sfc_layer_type == sfc_layer_varthick
            ! ====================================

      rf_volsfc_new = global_sum(TAREA(:,:,:)*(dz(1) + PSURF(:,:,newtime,:)/grav),   &
                             distrb_clinic,field_loc_center,MASK_BUDGT(:,:,1,:))
      rf_volsfc_cur = global_sum(TAREA(:,:,:)*(dz(1) + PSURF(:,:,curtime,:)/grav),   &
                             distrb_clinic,field_loc_center,MASK_BUDGT(:,:,1,:))

      rf_volume_total_cur = rf_volume_2_km + rf_volsfc_cur

      do n=1,nt

         rf_S(n) = rf_Svol(n)/rf_volume_total_cur

         if (nsteps_total == 1 .or. lrf_nonzero_newtime .or. lrf_conserveVT) then
           rf_conservation_factor = rf_S(n)
          !*** placeholder for lrf_nonzero_newtime; this is not correct
         else
          !*** control numerical instability when robert_alpha = 1
          rf_conservation_factor = p5*(rf_S(n)+rf_S_prev(n))
         endif

         rf_conserve_factor_new = rf_conservation_factor*robert_newtime
         rf_conserve_factor_cur = rf_conservation_factor*robert_curtime

         if (lrf_conserveVT) call step_RF_doc ('conserve', n, rf_conserve_factor_cur)

         !***  apply RF conservation adjustment to TRACERs at all vertical levels
         !$OMP PARALLEL DO PRIVATE(iblock,this_block,k)
         do iblock = 1,nblocks_clinic

          do k=1,km
            if (lrf_nonzero_newtime) &
            TRACER(:,:,k,n,newtime,iblock) =  &
            TRACER(:,:,k,n,newtime,iblock) - rf_conserve_factor_new*MASK_BUDGT(:,:,k,iblock)

            TRACER(:,:,k,n,curtime,iblock) =  &
            TRACER(:,:,k,n,curtime,iblock) - rf_conserve_factor_cur*MASK_BUDGT(:,:,k,iblock)
          enddo ! k

         end do ! block loop (iblock)
        !$OMP END PARALLEL DO

      enddo ! n

      if (lrf_prnt_conserve) call step_RF_diag ('final')

     !*** computation of tracer_mean_initial must precede ice formation
     if (lrf_budget_collect_tracer_mean_initial) &
     call diag_for_tracer_budgets (tracer_mean_initial,volume_t_initial,  &
                                   MASK_BUDGT, bgtarea_t_k, step_call = .true.)

      !***  temporary (?) diagnostics
      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)
      do iblock = 1,nblocks_clinic
        call step_diagnostics(iblock)
      end do ! block loop (iblock)
     !$OMP END PARALLEL DO


     if (check_time_flag(tavg_streams(budget_stream)%field_flag)) then
        !*** time to collect <S_n> for use in this budget interval
        !    note: this collects values computed in this timestep
        call diag_for_tracer_budgets_rf_Sterms (rf_S, rf_volume_total_cur, collect_Sn = .true.)
     endif

      !-----------------------------------------------------------------------
      ! after filtering, recompute and accumulate for time averaging
      !-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)
      do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

        !*** same as standard leapfrog
        FW_OLD(:,:,iblock) = FW(:,:,iblock)

        !*** form ice now, after RF adjustments have been made
        if (liceform .and. mix_pass /= 1) then
            
          !*** form ice from Robert-filtered TRACER(curtime)
          call ice_formation(TRACER(:,:,:,:,curtime,iblock),          &
                             PSURF(:,:,curtime,iblock),               &
                             STF(:,:,1,iblock) + SHF_QSW(:,:,iblock), &
                             iblock,this_block,lfw_as_salt_flx)
        endif
      end do ! block loop (iblock)
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)
      do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

        !*** form ice now, after RF adjustments have been made
        if (liceform .and. mix_pass /= 1) then
            
          !*** always form ice from TRACER(newtime)
          call ice_formation(TRACER(:,:,:,:,newtime,iblock),          &
                             PSURF(:,:,newtime,iblock),               &
                             STF(:,:,1,iblock) + SHF_QSW(:,:,iblock), &
                             iblock,this_block,lfw_as_salt_flx)
        endif
 
        !*** recalculate densities from averaged tracers
        do k = 1,km  

          !*** call state (newtime) every timestep -- is this ok?
!         if (lrf_nonzero_newtime) &
          call state(k,k,TRACER(:,:,k,1,newtime,iblock), TRACER(:,:,k,2,newtime,iblock), &
                         this_block, RHOOUT=RHO(:,:,k,newtime,iblock))

          call state(k,k,TRACER(:,:,k,1,curtime,iblock), TRACER(:,:,k,2,curtime,iblock), &
                     this_block, RHOOUT=RHO(:,:,k,curtime,iblock))
        enddo !k
 
        !*** accumulate TRACERs after RF curtime adjustment
        do k = 1,km  
          call tracer_accumulate_tavg(TRACER (:,:,k,:,oldtime,iblock),  &
                                      TRACER (:,:,k,:,curtime,iblock),  &
                                      RHO    (:,:,k  ,curtime,iblock),  &
                                      PSURF  (:,:    ,curtime,iblock),  &
                                      DH     (:,:            ,iblock),  &
                                      k,iblock)
        enddo !k

      end do ! block loop (iblock)
     !$OMP END PARALLEL DO

     !-----------------------------------------------------------------------
     !  extrapolate next guess for pressure from three known time levels,
     !   after filtering modifications to PSURF at curtime and newtime
     !-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic
       PGUESS(:,:,iblock) = c3*(PSURF(:,:,newtime,iblock)  -  &
                                PSURF(:,:,curtime,iblock)) +  &
                                PSURF(:,:,oldtime,iblock)
      end do
     !$OMP END PARALLEL DO

    !*** update time indices just like standard leapfrog
    tmptime = oldtime
    oldtime = curtime
    curtime = newtime
    newtime = tmptime

    !*** update rf tracer*volume terms
    if (lrf_nonzero_newtime) then
      !*** figure out something later...
        rf_Svol_avg = rf_Svol
    else
      !*** average rf_Svol and rf_Svol_prev for budget diagnostics
        rf_Svol_avg = p5*(rf_Svol+rf_Svol_prev)
        rf_S_avg    = p5*(rf_S+rf_S_prev)

      !*** then update rf_S_prev and rf_Svol_prev
      do n=1,nt
        rf_S_prev   (n) = rf_S(n)
        rf_Svol_prev(n) = rf_Svol(n)
      enddo
    endif

   !*** FW_OLD

   end subroutine step_RF

!***********************************************************************

!BOP
! !IROUTINE: step_RF_diag
! !INTERFACE:

 subroutine step_RF_diag (input_string)

! !DESCRIPTION:
!  Conservation diagnostic <T*volume>_pre-RF = <T*volume>_post-RF
!
! !REVISION HISTORY:
!  added September 2016 njn01

! !INPUT PARAMETERS:

   character (*), intent(in) :: input_string

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   character (char_len) :: control_string
   integer (POP_i4) :: k,n            ! loop indices

   control_string = trim(input_string)

   select case (control_string)

     case ('initial')

       !*** compute <T*volume> pre-Robert-Filtering
       do n=1,nt
         rf_trvol_initial(n) = c0 
         do k=2,km
           rf_trvol_initial(n) = rf_trvol_initial(n)  &
                       + global_sum(TRACER(:,:,k,n,curtime,:)*TAREA(:,:,:)*dz(k),  &
                         distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
         enddo ! k
         k=1
         rf_trvol_initial(n) = rf_trvol_initial(n)  &
                       + global_sum(TRACER(:,:,k,n,curtime,:)*(dz(1) + PSURF(:,:,curtime,:)/grav)*TAREA(:,:,:),  &
                         distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
       enddo ! n

     case ('final')

       !*** compute <T*volume> post-Robert-Filtering
       do n=1,nt
          rf_trvol_final(n) = c0 
          do k=2,km
            rf_trvol_final(n) = rf_trvol_final(n)  &
                        + global_sum(TRACER(:,:,k,n,curtime,:)*TAREA(:,:,:)*dz(k),  &
                          distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
          enddo ! k
          k=1
          rf_trvol_final(n) = rf_trvol_final(n)  &
                        + global_sum(TRACER(:,:,k,n,curtime,:)*(dz(1) + PSURF(:,:,curtime,:)/grav)*TAREA(:,:,:),  &
                          distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
       enddo ! n

       !*** report <T*volume>pre_RF - <T*volume>post_RF
       do n=1,nt
         if (my_task == master_task) write (stdout,*) ' n = ', n
         call document ('step_RF_diag', 'rf_trvol_initial(n)', rf_trvol_initial(n))
         call document ('step_RF_diag', 'rf_trvol_final(n)  ', rf_trvol_final(n))
         call document ('step_RF_diag', 'difference         ', rf_trvol_initial(n)-rf_trvol_final(n))
         call document ('step_RF_diag', ' ')
       enddo ! n

       call document ('step_RF_diag', ' ')

      case default

        !*** silent failure is ok

   end select

   end subroutine step_RF_diag

!***********************************************************************

!BOP
! !IROUTINE: step_RF_doc
! !INTERFACE:

 subroutine step_RF_doc (input_string, nn, rf_conserve_factor_cur )

! !DESCRIPTION:
!  Docuemnt sequence of Robert Filter in step_mod
!
! !REVISION HISTORY:
!  added September 2016 njn01

! !INPUT PARAMETERS:

   character (*), intent(in) :: input_string

   integer (POP_i4), intent(in), optional :: nn
   real    (POP_r8), intent(in), optional :: rf_conserve_factor_cur

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   character (char_len) :: control_string
   integer (POP_i4) :: k,n            ! loop indices

   control_string = trim(input_string)

   select case (control_string)

     case ('begin')

      if (my_task == master_task) write(stdout,*) ' '
      call document ('step_RF', '*** BEGIN lrobert_filter')
      call document ('step_RF', 'nsteps_total', nsteps_total)

     case ('END')

      call document ('step_mod', '==> END of step_mod cycle <==')
      if (my_task == master_task) write(stdout,*) ' '

     case ('conserve')

      if (nn == 1) then
        call document ('step', 'TEMP rf_conserve_factor_cur', rf_conserve_factor_cur)
      else if (nn == 2) then
        call document ('step', 'SALT rf_conserve_factor_cur', rf_conserve_factor_cur)
      else if (nn == 3) then
        call document ('step', 'IAGE rf_conserve_factor_cur', rf_conserve_factor_cur)
      endif
      if (my_task == master_task) write(stdout,*) ' '

      case default

        !*** silent failure is ok

   end select

   end subroutine step_RF_doc

!***********************************************************************

!BOP
! !IROUTINE: init_step
! !INTERFACE:

 subroutine init_step

! !DESCRIPTION:
!  This routine initializes timers and flags used in subroutine step.
!
! !REVISION HISTORY:
!  added 17 August 2007 njn01
!  modified for RF terms January 2016 njn01

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   real (POP_r8) :: rfthick

   integer (POP_i4) :: iblock, k

!-----------------------------------------------------------------------
!
!  initialize timers
!
!-----------------------------------------------------------------------

   call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)
   call get_timer(timer_baroclinic,'BAROCLINIC',1,distrb_clinic%nprocs)
   call get_timer(timer_barotropic,'BAROTROPIC',1,distrb_clinic%nprocs)
   call get_timer(timer_3dupdate,'3D-UPDATE',1,distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!
!  allocate Robert-filter work array
!
!-----------------------------------------------------------------------

   allocate (STORE_RF(nx_block,ny_block,km,nt,nblocks_clinic)) ; STORE_RF = c0
   allocate (VOLSUM_RF(nx_block,ny_block,nt,nblocks_clinic))   ; VOLSUM_RF = c0
   allocate (WORK1(nx_block,ny_block)) ; WORK1 = c0
   allocate (WORK2(nx_block,ny_block)) ; WORK2 = c0

   !*** NOTE:
   !    rf_S_prev and rf_Svol are defined and initialized in time_management.F90 and restart.F90
   allocate (LMASK_BUDGT(nx_block,ny_block,km,nblocks_clinic))
   allocate ( MASK_BUDGT(nx_block,ny_block,km,nblocks_clinic))

   allocate (bgtarea_t_k(km))         ; bgtarea_t_k=c0
   allocate (rf_S(nt))                ; rf_S       =c0  ! Robert filter tracer conservations term
   allocate (rf_Svol(nt))             ; rf_Svol    =c0  ! S*volume 
   allocate (rf_Svol_avg(nt))         ; rf_Svol_avg=c0  ! S*volume average two time levels
   allocate (rf_S_avg(nt))            ; rf_S_avg=c0    
   allocate (rf_sum(km,nt))           ; rf_sum=c0
   allocate (rf_trvol_initial(nt))    ; rf_trvol_initial=c0 
   allocate (rf_trvol_final  (nt))    ; rf_trvol_final=c0 

   lrf_budget_collect_S1                  = .false.
   lrf_budget_collect_tracer_mean_initial = .true.

   rf_volume_2_km = c0
   rf_volsfc_cur  = c0
   rf_volsfc_new  = c0

   !*** LMASK_BUDGT

   !$OMP PARALLEL DO PRIVATE(iblock,k)
   do iblock = 1,nblocks_clinic
     do k=1,km
       if (lrobert_filter) then
         LMASK_BUDGT(:,:,k,iblock) = (KMT(:,:,iblock) >= k .and.  MASK_SR(:,:,iblock) > 0)
       else
         LMASK_BUDGT(:,:,k,iblock) = .true.   ! b4b with nonRF version (partially coupled and avgfit)
       endif
     enddo ! k
   end do ! block loop (iblock)
   !$OMP END PARALLEL DO

     where (LMASK_BUDGT)  
       MASK_BUDGT = c1 
     elsewhere
       MASK_BUDGT = c0
     endwhere

     !*** time-invariant areas and volumes
     write(stdout,*) 'k, area_t_k(k), bgtarea_t_k(k)'

     do k=1,km
       bgtarea_t_k(k) = global_sum(TAREA,distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))
       rfthick = bgtarea_t_k(k)*dz(k)
       if (k .ge. 2) rf_volume_2_km = rf_volume_2_km + rfthick
       write(stdout,1099) k, area_t_k(k),  bgtarea_t_k(k)
     enddo ! k

     write(stdout,*) ' '
     write(stdout,1100)    'rf_volume_2_km', rf_volume_2_km

     !*** time-invariant surface-layer volume
     k=1
     rf_volsfc_invar = global_sum(TAREA(:,:,:)*(dz(k)),  &
                       distrb_clinic,field_loc_center,MASK_BUDGT(:,:,k,:))

1099 format (2x, i4, 2x, 2(1pe25.15))
1100 format (2x, a,  2x, 2(1pe25.15))

   end subroutine init_step
!***********************************************************************

!BOP
! !IROUTINE: step_diagnostics(iblock)
! !INTERFACE:

 subroutine step_diagnostics(iblock)

! !DESCRIPTION:
!  This routine consolidates extra diagnostics 
!
! !REVISION HISTORY:
!  added 11 May 2016 njn01

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: iblock
!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: n
   logical, save :: first_trip = .true.

!-----------------------------------------------------------------------
!
!  misc diagnostics
!
!-----------------------------------------------------------------------

   !*** temporarily, document lrf_nonzero_newtime
   if (first_trip) then
     if (iblock == 1) &
     write(stdout,*) 'lrf_nonzero_newtime = ', lrf_nonzero_newtime
   endif

   !*** collect TEMP at specified levels
   !*** collect TEMP at specified levels
   if (km >= 27) &
      call accumulate_tavg_field(TRACER(:,:,27,1,curtime,iblock),tavg_TEMP_27,iblock,27)
!jt   if (km >= 43) &
!jt      call accumulate_tavg_field(TRACER(:,:,43,1,curtime,iblock),tavg_TEMP_43,iblock,43)

   if (iblock == 1) first_trip = .false.

   end subroutine step_diagnostics

 end module step_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
