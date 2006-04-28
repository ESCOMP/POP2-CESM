!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing

!BOP
! !MODULE: forcing
!
! !DESCRIPTION:
!  This is the main driver module for all surface and interior
!  forcing.  It contains necessary forcing fields as well as
!  necessary routines for call proper initialization and
!  update routines for those fields.
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use constants
   use blocks
   use distribution
   use domain
   use grid
   use ice, only: salice, tfreez, FW_FREEZE
   use forcing_ws
   use forcing_shf
   use forcing_sfwf
   use forcing_pt_interior
   use forcing_s_interior
   use forcing_ap
   use forcing_coupled
   use forcing_tools
   use passive_tracers, only :      &
       set_passive_tracers_sflux,   &
       init_passive_tracers_sflux,  &
       init_passive_tracers_interior_restore
   use prognostic
   use tavg
   use time_management
   use exit_mod

   !*** ccsm
   use sw_absorption, only: set_chl
   use registry
   use shr_sys_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_forcing,        &
             set_surface_forcing, &
             set_combined_forcing,&
             tavg_forcing

! !PUBLIC DATA MEMBERS:

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), &
      public, target :: &
      SMF,  &!  surface momentum fluxes (wind stress)
      SMFT   !  surface momentum fluxes on T points if avail

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      public, target :: &
      STF,  &!  surface tracer fluxes
      TFW    ! tracer content in freshwater flux

   logical (log_kind), public :: &
      lsmft_avail   ! true if SMFT is an available field

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      public, target ::  &
      IFRAC,             &! ice fraction; not initialized in this routine
      U10_SQR,           &! 10m wind speed squared; not initialized in this routine
      ATM_PRESS,         &! atmospheric pressure forcing
      FW,FW_OLD           ! freshwater flux at T points (cm/s)
                          ! FW_OLD is at time n-1

!EOP
!BOC

   integer (int_kind) :: &
      tavg_SHF,          &! tavg_id for surface heat flux
      tavg_SHF_QSW,      &! tavg_id for short-wave solar heat flux
      tavg_SFWF,         &! tavg_id for surface freshwater flux
      tavg_TAUX,         &! tavg_id for wind stress in X direction
      tavg_TAUY,         &! tavg_id for wind stress in Y direction
      tavg_FW,           &! tavg_id for freshwater flux
      tavg_TFW_T,        &! tavg_id for T flux due to freshwater flux
      tavg_TFW_S          ! tavg_id for S flux due to freshwater flux


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_forcing
! !INTERFACE:

   subroutine init_forcing

! !DESCRIPTION:
!  Initializes forcing by calling a separate routines for
!  wind stress, heat flux, fresh water flux, passive tracer flux,
!  interior restoring, and atmospheric pressure.
!
! !REVISION HISTORY:
!  same as module

!-----------------------------------------------------------------------
!
!  write out header for forcing options to stdout.
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a15)') 'Forcing options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

!-----------------------------------------------------------------------
!
!  initialize forcing arrays
!
!-----------------------------------------------------------------------

   ATM_PRESS = c0
   FW        = c0
   FW_OLD    = c0
   SMF       = c0
   SMFT      = c0
   STF       = c0
   TFW       = c0

!-----------------------------------------------------------------------
!
!  call individual initialization routines
!
!-----------------------------------------------------------------------

   call init_ws(SMF,SMFT,lsmft_avail)

   !*** NOTE: with bulk NCEP forcing init_shf must be called before
   !***       init_sfwf

   call init_shf (STF)
   call init_sfwf(STF)
   call init_pt_interior
   call init_s_interior
   call init_ap(ATM_PRESS)
   if (nt > 2) then
      call init_passive_tracers_sflux
      call init_passive_tracers_interior_restore
   endif
   call init_coupled(SMF, SMFT, STF, SHF_QSW, lsmft_avail)

!-----------------------------------------------------------------------
!
!  error check for coupled option
!
!-----------------------------------------------------------------------

#ifndef coupled
   if (lcoupled) then
      call exit_POP(sigAbort, &
               'ERROR: code must be compiled with coupled ifdef option')
   endif
#endif

!-----------------------------------------------------------------------
!
!  check compatibility of partially-coupled option
!
!-----------------------------------------------------------------------

   if ( .not. lcoupled  .and.                           &
        ( shf_formulation  == 'partially-coupled' .or.  &
          sfwf_formulation == 'partially-coupled' ) ) then
     call exit_POP(sigAbort, &
              'ERROR: partially-coupled option is allowed only when coupled')
   endif

!-----------------------------------------------------------------------
!
!     check coupled compatibility with other forcing options
!
!-----------------------------------------------------------------------

   if (lcoupled) then
     if (ws_data_type /= 'none') then
       call exit_POP(sigAbort, &
                'ws_data_type must be set to none in coupled mode')
     endif
     if ( (shf_formulation  == 'partially-coupled' .and.  &
           sfwf_formulation /= 'partially-coupled') .or.  &
          (shf_formulation  /= 'partially-coupled' .and.  &
           sfwf_formulation == 'partially-coupled') ) then
          call exit_POP(sigAbort, &
                   'partially-coupled must be used for both shf and sfwf')
     endif
     if ( shf_formulation /= 'partially-coupled' .and.  &
          shf_data_type /= 'none') then
       call exit_POP(sigAbort, &
                'shf_data_type must be set to none or '/&
              &/ 'shf_formulation must be partially_coupled when lcoupled is true')
     endif
     if ( sfwf_formulation /= 'partially-coupled' .and.  &
          sfwf_data_type /= 'none') then
       call exit_POP(sigAbort, &
                'sfwf_data_type must be set to none or '/&
             &/ 'sfwf_formulation must be partially_coupled when lcoupled is true')
     endif


     if ( lcoupled .and. shf_formulation /= 'partially-coupled' ) then
       shf_num_comps = 1
       shf_comp_qsw  = 1

       allocate(SHF_COMP(nx_block,ny_block,max_blocks_clinic,shf_num_comps))
       SHF_COMP = c0
      endif



     if ( lcoupled .and. sfwf_formulation /= 'partially-coupled' &
          .and. sfc_layer_type == sfc_layer_varthick .and.       &
          .not. lfw_as_salt_flx .and. liceform ) then

       sfwf_num_comps = 1
       sfwf_comp_cpl  = 1
       tfw_num_comps  = 1
       tfw_comp_cpl   = 1

       allocate(SFWF_COMP(nx_block,ny_block,   max_blocks_clinic,sfwf_num_comps))
       allocate( TFW_COMP(nx_block,ny_block,nt,max_blocks_clinic, tfw_num_comps))

       SFWF_COMP = c0
       TFW_COMP  = c0
     endif

   endif !(lcoupled)

!-----------------------------------------------------------------------
!
!  define tavg diagnostic fields
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_SHF, 'SHF', 2,                                &
                          long_name='Total Surface Heat Flux, Including SW', &
                          missing_value=undefined_nf_r4,                     &
                          units='watt/m^2', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SHF_QSW, 'SHF_QSW', 2,                        &
                          long_name='Solar Short-Wave Heat Flux',            &
                          missing_value=undefined_nf_r4,                     &
                          units='watt/m^2', grid_loc='2110',                 &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_SFWF,'SFWF',2,                                   &
                          long_name='Virtual Salt Flux in FW Flux formulation', &
                          missing_value=undefined_nf_r4,                        &
                          units='kg/m^2/s', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')


   call define_tavg_field(tavg_TAUX,'TAUX',2,                         &
                          long_name='Windstress in grid-x direction', &
                          missing_value=undefined_nf_r4,              &
                          units='dyne/centimeter^2', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_TAUY,'TAUY',2,                         &
                          long_name='Windstress in grid-y direction', &
                          missing_value=undefined_nf_r4,              &
                          units='dyne/centimeter^2', grid_loc='2220', &
                          coordinates='ULONG ULAT time')

   call define_tavg_field(tavg_FW,'FW',2,                        &
                          long_name='Freshwater Flux',           &
                          missing_value=undefined_nf_r4,         &
                          units='centimeter/s', grid_loc='2110', &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_TFW_T,'TFW_T',2,                      &
                          long_name='T flux due to freshwater flux', &
                          missing_value=undefined_nf_r4,             &
                          units='watt/m^2', grid_loc='2110',         &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_TFW_S,'TFW_S',2,                      &
                          long_name='S flux due to freshwater flux (kg of salt/m^2/s)', &
                          missing_value=undefined_nf_r4,             &
                          units='kg/m^2/s', grid_loc='2110',         &
                          coordinates='TLONG TLAT time')

!-----------------------------------------------------------------------
!EOC

 end subroutine init_forcing

!***********************************************************************
!BOP
! !IROUTINE: set_surface_forcing
! !INTERFACE:

 subroutine set_surface_forcing

! !DESCRIPTION:
!  Calls surface forcing routines if necessary.
!  If forcing does not depend on the ocean state, then update
!     forcing if current time is greater than the appropriate
!     interpolation time or if it is the first step.
!  If forcing DOES depend on the ocean state, then call every
!     timestep.  interpolation check will be done within the set\_*
!     routine.
!  Interior restoring is assumed to take place every
!     timestep and is set in subroutine tracer\_update, but
!     updating the data fields must occur here outside
!     any block loops.
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

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      TFRZ               

!-----------------------------------------------------------------------
!
!  Get any interior restoring data and interpolate if necessary.
!
!-----------------------------------------------------------------------

   call get_pt_interior_data
   call get_s_interior_data

!-----------------------------------------------------------------------
!
!  Call individual forcing update routines.
!
!-----------------------------------------------------------------------

   if (lsmft_avail) then
      call set_ws(SMF,SMFT=SMFT)
   else
      call set_ws(SMF)
   endif

   call set_coupled_forcing(SMF,SMFT,STF,SHF_QSW,FW,TFW,IFRAC, &
        ATM_PRESS, U10_SQR)

   call set_chl   !  specify chlorophyll amount for sw absorption
                  !  if-test for chlorophyll is in subroutine set_chl

   !*** NOTE: with bulk NCEP and partially-coupled forcing 
   !***       set_shf must be called before set_sfwf

   call set_shf(STF)
   call set_sfwf(STF,FW,TFW)

   if ( shf_formulation  == 'partially-coupled' .or.  &
        sfwf_formulation == 'partially-coupled' ) then
      call set_combined_forcing(STF,FW,TFW)
   endif


!-----------------------------------------------------------------------
!
! apply solar diurnal cycle if chosen
!
!-----------------------------------------------------------------------

!  index_qsw = mod(nsteps_this_interval,nsteps_per_interval) + 1

!      SHF_QSW = diurnal_cycle_factor(index_qsw) &
!              * SHF_COMP(:,:,:,shf_comp_qsw)

      if ( lcoupled .and. sfwf_formulation /= 'partially-coupled'  &
           .and. sfc_layer_type == sfc_layer_varthick .and.        &
           .not. lfw_as_salt_flx .and. liceform ) then
        FW  = SFWF_COMP(:,:,:,  sfwf_comp_cpl)
        TFW =  TFW_COMP(:,:,:,:, tfw_comp_cpl)
      endif

      if ( sfc_layer_type == sfc_layer_varthick .and.   &
           .not. lfw_as_salt_flx .and. liceform ) then
        FW = FW + FW_FREEZE

        call tfreez(TFRZ,TRACER(:,:,1,2,curtime,:))

        TFW(:,:,1,:) = TFW(:,:,1,:) + FW_FREEZE(:,:,:)*TFRZ(:,:,:)
        TFW(:,:,2,:) = TFW(:,:,2,:) + FW_FREEZE(:,:,:)*salice
      endif


   call set_ap(ATM_PRESS)


   if (nt > 2)  &
      call set_passive_tracers_sflux(SMFT,IFRAC,ATM_PRESS,STF)


   if (ANY(SHF_QSW < c0)) then
      call exit_POP(sigAbort,'ERROR: SHF_QSW < c0 in set_surface_forcing')
   endif


!-----------------------------------------------------------------------
!EOC

 end subroutine set_surface_forcing

!***********************************************************************
!BOP
! !IROUTINE: tavg_forcing
! !INTERFACE:

 subroutine tavg_forcing

! !DESCRIPTION:
!  This routine accumulates tavg diagnostics related to surface
!  forcing.
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

   integer (int_kind) :: &
      iblock              ! block loop index

   type (block) ::       &
      this_block          ! block information for current block

   real (r8), dimension(nx_block,ny_block) :: &
      WORK                ! local temp space for tavg diagnostics

!-----------------------------------------------------------------------
!
!  compute and accumulate tavg forcing diagnostics
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block,WORK)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      if (tavg_requested(tavg_SHF)) then
         where (KMT(:,:,iblock) > 0)
            WORK = (STF(:,:,1,iblock)+SHF_QSW(:,:,iblock))/ &
                   hflux_factor ! W/m^2
         elsewhere
            WORK = c0
         end where

         call accumulate_tavg_field(WORK,tavg_SHF,iblock,1)
      endif

      if (tavg_requested(tavg_SHF_QSW)) then
         where (KMT(:,:,iblock) > 0)
            WORK = SHF_QSW(:,:,iblock)/hflux_factor ! W/m^2
         elsewhere
            WORK = c0
         end where

         call accumulate_tavg_field(WORK,tavg_SHF_QSW,iblock,1)
      endif

      if (tavg_requested(tavg_SFWF)) then
         if (sfc_layer_type == sfc_layer_varthick .and. &
             .not. lfw_as_salt_flx) then
            where (KMT(:,:,iblock) > 0)
               WORK = FW(:,:,iblock)*seconds_in_year*mpercm ! m/yr
            elsewhere
               WORK = c0
            end where
         else
            where (KMT(:,:,iblock) > 0) ! convert to kg(freshwater)/m^2/s
               WORK = STF(:,:,2,iblock)/salinity_factor
            elsewhere
               WORK = c0
            end where
         endif

         call accumulate_tavg_field(WORK,tavg_SFWF,iblock,1)
      endif

      if (tavg_requested(tavg_TAUX)) then
         call accumulate_tavg_field(SMF(:,:,1,iblock), &
                                    tavg_TAUX,iblock,1)
      endif

      if (tavg_requested(tavg_TAUY)) then
         call accumulate_tavg_field(SMF(:,:,2,iblock), &
                                    tavg_TAUY,iblock,1)
      endif

      if (tavg_requested(tavg_FW)) then
         call accumulate_tavg_field(FW(:,:,iblock), &
                                    tavg_FW,iblock,1)
      endif

      if (tavg_requested(tavg_TFW_T)) then
         call accumulate_tavg_field(TFW(:,:,1,iblock)/hflux_factor, &
                                    tavg_TFW_T,iblock,1)
      endif

      if (tavg_requested(tavg_TFW_S)) then
         call accumulate_tavg_field(TFW(:,:,2,iblock)*rho_sw*c10, &
                                    tavg_TFW_T,iblock,1)
      endif



   end do

   !$OMP END PARALLEL DO

   if (lcoupled) call tavg_coupled_forcing

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_forcing

!***********************************************************************

 end module forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
