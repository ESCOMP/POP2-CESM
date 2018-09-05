!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module hmix_gm_share

!BOP
! !MODULE: hmix_gm_share

! !DESCRIPTION:
!  This module contains some of the variables, code, and diagnostics
!  that are common to hmix_gm.F90 and hmix_gm_aniso.F90

! !USES:

   use domain_size, only: nt
   use kinds_mod
   use constants

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_gm_share
   public :: hdifft_gm_share

! !PUBLIC DATA MEMBERS:

   public :: UIT, VIT, WTOP_ISOP, WBOT_ISOP, diag_gm_bolus
   public :: tavg_Redi_TEND_TRACER

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  variables to save from one call to next
!
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:), allocatable :: &
      UIT, VIT,             &! work arrays for isopycnal mixing velocities
      WTOP_ISOP, WBOT_ISOP   ! vertical component of isopycnal velocities

   logical (log_kind) :: &
      diag_gm_bolus          ! true for diagnostic bolus velocity computation

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to diffusivities and
!  isopycnal velocities. Zonal and meridional refer here to logical
!  space only.
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_UISOP,        &   ! zonal      isopycnal velocity
      tavg_VISOP,        &   ! meridional isopycnal velocity
      tavg_WISOP,        &   ! vertical   isopycnal velocity
      tavg_ADVT_ISOP,    &   ! vertically-integrated T eddy-induced
                             !  advection tendency
      tavg_ADVS_ISOP,    &   ! vertically-integrated S eddy-induced
                             !  advection tendency
      tavg_VNT_ISOP,     &   ! heat flux tendency in grid-y direction
                             !  due to eddy-induced velocity
      tavg_VNS_ISOP          ! salt flux tendency in grid-y direction
                             !  due to eddy-induced velocity

   integer (int_kind), dimension(nt) :: &
      tavg_ISOP_ADV_TEND_TRACER,  &! tavg id for eddy-induced advective tendency of tracer
      tavg_Redi_TEND_TRACER        ! tavg id for Redi tendency of tracer

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_gm_share
! !INTERFACE:

   subroutine init_gm_share

! !DESCRIPTION:
!  Initialize variables that are common to hmix_gm.F90 and hmix_gm_aniso.F90
!  that have been migrated to this module.
!
! !REVISION HISTORY:
!  same as module

   use prognostic, only : tracer_d
   use blocks, only: nx_block, ny_block
   use domain, only: nblocks_clinic
   use registry, only: register_string
   use tavg, only: define_tavg_field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n                  ! tracer index

!-----------------------------------------------------------------------
!
!  allocate and initialize bolus velocity arrays if requested
!  define tavg fields related to bolus velocity
!
!-----------------------------------------------------------------------

   if ( diag_gm_bolus ) then

     call register_string ('diag_gm_bolus')

     allocate(WTOP_ISOP(nx_block,ny_block,nblocks_clinic), &
              WBOT_ISOP(nx_block,ny_block,nblocks_clinic), &
                    UIT(nx_block,ny_block,nblocks_clinic), &
                    VIT(nx_block,ny_block,nblocks_clinic))

     WTOP_ISOP = c0
     WBOT_ISOP = c0
     UIT       = c0
     VIT       = c0

     call define_tavg_field (tavg_UISOP, 'UISOP', 3,                &
      long_name='Bolus Velocity in grid-x direction (diagnostic)',  &
                   units='cm/s', grid_loc='3211',                   &
                   coordinates='ULONG TLAT z_t time')

     call define_tavg_field (tavg_VISOP, 'VISOP', 3,                &
      long_name='Bolus Velocity in grid-y direction (diagnostic)',  &
                   units='cm/s', grid_loc='3121',                   &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_WISOP, 'WISOP', 3,                &
      long_name='Vertical Bolus Velocity (diagnostic)',             &
                   units='cm/s', grid_loc='3112',                   &
                   coordinates='TLONG TLAT z_w time')

     call define_tavg_field (tavg_ADVT_ISOP, 'ADVT_ISOP', 2,                            &
      long_name='Vertically-Integrated T Eddy-Induced Advection Tendency (diagnostic)', &
                   units='cm degC/s', grid_loc='2110',                                  &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_ADVS_ISOP, 'ADVS_ISOP', 2,                            &
      long_name='Vertically-Integrated S Eddy-Induced Advection Tendency (diagnostic)', &
                   scale_factor=1000.0_r8,                                              &
                   units='cm gram/kilogram/s', grid_loc='2110',                         &
                   coordinates='TLONG TLAT time')

     call define_tavg_field (tavg_VNT_ISOP, 'VNT_ISOP', 3,                               &
      long_name='Heat Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   units='degC/s', grid_loc='3121',                                      &
                   coordinates='TLONG ULAT z_t time')

     call define_tavg_field (tavg_VNS_ISOP, 'VNS_ISOP', 3,                               &
      long_name='Salt Flux Tendency in grid-y Dir due to Eddy-Induced Vel (diagnostic)', &
                   scale_factor=1000.0_r8,                                               &
                   units='gram/kilogram/s', grid_loc='3121',                             &
                   coordinates='TLONG ULAT z_t time')

      do n = 1,nt
         call define_tavg_field(tavg_ISOP_ADV_TEND_TRACER(n),                &
                                'ISOP_ADV_TEND_' /&
                                        &/ trim(tracer_d(n)%short_name),3,   &
                                long_name='Eddy-induced advective tendency for ' /&
                                        &/trim(tracer_d(n)%short_name),      &
                                units=trim(tracer_d(n)%tend_units),          &
                                scale_factor=tracer_d(n)%scale_factor,       &
                                grid_loc='3111',                             &
                                coordinates='TLONG TLAT z_t time' )

         call define_tavg_field(tavg_Redi_TEND_TRACER(n),                    &
                                'Redi_TEND_' /&
                                        &/ trim(tracer_d(n)%short_name),3,   &
                                long_name='Redi tendency for ' /&
                                        &/trim(tracer_d(n)%short_name),      &
                                units=trim(tracer_d(n)%tend_units),          &
                                scale_factor=tracer_d(n)%scale_factor,       &
                                grid_loc='3111',                             &
                                coordinates='TLONG TLAT z_t time' )
      enddo

   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine init_gm_share

!***********************************************************************
!BOP
! !IROUTINE: hdifft_gm_share
! !INTERFACE:

   subroutine hdifft_gm_share (k, GTK, TMIX, U_ISOP, V_ISOP, FX, FY, FZTOP, &
        tavg_HDIFE_TRACER, tavg_HDIFN_TRACER, tavg_HDIFB_TRACER, this_block)

! !DESCRIPTION:
!  Compute diagnostics that are common to hmix_gm.F90 and hmix_gm_aniso.F90
!  that have been migrated to this module.
!
! !REVISION HISTORY:
!  same as module

   use domain_size, only: km
   use blocks, only: nx_block, ny_block, block
   use time_management, only: mix_pass
   use tavg, only: accumulate_tavg_now, accumulate_tavg_field
   use grid, only: dz, dzr, dz2r
   use grid, only: partial_bottom_cells, DZT
   use grid, only: KMT, HTE, HTN, TAREA_R

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: k  ! depth level index

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      GTK                 ! diffusion+bolus advection for nth tracer at level k

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TMIX                ! tracers at all vertical levels
                          !   at mixing time level

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      U_ISOP, V_ISOP      ! horizontal components of isopycnal velocities

   real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
      FX, FY,            &! fluxes across east, north faces
      FZTOP               ! vertical flux

   integer (int_kind), dimension(nt), intent(in) :: &
      tavg_HDIFE_TRACER, &! tavg id for east face diffusive flux of tracer
      tavg_HDIFN_TRACER, &! tavg id for north face diffusive flux of tracer
      tavg_HDIFB_TRACER   ! tavg id for bottom face diffusive flux of tracer

   type (block), intent(in) :: &
      this_block          ! block info for this sub block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,              &! dummy loop counters
      n,                &! tracer index
      bid                ! local block address for this sub block

   real (r8), dimension(nx_block,ny_block) :: &
      WORK1, WORK2, WORK3

!-----------------------------------------------------------------------
!
!  initialize various quantities
!
!-----------------------------------------------------------------------

   bid = this_block%local_id

!-----------------------------------------------------------------------

   if ( mix_pass /= 1 ) then
        if ( diag_gm_bolus ) then

            call accumulate_tavg_field (U_ISOP, tavg_UISOP, bid, k) 
            call accumulate_tavg_field (V_ISOP, tavg_VISOP, bid, k) 
            call accumulate_tavg_field (WTOP_ISOP(:,:,bid), tavg_WISOP,bid, k)

          if (accumulate_tavg_now(tavg_ADVT_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,1)  &
                      + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )  
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            if (partial_bottom_cells) then
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - DZT(i,j,k,bid) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            else
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            endif

            call accumulate_tavg_field (WORK1, tavg_ADVT_ISOP, bid, k)

          endif

           if (accumulate_tavg_now(tavg_ADVS_ISOP)) then

            WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=1, shift=1) )
            WORK2 = eoshift(WORK1, dim=1, shift=-1)
            WORK3 = WORK1 - WORK2

            WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,2)  &
                      + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
            WORK2 = eoshift(WORK1, dim=2, shift=-1)
            WORK3 = WORK3 + WORK1 - WORK2

            WORK1 = c0
            if (partial_bottom_cells) then
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - DZT(i,j,k,bid) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            else
             do j=this_block%jb,this_block%je
              do i=this_block%ib,this_block%ie
                if ( k <= KMT(i,j,bid) ) then
                  WORK1(i,j) = - dz(k) * TAREA_R(i,j,bid) * WORK3(i,j)
                endif
              enddo
             enddo
            endif

            call accumulate_tavg_field (WORK1, tavg_ADVS_ISOP, bid, k)

          endif

          if ( accumulate_tavg_now(tavg_VNT_ISOP)  .or.  &
               accumulate_tavg_now(tavg_VNS_ISOP) ) then

            WORK1 = p5 * V_ISOP * HTN(:,:,bid) * TAREA_R(:,:,bid) 

            if (accumulate_tavg_now(tavg_VNT_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,1)  &
                         + eoshift(TMIX(:,:,k,1), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNT_ISOP, bid, k) 
            endif

            if (accumulate_tavg_now(tavg_VNS_ISOP)) then
              WORK2 = WORK1 * (    TMIX(:,:,k,2)  &
                         + eoshift(TMIX(:,:,k,2), dim=2, shift=1) )
              call accumulate_tavg_field (WORK2, tavg_VNS_ISOP, bid, k)
            endif

          endif

          ! For each tracer, separate GTK into isop (bolus) and Redi terms, if either term is requested.

          do n = 1,nt
            if ( accumulate_tavg_now(tavg_ISOP_ADV_TEND_TRACER(n)) .or. &
                 accumulate_tavg_now(tavg_Redi_TEND_TRACER(n)) ) then

              ! compute ISOP_ADV_TEND_TRACER, storing into WORK3

              WORK1 = p5 * HTE(:,:,bid) * U_ISOP * ( TMIX(:,:,k,n)  &
                        + eoshift(TMIX(:,:,k,n), dim=1, shift=1) )
              WORK2 = eoshift(WORK1, dim=1, shift=-1)
              WORK3 = WORK1 - WORK2

              WORK1 = p5 * HTN(:,:,bid) * V_ISOP * ( TMIX(:,:,k,n)  &
                        + eoshift(TMIX(:,:,k,n), dim=2, shift=1) )  
              WORK2 = eoshift(WORK1, dim=2, shift=-1)
              WORK3 = WORK3 + WORK1 - WORK2

              WORK3 = -TAREA_R(:,:,bid) * WORK3

              if (k > 1) then
                WORK1 = WTOP_ISOP(:,:,bid) * (TMIX(:,:,k-1,n) + TMIX(:,:,k,n))
              else
                WORK1 = c0
              end if

              if (k < km) then
                WORK2 = WBOT_ISOP(:,:,bid) * (TMIX(:,:,k,n) + TMIX(:,:,k+1,n))
              else
                WORK2 = c0
              end if

              if (partial_bottom_cells) then
                WORK3 = WORK3 - p5 * (WORK1 - WORK2) / DZT(:,:,k,bid)
              else
                WORK3 = WORK3 - (WORK1 - WORK2) * dz2r(k)
              endif

              call accumulate_tavg_field(WORK3,tavg_ISOP_ADV_TEND_TRACER(n),bid,k)

              ! compute Redi_TEND_TRACER as a residual

              WORK3 = GTK(:,:,n) - WORK3

              call accumulate_tavg_field(WORK3,tavg_Redi_TEND_TRACER(n),bid,k)

            endif
          enddo

        endif ! bolus velocity option on

        do n = 1,nt
          if (accumulate_tavg_now(tavg_HDIFE_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FX(i,j,n)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FX(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFE_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFN_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FY(i,j,n)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FY(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFN_TRACER(n),bid,k)
          endif

          if (accumulate_tavg_now(tavg_HDIFB_TRACER(n))) then
            if (partial_bottom_cells) then
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FZTOP(i,j,n)/DZT(i,j,k,bid)*TAREA_R(i,j,bid)
               enddo
               enddo
            else
               do j=this_block%jb,this_block%je
               do i=this_block%ib,this_block%ie
                  WORK1(i,j) = FZTOP(i,j,n)*dzr(k)*TAREA_R(i,j,bid)
               enddo
               enddo
            endif
            call accumulate_tavg_field(WORK1,tavg_HDIFB_TRACER(n),bid,k)
          endif
        enddo
   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine hdifft_gm_share

end module hmix_gm_share

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
