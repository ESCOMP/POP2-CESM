!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module mcog

!BOP
! !MODULE: mcog
!
! !REVISION HISTORY:
! SVN:$Id$

! !USES

   use kinds_mod
   use communicate,   only: my_task, master_task
   use io_types,      only: nml_in, nml_filename, stdout
   use broadcast,     only: broadcast_scalar, broadcast_array
   use exit_mod,      only: sigAbort, exit_pop
   use constants,     only: c0, blank_fmt, ndelim_fmt
   use io_tools,      only: document
   use blocks,        only: nx_block, ny_block
   use domain,        only: nblocks_clinic
   use tavg,          only: define_tavg_field, accumulate_tavg_field

   implicit none
   private
   save

!-----------------------------------------------------------------------
! This module cotains variables, initialization, and tavg routines
! for the handling of shortwave heat flux over subsets of each ocean
! model cell. These subsets are referred to as columns. The module
! name mcog denotes Multi-Column Ocean Grid, a term from the CPT
! (Climate Process Team) project
!
!    Ocean Mixing Processes Associated with High Spatial Heterogeneity
!           in Sea Ice and the Implications for Climate Models
!
! that funded the initial implementation. The initial mcog implementation
! handled heat and freshwater fluxes over seperate columns and was
! coupled to the KPP vertical mixing scheme. The current implementation
! only handles shortwave heat flux over seperate columns and is not
! directly coupled to the ocean model's physics. Notes from the initial
! implementation are included below.
!
! The columns correspond to thickness categories of the CESM sea ice
! model, CICE, and the open ocean portion of the ocean cell. However,
! the primary portion of the code where this matters is in
! POP_CplIndices.F90, where the names of the fields passed to the
! ocean are mapped to column indices. The remainder of the mcog
! implementation is independent of the nature of the columns.
!
! The ocean receives, for each column, the variables
!
!    frac_n        cell fraction
!    fracr_n       fraction for radiative computations
!    fracr_qsw_n   fracr_n * shortwave heat flux
!
! The n suffix denotes the column index, and n ranges from 1 to the
! number of columns. The variables are time means from the previous
! ocean coupling interval. The variables fracr_n are potentially
! distinct from the variables frac_n because of the treatment of
! radiative fluxes in CESM.
!
! mcog aggregates these per column variables into bins. This is of
! use if the columnar resolution provided to the ocean model is not
! necessary in the ocean (e.g. if the sea ice model is configured
! with a large number of thickness categories). The mapping of columns
! to bins is implemented with the mcog_col_to_bin index array. Upon
! receiving the per column variables, mcog generates bin variants of
! the per column variables and then computes
!    qsw_n = fracr_qsw_n / fracr_n
! for each column and bin.
!
! Outside of the mcog internal code, only the binned variables should
! be used. The per column terms are intended to be used only as a
! debugging tool. They are available as tavg output only if lmcog_debug
! is set to .true.
!
! The code checks that the sum of fracr_qsw_n over columns and bins
! agrees with coupler aggregated shortwave heat flux. The code aborts if
! the mismatch exceeds mcog_dagg_qsw_abort_thres, which has a default
! value of 1.0e-10 W/m^2.
!
! The treatment of radiative fluxes in CESM leads to different time
! time lags in frac_n and fracr_n in the sea ice covered and open
! ocean columns. Because of this, the frac_n variables do not necessarily
! sum to 1. The mcog code handles this by adjusting the open ocean
! column fraction, which is column col_frac_adjust_ncol, to ensure
! that the adjusted fractions sum to 1. This adjustment is applied
! after per column and per bin fluxes are normalized by the column
! and bin fractions. The choice for this adjustment would likely need
! to change if the nature of the columns were to change.
!
! fracr_n is not used after this. For columns that arise only from
! sea ice thickness categories, the adjusted fracr_n fractions would
! be identical to the adjusted frac_n fractions.
!
! The subcoupling treatment that is applied to coupler aggregated
! shortwave heat flux (e.g. coszen normalization) is applied to the
! per column and per bin shortwave heat fluxes.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! The following comments are implementation notes from the the initial
! implementation of mcog. They are included for historical reasons
! and to document some of the issues that would arise if this
! implementation is extended to fluxes beyond shortwave heat flux.
! Note that many of the implementation details have changed since the
! initial implementation, so the details below do not necessarily
! apply to the current implementation.
!
! Overview comments for the MCOG CPT (Climate Process Team) project:
!
!    Ocean Mixing Processes Associated with High Spatial Heterogeneity
!           in Sea Ice and the Implications for Climate Models
!
! or otherwise known as MCOG (Multi-Column Ocean Grid) from the term
! used in the proposal description. All code changes for this project
! are denoted by "! MCOG +" and "! MCOG -" comments (for starting and
! ending of modified code).
!
! Here we give a commented overview of the changes made to implement
! MCOG for the ocean component. See the sea ice component and coupler
! for descriptions of fields sent from the sea ice component and ice
! fraction weighted averaging over the coupling interval, as well as
! averaging over the open ocean fraction for the appropriate fluxes.
! (Search for "Overview" in the modified code sub-directories to
! find descriptive comments).
!
! Here we start with the category ice fractions and various ice/ocn
! fluxes and stresses received by the ocean component.
!
! The category fields received by the ocean component from the coupler
! in MCOG are:
!
!            a_n              n_th category sea ice concentration
!            FW_n             n_th category fresh water flux
!            FQ_n             n_th category heat flux
!            FSW_n            n_th category shortwave flux
!            taux_n           n_th category ice/ocn zonal stress
!            tauy_n           n_th category ice/ocn meridional stress
!
! for {n=1,Ncat} categories. All fields are averaged over the coupling
! interval weighted by the time varying ice category fractions, and the
! open ocean (category n=0) fields of shortwave fluxes and ice/ocn
! stresses are computed and sent to the ocean.
!
! Please note that in the sea ice and coupler descriptions of the ice
! fraction weighting for fluxes, equations are given showing the full
! normalized averages. But actually, the sea ice component sends category
! fluxes to the coupler which are already multiplied by category ice
! fraction. The normalization is completed in the ocean component where
! the coupling interval averaged category ice fractions are divided into
! the ice fraction weighted fluxes received.
!
! Specifically, here are the averages:
!                         Nstp
!      a_n = [ a_nm ] = Sum ( a_nm ) / Nstp                            (1)
!                         m=1
!
!                         Nstp                         Ncat
!      A   = [ A_m ] = Sum ( A_m ) / Nstp      A_m = Sum a_nm          (1a)
!                         m=1                          n=1
!                                    Nstp
!      F_n = [ a_nm F_nm ] / a_n = Sum ( a_nm F_nm ) / (a_n Nstp)      (2)
!                                    m=1
!                                     Nstp
!  F_0 = [(1-A_m) Fatm_m] / (1-A) = Sum ((1-A_m) Fatm_m) / (1-A)Nstp   (3)
!                                     m=1
!
! As just said, the sea ice component multiplies each category and
! time step flux (F_nm) by the corresponding category and time step
! sea ice concentration (a_nm). The coupler then performs the usual
! averaging over Nstp, and sends the resulting fluxes to the ocean
! component (i.e. makes averages as in Eq. 1). For the open ocean category
! fluxes which are formed in the coupler, the coupler multiplies by the
! open ocean fraction (i.e. term 1-A_m in Eq.3) before averaging. The ocean
! component then completes the average by dividing with a_n from Eq 1 in
! Eq 2, and by (1-A) as in Eq 3 for category 0 fluxes, where A is the
! total ice fraction given by summing over the coupling interval A_m, which
! is the mth time step total ice fraction, as given in Eq 1a above.
!
! In the ocean, the KPP boundary layer parameterization requires as inputs
! surface forcing and profiles of temperature, salinity and tracers in the
! full ocean column. The temperature, salinity and tracer profiles are for
! the full ocean grid box. Specific surface forcings are the surface friction
! velocity (computed from surface stress), solar and non-solar buoyancy fluxes
! (evaluated from the surface shortwave flux and the sensible/latent heat
! fluxes plus the longwave flux, and finally the fresh water flux), and the
! kinematic surface tracer fluxes for both heat and virtual salt flux. The
! output of KPP are the full column diffusivity and viscosity, which are then
! input to the vertical diffusion solver to update the temperature, salinity
! and tracer profiles.
!
! We note some specifics about the category kinematic surface tracer fluxes
! for heat and virtual salt flux. Contributions to the fluxes for heat
! include atmosphere to ocean fluxes (sensible, latent, up and down longwave)
! which are segregated into the open ocean category (0), snow and ice runoff
! which is included in both open ocean and all cateogries (since there is
! no category specific partition that makes physical sense here), and the
! usual sea ice to ocean fluxes placed in the appropriate category. The
! virtual salt flux surface tracer flux segregates the atmosphere precipitation
! and ocean evaporation into the open ocean category forcing, with land
! and ice runoff placed again in every category, with the sea ice to ocean
! category melt fluxes placed in the appropriate category forcing.
!
! Finally, the coupler multiplies the atmosphere to ocean heat and water
! fluxes by the open ocean fraction, which must be divided back out when
! forming category specific (like open ocean) forcing fluxes.
!
! For the open ocean shortwave absorbed, we note that the ocean model
! typically uses a diurnal cycle partition of the received daily shortwave
! in its time integration. We take the open ocean shortwave absorbed from
! the coupler and partition it diurnally into the open ocean shortwave
! category forcing used by the ocean model. The remainder of the sea ice
! to ocean penetrating shortwave fluxes are placed in the appropriate
! category shortwave forcing.
!
! To apply MCOG to KPP, and minimize overhead of always doing all Ncat+1
! columns, we will run KPP over only those ice-covered columns for which
! {a_n > 0}. Surface forcing is defined using the category fluxes as just
! explained. These will include surface friction velocity u*_n, solar buoyancy
! flux BS_n, non-solar buoyancy flux BNS_n, and kinematic surface tracer fluxes.
! Then, KPP will be run up to N+1 times (depending on ice concentration, which
! sets 0 <= N <= Ncat) producing the multi-column diffusivities k_n and
! viscosities mu_n for {n=0,N}.
!
! We proceed by homogenizing k_n and mu_n for the grid-cell as:
!                      N                         N
!               k = Sum k_n a_n          mu = Sum mu_n a_n                  (4)
!                     n=0                       n=0
! and then run the vertical diffusion solver once to produce modified temperature,
! salinity and tracer profiles. Note that other forcing arrays used outside
! of KPP are also evaluated column by column and aggregrated as in Eq.4, in
! particular the boundary layer depth.
!
! In the implementation here, we have placed a category loop in KPP which runs
! from 0 (i.e. open ocean) through all five (currently, Ncat = 5) sea ice
! categories, and then finish with the originally forced single column ocean grid
! (termed SCOG). When MCOG is prognostic, the code uses the MCOG aggregated terms
! as in Eq 1 (and others as noted), while if it is diagnostic the original SCOG
! terms are used. In this way, we can compare diagnostically MCOG and SCOG from
! history field information regardless of whether MCOG is prognostic or not. We
! save both MCOG and SCOG diffusivities and viscosities (Eq. 1) to history file.
! We also save the individual category terms as well.
!
!   Bruce P. Briegleb   September 2011
!
!-----------------------------------------------------------------------

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_mcog, tavg_mcog

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public ::  &
      lmcog,                      &! .true. if mcog is on, i.e. select per columns fields are passed to the ocean
      lmcog_debug                  ! enable more tavg output if true

   integer, public ::  &
      mcog_ncols,                 &! number of mcog columns (set in POP_CplIndicesSet)
      mcog_nbins,                 &! number of mcog bins (set in init_mcog)
      col_frac_adjust_ncol         ! index of column whose fraction is adjusted, if any (set in POP_CplIndicesSet)

   real (r8), public :: &
      mcog_dagg_qsw_abort_thres    ! call abort if abs(dagg_qsw_raw) exceeds this threshold

   integer, dimension(:), allocatable, public ::  &
      mcog_col_to_bin              ! index mapping of mcog columns to mcog bins

   real (r8), dimension(:,:,:,:), allocatable, public :: &
      FRAC_COL,                   &! fraction of cell occupied by mcog column
      FRAC_BIN,                   &! fraction of cell occupied by mcog bin
      QSW_RAW_COL,                &! raw (directly from cpl) shortwave into each mcog column
      QSW_RAW_BIN,                &! raw (directly from cpl) shortwave into each mcog bin
      QSW_COL,                    &! shortwave into each mcog column, potentially modified by coszen factor
      QSW_BIN                      ! shortwave into each mcog bin, potentially modified by coszen factor

   real (r8), dimension(:,:,:), allocatable, public :: &
      QSW_RAW_COL_DAGG,           &! difference between QSW_RAW aggregated over columns and cpl aggregate
      QSW_RAW_BIN_DAGG             ! difference between QSW_RAW aggregated over columns and cpl aggregate

! !PRIVATE DATA MEMBERS:

   integer (int_kind), dimension(:), allocatable ::  &
      tavg_FRAC_COL,              &! tavg id for FRAC_COL
      tavg_FRAC_BIN,              &! tavg id for FRAC_BIN
      tavg_QSW_RAW_COL,           &! tavg id for QSW_RAW_COL
      tavg_QSW_RAW_BIN,           &! tavg id for QSW_RAW_BIN
      tavg_QSW_COL,               &! tavg id for QSW_COL
      tavg_QSW_BIN                 ! tavg id for QSW_BIN

   integer (int_kind) ::  &
      tavg_QSW_RAW_COL_DAGG,      &! tavg id for QSW_RAW_COL_DAGG
      tavg_QSW_RAW_BIN_DAGG        ! tavg id for QSW_RAW_BIN_DAGG

!EOP
!BOC

!EOC
!***********************************************************************

   contains

!***********************************************************************
!BOP
! !IROUTINE: init_mcog
! !INTERFACE:
 subroutine init_mcog

! !DESCRIPTION:
!  Initializes mcog variables.
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
      ncol, nbin,        &! loop indices
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      sname, lname        ! strings for define_tavg_field calls

   namelist /mcog_nml/ lmcog, lmcog_debug, mcog_dagg_qsw_abort_thres, mcog_col_to_bin

!-----------------------------------------------------------------------
!  tavg defaults
!-----------------------------------------------------------------------

   lmcog = .false.

   lmcog_debug = .false.

   mcog_dagg_qsw_abort_thres = 1.0e-10_r8

   allocate(mcog_col_to_bin(mcog_ncols))

   do ncol = 1, mcog_ncols
      mcog_col_to_bin(ncol) = ncol
   end do

!-----------------------------------------------------------------------
!  read namelist, echo its values, and broadcast results
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=mcog_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP (SigAbort, 'ERROR reading mcog_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' MCOG information'
      write(stdout,blank_fmt)
      write(stdout,*) ' mcog_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,mcog_nml)
      write(stdout,blank_fmt)
   endif

   call broadcast_scalar(lmcog, master_task)
   call broadcast_scalar(lmcog_debug, master_task)
   call broadcast_scalar(mcog_dagg_qsw_abort_thres, master_task)
   call broadcast_array(mcog_col_to_bin, master_task)

   if (.not. lmcog) return

!-----------------------------------------------------------------------
!  ensure that mcog_col_to_bin has valid values
!-----------------------------------------------------------------------

   do ncol = 1, mcog_ncols
      if ((mcog_col_to_bin(ncol) <= 0) .or. (mcog_col_to_bin(ncol) > mcog_ncols)) then
         call document('init_mcog', 'ncol', ncol)
         call document('init_mcog', 'mcog_col_to_bin(ncol)', mcog_col_to_bin(ncol))
         call exit_POP (SigAbort, 'FATAL ERROR: out of range mcog_col_to_bin(ncol) value')
      endif
   end do

   mcog_nbins = maxval(mcog_col_to_bin)
   call document('init_mcog', 'mcog_nbins', mcog_nbins)

!-----------------------------------------------------------------------
!  allocate and initialize mcog arrays
!-----------------------------------------------------------------------

   allocate(FRAC_BIN            (nx_block,ny_block,mcog_nbins,nblocks_clinic))
   allocate(QSW_RAW_BIN         (nx_block,ny_block,mcog_nbins,nblocks_clinic))
   allocate(QSW_BIN             (nx_block,ny_block,mcog_nbins,nblocks_clinic))

   FRAC_BIN    = c0
   QSW_RAW_BIN = c0
   QSW_BIN     = c0

   if (lmcog_debug) then
      allocate(FRAC_COL         (nx_block,ny_block,mcog_ncols,nblocks_clinic))
      allocate(QSW_RAW_COL      (nx_block,ny_block,mcog_ncols,nblocks_clinic))
      allocate(QSW_COL          (nx_block,ny_block,mcog_ncols,nblocks_clinic))
      allocate(QSW_RAW_COL_DAGG (nx_block,ny_block,nblocks_clinic))
      allocate(QSW_RAW_BIN_DAGG (nx_block,ny_block,nblocks_clinic))

      FRAC_COL         = c0
      QSW_RAW_COL      = c0
      QSW_COL          = c0
      QSW_RAW_COL_DAGG = c0
      QSW_RAW_BIN_DAGG = c0
   endif

!-----------------------------------------------------------------------
!  define tavg fields
!-----------------------------------------------------------------------

   allocate(tavg_FRAC_BIN (mcog_nbins))
   allocate(tavg_QSW_BIN  (mcog_nbins))

   do nbin = 1, mcog_nbins
      write(sname, '(a,i2.2)') 'FRAC_BIN_', nbin
      write(lname, '(a,i2.2)') 'fraction of ocean cell occupied by mcog bin ', nbin
      call define_tavg_field(tavg_FRAC_BIN(nbin), trim(sname), 2, &
                             long_name=trim(lname),               &
                             units='1', grid_loc='2110',          &
                             coordinates='TLONG TLAT time')

      write(sname, '(a,i2.2)') 'QSW_BIN_', nbin
      write(lname, '(a,i2.2)') 'net shortwave into mcog bin ', nbin
      call define_tavg_field(tavg_QSW_BIN(nbin), trim(sname), 2, &
                             long_name=trim(lname),              &
                             units='W m-2', grid_loc='2110',     &
                             coordinates='TLONG TLAT time')
   end do

   if (lmcog_debug) then
      allocate(tavg_QSW_RAW_BIN (mcog_nbins))
      allocate(tavg_FRAC_COL    (mcog_ncols))
      allocate(tavg_QSW_RAW_COL (mcog_ncols))
      allocate(tavg_QSW_COL     (mcog_ncols))

      do nbin = 1, mcog_nbins
         write(sname, '(a,i2.2)') 'QSW_RAW_BIN_', nbin
         write(lname, '(a,i2.2)') 'raw (from cpl) net shortwave into mcog bin ', nbin
         call define_tavg_field(tavg_QSW_RAW_BIN(nbin), trim(sname), 2, &
                                long_name=trim(lname),                  &
                                units='W m-2', grid_loc='2110',         &
                                coordinates='TLONG TLAT time')
      end do

      do ncol = 1, mcog_ncols
         write(sname, '(a,i2.2)') 'FRAC_COL_', ncol
         write(lname, '(a,i2.2)') 'fraction of ocean cell occupied by mcog column ', ncol
         call define_tavg_field(tavg_FRAC_COL(ncol), trim(sname), 2, &
                                long_name=trim(lname),               &
                                units='1', grid_loc='2110',          &
                                coordinates='TLONG TLAT time')

         write(sname, '(a,i2.2)') 'QSW_RAW_COL_', ncol
         write(lname, '(a,i2.2)') 'raw (from cpl) net shortwave into mcog column ', ncol
         call define_tavg_field(tavg_QSW_RAW_COL(ncol), trim(sname), 2, &
                                long_name=trim(lname),                  &
                                units='W m-2', grid_loc='2110',         &
                                coordinates='TLONG TLAT time')

         write(sname, '(a,i2.2)') 'QSW_COL_', ncol
         write(lname, '(a,i2.2)') 'net shortwave into mcog column ', ncol
         call define_tavg_field(tavg_QSW_COL(ncol), trim(sname), 2, &
                                long_name=trim(lname),              &
                                units='W m-2', grid_loc='2110',     &
                                coordinates='TLONG TLAT time')

      end do

      sname = 'QSW_RAW_COL_DAGG'
      lname = 'difference between QSW_RAW aggregated over columns and cpl aggregate'
      call define_tavg_field(tavg_QSW_RAW_COL_DAGG, trim(sname), 2, &
                             long_name=trim(lname),                 &
                             units='W m-2', grid_loc='2110',        &
                             coordinates='TLONG TLAT time')

      sname = 'QSW_RAW_BIN_DAGG'
      lname = 'difference between QSW_RAW aggregated over bins and cpl aggregate'
      call define_tavg_field(tavg_QSW_RAW_BIN_DAGG, trim(sname), 2, &
                             long_name=trim(lname),                 &
                             units='W m-2', grid_loc='2110',        &
                             coordinates='TLONG TLAT time')

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_mcog

!***********************************************************************
!BOP
! !IROUTINE: tavg_mcog
! !INTERFACE:

 subroutine tavg_mcog

! !DESCRIPTION:
!   calls tavg accumulation routines for mcog variables
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

   integer (int_kind) :: iblock, nbin, ncol ! loop indices

!-----------------------------------------------------------------------

   if (.not. lmcog) return

   do iblock = 1, nblocks_clinic
      do nbin = 1, mcog_nbins
         call accumulate_tavg_field(FRAC_BIN(:,:,nbin,iblock), tavg_FRAC_BIN(nbin), iblock, 1)
         call accumulate_tavg_field(QSW_BIN(:,:,nbin,iblock), tavg_QSW_BIN(nbin), iblock, 1)
      end do
   end do

   if (lmcog_debug) then
      do iblock = 1, nblocks_clinic
         do nbin = 1, mcog_nbins
            call accumulate_tavg_field(QSW_RAW_BIN(:,:,nbin,iblock), tavg_QSW_RAW_BIN(nbin), iblock, 1)
         end do

         do ncol = 1, mcog_ncols
            call accumulate_tavg_field(FRAC_COL(:,:,ncol,iblock), tavg_FRAC_COL(ncol), iblock, 1)
            call accumulate_tavg_field(QSW_RAW_COL(:,:,ncol,iblock), tavg_QSW_RAW_COL(ncol), iblock, 1)
            call accumulate_tavg_field(QSW_COL(:,:,ncol,iblock), tavg_QSW_COL(ncol), iblock, 1)
         end do

         call accumulate_tavg_field(QSW_RAW_COL_DAGG(:,:,iblock), tavg_QSW_RAW_COL_DAGG, iblock, 1)
         call accumulate_tavg_field(QSW_RAW_BIN_DAGG(:,:,iblock), tavg_QSW_RAW_BIN_DAGG, iblock, 1)
      end do
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine tavg_mcog

end module mcog

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
