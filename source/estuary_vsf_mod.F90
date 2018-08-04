!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module estuary_vsf_mod

! !DESCRIPTION:
!  This module contains the virtual salt flux version of the
!  estuary parameterization
!
!  Algorithm can be found in
!  1. Estuary parameterization:
!     Sun et al., A Box Model for Representing Estuarine Physical Processes in Earth System Models, Ocean. Model, submitted
!  2. Global salinity correction:
!     Tseng et al. (2016), Impacts of the Representation of Riverine Freshwater Input in the Community Earth System Model, Ocean Modell., 105, 71-86.
! !USES:

   use kinds_mod
   use constants,         only: c0, c1, salt_to_ppt, ppt_to_salt, fwmass_to_fwflux
   use communicate,       only: my_task, master_task
   use domain_size,       only: nt, max_blocks_clinic
   use domain,            only: nblocks_clinic, distrb_clinic, blocks_clinic
   use blocks,            only: block, get_block
   use exit_mod,          only: exit_POP, sigAbort
   use io_types,          only: stdout
   use blocks,            only: nx_block, ny_block
   use grid,              only: TAREA, KMT, RCALCT, zw, dz
   use forcing_fields,    only: ROFF_F
   use prognostic,        only: TRACER, curtime
   use global_reductions, only: global_sum, global_sum_prod
   use diagnostics,       only: ldiag_global, DIAG_TRACER_ROFF_VSF_2D, DIAG_TRACER_EXCHCIRC_2D

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_estuary
   public :: set_estuary_vsf_forcing
   public :: set_estuary_exch_circ
   public :: add_estuary_vsf_tracer_tend

   private :: estuary_box_model

! !PUBLIC DATA MEMBERS:

   logical, public :: lestuary_on, lvsf_river, lebm_on

   ! per-tracer correction for using local tracer concenctration in application of ROFF_F
   real(r8), dimension(nt), public :: vsf_river_correction

   ! mask for estuary model, 1 where it runs and 0 elsewhere
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), public :: MASK_ESTUARY

   ! Surface Virtual Salt Flux Associated with Rivers
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), public :: FLUX_ROFF_VSF_SRF

! !PRIVATE DATA MEMBERS:

   ! Maximum layer index to consider estuary contributions
   integer :: kmax_estuary = 15

   ! area where MASK_ESTUARY==1
   real (r8) :: area_vsf_correction

   ! Vertical Salt Flux Across Upper/Lower Layer Interface (From box model)
   REAL (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic) :: &
        FLUX_EXCH_INTRF

   ! Vertical salt fluxes across bottom of a cell interface
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic) ::&
        FLUX_ROFF_VSF_KM1, &
        FLUX_EXCH_KM1

! EBM parameters expanded to 2D fields
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
        & est_tide_amp,&
        & est_mouth_width, &
        & est_mouth_depth, &
        & est_length_a1,&
        & est_tidal_pump_a2, &
        & est_lower_depth_ratio, &
        & est_thick_upper, &
        & est_thick_lower

! EBM input and output fields
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
        & S_upper_layer, &
        & Q_river, Q_lower, Q_upper

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic) :: &
        & tracer_upper_layer, & !Upper-layer tracer value
        & tracer_lower_layer    !Lower-layer tracer value

!-----------------------------------------------------------------------
!
!  tavg ids for tavg diagnostics related to EBM
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &        !all-tracer exchange
         tavg_S_lower_layer,     &!Lower layer salinity
         tavg_S_upper_layer,     &!Upper layer salinity
         tavg_Q_lower,           &!Lower layer flow rate
         tavg_Q_upper             !Upper layer flow rate

   integer (int_kind), dimension(nt) :: &   !all-tracer exchange
         tavg_FLUX_ROFF_VSF_SRF, &!Surf virt salt flux from river
         tavg_FLUX_ROFF_VSF,     &!Vert salt flux across cell level interface
         tavg_FLUX_EXCH_INTRF,   &!Vert tracer flux across layer interface from EBM
         tavg_FLUX_EXCH           !Vert tracer flux across cell level interface from exchange circ.
!EOC
!***********************************************************************

   contains

!***********************************************************************

   subroutine init_estuary

!-----------------------------------------------------------------------
!
!     initialize estuary parameterization
!
!-----------------------------------------------------------------------0

   use constants,    only: blank_fmt, ndelim_fmt, field_loc_center, field_type_scalar
   use domain_size,  only: nx_global, ny_global
   use registry,     only: registry_match
   use forcing_sfwf, only: lfw_as_salt_flx
   use tavg,         only: define_tavg_field
   use grid,         only: RCALCT_OPEN_OCEAN
   use io,           only: data_set
   use io_types,     only: nml_in, nml_filename
   use io_types,     only: datafile, io_field_desc, io_dim
   use io_types,     only: construct_file, rec_type_dbl, construct_io_dim
   use io_types,     only: construct_io_field, destroy_io_field
   use broadcast,    only: broadcast_scalar
   use prognostic,   only: tracer_d

   integer (int_kind) :: &
      i,j,k,n,           &! level index
      nu,                &! unit for input dataset
      nml_error,         &! namelist error flag
      iblock

   type (block) :: this_block ! block information for current block

   !---- Varaibles for global EBM parameter input file
   type(datafile) :: ebm_param_file
   type(io_field_desc) :: io_tide_amp, io_W_h, io_H, io_a1, io_a2, &
        &                 io_h0, io_h_upper, io_h_lower
   type(io_dim) :: i_dim, j_dim
   character(char_len) :: io_field_name
   integer(i4) :: io_field_loc, io_field_type
   real (r8), allocatable, dimension(:,:,:), target :: TEMP_DATA   ! temporary data array

   !Set parameters to the same value for all estuaries (simple test)
   !Eventually replace this with an input file
   real(r8) :: &
        tide_amp, &           ! Tidal amplitude for estuary (m/s), used only if a2>0
        W_h, &                ! Estuary width (m)
        H, &                  ! Estuary mouth depth (m)
        a1, &                 ! Coefficient for estuary length (optimized for Geyer 2010)
        a2, &                 ! Coefficient for tidal pumping  (a2=0 for no tidal pumping)
        h0, &                 ! Ratio of lower-layer depth to H at estuary mouth
        h_upper, &            ! Thickness of upper layer of exchange flow
        h_lower               ! Thickness of lower layer of exchange flow

   character (char_len) ::       &
      estuary_option, &       ! estuary option (on or off)
      estuary_type, &         ! estuary type (vsf or vsf_ebm)
      ebm_param_option, &     ! estuary box model parameter type (internal or file)
      ebm_param_filename, &   ! estuary box model parameter file
      ebm_param_file_fmt      ! estuary box model parameter file format (nc,bin)

   namelist /estuary_nml/ &
        estuary_option, &                               ! 'on', 'off'
        estuary_type, &                                 ! 'vsf', 'vsf_ebm'
        ebm_param_option,     &                         ! 'internal', 'file'
        ebm_param_filename,       &                     ! full pathname for input file
        ebm_param_file_fmt, &                           ! 'nc,bin'
        tide_amp, W_h, H, a1, a2, h0, h_upper, h_lower  ! EBM parameters

!-----------------------------------------------------------------------

   ! Everything off by default
   estuary_option = 'off'
   estuary_type = 'vsf'
   ebm_param_option   = 'internal'
   ebm_param_filename = 'unknown-file'
   ebm_param_file_fmt = 'nc'
!-----------------------------------------------------------------------
!
!     set defaults for EBM parameters, then read them from namelist
!
!-----------------------------------------------------------------------
        tide_amp = 1._r8
        W_h      = 2000._r8
        H        = 10._r8
        a1       = 0.876_r8
        a2       = 0._r8
        h0       = 0.5_r8
        h_upper  = 10._r8
        h_lower  = 10._r8
   ! Look for namelist block
   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=estuary_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading estuary_nml')
   endif

   ! Echo input parameters
   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' Estuary VSF parameterization:'
       write(stdout,blank_fmt)
       write(stdout,*) ' Estuary VSF parametierzation namelist  '
       write(stdout,blank_fmt)
       write(stdout, estuary_nml)
       write(stdout,blank_fmt)
   endif

   call broadcast_scalar(estuary_option,    master_task)
   call broadcast_scalar(estuary_type,      master_task)
   call broadcast_scalar(ebm_param_option,  master_task)
   call broadcast_scalar(ebm_param_filename,master_task)
   call broadcast_scalar(ebm_param_file_fmt,master_task)
   call broadcast_scalar(tide_amp,          master_task)
   call broadcast_scalar(W_h,               master_task)
   call broadcast_scalar(H,                 master_task)
   call broadcast_scalar(a1,                master_task)
   call broadcast_scalar(a2,                master_task)
   call broadcast_scalar(h0,                master_task)
   call broadcast_scalar(h_upper,           master_task)
   call broadcast_scalar(h_lower,           master_task)

   if (estuary_option .EQ. 'on' ) then
      lestuary_on = .true.
   elseif (estuary_option .EQ. 'off' ) then
      lestuary_on = .false.
   else
      call exit_POP(sigAbort,'ERROR estuary_option unknown')
   endif

   if ( lestuary_on ) then
      MASK_ESTUARY(:,:,:) = RCALCT_OPEN_OCEAN(:,:,:)
      area_vsf_correction = global_sum(TAREA, distrb_clinic, field_loc_center, MASK_ESTUARY)
   else
      ! set MASK_ESTUARY to indicate that it doesn't run anywhere
      MASK_ESTUARY(:,:,:) = c0
      area_vsf_correction = c0
   endif

   if ( .NOT. lestuary_on ) then
      lvsf_river = .false.
      return
   endif

   if (.not. lfw_as_salt_flx) then
      call exit_POP(sigAbort, 'estuary_vsf is not implemented for lfw_as_salt_flx=.false.')
   endif

   if (estuary_type == 'vsf' ) then
      lvsf_river = .true.
      lebm_on = .false.
   elseif ( estuary_type == 'vsf_ebm' ) then
      lvsf_river = .true.
      lebm_on = .true.
   else
      call exit_POP(sigAbort,'ERROR estuary_type unknown')
   endif


   if (ebm_param_option == 'file') then

      allocate( TEMP_DATA(nx_block,ny_block,max_blocks_clinic))

      ! Set up the input file
      ebm_param_file = construct_file(ebm_param_file_fmt,&
           &                        full_name=trim(ebm_param_filename),&
           &                        record_length=rec_type_dbl,&
           &                        recl_words=nx_global*ny_global)
      call data_set(ebm_param_file,'open_read')

      i_dim = construct_io_dim('nx',nx_global)
      j_dim = construct_io_dim('ny',ny_global)
      io_field_loc = field_loc_center
      io_field_type = field_type_scalar

      ! Grab the data fields and copy them into the internal 2D arrays (est_* )
      TEMP_DATA = c0
      io_tide_amp = construct_io_field('tide_amp',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_tide_amp)
      call data_set(ebm_param_file,'read',io_tide_amp)
      call destroy_io_field(io_tide_amp)
      est_tide_amp = TEMP_DATA

      TEMP_DATA = c0
      io_W_h = construct_io_field('W_h',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_W_h)
      call data_set(ebm_param_file,'read',io_W_h)
      call destroy_io_field(io_W_h)
      est_mouth_width = TEMP_DATA

      TEMP_DATA = c0
      io_H = construct_io_field('H',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_H)
      call data_set(ebm_param_file,'read',io_H)
      call destroy_io_field(io_H)
      est_mouth_depth = TEMP_DATA

      TEMP_DATA = c0
      io_a1 = construct_io_field('a1',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_a1)
      call data_set(ebm_param_file,'read',io_a1)
      call destroy_io_field(io_a1)
      est_length_a1 = TEMP_DATA

      TEMP_DATA = c0
      io_a2 = construct_io_field('a2',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_a2)
      call data_set(ebm_param_file,'read',io_a2)
      call destroy_io_field(io_a2)
      est_tidal_pump_a2 = TEMP_DATA

      TEMP_DATA = c0
      io_h0 = construct_io_field('h0',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_h0)
      call data_set(ebm_param_file,'read',io_h0)
      call destroy_io_field(io_h0)
      est_lower_depth_ratio = TEMP_DATA

      TEMP_DATA = c0
      io_h_upper = construct_io_field('h_upper',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_h_upper)
      call data_set(ebm_param_file,'read',io_h_upper)
      call destroy_io_field(io_h_upper)
      est_thick_upper = TEMP_DATA

      TEMP_DATA = c0
      io_h_lower = construct_io_field('h_lower',dim1=i_dim,dim2=j_dim,&
           &        field_loc=io_field_loc,field_type=io_field_type,&
           &        d2d_array=TEMP_DATA(:,:,:))
      call data_set(ebm_param_file,'define',io_h_lower)
      call data_set(ebm_param_file,'read',io_h_lower)
      call destroy_io_field(io_h_lower)
      est_thick_lower = TEMP_DATA

   elseif (ebm_param_option == 'internal' ) then
      est_tide_amp = tide_amp
      est_mouth_width = W_h
      est_mouth_depth = H
      est_length_a1 = a1
      est_tidal_pump_a2 = a2
      est_lower_depth_ratio = h0
      est_thick_upper = h_upper*100._r8   !   m -> cm
      est_thick_lower = h_lower*100._r8   !   m -> cm
   else
      call exit_POP(sigAbort,'ERROR ebm_param_option unknown')
   endif

!-----------------------------------------------------------------------
!
!  define tavg fields related to EBM
!
!-----------------------------------------------------------------------
   call define_tavg_field(tavg_S_lower_layer,'S_lower_layer',2,                  &
                          long_name='Lower layer salinity in EBM',               &
                          units='gram/kilogram', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_S_upper_layer,'S_upper_layer',2,                  &
                          long_name='Upper layer salinity in EBM',               &
                          units='gram/kilogram', grid_loc='2110',                &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_Q_lower,'Q_lower',2,                              &
                          long_name='Lower layer flow rate in EBM',              &
                          units='m^3/s', grid_loc='2110',                        &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_Q_upper,'Q_upper',2,                              &
                          long_name='Upper layer flow rate in EBM',              &
                          units='m^3/s', grid_loc='2110',                        &
                          coordinates='TLONG TLAT time')

!QS: define the tavg field for all tracers
   call define_tavg_field(tavg_FLUX_ROFF_VSF_SRF(1),'T_FLUX_ROFF_VSF_SRF',2,     &
                    long_name='Surface Temperature Virtual Salt Flux Associated with Rivers',&
                          units='degC*cm/s', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_FLUX_ROFF_VSF_SRF(2),'S_FLUX_ROFF_VSF_SRF',2,     &
                    long_name='Surface Salt Virtual Salt Flux Associated with Rivers',&
                          units='g/kg*cm/s', grid_loc='2110',                    &
                          scale_factor=1000.0_r8,                                &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_FLUX_EXCH_INTRF(1),'T_FLUX_EXCH_INTRF',2,         &
         long_name='Vertical Temperature Flux Across Upper/Lower Layer Interface (From EBM)',&
                          units='degC*cm/s', grid_loc='2110',                    &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_FLUX_EXCH_INTRF(2),'S_FLUX_EXCH_INTRF',2,         &
         long_name='Vertical Salt Flux Across Upper/Lower Layer Interface (FromEBM)',&
                          units='g/kg*cm/s', grid_loc='2110',                    &
                          scale_factor=1000.0_r8,                                &
                          coordinates='TLONG TLAT time')

    call define_tavg_field(tavg_FLUX_ROFF_VSF(1),'T_FLUX_ROFF_VSF',3,             &
                          long_name='Vertical temperature fluxes across the cell interface from ROFF',  &
                          units='degC*cm/s', grid_loc='3112',                    &
                          coordinates='TLONG TLAT z_w time')

     call define_tavg_field(tavg_FLUX_ROFF_VSF(2),'S_FLUX_ROFF_VSF',3,             &
                          long_name='Vertical salt fluxes across the cell interface from ROFF',  &
                          units='g/kg*cm/s', grid_loc='3112',                    &
                          scale_factor=1000.0_r8,                                &
                          coordinates='TLONG TLAT z_w time')

     call define_tavg_field(tavg_FLUX_EXCH(1),'T_FLUX_EXCH',3,                     &
         long_name='Vertical temperature fluxes across the cell interface from exchange circ',&
                          units='degC*cm/s', grid_loc='3112',                    &
                          coordinates='TLONG TLAT z_w time')

     call define_tavg_field(tavg_FLUX_EXCH(2),'S_FLUX_EXCH',3,                     &
         long_name='Vertical salt fluxes across the cell interface from exchange circ',&
                          units='g/kg*cm/s', grid_loc='3112',                    &
                          scale_factor=1000.0_r8,                                &
                          coordinates='TLONG TLAT z_w time')

   do n=3,nt
     call define_tavg_field(tavg_FLUX_ROFF_VSF_SRF(n),                             &
                            trim(tracer_d(n)%short_name)//'_FLUX_ROFF_VSF_SRF',2,  &
                            long_name='Surface '//trim(tracer_d(n)%short_name)/&
                                      &/' Virtual Salt Flux Associatedwith Rivers',&
                            units=trim(tracer_d(n)%tend_units),                    &
                            scale_factor=tracer_d(n)%scale_factor,                 &
                            grid_loc='2110',                                       &
                            coordinates='TLONG TLAT time')

     call define_tavg_field(tavg_FLUX_EXCH_INTRF(n),                             &
                              trim(tracer_d(n)%short_name)//'_FLUX_EXCH_INTRF',2,&
                             long_name='Vertical '//trim(tracer_d(n)%short_name)/&
                                    &/ ' Flux Across Upper/Lower Layer Interface (From EBM)', &
                             units=trim(tracer_d(n)%tend_units),                 &
                             scale_factor=tracer_d(n)%scale_factor,              &
                             grid_loc='2110',                                    &
                             coordinates='TLONG TLAT time')

       call define_tavg_field(tavg_FLUX_ROFF_VSF(n),                               &
                             trim(tracer_d(n)%short_name)//'_FLUX_ROFF_VSF',3,   &
                          long_name='Vertical '//trim(tracer_d(n)%short_name)/&
                                  &/' fluxes across the cellinterface from ROFF',&
                          units=trim(tracer_d(n)%tend_units),                    &
                          scale_factor=tracer_d(n)%scale_factor,                 &
                          grid_loc='3112',                                       &
                          coordinates='TLONG TLAT z_w time')

       call define_tavg_field(tavg_FLUX_EXCH(n),                                   &
                              trim(tracer_d(n)%short_name)//'_FLUX_EXCH',3,      &
                             long_name='Vertical '//trim(tracer_d(n)%short_name)/&
                                    &/ ' fluxes across the cell interface from exchange circ', &
                             units=trim(tracer_d(n)%tend_units),                 &
                             scale_factor=tracer_d(n)%scale_factor,              &
                             grid_loc='3112',                                    &
                             coordinates='TLONG TLAT z_w time')
   enddo

!-----------------------------------------------------------------------
!
!  ensure that
!  1) est_thick_upper and est_thick_lower are non-zero at the same points
!  2) est_thick_upper is at least as large as the top layer thickness
!
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)
     do j=1,ny_block
       do i=1,nx_block
         if ((est_thick_upper(i,j,iblock) == c0) .neqv. (est_thick_lower(i,j,iblock) == c0)) then
           write(stdout,*) 'i_glob=', this_block%i_glob(i), ', j_glob=', this_block%j_glob(j)
           write(stdout,*) 'est_thick_upper=', est_thick_upper(i,j,iblock)
           write(stdout,*) 'est_thick_lower=', est_thick_lower(i,j,iblock)
           call exit_POP(sigAbort, 'est_thick_upper and est_thick_lower must be non-zero at the same locations')
         endif
         if ((est_thick_upper(i,j,iblock) > c0) .and. (est_thick_upper(i,j,iblock) < dz(1))) then
           write(stdout,*) 'i_glob=', this_block%i_glob(i), ', j_glob=', this_block%j_glob(j)
           write(stdout,*) 'est_thick_upper=', est_thick_upper(i,j,iblock), ', dz(1)=', dz(1)
           call exit_POP(sigAbort, 'est_thick_upper must be at least as large as dz(1)')
         endif
       end do
     end do
   enddo

   end subroutine init_estuary

!***********************************************************************
!BOP
! !IROUTINE: set_estuary_vsf_forcing
! !INTERFACE:

 subroutine set_estuary_vsf_forcing

! !DESCRIPTION:
!  This routine calucates the tracer flux through the sea surface for
!  the virtual flux forcing option. It uses the local surface tracer
!  concentration.
! !REVISION HISTORY:
!  same as module

   use passive_tracers, only: tracer_ref_val
   use forcing_fields,  only: lhas_vflux, vflux_tracer_cnt
   use constants,       only: ocn_ref_salinity, field_loc_center
   use prognostic,      only: tracer_d

! !INPUT/OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: iblock ! block loop index
   integer :: i, j   ! loop indices
   integer :: n      ! tracer loop index
   integer :: n2     ! tracer index for vflux tracers

   real(r8) :: GSUM_FIELD(nx_block,ny_block,vflux_tracer_cnt,nblocks_clinic)
   real(r8) :: gsum_vsf_sloc(vflux_tracer_cnt)
   real(r8) :: glo_ref_val, gsum_vsf_ref1

   type (block) :: this_block ! block information for current block

!-----------------------------------------------------------------------

   ! ensure that ROFF_F == 0 where est_thick_upper == 0 and MASK_ESTUARY == 1
   do iblock=1,nblocks_clinic
     this_block = get_block(blocks_clinic(iblock),iblock)
     do j=1,ny_block
       do i=1,nx_block
         if ((est_thick_upper(i,j,iblock) == c0) .and. (MASK_ESTUARY(i,j,iblock) == c1) .and. (ROFF_F(i,j,iblock) /= c0)) then
           write(stdout,*) 'i_glob=', this_block%i_glob(i), ', j_glob=', this_block%j_glob(j)
           call exit_POP(sigAbort, 'est_thick_upper == 0, MASK_ESTUARY == 1, but ROFF_F /= 0')
         endif
       end do
     end do
   enddo

   ! Compute surface virtual tracer flux using local tracer value
   ! temperature (n=1) is effectively skipped because lhas_vflux(1)==.false.
   !$OMP PARALLEL DO PRIVATE(iblock,n,n2)
   do iblock=1,nblocks_clinic
     n2 = 1
     do n = 1, nt
       if (lhas_vflux(n)) then
         ! fwmass_to_fwflux = 0.1_r8, ROFF_F>0
         ! VSF at surface = -(ROFF*fwmass_to_fwflux)*Su
         ! But set FLUX_ROFF_VSF_SRF>0 at surface (postive upward)
         FLUX_ROFF_VSF_SRF(:,:,n,iblock) = fwmass_to_fwflux*ROFF_F(:,:,iblock)* &
                                           TRACER(:,:,1,n,curtime,iblock)*      &
                                           MASK_ESTUARY(:,:,iblock)
         GSUM_FIELD(:,:,n2,iblock) = TAREA(:,:,iblock)*FLUX_ROFF_VSF_SRF(:,:,n,iblock)
         n2 = n2 + 1
       endif
     enddo
   enddo
   !$OMP END PARALLEL DO

   ! global area-weighted sum of vsf for tracer with global reference value == 1
   gsum_vsf_ref1 = -fwmass_to_fwflux*global_sum_prod(ROFF_F,TAREA(:,:,:),distrb_clinic,&
        & field_loc_center, MASK_ESTUARY)

   ! global area-weighted sum of vsf using local tracer values
   gsum_vsf_sloc(:) = -global_sum(GSUM_FIELD,distrb_clinic,field_loc_center)

   n2 = 1
   do n = 1, nt
     if (lhas_vflux(n)) then
       if (n == 2) then
         ! model salinity units are msu, ocn_ref_salinity units are psu
         glo_ref_val = ppt_to_salt * ocn_ref_salinity
       else
         glo_ref_val = tracer_ref_val(n)
       end if

       vsf_river_correction(n) = (gsum_vsf_ref1*glo_ref_val - gsum_vsf_sloc(n2))/area_vsf_correction
       if (ldiag_global) then
         if ( my_task == master_task) then
           if (n2 == 1) then
             write(stdout,*) 'ESTUARY FORCING (in estuary_vsf_mod)'
           endif
           write(stdout,*)   'tracer = ' // trim(tracer_d(n)%short_name)
           write(stdout,11)  'vsf_tref   (global sum)          ',gsum_vsf_ref1*glo_ref_val
           write(stdout,11)  'vsf_tlocal (global sum)          ',gsum_vsf_sloc(n2)
           write(stdout,11)  'vsf_corr (sum_sref-sum_sloc)/area',vsf_river_correction(n)
         endif
       endif
       n2 = n2 + 1
     endif
   end do

11 format(A,1(1pe11.4,2x))

 end subroutine set_estuary_vsf_forcing

!***********************************************************************
!BOP
! !IROUTINE: set_estuary_exch_circ
! !INTERFACE:

 subroutine set_estuary_exch_circ

   ! !DESCRIPTION:
   !  This routine calculates the salinity flux through the interface
   !  separating the lower and upper layers of the exchange circulation
   !
   ! !REVISION HISTORY:
   !  same as module

   ! !INPUT/OUTPUT PARAMETERS:

   ! !INPUT PARAMETERS:

   !EOP
   !BOC
   !-----------------------------------------------------------------------
   !
   !  local variables
   !
   !-----------------------------------------------------------------------

   real(r8) :: z_layer_top, z_layer_bot    ! depths of top and bottom of lower layer
   real(r8) :: z_cell_top, z_cell_bot    ! depths of top and bottom of grid cell
   real(r8) :: wgt
   integer  :: i,j,n,k,iblock
   real(r8) :: unit_cvrt ! converting factor for the tracer calculation
   real(r8) :: WORK0    ! temporary work space

   !-----------------------------------------------------------------------

   tracer_upper_layer = c0
   tracer_lower_layer = c0
   Q_river = fwmass_to_fwflux*ROFF_F*TAREA*1.0e-6_r8   ! kg/m^2/s -> m^3/s
   FLUX_EXCH_INTRF = c0

   do iblock=1,nblocks_clinic

      do j=1,ny_block
         do i=1,nx_block

            if ( est_thick_lower(i,j,iblock) == c0 .OR. &
                 & Q_river(i,j,iblock) .LE. 0._r8 ) cycle     ! no estuary at this point

            z_layer_top = est_thick_upper(i,j,iblock)
            z_layer_bot = est_thick_upper(i,j,iblock)+est_thick_lower(i,j,iblock)

            ! Compute lower layer averaged salinity
            do k=1,kmax_estuary
               if ( k == 1 ) then
                  z_cell_top = 0._r8
               else
                  z_cell_top = zw(k-1)
               endif
               z_cell_bot = zw(k)

!               if ( z_cell_bot < z_layer_top ) cycle      ! still above lower layer
               if ( z_cell_top > z_layer_bot ) exit       ! now below lower layer

               do n=1,nt

                 IF(n .EQ. 2) THEN ! salinity
                   unit_cvrt=salt_to_ppt
                 ELSE
                   unit_cvrt=c1
                 endif

                 ! Fraction of the upper layer in this cell
                 if ( z_cell_top < z_layer_top ) THEN ! upper layer
                   wgt = ( min(z_cell_bot,z_layer_top) - z_cell_top ) / est_thick_upper(i,j,iblock)

                   tracer_upper_layer(i,j,n,iblock) = tracer_upper_layer(i,j,n,iblock) + &
                                           wgt*TRACER(i,j,k,n,curtime,iblock)*unit_cvrt
                 endif

                 ! Fraction of the lower layer in this cell
                 if ( z_cell_bot > z_layer_top ) THEN ! lower layer
                   wgt = ( min(z_cell_bot,z_layer_bot) - max(z_cell_top,z_layer_top) ) &
                      & / est_thick_lower(i,j,iblock)

                   tracer_lower_layer(i,j,n,iblock) = tracer_lower_layer(i,j,n,iblock) + &
                                           wgt*TRACER(i,j,k,n,curtime,iblock)*unit_cvrt
                 endif

               enddo ! tracer n-loop
            enddo

            ! Run the estuary box model to get the exchange flow volume flux
            call estuary_box_model(Q_river(i,j,iblock),                                       &
                 &                 est_tide_amp(i,j,iblock),tracer_lower_layer(i,j,2,iblock), &
                 &                 est_mouth_width(i,j,iblock), est_mouth_depth(i,j,iblock),  &
                 &                 est_length_a1(i,j,iblock), est_tidal_pump_a2(i,j,iblock),  &
                 &                 est_lower_depth_ratio(i,j,iblock), Q_upper(i,j,iblock),    &
                 &                 Q_lower(i,j,iblock), S_upper_layer(i,j,iblock))

            !  Translate the exchange flow volume flux into a vertical salt flux
            ! NB Q_lower is negative
            do n=1,nt
              if (n .EQ. 2) THEN
                WORK0=(tracer_lower_layer(i,j,n,iblock) - S_upper_layer(i,j,iblock))*ppt_to_salt! psu -> msu
              ELSE
                WORK0=(tracer_lower_layer(i,j,n,iblock) - tracer_upper_layer(i,j,n,iblock))       ! no tracer unit conversion needed
              endif
              FLUX_EXCH_INTRF(i,j,n,iblock) =-Q_lower(i,j,iblock)*(1.e6_r8)* &          ! m^3/s -> cm^3/s
                 &   WORK0 / TAREA(i,j,iblock)*RCALCT(i,j,iblock)
           enddo ! tracer n loop

         enddo ! i: nx_block loop
      enddo ! j: ny_block loop
   enddo ! iblock: nblocks_clinic loop

 end subroutine set_estuary_exch_circ

!***********************************************************************
!BOP
! !IROUTINE: add_estuary_vsf_tracer_tend
! !INTERFACE:

 subroutine add_estuary_vsf_tracer_tend(T_SOURCE, k, this_block)

! !DESCRIPTION:
!  This routine calculates the salinity flux through the bottom of a
!  tracer cell due to the vertically distributed surface virtual
!  salt flux and due to the estuary exchange flow. The salinity tendency
!  is computed as the divergence of the fluxes across the cell.
!  Only salinity is affected. Transport of other tracers by the exchange
!  flow could be easily added later.
!
!  This routine must be called in level-order k=1,2,3...
!
! !REVISION HISTORY:
!  same as module

   use blocks,         only: block
   use grid,           only: DZT, partial_bottom_cells
   use forcing_fields, only: STF_RIV, lhas_vflux, lhas_riv_flux
   use tavg,           only: accumulate_tavg_field

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      T_SOURCE     ! source terms for all tracers (to avoid copies)
                   ! contribution to tendency only added to salinity
                   ! (tracer index 2) for now.

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block info for this block

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,n,bid

   real (r8), dimension(nx_block,ny_block) :: &
      VSF_SRF_p_STF_RIV ,&! ROFF_VSF + STF_RIV
      FLUX_ROFF_VSF     ,&! fluxes across bottom of a cell interface
      FLUX_EXCH         ,&! fluxes across bottom of a cell interface
      WORK1             ,&! temporary work space
      DZ_loc              ! copy of DZ, whether partial_bottom_cells is .true. or .false.
!-----------------------------------------------------------------------
!EOC

      ! Skip all of this if we are beyond the maximum allowable depth
      ! for estuary fluxes
      if ( k > kmax_estuary ) then
         return
      endif

      bid = this_block%local_id

      if ( partial_bottom_cells ) then
        DZ_loc=DZT(:,:,k,bid)
      else
        DZ_loc=DZ(k)
      endif

      if (lvsf_river) then
        do n=1,nt
          if (lhas_vflux(n) .or. lhas_riv_flux(n)) then
            if (lhas_vflux(n)) then
              VSF_SRF_p_STF_RIV(:,:) = FLUX_ROFF_VSF_SRF(:,:,n,bid) ! Virtual salt flux from runoff at sea surface
                                                                    !   masking has been applied in set_estuary_vsf_forcing
            else
              VSF_SRF_p_STF_RIV(:,:) = c0
            endif

            ! add STF_RIV
            if (lhas_riv_flux(n)) then
              ! ensure that STF_RIV == 0 where est_thick_upper == 0 and MASK_ESTUARY == 1
              do j=1,ny_block
                do i=1,nx_block
                  if ((est_thick_upper(i,j,bid) == c0) .and. (MASK_ESTUARY(i,j,bid) == c1) .and. (STF_RIV(i,j,n,bid) /= c0)) then
                    write(stdout,*) 'i_glob=', this_block%i_glob(i), ', j_glob=', this_block%j_glob(j), ', n=', n
                    call exit_POP(sigAbort, 'est_thick_upper == 0, MASK_ESTUARY == 1, but STF_RIV(n) /= 0')
                  endif
                end do
              end do

              ! note that FLUX_ROFF_VSF_SRF is positive up and STF_RIV is positive down
              VSF_SRF_p_STF_RIV(:,:) = VSF_SRF_p_STF_RIV(:,:) - MASK_ESTUARY(:,:,bid)*STF_RIV(:,:,n,bid)
            endif

            !  Set fluxes at top of top cell
            if ( k == 1 ) then
              FLUX_ROFF_VSF_KM1(:,:,n,bid) = VSF_SRF_p_STF_RIV(:,:) ! Virtual salt flux from runoff at sea surface
            endif

            ! Set virtual salt flux at bottom of cell
            ! NB zw(k) always > 0, so should not get divide by zero if est_thick_upper = 0
            do j=1,ny_block
              do i=1,nx_block
                if ( zw(k) <= est_thick_upper(i,j,bid) ) then
                  FLUX_ROFF_VSF(i,j) = VSF_SRF_p_STF_RIV(i,j)*     &
                           &        (c1 - zw(k)/est_thick_upper(i,j,bid))
                else
                  FLUX_ROFF_VSF(i,j) = c0
                endif
              enddo
            enddo

            ! Make sure there is no flux through the sea floor
            where(k >= KMT(:,:,bid) .or. MASK_ESTUARY(:,:,bid) == c0)
              FLUX_ROFF_VSF(:,:) = c0
            endwhere

            ! Compute contribution to tendency from divergence of flux
            ! Surface Virtual salt flux contribution
            ! Positive upward dS/dt=-(F_{k-1}-F_{k})/dz
            WORK1(:,:) = RCALCT(:,:,bid)*(FLUX_ROFF_VSF(:,:) - FLUX_ROFF_VSF_KM1(:,:,n,bid))/DZ_loc
            T_SOURCE(:,:,n) = T_SOURCE(:,:,n) + WORK1(:,:)

            ! compute diagnostics
            if (ldiag_global) then
              where (k <= KMT(:,:,bid))            &
                DIAG_TRACER_ROFF_VSF_2D(:,:,n,bid) = &
                DIAG_TRACER_ROFF_VSF_2D(:,:,n,bid) + &
                WORK1(:,:)*DZ_loc
            endif

            if ((k == 1) .and. (lhas_vflux(n))) then
              call accumulate_tavg_field(FLUX_ROFF_VSF_SRF(:,:,n,bid), &
                                         tavg_FLUX_ROFF_VSF_SRF(n),bid,1)
            endif
            call accumulate_tavg_field(FLUX_ROFF_VSF_KM1(:,:,n,bid),tavg_FLUX_ROFF_VSF(n),bid,k)

            ! Move the bottom vsf to the top of the cell
            FLUX_ROFF_VSF_KM1(:,:,n,bid) = FLUX_ROFF_VSF(:,:)
          endif
        enddo
      endif ! lvsf_river

      if (lebm_on) then
        do n=1,nt
          !  Set fluxes at top of top cell
          if ( k == 1 ) then
            FLUX_EXCH_KM1(:,:,n,bid) = c0                          ! Vertical salt flux across interface between
          endif
          ! Set exchange flow flux at bottom of cell
          do j=1,ny_block
            do i=1,nx_block
              ! Guard against divide by zero
              if ( est_thick_lower(i,j,bid) == c0 ) then
                 FLUX_EXCH(i,j) = c0
                 cycle
              endif

              if ( zw(k) <= est_thick_upper(i,j,bid) ) then
                 FLUX_EXCH(i,j) = FLUX_EXCH_INTRF(i,j,n,bid)*zw(k)/est_thick_upper(i,j,bid)
              elseif ( zw(k) <= est_thick_upper(i,j,bid)+est_thick_lower(i,j,bid) ) then
                 FLUX_EXCH(i,j) =  FLUX_EXCH_INTRF(i,j,n,bid)*   &
                      & (est_thick_upper(i,j,bid)+est_thick_lower(i,j,bid)-zw(k))/est_thick_lower(i,j,bid)
              else
                 FLUX_EXCH(i,j) = c0
              endif

            enddo
          enddo

          ! Make sure there is no flux through the sea floor
          where(k >= KMT(:,:,bid) .or. MASK_ESTUARY(:,:,bid) == c0)
            FLUX_EXCH(:,:) = c0
          endwhere

          ! Exchange flow contribution
          ! Positive upward dS/dt=-(F_{k-1}-F_{k})/dz
          WORK1(:,:) = RCALCT(:,:,bid)*(FLUX_EXCH(:,:) - FLUX_EXCH_KM1(:,:,n,bid))/DZ_loc
          T_SOURCE(:,:,n) = T_SOURCE(:,:,n) + WORK1(:,:)

          ! compute diagnostics
          if (ldiag_global) then
            where (k <= KMT(:,:,bid))              &
              DIAG_TRACER_EXCHCIRC_2D(:,:,n,bid) = &
              DIAG_TRACER_EXCHCIRC_2D(:,:,n,bid) + &
              WORK1(:,:)*DZ_loc
          endif

          if (k == 1)  then
            call accumulate_tavg_field(FLUX_EXCH_INTRF(:,:,n,bid)*MASK_ESTUARY(:,:,bid),  &
                                      tavg_FLUX_EXCH_INTRF(n),bid,1)
            if (n == 2) then
              call accumulate_tavg_field(tracer_lower_layer(:,:,n,bid), &
                                        tavg_S_lower_layer,bid,1)
              call accumulate_tavg_field(S_upper_layer(:,:,bid), &
                                        tavg_S_upper_layer,bid,1)
              call accumulate_tavg_field(Q_lower(:,:,bid), &
                                        tavg_Q_lower,bid,1)
              call accumulate_tavg_field(Q_upper(:,:,bid), &
                                      tavg_Q_upper,bid,1)
            endif
          endif
          call accumulate_tavg_field(FLUX_EXCH_KM1(:,:,n,bid),tavg_FLUX_EXCH(n),bid,k)

          ! Move lower face flux to upper face for next call
          FLUX_EXCH_KM1(:,:,n,bid) = FLUX_EXCH(:,:)
        enddo
      endif ! lebm_on

 end subroutine add_estuary_vsf_tracer_tend

!***********************************************************************

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!No modifications have been made for this application (MMWhitney 7/22/15)

subroutine estuary_box_model (Q_r,tide_amp,S_l,W_h,H,a1,a2,h0,   &
                              Q_u,Q_l,S_u)
! !DESCRIPTION:
!  This subroutine is a test version of estaury Box Model v2_4_FORTRA
!  EBM v2.4 FORTRAN includes river discharge, tides.
!  EBM v2.4 FORTRAN does not include wind mixing, thermal expension
!           and precipitation/evaporation.
!-----------------------------------------------------------------------
! Estuary Box:
!
!             ^ z
!             |____________________         _______
!             |                    |               |
! Q_u,S_u <---+--  upper layer   <-+-- Q_r         |
!             |--------------------|        -+-    | H
!    Q_l,S_l -+->  lower layer     |         | h_l |
!          ___|____________________|        _|_____|
!         x   0                   -LE
!
! NOTE: the negative lower layer volume flux leaves POP, positive one
!      flows into POP.
!-----------------------------------------------------------------------
! Governing Equations:
!   The Estuary Box Model (EBM) is assumed steady state, flat bottom and
! flat surface with straight channel and rectangular cross-section.
!   The EBM is built basic on three global conservations, water mass
! conservation, water volume conservation and potential energy
! conservation. Combine water mass and volume conservation with linear
! equation of state, the salinity flux conservation can be obtained.
! >>> Water mass conservation:
!  rho_r*Q_r = rho_u*Q_u + rho_l*Q_l  + (rho_l - rho_u)*a2*Q_ut/2
! >>> Water volume conservation:
!  Q_r = Q_u +Q_l
! >>> Potential energy conservation (PEF: potential energy flux):
!  PEF_AD + PEF_HD + PEF_VD + PEF_LF = 0
! --> PEF_AD=1/2*g*rho_l*Q_l*h+1/2*g*(rho_u*Q_u-rho_r*Q_r)*(H+h)
! --> PEF_HD=1/2*g*a2*(rho_l-rho_u)Q_ut*(H+h)
! --> PEF_VD=-1/2*g*KV*a1*LE*(rho_l-rho_u)*(rho_l+rho_u-2*rho_r)/(rho_l-rho_r)*W
! --> PEF_LF=-1/4*g*((rho_u^2+2*rho_l*rho_r-2*rho_u*rho_r)*(H-h)-H*rho_r^2+h*rho_l^2)/(rho_l-rho_r)*Q_l*h
!
!-----------------------------------------------------------------------
! QSun, MMWhitney, 2015.04.20

      use constants, only: grav, rho_fw, pi

! input variables

      real (r8), intent(in) :: &
       Q_r,      & !River discharge (m**3/s)
       tide_amp, & !Averaged tidal amplitude at estuary mouth (m)
       S_l,      & !Salinity at estuary lower layer (ppt)
       W_h,      & !Estuary head width (m)
       H,        & !Estuary averaged depth (m)
       a2,       & !A constant (-), of tidal diffusion
       a1,       & !A constant (-), of estuarine mixing length
       h0          !A constant (-), of ratio of geometry: h_l/H

! output variables

      real (r8), intent(out) :: &
       Q_u,      & !Upper layer volume flux (m**3/s, positive)
       Q_l,      & !Lower layer volume flux (m**3/s, negative)
       S_u     !Salinity at estuary upper layer (ppt)

      real (r8) :: &
       ERR_EBMv24  !Closing error of EBMv2.4

! local variables

      real (r8) ::  &
       rho_r,    & !River water density (kg/m**3)
       rho_l,    & !Lower layer inflow density / ocean water density (kg/m**3)
       rho_u,    & !Upper layer inflow density / ocean water density (kg/m**3)
       u_t,      & !Tidal current amplitude near bottom (m/s)
       u_r,      & !Riverine velocity at head of EBM (m/s)
       u_l,      & !Estuarine inflow velocity at mouth of EBM (m/s)
       u_u,      & !Estuarine outflow velocity at mouth of EBM (m/s)
       u_bar,    & !Layer averaged net velocity in EBM (m/s)
       c_wave,   & !Densimetric wave phase speed (m/s)
       ur0,      & !Densimetric riverine Froude number (-)
       ut0,      & !Densimetric tidal current Froude number (-)
       ul0,      & !Densimetric lower layer inflow Froude number (-)
       uu0,      & !Densimetric upper layer outflow Froude number (-)
       R0,       & !Layer densimetric riverine Froude number (-)
       T0,       & !Layer densimetric tidal current Froude number (-)
       h_l,      & !Lower layer water depth of EBM (m)
       a,b,c,d,  & !parameters to solve the cubic equation
       AD,HD,VD,LF  !Different mechanism of salt mass transport in EBM (kg/s)

      integer (int_kind) :: i,n                !Indexes
      integer (int_kind), dimension(3) :: mask !Arrow

      real (r8), dimension(3,2) ::  &
       roots      !roots of 3rd order polynomial:
                  !Dimension one: real part of the root,
                  !Dimension two: imaginary part of the root.

      real (r8) :: &
       g,          &!Gravity in mks (m/s**2)
       rho_0,      &!water reference density (kg/m**3)
       beta_S,     &!Saline contraction coefficient (1/ppt)
       Sc           !Schmidt number (-)
!-----------------------------------------------------------------------
!
!  constants
!
!-----------------------------------------------------------------------
      g      = grav/100.0_r8
      rho_0  = rho_fw*1000.0_r8
      beta_S = 7.7e-4_r8
      Sc     = 2.2_r8
! Executable Statements

!  skip computations if S_l <= 0, using limit as S_l->0
      if (S_l <= 0) then
        Q_u = Q_r
        Q_l = c0
        S_u = c0
        return
      end if

  !----Water density----
      rho_r = rho_0
      rho_l = rho_0*(1.0_r8+beta_S*S_l)
  !----River and tide----
      u_t = -tide_amp*sqrt(g/H) ! Tidal velocity (toward river)
      u_r = Q_r/(W_h*H*(1.0_r8-h0))    ! Riverine velocity at head of upper layer
      c_wave = sqrt(beta_S*S_l*g*H)  ! Densimetric wave phase speed (m/s)

  !----governing equations----
  !solve upper outflow: Q_l
      ur0=u_r/c_wave
      ut0=u_t/c_wave
      R0=ur0*(1.0_r8-h0)
      T0=ut0*(1.0_r8-h0)/pi

      a =-h0**3.0_r8

      b = 2.0_r8*h0**2.0_r8*((2.0_r8-h0)*R0-a2*T0)

      c = 0.096_r8*a1*h0*(Sc**2.0_r8*R0)**(-1.0_r8/3.0_r8)*R0      &
         -h0*((2.0_r8-h0)*R0*(R0-2.0_r8*a2*T0)+a2**2.0_r8*T0**2.0_r8)

      d =-0.048_r8*a1*(Sc**2.0_r8*R0)**(-1.0_r8/3.0_r8)   &
          *R0*(R0-2.0_r8*a2*T0)

      call cubsolve(b/a,c/a,d/a,roots)
  !----find the solution of dimensionless lower layer inflow
      mask=0
      n=0
      do i=1,3
        if (roots(i,1) .LT. 0.0_r8 .AND. roots(i,2) .EQ. 0.0_r8) then
          mask(i)=1
          n=n+1
        end if
      end do

      if (n .eq. 0) then
        write(stdout,100)roots(1,1),roots(1,2),roots(2,1),roots(2,2),roots(3,1),roots(3,2)
100     format('EBM error NO solution: ',6(1pe11.4,1x))
        ul0=0.0_r8
      else if (n .eq. 1) then
        ul0=SUM(roots(1:3,1)*mask)
      else if (n .gt. 1) then
        write(stdout,101)roots(1,1),roots(1,2),roots(2,1),roots(2,2),roots(3,1),roots(3,2)
101     format('EBM error Multiple solution: ',6(1pe11.4,1x))
        ul0=0.0_r8
      end if

  !----solving for the upper layer salinity and volume flux at EBM mouth
      uu0=R0/(1.0_r8-h0)-h0/(1.0_r8-h0)*ul0
      S_u=(-S_l*ul0*h0-S_l*a2*T0)/(R0-ul0*h0-a2*T0)
      Q_l=ul0*h0*H*W_h*c_wave
      Q_u=uu0*(1.0_r8-h0)*H*W_h*c_wave
!      write(stdout,102)Q_r,Q_u,Q_l,S_u,S_l
!102   format('EBM--->POP original: ',5(1pe11.4,1x))

  !----check if the solution of EBMv24 is colsed
      u_l=ul0*c_wave
      u_u=uu0*c_wave
      u_bar=Q_r/(W_h*H)
      h_l=H*h0
      rho_u=rho_0*(1.0_r8+beta_S*S_u)

      AD= 0.5_r8*g*rho_l*u_l*h_l**2_r8   &
         +0.5_r8*g*(rho_u*u_u-rho_r*u_r)*(H**2.0_r8-h_l**2.0_r8)

      HD=-0.5_r8*a2*g*(rho_l-rho_u)*(H**2.0_r8-h_l**2.0_r8)*u_t/pi

      VD=-0.5_r8*g*(rho_u-rho_l)                                    &
          *(rho_l+rho_u-2.0_r8*rho_r)/(rho_l-rho_r)                 &
          *0.024_r8*a1*H**2.0_r8                                    &
          *(c_wave**4.0_r8/(u_bar*Sc**2.0_r8))**(1.0_r8/3.0_r8)

      LF=-0.25_r8*g*Q_l/W_h                                         &
          *( (rho_u**2.0_r8+2.0_r8*rho_l*rho_r-2.0_r8*rho_u*rho_r)  &
            *(H-h_l)                                                &
         -rho_r**2.0_r8*H+rho_l**2.0_r8*h_l)/(rho_l-rho_r)

      ERR_EBMv24=AD-HD-VD-LF
!      write(stdout,103)AD,HD,VD,LF,ERR_EBMv24
!103   format('Salt budget in EBM: ',5(1pe11.4,1x))

! effective salinity for POP using, the Q_l is a negative real value
      S_u=-Q_l*S_l/Q_u
!      write(stdout,104)Q_r,Q_u,Q_l,S_u,S_l
!104   format('EBM--->POP adjusted: ',5(1pe11.4,1x))

end subroutine estuary_box_model

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!No modifications have been made for this application (MMWhitney 7/22/15)

subroutine cubsolve(a,b,c,roots)

! !DESCRIPTION:
!  This subroutine is to solve the 3rd order polynomial analytically.
!  Source of code: https://web.archive.org/web/20141219193637/http://www.helioscorner.com/
!                  http://www.helioscorner.com
! QSun, MMWhitney, 2014.06.18

!     use constants, only: pi

! input variables

      real (r8), intent(in) :: &
       a,b,c     !parameters to solve the cubic equation

! output variables

      real (r8), dimension(3,2), intent(out) :: &
       roots     !roots of 3rd order polynomial:
                 !Column one: real part of the root,
                 !Column two: imaginary part of the root.

! local variables

      real (r8) :: &
       Q,      & !Parameter to solve the cubic equation
       R,      & !Parameter to solve the cubic equation
       Rsqu,   & !Parameter to solve the cubic equation
       Qcub,   & !Parameter to solve the cubic equation
       SQ,     & !Parameter to solve the cubic equation
       theta,  & !Parameter to solve the cubic equation
       X,      & !Parameter to solve the cubic equation
       Y,      & !Parameter to solve the cubic equation
       XY        !Parameter to solve the cubic equation

      real (r8), parameter :: &
       pi = 3.1415927_r8    !Circumference

! Executable Statements
      Q = (a**2.0_r8-3.0_r8*b)/9.0_r8
      R = (2.0_r8*a**3.0_r8-9.0_r8*a*b+27.0_r8*c)/54.0_r8
      Rsqu = R**2.0_r8
      Qcub = Q**3.0_r8

      if (Rsqu .LT. Qcub) then !There are three real roots
        theta=acos(R/sqrt(Qcub))
        SQ=sqrt(Q)
        roots(1,1) = -2.0_r8*SQ*cos(theta/3.0_r8)-a/3.0_r8
        roots(2,1) = -2.0_r8*SQ*cos((theta+2.0_r8*pi)/3.0_r8)-a/3.0_r8
        roots(3,1) = -2.0_r8*SQ*cos((theta-2.0_r8*pi)/3.0_r8)-a/3.0_r8
        roots(1,2) = 0.0_r8
        roots(2,2) = 0.0_r8
        roots(3,2) = 0.0_r8
        return
      endif

      !Otherwise, there are one real root and two conjugate complex roots
      X=-(abs(R)+sqrt(Rsqu-Qcub))**(1.0_r8/3.0_r8)
      if (R .LT. 0.0_r8) then
        X = -X
      endif

      if (X .EQ. 0.0_r8) then
        Y = 0.0_r8
      else
        Y = Q/X
      endif

      XY = X+Y
      roots(1,1) = XY-a/3.0_r8
      roots(1,2) = 0.0_r8

      roots(2,1) = -0.5_r8*XY-a/3.0_r8
      roots(3,1) = -0.5_r8*XY-a/3.0_r8
      roots(2,2) = sqrt(3.0_r8)*(X-Y)/2.0_r8
      roots(3,2) = -sqrt(3.0_r8)*(X-Y)/2.0_r8

end subroutine cubsolve

end module estuary_vsf_mod


