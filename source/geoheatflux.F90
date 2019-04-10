!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  module geoheatflux

!BOP
! !MODULE: geoheatflux
! !DESCRIPTION:
!  This module contains geothermal-heat-flux-related variables and methods.
!
! !REVISION HISTORY:
!  SVN:$Id: geoheatflux.F90 $

! !USES:

   use kinds_mod
   use blocks
   use broadcast
   use constants
   use diagnostics
   use distribution
   use domain
   use exit_mod
   use forcing_shf, only: MASK_SR
   use global_reductions
   use grid
   use io
   use io_tools
   use registry
   use tavg
   use time_management

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_geoheatflux
   public :: geoheatflux_accumulate_tavg_once

! !PUBLIC DATA MEMBERS:

   public :: GEOHFLUX
   public :: geoheatflux_const
   public :: geoheatflux_depth
   public :: geo_mean
   public :: lgeoheatflux
   public :: tavg_GEOHFLUX
   public :: tavg_GEOHFLUX_3D

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  geothermal (bottom) heat flux
!
!-----------------------------------------------------------------------

   integer (int_kind)   :: geoheatflux_itype
   integer (int_kind)   :: geoheatflux_itype_zero          = 0
   integer (int_kind)   :: geoheatflux_itype_nonzero_const = 1
   integer (int_kind)   :: geoheatflux_itype_spatial       = 2
   integer (int_kind)   :: geoheatflux_itype_error         = -1000

   logical (log_kind)   :: lgeoheatflux ! true if geothermal heat flux /= 0

   real (r8), allocatable, dimension (:,:,:) ::  &
      GEOHFLUX            ! spatial array of geothermal heat flux

   real (r8)  :: geoheatflux_const  ! geothermal heat flux (W/m2)
   real (r8)  :: geoheatflux_depth  ! depth below which geothermal heat applied

!-----------------------------------------------------------------------
!
!  geothermal heat flux budgets
!
!-----------------------------------------------------------------------

   real (r8), allocatable, dimension(:,:,:,:)  ::  MASK_TRBUDGET_GEO

!-----------------------------------------------------------------------
!
!  geothermal heat flux diagnostics
!
!-----------------------------------------------------------------------

   real (r8) :: geo_mean

   integer (int_kind) :: tavg_GEOHFLUX    ! tavg id for geothermal heat flux
   integer (int_kind) :: tavg_GEOHFLUX_3D ! tavg id for geothermal heat flux*3D MASK

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_geoheatflux
! !INTERFACE:

 subroutine init_geoheatflux

! !DESCRIPTION:
!  Initialization for geothermal heat flux
!
! !REVISION HISTORY:
!  njn01 2016

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock    ! block index
   integer (int_kind) :: i, j      ! loop indices
   integer (int_kind) :: k         ! vertical level index
   integer (int_kind) :: n         ! dummy loop index
   integer (int_kind) :: ncol      ! number of columns
   integer (int_kind) :: nml_error ! namelist i/o error flag
   integer (int_kind) :: nu        ! i/o unit

   character (char_len) :: geoheatflux_choice ! input choice for desired geothermal heat-flux parameterization
   character (char_len) :: geo_string

   namelist /geoheatflux_nml/ geoheatflux_choice, geoheatflux_const, geoheatflux_depth


!-----------------------------------------------------------------------
!
!  set defaults, then read input namelist 
!
!-----------------------------------------------------------------------

   geoheatflux_choice = 'unknown_geoheatflux_choice'
   geoheatflux_const  = c0
   geoheatflux_depth  = 1000.00e2_r8

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=geoheatflux_nml, iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)

   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading geoheatflux_nml') 
   endif

   if (my_task == master_task) then
     write(stdout,blank_fmt)
     write(stdout,ndelim_fmt)
     write(stdout,blank_fmt)
     write(stdout,'(a27)') 'Geothermal Heat Flux options'
     write(stdout,blank_fmt)
     write(stdout,geoheatflux_nml)
     write(stdout,delim_fmt)
   endif

   call broadcast_scalar(geoheatflux_choice, master_task)
   call broadcast_scalar(geoheatflux_const,  master_task)
   call broadcast_scalar(geoheatflux_depth,  master_task)

!-----------------------------------------------------------------------
!
!  define geothermal heat flux and logical controls
!
!-----------------------------------------------------------------------

   lgeoheatflux = .false.

   if (trim(geoheatflux_choice) == 'const') then
     allocate ( GEOHFLUX(1,1,1))

     if (geoheatflux_const == c0) then
       geoheatflux_itype = geoheatflux_itype_zero
       lgeoheatflux = .false.
       geo_string   = 'zero' 
     else
       geoheatflux_itype = geoheatflux_itype_nonzero_const
       lgeoheatflux = .true.
       geo_string   = 'non-zero constant' 

       call register_string ('geoheatflux_choice_nonzero_const')

       !*** convert from W/m^2 to cgs:    
       geoheatflux_const = geoheatflux_const * hflux_factor
     endif

   else if (trim(geoheatflux_choice) == 'spatial') then
     allocate ( GEOHFLUX(nx_block,ny_block,nblocks_clinic))

     geoheatflux_itype = geoheatflux_itype_spatial
     lgeoheatflux = .true.
     geo_string   = 'spatial'

     call register_string ('geoheatflux_choice_spatial')

     !*** compute GEOHFLUX
     call form_geoheatflux

     !*** convert from mW/m^2 to cgs:    
     GEOHFLUX = GEOHFLUX*hflux_factor/c1000

   else
      geoheatflux_itype = geoheatflux_itype_error
      lgeoheatflux = .false.
      geo_string   = 'error condition -- model will stop'
   endif

!-----------------------------------------------------------------------
!
!  geothermal heat flux budget mask
!
!-----------------------------------------------------------------------
   allocate (MASK_TRBUDGET_GEO(nx_block,ny_block,km,nblocks_clinic))
   MASK_TRBUDGET_GEO = c0

  !$OMP PARALLEL DO PRIVATE(iblock,k)
   do iblock = 1,nblocks_clinic
     do k=1,km
      if (zw(k) >= geoheatflux_depth) then
      do j=1,ny_block
      do i=1,nx_block
        if (k == KMT(i,j,iblock)) then
         !write(stdout,1313) '(init_geoheatflux) BINGO! iblock, i,j,k, KMT(i,j,iblock) = ', &
         !                                              iblock, i,j,k, KMT(i,j,iblock)
!########### DEBUG
!         MASK_TRBUDGET_GEO(i,j,k,iblock) = c1
          MASK_TRBUDGET_GEO(i,j,k,iblock) = RCALCT_OPEN_OCEAN_3D(i,j,k,iblock)
!########### DEBUG
        endif ! k KTMI
      enddo ! j
      enddo ! i
      endif ! zw
     enddo ! k
   end do ! block loop (iblock)
   !$OMP END PARALLEL DO
1313 format(1x,a, 5(1x,i2.2))

!-----------------------------------------------------------------------
!
!  document geothermal heat flux options
!
!-----------------------------------------------------------------------

   call document ('init_geoheatflux', 'geothermal heat flux option', geo_string)

   if (geoheatflux_itype == geoheatflux_itype_zero) then
     call document ('init_geoheatflux', 'geothermal (bottom) heat flux is disabled')
   elseif (geoheatflux_itype == geoheatflux_itype_nonzero_const) then
     call document ('init_geoheatflux', 'geothermal (bottom) heat flux is enabled')
     call document ('init_geoheatflux', 'constant geothermal (bottom) heat flux (cgs)',  &
                    geoheatflux_const)
   elseif (geoheatflux_itype == geoheatflux_itype_spatial) then
     call document ('init_geoheatflux', 'Spatially varying geothermal heat flux')
   endif

   if (lgeoheatflux) then
     call document ('init_geoheatflux','geothermal heat flux applied below (cm)',geoheatflux_depth)
   endif
  
!-----------------------------------------------------------------------
!
!  perform some error and consistency checks
!
!-----------------------------------------------------------------------

  if (geoheatflux_itype == geoheatflux_itype_error) then
    geo_string = trim(geoheatflux_choice)
    write(stdout,*) 'geoheatflux_choice = ', geo_string
    call document ('init_geoheatflux', 'geothermal heat flux option', geo_string)
    call exit_POP(sigAbort, '(init_geoheatflux) unknown geoheatflux_choice')
  endif

  if (geoheatflux_itype == geoheatflux_itype_spatial .and. partial_bottom_cells) then
    geo_string = trim(geoheatflux_choice)
    write(stdout,*) 'geoheatflux_choice = ', geo_string
    call document ('init_geoheatflux', 'geothermal heat flux option', geo_string)
    geo_string = 'spatial geothermal heat flux not yet implemented with partial bottom cells'
    call exit_POP(sigAbort, '(init_geoheatflux) unknown geoheatflux_choice')
  endif


!-----------------------------------------------------------------------
!
!  define fields for accumulating tavg diagnostics
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_GEOHFLUX,'GEOHFLUX',2,   &
        long_name='geothermal (bottom) heat flux data', &
        tavg_method=tavg_method_constant,               &
        units='mW/s^2', grid_loc='2110',                &
        coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_GEOHFLUX_3D,'GEOHFLUX',3,    &
        long_name='geothermal (bottom) heat flux by level', &
        tavg_method=tavg_method_constant,                   &
        units='mW/s^2', grid_loc='3111',                    &
        coordinates='TLONG TLAT z_t time')

!-----------------------------------------------------------------------
!EOC

 call flushm(stdout)

 end subroutine init_geoheatflux

!***********************************************************************
!BOP
! !IROUTINE: form_geoheatflux
! !INTERFACE:

 subroutine form_geoheatflux

! !DESCRIPTION:
!     This subroutine was contributed to CESM by Andreas Schmittner and
!     then lightly modified ("POP-ified") by njn01 Nov 2016
!   
!     The following comments are from Schmittner's code:
!     source file: /raid24/aschmitt/UVic2.9/MOBI1.7/updates/bhf.F
!     written by Alex Hoffman
!     put into 2.9 by Andreas Schmittner
!     expects latitude and longitude in degrees!
!     p is longitude and q is colatitude

! !REVISION HISTORY:
!  njn01 2016

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

      real (r8) ::  anm(0:12,0:12), bnm(0:12,0:12), p, con, gcon
      real (r8) ::  q, x, krt, sum, pprime, pp, qq
      real (r8) ::  xt, yt
      real (r8) ::  cosyt, sinyt
      real (r8) ::  twotothen, sinyttothem, twonp1, qq_term

      integer (int_kind) ::  i, j, iblock
      integer (int_kind) ::  h, m, n, t

!-----------------------------------------------------------------------
!
!  use coefficients from  Pollack, Henry N and Suzanne Hurter,
!   _Heat Flow from the Earth's !Interior: Analysis of the Global Data Set_,
!   Reviews of Geophysics, 31, 267-280, 1993.
!
!-----------------------------------------------------------------------

      anm(0,0)  = 86.674_r8
      anm(1,0)  =-12.999_r8
      anm(1,1)  = -2.689_r8
      bnm(1,1)  =-10.417_r8
      anm(2,0)  = -1.917_r8
      anm(2,1)  =  4.578_r8
      bnm(2,1)  =  1.022_r8
      anm(2,2)  =-14.076_r8
      bnm(2,2)  =  6.507_r8
      anm(3,0)  =  7.122_r8
      anm(3,1)  = -2.934_r8
      bnm(3,1)  =  3.555_r8
      anm(3,2)  =  7.232_r8
      bnm(3,2)  = -3.295_r8
      anm(3,3)  = 10.299_r8
      bnm(3,3)  =  4.646_r8
      anm(4,0)  = -3.511_r8
      anm(4,1)  =  2.778_r8
      bnm(4,1)  = -1.873_r8
      anm(4,2)  =  1.728_r8
      bnm(4,2)  = -2.546_r8
      anm(4,3)  = -4.822_r8
      bnm(4,3)  =   .486_r8
      anm(4,4)  =  4.408_r8
      bnm(4,4)  =-17.946_r8
      anm(5,0)  =  5.316_r8
      anm(5,1)  = -1.984_r8
      bnm(5,1)  = -2.642_r8
      anm(5,2)  =  2.167_r8
      bnm(5,2)  =  3.835_r8
      anm(5,3)  =  4.570_r8
      bnm(5,3)  = -6.087_r8
      anm(5,4)  = -8.353_r8
      bnm(5,4)  = 10.283_r8
      anm(5,5)  = -6.896_r8
      bnm(5,5)  = -4.199_r8
      anm(6,0)  = -5.204_r8
      anm(6,1)  =  2.795_r8
      bnm(6,1)  =  3.162_r8
      anm(6,2)  =  2.065_r8
      bnm(6,2)  = -2.889_r8
      anm(6,3)  = -2.740_r8
      bnm(6,3)  =  -.252_r8
      anm(6,4)  =  -.012_r8
      bnm(6,4)  = -1.897_r8
      anm(6,5)  =   .637_r8
      bnm(6,5)  =   .476_r8
      anm(6,6)  =  3.739_r8
      bnm(6,6)  =  7.849_r8
      anm(7,0)  =  2.010_r8
      anm(7,1)  =   .912_r8
      bnm(7,1)  =   .116_r8
      anm(7,2)  = -6.044_r8
      bnm(7,2)  =  -.179_r8
      anm(7,3)  =  4.999_r8
      bnm(7,3)  =  -.123_r8
      anm(7,4)  = -1.605_r8
      bnm(7,4)  = -3.721_r8
      anm(7,5)  =  -.334_r8
      bnm(7,5)  =  3.466_r8
      anm(7,6)  = -4.111_r8
      bnm(7,6)  =  -.639_r8
      anm(7,7)  =  4.126_r8
      bnm(7,7)  = -1.659_r8
      anm(8,0)  =  2.621_r8
      anm(8,1)  = -1.376_r8
      bnm(8,1)  =  1.795_r8
      anm(8,2)  =  7.201_r8
      bnm(8,2)  =  1.436_r8
      anm(8,3)  = -1.947_r8
      bnm(8,3)  =   .679_r8
      anm(8,4)  =   .204_r8
      bnm(8,4)  =  1.171_r8
      anm(8,5)  =  1.851_r8
      bnm(8,5)  =  1.771_r8
      anm(8,6)  =  3.579_r8
      bnm(8,6)  =  -.250_r8
      anm(8,7)  =  1.886_r8
      bnm(8,7)  =  4.903_r8
      anm(8,8)  = -5.285_r8
      bnm(8,8)  = -4.412_r8
      anm(9,0)  =  -.211_r8
      anm(9,1)  =  3.140_r8
      bnm(9,1)  =   .886_r8
      anm(9,2)  =  -.360_r8
      bnm(9,2)  = -3.894_r8
      anm(9,3)  = -3.004_r8
      bnm(9,3)  = -2.056_r8
      anm(9,4)  =  1.947_r8
      bnm(9,4)  = -2.511_r8
      anm(9,5)  =   .328_r8
      bnm(9,5)  = -3.064_r8
      anm(9,6)  =  1.030_r8
      bnm(9,6)  =  -.745_r8
      anm(9,7)  = -4.117_r8
      bnm(9,7)  = -3.888_r8
      anm(9,8)  =  6.529_r8
      bnm(9,8)  =  3.889_r8
      anm(9,9)  = -4.084_r8
      bnm(9,9)  =  -.082_r8
      anm(10,0) =  2.735_r8
      anm(10,1) = -1.624_r8
      bnm(10,1) = -1.998_r8
      anm(10,2) = -1.309_r8
      bnm(10,2) =  1.333_r8
      anm(10,3) =  4.576_r8
      bnm(10,3) =   .641_r8
      anm(10,4) = -4.506_r8
      bnm(10,4) =   .927_r8
      anm(10,5) =  -.363_r8
      bnm(10,5) =  -.927_r8
      anm(10,6) = -4.528_r8
      bnm(10,6) = -1.353_r8
      anm(10,7) =  -.952_r8
      bnm(10,7) =  1.810_r8
      anm(10,8) = -1.104_r8
      bnm(10,8) =  -.739_r8
      anm(10,9) =   .129_r8
      bnm(10,9) =   .644_r8
      anm(10,10)=  4.164_r8
      bnm(10,10)= -3.463_r8
      anm(11,0) = -1.708_r8
      anm(11,1) =   .429_r8
      bnm(11,1) =  2.902_r8
      anm(11,2) =  2.106_r8
      bnm(11,2) =   .915_r8
      anm(11,3) = -5.078_r8
      bnm(11,3) =   .595_r8
      anm(11,4) =  3.441_r8
      bnm(11,4) =   .907_r8
      anm(11,5) =   .784_r8
      bnm(11,5) =  2.762_r8
      anm(11,6) =   .158_r8
      bnm(11,6) =   .782_r8
      anm(11,7) =  -.377_r8
      bnm(11,7) =  -.355_r8
      anm(11,8) =  -.818_r8
      bnm(11,8) =  1.851_r8
      anm(11,9) =  3.654_r8
      bnm(11,9) =  1.336_r8
      anm(11,10)= -1.765_r8
      bnm(11,10)=  4.245_r8
      anm(11,11)=  -.505_r8
      bnm(11,11)= -3.520_r8
      anm(12,0) =  1.003_r8
      anm(12,1) =  -.689_r8
      bnm(12,1) = -1.476_r8
      anm(12,2) = -2.359_r8
      bnm(12,2) =  -.066_r8
      anm(12,3) =  3.863_r8
      bnm(12,3) =   .504_r8
      anm(12,4) =   .793_r8
      bnm(12,4) = -1.034_r8
      anm(12,5) = -1.761_r8
      bnm(12,5) =  -.267_r8
      anm(12,6) =  2.439_r8
      bnm(12,6) = -2.484_r8
      anm(12,7) = -2.080_r8
      bnm(12,7) =  3.714_r8
      anm(12,8) =  2.237_r8
      bnm(12,8) =   .809_r8
      anm(12,9) =   .289_r8
      bnm(12,9) =  -.838_r8
      anm(12,10)=  1.516_r8
      bnm(12,10)= -4.821_r8
      anm(12,11)=  4.114_r8
      bnm(12,11)=  -.533_r8
      anm(12,12)= -3.033_r8
      bnm(12,12)=  2.175_r8

!-----------------------------------------------------------------------
!
!  loop over all points in this block
!
!-----------------------------------------------------------------------

 do iblock = 1,nblocks_clinic

   do i=1,nx_block
   do j=1,ny_block

!-----------------------------------------------------------------------
!
!  use functions as defined in Appendix 1 of Hamza et al. 2007
!
!-----------------------------------------------------------------------

      !*** co-latitude and longitude in degrees
      yt=(90._r8-TLATD(i,j,iblock))/radian
      xt=TLON(i,j,iblock)

      cosyt = cos(yt)
      sinyt = sin(yt)
      qq=c0

      do n=0,12
         twotothen=2**n
         twonp1   =2*n+1._r8
         do m=0,n
            sinyttothem = sinyt**m
            qq_term = anm(n,m)*cos(m*xt)+bnm(n,m)*sin(m*xt)
            sum=c0
            do t=0,int((n-m)/2)

               sum = ((((-1)**t)*factorial(2*n-2*t))/(factorial(t)*  &
               factorial(n-t)*factorial(n-m-2*t))*(cosyt**(n-m-2*t)  &
                    ) )+sum

            enddo
            pprime = (sinyttothem/twotothen)*sum 
            if (m.eq.0) then
              krt = sqrt((factorial(n+m)/factorial(n-m))/(twonp1))
            else
               h=2._r8
              krt = sqrt(p5*(factorial(n+m)/factorial(n-m))/(twonp1))
            endif
            pp = (pprime)/(krt)

            qq = (qq_term*pp)+qq

         enddo
      enddo

!-----------------------------------------------------------------------
!
!  Geothermal heat flux in mW/m^2  (convert to cgs in calling routine)
!
!-----------------------------------------------------------------------

      GEOHFLUX(i,j,iblock) = qq

     enddo ! j
     enddo ! i

 enddo ! iblock

 end subroutine form_geoheatflux

!***********************************************************************
!BOP
! !IROUTINE: factorial
! !INTERFACE:

 real function factorial(n)

! !DESCRIPTION:
!  This function was contributed by Andreas Schmittner and
!  lightly modified/"POP-ified" by njn01

! !REVISION HISTORY:
!  njn01 2016

! !INPUT PARAMETERS:
   integer (int_kind), intent(in) :: n

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::  i
   real (r8) :: f ! factorial function value

      f=1._r8
      do i=1,n
         f=i*f
      enddo
      factorial = f

 end function factorial

!***********************************************************************
!BOP
! !IROUTINE: geoheatflux_accumulate_tavg_once
! !INTERFACE:

   subroutine geoheatflux_accumulate_tavg_once(MASK_TRBUDGET)

! !DESCRIPTION:
!  There is no good place to accumulate one-time ("once") quantities
!   in the tavg "once" output stream. The timing of the various time
!   flags is complicated and confusing; the variables cannot be
!   accumulated too early nor too late. This subroutine is intended
!   to be a place-holder for those variables that are typically
!   added to the "once" tavg stream.
!
! !REVISION HISTORY:
!  added 14 November 2016 njn01

 !INPUT PARAMETERS:

   real(r8),dimension(nx_block,ny_block,km,nblocks_clinic), intent(in) :: MASK_TRBUDGET

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  local or common variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: iblock
   integer (int_kind) :: k
  
   real (r8) :: geo_volume

   real (r8), allocatable :: WORK(:,:,:)

!-----------------------------------------------------------------------
!
!  miscellaneous one-time output variables
!
!-----------------------------------------------------------------------

   allocate (WORK(nx_block,ny_block,nblocks_clinic))
   WORK = c0

   if (registry_match('geoheatflux_choice_spatial')) then

     !*** compute total geo volume
     geo_volume = c0
     do k=1,km
       geo_volume = geo_volume + global_sum(MASK_TRBUDGET_GEO(:,:,k,:)*TAREA(:,:,1:nblocks_clinic)*dz(k), &
                                 distrb_clinic, field_loc_center)
     enddo

     do iblock = 1,nblocks_clinic
       call accumulate_tavg_field(GEOHFLUX(:,:,iblock)/hflux_factor,tavg_GEOHFLUX,iblock,1)

       !*** compute vertical integratal of GEOHFLUX 
       do k=1,km
         WORK(:,:,iblock) = WORK(:,:,iblock) +   &
                            GEOHFLUX(:,:,iblock)*MASK_TRBUDGET_GEO(:,:,k,iblock)*dz(k)
       enddo ! k
     enddo ! iblock

     WORK = WORK*TAREA(:,:,1:nblocks_clinic)

     !*** divide by total geo volume and convert to MKS
     geo_mean = global_sum(WORK,distrb_clinic,field_loc_center)
     geo_mean = geo_mean/geo_volume/hflux_factor

   else
     !### DEBUG fix this!
!########## DEBUG
     !*** compute total geo volume
     geo_volume = c0
     do k=1,km
       geo_volume = geo_volume + global_sum(MASK_TRBUDGET_GEO(:,:,k,:)*TAREA(:,:,1:nblocks_clinic)*dz(k), &
                                 distrb_clinic, field_loc_center)
     enddo
     call document('geoheatflux_accumulate_tavg_once', ' geo_volume', geo_volume)
     call document('geoheatflux_accumulate_tavg_once', ' volume_t  ' ,volume_t)
!########## DEBUG
     geo_mean = geoheatflux_const/hflux_factor
   endif

   call document('geoheatflux_accumulate_tavg_once', ' geo_mean (mks)' ,geo_mean)

   deallocate (WORK)

   end subroutine geoheatflux_accumulate_tavg_once
!***********************************************************************

 end module geoheatflux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
