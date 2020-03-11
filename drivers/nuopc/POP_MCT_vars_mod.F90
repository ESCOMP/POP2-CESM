module POP_MCT_vars_mod

  use mct_mod
  use kinds_mod
  use blocks      , only : get_block, block
  use domain      , only : nblocks_clinic, blocks_clinic
  use communicate , only : my_task, master_task
  use grid        , only : TLOND, TLATD, TAREA, KMT
  use domain_size , only : km, nx_global, ny_global
  use io_types    , only : stdout
  use constants   , only : radius

  implicit none
  private

  public  :: pop_mct_init
  private :: set_gsmap
  private :: set_domain

  integer(int_kind) , public :: POP_MCT_OCNID
  type(mct_gsMap)   , public :: POP_MCT_gsMap_o     ! 2d, points to cdata
  type(mct_gGrid)   , public :: POP_MCT_dom_o       ! 2d, points to cdata
  type(mct_gsMap)   , public :: POP_MCT_gsMap3d_o   ! for 3d streams, local
  type(mct_gGrid)   , public :: POP_MCT_dom3d_o     ! for 3d streams, local

!=======================================================================
contains
!=======================================================================

  subroutine pop_mct_init(ocnid_in, mpicom_in)

    ! input/output variables
    integer, intent(in) :: ocnid_in
    integer, intent(in) :: mpicom_in

    ! local variables
    integer :: lsize
    !-----------------------------------------------------------------------

    POP_MCT_OCNID = ocnid_in

    call set_gsmap( mpicom_in, POP_MCT_OCNID, POP_MCT_gsmap_o, POP_MCT_gsmap3d_o)
    lsize = mct_gsMap_lsize(POP_MCT_gsmap_o, mpicom_in)

    call set_domain( lsize   , POP_MCT_gsmap_o  , POP_MCT_dom_o)
    call set_domain( lsize*km, POP_MCT_gsmap3d_o, POP_MCT_dom3d_o)

  end subroutine pop_mct_init

  !=======================================================================

  subroutine set_gsmap( mpicom_ocn, OCNID, gsMap_ocn, gsMap3d_ocn )

    ! Set the module level gsmap variables pop decomposition

    ! input/output variables
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn
    type(mct_gsMap), intent(inout) :: gsMap3d_ocn

    !  local variables
    integer,allocatable :: gindex(:)
    integer (int_kind)  :: i,j, k, n, iblock
    integer (int_kind)  :: lsize, gsize
    integer (int_kind)  :: ier
    type (block)        :: this_block ! block information for current block
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Build the POP grid numbering for MCT
    !  NOTE:  Numbering scheme is: West to East and South to North starting
    !  at the south pole.  Should be the same as what's used in SCRIP
    !-----------------------------------------------------------------------

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n=n+1
          enddo
       enddo
    enddo
    lsize = n

    !--- 2d ---
    gsize = nx_global*ny_global
    allocate(gindex(lsize),stat=ier)

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n=n+1
             gindex(n) = (this_block%j_glob(j)-1)*(nx_global) + this_block%i_glob(i) 
          enddo
       enddo
    enddo

    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )
    deallocate(gindex)

    !--- 3d ---
    gsize = nx_global*ny_global*km
    lsize = lsize*km
    allocate(gindex(lsize),stat=ier)

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do k = 1,km
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n=n+1
                gindex(n) = (k-1)*(nx_global*ny_global) + &
                     (this_block%j_glob(j)-1)*(nx_global) + this_block%i_glob(i) 
             enddo
          enddo
       enddo
    enddo

    call mct_gsMap_init( gsMap3d_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )
    deallocate(gindex)

  end subroutine set_gsmap

  !=======================================================================

  subroutine set_domain( lsize, gsMap_o, dom_o)

    !-------------------------------------------------------------------
    ! create mct domain type, lat/lon in degrees, area in radians^2, 
    ! mask is 1 (ocean), 0 (non-ocean)
    !-------------------------------------------------------------------

    ! input/output parameters
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o     

    !  local variables
    integer, pointer            :: idata(:)
    real(r8), pointer           :: data(:)
    integer (int_kind)          :: i,j, k, n, iblock, klev, ier
    type (block)                :: this_block ! block information for current block
    character(len=*), parameter :: subname = "ocn_domain_mct"
    !-------------------------------------------------------------------

    ! initialize mct domain type

    call mct_gGrid_init( GGrid=dom_o, &
         CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsize )

    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))

    !-------------------------------------------------------------------
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_o, my_task, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)

    !-------------------------------------------------------------------
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !-------------------------------------------------------------------

    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize) 

    !-------------------------------------------------------------------
    ! Fill in correct values for domain components
    !-------------------------------------------------------------------

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n=n+1
          enddo
       enddo
    enddo

    ! check that the lsize is an exact multiple of 2d size
    if (n > 0) then
       if (mod(lsize,n) /= 0) then
          write(stdout,*) trim(subname),' ERROR: 2d/3d size ',lsize,n
          call shr_sys_abort( SubName//":: 2d/3d size mismatch")
       else 
          klev = lsize/n
       endif
    else
       klev = 0
    endif

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do k = 1,klev
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n=n+1
                data(n) = TLOND(i,j,iblock)
             enddo
          enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize) 

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do k = 1,klev
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n=n+1
                data(n) = TLATD(i,j,iblock)
             enddo
          enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize) 

    n=0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do k = 1,klev
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n=n+1
                data(n) = TAREA(i,j,iblock)/(radius*radius)
             enddo
          enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"area",data,lsize) 

    n=0 
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do k = 1,klev
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n=n+1
                data(n) = float(KMT(i,j,iblock)) 
                if (data(n) > 1.0_r8) data(n) = 1.0_r8
             enddo
          enddo
       enddo
    enddo
    call mct_gGrid_importRattr(dom_o,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine set_domain

end module POP_MCT_vars_mod
