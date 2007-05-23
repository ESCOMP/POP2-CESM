!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module diag_bsf 

!BOP
! !MODULE: diag_bsf
!
! !DESCRIPTION:
!  This module contains code to compute the barotropic stream function

! !REVISION HISTORY:
!     SVN:$Id$

! !USES:
   use kinds_mod 
   use broadcast
   use boundary
   use domain_size
   use domain
   use constants
   use exit_mod
   use global_reductions
   use gather_scatter
   use grid
   use io
   use shr_sys_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:
   public :: init_diag_bsf,         &
             pcg_diag_bsf_solver

!EOP
!BOC
   integer (int_kind), parameter :: &
      mxisle = 200  ! max number of islands

   real (r8) ::                     &
      r_norm                         ! normalization constant

   integer (int_kind) ::            &
      nisle                          ! number of islands

   integer (int_kind), dimension(mxisle) ::  &
      npts_is                        ! number of island boundary points

   integer (int_kind), dimension(nx_block,ny_block,max_blocks_tropic) ::  &
      ISMASK_B                       ! island mask on barotropic decomposition

!     coefficients for the 9-point operator

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      BSFC, BSFN, BSFS, BSFE, BSFW, BSFNE, BSFSE, BSFNW, BSFSW

   real (r8),dimension(nx_block,ny_block,max_blocks_tropic)  ::   &
      BSF_AN,     &
      BSF_AE,     &
      BSF_ANE,    &
      BSF_A0,     &
      BSF_RCALCT_B    

   real (r8), dimension(mxisle) ::  &
      aislandr                       ! island diagonal coefficient


!EOC
!***********************************************************************

 contains

!***********************************************************************

!BOP
! !IROUTINE: init_diag_bsf
! !INTERFACE:
!
 subroutine init_diag_bsf 

!
! !DESCRIPTION:
!  This subroutine initializes variables associated with the computation of
!  the barotropic stream function

! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     initialize the necessary fields for the diagnostic barotropic
!     streamfunction calculation
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      isle,               &
      iblock,             &
      i,j

   integer (int_kind), dimension(:,:,:), allocatable ::  &
      ISMASK,                        &! island mask on baroclinic decomposition
      WORK1

   logical (log_kind), dimension(:,:,:), allocatable ::  &
      MASKI

   logical (log_kind) ::  &
      maskt


!     coefficients for the 9-point operator

   real (r8), dimension(:,:,:), allocatable :: &
      BSFC, BSFN, BSFS, BSFE, BSFW, BSFNE, BSFSE, BSFNW, BSFSW,  &
      WORK0,WORK2,RCALC_TMP

   real (r8) ::         &
      xne,xse,xnw,xsw,  &! contribution to coefficients from x,y
      yne,yse,ynw,ysw,  &!   components of divergence
      ase,anw,asw


   if ( my_task == master_task ) then
     write (stdout,*) ' '
     write (stdout,*) 'Initializing diagnostic BSF variables ....'
     call shr_sys_flush(stdout)
   endif


   allocate(WORK0    (nx_block,ny_block,nblocks_clinic), &
            WORK2    (nx_block,ny_block,nblocks_clinic), &
            BSFC     (nx_block,ny_block,nblocks_clinic), &
            BSFN     (nx_block,ny_block,nblocks_clinic), &
            BSFNE    (nx_block,ny_block,nblocks_clinic), &
            BSFE     (nx_block,ny_block,nblocks_clinic), &
            BSFSE    (nx_block,ny_block,nblocks_clinic), &
            BSFS     (nx_block,ny_block,nblocks_clinic), &
            BSFSW    (nx_block,ny_block,nblocks_clinic), &
            BSFW     (nx_block,ny_block,nblocks_clinic), &
            BSFNW    (nx_block,ny_block,nblocks_clinic), &
            RCALC_TMP(nx_block,ny_block,nblocks_clinic))
   allocate(ISMASK   (nx_block,ny_block,nblocks_clinic), &
            WORK1    (nx_block,ny_block,nblocks_clinic))
   allocate(MASKI    (nx_block,ny_block,nblocks_clinic))

   r_norm = c0

   ISMASK_B = 0
   WORK1 = 0 
   WORK2 = p25*merge( 1, 0, KMU > 0 )   ! WORK2 = 0.25*(0,1) on (land,ocean)

!-----------------------------------------------------------------------
!
!     9-point coefficients
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic

      WORK0(:,:,iblock) = c0
      BSFC (:,:,iblock) = c0
      BSFN (:,:,iblock) = c0
      BSFNE(:,:,iblock) = c0
      BSFE (:,:,iblock) = c0
      BSFSE(:,:,iblock) = c0
      BSFS (:,:,iblock) = c0
      BSFSW(:,:,iblock) = c0
      BSFW (:,:,iblock) = c0
      BSFNW(:,:,iblock) = c0
      RCALC_TMP(:,:,iblock) = c0

      do j=2,ny_block
      do i=2,nx_block

       xne = WORK2(i  ,j  ,iblock)*DXUR(i  ,j  ,iblock)*DYU(i  ,j  ,iblock)
       xse = WORK2(i  ,j-1,iblock)*DXUR(i  ,j-1,iblock)*DYU(i  ,j-1,iblock)
       xnw = WORK2(i-1,j  ,iblock)*DXUR(i-1,j  ,iblock)*DYU(i-1,j  ,iblock)
       xsw = WORK2(i-1,j-1,iblock)*DXUR(i-1,j-1,iblock)*DYU(i-1,j-1,iblock)

       yne = WORK2(i  ,j  ,iblock)*DYUR(i  ,j  ,iblock)*DXU(i  ,j  ,iblock)
       yse = WORK2(i  ,j-1,iblock)*DYUR(i  ,j-1,iblock)*DXU(i  ,j-1,iblock)
       ynw = WORK2(i-1,j  ,iblock)*DYUR(i-1,j  ,iblock)*DXU(i-1,j  ,iblock)
       ysw = WORK2(i-1,j-1,iblock)*DYUR(i-1,j-1,iblock)*DXU(i-1,j-1,iblock)

       BSFNE(i,j,iblock) = xne + yne
       BSFSE(i,j,iblock) = xse + yse
       BSFNW(i,j,iblock) = xnw + ynw
       BSFSW(i,j,iblock) = xsw + ysw

       BSFE(i,j,iblock)  = xne + xse - yne - yse
       BSFW(i,j,iblock)  = xnw + xsw - ynw - ysw
       BSFN(i,j,iblock)  = yne + ynw - xne - xnw
       BSFS(i,j,iblock)  = yse + ysw - xse - xsw

       !*** diagnoal coefficient
       BSFC(i,j,iblock)  = -(BSFNE(i,j,iblock) + BSFSE(i,j,iblock)   &
                           + BSFNW(i,j,iblock) + BSFSW(i,j,iblock))   

       WORK0    (i,j,iblock) = TAREA(i,j,iblock)**2

      enddo ! i
      enddo ! j
   enddo  ! iblock

   call update_ghost_cells (BSFN, bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFNE,bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFE, bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFSE,bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFS, bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFSW,bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFW, bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call update_ghost_cells (BSFNW,bndy_clinic,field_loc_center,  &
                            field_type_scalar)
   call redistribute_blocks(BSF_AN,  distrb_tropic, BSFN,  distrb_clinic)
   call redistribute_blocks(BSF_AE,  distrb_tropic, BSFE,  distrb_clinic)
   call redistribute_blocks(BSF_ANE, distrb_tropic, BSFNE, distrb_clinic)
   call redistribute_blocks(BSF_A0,  distrb_tropic, BSFC , distrb_clinic)

   r_norm = c1/global_sum(WORK0,distrb_clinic,field_loc_center,RCALCT)

   call redistribute_blocks(BSF_RCALCT_B, distrb_tropic, RCALCT, distrb_clinic)

!-----------------------------------------------------------------------
!
!  calculate island mask and redistribute ISMASK after it has been computed
!
!-----------------------------------------------------------------------

   call island_mask_diag_bsf (ISMASK)

   call redistribute_blocks(ISMASK_B, distrb_tropic, ISMASK, distrb_clinic)

!-----------------------------------------------------------------------
!  calculate island coefficients
!-----------------------------------------------------------------------

   do isle=1,mxisle
     aislandr(isle) = c0
   enddo

   if ( my_task == master_task ) then
     write (stdout,*) ' '
     write (stdout,*) ' island coefficients: '
     call shr_sys_flush(stdout)
   endif

   do isle=1,nisle

     WORK0 = c0
     MASKI(:,:,:) = (ISMASK(:,:,:)  ==  isle)

     do iblock = 1,nblocks_clinic

      do j=2,ny_block-1
      do i=2,nx_block-1

        !*** north faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i,j+1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = BSFN(i,j,iblock)

        !*** northeast faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i+1,j+1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFNE(i,j,iblock) 

        !*** northwest faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i-1,j+1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFNW(i,j,iblock)

        !*** south faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i,j-1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFS(i,j,iblock)

        !*** southeast faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i+1,j-1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFSE(i,j,iblock)

        !*** southwest faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i-1,j-1,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFSW(i,j,iblock)

        !*** east faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i+1,j,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFE(i,j,iblock)

        !*** west faces
        maskt = MASKI(i,j,iblock) .and. (ISMASK(i-1,j,iblock) == 0)
        if (maskt) WORK0(i,j,iblock) = WORK0(i,j,iblock) + BSFW(i,j,iblock)

      end do  ! i
      end do  ! j
   end do  ! iblock

   !*** diagonal coefficient on island
   aislandr(isle) = -c1/global_sum(WORK0, distrb_clinic, field_loc_center)

   if ( my_task == master_task ) then 
     write (stdout,*) ' island ', isle, ' coefficient = ', c1/aislandr(isle)
     call shr_sys_flush(stdout)
   endif

 enddo ! isle

 deallocate(WORK0,WORK2,BSFC,BSFN,BSFNE,BSFE,BSFSE,BSFS,BSFSW,  &
            BSFW,BSFNW,RCALC_TMP)
 deallocate(ISMASK,WORK1)
 deallocate(MASKI)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_diag_bsf 

!***********************************************************************

!BOP
! !IROUTINE: pcg_diag_bsf
! !INTERFACE:

 subroutine pcg_diag_bsf ( X, B )

! !DESCRIPTION:
!
!  Conjugate-gradient solver for diagnosed barotropic streamfunction 
!  in free-surface formulation.
!  Based upon pcg in module solvers.F90

! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8),dimension(nx_block,ny_block,max_blocks_tropic),intent(in) ::   &
     B

! !INPUT/OUTPUT PARAMETERS:

   real (r8),dimension(nx_block,ny_block,max_blocks_tropic),intent(inout) ::  &
     X

!EOP
!BOC

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: mxscan =1000, ncheck = 100 

   real (r8) ::  &
      err = 1.0e-6_r8

   integer (int_kind) ::  &
      m, isle, mscan, iblock

   real (r8) ::  &
      eta0, eta1, eta2, rr, risle, rms_residual

   logical (log_kind), dimension(nx_block,ny_block,max_blocks_tropic) ::  &
      CALCZ

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) ::  &
      R, S, Q, WORK1, WORK2, WORK3, BSFCR

   type (block) ::      &
      this_block         ! block information for current block

!-----------------------------------------------------------------------
!
!     initialization 
!
!-----------------------------------------------------------------------

!######### debug? ##########
      call update_ghost_cells (X,bndy_tropic,field_loc_center,  &
                               field_type_scalar)
!##### end debug? ##########
      R     = c0  
      Q     = c0
      WORK1 = c0
      WORK2 = c0

      CALCZ = (ISMASK_B == 0)  ! mask for interior ocean points

!-----------------------------------------------------------------------
!
!     diagonal preconditioner
!
!-----------------------------------------------------------------------

      where (BSF_A0 /= c0) 
        BSFCR = c1/BSF_A0
      elsewhere
        BSFCR = c0
      end where

      where ( .not. CALCZ )  BSFCR = c0

      do isle=1,nisle
        where ( abs(ISMASK_B) == isle )  BSFCR = aislandr(isle)
      enddo

!-----------------------------------------------------------------------
!
!     compute initial residual
!
!-----------------------------------------------------------------------

!########### debug? ##########
      call update_ghost_cells (X,bndy_tropic,field_loc_center,  &
                               field_type_scalar)
!########### debug? ##########

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)
   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)

      call btrop_operator_diag(WORK1,X,this_block,iblock)
      WORK1(:,:,iblock) = B(:,:,iblock) - WORK1(:,:,iblock)

      where (CALCZ(:,:,iblock))
        R    (:,:,iblock) = WORK1(:,:,iblock)
        WORK2(:,:,iblock) = WORK1(:,:,iblock)
      endwhere
   end do ! block loop
   !$OMP END PARALLEL DO

   do isle=1,nisle                           ! island loop
     WORK3 = merge(c1, c0, ISMASK_B == isle)
     risle = global_sum( WORK1*WORK3,distrb_tropic,field_loc_center)
     where ( abs(ISMASK_B) == isle )  R = risle
     where (     ISMASK_B  == isle )  WORK2 = risle/npts_is(isle)

   enddo

   rms_residual = sqrt(global_sum(WORK2*R,distrb_tropic,  &
                                  field_loc_center, BSF_RCALCT_B)*r_norm)
   if ( my_task == master_task ) then
     write (stdout,*) ' '
     write (stdout,*) ' convergence info from pcg_diag_bsf: '
     write (stdout,*) ' iter = ', 0, ' rms_resid = ', rms_residual
     call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!
!     initialize fields and scalars
!
!-----------------------------------------------------------------------

      call update_ghost_cells (R,bndy_tropic,field_loc_center,  &
                               field_type_scalar)
      S     = c0
      eta0  = c1
      mscan = mxscan

!-----------------------------------------------------------------------
!
!     iterate
!
!-----------------------------------------------------------------------

   iter_loop: do m=1,mxscan

     !$OMP PARALLEL DO PRIVATE(iblock,this_block)
     do iblock=1,nblocks_tropic
        this_block = get_block(blocks_tropic(iblock),iblock)

        !*** use diagnoal preconditioner
        WORK1(:,:,iblock) = R(:,:,iblock) * BSFCR (:,:,iblock)  

!-----------------------------------------------------------------------
!
!       update conjugate direction vector s
!
!-----------------------------------------------------------------------

        WORK2(:,:,iblock) = R(:,:,iblock)

        do isle=1,nisle   ! island loop
          where ( ISMASK_B(:,:,iblock)  == isle )  &
                  WORK2(:,:,iblock) = R(:,:,iblock)/npts_is(isle)
        enddo

     end do ! block loop
     !$OMP END PARALLEL DO

     !*** (r,(PC)r)
     eta1 = global_sum( WORK2*WORK1, distrb_tropic,  &
                        field_loc_center,BSF_RCALCT_B ) 
     eta2 = eta1/eta0

     call update_ghost_cells (WORK1,bndy_tropic,field_loc_center,  &
                              field_type_scalar)
     call update_ghost_cells (S    ,bndy_tropic,field_loc_center,  &
                              field_type_scalar)


     !$OMP PARALLEL DO PRIVATE(iblock,this_block)
     do iblock=1,nblocks_tropic
        this_block = get_block(blocks_tropic(iblock),iblock)

        S(:,:,iblock) =  WORK1(:,:,iblock) + S(:,:,iblock)*eta2

!-----------------------------------------------------------------------
!
!       compute As
!
!-----------------------------------------------------------------------

        call btrop_operator_diag ( WORK1,S,this_block, iblock)
     end do ! block loop
     !$OMP END PARALLEL DO


     where (CALCZ)
       Q     = WORK1
       WORK2 = WORK1
     endwhere

     do isle=1,nisle   ! island loop
       WORK3 = merge(c1, c0, ISMASK_B == isle)
       risle = global_sum( WORK1, distrb_tropic, field_loc_center,WORK3 )
       where ( abs(ISMASK_B) == isle ) Q     = risle
       where (     ISMASK_B  == isle ) WORK2 = risle/npts_is(isle)
     enddo

!-----------------------------------------------------------------------
!
!       compute next solution and residual
!
!-----------------------------------------------------------------------

        call update_ghost_cells (Q,bndy_tropic, field_loc_center,  &
                                 field_type_scalar)
        call update_ghost_cells (S,bndy_tropic, field_loc_center,  &
                                 field_type_scalar)
        eta0 = eta1
        eta1 = eta0/global_sum( WORK2*S, distrb_tropic,  &
                                field_loc_center, BSF_RCALCT_B )

     !$OMP PARALLEL DO PRIVATE(iblock,this_block)
     do iblock=1,nblocks_tropic
        this_block = get_block(blocks_tropic(iblock),iblock)

        X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)
        R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)

        if (mod(m,ncheck) == 0) then
          call btrop_operator_diag ( WORK1, X, this_block, iblock)

          WORK1(:,:,iblock) = B(:,:,iblock) - WORK1(:,:,iblock)
          where (CALCZ(:,:,iblock))
            R(:,:,iblock)     = WORK1(:,:,iblock)
            WORK2(:,:,iblock) = WORK1(:,:,iblock)
          endwhere
        endif ! mod(m,ncheck)

     end do ! block loop
     !$OMP END PARALLEL DO

     if (mod(m,ncheck) == 0) then
       do isle=1,nisle   ! island loop
         WORK3 = merge(c1, c0, ISMASK_B == isle)
         risle = global_sum( WORK1, distrb_tropic,  &
                             field_loc_center, WORK3 )
         where ( abs(ISMASK_B) == isle ) R     = risle
         where (     ISMASK_B  == isle ) WORK2 = risle/npts_is(isle)
       enddo

!-----------------------------------------------------------------------
!
!       test for convergence
!
!-----------------------------------------------------------------------
          rr = global_sum( WORK2*R, distrb_tropic,  &
                           field_loc_center,BSF_RCALCT_B )   ! (r,r)

          rms_residual = sqrt(rr*r_norm)
 
          if ( my_task == master_task ) then
            write (stdout,*) ' iter = ', m, ' rms_resid = ', rms_residual 
          endif

          if (rms_residual < err) then
            mscan = m
            exit iter_loop
          endif

        endif ! mod(m,ncheck)

      enddo iter_loop

      if ( mscan == mxscan ) then
        if ( my_task == master_task ) then
          write (stdout,*)                                    &
            ' WARNING (pcg_diag_bsf): No convergence after ', &
                                      mxscan, ' scans'
        endif
      endif

!EOC
!-----------------------------------------------------------------------

 end subroutine pcg_diag_bsf

!***********************************************************************

!BOP
! !IROUTINE: island_mask_diag_bsf
! !INTERFACE:

 subroutine island_mask_diag_bsf (ISMASK)

!
! !DESCRIPTION:
!  This subroutine finds island mask for the barotropic streamfunction
!  diagnostic computation

! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:
   integer (int_kind), dimension(:,:,:), intent(inout) ::  &
      ISMASK                         ! island mask on baroclinic decomposition
  
!EOP
!BOC
  
!-----------------------------------------------------------------------
!
!     output variables: nisle, npts_is, ISMASK
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i, j, n, imax, error_code
   integer (int_kind), dimension(:,:), allocatable ::  &
     KALL, KINT, I1, I2, ISMASK_G, WORK1, WORK2

   logical (log_kind), dimension(:,:), allocatable ::  &
     LAND, LANDLEFT 

   allocate(KALL     (nx_global,ny_global), &
            KINT     (nx_global,ny_global), &
            I1       (nx_global,ny_global), &
            I2       (nx_global,ny_global), &
            ISMASK_G (nx_global,ny_global), &
            WORK1    (nx_global,ny_global), &
            WORK2    (nx_global,ny_global)  )

   allocate(LAND     (nx_global,ny_global), &
            LANDLEFT (nx_global,ny_global)  )

   KALL = 0; KINT = 0;  I1 = 0;  I2 = 0;  ISMASK_G = 0;  WORK1 = 0;  WORK2 = 0

   call gather_global (WORK1, KMT, master_task,distrb_clinic)

   error_code = 0

   if ( my_task == master_task ) then

     LAND = .false.
     where ( WORK1 == 0 )  LAND = .true.

     do j=1,ny_global
       do i=1,nx_global
         I1(i,j) = i
         I2(i,j) = j
       enddo
     enddo

     where ( .not. LAND )
       KINT = nx_global * ny_global + 1
     elsewhere
       KINT = I1 + (I2-1) * nx_global
     endwhere

     KALL = KINT

     imax = nx_global + 2 * ny_global

     do i=1,imax

       WORK1 = cshift(KALL, +1, 2)
       WORK1(:,ny_global) = KALL(:,ny_global)
       WORK2 = cshift(KALL, -1, 2)
       WORK2(:,1) = KALL(:,1)
       KALL = min(WORK1, WORK2, KALL)

       WORK1 = cshift(KALL, +1, 2)
       WORK1(:,ny_global) = KALL(:,ny_global)
       WORK2 = cshift(KALL, -1, 2)
       WORK2(:,1) = KALL(:,1)
       KALL = min(WORK1, WORK2, KALL)

       WORK1 = cshift(KALL, +1, 1)
       WORK2 = cshift(KALL, -1, 1)
       KALL  = min(WORK1, WORK2, KALL)

       WORK1 = cshift(KALL, +1, 1)
       WORK2 = cshift(KALL, -1, 1)
       KALL  = min(WORK1, WORK2, KALL)

       where (LAND)  KINT = KALL
       KALL = KINT

     enddo

     ISMASK_G = 0
     LANDLEFT = LAND

     where ( KINT == KINT(1,ny_global) )  ! north continent
       ISMASK_G = 999                      ! NOTE: nisle must be less than 999
       LANDLEFT = .false.
     endwhere

     n = 0
     do i=1,nx_global
       do j=1,ny_global
         if (LANDLEFT(i,j)) then
           n = n + 1 
           where ( KINT == KINT(i,j) )
             ISMASK_G = n
             LANDLEFT = .false.
           endwhere
         endif
       enddo
     enddo

     nisle = n

     write (stdout,*) ' '
     write (stdout,*) ' number of islands = ', nisle

     if ( nisle > mxisle )  error_code = 1 

   endif        ! if master_task

   call broadcast_scalar (error_code, master_task)

   if ( error_code /= 0) then
   call exit_POP(sigAbort, &
                 'ERROR (island_mask_diag_bsf): nisle > mxisle')
   endif

   if ( my_task == master_task ) then

     WORK1 = cshift(ISMASK_G, +1, 1)
     WORK2 = cshift(ISMASK_G, -1, 1)
     KALL  = max(ISMASK_G, WORK1, WORK2)

     WORK1 = cshift(ISMASK_G, +1, 2)
     WORK2 = cshift(ISMASK_G, -1, 2)
     KALL  = max(KALL, WORK1, WORK2)

     WORK1 = cshift(ISMASK_G, +1, 2)
     WORK2 = cshift(WORK1,    +1, 1)
     KALL  = max(KALL, WORK2)

     WORK1 = cshift(ISMASK_G, +1, 2)
     WORK2 = cshift(WORK1,    -1, 1)
     KALL  = max(KALL, WORK2)

     WORK1 = cshift(ISMASK_G, -1, 2)  
     WORK2 = cshift(WORK1,    +1, 1)
     KALL  = max(KALL, WORK2)
 
     WORK1 = cshift(ISMASK_G, -1, 2)
     WORK2 = cshift(WORK1,    -1, 1)
     KALL  = max(KALL, WORK2)
  
     where (.not.LAND)
       ISMASK_G =  -KALL
     endwhere
     ISMASK_G = -ISMASK_G

     write (stdout,*) ' '
     do i=1,nisle
       npts_is(i) = count(ISMASK_G == i)
       write (stdout,*) ' island ', i,' # points on boundary = ', npts_is(i)  
     enddo

   endif        ! if master_task

   call scatter_global(ISMASK, ISMASK_G, master_task, distrb_clinic, &
        field_loc_NEcorner, field_type_scalar)

   call update_ghost_cells(ISMASK, bndy_clinic, &
                           field_loc_center, field_type_scalar)


   call broadcast_scalar ( nisle,   master_task )
   call broadcast_array  ( npts_is, master_task )

   deallocate(KALL,KINT,I1,I2,ISMASK_G,WORK1,WORK2)
   deallocate(LAND,LANDLEFT)
!-----------------------------------------------------------------------
!EOC

 end subroutine island_mask_diag_bsf


!***********************************************************************

!BOP
! !IROUTINE: btrop_operator_diag
! !INTERFACE:

 subroutine btrop_operator_diag(AX,X,this_block,bid)

! !DESCRIPTION:
!  This routine applies the nine-point stencil operator for the
!  barotropic solver, using weights local to this diagnostics
!  module.  It takes advantage of some 9pt weights being 
!  shifted versions of others.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: & 
      X                  ! array to be operated on 

   type (block), intent(in) :: &
      this_block             ! block info for this block

   integer (int_kind), intent(in) :: &
      bid                    ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX                     ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j                ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = BSF_A0 (i  ,j  ,bid)*X(i  ,j  ,bid) + &
                    BSF_AN (i  ,j  ,bid)*X(i  ,j+1,bid) + &
                    BSF_AN (i  ,j-1,bid)*X(i  ,j-1,bid) + &
                    BSF_AE (i  ,j  ,bid)*X(i+1,j  ,bid) + &
                    BSF_AE (i-1,j  ,bid)*X(i-1,j  ,bid) + &
                    BSF_ANE(i  ,j  ,bid)*X(i+1,j+1,bid) + &
                    BSF_ANE(i  ,j-1,bid)*X(i+1,j-1,bid) + &
                    BSF_ANE(i-1,j  ,bid)*X(i-1,j+1,bid) + &
                    BSF_ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator_diag

!BOP
! !IROUTINE: pcg_diag_bsf_solver
! !INTERFACE:

 subroutine pcg_diag_bsf_solver(SOLN, RHS)

! !DESCRIPTION:
!  Solves the elliptic equation for bsf computation.
!  First redistributes arrays to the barotropic block distribution for
!  better performance of the solver.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(in) :: &
      RHS                  ! right-hand-side of linear system
                           !  for blocks in baroclinic distribution

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
     intent(inout) :: &
     SOLN              ! on input,  initial guess
                       ! on output, final solution for sfc pressure

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      S_TROPIC,         &! surface pressure in barotropic distribution
      RHS_TROPIC         ! right hand side  in barotropic distribution

!-----------------------------------------------------------------------
!
!  switch to the barotropic distribution for iterative solvers
!
!-----------------------------------------------------------------------

   call redistribute_blocks(S_TROPIC,   distrb_tropic, &
                            SOLN,       distrb_clinic)
   call redistribute_blocks(RHS_TROPIC, distrb_tropic, &
                            RHS,        distrb_clinic)

!-----------------------------------------------------------------------
!
!  call bsf diagnostic solver
!
!-----------------------------------------------------------------------

   call pcg_diag_bsf(S_TROPIC,RHS_TROPIC)      ! precond conjg grad solver

!-----------------------------------------------------------------------
!
!  switch solution back to the baroclinic distribution
!
!-----------------------------------------------------------------------

   call redistribute_blocks(SOLN,     distrb_clinic, &
                            S_TROPIC, distrb_tropic)

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg_diag_bsf_solver
!***********************************************************************

 end module diag_bsf

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
