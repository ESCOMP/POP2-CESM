!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module passive_tracers

!BOP
! !MODULE: passive_tracers
! !DESCRIPTION:
!  This module provides support for passive tracers.
!  The base model calls subroutines in this module which then call
!     subroutines in individual passive tracer modules.

! !USES:

   use kinds_mod, only: r8, int_kind, log_kind, char_len
   use blocks, only: block, nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km, nt
   use communicate, only: my_task, master_task
   use broadcast, only: broadcast_scalar
   use prognostic, only: TRACER, tracer_d, curtime
   use forcing_shf, only: SHF_QSW
   use io_types, only: stdout, nml_in, nml_filename,io_field_desc
   use exit_mod, only: sigAbort, exit_pop
   use shr_sys_mod
   use timers, only : timer_start, timer_stop
   use tavg, only: define_tavg_field
   use constants, only : delim_fmt, char_blank, salt_to_ppt, ocn_ref_salinity, &
      ppt_to_salt, c0
   use registry, only: register_string, registry_match

   use cfc11_mod, only :           &
      cfc11_tracer_cnt,            &
      cfc11_ind_begin,             &
      cfc11_ind_end,               &
      cfc11_tracer_names,          &
      cfc11_init,                  &
      cfc11_init_sflux,            &
      cfc11_set_sflux,             &
      cfc11_ind2name,              &
      cfc11_name2ind,              &
      cfc11_sflux_timer,           &
      cfc11_tavg,                  &
      cfc11_tavg_forcing,          &
      cfc11_tracer_field_info


   use iage_mod, only :           &
      iage_tracer_cnt,            &
      iage_ind_begin,             &
      iage_ind_end,               &
      iage_tracer_names,          &
      iage_init,                  &
      iage_set_interior,          &
      iage_reset,                 &
      iage_ind2name,              &
      iage_name2ind,              &
      iage_tavg


   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public ::                                 &
      init_passive_tracers,                  &
      set_passive_tracers_interior,          &
      set_passive_tracers_sflux,             &
      init_passive_tracers_sflux,            &
      init_passive_tracers_interior_restore, &
      reset_passive_tracers,                 &
      write_restart_passive_tracers,         &
      tavg_passive_tracers,                  &
      tavg_passive_tracers_sflux,            &
      passive_tracer_ind2name,               &
      tracer_ref_val

!-----------------------------------------------------------------------
!  variables for automatically generated tavg passive-tracer fields
!-----------------------------------------------------------------------

   integer (int_kind), parameter, public :: & 
      num_auto_gen_tr = 11,         &!number of auto-generated passive-tracer fields
      nt_passive      = nt-2         !number of passive tracers

   integer (int_kind), dimension (num_auto_gen_tr*nt_passive), public ::  &
      tavg_PASSIVE
!-----------------------------------------------------------------------
!  logical variables that denote if a passive tracer module is on
!-----------------------------------------------------------------------

   logical (kind=log_kind) ::  &
      ecosys_on, cfc11_on, mchl_on, iage_on, tracegas_on

   namelist /passive_tracers_on_nml/  &
      ecosys_on, cfc11_on, mchl_on, iage_on, tracegas_on

!EOP
!BOC

!EOC
!***********************************************************************

 contains

!***********************************************************************

 subroutine init_passive_tracers


   integer (int_kind) :: cumulative_nt, n

   character (char_len) :: tracer_name

   !-----------------------------------------------------------------
   !   register init_passive_tracers
   !-----------------------------------------------------------------
 
   call register_string('init_passive_tracers')


   ecosys_on    = .false.
   cfc11_on     = .false.
   mchl_on      = .false.
   iage_on      = .false.
   tracegas_on  = .false.


   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old')
   10    continue  !*** keep reading until find right namelist
      read(nml_in, nml=passive_tracers_on_nml, err=10)
      close(nml_in)
   end if
   
   if (my_task == master_task) then
     write(stdout,*) ' '
     write(stdout,*) ' Document Namelist Parameters:'
     write(stdout,*) ' ============================ '
     write(stdout,*) ' '
     write(stdout, passive_tracers_on_nml)  
     write(stdout,*) ' '
     call shr_sys_flush (stdout)
   endif

   call broadcast_scalar(ecosys_on, master_task)
   call broadcast_scalar(cfc11_on,  master_task)
   call broadcast_scalar(mchl_on,  master_task)
   call broadcast_scalar(iage_on,  master_task)
   call broadcast_scalar(tracegas_on,  master_task)



   !-----------------------------------------------------------------
   !  check for modules that require the flux coupler
   !-----------------------------------------------------------------

   if (cfc11_on .and. .not. registry_match('lcoupled')) then
     call exit_POP(sigAbort,'cfc module requires the flux coupler')
   end if



   !-----------------------------------------------------------------
   !  set up indices for passive tracer modules that are on
   !-----------------------------------------------------------------

   cumulative_nt = 2
      

   if (cfc11_on) then
     call set_tracer_indices('CFC11', cfc11_tracer_cnt, cumulative_nt,  &
                             cfc11_ind_begin, cfc11_ind_end)
   end if


   if (iage_on) then
     call set_tracer_indices('IAGE', iage_tracer_cnt, cumulative_nt,  &
                             iage_ind_begin, iage_ind_end)
   end if


   if (cumulative_nt /= nt) then
     
     call exit_POP(sigAbort, &
      'ERROR in init_passive_tracers: declared nt does not match cumulative nt')
   end if

      !-----------------------------------------------------------------
      !  print out tracer names from tracer modules that are on
      !-----------------------------------------------------------------

   if (my_task == master_task) then
     write(stdout,delim_fmt)
     write(stdout,*) 'TRACER INDEX    TRACER NAME'
     write(stdout,1010) 1, 'TEMP'
     write(stdout,1010) 2, 'SALT'
     call shr_sys_flush (stdout)
     do n = 3, nt
        call passive_tracer_ind2name(n, tracer_name)
        write(stdout,1010) n, TRIM(tracer_name)
        call shr_sys_flush (stdout)
     enddo
     write(stdout,delim_fmt)
     call shr_sys_flush (stdout)
   end if


      !-----------------------------------------------------------------
      !  CFC11 block 
      !-----------------------------------------------------------------

   if (cfc11_on) then
      if (my_task == master_task) then
         write(stdout,delim_fmt)
         write(stdout,*) 'cfc11_ind_begin = ', cfc11_ind_begin
         write(stdout,*) 'cfc11_ind_end   = ', cfc11_ind_end
         write(stdout,delim_fmt)
      end if

      call cfc11_init(TRACER(:,:,:,cfc11_ind_begin:cfc11_ind_end,:,:))

   end if


      !-----------------------------------------------------------------
      !  Ideal Age (IAGE) block 
      !-----------------------------------------------------------------

   if (iage_on) then

      if (my_task == master_task) then
         write(stdout,delim_fmt)
         write(stdout,*) 'iage_ind_begin = ', iage_ind_begin
         write(stdout,*) 'iage_ind_end   = ', iage_ind_end
         write(stdout,delim_fmt)
         call shr_sys_flush(stdout)
      end if

      call iage_init(TRACER(:,:,:,iage_ind_begin:iage_ind_end,:,:))

   end if



 1010 format(5X,I2,10X,A)

!-----------------------------------------------------------------------

   end subroutine init_passive_tracers

!***********************************************************************

   subroutine set_passive_tracers_interior(k, this_block, TRACER_SOURCE)

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: k  ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block information for this block

!-----------------------------------------------------------------------
!     input/output variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) ::   &
      TRACER_SOURCE

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      bid                  ! local block address for this block


   bid = this_block%local_id

   if (iage_on) then
      call iage_set_interior(k,                                    &
         TRACER_SOURCE (:,:,iage_ind_begin:iage_ind_end) )
   end if


   end subroutine set_passive_tracers_interior

!***********************************************************************

   subroutine set_passive_tracers_sflux(SMFT,ICE_FRAC,PRESS,STF)


!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,2,max_blocks_clinic), intent(in) :: &
      SMFT   ! 'zonal' & 'merdional' surface velocity fluxes (cm^2/s^2)

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) ::   &
      ICE_FRAC, & ! sea ice fraction (non-dimensional)  
      PRESS       ! sea level atmospheric pressure (Pascals)

!-----------------------------------------------------------------------
!     output variables
!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), intent(inout) :: &
      STF   ! surface fluxes for tracers

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
   real (r8)          :: ref_val
   integer (int_kind) :: n


   if (iage_on) then
!     iage has no surface flux
   end if

!-----------------------------------------------------------------------
!     add virtual fluxes for tracers that specify a non-zero ref_val
!-----------------------------------------------------------------------

      do n=3,nt
        ref_val = tracer_ref_val(n)
        if (ref_val /= c0) STF(:,:,n,:) = STF(:,:,n,:) + &
          (ref_val/(ocn_ref_salinity*ppt_to_salt))*STF(:,:,2,:)
      end do

   end subroutine set_passive_tracers_sflux

!***********************************************************************

   subroutine init_passive_tracers_sflux

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   if (ecosys_on) then
      !call ecosys_init_sflux
   end if

   if (cfc11_on) then
      call cfc11_init_sflux
   end if

   if (mchl_on) then
      !call mchl_init_sflux
   end if

   if (tracegas_on) then
      !call tracegas_init_sflux
   end if

   end subroutine init_passive_tracers_sflux

!***********************************************************************

   subroutine init_passive_tracers_interior_restore

!-----------------------------------------------------------------------
!  assume actual restoring and updating occurs in set_passive_tracers_interior
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   if (ecosys_on) then
      !call ecosys_init_interior_restore
   end if

   end subroutine init_passive_tracers_interior_restore

!***********************************************************************

   subroutine write_restart_passive_tracers(restart_file, action)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   use io_types, only: datafile

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file
   character(*), intent(in) :: action

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   if (ecosys_on) then
      !call ecosys_write_restart(restart_file, action)
   end if

   end subroutine write_restart_passive_tracers

!***********************************************************************

   subroutine reset_passive_tracers(TRACER_NEW, this_block, iblock)

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

   real(r8), dimension(nx_block,ny_block,km,nt), intent(inout) :: &
         TRACER_NEW      ! all tracers at new time for a given block

   integer(int_kind), intent(in) :: &
         iblock          ! block index

   type (block), intent(in) :: &
      this_block   ! block information for this block

!-----------------------------------------------------------------------
!     ecosys does not reset
!     cfc11  does not reset
!     mchl   does not reset
!
!     iage   does reset
!-----------------------------------------------------------------------

   if (iage_on) then
      call iage_reset(  &
         TRACER_NEW(:,:,:,iage_ind_begin:iage_ind_end) )
   end if

   end subroutine reset_passive_tracers


   subroutine tavg_passive_tracers(bid, k)

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: k, bid  ! vertical level index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


   if (iage_on) then
      call iage_tavg(bid, k,   &
         TRACER(:,:,k,iage_ind_begin:iage_ind_end,curtime,bid))
   end if


   end subroutine tavg_passive_tracers
!***********************************************************************


   subroutine tavg_passive_tracers_sflux(STF)

!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------

  REAL(r8), DIMENSION(nx_block,ny_block,nt,max_blocks_clinic), &
      INTENT(IN) :: STF

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!  if (cfc11_on) then
!     call cfc11_tavg_forcing(STF(:,:,cfc11_ind_begin:cfc11_ind_end,:))
!  end if


   end subroutine tavg_passive_tracers_sflux


   subroutine set_tracer_indices(module_string, module_nt,  &
        cumulative_nt, ind_begin, ind_end)
      
!-----------------------------------------------------------------------
!     set the index bounds of a single passive tracer module
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
!     input variables
!-----------------------------------------------------------------------
          
      character (*), intent(in) ::  &
        module_string
          
      integer (kind=int_kind), intent(in) ::  &
        module_nt
          
!-----------------------------------------------------------------------
!     input/output variables
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(inout) ::  &
        cumulative_nt
      
      integer (kind=int_kind), intent(out) ::  &
        ind_begin  &
      , ind_end
      
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
          
      character (char_len) ::  &
        error_string

!-----------------------------------------------------------------------

      ind_begin = cumulative_nt + 1
      ind_end = ind_begin + module_nt - 1
      cumulative_nt = ind_end
      if (my_task == master_task) then
        write(stdout,delim_fmt)
        write(stdout,*) module_string // ' ind_begin = ', ind_begin
        write(stdout,*) module_string // ' ind_end   = ', ind_end
        write(stdout,delim_fmt)
      end if
      if (cumulative_nt > nt) then
         error_string = 'nt too small for module ' // module_string
         call exit_POP(sigAbort, error_string)
      end if


!-----------------------------------------------------------------------

      end subroutine set_tracer_indices

!***********************************************************************

      subroutine passive_tracer_ind2name(ind, name)

!-----------------------------------------------------------------------
!     convert a TRACER index into a tracer name
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------

      integer(kind=int_kind), intent(in)   :: ind
      character(len=char_len), intent(out) :: name

!-----------------------------------------------------------------------

      if (ecosys_on) then
        !if (ind >= ecosys_ind_begin .and. ind <= ecosys_ind_end) then
        !  name = ecosys_tracer_names(ind-(ecosys_ind_begin-1))
        !  return
        !end if
      end if

      if (cfc11_on) then
        if (ind >= cfc11_ind_begin .and. ind <= cfc11_ind_end) then
          name = cfc11_tracer_names(ind-(cfc11_ind_begin-1))
          return
        end if
      end if

      if (mchl_on) then
        !if (ind >= mchl_ind_begin .and. ind <= mchl_ind_end) then
        !  name = mchl_tracer_names(ind-(mchl_ind_begin-1))
        !  return
        !end if
      end if

      if (iage_on) then
        if (ind >= iage_ind_begin .and. ind <= iage_ind_end) then
          name = iage_tracer_names(ind-(iage_ind_begin-1))
          return
        end if
      end if

      if (tracegas_on) then
        !if (ind >= tracegas_ind_begin .and. ind <= tracegas_ind_end) then
        !  name = tracegas_tracer_names(ind-(tracegas_ind_begin-1))
        !  return
        !end if
      end if

      call exit_POP(sigAbort, 'internal error in passive_tracer_ind2name')

      end subroutine passive_tracer_ind2name

  function tracer_ref_val(ind)

!-----------------------------------------------------------------------
! return reference value for tracer with global tracer index ind
! this is used in virtual flux computations
!
! hooks are only put in for modules that use virtual fluxes
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  result & argument declarations
!-----------------------------------------------------------------------

   real(r8) :: tracer_ref_val

   integer(int_kind), intent(in) :: ind

!-----------------------------------------------------------------------
!  default value for reference value is 0
!-----------------------------------------------------------------------

   tracer_ref_val = c0

!-----------------------------------------------------------------------

  end function tracer_ref_val


!***********************************************************************

 end module passive_tracers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
