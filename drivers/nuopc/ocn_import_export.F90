module ocn_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use POP_KindsMod,          only: POP_i4, POP_r8
  use POP_ErrorMod,          only: POP_ErrorSet, POP_Success
  use POP_FieldMod,          only: POP_fieldKindScalar
  use POP_GridHorzMod,       only: POP_gridHorzLocCenter
  use POP_HaloMod,           only: POP_HaloUpdate
  use kinds_mod,             only: int_kind, r8
  use shr_kind_mod,          only: cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_cal_mod,           only: shr_cal_date2ymd
  use shr_sys_mod,           only: shr_sys_flush, shr_sys_abort
  use shr_const_mod,         only: shr_const_spval
  use shr_string_mod,        only: shr_string_listGetNum, shr_string_listGetName
  use shr_mpi_mod,           only: shr_mpi_min, shr_mpi_max
  use shr_ndep_mod,          only: shr_ndep_readnl
  use communicate,           only: my_task, master_task
  use ocn_communicator,      only: mpi_communicator_ocn
  use domain,                only: distrb_clinic, POP_haloClinic
  use forcing_sfwf,          only: lsend_precip_fact, precip_fact
  use forcing_shf,           only: SHF_QSW
  use forcing_fields,        only: EVAP_F, PREC_F, SNOW_F, MELT_F, ROFF_F, IOFF_F, SALT_F
  use forcing_fields,        only: SENH_F, LWUP_F, LWDN_F, MELTH_F
  use forcing_fields,        only: ATM_CO2_PROG_nf_ind, ATM_CO2_DIAG_nf_ind
  use forcing_fields,        only: ATM_NHx_nf_ind, ATM_NOy_nf_ind
  use forcing_fields,        only: IFRAC, U10_SQR, ATM_PRESS
  use forcing_fields,        only: LAMULT, USTOKES, VSTOKES, LASL
  use forcing_fields,        only: ATM_FINE_DUST_FLUX, ATM_COARSE_DUST_FLUX, SEAICE_DUST_FLUX
  use forcing_fields,        only: ATM_BLACK_CARBON_FLUX, SEAICE_BLACK_CARBON_FLUX
  use mcog,                  only: lmcog, mcog_ncols, lmcog_flds_sent, import_mcog
  use forcing_coupled,       only: update_ghost_cells_coupler_fluxes, rotate_wind_stress
  use ice,                   only: QFLUX, QICE, AQICE, tlast_ice
  use global_reductions,     only: global_sum_prod
  use io_tools,              only: document
  use named_field_mod,       only: named_field_register, named_field_get_index, named_field_set, named_field_get
  use vmix_kpp,              only: KPP_HBLT      ! ocn -> wav, bounadry layer depth
  use grid,                  only: KMT, TAREA
  use nuopc_shr_methods,     only: chkerr
  use ecosys_forcing_mod,    only: ldriver_has_ndep, ldriver_has_atm_co2_diag, ldriver_has_atm_co2_prog
  use constants
  use blocks
  use exit_mod
  use prognostic
  use time_management

  implicit none
  public

  public  :: ocn_advertise_fields
  public  :: ocn_realize_fields
  public  :: ocn_import
  public  :: ocn_export
  public  :: pop_sum_buffer

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_FldChk

  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToOcn_num = 0
  integer                  :: fldsFrOcn_num = 0
  type (fld_list_type)     :: fldsToOcn(fldsMax)
  type (fld_list_type)     :: fldsFrOcn(fldsMax)

  ! area correction factors for fluxes send and received from mediator
  real(r8), allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
  real(r8), allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr

  logical :: ocn2glc_coupling
  integer :: num_ocn2glc_levels
  integer, allocatable :: ocn2glc_levels(:)

  ! accumulated sum of send buffer quantities for averaging before being sent
  real(r8) :: sbuff_sum_u    (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_v    (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_t    (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_s    (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_dhdx (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_dhdy (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_bld  (nx_block,ny_block,max_blocks_clinic)
  real(r8) :: sbuff_sum_co2  (nx_block,ny_block,max_blocks_clinic)
  real(r8), allocatable :: sbuff_sum_t_depth (:,:,:,:)
  real(r8), allocatable :: sbuff_sum_s_depth (:,:,:,:)

  ! tlast_coupled is incremented by delt every time pop_sum_buffer is called
  ! tlast_coupled is reset to 0 when ocn_export is called
  real (r8) :: tlast_coupled


  integer     , parameter :: dbug = 1        ! i/o debug messages
  character(*), parameter :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    character(CS) :: cname
    integer       :: ice_ncat
    logical       :: flds_i2o_per_cat  ! .true. => select per ocn thickness category
    logical       :: flds_co2a
    logical       :: flds_co2b
    logical       :: flds_co2c
    integer       :: ndep_nflds
    character(len=*), parameter :: subname='(ocn_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! pop has not yet initialized my_task and master_task so do it here
    master_task = 0
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=my_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------
    ! optional input from cice columns due to ice thickness categories
    !-----------------
    ! Note that flds_i2o_per_cat is set by the env_run.xml variable CPL_I2O_PER_CAT
    ! This xml variable is set by the POP build-namelist
    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lmcog_flds_sent
    if (my_task == master_task) then
       write(stdout,*) 'lmcog_flds_sent = ',lmcog_flds_sent
    end if

    ! Note that ice_ncat is set by the env_run.xml variable ICE_NCAT which is set
    ! by the ice component (default is 1)
    call NUOPC_CompAttributeGet(gcomp, name='ice_ncat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ice_ncat
    if (my_task == master_task) then
       write(stdout,*) 'ice_ncat = ',ice_ncat
    end if

    if (lmcog_flds_sent) then
       mcog_ncols = ice_ncat+1
    else
       mcog_ncols = 1
    end if

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToOcn_num, fldsToOcn, trim(flds_scalar_name))

    ! from ice
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_melth')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_meltw')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_salt')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_bcpho')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_bcphi')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_flxdst')
    if (lmcog_flds_sent) then
       ! this implementation only handles columns due to ice thickness categories
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sf_afrac')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sf_afracr')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet_afracr' )
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_swpen_ifrac_n', ungridded_lbound=1, ungridded_ubound=ice_ncat)
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac_n'        , ungridded_lbound=1, ungridded_ubound=ice_ncat)
    endif

    ! from river
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofl')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofi')

    ! from mediator
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'So_duu10n')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_tauy')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_taux')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lat')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_sen')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lwup')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_evap')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet')

    ! from wave
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_lamult')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_ustokes')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_vstokes')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sw_hstokes')

    ! from atmosphere
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_pslv')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_lwdn')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_snow')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_rain')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    ! from atm co2 fields
    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    if (my_task == master_task) then
       write(stdout,'(a)') trim(subname)//'flds_co2a = '// trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    if (my_task == master_task) then
       write(stdout,'(a)') trim(subname)//'flds_co2b = '// trim(cvalue)
    end if

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    if (my_task == master_task) then
       write(stdout,'(a)') trim(subname)//'flds_co2c = '// trim(cvalue)
    end if

    if (flds_co2a .or. flds_co2c) then
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2diag')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2prog')
       ldriver_has_atm_co2_prog = .true.
       ldriver_has_atm_co2_diag = .true.
    else
       ldriver_has_atm_co2_prog = .false.
       ldriver_has_atm_co2_diag = .false.
    end if

    ! Determine if will get nitrogen deposition from atm
    call shr_ndep_readnl("drv_flds_in", ndep_nflds)
    if (ndep_nflds > 0) then
       ldriver_has_ndep = .true.
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_ndep'  , ungridded_lbound=1, ungridded_ubound=2)
    else
       ldriver_has_ndep = .false.
    end if

    ! optional per thickness category fields
    do n = 1,fldsToOcn_num
       call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    ! Determine if ocn is sending temperature and salinity data to glc
    call NUOPC_CompAttributeGet(gcomp, name="ocn2glc_coupling", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn2glc_coupling
    if (my_task == master_task) then
       write(stdout,'(a,L1)') trim(subname) // 'ocn2glc coupling is ',ocn2glc_coupling
    end if

    ! Determine number of ocean levels and ocean level indices
    if (ocn2glc_coupling) then
       call NUOPC_CompAttributeGet(gcomp, name="ocn2glc_levels", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       num_ocn2glc_levels = shr_string_listGetNum(cvalue)
       allocate(ocn2glc_levels(num_ocn2glc_levels))
       do n = 1,num_ocn2glc_levels
          call shr_string_listGetName(cvalue, n, cname, rc)
          read(cname,*) ocn2glc_levels(n)
       end do
       if (my_task == master_task) then
          write(stdout,'(a,i0)') trim(subname)//' number of ocean levels sent to glc = ',num_ocn2glc_levels
          write(stdout,*)' ',trim(subname)//' ocean level indices are ',ocn2glc_levels
       end if
    end if

    call fldlist_add(fldsFrOcn_num, fldsFrOcn, trim(flds_scalar_name))
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_omask')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_t')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_u')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_v')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_s')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdx')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdy')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_bldepth')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Fioo_q')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Faoo_fco2_ocn')
    if (ocn2glc_coupling) then
       call fldlist_add(fldsFrOcn_num , fldsFrOcn, 'So_t_depth', &
            ungridded_lbound=1, ungridded_ubound=num_ocn2glc_levels)
       call fldlist_add(fldsFrOcn_num , fldsFrOcn, 'So_s_depth', &
            ungridded_lbound=1, ungridded_ubound=num_ocn2glc_levels)
    end if

    do n = 1,fldsFrOcn_num
       call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_advertise_fields

  !==============================================================================
  subroutine ocn_realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)  , intent(in)  :: mesh
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: flds_scalar_num
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Field)      :: lfield
    integer               :: spatialDim
    integer               :: numOwnedElements
    integer               :: i,j,iblock,n
    type(block)           :: this_block         ! block information for current block
    real(r8), allocatable :: mesh_areas(:)
    real(r8), allocatable :: model_areas(:)
    real(r8), pointer     :: dataptr(:)
    integer               :: num_ocn
    real(r8)              :: max_mod2med_areacor
    real(r8)              :: max_med2mod_areacor
    real(r8)              :: min_mod2med_areacor
    real(r8)              :: min_med2mod_areacor
    real(r8)              :: max_mod2med_areacor_glob
    real(r8)              :: max_med2mod_areacor_glob
    real(r8)              :: min_mod2med_areacor_glob
    real(r8)              :: min_med2mod_areacor_glob
    real(r8), pointer     :: ownedElemCoords(:)
    real(r8), pointer     :: latModel(:), latMesh(:)
    real(r8), pointer     :: lonModel(:), lonMesh(:)
    real(r8)              :: diff_lon
    real(r8)              :: diff_lat
    type(ESMF_StateItem_Flag) :: itemflag
    character(len=*), parameter :: subname='(ocn_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Get import and export states
    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Error checks before realizing fields
    if (lmcog) then
       if ( (.not. State_FldChk(importState, 'Si_ifrac_n'         )) .and. &
            (.not. State_FldChk(importState, 'Fioi_swpen_ifrac_n' )) .and. &
            (.not. State_FldChk(importState, 'Foxx_swnet_afracr'  )) .and. &
            (.not. State_FldChk(importState, 'Sf_afrac'           )) .and. &
            (.not. State_FldChk(importState, 'Sf_afracr'          ))) then

          if (my_task == master_task) then
             write(stdout,*) ' Query for Si_ifrac_n         in import state is ',&
                  State_FldChk(importState, 'Si_ifrac_n')
             write(stdout,*) ' Query for Fioi_swpen_ifrac_n in import state is ',&
                  State_FldChk(importState, 'Foxx_swnet_afracr')
             write(stdout,*) ' Query for Foxx_swnet_afracr  in import state is ',&
                  State_FldChk(importState, 'Foxx_swnet_afracr')
             write(stdout,*) ' Query for Sf_afrac           in import state is ',&
                  State_FldChk(importState, 'Sf_afrac')
             write(stdout,*) ' Query for Sf_afracr          in import state is ',&
                  State_FldChk(importState, 'Sf_afracr')
             write(stdout,*) ' Aborting: all above import fields must be in import state if lmcog is true'
          end if
          call shr_sys_abort(trim(subname)//": all import fields not set if lmcog is .true.")
       end if
    end if

    ! Realize import and export states
    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrOcn, &
         numflds=fldsFrOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':POP_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToOcn, &
         numflds=fldsToOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':POP_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    allocate(lonMesh(numOwnedElements))
    allocate(latMesh(numOwnedElements))
    allocate(lonModel(numOwnedElements))
    allocate(latModel(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

    ! Compare mesh lats/lons to model generated lats/lons
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n+1
             lonModel(n) = tlond(i,j,iblock)
             latModel(n) = tlatd(i,j,iblock)
             diff_lon = abs(lonMesh(n) - lonModel(n))
             if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360.) > 1.e-1) .or.&
                  (diff_lon > 1.e-3 .and. diff_lon < c1) ) then
                write(6,'(a,i6,2(f21.13,3x),d21.5)') &
                     'ERROR: POP n, lonMesh, lonModel, diff_lon = ',&
                     n,lonMesh(n),lonModel(n), diff_lon
                call shr_sys_abort()
             end if
             if (abs(latMesh(n) - latModel(n)) > 1.e-1) then
                write(6,'(a,i6,2(f21.13,3x),d21.5)') &
                     'ERROR: POP n, latMesh, latModel, diff_lat = ', &
                     n,latMesh(n),latModel(n), abs(latMesh(n)-latModel(n))
                call shr_sys_abort()
             end if
          enddo
       enddo
    enddo

    ! Determine mesh areas used in regridding
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)
    call ESMF_FieldDestroy(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine flux correction factors (module variables)
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    mod2med_areacor(:) = 1._r8
    med2mod_areacor(:) = 1._r8
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n+1
             model_areas(n) = tarea(i,j,iblock)/(radius*radius)
             mod2med_areacor(n) = model_areas(n) / mesh_areas(n)
             med2mod_areacor(n) = mesh_areas(n) / model_areas(n)
          end do
       end do
    end do
    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, mpi_communicator_ocn)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, mpi_communicator_ocn)

    if (my_task == master_task) then
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'POP'
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'POP'
    end if

    deallocate(model_areas)
    deallocate(mesh_areas)

    call ESMF_StateGet(importState, 'Faxa_ndep', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call ESMF_LogWrite(subname//' Faxa_ndep is in import state', ESMF_LOGMSG_INFO)
    end if

  end subroutine ocn_realize_fields

  !==============================================================================
  subroutine ocn_import( importState, flds_scalar_name, ldiag_cpl, errorCode, rc )

    !-----------------------------------------------------------------------
    ! swnet  -- net short-wave heat flux                 (W/m2   )
    ! lwup   -- longwave radiation (up)                  (W/m2   )
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)   , intent(in)  :: importState
    character(len=*)   , intent(in)  :: flds_scalar_name
    logical (log_kind) , intent(inout)  :: ldiag_cpl
    integer            , intent(out) :: errorCode
    integer            , intent(out) :: rc

    ! local variables
    type (block)         :: this_block ! local block info
    character (char_len) :: label,  message
    real (r8)            :: work1(nx_block,ny_block,max_blocks_clinic)
    real (r8)            :: work2(nx_block,ny_block,max_blocks_clinic)  ! local work space
    integer (int_kind)   :: i,j,k,n,ncol,iblock,nfld,nf
    real (r8)            :: m2percm2, gsum
    real (r8), pointer   :: dataptr1d(:)
    real (r8), pointer   :: dataptr2d(:,:)
    ! from mediator (virtual ocn)
    real (r8), pointer   :: foxx_swnet(:)
    real (r8), pointer   :: foxx_swnet_afracr(:)
    real (r8), pointer   :: foxx_taux(:)
    real (r8), pointer   :: foxx_tauy(:)
    real (r8), pointer   :: foxx_lwup(:)
    real (r8), pointer   :: foxx_sen(:)
    real (r8), pointer   :: foxx_evap(:)
    real (r8), pointer   :: foxx_rofl(:)
    real (r8), pointer   :: foxx_rofi(:)
    real (r8), pointer   :: so_duu10n(:)
    ! from atm
    real (r8), pointer   :: sa_pslv(:)
    real (r8), pointer   :: faxa_rain(:)
    real (r8), pointer   :: faxa_snow(:)
    real (r8), pointer   :: faxa_lwdn(:)
    real (r8), pointer   :: faxa_dstwet(:,:)
    real (r8), pointer   :: faxa_dstdry(:,:)
    real (r8), pointer   :: faxa_bcph(:,:)
    ! from ice
    real (r8), pointer   :: sf_afrac(:)
    real (r8), pointer   :: sf_afracr(:)
    real (r8), pointer   :: si_ifrac(:)
    real (r8), pointer   :: si_ifrac_n(:,:)
    real (r8), pointer   :: fioi_flxdst(:)
    real (r8), pointer   :: fioi_bcpho(:)
    real (r8), pointer   :: fioi_bcphi(:)
    real (r8), pointer   :: fioi_meltw(:)
    real (r8), pointer   :: fioi_melth(:)
    real (r8), pointer   :: fioi_salt(:)
    real (r8), pointer   :: fioi_swpen_ifrac_n(:,:)
    real (r8)            :: frac_col_1pt(mcog_ncols)
    real (r8)            :: fracr_col_1pt(mcog_ncols)
    real (r8)            :: qsw_fracr_col_1pt(mcog_ncols)
    ! from wave
    real (r8), pointer   :: sw_lamult(:)
    real (r8), pointer   :: sw_ustokes(:)
    real (r8), pointer   :: sw_vstokes(:)
    real (r8), pointer   :: sw_hstokes(:)
    !
    integer (int_kind)   :: fieldCount
    character (char_len) :: fldname
    type(ESMF_StateItem_Flag) :: itemflag
#ifdef _HIRES
    real (r8)            :: qsw_eps = -1.e-3_r8
#else
    real (r8)            :: qsw_eps = 0._r8
#endif
    character (char_len), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_import)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    !  zero out padded cells
    !-----------------------------------------------------------------------

    work1 = c0
    work2 = c0

    ! NOTE: if there are code changes associated with changing the names or
    !       the number of fluxes received from the coupler, then subroutine
    !       update_ghost_cells_coupler_fluxes will need to be modified also

    ! NOTE : RCALCT acts as a KMT mask, its 1 if KMT>=1 and 0 otherwise

    !-----------------------------------------------------------------------
    ! from mediator (virtual ocean)
    !-----------------------------------------------------------------------

    !  unpack and distribute wind stress, then convert to correct units
    !  and rotate components to local coordinates

    ! zonal wind and meridional wind stress  (W/m2)
    call state_getfldptr(importState, 'Foxx_taux', foxx_taux, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Foxx_tauy', foxx_tauy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n + 1
             work1(i,j,iblock)  = foxx_taux(n) * med2mod_areacor(n)
             work2(i,j,iblock)  = foxx_tauy(n) * med2mod_areacor(n)
          end do
       end do
    end do
    ! rotate true zonal/meridional wind stress into local coordinates,
    ! convert to dyne/cm**2, and shift SMFT to U grid
    ! halo updates are performed in subroutine rotate_wind_stress,
    ! following the rotation
    call rotate_wind_stress(work1, work2)

    ! evaporation flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_evap', foxx_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! sensible heat flux (W/m2)
    call state_getfldptr(importState, 'Foxx_sen', foxx_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! long wave up flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_lwup', foxx_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! shortwave net flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_swnet', foxx_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! 10m wind speed squared (m^2/s^2)
    call state_getfldptr(importState, 'So_duu10n', so_duu10n, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n + 1
             EVAP_F(i,j,iblock) = foxx_evap(n) * med2mod_areacor(n)
             SENH_F(i,j,iblock) = foxx_sen(n) * med2mod_areacor(n)
             LWUP_F(i,j,iblock) = foxx_lwup(n) * med2mod_areacor(n)
             !  convert from W/m**2
             SHF_QSW(i,j,iblock) = foxx_swnet(n) * med2mod_areacor(n) * RCALCT(i,j,iblock)*hflux_factor
             ! convert from m**2/s**2 to cm**2/s**2
             U10_SQR(i,j,iblock) = cmperm * cmperm * so_duu10n(n) * RCALCT(i,j,iblock)
         end do
       end do
    end do
    if (ANY(SHF_QSW < qsw_eps)) then
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                write(6,*)'ERROR: j,i,shf_qsw = ',&
                     this_block%j_glob(j),this_block%i_glob(i),SHF_QSW(i,j,iblock)
             enddo
          enddo
       enddo
       call shr_sys_abort('(set_surface_forcing) ERROR: SHF_QSW < qsw_eps in set_surface_forcing')
    endif

    !-----------------------------------------------------------------------
    ! from atmosphere
    !-----------------------------------------------------------------------

    ! sea-level pressure (Pa)
    call state_getfldptr(importState, 'Sa_pslv', sa_pslv, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow
    call state_getfldptr(importState, 'Faxa_snow', faxa_snow, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! rain
    call state_getfldptr(importState, 'Faxa_rain', faxa_rain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! longwave radiation (down) (W/m2)
    call state_getfldptr(importState, 'Faxa_lwdn', faxa_lwdn, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! wet dust flux from atm (convert from from MKS (kg/m^2/s) to CGS (g/cm^2/s))
    call state_getfldptr(importState, 'Faxa_dstwet', faxa_dstwet, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! dry dust flux from atm (convert from from MKS (kg/m^2/s) to CGS (g/cm^2/s))
    call state_getfldptr(importState, 'Faxa_dstdry', faxa_dstdry, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! black carbon flux from atm (convert from MKS (kg/m^2/s) to CGS (g/cm^2/s))
    call state_getfldptr(importState, 'Faxa_bcph', faxa_bcph, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n+1
             ATM_PRESS(i,j,iblock) = c10 * sa_pslv(n) * RCALCT(i,j,iblock) ! convert from Pa to dynes/cm**2
             SNOW_F(i,j,iblock) = faxa_snow(n) * med2mod_areacor(n)
             PREC_F(i,j,iblock) = (faxa_rain(n) + faxa_snow(n)) * med2mod_areacor(n) ! rain + snow
             LWDN_F(i,j,iblock) = faxa_lwdn(n) * med2mod_areacor(n)
             ATM_FINE_DUST_FLUX(i,j,iblock)   = (faxa_dstwet(1,n) + faxa_dstdry(1,n)) * med2mod_areacor(n)
             ATM_COARSE_DUST_FLUX(i,j,iblock) = (faxa_dstwet(2,n) + faxa_dstdry(2,n) + &
                                                 faxa_dstwet(3,n) + faxa_dstdry(3,n) + &
                                                 faxa_dstwet(4,n) + faxa_dstdry(4,n)) * med2mod_areacor(n)
             ATM_BLACK_CARBON_FLUX(i,j,iblock) = (faxa_bcph(1,n) + faxa_bcph(2,n) + faxa_bcph(3,n)) * med2mod_areacor(n)
             ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s)
             ATM_FINE_DUST_FLUX(i,j,iblock)    = ATM_FINE_DUST_FLUX(i,j,iblock)    * 0.1_r8 * RCALCT(i,j,iblock)
             ATM_COARSE_DUST_FLUX(i,j,iblock)  = ATM_COARSE_DUST_FLUX(i,j,iblock)  * 0.1_r8 * RCALCT(i,j,iblock)
             ATM_BLACK_CARBON_FLUX(i,j,iblock) = ATM_BLACK_CARBON_FLUX(i,j,iblock) * 0.1_r8 * RCALCT(i,j,iblock)
          end do
       end do
    end do

    !-----------------------------------------------------------------------
    ! from sea-ice
    !-----------------------------------------------------------------------

    ! ice fraction
    call state_getfldptr(importState, 'Si_ifrac', si_ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! heat flux from sea ice snow & ice melt (W/m2)
    call state_getfldptr(importState, 'Fioi_melth', fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! salt from sea ice (kg(salt)/m2/s)
    call state_getfldptr(importState, 'Fioi_salt', fioi_salt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow melt flux from sea ice (kg/m2/s)
    call state_getfldptr(importState, 'Fioi_meltw', fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! dust flux from sea ice
    call state_getfldptr(importState, 'Fioi_flxdst', fioi_flxdst, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! black carbon flux from sea ice
    call state_getfldptr(importState, 'Fioi_bcpho', fioi_bcpho, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Fioi_bcphi', fioi_bcphi, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n + 1
             IFRAC(i,j,iblock) = si_ifrac(n) * RCALCT(i,j,iblock)
             MELTH_F(i,j,iblock) = fioi_melth(n) * med2mod_areacor(n)
             MELT_F(i,j,iblock) = fioi_meltw(n) * med2mod_areacor(n)
             SALT_F(i,j,iblock) = fioi_salt(n) * med2mod_areacor(n)
             SEAICE_DUST_FLUX(i,j,iblock) = fioi_flxdst(n)
             ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s) and apply flux area correction factor
             SEAICE_DUST_FLUX(i,j,iblock) = SEAICE_DUST_FLUX(i,j,iblock) * &
                                            0.1_r8 * RCALCT(i,j,iblock) * med2mod_areacor(n)
             SEAICE_BLACK_CARBON_FLUX(i,j,iblock) = fioi_bcpho(n) + fioi_bcphi(n)
             ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s) and apply flux area correction factor
             SEAICE_BLACK_CARBON_FLUX(i,j,iblock) = SEAICE_BLACK_CARBON_FLUX(i,j,iblock) * &
                                                    0.1_r8 * RCALCT(i,j,iblock) * med2mod_areacor(n)
          end do
       end do
    end do

    ! net sw to ocn needed below
    call state_getfldptr(importState, 'Foxx_swnet', foxx_swnet, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,size(Foxx_swnet)
       foxx_swnet(n) = foxx_swnet(n) * med2mod_areacor(n)
    end do

    if (lmcog) then
       ! extract fields for each column and pass to import_mcog

       call state_getfldptr(importState, 'Sf_afrac', sf_afrac, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Sf_afracr', sf_afracr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_swnet_afracr', Foxx_swnet_afracr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(Foxx_swnet_afracr)
          foxx_swnet_afracr(n) = Foxx_swnet_afracr(n) * med2mod_areacor(n)
       end do
       call state_getfldptr(importState, 'Si_ifrac_n', Si_ifrac_n, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Fioi_swpen_ifrac_n', Fioi_swpen_ifrac_n, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,size(Fioi_swpen_ifrac_n, dim=2)
          fioi_swpen_ifrac_n(:,n) = Fioi_swpen_ifrac_n(:,n) * med2mod_areacor(n) !
       end do

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n = n + 1
                frac_col_1pt(1)  = max(c0, min(c1, Sf_afrac(n)))
                fracr_col_1pt(1) = max(c0, min(c1, Sf_afracr(n)))
                qsw_fracr_col_1pt(1) = Foxx_swnet_afracr(n)
                do ncol = 2,mcog_ncols ! mcog_ncols is the same as ice_ncat+1
                   frac_col_1pt(ncol)  = max(c0, min(c1, Si_ifrac_n(ncol-1,n)))
                   fracr_col_1pt(ncol) = max(c0, min(c1, Si_ifrac_n(ncol-1,n)))
                   qsw_fracr_col_1pt(ncol) = Fioi_swpen_ifrac_n(ncol-1,n)
                end do
                call import_mcog(frac_col_1pt, fracr_col_1pt, qsw_fracr_col_1pt, Foxx_swnet(n), iblock, i, j)
             end do
          end do
       end do

    else

       ! if mcog is off, fill its arrays with cpl aggregated full cell means
       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n = n + 1
                ncol = 1
                frac_col_1pt(ncol)  = c1
                fracr_col_1pt(ncol) = c1
                qsw_fracr_col_1pt(ncol) = Foxx_swnet(n)
                call import_mcog(frac_col_1pt, fracr_col_1pt, qsw_fracr_col_1pt, Foxx_swnet(n), iblock, i, j)
             enddo
          enddo
       enddo

    endif ! if (lmcog) then

    !-----------------------------------------------------------------------
    ! from wave
    !-----------------------------------------------------------------------

    ! enhancement factor (unitless)
    call state_getfldptr(importState, 'Sw_lamult', sw_lamult, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Stokes drift (m/s)
    call state_getfldptr(importState, 'Sw_ustokes', sw_ustokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! Stokes drift (m/s)
    call state_getfldptr(importState, 'Sw_vstokes', sw_vstokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! surface layer Lanmguir number (unitless)
    call state_getfldptr(importState, 'Sw_hstokes', sw_hstokes, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n+1
             if (IFRAC(i,j,iblock) <= 0.05_r8) then
                LAMULT(i,j,iblock) = sw_lamult(n)  * RCALCT(i,j,iblock)
                LASL(i,j,iblock)   = sw_hstokes(n) * RCALCT(i,j,iblock)
             else
                LAMULT(i,j,iblock) =  c1
                LASL(i,j,iblock)   = -c1
             end if
             USTOKES(i,j,iblock) = sw_ustokes(n) * RCALCT(i,j,iblock)
             VSTOKES(i,j,iblock) = sw_vstokes(n) * RCALCT(i,j,iblock)
          end do
       end do
    end do

    !-----------------------------------------------------------------------
    ! from river
    !-----------------------------------------------------------------------

    ! liquid runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofl', foxx_rofl, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ice runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofi', foxx_rofi, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n+1
             ROFF_F(i,j,iblock) = foxx_rofl(n) * med2mod_areacor(n)
             IOFF_F(i,j,iblock) = foxx_rofi(n) * med2mod_areacor(n)
          end do
       end do
    end do

    !-----------------------------------------------------------------------
    ! update ghost cells for fluxes received from the mediator
    !-----------------------------------------------------------------------

    call update_ghost_cells_coupler_fluxes(errorCode)

    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, &
            'ocn_import: error in update_ghost_cells_coupler_fluxes')
       return
    endif

    !-----------------------------------------------------------------------
    ! CO2 from atm
    !-----------------------------------------------------------------------

    ! co2prog-- bottom atm level prognostic co2
    call ESMF_StateGet(importState, 'Sa_co2prog', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getimport(importState, 'Sa_co2prog', work1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call POP_HaloUpdate(work1,POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating PROG CO2 halo')
          return
       endif
       call named_field_set(ATM_CO2_PROG_nf_ind, work1)
    endif

    ! co2diag-- bottom atm level diagnostic co2
    call ESMF_StateGet(importState, 'Sa_co2diag', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getimport(importState, 'Sa_co2diag', work1, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call POP_HaloUpdate(work1,POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating DIAG CO2 halo')
          return
       endif
       call named_field_set(ATM_CO2_DIAG_nf_ind, work1)
    endif

    call ESMF_StateGet(importState, 'Faxa_ndep', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getfldptr(importState, 'Faxa_ndep', dataptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Note that nhx is ungridded_index=1, nhy is ungridded_index = 2
       ! Note - the input units are kgN/m2/s to nmolN/cm2/s
       ! TODO: Keith has pointed out might want to use 14.007_r8 instead of 14.0_r8 for more
       ! consistency when bringing in N isotopes into the code

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                n = n+1
                work1(i,j,iblock) = dataptr2d(1,n) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8) * med2mod_areacor(n)
             end do
          end do
       end do
       call POP_HaloUpdate(work1, POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating DIAG NHx halo')
          return
       endif
       call named_field_set(ATM_NHx_nf_ind, work1)

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                n = n+1
                work1(i,j,iblock) = dataptr2d(2,n) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8) * med2mod_areacor(n)
             end do
          end do
       end do
       call POP_HaloUpdate(work1, POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating DIAG NOy halo')
          return
       endif
       call named_field_set(ATM_NOy_nf_ind, work1)
    endif

    !-----------------------------------------------------------------------
    !  diagnostics
    !-----------------------------------------------------------------------

    if (ldiag_cpl) then
       write(message,'(6a,1x,5a)')  &
            ' Global averages of fluxes received from cpl at ',  &
            cyear,'/',cmonth ,'/',cday,  chour,':',cminute,':',csecond
       call document ('pop_recv_from_coupler', trim(message))

       ! Determine all field names in import state
       call ESMF_StateGet(importState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in import state
       ! from atm - black carbon deposition fluxes (3)
       ! (1) => Faxa_bcphidry, (2) => Faxa_bcphodry, (3) => Faxa_bcphiwet
       ! from atm - organic carbon deposition fluxes (3)
       ! (1) => Faxa_dstwet1, (2) => Faxa_dstwet2, (3) => Faxa_dstwet3, (4) => Faxa_dstwet4
       ! from atm - dry dust deposition frluxes (4 sizes)
       ! (1) => Faxa_dstdry1, (2) => Faxa_dstdry2, (3) => Faxa_dstdry3, (4) => Faxa_dstdry4

       m2percm2  = mpercm*mpercm
       do nfld = 1, fieldCount
          if (fieldNameList(nfld) == trim(flds_scalar_name)) then
             CYCLE
          else if (fieldNameList(nfld) == 'Faxa_bcph') then
             call state_getfldptr(importState, 'Faxa_bcph', dataPtr2d, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do nf = 1,3
                if (nf == 1) fldname = 'Faxa_bcphidry'
                if (nf == 2) fldname = 'Faxa_bcphodry'
                if (nf == 3) fldname = 'Faxa_bcphiwet'
                n = 0
                do iblock = 1, nblocks_clinic
                   this_block = get_block(blocks_clinic(iblock),iblock)
                   do j = this_block%jb,this_block%je
                      do i = this_block%ib,this_block%ie
                         n = n + 1
                         work1(i,j,iblock) = dataPtr2d(nf,n)
                      end do
                   end do
                end do
                gsum = global_sum_prod(work1, TAREA, distrb_clinic,  field_loc_center, RCALCT) * m2percm2
                if (my_task == master_task) then
                   write(stdout,1100)'ocn','recv', trim(fldname), gsum
                   call shr_sys_flush(stdout)
                endif
             end do
          else if (fieldNameList(nfld) == 'Faxa_dstwet') then
             call state_getfldptr(importState, 'Faxa_dstwet', dataPtr2d, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do nf = 1,4
                if (nf == 1) fldname = 'Faxa_dstwet1'
                if (nf == 2) fldname = 'Faxa_dstwet2'
                if (nf == 3) fldname = 'Faxa_dstwet3'
                if (nf == 4) fldname = 'Faxa_dstwet4'
                n = 0
                do iblock = 1, nblocks_clinic
                   this_block = get_block(blocks_clinic(iblock),iblock)
                   do j = this_block%jb,this_block%je
                      do i = this_block%ib,this_block%ie
                         n = n + 1
                         work1(i,j,iblock) = dataPtr2d(nf,n)
                      end do
                   end do
                end do
                gsum = global_sum_prod(work1, TAREA, distrb_clinic,  field_loc_center, RCALCT) * m2percm2
                if (my_task == master_task) then
                   write(stdout,1100)'ocn','recv', trim(fldname), gsum
                   call shr_sys_flush(stdout)
                endif
             end do
          else if (fieldNameList(nfld) == 'Faxa_dstdry') then
             call state_getfldptr(importState, 'Faxa_dstdry', dataPtr2d, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             do nf = 1,4
                if (nf == 1) fldname = 'Faxa_dstdry1'
                if (nf == 2) fldname = 'Faxa_dstdry2'
                if (nf == 3) fldname = 'Faxa_dstdry3'
                if (nf == 4) fldname = 'Faxa_dstdry4'
                n = 0
                do iblock = 1, nblocks_clinic
                   this_block = get_block(blocks_clinic(iblock),iblock)
                   do j = this_block%jb,this_block%je
                      do i = this_block%ib,this_block%ie
                         n = n + 1
                         work1(i,j,iblock) = dataPtr2d(nf,n)
                      end do
                   end do
                end do
                gsum = global_sum_prod(work1, TAREA, distrb_clinic,  field_loc_center, RCALCT) * m2percm2
                if (my_task == master_task) then
                   write(stdout,1100)'ocn','recv', trim(fldname), gsum
                   call shr_sys_flush(stdout)
                endif
             end do
          else if (fieldNameList(nfld) == 'Fioi_swpen_ifrac_n' .or. fieldNameList(nfld) == 'Si_ifrac_n') then
             ! do nothing for now
          else
             call state_getfldptr(importState, trim(fieldNameList(nfld)), dataPtr1d, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             n = 0
             do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock),iblock)
                do j = this_block%jb,this_block%je
                   do i = this_block%ib,this_block%ie
                      n = n + 1
                      work1(i,j,iblock) = dataPtr1d(n)
                   end do
                end do
             end do
             gsum = global_sum_prod(work1, TAREA, distrb_clinic,  field_loc_center, RCALCT) * m2percm2
             if (my_task == master_task) then
                write(stdout,1100)'ocn','recv', trim(fieldNameList(nfld)), gsum
                call shr_sys_flush(stdout)
             endif
          end if
       end do
    end if
1100 format ('comm_diag ', a3, 1x, a4, 1x, a15, 1x, es26.19:, 1x, a6)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_import

  !==============================================================================
  subroutine ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)

    !-----------------------------------------------------------------------
    ! Create export state
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                 :: exportState
    character(len=*)   , intent(in)  :: flds_scalar_name
    logical (log_kind) , intent(inout)  :: ldiag_cpl
    integer (POP_i4)   , intent(out) :: errorCode  ! pop error code
    integer            , intent(out) :: rc         ! returned error code

    ! local variables
    type (block)         :: this_block ! local block info
    integer (int_kind)   :: n,i,j,k,iblock,nfld,lev
    character (char_len) :: label
    real (r8)            :: work1(nx_block,ny_block)
    real (r8)            :: work2(nx_block,ny_block)
    real (r8)            :: work3(nx_block,ny_block)
    real (r8)            :: work4(nx_block,ny_block)
    real (r8)            :: worka(nx_block,ny_block,max_blocks_clinic)
    real (r8)            :: m2percm2
    real (r8)            :: gsum
    real (r8), pointer   :: dataptr1(:)
    real (r8), pointer   :: dataptr2(:)
    real (r8), pointer   :: dataptr2d(:,:)
    integer (int_kind)   :: fieldCount
    character (char_len), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_export)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    ! ocean mask
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_omask', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dataptr1(:) = c0
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n+1
             dataptr1(n) = float(KMT(i,j,iblock))
             if (dataptr1(n) > 1.0_r8) dataptr1(n) = 1.0_r8
          enddo
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! interpolate onto T-grid points and rotate on T grid
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_u', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_v', dataPtr2, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    dataptr2(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)

       call ugrid_to_tgrid(work3, sbuff_sum_u(:,:,iblock),iblock)
       call ugrid_to_tgrid(work4, sbuff_sum_v(:,:,iblock),iblock)

       work1 = (work3*cos(ANGLET(:,:,iblock)) + work4*sin(-ANGLET(:,:,iblock))) * mpercm/tlast_coupled
       work2 = (work4*cos(ANGLET(:,:,iblock)) - work3*sin(-ANGLET(:,:,iblock))) * mpercm/tlast_coupled

       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = work1(i,j)
             dataptr2(n) = work2(i,j)
          enddo
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! convert and pack surface temperature
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_t', dataptr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dataptr1(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n)= sbuff_sum_t(i,j,iblock)/tlast_coupled + T0_Kelvin
          enddo
       enddo
    enddo
    if (ocn2glc_coupling) then
       call state_getfldptr(exportState, 'So_t_depth', dataptr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr2d(:,:) = c0
       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             do lev = 1,num_ocn2glc_levels
                if (KMT(i,j,iblock) >= ocn2glc_levels(lev)) then
                   dataptr2d(lev,n) = sbuff_sum_t_depth(i,j,iblock,lev)/tlast_coupled + T0_Kelvin
                end if
             end do
          enddo
          enddo
       enddo
    endif

    !-----------------------------------------------------------------------
    ! convert and pack salinity
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_s', dataptr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = sbuff_sum_s(i,j,iblock) * salt_to_ppt / tlast_coupled
          enddo
       enddo
    enddo
    if (ocn2glc_coupling) then
       call state_getfldptr(exportState, 'So_s_depth', dataptr2d, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr2d(:,:) = c0
       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             do lev = 1,num_ocn2glc_levels
                if (KMT(i,j,iblock) >= ocn2glc_levels(lev)) then
                   dataptr2d(lev,n) = sbuff_sum_s_depth(i,j,iblock,lev)*salt_to_ppt/tlast_coupled
                end if
             end do
          end do
          end do
       end do
    end if

    !-----------------------------------------------------------------------
    ! convert and pack boundary layer depth
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_bldepth', dataptr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = sbuff_sum_bld(i,j,iblock)/100./tlast_coupled
          enddo
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! interpolate onto T-grid points, then rotate on T grid
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_dhdx', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_dhdy', dataPtr2, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)

       call ugrid_to_tgrid(work3, sbuff_sum_dhdx(:,:,iblock), iblock)
       call ugrid_to_tgrid(work4, sbuff_sum_dhdy(:,:,iblock), iblock)

       work1 = (work3*cos(ANGLET(:,:,iblock)) + work4*sin(-ANGLET(:,:,iblock))) /grav /tlast_coupled
       work2 = (work4*cos(ANGLET(:,:,iblock)) - work3*sin(-ANGLET(:,:,iblock))) /grav /tlast_coupled

       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = work1(i,j)
             dataptr2(n) = work2(i,j)
          enddo
       enddo
    enddo

    !-----------------------------------------------------------------------
    ! pack heat flux due to freezing/melting (W/m^2)
    ! QFLUX computation and units conversion occurs in ice.F
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'Fioo_q', dataPtr1, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    dataptr1(:) = shr_const_spval
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = QFLUX(i,j,iblock) * mod2med_areacor(n)
          enddo
       enddo
    enddo

    tlast_ice = c0
    AQICE     = c0
    QICE      = c0

    !-----------------------------------------------------------------------
    ! pack co2 flux, if requested (kg CO2/m^2/s)
    ! units conversion occurs where co2 flux is computed
    !-----------------------------------------------------------------------

    if ( State_FldChk(exportState, 'Faoo_fco2_ocn')) then
       call state_getfldptr(exportState, 'Faoo_fco2_ocn', dataPtr1, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dataptr1(:) = shr_const_spval

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n = n + 1
                dataptr1(n) = (sbuff_sum_co2(i,j,iblock)/tlast_coupled) * mod2med_areacor(n)
             enddo
          enddo
       enddo
    endif

    !-----------------------------------------------------------------------
    ! diagnostics
    !-----------------------------------------------------------------------

    if (ldiag_cpl) then
       call ccsm_char_date_and_time
       if (my_task == master_task) then
          write(stdout,*)'pop_send_to_mediator'
       end if

       ! Determine all field names in export state
       call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in export state
       m2percm2  = mpercm*mpercm
       do nfld = 1, fieldCount
          if (trim(fieldNameList(nfld)) /= flds_scalar_name) then
             call state_getfldptr(exportState, trim(fieldNameList(nfld)), dataptr1, rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             n = 0
             do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock),iblock)
                do j = this_block%jb,this_block%je
                   do i = this_block%ib,this_block%ie
                      n = n + 1
                      worka(i,j,iblock) = dataptr1(n)
                   enddo
                enddo
             enddo

             call POP_HaloUpdate(worka, POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
                  errorCode, fillValue = 0.0_POP_r8)
             if (errorCode /= POP_Success) then
                call POP_ErrorSet(errorCode, 'ocn_export: error updating halo for state')
                rc = ESMF_FAILURE
                return
             endif

             gsum = global_sum_prod(worka, TAREA, distrb_clinic, field_loc_center, RCALCT) * m2percm2

             if (my_task == master_task) then
                write(stdout,1100)'ocn','send', trim(fieldNameList(nfld)), gsum
                call shr_sys_flush(stdout)
             endif
          end if
       end do
    end if
1100 format ('comm_diag ', a3, 1x, a4, 1x, a15, 1x, es26.19:, 1x, a6)

    tlast_coupled = c0

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_export

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC, only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(ocn_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                                         ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                                         ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                                         gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(fldlist_realize:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !==============================================================================
  subroutine state_getimport(state, fldname, output, areacor, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    real (r8)           , intent(inout) :: output(:,:,:)
    real(r8) , optional , intent(in)    :: areacor(:)
    integer             , intent(out)   :: rc

    ! local variables
    type(block)       :: this_block         ! block information for current block
    integer           :: i, j, iblock, n   ! incides
    real(r8), pointer :: dataPtr1d(:)
    character(len=*), parameter :: subname='(import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call state_getfldptr(state, trim(fldname), dataptr1d, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       if (present(areacor)) then
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                n = n + 1
                output(i,j,iblock) = dataPtr1d(n) * areacor(n)
             end do
          end do
       else
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                n = n + 1
                output(i,j,iblock) = dataPtr1d(n)
             end do
          end do
       end if
    end do

  end subroutine state_getimport

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 1d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(r8), pointer , intent(inout)  :: fldptr(:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_1d

  !===============================================================================
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 2d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(r8), pointer , intent(inout)  :: fldptr(:,:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_2d

  !===============================================================================
  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)
    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

  !===============================================================================

  subroutine pop_sum_buffer(exportState, rc)

    ! ----------------------------------------------
    ! Accumulates sums for averaging fields to be sent to mediator
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)     :: exportState
    integer, intent(out) :: rc

    ! local variables
    type (block)             :: this_block ! local block info
    real (r8)                :: work(nx_block,ny_block,max_blocks_clinic) ! local work arrays
    real (r8)                :: delt                                      ! time interval since last step
    real (r8)                :: delt_last                                 ! time interval for previous step
    integer (int_kind)       :: iblock                                    ! block index
    integer (int_kind)       :: sflux_co2_nf_ind = 0                      ! named field index of fco2
    integer (int_kind)       :: i,j,n,lev                                 ! indices
    logical (log_kind), save :: first = .true.                            ! only true for first call

    !-----------------------------------------------------------------------
    ! zero buffer if this is the first time after a coupling interval
    !-----------------------------------------------------------------------

    if (tlast_coupled == c0) then
       sbuff_sum_u    (:,:,:) = c0
       sbuff_sum_v    (:,:,:) = c0
       sbuff_sum_t    (:,:,:) = c0
       sbuff_sum_s    (:,:,:) = c0
       sbuff_sum_dhdx (:,:,:) = c0
       sbuff_sum_dhdy (:,:,:) = c0
       sbuff_sum_bld  (:,:,:) = c0
       sbuff_sum_co2  (:,:,:) = c0
       if (.not. allocated(sbuff_sum_t_depth)) then
          allocate(sbuff_sum_t_depth (nx_block,ny_block,max_blocks_clinic,num_ocn2glc_levels))
       end if
       sbuff_sum_t_depth(:,:,:,:) = c0
       if (.not. allocated(sbuff_sum_s_depth)) then
          allocate(sbuff_sum_s_depth (nx_block,ny_block,max_blocks_clinic,num_ocn2glc_levels))
       end if
       sbuff_sum_s_depth(:,:,:,:) = c0
    end if

    work = c0

    !-----------------------------------------------------------------------
    ! update time since last coupling
    !-----------------------------------------------------------------------

    if (avg_ts .or. back_to_back) then
       delt = p5*dtt
    else
       delt =    dtt
    endif
    tlast_coupled = tlast_coupled + delt

    !-----------------------------------------------------------------------
    ! allow for fco2 field to not be registered on first call
    !    because init_forcing is called before init_passive_tracers
    ! use weight from previous timestep because flux used here is that
    !    computed during the previous timestep
    !-----------------------------------------------------------------------

    if ( State_FldChk(exportState, 'Faoo_fco2_ocn')) then
       if (sflux_co2_nf_ind == 0) then
          call named_field_get_index('SFLUX_CO2', sflux_co2_nf_ind, exit_on_err=.not. first)
       endif

       if (avg_ts .or. back_to_back) then
          delt_last = p5*dtt
       else
          delt_last =    dtt
       endif
    endif

    !-----------------------------------------------------------------------
    !  accumulate sums of U,V,T,S and GRADP
    !  accumulate sum of co2 flux, if requested
    !     implicitly use zero flux if fco2 field not registered yet
    !  ice formation flux is handled separately in ice routine
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       sbuff_sum_u    (:,:,iblock) = sbuff_sum_u   (:,:,iblock) + delt * UVEL(:,:,1,curtime,iblock)
       sbuff_sum_v    (:,:,iblock) = sbuff_sum_v   (:,:,iblock) + delt * VVEL(:,:,1,curtime,iblock)
       sbuff_sum_t    (:,:,iblock) = sbuff_sum_t   (:,:,iblock) + delt * TRACER(:,:,1,1,curtime,iblock)
       sbuff_sum_s    (:,:,iblock) = sbuff_sum_s   (:,:,iblock) + delt * TRACER(:,:,1,2,curtime,iblock)
       sbuff_sum_dhdx (:,:,iblock) = sbuff_sum_dhdx(:,:,iblock) + delt * GRADPX(:,:,curtime,iblock)
       sbuff_sum_dhdy (:,:,iblock) = sbuff_sum_dhdy(:,:,iblock) + delt * GRADPY(:,:,curtime,iblock)
       sbuff_sum_bld  (:,:,iblock) = sbuff_sum_bld (:,:,iblock) + delt * KPP_HBLT(:,:,iblock)
    end do
    if (ocn2glc_coupling) then
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j = this_block%jb,this_block%je
             do i = this_block%ib,this_block%ie
                do n = 1,num_ocn2glc_levels
                   lev = ocn2glc_levels(n)
                   if (KMT(i,j,iblock) >= lev) then
                      sbuff_sum_t_depth(i,j,iblock,n) = sbuff_sum_t_depth(i,j,iblock,n) + &
                           delt * TRACER(i,j,lev,1,curtime,iblock)
                      sbuff_sum_s_depth(i,j,iblock,n) = sbuff_sum_s_depth(i,j,iblock,n) + &
                           delt * TRACER(i,j,lev,2,curtime,iblock)
                   end if
                end do
             end do
          end do
       end do
    end if

    if ( State_FldChk(exportState, 'Faoo_fco2_ocn') .and. sflux_co2_nf_ind > 0) then
       do iblock = 1, nblocks_clinic
          call named_field_get(sflux_co2_nf_ind, iblock, work(:,:,iblock))
          sbuff_sum_co2(:,:,iblock) = sbuff_sum_co2(:,:,iblock) + delt_last*work(:,:,iblock)
       enddo
    end if

    first = .false.

  end subroutine pop_sum_buffer

end module ocn_import_export
