module ocn_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use POP_KindsMod,          only: POP_i4, POP_r8
  use POP_ErrorMod,          only: POP_ErrorSet, POP_Success
  use POP_FieldMod,          only: POP_fieldKindScalar
  use POP_GridHorzMod,       only: POP_gridHorzLocCenter
  use POP_HaloMod,           only: POP_HaloUpdate
  use med_constants_mod,     only: CS,CL
  use kinds_mod,             only: int_kind, r8
  use shr_cal_mod,           only: shr_cal_date2ymd
  use shr_sys_mod,           only: shr_sys_flush, shr_sys_abort
  use communicate,           only: my_task, master_task
  use domain,                only: distrb_clinic, POP_haloClinic
  use forcing_shf,           only: SHF_QSW
  use forcing_sfwf,          only: lsend_precip_fact, precip_fact
  use forcing_fields,        only: EVAP_F, PREC_F, SNOW_F, MELT_F, ROFF_F, IOFF_F
  use forcing_fields,        only: SALT_F
  use forcing_fields,        only: SENH_F, LWUP_F, LWDN_F, MELTH_F
  use forcing_fields,        only: ATM_CO2_PROG_nf_ind, ATM_CO2_DIAG_nf_ind
  use forcing_fields,        only: ATM_NHx_nf_ind, ATM_NOy_nf_ind
  use forcing_fields,        only: IFRAC, U10_SQR, ATM_PRESS
  use forcing_fields,        only: LAMULT, USTOKES, VSTOKES
  use forcing_fields,        only: ATM_FINE_DUST_FLUX, ATM_COARSE_DUST_FLUX, SEAICE_DUST_FLUX
  use forcing_fields,        only: ATM_BLACK_CARBON_FLUX, SEAICE_BLACK_CARBON_FLUX
  use mcog,                  only: lmcog, mcog_ncols, import_mcog
  use forcing_coupled,       only: update_ghost_cells_coupler_fluxes, rotate_wind_stress
  use ice,                   only: QFLUX, QICE, AQICE, tlast_ice
  use global_reductions,     only: global_sum_prod
  use io_tools,              only: document
  use named_field_mod,       only: named_field_register, named_field_get_index, named_field_set, named_field_get
  use vmix_kpp,              only: KPP_HBLT      ! ocn -> wav, bounadry layer depth
  use shr_nuopc_scalars_mod, only: flds_scalar_name
  use shr_nuopc_scalars_mod, only: flds_scalar_num
  use shr_nuopc_methods_mod, only: shr_nuopc_methods_chkerr
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

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr

  ! accumulated sum of send buffer quantities for averaging before being sent
  real (r8) :: sbuff_sum_u    (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_v    (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_t    (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_s    (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_dhdx (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_dhdy (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_bld  (nx_block,ny_block,max_blocks_clinic)
  real (r8) :: sbuff_sum_co2  (nx_block,ny_block,max_blocks_clinic)

  ! tlast_coupled is incremented by delt every time pop_sum_buffer is called
  ! tlast_coupled is reset to 0 when ocn_export is called
  real (r8) :: tlast_coupled

  integer     , parameter :: dbug = 1        ! i/o debug messages
  character(*), parameter :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ocn_advertise_fields(gcomp, importState, exportState, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    integer, intent(out) :: rc

    ! local variables
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    logical       :: flds_i2o_per_cat  ! .true. => select per ocn thickness category
    integer       :: dbrc
    character(len=*), parameter :: subname='(ocn_import_export:ocn_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_i2o_per_cat
    call ESMF_LogWrite('flds_i2o_per_cat = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

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

    ! from atmosphere
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_pslv')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2prog')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_co2diag')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_lwdn')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_snow')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_rain')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_bcphidry')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_bcphodry')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_bcphiwet')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_ocphidry')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_ocphodry')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_ocphiwet')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry1')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry2')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry3')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstdry4')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet1')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet2')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet3')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_dstwet4')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_nhx')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_noy')

    ! optional per thickness category fields
    ! this implementation only handles columns due to ice thickness categories

    if (flds_i2o_per_cat) then
       mcog_ncols = 6 ! TODO: hard-wire ice_ncat to 5 for now
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sf_afrac')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sf_afracr')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet_afracr')
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_swpen_ifrac_n', ungridded_lbound=1, ungridded_ubound=mcog_ncols-1)
       call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac_n'        , ungridded_lbound=1, ungridded_ubound=mcog_ncols-1)
    else
       mcog_ncols = 1
    endif

    do n = 1,fldsToOcn_num
       call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

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

    do n = 1,fldsFrOcn_num
       call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ocn_advertise_fields

!==============================================================================

  subroutine ocn_realize_fields(gcomp, mesh, rc)

    ! input/output variables
    type(ESMF_GridComp)           :: gcomp
    type(ESMF_Mesh) , intent(in)  :: mesh
    integer         , intent(out) :: rc

    ! local variables
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    character(len=*), parameter :: subname='(ocn_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Only mesh is supported for now

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrOcn, &
         numflds=fldsFrOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':POP_Export',&
         mesh=mesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToOcn, &
         numflds=fldsToOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':POP_Import',&
         mesh=mesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ocn_realize_fields

  !==============================================================================

  subroutine ocn_import( importState, ldiag_cpl, errorCode, rc )

    !-----------------------------------------------------------------------
    ! swnet  -- net short-wave heat flux                 (W/m2   )
    ! lwup   -- longwave radiation (up)                  (W/m2   )
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)   , intent(in)  :: importState
    logical (log_kind) , intent(in)  :: ldiag_cpl
    integer            , intent(out) :: errorCode
    integer            , intent(out) :: rc

    ! local variables
    type (block)         :: this_block ! local block info
    character (char_len) :: label,  message
    real (r8)            :: work1(nx_block,ny_block,max_blocks_clinic)
    real (r8)            :: work2(nx_block,ny_block,max_blocks_clinic)  ! local work space
    integer (int_kind)   :: i,j,k,n,ncol,iblock,nfld
    real (r8)            :: frac_col_1pt(mcog_ncols)
    real (r8)            :: fracr_col_1pt(mcog_ncols)
    real (r8)            :: qsw_fracr_col_1pt(mcog_ncols)
    real (r8)            :: m2percm2, gsum
    real (r8), pointer   :: Foxx_swnet(:)
    real (r8), pointer   :: Foxx_swnet_afracr(:)
    real (r8), pointer   :: Sf_afrac(:)
    real (r8), pointer   :: Sf_afracr(:)
    real (r8), pointer   :: Si_ifrac_n(:,:)
    real (r8), pointer   :: Fioi_swpen_ifrac(:,:)
    integer (int_kind)   :: fieldCount
    integer              :: dbrc
    character (char_len), allocatable :: fieldNameList(:)
    type(ESMF_StateItem_Flag) :: itemflag
    character(len=*), parameter :: subname='(ocn_import_export:ocn_import)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

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

    ! zonal wind stress  (W/m2)
    call state_getimport(importState, 'Foxx_taux', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! meridonal wind stress (W/m2)
    call state_getimport(importState, 'Foxx_tauy', work2, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! rotate true zonal/meridional wind stress into local coordinates,
    ! convert to dyne/cm**2, and shift SMFT to U grid
    ! halo updates are performed in subroutine rotate_wind_stress,
    ! following the rotation

    call rotate_wind_stress(work1, work2)

    ! evaporation flux (kg/m2/s)
    call state_getimport(importState, 'Foxx_evap', EVAP_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! sensible heat flux (W/m2)
    call state_getimport(importState, 'Foxx_sen', SENH_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Foxx_lwup', LWUP_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Foxx_swnet', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    SHF_QSW(:,:,:) = work1(:,:,:) * RCALCT(:,:,:)*hflux_factor  !  convert from W/m**2

    ! 10m wind speed squared (m^2/s^2)
    call state_getimport(importState, 'So_duu10n', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    U10_SQR(:,:,:)  = cmperm * cmperm * work1(:,:,:) * RCALCT(:,:,:) ! convert from m**2/s**2 to cm**2/s**2

    !-----------------------------------------------------------------------
    ! from atmosphere
    !-----------------------------------------------------------------------

    ! sea-level pressure (Pa)
    call state_getimport(importState, 'Sa_pslv', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ATM_PRESS(:,:,:) = c10 * work1(:,:,:) * RCALCT(:,:,:) ! convert from Pa to dynes/cm**2

    ! water flux due to snow (kg/m2/s)
    call state_getimport(importState, 'Faxa_snow', SNOW_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! water flux due to rain (kg/m2/s)
    call state_getimport(importState, 'Faxa_rain', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    PREC_F(:,:,:) = work1(:,:,:) + SNOW_F(:,:,:) ! rain + snow

    ! longwave radiation (down)   (W/m2)
    call state_getimport(importState, 'Faxa_lwdn', LWDN_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Faxa_dstwet1', output=work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry1', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ATM_FINE_DUST_FLUX(:,:,:) = 0.1_r8 * RCALCT(:,:,:) * work1(:,:,:) ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s)

    call state_getimport(importState, 'Faxa_dstwet2', output=work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry2', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstwet3', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry3', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstwet4', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Faxa_dstdry4', output=work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ATM_COARSE_DUST_FLUX(:,:,:) = 0.1_r8 * RCALCT(:,:,:) * work1(:,:,:) ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s)

    !-----------------------------------------------------------------------
    ! from sea-ice
    !-----------------------------------------------------------------------

    ! ice fraction
    call state_getimport(importState, 'Si_ifrac', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    IFRAC(:,:,:) = work1(:,:,:) * RCALCT(:,:,:)

    ! snow melt flux (kg/m2/s)
    call state_getimport(importState, 'Fioi_meltw', MELT_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! heat flux from snow & ice melt (W/m2)
    call state_getimport(importState, 'Fioi_melth', MELTH_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! salt from ice (kg(salt)/m2/s)
    call state_getimport(importState, 'Fioi_salt', SALT_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Fioi_flxdst', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    SEAICE_DUST_FLUX(:,:,:) = 0.1_r8 * RCALCT(:,:,:) * work1(:,:,:) ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s)

    call state_getimport(importState, 'Fioi_bcpho', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(importState, 'Fioi_bcphi', work1, do_sum=.true., rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    SEAICE_BLACK_CARBON_FLUX(:,:,:) = 0.1_r8 * RCALCT(:,:,:) * work1(:,:,:) ! convert from MKS (kg/m^2/s) to CGS (g/cm^2/s)

    !  optional fields per mcog column

    call state_getfldptr(importState, 'Foxx_swnet', Foxx_swnet, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (lmcog) then

       if ( State_FldChk(importState, 'Si_ifrac_n')        .and. &
            State_FldChk(importState, 'Fioi_swpen_ifrac')  .and. &
            State_FldChk(importState, 'Foxx_swnet_afracr') .and. &
            State_FldChk(importState, 'Sf_afrac')          .and. &
            State_FldChk(importState, 'Sf_afracr')) then

          call state_getfldptr(importState, 'Si_ifrac_n', Si_ifrac_n, rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call state_getfldptr(importState, 'Fioi_swpen_ifrac', Fioi_swpen_ifrac, rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call state_getfldptr(importState, 'Foxx_swnet_afracr', Foxx_swnet_afracr, rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call state_getfldptr(importState, 'Sf_afrac', Sf_afrac, rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call state_getfldptr(importState, 'Sf_afracr', Sf_afracr, rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          n = 0
          do iblock = 1, nblocks_clinic
             this_block = get_block(blocks_clinic(iblock),iblock)
             do j=this_block%jb,this_block%je
                do i=this_block%ib,this_block%ie
                   n = n + 1
                   ! extract fields for each column and pass to import_mcog
                   do ncol = 1, mcog_ncols
                      if (ncol == 1) then
                         frac_col_1pt(ncol)  = max(c0, min(c1, Sf_afrac(n)))
                         fracr_col_1pt(ncol) = max(c0, min(c1, Sf_afracr(n)))
                         qsw_fracr_col_1pt(ncol) = Foxx_swnet_afracr(n)
                      else
                         frac_col_1pt(ncol)  = max(c0, min(c1, Si_ifrac_n(n,ncol)))
                         fracr_col_1pt(ncol) = max(c0, min(c1, Si_ifrac_n(n,ncol)))
                         qsw_fracr_col_1pt(ncol) = Fioi_swpen_ifrac(n,ncol)
                      end if
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
                enddo ! do i
             enddo ! do j
          enddo ! do iblock = 1, nblocks_clinic

       end if

    endif ! if (lmcog) then

    !-----------------------------------------------------------------------
    ! from wave
    !-----------------------------------------------------------------------

    call state_getimport(importState, 'Sw_lamult', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    where (IFRAC <= 0.05_r8)
       LAMULT(:,:,:) = work1 * RCALCT(:,:,:) ! import enhancement factor (unitless)
    elsewhere
       LAMULT(:,:,:) = c1
    end where

    call state_getimport(importState, 'Sw_ustokes', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    USTOKES(:,:,:) = work1(:,:,:) * RCALCT(:,:,:) ! Stokes drift (m/s)

    call state_getimport(importState, 'Sw_vstokes', work1, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    VSTOKES(:,:,:) = work1(:,:,:) * RCALCT(:,:,:) ! Stokes drift (m/s)

    !-----------------------------------------------------------------------
    ! from river
    !-----------------------------------------------------------------------

    ! liquid runoff flux (kg/m2/s)
    call state_getimport(importState, 'Foxx_rofl', ROFF_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ice runoff flux (kg/m2/s)
    call state_getimport(importState, 'Foxx_rofi', IOFF_F, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call POP_HaloUpdate(work1,POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating DIAG CO2 halo')
          return
       endif
       call named_field_set(ATM_CO2_DIAG_nf_ind, work1)
    endif

    call ESMF_StateGet(importState, 'Faxa_nhx', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getimport(importState, 'Faxa_nhx', work1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Note - the input units are kgN/m2/s to nmolN/cm2/s
       ! TODO: Keith has pointed out might want to use 14.007_r8 instead of 14.0_r8 for more
       ! consistency when bringing in N isotopes into the code
       work1(:,:,:) = work1(:,:,:) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8)
       call POP_HaloUpdate(work1,POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
            errorCode, fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call POP_ErrorSet(errorCode, 'ocn_import_import: error updating DIAG NHx halo')
          return
       endif
       call named_field_set(ATM_NHx_nf_ind, work1)
    endif

    call ESMF_StateGet(importState, 'Faxa_noy', itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

       call state_getimport(importState, 'Faxa_noy', work1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Note - the input units are kgN/m2/s to nmolN/cm2/s
       ! TODO: Keith has pointed out might want to use 14.007_r8 instead of 14.0_r8 for more
       ! consistency when bringing in N isotopes into the code
       work1(:,:,:) = work1(:,:,:) * (1.0e-1_r8 * (c1/14.0_r8) * 1.0e9_r8)
       call POP_HaloUpdate(work1,POP_haloClinic, POP_gridHorzLocCenter, POP_fieldKindScalar, &
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
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(importState, itemNameList=fieldNameList, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in import state
       m2percm2  = mpercm*mpercm
       do nfld = 1, fieldCount
          if (trim(fieldNameList(nfld)) /= flds_scalar_name) then
             call state_getimport(importState, trim(fieldNameList(nfld)), work1, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             gsum = global_sum_prod(work1, TAREA, distrb_clinic,  field_loc_center, RCALCT) * m2percm2
             
             if (my_task == master_task) then
                write(stdout,1100)'ocn','recv', trim(fieldNameList(nfld)), gsum
                call shr_sys_flush(stdout)
             endif
          end if
       end do

    end if
1100 format ('comm_diag ', a3, 1x, a4, 1x, a15, 1x, es26.19:, 1x, a6)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ocn_import

!==============================================================================

  subroutine ocn_export(exportState, ldiag_cpl, errorCode, rc)

    !-----------------------------------------------------------------------
    ! Create export state
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                 :: exportState
    logical (log_kind) , intent(in)  :: ldiag_cpl
    integer (POP_i4)   , intent(out) :: errorCode  ! pop error code
    integer            , intent(out) :: rc         ! returned error code

    ! local variables
    type (block)         :: this_block ! local block info
    integer (int_kind)   :: n, i,j,k,iblock,nfld
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
    integer (int_kind)   :: fieldCount
    integer              :: dbrc
    character (char_len), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_export)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !-----------------------------------------------------------------------
    ! ocean mask
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_omask', dataPtr1, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    n=0
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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_v', dataPtr2, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

    !-----------------------------------------------------------------------
    ! convert and pack salinity
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_s', dataptr1, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

    !-----------------------------------------------------------------------
    ! convert and pack boundary layer depth
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_bldepth', dataptr1, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'So_dhdy', dataPtr2, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n + 1
             dataptr1(n) = QFLUX(i,j,iblock)
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
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock),iblock)
          do j=this_block%jb,this_block%je
             do i=this_block%ib,this_block%ie
                n = n + 1
                dataptr1(n) = sbuff_sum_co2(i,j,iblock)/tlast_coupled
             enddo
          enddo
       enddo
    endif

    !-----------------------------------------------------------------------
    ! diagnostics
    !-----------------------------------------------------------------------

    if (ldiag_cpl) then
       call ccsm_char_date_and_time
       write(stdout,*)'pop_send_to_coupler'

       ! Determine all field names in export state
       call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(fieldNameList(fieldCount))
       call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! loop over all fields in export state
       m2percm2  = mpercm*mpercm
       do nfld = 1, fieldCount
          if (trim(fieldNameList(nfld)) /= flds_scalar_name) then
             call state_getfldptr(exportState, trim(fieldNameList(nfld)), dataptr1, rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ocn_export

!==============================================================================

  subroutine state_getimport(state, fldname, output, do_sum, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)    :: state
    character(len=*)  , intent(in)    :: fldname
    real (r8)         , intent(inout) :: output(:,:,:)
    logical, optional , intent(in)    :: do_sum
    integer           , intent(out)   :: rc

    ! local variables
    type(block)       :: this_block         ! block information for current block
    integer           :: i, j, iblock, n   ! incides
    real(r8), pointer :: dataPtr1d(:)
    character(len=*), parameter :: subname='(ice_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call state_getfldptr(state, trim(fldname), dataptr1d, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine output array
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j = this_block%jb,this_block%je
          do i = this_block%ib,this_block%ie
             n = n + 1
             if (present(do_sum)) then
                output(i,j,iblock)  = output(i,j,iblock) + dataPtr1d(n)
             else
                output(i,j,iblock)  = dataPtr1d(n)
             end if
          end do
       end do
    end do

  end subroutine state_getimport

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
    use ESMF , only : ESMF_VM

    use shr_nuopc_methods_mod, only : shr_nuopc_methods_State_SetScalar

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
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                                         ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                                         ungriddedUbound=(/fldlist(n)%ungridded_ubound/), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
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
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    real (r8)                :: work(nx_block,ny_block,max_blocks_clinic) ! local work arrays
    real (r8)                :: delt                                      ! time interval since last step
    real (r8)                :: delt_last                                 ! time interval for previous step
    integer (int_kind)       :: iblock                                    ! block index
    integer (int_kind)       :: sflux_co2_nf_ind = 0                      ! named field index of fco2
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

    if ( State_FldChk(exportState, 'Faoo_fco2_ocn') .and. sflux_co2_nf_ind > 0) then
       do iblock = 1, nblocks_clinic
          call named_field_get(sflux_co2_nf_ind, iblock, work(:,:,iblock))
          sbuff_sum_co2(:,:,iblock) = sbuff_sum_co2(:,:,iblock) + delt_last*work(:,:,iblock)
       enddo
    end if

    first = .false.

  end subroutine pop_sum_buffer

end module ocn_import_export
