module ocn_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for POP
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_CompSetClock
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_CheckImport    => label_CheckImport
  use NUOPC_Model           , only : model_label_SetClock       => label_SetClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use constants             , only : c0, blank_fmt, ndelim_fmt
  use POP_IOUnitsMod        , only : POP_IOUnitsFlush, POP_stdout, inst_suffix, inst_index, inst_name
  use POP_ErrorMod          , only : POP_ErrorSet, POP_Success, POP_ErrorPrint
  use shr_file_mod          , only : shr_file_getLogUnit, shr_file_setLogUnit
  use shr_cal_mod           , only : shr_cal_date2ymd, shr_cal_ymd2date
  use shr_sys_mod           , only : shr_sys_abort
  use shr_kind_mod          , only : cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_orb_mod           , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use POP_KindsMod          , only : POP_i4
  use kinds_mod             , only : int_kind, log_kind, char_len, r8
  use ocn_communicator      , only : mpi_communicator_ocn
  use communicate           , only : my_task, master_task, exit_message_environment
  use blocks                , only : block, get_block, nblocks_x, nblocks_y
  use exit_mod              , only : exit_POP, sigAbort, flushm
  use POP_DistributionMod   , only : POP_DistributionGet
  use domain                , only : POP_distrbclinic, distrb_clinic, POP_haloClinic, nblocks_clinic, blocks_clinic
  use domain_size           , only : nx_global, ny_global
  use forcing_fields        , only : ATM_CO2_PROG_nf_ind, ATM_CO2_DIAG_nf_ind, ATM_NHx_nf_ind, ATM_NOy_nf_ind, lsmft_avail
  use forcing_shf           , only : SHF_QSW
  use forcing_sfwf          , only : lsend_precip_fact, precip_fact
  use forcing_coupled       , only : ncouple_per_day, pop_set_coupled_forcing
  use forcing_coupled       , only : orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp
  use grid                  , only : TLAT, TLON, TAREA, KMT
  use io_tools              , only : document
  use io_types              , only : stdout
  use named_field_mod       , only : named_field_register, named_field_get_index, named_field_set, named_field_get
  use timers                , only : get_timer, timer_start, timer_stop, timer_print_all
  use diagnostics           , only : check_KE
  use output                , only : output_driver
  use step_mod              , only : step
  use time_management       , only : iyear0, imonth0, iday0, runid, init_time_flag, nsteps_run
  use time_management       , only : iyear, imonth, iday, ihour, iminute, isecond, document_time_flags
  use time_management       , only : cyear, cmonth, cday, seconds_in_day, seconds_in_hour, seconds_in_minute
  use time_management       , only : check_time_flag, override_time_flag, ccsm_char_date_and_time
  use time_management       , only : access_time_flag
  use registry              , only : registry_match, register_string
  use ecosys_forcing_mod    , only : ldriver_has_ndep, ldriver_has_atm_co2_diag, ldriver_has_atm_co2_prog
  use initial               , only : pop_init_phase1, pop_init_phase2
  use POP_MCT_vars_mod      , only : pop_mct_init
  use perf_mod              , only : t_startf, t_stopf
  use ocn_import_export     , only : ocn_advertise_fields, ocn_realize_fields
  use ocn_import_export     , only : ocn_import, ocn_export, pop_sum_buffer, tlast_coupled
  use ocn_shr_methods       , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use ocn_shr_methods       , only : set_component_logging, get_component_instance, log_clock_advance

  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private                              ! By default make data private

  public  :: SetServices

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  character(len=CL)   :: flds_scalar_name = ''
  integer             :: flds_scalar_num = 0
  integer             :: flds_scalar_index_nx = 0
  integer             :: flds_scalar_index_ny = 0
  integer             :: flds_scalar_index_precip_factor = 0

  integer (int_kind)  :: stop_now    ! flag id for stop_now flag
  integer (int_kind)  :: cpl_ts      ! flag id for coupled timestep flag
  integer (int_kind)  :: timer_total ! timer for entire run phase

  integer, parameter  :: dbug = 1
  logical (log_kind)  :: ldiag_cpl = .false.
  integer (int_kind)  :: cpl_write_restart   ! flag id for write restart
  integer (int_kind)  :: cpl_write_history   ! flag id for write history
  integer (int_kind)  :: cpl_write_tavg      ! flag id for write tavg
  integer (int_kind)  :: cpl_diag_global     ! flag id for computing diagnostics
  integer (int_kind)  :: cpl_diag_transp     ! flag id for computing diagnostics
  character(char_len) :: runtype

  character(*), parameter :: u_FILE_u = &
       __FILE__

  character(len=CL)      :: attribute_orb_mode        ! attribute - orbital mode
  integer                :: attribute_orb_iyear       ! attribute - orbital year
  integer                :: attribute_orb_iyear_align ! attribute - associated with model year
  real(R8)               :: attribute_orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: attribute_orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: attribute_orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

!=======================================================================
contains
!=======================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    character(len=*),parameter  :: subname='ocn_comp_nuopc:(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
         specRoutine=ModelCheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    use NUOPC, only : NUOPC_isConnected

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: logmsg
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname='ocn_comp_nuopc:(InitializeAdvertise) '
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxPrecipFactor", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_precip_factor
       write(logmsg,*) flds_scalar_index_precip_factor
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_precip_factor = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxPrecipFactor')
    endif

    call ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    !-----------------------------------------------------------------------
    !  first initializaiton phase of pop2
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------

    use shr_const_mod      , only: shr_const_pi  
    use constants          , only: radius

    ! Initialize POP

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !  local variables
    type (block)            :: this_block          ! block information for current block
    integer                 :: iblock
    type(ESMF_VM)           :: vm
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_Array)        :: EArray
    integer                 :: spatialDim
    integer                 :: numOwnedElements
    real(R8), pointer       :: ownedElemCoords(:)
    real(R8)                :: lon, lat, area
    real(R8)                :: diff_lon, diff_lat, diff_area
    real(R8), pointer       :: areaMesh(:)  
    real(R8), allocatable   :: lonMesh(:), latMesh(:)
    integer , allocatable   :: gindex_ocn(:)
    integer , allocatable   :: gindex_elim(:)
    integer , allocatable   :: gindex(:)
    integer                 :: globalID
    character(CL)           :: cvalue
    integer                 :: num_elim_global
    integer                 :: num_elim_local
    integer                 :: num_elim
    integer                 :: num_ocn
    integer                 :: num_elim_gcells ! local number of eliminated gridcells
    integer                 :: num_elim_blocks ! local number of eliminated blocks
    integer                 :: num_total_blocks
    integer                 :: my_elim_start
    integer                 :: my_elim_end
    integer(int_kind)       :: lsize
    integer(int_kind)       :: shrlogunit      ! old values
    integer(int_kind)       :: nThreads
    integer(int_kind)       :: npes
    integer(int_kind)       :: iam
    character(len=32)       :: starttype
    integer                 :: n,i,j,iblk,jblk,ig,jg
    integer                 :: lbnum
    integer(POP_i4)         :: errorCode       ! error code
    integer                 :: lmpicom
    integer(POP_i4) , pointer    :: blockLocation(:)
    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(InitializeRealize)"
#ifdef _OPENMP
    integer, external :: omp_get_max_threads  ! max threads that can execute concurrently in a single parallel region
#endif
    type(ESMF_Field)  :: areaField
    type(ESMF_Array)  :: areaArray
    real(R8), pointer :: areaptr(:)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = POP_Success

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

#ifdef _OPENMP
    nThreads = omp_get_max_threads()
#endif

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=iam, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    mpi_communicator_ocn = lmpicom

    ! reset shr logging to my log file
    if (iam == 0) then
       call set_component_logging(gcomp, .true., stdout, shrlogunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call set_component_logging(gcomp, .false., stdout, shrlogunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

#if (defined _MEMTRACE)
    if (iam == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','InitializeRealize:start::',lbnum)
    endif
#endif

    !-----------------------------------------------------------------------
    ! initialize the model run
    !-----------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) runid

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (trim(starttype) == trim('startup')) then
       runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
       runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
       runtype = "continue"
    else
       write(stdout,*) 'ocn_comp_nuopc: ERROR: unknown starttype'
       call exit_POP(sigAbort,' ocn_comp_nuopc: ERROR: unknown starttype')
    end if

    ! TODO: Set model_doi_url

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    inst_name = "OCN"//trim(inst_suffix)
    write(6,*)'DEBUG: inst_index,inst_suffix,inst_name = ',inst_index,trim(inst_suffix),trim(inst_name)

    !-----------------------------------------------------------------------
    !  first initializaiton phase of pop2
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------

    call t_startf ('pop_init1')

    call pop_init_phase1(errorCode, nuopc_cap=.true.)

    if (errorCode /= POP_Success) then
       call ESMF_LogWrite(trim(subname)//'POP_Initialize1: error in pop_init_phase1',ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    endif

    call t_stopf ('pop_init1')

    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

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
    allocate(gindex_ocn(lsize))

    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n=n+1
             gindex_ocn(n) = (this_block%j_glob(j)-1)*(nx_global) + this_block%i_glob(i)
          enddo
       enddo
    enddo

    ! Determine total number of eliminated blocks globally
    globalID = 0
    num_elim_global = 0  ! number of eliminated blocks
    do jblk = 1,nblocks_y
       do iblk = 1,nblocks_x
          globalID = globalID + 1
          call POP_DistributionGet(POP_distrbclinic, errorCode, blockLocation=blockLocation)
          if (blockLocation(globalID) == 0) then
             num_elim_global = num_elim_global + 1
          end if
       end do
    end do

    if (num_elim_global > 0) then

       ! Distribute the eliminated blocks in a round robin fashion amoung processors
       num_elim_local = num_elim_global / npes
       my_elim_start = num_elim_local*iam + min(iam, mod(num_elim_global, npes)) + 1
       if (iam < mod(num_elim_global, npes)) then
          num_elim_local = num_elim_local + 1
       end if
       my_elim_end = my_elim_start + num_elim_local - 1

       ! Determine the number of eliminated gridcells locally
       globalID = 0
       num_elim_blocks = 0  ! local number of eliminated blocks
       num_elim_gcells = 0
       do jblk = 1,nblocks_y
          do iblk= 1,nblocks_x
             globalID = globalID + 1
             call POP_DistributionGet(POP_distrbclinic, errorCode, blockLocation=blockLocation)
             if (blockLocation(globalID) == 0) then
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   this_block = get_block(globalID, globalID)
                   num_elim_gcells = num_elim_gcells + &
                        (this_block%je-this_block%jb+1) * (this_block%ie-this_block%ib+1)
                end if
             end if
          end do
       end do
       allocate(gindex_elim(num_elim_gcells))

       ! Determine the global index space of the eliminated gridcells
       globalID = 0
       num_elim_gcells = 0  ! local number of eliminated gridcells
       num_elim_blocks = 0  ! local number of eliminated blocks
       do jblk = 1,nblocks_y
          do iblk = 1,nblocks_x
             globalID = globalID + 1
             call POP_DistributionGet(POP_distrbclinic, errorCode, blockLocation=blockLocation)
             if (blockLocation(globalID) == 0) then
                this_block = get_block(globalID, globalID)
                num_elim_blocks = num_elim_blocks + 1
                if (num_elim_blocks >= my_elim_start .and. num_elim_blocks <= my_elim_end) then
                   do j = this_block%jb,this_block%je
                      do i = this_block%ib,this_block%ie
                         num_elim_gcells = num_elim_gcells + 1
                         ig = this_block%i_glob(i)
                         jg = this_block%j_glob(j)
                         gindex_elim(num_elim_gcells) = (jg-1)*nx_global + ig
                      end do
                   end do
                end if
             end if
          end do
       end do

       ! create a global index that includes both active and eliminated gridcells
       num_ocn  = size(gindex_ocn)
       num_elim = size(gindex_elim)
       allocate(gindex(num_elim + num_ocn))
       do n = 1,num_ocn
          gindex(n) = gindex_ocn(n)
       end do
       do n = num_ocn+1,num_ocn+num_elim
          gindex(n) = gindex_elim(n-num_ocn)
       end do

       deallocate(gindex_elim)

    else

       ! No eliminated land blocks
       num_ocn = size(gindex_ocn)
       allocate(gindex(num_ocn))
       do n = 1,num_ocn
          gindex(n) = gindex_ocn(n)
       end do

    end if

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------

    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the POP mesh
    !---------------------------------------------------------------------------

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(stdout,*)'mesh file for pop domain is ',trim(cvalue)
    end if

    !---------------------------------------------------------------------------
    ! Error checking
    !---------------------------------------------------------------------------

    ! obtain mesh lats and lons
    call ESMF_MeshGet(Emesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(lonMesh(numOwnedElements), latMesh(numOwnedElements))
    call ESMF_MeshGet(Emesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do

    ! obtain mesh area
    ! areaField = ESMF_FieldCreate(Emesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call ESMF_FieldGet(areaField, array=areaArray, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call ESMF_MeshGet(Emesh, elemAreaArray=areaArray, rc=rc)  !<== crashes if this is added
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call ESMF_ArrayGet(areaArray, farrayptr=areaMesh, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    ! error check differences between internally generated lons and those read in
    n = 0
    do iblock = 1, nblocks_clinic
       this_block = get_block(blocks_clinic(iblock),iblock)
       do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
             n = n+1
             lon  = TLON(i,j,iblock) * 180._R8/shr_const_pi
             lat  = TLAT(i,j,iblock) * 180._R8/shr_const_pi
             area = TAREA(i,j,iblock) / (radius*radius) 
             diff_lon  = abs(lonMesh(n) - lon)
             diff_lat  = abs(latMesh(n) - lat)
             !diff_area = abs(areaMesh(n) - area)
             if (diff_lon  > 1.e-10 ) write(6,100) n,lonMesh(n),lon,diff_lon
             if (diff_lat  > 1.e-10 ) write(6,101) n,latMesh(n),lat,diff_lat
             !if (diff_area > 1.e-12) write(6,102) n,areaMesh(n),area,diff_area
          end do
       end do
    end do
100 format('WARNING: POP n, lonmesh(n) , lon , diff_lon  = ',i6,2(f21.13,3x),d21.5)
101 format('WARNING: POP n, latmesh(n) , lat , diff_lat  = ',i6,2(f21.13,3x),d21.5)
102 format('WARNING: POP n, areamesh(n), area, diff_area = ',i6,2(f21.13,3x),d21.5)

    ! deallocate memory
    deallocate(ownedElemCoords, lonMesh, latMesh)

    !call ESMF_FieldDestroy(areaField)

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ocn_realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Initialize MCT gsmaps and domains
    !-----------------------------------------------------------------

    ! call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) ocnid  ! convert from string to integer

    ! call pop_mct_init(ocnid, mpi_communicator_ocn)
    ! if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine DataInitialize(gcomp, rc)

    !-----------------------------------------------------------------------
    !  second initializaiton phase of pop2
    !  - initialize passive tracer modules -- after call init_forcing_coupled
    !  - initialize vertical mixing variables
    !  - initialize niw driven mixing
    !  - initialize geothermal heat flux
    !  - initialize horizontal mixing variables
    !  - initialize overflow regional values
    !  - initialize advection variables
    !  - initialize shortwave absorption
    !  - initialize estuary parameterization
    !  - partial coupling forcing initialization
    !  - initialize overflows output diagnostics filename
    !  - initialize output; subroutine init_output calls
    !  - initialize global budget diagnostics
    !  - initialize step timers
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)          :: clock
    type(ESMF_State)          :: importState
    type(ESMF_State)          :: exportState
    type(ESMF_StateItem_Flag) :: itemType
    type(ESMF_StateItem_Flag) :: itemType1
    type(ESMF_StateItem_Flag) :: itemType2
    type(ESMF_TimeInterval)   :: timeStep        ! Model timestep
    type(ESMF_Time)           :: starttime
    character(CL)             :: cvalue
    integer(int_kind)         :: ocn_cpl_dt
    integer(int_kind)         :: pop_cpl_dt
    integer(int_kind)         :: start_ymd
    integer(int_kind)         :: start_tod
    integer(int_kind)         :: start_year
    integer(int_kind)         :: start_day
    integer(int_kind)         :: start_month
    integer(int_kind)         :: start_hour
    integer(POP_i4)           :: errorCode       ! error code
    integer(int_kind)         :: shrlogunit      ! old values
    integer                   :: ocnid
    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(DataInitialize)"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (stdout)

    !-----------------------------------------------------------------------
    ! register non-standard incoming fields
    !-----------------------------------------------------------------------

    ! query the Component for its importState, exportState and clock
    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(importState, 'Sa_co2prog', itemType, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ldriver_has_atm_co2_prog = (itemType /= ESMF_STATEITEM_NOTFOUND)

    call ESMF_StateGet(importState, 'Sa_co2diag', itemType, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ldriver_has_atm_co2_diag = (itemType /= ESMF_STATEITEM_NOTFOUND)

    call ESMF_StateGet(importState, 'Faxa_nhx', itemType1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(importState, 'Faxa_noy', itemType2, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ldriver_has_ndep = ((itemType1 /= ESMF_STATEITEM_NOTFOUND) .or. (itemType2 /= ESMF_STATEITEM_NOTFOUND))

    if (ldriver_has_atm_co2_prog) then
       call named_field_register('ATM_CO2_PROG', ATM_CO2_PROG_nf_ind)
    endif
    if (ldriver_has_atm_co2_diag) then
       call named_field_register('ATM_CO2_DIAG', ATM_CO2_DIAG_nf_ind)
    endif
    if (ldriver_has_ndep) then
       if (my_task == master_task) write(stdout,'(" using ATM_NHx and ATM_NOy from mediator")')
       call named_field_register('ATM_NHx', ATM_NHx_nf_ind)
       call named_field_register('ATM_NOy', ATM_NOy_nf_ind)
    endif

    call register_string('pop_init_coupled')
    call flushm (stdout)

    !-----------------------------------------------------------------------
    ! second initialization phase of pop2
    !-----------------------------------------------------------------------

    call t_startf ('pop_init2')

    call pop_init_phase2(errorCode)
    if (errorCode /= POP_Success) then
       call POP_ErrorSet(errorCode, 'POP_Initialize2: error in pop_init_phase2')
       return
    endif

    if (errorCode /= POP_Success) then
       call ESMF_LogWrite(trim(subname)//'POP_Initialize2: error in pop_init_phase2',ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    endif

    ! initialize driver-level flags and timers
    call access_time_flag ('stop_now', stop_now)
    call access_time_flag ('coupled_ts', cpl_ts)
    call get_timer(timer_total,'TOTAL', 1, distrb_clinic%nprocs)

    call t_stopf ('pop_init2')

    !-----------------------------------------------------------------------
    !  initialize time-stamp information
    !-----------------------------------------------------------------------

    call ccsm_char_date_and_time

    !-----------------------------------------------------------------------
    ! check for consistency of pop and sync clock initial time
    !-----------------------------------------------------------------------

    if (runtype == 'initial') then

       call ESMF_ClockGet( clock, startTime=startTime, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_TimeGet( startTime, yy=start_year, mm=start_month, dd=start_day, s=start_tod, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_cal_ymd2date(start_year,start_month,start_day,start_ymd)

       if (iyear0 /= start_year) then
          call document ('DataInitialize', 'iyear0      ', iyear0)
          call document ('DataInitialize', 'imonth0     ', imonth0)
          call document ('DataInitialize', 'iday0       ', iday0)
          call document ('DataInitialize', 'start_year  ', start_year)
          call document ('DataInitialize', 'start_month ', start_month)
          call document ('DataInitialize', 'start_day   ', start_day)

          ! skip exit_POP if pop is 1 day ahead across a year boundary
          if (.not. (iyear0 == start_year + 1 .and. &
                    (imonth0 == 1 .and. start_month == 12) .and. &
                    (iday0 == 1   .and. start_day == 31))) then
             call exit_POP(sigAbort,' iyear0 does not match start_year')
          endif

       else if (imonth0 /= start_month) then
          call document ('DataInitialize', 'imonth0     ', imonth0)
          call document ('DataInitialize', 'iday0       ', iday0)
          call document ('DataInitialize', 'start_month ', start_month)
          call document ('DataInitialize', 'start_day   ', start_day)
          ! skip exit_POP if pop is 1 day ahead across a month boundary
          !   this conditional doesn't confirm that start_day is the last day of the month,
          !   only that iday0 is the first day of the month
          if (.not. (imonth0 == start_month + 1 .and. iday0 == 1)) then
             call exit_POP(sigAbort,' imonth0 does not match start_month')
          endif

       else if (iday0 /= start_day) then
          call document ('DataInitialize', 'iday0       ', iday0)
          call document ('DataInitialize', 'start_day   ', start_day)

          ! skip exit_POP if pop is 1 day ahead
          if (.not. (iday0 == start_day + 1)) then
             call exit_POP(sigAbort,' iday0 does not match start_day')
          endif
       end if
    end if

    !-----------------------------------------------------------------
    ! Initialize MCT gsmaps and domains
    !-----------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocnid  ! convert from string to integer

    call pop_mct_init(ocnid, mpi_communicator_ocn)
    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    ! Initialize flags and shortwave absorption profile
    ! Note that these cpl_write_xxx flags have no freqency options
    ! set; therefore, they will retain a default value of .false.
    ! unless they are explicitly set .true.  at the appropriate times
    !-----------------------------------------------------------------------

    call init_time_flag('cpl_write_restart',cpl_write_restart, owner = 'DataInitialize')
    call init_time_flag('cpl_write_history',cpl_write_history, owner = 'DataInitialize')
    call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg,    owner = 'DataInitialize')
    call init_time_flag('cpl_diag_global'  ,cpl_diag_global,   owner = 'DataInitialize')
    call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp,   owner = 'DataInitialize')

    lsmft_avail = .true.
    tlast_coupled = c0

    !-----------------------------------------------------------------------
    ! initialize necessary coupling info
    !-----------------------------------------------------------------------

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet( timeStep, s=ocn_cpl_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    pop_cpl_dt = seconds_in_day / ncouple_per_day

    if (pop_cpl_dt /= ocn_cpl_dt) then
       write(stdout,*)'pop_cpl_dt= ',pop_cpl_dt,' ocn_cpl_dt= ',ocn_cpl_dt
       call exit_POP(sigAbort,'ERROR pop_cpl_dt and ocn_cpl_dt must be identical')
    end if

    !-----------------------------------------------------------------------
    ! send export state
    !-----------------------------------------------------------------------

    call pop_sum_buffer(exportState, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (errorCode /= POP_Success) then
       call POP_ErrorPrint(errorCode)
       call exit_POP(sigAbort, 'ERROR in ocn_export')
    endif

    call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( lsend_precip_fact ) then
       if ( flds_scalar_index_precip_factor > 0 ) then
          call State_SetScalar(precip_fact, flds_scalar_index_precip_factor, exportState, &
               flds_scalar_name, flds_scalar_num, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_sys_abort('flds_scalar_index_precip_factor must be > 0 when lsend_precip_fact is .true.')
       end if
    end if

#if (defined _MEMTRACE)
    if (iam  == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','DataInitialize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    !-----------------------------------------------------------------------
    ! Document orbital parameters
    !-----------------------------------------------------------------------

    if (registry_match('qsw_distrb_iopt_cosz')) then

       call pop_orbital_init(gcomp, stdout, my_task==master_task, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call pop_orbital_update(clock, stdout, my_task==master_task, orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       write(stdout,*) ' '
       call document ('DataInitialize', 'orb_eccen   ',  orb_eccen)
       call document ('DataInitialize', 'orb_mvelpp  ',  orb_mvelpp)
       call document ('DataInitialize', 'orb_lambm0  ',  orb_lambm0)
       call document ('DataInitialize', 'orb_obliqr  ',  orb_obliqr)
    endif

    !-----------------------------------------------------------------------
    ! check whether all Fields in the exportState are "Updated"
    !-----------------------------------------------------------------------

    if (NUOPC_IsUpdated(exportState)) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

      call ESMF_LogWrite("POP - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    !-----------------------------------------------------------------------
    ! Now document all time flags, last step of pop initialization
    !-----------------------------------------------------------------------

    call document_time_flags

    !-----------------------------------------------------------------------
    ! Output delimiter to log file
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,'(" End of initialization")')
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       call POP_IOUnitsFlush(POP_stdout)
       call POP_IOUnitsFlush(stdout)
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

  end subroutine DataInitialize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Run POP for a coupling interval
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    !  local variables
    type(ESMF_State)             :: importState
    type(ESMF_State)             :: exportState
    type(ESMF_Clock)             :: clock
    type(ESMF_Alarm)             :: alarm
    type(ESMF_Time)              :: currTime
    type(ESMF_Time)              :: nextTime
    character(ESMF_MAXSTR)       :: cvalue
    integer(int_kind)            :: errorCode  ! error flag
    integer(int_kind)            :: ymd        ! POP2 current date (YYYYMMDD)
    integer(int_kind)            :: tod        ! POP2 current time of day (sec)
    integer(int_kind)            :: ymd_sync   ! Sync clock current date (YYYYMMDD)
    integer(int_kind)            :: tod_sync   ! Sync clcok current time of day (sec)
    integer                      :: lbnum
    integer                      :: yr_sync
    integer                      :: mon_sync
    integer                      :: day_sync
    integer                      :: shrlogunit ! old values
    character(char_len)          :: message
    logical                      :: first_time = .true.
    character(len=*), parameter  :: subname = "ocn_comp_nuopc: (ModelAdvance)"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = POP_Success

    !-----------------------------------------------------------------------
    ! skip first coupling interval for an initial run
    !-----------------------------------------------------------------------

    ! NOTE: pop starts one coupling interval ahead of the rest of the system
    ! so to have it be in sync with the rest of the system, simply skip the first
    ! coupling interval for a initial run only

    if (first_time) then
       first_time = .false.
       if (runtype == 'initial') then
          if (my_task == master_task) then
             write(stdout,*)'Returning at first coupling interval'
          end if
          RETURN
       end if
    end if


#if (defined _MEMTRACE)
    if(my_task == 0 ) then
       lbnum=1
       call memmon_dump_fort('memmon.out',subname//':start::',lbnum)
    endif
#endif

    !--------------------------------------------------------------------
    ! check that pop internal clock is in sync with ESMF clock
    !--------------------------------------------------------------------

    ! NOTE: that in nuopc - the ESMF clock is updated at the end of the timestep
    ! whereas in cpl7 it was updated at the beginning - so need to have the check
    ! at the beginning of the time loop BEFORE pop updates its time step

    ! pop clock
    ymd = iyear*10000 + imonth*100 + iday
    tod = ihour*seconds_in_hour + iminute*seconds_in_minute + isecond

    ! model clock
    call NUOPC_ModelGet(gcomp, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)

    ! check
    if ( (ymd /= ymd_sync) .or. (tod /= tod_sync) ) then
       write(stdout,*)' pop2 ymd=',ymd     ,'  pop2 tod= ',tod
       write(stdout,*)' sync ymd=',ymd_sync,'  sync tod= ',tod_sync
       write(stdout,*)' Internal pop2 clock not in sync with Sync Clock'
       call ESMF_LogWrite(subname//" Internal POP clock not in sync with ESMF model clock", ESMF_LOGMSG_INFO)
       !call shr_sys_abort(subName// ":: Internal POP clock not in sync with ESMF model Clock")
    end if

    !-----------------------------------------------------------------------
    !  start up the main timer
    !-----------------------------------------------------------------------

    call timer_start(timer_total)

    !-----------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (stdout)

    if (ldiag_cpl) then
       call register_string ('info_debug_ge2')
    endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1 .and. master_task) then
       call log_clock_advance(clock, 'POP', stdout, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !----------------------------------------------------------------------------
    ! restart flag (rstwr) will assume only an eod restart for now
    !----------------------------------------------------------------------------

    ! Note this logic triggers off of the component clock rather than the internal pop time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call override_time_flag(cpl_write_restart, value=.true.)
       call ccsm_char_date_and_time ! set time_management module vars cyear, cmonth, ..
       write(message,'(6a)') 'driver requests restart file at eod  ', cyear,'/',cmonth,'/',cday
       call document ('ModelAdvance:', message)
    endif

    !-----------------------------------------------------------------------
    ! advance the model in time over coupling interval
    ! write restart dumps and archiving
    !-----------------------------------------------------------------------

    ! Note that all ocean time flags are evaluated each timestep in time_manager
    ! tlast_coupled is set to zero at the end of ocn_export

    advance: do

       ! -----
       ! obtain import state data
       ! -----
       if (check_time_flag(cpl_ts) .or. nsteps_run == 0) then

          call ocn_import(importState, flds_scalar_name, ldiag_cpl, errorCode, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (errorCode /= POP_Success) then
             call POP_ErrorPrint(errorCode)
             call exit_POP(sigAbort, 'ERROR in step')
          endif

          ! update orbital parameters

          call pop_orbital_update(clock, stdout, my_task==master_task, orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call pop_set_coupled_forcing
       end if

       ! -----
       ! advance pop
       ! -----

       call step(errorCode)

       if (errorCode /= POP_Success) then
          call POP_ErrorPrint(errorCode)
          call exit_POP(sigAbort, 'ERROR in step')
       endif

       if (check_KE(100.0_r8)) then
          !*** exit if energy is blowing
          call output_driver
          call exit_POP(sigAbort,'ERROR: k.e. > 100 ')
       endif
       call output_driver()

       ! -----
       ! create export state
       ! -----
       call pop_sum_buffer(exportState, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (check_time_flag(cpl_ts)) then

          call ocn_export(exportState, flds_scalar_name, ldiag_cpl, errorCode, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          if (errorCode /= POP_Success) then
             call POP_ErrorPrint(errorCode)
             call exit_POP(sigAbort, 'ERROR in ocn_export')
          endif

          if ( lsend_precip_fact ) then
             if ( flds_scalar_index_precip_factor > 0 ) then
                call State_SetScalar(precip_fact, flds_scalar_index_precip_factor, exportState, &
                     flds_scalar_name, flds_scalar_num, rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call shr_sys_abort('flds_scalar_index_precip_factor must be > 0 when lsend_precip_fact is .true.')
             end if
          end if

          exit advance
       end if

    enddo advance

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

    call timer_stop(timer_total)

#if (defined _MEMTRACE)
    if(my_task == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out',subname//':end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelSetRunClock)'
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelCheckImport(model, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)  :: clock
    type(ESMF_Time)   :: currTime
    integer(int_kind) :: yy  ! current date (YYYYMMDD)
    integer(int_kind) :: mon ! current month 
    integer(int_kind) :: day ! current day
    integer(int_kind) :: tod ! current time of day (sec)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
      
    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call ESMF_TimeGet(currTime, yy=yy, mm=mon, dd=day, s=tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) then
       write(stdout,*)' CheckImport pop2 year = ',yy
       write(stdout,*)' CheckImport pop2 mon  = ',mon
       write(stdout,*)' CheckImport pop2 day  = ',day
       write(stdout,*)' CheckImport pop2 tod  = ',tod
    end if

  end subroutine ModelCheckImport

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    !--------------------------------
    ! POP finalization that shuts down POP gracefully (we hope).
    ! Exits the message environment and checks for successful execution.
    ! --------------------------------

    use output          , only : final_output
    use passive_tracers , only : passive_tracers_timer_print_all

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer (POP_i4) :: errorCode         ! error code
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelFinalize)'
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call POP_ErrorPrint(errorCode, printTask=master_task)

    !  close any open files
    call final_output

    !  clear any open displays and print all timers with statistics
    call timer_print_all(stats=.true.)
    call passive_tracers_timer_print_all(stats=.true.)

    !  write final message to pop output log

    if (my_task == master_task) then
      write(stdout,*) '==================='
      write(stdout,*) 'completed POP_Final'
      write(stdout,*) '==================='
    endif

    !  exit the communication environment
    call exit_message_environment(ErrorCode)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine pop_orbital_init(gcomp, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask 
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(pop_orbital_init)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_mode

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_iyear

    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_iyear_align

    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_obliq

    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_eccen

    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) attribute_orb_mvelp

    ! Error checks
    if (trim(attribute_orb_mode) == trim(orb_fixed_year)) then
       attribute_orb_obliq = SHR_ORB_UNDEF_REAL
       attribute_orb_eccen = SHR_ORB_UNDEF_REAL
       attribute_orb_mvelp = SHR_ORB_UNDEF_REAL
       if (attribute_orb_iyear == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(attribute_orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',attribute_orb_iyear
             write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(attribute_orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(attribute_orb_mode) == trim(orb_variable_year)) then
       attribute_orb_obliq = SHR_ORB_UNDEF_REAL
       attribute_orb_eccen = SHR_ORB_UNDEF_REAL
       attribute_orb_mvelp = SHR_ORB_UNDEF_REAL
       if (attribute_orb_iyear == SHR_ORB_UNDEF_INT .or. attribute_orb_iyear_align == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(attribute_orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',attribute_orb_iyear, attribute_orb_iyear_align
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(attribute_orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    elseif (trim(attribute_orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       attribute_orb_iyear = SHR_ORB_UNDEF_INT
       attribute_orb_iyear_align = SHR_ORB_UNDEF_INT
       if ( attribute_orb_eccen == SHR_ORB_UNDEF_REAL .or. &
            attribute_orb_obliq == SHR_ORB_UNDEF_REAL .or. &
            attribute_orb_mvelp == SHR_ORB_UNDEF_REAL) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(attribute_orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',attribute_orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',attribute_orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',attribute_orb_mvelp
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(attribute_orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(attribute_orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif

  end subroutine pop_orbital_init

  !===============================================================================

  subroutine pop_orbital_update(clock, logunit,  mastertask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings 
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit 
    logical          , intent(in)    :: mastertask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time 
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    logical           :: first_time = .true.
    character(len=*) , parameter :: subname = "(pop_orbital_update)"
    !-------------------------------------------

    if (trim(attribute_orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = attribute_orb_iyear + (year - attribute_orb_iyear_align)
       lprint = mastertask
    else
       orb_year = attribute_orb_iyear 
       if (first_time) then
          lprint = mastertask
          first_time = .false.
       else
          lprint = .false.
       end if
    end if

    eccen = attribute_orb_eccen
    call shr_orb_params(orb_year, eccen, attribute_orb_obliq, attribute_orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( orb_eccen  == SHR_ORB_UNDEF_REAL .or. orb_obliqr == SHR_ORB_UNDEF_REAL .or. &
         orb_mvelpp == SHR_ORB_UNDEF_REAL .or. orb_lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine pop_orbital_update

end module ocn_comp_nuopc
