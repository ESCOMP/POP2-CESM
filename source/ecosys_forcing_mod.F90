module ecosys_forcing_mod

  ! !DESCRIPTION:
  !  This module sets up data types to keep track of the surface forcing data
  !  being sent to MARBL. Data may be passed through to MARBL from the flux
  !  coupler, read from a file, or set to a constant value based on namelist
  !  variables. Additionally, this module handles some interior forcing --
  !  namely tracer restoring.

  use constants, only : c0, c1

  use kinds_mod, only : r8, int_kind, log_kind, char_len, char_len_long

  use io_types, only : stdout

  use io_tools, only : document

  use exit_mod, only : sigAbort, exit_POP

  use passive_tracer_tools, only : tracer_read
  use passive_tracer_tools, only : forcing_monthly_every_ts

  use communicate, only : my_task, master_task

  use blocks, only : block, nx_block, ny_block

  use domain, only : nblocks_clinic

  use domain_size, only : km, max_blocks_clinic

  use timers, only : timer_start
  use timers, only : timer_stop
  use timers, only : get_timer

  use strdata_interface_mod     , only : strdata_input_type

  use ecosys_tracers_and_saved_state_mod, only : marbl_tracer_cnt

  implicit none
  private

  public :: ecosys_forcing_init
  public :: ecosys_forcing_set_surface_time_varying_forcing_data
  public :: ecosys_forcing_set_interior_time_varying_forcing_data
  public :: ecosys_forcing_tracer_ref_val

  !*****************************************************************************

  type, private :: forcing_constant_type
     real(kind=r8) :: field_constant           ! constant value for field_source
   contains
      procedure :: initialize  => forcing_constant_init
  end type forcing_constant_type


  !*****************************************************************************

  type, private :: forcing_driver_type
     character(char_len) :: driver_varname     ! Name of variable in coupler
   contains
     procedure :: initialize  => forcing_driver_init
  end type forcing_driver_type

  !*****************************************************************************

  type, private :: forcing_named_field_type
     character(char_len) :: field_name     ! Name of variable in named_field
     integer             :: field_ind
   contains
     procedure :: initialize  => forcing_named_field_init
  end type forcing_named_field_type

  !*****************************************************************************

  type, private :: forcing_file_type
     character (char_len)      :: filename
     character (char_len)      :: file_varname
     character (char_len)      :: temporal      ! temporarily to support current I/O routines
     integer   (kind=int_kind) :: year_first
     integer   (kind=int_kind) :: year_last
     integer   (kind=int_kind) :: year_align
     integer   (kind=int_kind) :: date
     integer   (kind=int_kind) :: time
   contains
     procedure :: initialize  => forcing_file_init
  end type forcing_file_type

  !*****************************************************************************

  type, private :: forcing_monthly_calendar_type
     type (forcing_monthly_every_ts), pointer :: forcing_calendar_name
   contains
     procedure :: initialize  => forcing_monthly_calendar_init
  end type forcing_monthly_calendar_type

  !*****************************************************************************

  type, private :: forcing_fields_metadata_type
     character(char_len)                  :: marbl_varname
     character(char_len)                  :: field_units      ! field data units, not the file (driver must do unit conversion)
     character(char_len)                  :: field_source     ! see valid_field_source in forcing_field_metadata_set
     character(char_len)                  :: temporal_interp  ! information on interpolation scheme used to populate field data
     real(kind=r8)                        :: unit_conv_factor ! unit conversion factor, incorporates scale_factor
     type (forcing_constant_type)         :: field_constant_info
     type (forcing_driver_type)           :: field_driver_info
     type (forcing_named_field_type)      :: field_named_info
     type (forcing_file_type)             :: field_file_info
     type (forcing_monthly_calendar_type) :: field_monthly_calendar_info
   contains
     procedure :: set  => forcing_field_metadata_set
  end type forcing_fields_metadata_type

  !*****************************************************************************

  type, private :: forcing_fields_type
    ! This is modelled after the marbl_forcing_fields_type, and is used to
    ! pass forcing data from POP to MARBL
    type(forcing_fields_metadata_type) :: metadata
    real(r8), allocatable :: field_0d(:,:,:)    ! (nx, ny, bid)
    real(r8), allocatable :: field_1d(:,:,:,:)  ! (nx, ny, [any dim], bid)
  contains
    procedure, public :: add_forcing_field => forcing_fields_add
  end type forcing_fields_type

  type(forcing_fields_type), dimension(:), allocatable, public :: surface_forcing_fields
  type(forcing_fields_type), dimension(:), allocatable, public :: interior_forcing_fields

  !---------------------------------------------------------------------
  !  Variables read in via &ecosys_forcing_data_nml
  !---------------------------------------------------------------------

  character(char_len) :: dust_flux_source             ! option for atmospheric dust deposition
  type(tracer_read)   :: dust_flux_input             ! namelist input for dust_flux
  character(char_len) :: iron_flux_source             ! option for atmospheric iron deposition
  type(tracer_read)   :: iron_flux_input             ! namelist input for iron_flux
  type(tracer_read)   :: fesedflux_input                    ! namelist input for iron_flux
  logical(log_kind)   :: lignore_driver_ndep
  character(char_len) :: ndep_data_type               ! type of ndep forcing
  type(tracer_read)   :: nox_flux_monthly_input      ! namelist input for nox_flux_monthly
  type(tracer_read)   :: nhy_flux_monthly_input      ! namelist input for nhy_flux_monthly
  integer(int_kind)   :: ndep_shr_stream_year_first   ! first year in stream to use
  integer(int_kind)   :: ndep_shr_stream_year_last    ! last year in stream to use
  integer(int_kind)   :: ndep_shr_stream_year_align   ! align ndep_shr_stream_year_first with this model year
  character(char_len) :: ndep_shr_stream_file         ! file containing domain and input data
  real(r8)            :: ndep_shr_stream_scale_factor ! unit conversion factor
  type(tracer_read)   :: din_riv_flux_input          ! namelist input for din_riv_flux
  type(tracer_read)   :: dip_riv_flux_input          ! namelist input for dip_riv_flux
  type(tracer_read)   :: don_riv_flux_input          ! namelist input for don_riv_flux
  type(tracer_read)   :: dop_riv_flux_input          ! namelist input for dop_riv_flux
  type(tracer_read)   :: dsi_riv_flux_input          ! namelist input for dsi_riv_flux
  type(tracer_read)   :: dfe_riv_flux_input          ! namelist input for dfe_riv_flux
  type(tracer_read)   :: dic_riv_flux_input          ! namelist input for dic_riv_flux
  type(tracer_read)   :: alk_riv_flux_input          ! namelist input for alk_riv_flux
  type(tracer_read)   :: doc_riv_flux_input          ! namelist input for doc_riv_flux
  character(char_len) :: gas_flux_forcing_opt        ! option for forcing gas fluxes
  character(char_len) :: gas_flux_forcing_file        ! file containing gas flux forcing fields
  type(tracer_read)   :: gas_flux_fice               ! ice fraction for gas fluxes
  type(tracer_read)   :: gas_flux_ws                 ! wind speed for gas fluxes
  type(tracer_read)   :: gas_flux_ap                 ! atmospheric pressure for gas fluxes
  character(char_len) :: atm_co2_opt                 ! option for atmospheric co2 concentration
  real(r8)            :: atm_co2_const                ! value of atmospheric co2 (ppm, dry-air, 1 atm)
  character(char_len) :: atm_alt_co2_opt             ! option for atmospheric alternative CO2
  real(r8)            :: atm_alt_co2_const            ! value of atmospheric alternative co2 (ppm, dry-air, 1 atm)
  logical(log_kind)   :: liron_patch                  ! flag for iron patch fertilization
  character(char_len) :: iron_patch_flux_filename     ! file containing name of iron patch file
  integer(int_kind)   :: iron_patch_month             ! integer month to add patch flux
  integer(int_kind)   :: ciso_atm_model_year            ! arbitrary model year
  integer(int_kind)   :: ciso_atm_data_year             ! year in atmospheric ciso data that corresponds to ciso_atm_model_year
  integer(int_kind)   :: ciso_atm_d13c_data_nbval       ! number of values in ciso_atm_d13c_filename
  integer(int_kind)   :: ciso_atm_d14c_data_nbval       ! number of values in ciso_atm_d14c_filename
  real(r8), allocatable :: ciso_atm_d13c_data(:)          ! atmospheric D13C values in datafile
  real(r8), allocatable :: ciso_atm_d13c_data_yr(:)       ! date of atmospheric D13C values in datafile
  real(r8), allocatable :: ciso_atm_d14c_data(:,:)        ! atmospheric D14C values in datafile (sh, eq, nh, in permil)
  real(r8), allocatable :: ciso_atm_d14c_data_yr(:,:)     ! date of atmospheric D14C values in datafile (sh, eq, nh)
  real(r8)            :: ciso_atm_d13c_const            ! atmospheric D13C constant [permil]
  real(r8)            :: ciso_atm_d14c_const            ! atmospheric D14C constant [permil]
  character(char_len) :: ciso_atm_d13c_opt              ! option for CO2 and D13C varying or constant forcing
  character(char_len) :: ciso_atm_d13c_filename         ! filenames for varying atm D13C
  character(char_len) :: ciso_atm_d14c_opt              ! option for CO2 and D13C varying or constant forcing
  character(char_len) :: ciso_atm_d14c_filename(3)      ! filenames for varying atm D14C (one each for NH, SH, EQ)

  ! For tracer restoring, we need three things:
  ! (1) list of tracers to apply restoring to
  character(char_len), dimension(marbl_tracer_cnt) :: restoreable_tracer_names
  ! (2) List of files containing the restoring fields
  character(char_len), dimension(marbl_tracer_cnt) :: restore_data_filenames
  ! (3) List of the name of the variable in (2) that we will restore (1) towards
  character(char_len), dimension(marbl_tracer_cnt) :: restore_data_file_varnames

  ! Also need to know what time scale to restore on
  character(char_len)                              :: restore_inv_tau_opt
  real(r8)                                         :: restore_inv_tau_const

  ! Is NDEP available from the driver?
  logical(log_kind), public :: ldriver_has_ndep

  !-----------------------------------------------------------------------
  !  input surface forcing
  !-----------------------------------------------------------------------

  ! These variables are necessary because of the way POP uses pointers store
  ! data that is read in via the forcing_monthly_every_ts data type. POP needs
  ! to allocate memory to store data read in from each file, and for each file
  ! *_file_loc%DATA points to that memory.
  type(forcing_monthly_every_ts), target :: dust_flux_file_loc
  type(forcing_monthly_every_ts), target :: iron_flux_file_loc
  type(forcing_monthly_every_ts), target :: fice_file_loc
  type(forcing_monthly_every_ts), target :: xkw_file_loc
  type(forcing_monthly_every_ts), target :: ap_file_loc
  type(forcing_monthly_every_ts), target :: nox_flux_monthly_file_loc
  type(forcing_monthly_every_ts), target :: nhy_flux_monthly_file_loc
  type(forcing_monthly_every_ts), target :: din_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: dip_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: don_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: dop_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: dsi_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: dfe_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: dic_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: alk_riv_flux_file_loc
  type(forcing_monthly_every_ts), target :: doc_riv_flux_file_loc

  integer (int_kind), parameter :: shr_stream_var_cnt    = 2 ! number of variables in ndep shr_stream
  integer (int_kind), parameter :: shr_stream_no_ind     = 1 ! index for NO forcing
  integer (int_kind), parameter :: shr_stream_nh_ind     = 2 ! index for NH forcing
  type    (strdata_input_type)  :: strdata_inputlist(shr_stream_var_cnt)  ! FIXME - need to make this more flexible

  ! Surface and interior forcing fields
  real(r8), target :: iron_patch_flux(nx_block, ny_block, max_blocks_clinic)
  real(r8)         :: dust_flux_in(nx_block, ny_block, max_blocks_clinic)


  !  ciso_data_ind_d13c is the index for the D13C data for the current timestep
  !  Note that ciso_data_ind_d13c is always less than ciso_atm_d13c_data_nbval.
  !  To enable OpenMP parallelism, duplicating data_ind for each block

  integer (int_kind), dimension(max_blocks_clinic) :: ciso_data_ind_d13c = -1 ! data index for D13C data
  integer (int_kind), dimension(max_blocks_clinic) :: ciso_data_ind_d14c = -1 ! data index for D14C data

  integer (int_kind) :: ecosys_pre_sflux_timer
  integer (int_kind) :: ecosys_shr_strdata_advance_timer

  ! Some surface forcing fields need special treatment, so we store indices
  integer(int_kind) :: dust_ind     = 0, &
                       Fe_ind       = 0, &
                       bc_ind       = 0, &
                       nox_ind      = 0, &
                       nhy_ind      = 0, &
                       xco2_ind     = 0, &
                       mask_ind     = 0, &
                       ifrac_ind    = 0, &
                       ap_ind       = 0, &
                       sst_ind      = 0, &
                       sss_ind      = 0, &
                       u10sqr_ind   = 0, &
                       d13c_ind     = 0, &
                       d14c_ind     = 0, &
                       d14c_glo_ind = 0

  ! We also need to track the indices of all the interior forcing fields
  integer(int_kind), public :: dustflux_ind       = 0, &
                       PAR_col_frac_ind   = 0, &
                       surf_shortwave_ind = 0, &
                       temperature_ind    = 0, &
                       salinity_ind       = 0, &
                       pressure_ind       = 0, &
                       fesedflux_ind      = 0

  !-----------------------------------------------------------------------
  ! Other private variables
  !-----------------------------------------------------------------------

  ! Virtual fluxes
  logical(log_kind), dimension(marbl_tracer_cnt) :: vflux_flag           ! which tracers get virtual fluxes applied
  real(r8), dimension(marbl_tracer_cnt) :: surf_avg                      ! average surface tracer values

  real(r8) :: iron_frac_in_dust
  real(r8) :: iron_frac_in_bc

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine ecosys_forcing_init(ciso_on, land_mask,                          &
                                 fe_frac_dust,fe_frac_bc, surface_forcings,   &
                                 interior_forcings, forcing_nml)

    use ecosys_tracers_and_saved_state_mod, only : set_defaults_tracer_read
    use ecosys_tracers_and_saved_state_mod, only : dic_ind
    use ecosys_tracers_and_saved_state_mod, only : dic_alt_co2_ind
    use ecosys_tracers_and_saved_state_mod, only : alk_ind
    use ecosys_tracers_and_saved_state_mod, only : di13c_ind
    use ecosys_tracers_and_saved_state_mod, only : di14c_ind

    use marbl_namelist_mod, only : marbl_nl_buffer_size
    use marbl_interface_types, only : marbl_forcing_fields_type

    use constants, only : delim_fmt, char_blank, ndelim_fmt

    use domain, only : distrb_clinic

    use mcog, only : mcog_nbins

    logical,                         intent(in)    :: ciso_on
    logical,                         intent(in)    :: land_mask(:,:,:)
    real(r8),                        intent(in)    :: fe_frac_dust
    real(r8),                        intent(in)    :: fe_frac_bc
    type(marbl_forcing_fields_type), intent(in)    :: surface_forcings(:)
    type(marbl_forcing_fields_type), intent(in)    :: interior_forcings(:)
    character(marbl_nl_buffer_size), intent(in)    :: forcing_nml

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter  :: subname = 'ecosys_forcing_mod:ecosys_forcing_init'
    character(len=char_len)  :: err_msg
    character(char_len)      :: marbl_varname, tracer_name, units
    integer (int_kind)       :: nml_error                  ! error flag for nml read
    character(char_len_long) :: ioerror_msg
    integer                  :: m, n
    type(forcing_monthly_every_ts), pointer :: file_details
    logical                  :: var_processed

    ! Virtual fluxes set in namelist
    real(r8) :: surf_avg_dic_const
    real(r8) :: surf_avg_alk_const
    real(r8) :: surf_avg_di13c_const
    real(r8) :: surf_avg_di14c_const

    namelist /ecosys_forcing_data_nml/                                        &
         dust_flux_source, dust_flux_input,                                   &
         iron_flux_source, iron_flux_input,                                   &
         fesedflux_input,                                                     &
         lignore_driver_ndep, ndep_data_type,                                 &
         nox_flux_monthly_input, nhy_flux_monthly_input,                      &
         ndep_shr_stream_year_first, ndep_shr_stream_year_last,               &
         ndep_shr_stream_year_align, ndep_shr_stream_file,                    &
         ndep_shr_stream_scale_factor, din_riv_flux_input,                    &
         dip_riv_flux_input, don_riv_flux_input, dop_riv_flux_input,          &
         dsi_riv_flux_input, dfe_riv_flux_input, dic_riv_flux_input,          &
         alk_riv_flux_input, doc_riv_flux_input, gas_flux_forcing_opt,        &
         gas_flux_forcing_file, gas_flux_fice, gas_flux_ws, gas_flux_ap,      &
         atm_co2_opt, atm_co2_const, atm_alt_co2_opt, atm_alt_co2_const,      &
         liron_patch, iron_patch_flux_filename, iron_patch_month,             &
         ciso_atm_d13c_opt, ciso_atm_d13c_const, ciso_atm_d13c_filename,      &
         ciso_atm_d14c_opt, ciso_atm_d14c_const, ciso_atm_d14c_filename,      &
         ciso_atm_model_year, ciso_atm_data_year, restore_data_filenames,     &
         restore_data_file_varnames, restoreable_tracer_names,                &
         restore_inv_tau_opt, restore_inv_tau_const, surf_avg_alk_const,      &
         surf_avg_dic_const, surf_avg_di13c_const, surf_avg_di14c_const

    !-----------------------------------------------------------------------
    !  Set module variables from intent(in)
    !-----------------------------------------------------------------------

    iron_frac_in_dust = fe_frac_dust
    iron_frac_in_bc   = fe_frac_bc

    !-----------------------------------------------------------------------
    !  &ecosys_forcing_data_nml
    !-----------------------------------------------------------------------

    gas_flux_forcing_opt  = 'drv'
    gas_flux_forcing_file = 'unknown'
    call set_defaults_tracer_read(gas_flux_fice, file_varname='FICE')
    call set_defaults_tracer_read(gas_flux_ws, file_varname='XKW')
    call set_defaults_tracer_read(gas_flux_ap, file_varname='P')
    dust_flux_source             = 'monthly-calendar'
    call set_defaults_tracer_read(dust_flux_input, file_varname='dust_flux')
    iron_flux_source             = 'monthly-calendar'
    call set_defaults_tracer_read(iron_flux_input, file_varname='iron_flux')
    call set_defaults_tracer_read(fesedflux_input, file_varname='FESEDFLUXIN')
    lignore_driver_ndep = .false.
    ndep_data_type = 'monthly-calendar'
    call set_defaults_tracer_read(nox_flux_monthly_input, file_varname='nox_flux')
    call set_defaults_tracer_read(nhy_flux_monthly_input, file_varname='nhy_flux')
    ndep_shr_stream_year_first = 1
    ndep_shr_stream_year_last  = 1
    ndep_shr_stream_year_align = 1
    ndep_shr_stream_file       = 'unknown'
    ndep_shr_stream_scale_factor = c1
    call set_defaults_tracer_read(din_riv_flux_input, file_varname='din_riv_flux')
    call set_defaults_tracer_read(dip_riv_flux_input, file_varname='dip_riv_flux')
    call set_defaults_tracer_read(don_riv_flux_input, file_varname='don_riv_flux')
    call set_defaults_tracer_read(dop_riv_flux_input, file_varname='dop_riv_flux')
    call set_defaults_tracer_read(dsi_riv_flux_input, file_varname='dsi_riv_flux')
    call set_defaults_tracer_read(dfe_riv_flux_input, file_varname='dfe_riv_flux')
    call set_defaults_tracer_read(dic_riv_flux_input, file_varname='dic_riv_flux')
    call set_defaults_tracer_read(alk_riv_flux_input, file_varname='alk_riv_flux')
    call set_defaults_tracer_read(doc_riv_flux_input, file_varname='doc_riv_flux')
    liron_patch              = .false.
    iron_patch_flux_filename = 'unknown_iron_patch_filename'
    iron_patch_month         = 1
    atm_co2_opt   = 'const'
    atm_co2_const = 280.0_r8
    atm_alt_co2_opt   = 'const'
    atm_alt_co2_const = 280.0_r8
    ciso_atm_d13c_opt                       = 'const'
    ciso_atm_d13c_const                     = -6.379_r8
    ciso_atm_d13c_filename                  = 'unknown'
    ciso_atm_d14c_opt                       = 'const'
    ciso_atm_d14c_const                     = 0.0_r8
    ciso_atm_d14c_filename(1)               = 'unknown'
    ciso_atm_d14c_filename(2)               = 'unknown'
    ciso_atm_d14c_filename(3)               = 'unknown'
    ciso_atm_model_year                     = 1
    ciso_atm_data_year                      = 1

    restore_data_filenames     = ''
    restore_data_file_varnames = ''
    restoreable_tracer_names   = ''
    restore_inv_tau_opt   = 'const'
    restore_inv_tau_const = 2.3e-6_r8

    surf_avg_alk_const   = 2225.0_r8
    surf_avg_dic_const   = 1944.0_r8
    surf_avg_di13c_const = 1944.0_r8
    surf_avg_di14c_const = 1944.0_r8

    read(forcing_nml, nml=ecosys_forcing_data_nml, iostat=nml_error, iomsg=ioerror_msg)
    if (nml_error /= 0) then
       write(stdout, *) subname, ": process ", my_task, ": namelist read error: ", nml_error, " : ", ioerror_msg
       call exit_POP(sigAbort, 'ERROR reading ecosys_forcing_data_nml from buffer.')
    end if

    if (my_task == master_task) then
       write(stdout,*)
       write(stdout,ndelim_fmt)
       write(stdout,*)
       write(stdout,*) ' ecosys_forcing_data_nml namelist settings:'
       write(stdout,*)
       write(stdout,ecosys_forcing_data_nml)
       write(stdout,*)
       write(stdout,delim_fmt)
    endif

    ! Consistency check
    if ((.not. lignore_driver_ndep) .and. (ldriver_has_ndep) .and. &
        (ndep_data_type.ne.'driver')) then
        ! If either named field index is greater than 0, then either
        ! ndep_data_type must be driver or lignore_driver_ndep must be .true.
        call document(subname, "ERROR: coupler is providing nitrogen deposition but ndep_data_type is not 'driver'")
        call document(subname, "To use the coupler values, set ndep_data_type = 'driver'")
        call document(subname, 'Otherwise set lignore_driver_ndep = .true.')
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    if ((ndep_data_type.eq.'driver') .and. (.not. ldriver_has_ndep)) then
        call document(subname, "ERROR: ndep_data_type is 'driver' but coupler is not providing nitrogen deposition!")
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    if (ciso_on) then
       call ciso_init_atm_D13_D14()
    end if

    ! Set surf_avg for all tracers
    surf_avg = c0
    vflux_flag(:) = .false.
    if (any((/dic_ind, dic_alt_co2_ind, alk_ind/).eq.0)) then
      call document(subname, 'dic_ind, alk_ind, and dic_alt_co2_ind must be non-zero')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    surf_avg(dic_ind) = surf_avg_dic_const
    vflux_flag(dic_ind) = .true.

    surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
    vflux_flag(dic_alt_co2_ind) = .true.

    surf_avg(alk_ind) = surf_avg_alk_const
    vflux_flag(alk_ind) = .true.
    
    if (ciso_on) then
       if (any((/di13c_ind, di14c_ind/).eq.0)) then
         call document(subname, 'di13c_ind and di14c_ind must be non-zero')
         call exit_POP(sigAbort, 'Stopping in ' // subname)
       end if
       surf_avg(di13c_ind) = surf_avg_di13c_const
       vflux_flag(di13c_ind) = .true.
       surf_avg(di14c_ind) = surf_avg_di14c_const
       vflux_flag(di14c_ind) = .true.
    end if

    call get_timer(ecosys_pre_sflux_timer    , 'ECOSYS_PRE_SFLUX'  , 1              , distrb_clinic%nprocs)
    if (ndep_data_type == 'shr_stream') then
       call get_timer(ecosys_shr_strdata_advance_timer, 'ecosys_shr_strdata_advance', 1, distrb_clinic%nprocs)
    endif

    !--------------------------------------------------------------------------
    !  Surface forcing
    !--------------------------------------------------------------------------

    ! Set up where surface forcing data will come from
    fice_file_loc%input             = gas_flux_fice
    xkw_file_loc%input              = gas_flux_ws
    ap_file_loc%input               = gas_flux_ap
    dust_flux_file_loc%input        = dust_flux_input
    iron_flux_file_loc%input        = iron_flux_input
    nox_flux_monthly_file_loc%input = nox_flux_monthly_input
    nhy_flux_monthly_file_loc%input = nhy_flux_monthly_input
    din_riv_flux_file_loc%input     = din_riv_flux_input
    dip_riv_flux_file_loc%input     = dip_riv_flux_input
    don_riv_flux_file_loc%input     = don_riv_flux_input
    dop_riv_flux_file_loc%input     = dop_riv_flux_input
    dsi_riv_flux_file_loc%input     = dsi_riv_flux_input
    dfe_riv_flux_file_loc%input     = dfe_riv_flux_input
    dic_riv_flux_file_loc%input     = dic_riv_flux_input
    alk_riv_flux_file_loc%input     = alk_riv_flux_input
    doc_riv_flux_file_loc%input     = doc_riv_flux_input

    allocate(surface_forcing_fields(size(surface_forcings)))
    do n=1,size(surface_forcing_fields)
      marbl_varname = surface_forcings(n)%metadata%varname
      units         = surface_forcings(n)%metadata%field_units
      select case (trim(surface_forcings(n)%metadata%varname))
        case ('surface_mask')
          mask_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='SURFACE_MASK', id=n)

        case ('d13c')
          d13c_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='D13C', id=n)

        case ('d14c')
          d14c_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='D14C', id=n)

        case ('d14c_gloavg')
          d14c_glo_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='D14C_GLOAVG', id=n)

        case ('u10_sqr')
          u10sqr_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='U10_SQR', id=n)

        case ('sst')
          sst_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='SST', id=n)

        case ('sss')
          sss_ind = n
          call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                               marbl_varname=marbl_varname, field_units=units,      &
                               driver_varname='SSS', id=n)

        case ('xco2')
          xco2_ind = n
          if (trim(atm_co2_opt).eq.'const') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='const', &
                                 marbl_varname=marbl_varname, field_units=units,   &
                                 field_constant=atm_co2_const, id=n)
          else if (trim(atm_co2_opt).eq.'drv_prog') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='named_field', &
                                 marbl_varname=marbl_varname, field_units=units,         &
                                 named_field='ATM_CO2_PROG', id=n)
          else if (trim(atm_co2_opt).eq.'drv_diag') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='named_field', &
                                 marbl_varname=marbl_varname, field_units=units,         &
                                 named_field='ATM_CO2_DIAG', id=n)
          else
            write(err_msg, "(A,1X,A)") trim(atm_co2_opt),                     &
                 'is not a valid option for atm_co2_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('xco2_alt_co2')
          if (trim(atm_alt_co2_opt).eq.'const') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='const', &
                                 marbl_varname=marbl_varname, field_units=units,   &
                                 field_constant=atm_alt_co2_const, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(atm_alt_co2_opt),                 &
                 'is not a valid option for atm_alt_co2_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Ice Fraction')
          ifrac_ind = n
          if (trim(gas_flux_forcing_opt).eq.'drv') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                                 marbl_varname=marbl_varname, field_units=units,      &
                                 driver_varname='ICE Fraction', id=n)
          else if (trim(gas_flux_forcing_opt).eq.'file') then
            file_details => fice_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(gas_flux_forcing_opt),            &
                 'is not a valid option for gas_flux_forcing_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Atmospheric Pressure')
          ap_ind = n
          if (trim(gas_flux_forcing_opt).eq.'drv') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                                 marbl_varname=marbl_varname, field_units=units,      &
                                 driver_varname='AP_FILE_INPUT', id=n)
          else if (trim(gas_flux_forcing_opt).eq.'file') then
            file_details => ap_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(gas_flux_forcing_opt),            &
                 'is not a valid option for gas_flux_forcing_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Dust Flux')
          dust_ind = n
          if (trim(dust_flux_source).eq.'driver') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                                 marbl_varname=marbl_varname, field_units=units,      &
                                 driver_varname='DUST_FLUX', id=n)
          else if (trim(dust_flux_source).eq.'monthly-calendar') then
            file_details => dust_flux_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(dust_flux_source),                &
                 'is not a valid option for dust_flux_source'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Iron Flux')
          if (trim(iron_flux_source).eq.'driver-derived') then
            bc_ind = n
            call surface_forcing_fields(n)%add_forcing_field(field_source='internal', &
                                 marbl_varname=marbl_varname, field_units=units,      &
                                 driver_varname='BLACK_CARBON_FLUX', id=n)
          else if (trim(iron_flux_source).eq.'monthly-calendar') then
            Fe_ind = n
            file_details => iron_flux_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(iron_flux_source),                &
                 'is not a valid option for iron_flux_source'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('NOx Flux')
          nox_ind = n
          if (trim(ndep_data_type).eq.'shr_stream') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='shr_stream', &
                                 marbl_varname=marbl_varname, field_units=units,  &
                                 unit_conv_factor=ndep_shr_stream_scale_factor,   &
                                 file_varname='NOy_deposition',                   &
                                 year_first = ndep_shr_stream_year_first,         &
                                 year_last = ndep_shr_stream_year_last,           &
                                 year_align = ndep_shr_stream_year_align,         &
                                 filename = ndep_shr_stream_file, id=n)
          else if (trim(ndep_data_type).eq.'monthly-calendar') then
            file_details => nox_flux_monthly_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else if (trim(ndep_data_type).eq.'driver') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='named_field', &
                                 marbl_varname=marbl_varname, field_units=units,         &
                                 named_field='ATM_NOy', id=n)
          else
            write(err_msg, "(A,1X,A)") trim(ndep_data_type),                  &
                 'is not a valid option for ndep_data_type'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('NHy Flux')
          nhy_ind = n
          if (trim(ndep_data_type).eq.'shr_stream') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='shr_stream', &
                                 marbl_varname=marbl_varname, field_units=units,  &
                                 unit_conv_factor=ndep_shr_stream_scale_factor,   &
                                 file_varname='NHx_deposition',                   &
                                 year_first = ndep_shr_stream_year_first,         &
                                 year_last = ndep_shr_stream_year_last,           &
                                 year_align = ndep_shr_stream_year_align,         &
                                 filename = ndep_shr_stream_file, id=n)
          else if (trim(ndep_data_type).eq.'monthly-calendar') then
            file_details => nhy_flux_monthly_file_loc
            call init_monthly_surface_forcing_metadata(file_details)
            call surface_forcing_fields(n)%add_forcing_field(                    &
                                 field_source='POP monthly calendar',            &
                                 marbl_varname=marbl_varname, field_units=units, &
                                 forcing_calendar_name=file_details, id=n)
          else if (trim(ndep_data_type).eq.'driver') then
            call surface_forcing_fields(n)%add_forcing_field(field_source='named_field', &
                                 marbl_varname=marbl_varname, field_units=units,         &
                                 named_field='ATM_NHx', id=n)
          else
            write(err_msg, "(A,1X,A)") trim(ndep_data_type),                  &
                 'is not a valid option for ndep_data_type'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('DIN River Flux')
          file_details => din_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DIP River Flux')
          file_details => dip_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DON River Flux')
          file_details => don_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DOP River Flux')
          file_details => dop_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DSi River Flux')
          file_details => dsi_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DFe River Flux')
          file_details => dfe_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DIC River Flux')
          file_details => dic_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('ALK River Flux')
          file_details => alk_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case ('DOC River Flux')
          file_details => doc_riv_flux_file_loc
          call init_monthly_surface_forcing_metadata(file_details)
          call surface_forcing_fields(n)%add_forcing_field(                    &
                               field_source='POP monthly calendar',            &
                               marbl_varname=marbl_varname, field_units=units, &
                               forcing_calendar_name=file_details, id=n)

        case DEFAULT
          write(err_msg, "(A,1X,A)") trim(surface_forcings(n)%metadata%varname), &
                         'is not a valid surface forcing field name.'
          call document(subname, err_msg)
          call exit_POP(sigAbort, 'Stopping in ' // subname)
      end select

      ! All surface forcing fields are 0d; if a 1d field is introduced later,
      ! move this allocate into the select case
      allocate(surface_forcing_fields(n)%field_0d(nx_block, ny_block, nblocks_clinic))

      ! Zero out forcing field. If a 1d field is introduced later, check to see
      ! which of field_0d and field_1d is allocated.
      surface_forcing_fields(n)%field_0d = c0
    end do

    !--------------------------------------------------------------------------
    !  Interior forcing
    !--------------------------------------------------------------------------

    allocate(interior_forcing_fields(size(interior_forcings)))

    do n=1,size(interior_forcing_fields)
      marbl_varname = interior_forcings(n)%metadata%varname
      units = interior_forcings(n)%metadata%field_units

      var_processed = .false.
      ! Check to see if this forcing field is tracer restoring
      if (index(marbl_varname,'Restoring Field').gt.0) then
        tracer_name = trim(marbl_varname(1:scan(marbl_varname,' ')))
        do m=1,marbl_tracer_cnt
          if (trim(tracer_name).eq.trim(restoreable_tracer_names(m))) then
            ! Check to make sure restore_data_filenames and
            ! restore_data_file_varnames have both been provided by namelist
            if (len_trim(restore_data_filenames(m)).eq.0) then
              write(err_msg, "(3A)") "No file provided to read restoring ",   &
                                     "field for ", trim(tracer_name)
              call document(subname, err_msg)
              call exit_POP(sigAbort, 'Stopping in ' // subname)
            end if
            if (len_trim(restore_data_file_varnames(m)).eq.0) then
              write(err_msg, "(3A)") "No variable name provided to read ",    &
                                     "restoring field for ", trim(tracer_name)
              call document(subname, err_msg)
              call exit_POP(sigAbort, 'Stopping in ' // subname)
            end if
            if (my_task.eq.master_task) then
              write(stdout, "(6A)") "Will restore ", trim(tracer_name),       &
                            " with ", trim(restore_data_file_varnames(m)),    &
                            " from ", trim(restore_data_filenames(m))
            end if
            call interior_forcing_fields(n)%add_forcing_field(                &
                       field_source='file_time_invariant',                    &
                       marbl_varname=marbl_varname, field_units=units,        &
                       filename=restore_data_filenames(m),                    &
                       file_varname=restore_data_file_varnames(m),            &
                       id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
            var_processed = .true.
            exit
          end if
        end do
      end if

      ! Check to see if this forcing field is a restoring time scale
      if (index(marbl_varname,'Restoring Inverse Timescale').gt.0) then
        tracer_name = trim(marbl_varname(1:scan(marbl_varname,' ')))
        select case (trim(restore_inv_tau_opt))
          case('const')
            call interior_forcing_fields(n)%add_forcing_field(                &
                       field_source='const',                                  &
                       marbl_varname=marbl_varname, field_units=units,        &
                       field_constant = restore_inv_tau_const,                &
                       id=n)
          ! case('shr_stream')
          ! NOT SUPPORTED YET
          ! will require additional namelist variables, and we can consider
          ! reading in one file per tracer instead of using the same mask
          ! for all restoring fields
          case DEFAULT
            write(err_msg, "(A,1X,A)") trim(restore_inv_tau_opt),             &
                 'is not a valid option for restore_inv_tau_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
        end select
        allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
        var_processed = .true.
      end if

      if (.not.var_processed) then
        select case (trim(interior_forcings(n)%metadata%varname))
          case ('Dust Flux')
            dustflux_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='dust_flux', id=n)
            allocate(interior_forcing_fields(n)%field_0d(nx_block, ny_block, nblocks_clinic))
          case ('PAR Column Fraction')
            PAR_col_frac_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='PAR_col_frac', id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, mcog_nbins, nblocks_clinic))
          case ('Surface Shortwave')
            surf_shortwave_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='surf_shortwave', id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, mcog_nbins, nblocks_clinic))
          case ('Temperature')
            temperature_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='temperature', id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
          case ('Salinity')
            salinity_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='salinity', id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
          case ('Pressure')
            pressure_ind = n
            call interior_forcing_fields(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='pressure', id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
          case ('Iron Sediment Flux')
            fesedflux_ind = n
            call interior_forcing_fields(n)%add_forcing_field(                &
                          field_source='file_time_invariant',                 &
                          marbl_varname=marbl_varname, field_units=units,     &
                          filename=fesedflux_input%filename,                  &
                          file_varname=fesedflux_input%file_varname,          &
                          id=n)
            allocate(interior_forcing_fields(n)%field_1d(nx_block, ny_block, km, nblocks_clinic))
          case DEFAULT
            write(err_msg, "(A,1X,A)") trim(interior_forcings(n)%metadata%varname), &
                           'is not a valid interior forcing field name.'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
        end select
      end if

      ! Zero out field
      if (allocated(interior_forcing_fields(n)%field_0d)) then
        interior_forcing_fields(n)%field_0d = c0
      else
        interior_forcing_fields(n)%field_1d = c0
      end if
    end do

    if ((bc_ind.ne.0).and.(dust_ind.eq.0)) then
      call document(subname, "If deriving iron flux, must provide dust flux!")
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    call set_time_invariant_forcing_data(surface_forcing_fields)
    call set_time_invariant_forcing_data(interior_forcing_fields)
    call read_monthly_calendar_forcing_data(surface_forcing_fields, land_mask)
    call forcing_init_post_processing(land_mask)
    call mask_time_invariant_forcing_data(surface_forcing_fields, land_mask)
    call mask_time_invariant_forcing_data(interior_forcing_fields, land_mask)

  end subroutine ecosys_forcing_init

  !*****************************************************************************

  subroutine set_time_invariant_forcing_data(forcing_fields_in)

    use passive_tracer_tools, only : read_field

    type(forcing_fields_type), intent(inout) :: forcing_fields_in(:)

    integer :: n

    !--------------------------------------------------------------------------
    !  Set any constant forcing fields and
    !  read any forcing fields that are time invariant
    !
    !  For POP monthly calendar surface forcing fields, read all
    !  12 months of data
    !--------------------------------------------------------------------------

    do n=1, size(forcing_fields_in)
      associate (forcing_field =>forcing_fields_in(n)%metadata)
        select case (trim(forcing_field%field_source))
          case ('const','zero')
            if (allocated(forcing_fields_in(n)%field_0d)) then
              forcing_fields_in(n)%field_0d =                                 &
                             forcing_field%field_constant_info%field_constant
            else
              forcing_fields_in(n)%field_1d =                                 &
                             forcing_field%field_constant_info%field_constant
            end if
          case ('file_time_invariant')
            if (allocated(forcing_fields_in(n)%field_0d)) then
              call read_field('nc', forcing_field%field_file_info%filename,   &
                              forcing_field%field_file_info%file_varname,     &
                              forcing_fields_in(n)%field_0d)
            else
              call read_field('nc', forcing_field%field_file_info%filename,   &
                              forcing_field%field_file_info%file_varname,     &
                              forcing_fields_in(n)%field_1d)
            end if
        end select

      end associate
    end do

  end subroutine set_time_invariant_forcing_data

  !*****************************************************************************

  subroutine read_monthly_calendar_forcing_data(forcing_fields_in, land_mask)

    use passive_tracer_tools, only : read_field
    use forcing_tools       , only : find_forcing_times

    type(forcing_fields_type), intent(inout) :: forcing_fields_in(:)
    logical,                   intent(in)    :: land_mask(:,:,:)

    type(forcing_monthly_every_ts), pointer :: file
    real(r8), allocatable, target :: work_read(:,:,:,:)
    integer :: m, n, iblock

    !--------------------------------------------------------------------------
    !  For POP monthly calendar surface forcing fields, read all
    !  12 months of data
    !--------------------------------------------------------------------------

    do n=1, size(forcing_fields_in)
      associate (forcing_field =>forcing_fields_in(n)%metadata)
        if (trim(forcing_field%field_source).eq.'POP monthly calendar') then
           file => forcing_field%field_monthly_calendar_info%forcing_calendar_name

           allocate(work_read(nx_block, ny_block, 12, max_blocks_clinic))  
           if (trim(file%input%filename) == 'unknown') then
              file%input%filename = gas_flux_forcing_file
           end if
           if (trim(file%input%filename) /= 'none') then
              allocate(file%data(nx_block, ny_block, max_blocks_clinic, 1, 12))
              call read_field(file%input%file_fmt, file%input%filename, file%input%file_varname, work_read)
              do iblock=1, nblocks_clinic
                 do m=1, 12
                    file%data(:, :, iblock, 1, m) = work_read(:, :, m, iblock)
                    where (.not. land_mask(:, :, iblock)) file%data(:, :, iblock, 1, m) = c0
                    file%data(:, :, iblock, 1, m) = file%data(:, :, iblock, 1, m) * file%input%scale_factor
                 end do
              end do
              call find_forcing_times(     &
                   file%data_time, file%data_inc, file%interp_type, file%data_next, &
                   file%data_time_min_loc, file%data_update, file%data_type)
              file%has_data = .true.
           else
              file%has_data = .false.
           endif

           !  load iron PATCH flux fields (if required)
           !  assume patch file has same normalization and format as deposition file
           if (n == Fe_ind .and. liron_patch) then
              call read_field(file%input%file_fmt, file%input%filename, iron_patch_flux_filename, iron_patch_flux)
              do iblock=1, nblocks_clinic
                 do m=1, 12
                    where (.not. land_mask(:, :, iblock)) iron_patch_flux(:, :, iblock) = c0
                    file%data(:, :, iblock, 1, m) = iron_patch_flux(:, :, iblock) * file%input%scale_factor
                 end do
              end do
           end if

           deallocate(work_read)

        end if
      end associate
    end do

  end subroutine read_monthly_calendar_forcing_data

  !*****************************************************************************

  subroutine forcing_init_post_processing(land_mask)

    use passive_tracer_tools, only : read_field
    use grid                , only : KMT

    logical, intent(in)  :: land_mask(:,:,:)

    integer :: i, j, iblock, k, n
    real (r8) :: subsurf_fesed      ! sum of subsurface fesed values

    !-----------------------------------------------------------------------
    !  subsurf_fesed adjustment
    !-----------------------------------------------------------------------

    associate (fesedflux => interior_forcing_fields(fesedflux_ind)%field_1d)
    do iblock=1,nblocks_clinic
      do j=1, ny_block
        do i=1, nx_block
          if (KMT(i, j, iblock) > 0 .and. KMT(i, j, iblock) < km) then
            subsurf_fesed = c0
            do k=KMT(i, j, iblock)+1, km
              subsurf_fesed = subsurf_fesed + fesedflux(i, j, k, iblock)
            enddo
            fesedflux(i, j, KMT(i, j, iblock), iblock) = fesedflux(i, j, KMT(i, j, iblock), iblock) + subsurf_fesed
          endif
        enddo
      enddo

      do k = 1, km
        where (land_mask(:, :, iblock) .and. (k.le.KMT(:, :, iblock)))
          fesedflux(:, :, k, iblock) = fesedflux(:, :, k, iblock) * fesedflux_input%scale_factor
        end where
      enddo
    enddo
    end associate

  end subroutine forcing_init_post_processing

  !*****************************************************************************

  subroutine mask_time_invariant_forcing_data(forcing_fields_in, land_mask)

    use grid, only : KMT

    type(forcing_fields_type), intent(inout) :: forcing_fields_in(:)
    logical,                   intent(in)    :: land_mask(:,:,:)

    integer :: k, iblock, n

    do n=1,size(forcing_fields_in)
      associate (forcing_field =>forcing_fields_in(n)%metadata)
        if ((trim(forcing_field%field_source).eq.'const') .or.                &
            (trim(forcing_field%field_source).eq.'zero')  .or.                &
            (trim(forcing_field%field_source).eq.'file_time_invariant')) then
          do iblock=1,nblocks_clinic
            if (allocated(forcing_fields_in(n)%field_0d)) then
              where (.not.land_mask(:, :, iblock))
                forcing_fields_in(n)%field_0d(:,:,iblock) = c0
              endwhere
            else
              do k=1,km
                where (.not.land_mask(:, :, iblock) .or. (k.gt.KMT(:, :, iblock)))
                  forcing_fields_in(n)%field_1d(:,:,k,iblock) = c0
                end where
              end do
            end if
          end do
        end if
      end associate
    end do
    
  end subroutine mask_time_invariant_forcing_data

  !*****************************************************************************

  subroutine ecosys_forcing_set_interior_time_varying_forcing_data(           &
                                         FRACR_BIN, QSW_RAW_BIN, QSW_BIN,     &
                                         temperature, salinity, pressure,     &
                                         ecosys_qsw_distrb_const, bid)
  
    use mcog, only : mcog_nbins

    real(r8), dimension(nx_block, ny_block, mcog_nbins), intent(in) :: FRACR_BIN
    real(r8), dimension(nx_block, ny_block, mcog_nbins), intent(in) :: QSW_RAW_BIN
    real(r8), dimension(nx_block, ny_block, mcog_nbins), intent(in) :: QSW_BIN
    real(r8), dimension(nx_block, ny_block, km),         intent(in) :: temperature
    real(r8), dimension(nx_block, ny_block, km),         intent(in) :: salinity
    real(r8), dimension(nx_block, ny_block, km),         intent(in) :: pressure
    logical, intent(in) :: ecosys_qsw_distrb_const
    integer, intent(in) :: bid

    logical(log_kind) :: first_call = .true.
    integer :: index

    if (first_call) then
      ! If in the future interior forcings come from shr_stream we will need
      ! the first_call block
      first_call = .false.
    end if

    do index= 1,size(interior_forcing_fields)
      select case (interior_forcing_fields(index)%metadata%field_source)
        case('internal')
          if (index.eq.dustflux_ind) then
            interior_forcing_fields(index)%field_0d(:,:,bid) = dust_flux_in(:, :, bid)
          else if (index.eq.PAR_col_frac_ind) then
            interior_forcing_fields(index)%field_1d(:,:,:,bid) = FRACR_BIN
          else if (index.eq.surf_shortwave_ind) then
            if (ecosys_qsw_distrb_const) then
              interior_forcing_fields(index)%field_1d(:,:,:,bid) = QSW_RAW_BIN
            else
              interior_forcing_fields(index)%field_1d(:,:,:,bid) = QSW_BIN
            end if
          else if (index.eq.temperature_ind) then
            interior_forcing_fields(index)%field_1d(:,:,:,bid) = temperature
          else if (index.eq.salinity_ind) then
            interior_forcing_fields(index)%field_1d(:,:,:,bid) = salinity
          else if (index.eq.pressure_ind) then
            interior_forcing_fields(index)%field_1d(:,:,:,bid) = pressure
          end if
      end select
    end do

    call adjust_interior_time_varying_data()

  end subroutine ecosys_forcing_set_interior_time_varying_forcing_data

  !***********************************************************************

  subroutine adjust_interior_time_varying_data()
    ! This subroutine is empty because there are no interior forcing fields
    ! that need to be modified before being passed to MARBL
  end subroutine adjust_interior_time_varying_data

  !*****************************************************************************

  subroutine ecosys_forcing_set_surface_time_varying_forcing_data( &
       ciso_on,                               &
       land_mask,                             &
       u10_sqr,                               &
       ifrac,                                 &
       press,                                 &
       dust_flux,                             &
       black_carbon_flux,                     &
       sst,                                   &
       sss)

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.

    use POP_HaloMod           , only : POP_HaloUpdate 
    use POP_GridHorzMod       , only : POP_gridHorzLocCenter 
    use POP_CommMod           , only : POP_communicator 
    use POP_FieldMod          , only : POP_fieldKindScalar
    use POP_ErrorMod          , only : POP_Success
    use domain                , only : POP_haloClinic
    use domain                , only : blocks_clinic
    use blocks                , only : get_block
    use constants             , only : field_loc_center
    use constants             , only : field_type_scalar
    use constants             , only : xkw_coeff
    use forcing_tools         , only : interpolate_forcing
    use forcing_tools         , only : update_forcing_data
    use named_field_mod       , only : named_field_get
    use time_management       , only : isecond
    use time_management       , only : iminute
    use time_management       , only : ihour
    use time_management       , only : iday
    use time_management       , only : imonth
    use time_management       , only : iyear
    use time_management       , only : thour00
    use strdata_interface_mod , only : POP_strdata_advance 
    use strdata_interface_mod , only : POP_strdata_create
    use passive_tracer_tools  , only : read_field
    use marbl_constants_mod   , only : molw_Fe


    implicit none

    logical,   intent(in)  :: ciso_on
    logical,   intent(in)  :: land_mask            (nx_block,ny_block,max_blocks_clinic)
    real (r8), intent(in)  :: u10_sqr              (nx_block,ny_block,max_blocks_clinic) ! 10m wind speed squared (cm/s)**2
    real (r8), intent(in)  :: ifrac                (nx_block,ny_block,max_blocks_clinic) ! sea ice fraction (non-dimensional)
    real (r8), intent(in)  :: press                (nx_block,ny_block,max_blocks_clinic) ! sea level atmospheric pressure (dyne/cm**2)
    real (r8), intent(in)  :: dust_flux            (nx_block,ny_block,max_blocks_clinic) ! dust flux (g/cm**2/s)
    real (r8), intent(in)  :: black_carbon_flux    (nx_block,ny_block,max_blocks_clinic) ! black carbon flux (g/cm**2/s)
    real (r8), intent(in)  :: sst                  (nx_block,ny_block,max_blocks_clinic) ! sea surface temperature (c)
    real (r8), intent(in)  :: sss                  (nx_block,ny_block,max_blocks_clinic) ! sea surface salinity (psu)    

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (*), parameter       :: subname = 'ecosys_forcing_mod:ecosys_forcing_set_surface_time_varying_forcing_data'
    logical   (log_kind)           :: first_call = .true.
    type      (block)              :: this_block                                            ! block info for the current block
    integer   (int_kind)           :: index                                                 ! field index
    integer   (int_kind)           :: i, j, iblock, n                                       ! loop indices
    integer   (int_kind)           :: errorCode                                             ! errorCode from HaloUpdate call
    integer   (int_kind)           :: tracer_bndy_loc(1)                                    ! location   for ghost tracer_bndy_type cell updates
    integer   (int_kind)           :: tracer_bndy_type(1)                                   ! field type for ghost tracer_bndy_type cell updates
    character (char_len)           :: tracer_data_names(1)                                  ! short names for input data fields
    real      (r8)                 :: interp_work(nx_block, ny_block, max_blocks_clinic, 1) ! temp array for interpolate_forcing output
    real      (r8)                 :: shr_stream(nx_block, ny_block, max_blocks_clinic)
    real      (r8)                 :: d13c(nx_block, ny_block, max_blocks_clinic)           ! atm 13co2 value
    real      (r8)                 :: d14c(nx_block, ny_block, max_blocks_clinic)           ! atm 14co2 value
    real      (r8)                 :: d14c_glo_avg                                          ! global average D14C over the ocean, computed from current D14C field
    type(forcing_monthly_every_ts), pointer :: file
    integer   (int_kind)          :: stream_index
    integer   (int_kind)          :: nf_ind
    !-----------------------------------------------------------------------

    call timer_start(ecosys_pre_sflux_timer)

    !-----------------------------------------------------------------------
    ! Update carbon isotope atmosphere deltas if appropriate
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ciso_update_atm_D13C_D14C(land_mask, d13c, d14c, d14c_glo_avg)
    end if

    !-----------------------------------------------------------------------
    ! Initialize strdata_inputlist data type (only once)
    !-----------------------------------------------------------------------

    if (first_call) then
       do index = 1, size(surface_forcing_fields)
       associate(fields => surface_forcing_fields(index)%metadata)
          if (fields%field_source.eq.'shr_stream') then
             if (trim(ndep_data_type) == 'shr_stream') then
                stream_index = 0
                if (index == nox_ind) stream_index = shr_stream_no_ind
                if (index == nhy_ind) stream_index = shr_stream_nh_ind

                if (stream_index /= 0) then
                   strdata_inputlist(stream_index)%timer_label = 'marbl_file'
                   strdata_inputlist(stream_index)%year_first  = fields%field_file_info%year_first
                   strdata_inputlist(stream_index)%year_last   = fields%field_file_info%year_last
                   strdata_inputlist(stream_index)%year_align  = fields%field_file_info%year_align
                   strdata_inputlist(stream_index)%file_name   = fields%field_file_info%filename
                   strdata_inputlist(stream_index)%field_list  = fields%field_file_info%file_varname

                   call POP_strdata_create(strdata_inputlist(stream_index)) !FIXME - need more general scheme
                end if
             end if
          end if
       end associate
       end do
       first_call = .false.
    end if ! first_call

    if (trim(ndep_data_type) == 'shr_stream') then
       !FIXME - this should be moved into select for the case file
       ! call timer_start(ecosys_shr_strdata_advance_timer) FIXME - does not work
       do n = 1, shr_stream_var_cnt
          strdata_inputlist(n)%date = iyear*10000 + imonth*100 + iday
          strdata_inputlist(n)%time = isecond + 60 * (iminute + 60 * ihour)
          call POP_strdata_advance(strdata_inputlist(n))
       end do
       !call timer_stop(ecosys_shr_strdata_advance_timer) FIXME - does not work
    end if

    !-----------------------------------------------------------------------
    !  loop through forcing fields
    !-----------------------------------------------------------------------

    do index = 1, size(surface_forcing_fields)
    associate(fields => surface_forcing_fields(index)%metadata)

       select case (fields%field_source)

       !------------------------------------
       case ('POP monthly calendar')
       !------------------------------------

          file => fields%field_monthly_calendar_info%forcing_calendar_name

          if (file%has_data) then
             if (thour00 >= file%data_update) then
                tracer_data_names(1) = file%input%file_varname
                tracer_bndy_loc(1)   = field_loc_center
                tracer_bndy_type(1)  = field_type_scalar
                call update_forcing_data(                                &
                     forcing_time         = file%data_time,              &
                     forcing_time_min_loc = file%data_time_min_loc,      &
                     forcing_interp_type  = file%interp_type,            &
                     forcing_data_next    = file%data_next,              &
                     forcing_data_update  = file%data_update,            &
                     forcing_data_type    = file%data_type,              &
                     forcing_data_inc     = file%data_inc,               &
                     field                = file%data(:, :, :, :, 1:12), &
                     forcing_data_rescale = file%data_renorm,            &
                     forcing_data_label   = fields%marbl_varname, &
                     forcing_data_names   = tracer_data_names,           &
                     forcing_bndy_loc     = tracer_bndy_loc,             &
                     forcing_bndy_type    = tracer_bndy_type,            &
                     forcing_infile       = file%filename,               &
                     forcing_infile_fmt   = file%input%file_fmt)
             endif

             call interpolate_forcing(                                &
                  interp               = interp_work,                 &
                  field                = file%data(:, :, :, :, 1:12), &
                  forcing_time         = file%data_time,              &
                  forcing_interp_type  = file%interp_type,            &
                  forcing_time_min_loc = file%data_time_min_loc,      &
                  forcing_interp_freq  = file%interp_freq,            &
                  forcing_interp_inc   = file%interp_inc,             &
                  forcing_interp_next  = file%interp_next,            &
                  forcing_interp_last  = file%interp_last,            &
                  nsteps_run_check     = 0)

             surface_forcing_fields(index)%field_0d = interp_work(:, :, :, 1)
          endif

       !------------------------------------
       case ('named_field')
       !------------------------------------

       do iblock = 1,nblocks_clinic
          call named_field_get(fields%field_named_info%field_ind, iblock, &
                               surface_forcing_fields(index)%field_0d(:,:,iblock))
       end do
       !------------------------------------
       case ('internal')
       !------------------------------------

          do iblock = 1,nblocks_clinic

             if (index == mask_ind) then
                where(land_mask(:,:,iblock))
                  surface_forcing_fields(index)%field_0d(:,:,iblock) = c1
                elsewhere
                  surface_forcing_fields(index)%field_0d(:,:,iblock) = c0
                end where

             else if (index == ifrac_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = ifrac(:,:,iblock)

             else if (index == ap_ind) then
                !  assume PRESS is in cgs units (dyne/cm**2) since that is what is
                !    required for pressure forcing in barotropic
                !  want units to be atmospheres
                !  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
                !  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
                surface_forcing_fields(index)%field_0d(:,:,iblock) = press(:,:,iblock) / 101.325e+4_r8

             else if (index == sst_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = sst(:,:,iblock)

             else if (index == sss_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = sss(:,:,iblock)

             else if (index == bc_ind) then
                ! compute iron_flux in gFe/cm^2/s, then convert to nmolFe/cm^2/s
                surface_forcing_fields(index)%field_0d(:,:,iblock) =          &
                     (1.0e9_r8 / molw_Fe) *                                   &
                     ((dust_flux(:,:,iblock) * 0.98_r8) * iron_frac_in_dust + &
                      black_carbon_flux(:,:,iblock) * iron_frac_in_bc)

             else if (index == u10sqr_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = u10_sqr(:,:,iblock)

             else if (index == dust_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = dust_flux(:,:,iblock)

             else if (index == d13c_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = d13c(:,:,iblock)

             else if (index == d14c_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = d14c(:,:,iblock)
                
             else if (index == d14c_glo_ind) then
                surface_forcing_fields(index)%field_0d(:,:,iblock) = d14c_glo_avg

             end if  ! index

          end do

       !------------------------------------
       case ('shr_stream')
       !------------------------------------

          ! FIXME - move stream_index in forcing_field_file_type 
          stream_index = 0
          if (index == nox_ind) stream_index = shr_stream_no_ind
          if (index == nhy_ind) stream_index = shr_stream_nh_ind

          if (stream_index /= 0) then
             n = 0
             do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock), iblock)
                do j = this_block%jb, this_block%je
                   do i = this_block%ib, this_block%ie
                      n = n + 1
                      ! Note that each stream currently is assumed to have only 1 field in its
                      ! attribute vector
                      shr_stream(i, j, iblock) = strdata_inputlist(stream_index)%sdat%avs(1)%rAttr(1, n)
                   enddo
                enddo
             enddo

             call POP_HaloUpdate(shr_stream, POP_haloClinic, &
                  POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
             if (errorCode /= POP_Success) then
                call document(subname, 'error updating halo for Ndep fields')
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             endif
             
             do iblock = 1, nblocks_clinic
                where (land_mask(:, :, iblock))
                   surface_forcing_fields(index)%field_0d(:,:,iblock) =       &
                                                     shr_stream(:, :, iblock)
                endwhere
             enddo
          end if
             
       end select ! file, constant, driver, shr_stream

       if (fields%unit_conv_factor /= c1) then
          do iblock = 1, nblocks_clinic
             where (land_mask(:, :, iblock))
                surface_forcing_fields(index)%field_0d(:,:,iblock) =          &
                       fields%unit_conv_factor*                        &
                       surface_forcing_fields(index)%field_0d(:,:,iblock)
             endwhere
          enddo
       end if
          
    end associate
    end do  ! index

    call adjust_surface_time_varying_data()

    call timer_stop(ecosys_pre_sflux_timer)

  end subroutine ecosys_forcing_set_surface_time_varying_forcing_data

  !***********************************************************************

  subroutine adjust_surface_time_varying_data()

    use time_management, only : imonth
    use marbl_parms    , only : parm_Fe_bioavail

    integer :: index, iblock

    !-----------------------------------------------------------------------
    ! Some surface forcing fields need to be modified before being
    ! sent to MARBL
    !-----------------------------------------------------------------------

    do iblock = 1,nblocks_clinic

       ! Reduce surface dust flux due to assumed instant surface dissolution
       index = dust_ind
       if (index.gt.0) then
         surface_forcing_fields(index)%field_0d(:,:,iblock) =                 &
                 surface_forcing_fields(index)%field_0d(:,:,iblock) * 0.98_r8
       end if

       ! Add iron patch (if available)
       ! Apply bioavail scaling
       index = Fe_ind
       if (index.gt.0) then
         if (liron_patch .and. imonth == iron_patch_month) then
           surface_forcing_fields(index)%field_0d(:,:,iblock) =               &
                         surface_forcing_fields(index)%field_0d(:,:,iblock) + &
                         iron_patch_flux(:,:,iblock)
         endif
         surface_forcing_fields(index)%field_0d(:,:,iblock) = surface_forcing_fields(index)%field_0d(:,:,iblock) * parm_Fe_bioavail
       endif

       ! Add iron patch (if available)
       index = bc_ind
       if (index.gt.0) then
         if (liron_patch .and. imonth == iron_patch_month) then
           surface_forcing_fields(index)%field_0d(:,:,iblock) =               &
                         surface_forcing_fields(index)%field_0d(:,:,iblock) + &
                         iron_patch_flux(:,:,iblock) * parm_Fe_bioavail
         endif
       endif

       index = ifrac_ind
       if (index.gt.0) then
         if (surface_forcing_fields(index)%metadata%field_source == 'internal') then
           where (surface_forcing_fields(index)%field_0d(:,:,iblock) < c0)    &
                surface_forcing_fields(index)%field_0d(:,:,iblock) = c0
           where (surface_forcing_fields(index)%field_0d(:,:,iblock) > c1)    &
                surface_forcing_fields(index)%field_0d(:,:,iblock) = c1
         else
           ! Apply OCMIP ice fraction mask when input is from a file.
           where (surface_forcing_fields(index)%field_0d(:,:,iblock) < 0.2000_r8) &
                surface_forcing_fields(index)%field_0d(:,:,iblock) = 0.2000_r8
           where (surface_forcing_fields(index)%field_0d(:,:,iblock) > 0.9999_r8) &
                surface_forcing_fields(index)%field_0d(:,:,iblock) = 0.9999_r8
         end if
       end if

    end do

    index = dust_ind
    if (index.gt.0) then
      dust_flux_in = surface_forcing_fields(index)%field_0d
    end if

  end subroutine adjust_surface_time_varying_data

  !***********************************************************************

  subroutine ciso_update_atm_D13C_D14C (land_mask, D13C, D14C, D14C_glo_avg)

    ! Updates module variables D13C and D14C (for atmospheric ratios)

    use grid              , only : TAREA
    use domain            , only : blocks_clinic
    use domain            , only : distrb_clinic
    use blocks            , only : get_block
    use global_reductions , only : global_sum

    implicit none

    logical,   intent(in)  :: land_mask(nx_block, ny_block, max_blocks_clinic)
    real (r8), intent(out) :: D13C(nx_block, ny_block, max_blocks_clinic)  ! atm 13co2 value
    real (r8), intent(out) :: D14C(nx_block, ny_block, max_blocks_clinic)  ! atm 14co2 value
    real (r8), intent(out) :: D14C_glo_avg  ! global average D14C over the ocean, computed from current D14C field

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_update_atm_D13C_D14C'
    type (block) :: &
         this_block      ! block info for the current block

    real (r8), dimension(nx_block,ny_block) :: &
         work1, &! local work space
         tfact   ! factor for normalizing sums

    integer (int_kind) :: &
         ib,ie,jb,je, &
         iblock  ! index for looping over blocks

    real (r8), dimension(max_blocks_clinic) :: &
         d14c_local_sums, & ! array for holding block sums when calculating global D14C
         tarea_local_sums   ! array for holding block sums of TAREA when calculating global D14C

    real (r8) :: &
         d14c_sum_tmp,  & ! temp for local sum of D14C
         tarea_sum_tmp    ! temp for local sum of TAREA
    !-----------------------------------------------------------------------

    work1(:,:) = c0
    d14c_local_sums(:)  = c0
    tarea_local_sums(:) = c0
    
    !-----------------------------------------------------------------------
    ! Loop over blocks
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic

       !-----------------------------------------------------------------------
       !  Set D13C (constant or from files read in _init) and put on global grid
       !-----------------------------------------------------------------------

       select case (ciso_atm_d13c_opt)
       case ('const')
          D13C(:,:,:) = ciso_atm_d13c_const
       case ('file')
          call ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c(iblock), D13C(:,:,iblock))
       case default
          call document(subname, 'unknown ciso_atm_d13c_opt in ecosys_ciso_set_sflux')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
       end select

       !-----------------------------------------------------------------------
       !  Set D14C (constant or from files read in _init) and put on global grid
       !-----------------------------------------------------------------------

       select case (ciso_atm_d14c_opt)
       case ('const')
          D14C(:,:,:) = ciso_atm_d14c_const
       case ('file')
          call ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c(iblock), D14C(:,:,iblock))
       case default
          call document(subname, 'unknown ciso_atm_d14c_opt in ecosys_ciso_set_sflux')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
       end select

       !-----------------------------------------------------------------------
       ! Save local D14C field for making global mean after end of iblock loop
       !-----------------------------------------------------------------------

       this_block = get_block(blocks_clinic(iblock),iblock)
       ib = this_block%ib
       ie = this_block%ie
       jb = this_block%jb
       je = this_block%je

       where (land_mask(:,:,iblock))
          tfact(:,:) = TAREA(:,:,iblock)
       elsewhere
          tfact(:,:) = 0.0_r8
       endwhere

       work1(:,:) = D14C(:,:,iblock) * tfact(:,:)
       d14c_local_sums(iblock)  = sum(work1(ib:ie,jb:je))
       tarea_local_sums(iblock) = sum(tfact(ib:ie,jb:je))

    end do

    !-----------------------------------------------------------------------
    ! Compute D14C making global mean
    !-----------------------------------------------------------------------

    d14c_sum_tmp  = sum(d14c_local_sums)
    tarea_sum_tmp = sum(tarea_local_sums)

    d14c_glo_avg  = global_sum(d14c_sum_tmp ,distrb_clinic) / global_sum(tarea_sum_tmp,distrb_clinic)

  end subroutine ciso_update_atm_D13C_D14C

  !***********************************************************************

  subroutine ciso_comp_varying_D13C(iblock, ciso_data_ind_d13c, D13C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of d13c when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate data values to current model time step
    !  2. Spatial patern of D13Cis the same everywhere (90 S - 90 N)
    !-----------------------------------------------------------------------

    use time_management , only : days_in_year
    use time_management , only : frac_day
    use time_management , only : iday_of_year
    use time_management , only : iyear
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt

    implicit none

    ! note that ciso_data_ind_d13c is always strictly less than the length
    ! of the data and is initialized to -1 before the first call

    integer (int_kind) , intent(in)  :: iblock                  ! block index
    integer (int_kind) , intent(out) :: ciso_data_ind_d13c      ! inex for the data for current timestep,
    real (r8)          , intent(out) :: D13C(nx_block,ny_block) ! atmospheric D13C (permil)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_comp_varying_D13C'

    integer (int_kind) :: &
         i, j              ! loop indices

    real (r8) :: &
         model_date,     & ! date of current model timestep
         mapped_date,    & ! model_date mapped to data timeline
         weight            ! weighting for temporal interpolation
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Generate mapped_date and check to see if it is too large.
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year

    if (mapped_date >= ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval)) then
       call document(subname, 'ciso: Model date maps to date after end of D13C data in file.')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    !--------------------------------------------------------------------------------------------------------------
    !  Set atmospheric D13C to first value in record for years before record begins
    !--------------------------------------------------------------------------------------------------------------

    if (mapped_date < ciso_atm_d13c_data_yr(1)) then
       D13C = ciso_atm_d13c_data(1)
       ciso_data_ind_d13c = 1
       if(my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Mapped date less than start of D13C data --> using first value in D13C data file'
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
       return
    endif

    !-----------------------------------------------------------------------
    !  On first time step, perform linear search to find data_ind_d13c
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d13c == -1) then
       do ciso_data_ind_d13c = ciso_atm_d13c_data_nbval-1,1,-1
          if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) exit
       end do
    endif

    !-----------------------------------------------------------------------
    !  See if ciso_data_ind_d13c needs to be updated,
    !  but do not set it to atm_d13c_data_nbval.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d13c < ciso_atm_d13c_data_nbval-1) then
       if (mapped_date >= ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1)) ciso_data_ind_d13c = ciso_data_ind_d13c + 1
    endif

    !-----------------------------------------------------------------------
    !  Generate hemisphere values for current time step.
    !-----------------------------------------------------------------------

    weight = (mapped_date - ciso_atm_d13c_data_yr(ciso_data_ind_d13c)) &
         / (ciso_atm_d13c_data_yr(ciso_data_ind_d13c+1) - ciso_atm_d13c_data_yr(ciso_data_ind_d13c))

    D13C = weight * ciso_atm_d13c_data(ciso_data_ind_d13c+1) + (c1-weight) * ciso_atm_d13c_data(ciso_data_ind_d13c)

  end subroutine ciso_comp_varying_D13C

  !***********************************************************************

  subroutine ciso_comp_varying_D14C(iblock, ciso_data_ind_d14c, D14C)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Compute atmospheric mole fractions of CO2 when temporarily
    !  varying data is read from files
    !  1. Linearly interpolate hemispheric values to current time step
    !  2. Make global field of D14C, determined by:
    !   -Northern Hemisphere value is used for 20N - 90 N
    !   -Southern Hemisphere value is used for 20 S - 90 S
    !   -Equator value is used for 20 S- 20 N

    use time_management, only : days_in_year
    use time_management, only : frac_day
    use time_management, only : iday_of_year
    use time_management, only : iyear
    use grid           , only : TLATD
    
    implicit none

    !  note that data_ind is always strictly less than the length of D14C data
    !  and is initialized to -1 before the first call

    integer (int_kind) , intent(in)  :: iblock                  ! block index
    integer (int_kind) , intent(out) :: ciso_data_ind_d14c      ! data_ind_d14c is the index into data for current timestep,
    real (r8)          , intent(out) :: D14C(nx_block,ny_block) ! atmospheric delta C14 in permil on global grid

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_comp_varying_D14C'
    integer (int_kind) :: &
         i, j, il        ! loop indices

    real (r8) :: &
         model_date,      & ! date of current model timestep
         mapped_date,     & ! model_date mapped to data timeline
         weight,          & ! weighting for temporal interpolation
         d14c_curr_sh,    & ! current atmospheric D14C value for SH (interpolated from data to model date)
         d14c_curr_nh,    & ! current atmospheric D14C value for NH (interpolated from data to model date)
         d14c_curr_eq       ! current atmospheric D14C value for EQ (interpolated from data to model date)
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Generate mapped_date and check to see if it is too large.
    !-----------------------------------------------------------------------

    model_date = iyear + (iday_of_year-1+frac_day)/days_in_year
    mapped_date = model_date - ciso_atm_model_year + ciso_atm_data_year
    do il=1,3
       if (mapped_date >= ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,il)) then
          call document(subname, 'ciso: model date maps to date after end of D14C data in files.')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
       endif
    enddo

    !--------------------------------------------------------------------------------------------------------------
    !  Set atmospheric D14C concentrations to zero before D14C record begins
    !--------------------------------------------------------------------------------------------------------------

    if (mapped_date < ciso_atm_d14c_data_yr(1,1)) then
       D14C = c0
       ciso_data_ind_d14c = 1
       if(my_task == master_task) then
          write(stdout,*)'ciso: Model date less than start of D14C data --> D14C=0'
       endif
       return
    endif

    !-----------------------------------------------------------------------
    !  On first time step, perform linear search to find data_ind_d14c.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d14c == -1) then
       do ciso_data_ind_d14c = ciso_atm_d14c_data_nbval-1,1,-1
          if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) exit
       end do
    endif

    !-----------------------------------------------------------------------
    !  See if data_ind_d14c need to be updated,
    !  but do not set it to atm_co2_data_nbval.
    !-----------------------------------------------------------------------

    if (ciso_data_ind_d14c < ciso_atm_d14c_data_nbval-1) then
       if (mapped_date >= ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1))  then
          ciso_data_ind_d14c = ciso_data_ind_d14c + 1
       endif
    endif
    !
    !-----------------------------------------------------------------------
    !  Generate hemisphere values for current time step.
    !-----------------------------------------------------------------------

    weight = (mapped_date - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1)) &
         / (ciso_atm_d14c_data_yr(ciso_data_ind_d14c+1,1) - ciso_atm_d14c_data_yr(ciso_data_ind_d14c,1))

    d14c_curr_sh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,1) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,1)
    d14c_curr_eq = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,2) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,2)
    d14c_curr_nh = weight * ciso_atm_d14c_data(ciso_data_ind_d14c+1,3) + &
              (c1-weight) * ciso_atm_d14c_data(ciso_data_ind_d14c,3)

    !-----------------------------------------------------------------------
    !  Merge hemisphere values for D14C
    !      -Northern Hemisphere value is used for >20N - 90 N
    !      -Southern Hemisphere value is used for >20 S - 90 S
    !      -Equatorial value is used for 20 S to 20 N
    !-----------------------------------------------------------------------

    do j = 1, ny_block
       do i = 1, nx_block
          if (TLATD(i,j,iblock) < -20.0_r8) then
             D14C(i,j) = d14c_curr_sh
          else if (TLATD(i,j,iblock) > 20.0_r8) then
             D14C(i,j) = d14c_curr_nh
          else
             D14C(i,j) = d14c_curr_eq
          endif
       end do
    end do

  end subroutine ciso_comp_varying_D14C

  !*****************************************************************************

  subroutine ciso_init_atm_D13_D14

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Initialize surface flux computations for the ecosys_ciso tracer module.
    !  Includes reading CO2 and D13C and D14C data from file if option file is used
    !---------------------------------------------------------------------

    use io_types        , only : stdout
    use constants       , only : blank_fmt      
    use constants       , only : delim_fmt      
    use constants       , only : ndelim_fmt     

    implicit none

    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_init_atm_D13_D14'

    !-------------------------------------------------------------------------
    !     Set D13C data source
    !-------------------------------------------------------------------------

    select case (ciso_atm_d13c_opt)
    case ('const')
       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Using constant D13C values of ',ciso_atm_d13c_const
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
    case('file')
       call ciso_read_atm_D13C_data ! READ in D13C data from file
    case default
       call document(subname, 'unknown ciso_atm_d13c_opt in ecosys_ciso_init_atm_d13_d14')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    end select

    !-------------------------------------------------------------------------
    !     Set D14C data source
    !-------------------------------------------------------------------------

    select case (ciso_atm_d14c_opt)
    case ('const')
       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
          write(stdout,*)'ciso: Using constant D14C values of ',ciso_atm_d14c_const
          write(stdout,blank_fmt)
          write(stdout,ndelim_fmt)
          write(stdout,blank_fmt)
       endif
    case('file')
       call ciso_read_atm_D14C_data ! READ in D14C data from files
    case default
       call document(subname, 'unknown ciso_atm_d14c_opt in ecosys_ciso_init_atm_d13_d14')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    end select

  end subroutine ciso_init_atm_D13_D14

  !*****************************************************************************

  subroutine ciso_read_atm_D13C_data()

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Read atmospheric D13C [permil] data from file
    !
    !  Have the master_task do the following :
    !     1) get length of data
    !     2) allocate memory for data
    !     3) read in data, checking for consistent lengths
    !  Then, outside master_task conditional
    !     1) broadcast length of data
    !     2) have non-mastertasks allocate memory for data
    !     3) broadcast data
    !-----------------------------------------------------------------------

    use broadcast       , only : broadcast_array
    use broadcast       , only : broadcast_scalar
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt
    use io_types        , only : nml_in

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_forcing_mod:ciso_read_atm_D13C_data'
    integer (int_kind) ::    &
         stat,                 &  ! i/o status code
         irec,                 &  ! counter for looping
         skiplines,            &  ! number of comment lines at beginning of ascii file
         il                       ! looping index
    character (char_len) :: &
         sglchr                   ! variable to read characters from file into
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !     READ in D13C data from file
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*)'ciso: Using varying D13C values from file ',trim(ciso_atm_d13c_filename)
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       open (nml_in, file=ciso_atm_d13c_filename, status='old',iostat=stat)
       if (stat /= 0) then
          write(stdout,fmt=*) 'open failed'
          go to 99
       endif
       read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d13c_data_nbval
       if (stat /= 0) then
          write(stdout,fmt=*) '1st line read failed'
          go to 99
       endif
       allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
       allocate(ciso_atm_d13c_data   (ciso_atm_d13c_data_nbval))
       do irec=1,skiplines
          read(nml_in,FMT=*,iostat=stat) sglchr
          if (stat /= 0) then
             write(stdout,fmt=*) 'skipline read failed'
             go to 99
          endif
       enddo
       do irec=1,ciso_atm_d13c_data_nbval
          read(nml_in,FMT=*,iostat=stat) ciso_atm_d13c_data_yr(irec), ciso_atm_d13c_data(irec)
          if (stat /= 0) then
             write(stdout,fmt=*) 'data read failed'
             go to 99
          endif
       enddo
       close(nml_in)
    endif

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'Stopping in ' // subname)

    !---------------------------------------------------------------------
    !     Need to allocate and broadcast the variables to other tasks beside master-task
    !---------------------------------------------------------------------

    call broadcast_scalar(ciso_atm_d13c_data_nbval,master_task)

    if (my_task /= master_task) then
       allocate(ciso_atm_d13c_data_yr(ciso_atm_d13c_data_nbval))
       allocate(ciso_atm_d13c_data(ciso_atm_d13c_data_nbval))
    endif

    call broadcast_array(ciso_atm_d13c_data   , master_task)
    call broadcast_array(ciso_atm_d13c_data_yr, master_task)

  end subroutine ciso_read_atm_D13C_data

  !*****************************************************************************

  subroutine ciso_read_atm_D14C_data

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Read atmospheric D14C data from file
    !
    !  Have the master_task do the following :
    !     1) get length of data
    !     2) allocate memory for data
    !     3) read in data, checking for consistent lengths
    !  Then, outside master_task conditional
    !     1) broadcast length of data
    !     2) have non-mastertasks allocate memory for data
    !     3) broadcast data
    !-----------------------------------------------------------------------

    use broadcast       , only : broadcast_array
    use broadcast       , only : broadcast_scalar
    use constants       , only : blank_fmt
    use constants       , only : delim_fmt
    use constants       , only : ndelim_fmt
    use io_types        , only : nml_in

    implicit none

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_forcing_mod:ciso_read_atm_D14C_data'

    integer (int_kind) ::      &
         stat,                   &  ! i/o status code
         irec,                   &  ! counter for looping
         skiplines,              &  ! number of comment lines at beginning of ascii file
         il                         ! looping index

    character (char_len) ::  &
         sglchr                     ! variable to read characters from file into

    integer (int_kind) :: &
         ciso_atm_d14c_data_nbval_tmp

    logical (log_kind) :: &
         nbval_mismatch

    !-----------------------------------------------------------------------
    !     ensure that three datafiles have same number of entries
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       write(stdout,*)'ciso DIC14 calculation: Using varying C14 values from files'
       do il=1,3
          write(stdout,*) trim(ciso_atm_d14c_filename(il))
       enddo
       nbval_mismatch = .false.
       do il=1,3
          open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
          if (stat /= 0) then
             write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
          if (stat /= 0) then
             write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          close(nml_in)
          if (il == 1) then
             ciso_atm_d14c_data_nbval = ciso_atm_d14c_data_nbval_tmp
          else
             if (ciso_atm_d14c_data_nbval /= ciso_atm_d14c_data_nbval_tmp) nbval_mismatch = .true.
          endif
       enddo
    endif

    call broadcast_scalar(nbval_mismatch, master_task)
    if (nbval_mismatch) then
       call document(subname, 'D14C data files must all have the same number of values')
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    call broadcast_scalar(ciso_atm_d14c_data_nbval, master_task)
    allocate(ciso_atm_d14c_data_yr(ciso_atm_d14c_data_nbval,3))
    allocate(ciso_atm_d14c_data   (ciso_atm_d14c_data_nbval,3))

    !-----------------------------------------------------------------------
    !     READ in C14 data from files - three files, for SH, EQ, NH
    !-----------------------------------------------------------------------

    if (my_task == master_task) then
       do il=1,3
          open (nml_in,file=ciso_atm_d14c_filename(il),status='old',iostat=stat)
          if (stat /= 0) then
             write(stdout,*) 'open failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          read(nml_in,FMT=*,iostat=stat) skiplines,ciso_atm_d14c_data_nbval_tmp
          if (stat /= 0) then
             write(stdout,*) '1st line read failed for ', trim(ciso_atm_d14c_filename(il))
             go to 99
          endif
          do irec=1,skiplines
             read(nml_in,FMT=*,iostat=stat) sglchr
             if (stat /= 0) then
                write(stdout,fmt=*) 'skipline read failed for ', trim(ciso_atm_d14c_filename(il))
                go to 99
             endif
          enddo
          do irec=1,ciso_atm_d14c_data_nbval
             read(nml_in,FMT=*,iostat=stat) ciso_atm_d14c_data_yr(irec,il), ciso_atm_d14c_data(irec,il)
             if (stat /= 0) then
                write(stdout,fmt=*) 'data read failed for ', trim(ciso_atm_d14c_filename(il))
                go to 99
             endif
          enddo
          close(nml_in)
       enddo
    endif

99  call broadcast_scalar(stat, master_task)
    if (stat /= 0) call exit_POP(sigAbort, 'Stopping in ' // subname)

    !---------------------------------------------------------------------
    ! Broadcast the variables to other tasks beside master_task
    !---------------------------------------------------------------------

    call broadcast_array(ciso_atm_d14c_data   , master_task)
    call broadcast_array(ciso_atm_d14c_data_yr, master_task)

  end subroutine ciso_read_atm_D14C_data

  !*****************************************************************************

  subroutine forcing_constant_init(this, field_constant)

    class(forcing_constant_type), intent(inout) :: this
    real(kind=r8),                intent(in)    :: field_constant

    this%field_constant = field_constant

  end subroutine forcing_constant_init

  !*****************************************************************************

  subroutine forcing_driver_init(this, driver_varname)

    class(forcing_driver_type), intent(inout) :: this
    character(len=*),           intent(in)    :: driver_varname

    this%driver_varname = driver_varname

  end subroutine forcing_driver_init

  !*****************************************************************************

  subroutine forcing_named_field_init(this, field_name)

    use named_field_mod       , only : named_field_get_index

    class(forcing_named_field_type), intent(inout) :: this
    character(len=*),                intent(in)    :: field_name

    this%field_name = field_name
    call named_field_get_index(field_name, this%field_ind) 

  end subroutine forcing_named_field_init

  !*****************************************************************************

  subroutine forcing_file_init(this, filename, file_varname, temporal, year_first, &
                                     year_last, year_align, date, time)


    class(forcing_file_type),         intent(inout) :: this
    character(len=*),                 intent(in)    :: filename
    character(len=*),                 intent(in)    :: file_varname
    character(len=*),       optional, intent(in)    :: temporal
    integer(kind=int_kind), optional, intent(in)    :: year_first
    integer(kind=int_kind), optional, intent(in)    :: year_last
    integer(kind=int_kind), optional, intent(in)    :: year_align
    integer(kind=int_kind), optional, intent(in)    :: date
    integer(kind=int_kind), optional, intent(in)    :: time

    this%filename     = filename
    this%file_varname = file_varname
    if (present(temporal  )) this%temporal   = temporal
    if (present(year_first)) this%year_first = year_first
    if (present(year_last )) this%year_last  = year_last
    if (present(year_align)) this%year_align = year_align
    if (present(date      )) this%date       = date
    if (present(time      )) this%time       = time

  end subroutine forcing_file_init

  !*****************************************************************************

  subroutine forcing_monthly_calendar_init(this, forcing_calendar_name)

    class(forcing_monthly_calendar_type)   , intent(inout) :: this
    type (forcing_monthly_every_ts), target, intent(in)    :: forcing_calendar_name

    this%forcing_calendar_name => forcing_calendar_name

  end subroutine forcing_monthly_calendar_init

  !*****************************************************************************

  subroutine forcing_field_metadata_set(this, &
       num_elements,                          &
       field_source,                          &
       marbl_varname,                         &
       field_units,                           &
       unit_conv_factor, temporal_interp,     &
       field_constant,                        &
       driver_varname,                        &
       named_field,                           &
       filename,                              &
       file_varname,                          &
       temporal,                              &
       year_first, year_last, year_align,     &
       date,                                  &
       time,                                  &
       forcing_calendar_name)

    use io_tools        , only : document

    class(forcing_fields_metadata_type), intent(inout) :: this
    integer (kind=int_kind),             intent(in)    :: num_elements
    character (len=*),                   intent(in)    :: field_source  ! must  have valid_field_source value)
    character (len=*),                   intent(in)    :: marbl_varname ! required
    character (len=*),                   intent(in)    :: field_units
    real(kind=r8),           optional,   intent(in)    :: unit_conv_factor
    character (len=*),       optional,   intent(in)    :: temporal_interp
    real(kind=r8),           optional,   intent(in)    :: field_constant
    character (len=*),       optional,   intent(in)    :: driver_varname
    character (len=*),       optional,   intent(in)    :: named_field
    character (len=*),       optional,   intent(in)    :: filename
    character (len=*),       optional,   intent(in)    :: file_varname
    character (len=*),       optional,   intent(in)    :: temporal
    integer (kind=int_kind), optional,   intent(in)    :: year_first
    integer (kind=int_kind), optional,   intent(in)    :: year_last
    integer (kind=int_kind), optional,   intent(in)    :: year_align
    integer (kind=int_kind), optional,   intent(in)    :: date
    integer (kind=int_kind), optional,   intent(in)    :: time
    type(forcing_monthly_every_ts), optional, target, intent(in) :: forcing_calendar_name

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=char_len), dimension(7) :: valid_field_sources
    integer(kind=int_kind)  :: n
    logical(log_kind)       :: has_valid_source
    logical(log_kind)       :: has_valid_inputs
    character(*), parameter :: subname = 'ecosys_forcing_mod:forcing_field_metadata_set'
    character(len=char_len) :: log_message
    !-----------------------------------------------------------------------

    valid_field_sources(1) = 'const'
    valid_field_sources(2) = 'zero'
    valid_field_sources(3) = 'internal'
    valid_field_sources(4) = 'named_field'
    valid_field_sources(5) = 'shr_stream'
    valid_field_sources(6) = 'file_time_invariant'
    valid_field_sources(7) = 'POP monthly calendar'

    ! check for valid source
    has_valid_source = .false.
    do n = 1,size(valid_field_sources)
       if (trim(field_source) .eq. trim(valid_field_sources(n))) has_valid_source = .true.
    enddo
    if (.not. has_valid_source) then
       write(log_message,"(4A)") trim(field_source),                          &
                                 " is not a valid source for reading the ",   &
                                 trim(marbl_varname), " forcing field"

       call document(subname, log_message)
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    endif

    ! required variables for all forcing field sources
    this%field_source  = trim(field_source)
    this%marbl_varname = marbl_varname
    this%field_units   = field_units

    ! optional variables
    this%unit_conv_factor = c1
    if (present(unit_conv_factor)) this%unit_conv_factor = unit_conv_factor

    this%temporal_interp  = ''
    if (present(temporal_interp )) this%temporal_interp  = temporal_interp


    ! optional variables for forcing field type

    ! each forcing type has its own requirements - if we check here, then the
    ! separate type inits can have fewer optional arguments

    has_valid_inputs = .true.

    select case (trim(field_source))

    case('const')
       if (.not.present(field_constant)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          write(log_message,"(2A)") "Adding constant forcing_field_type for ", &
                                    trim(this%marbl_varname)
          call document(subname, log_message)
          call forcing_constant_init(this%field_constant_info, field_constant)
       endif

    case('zero')
       write(log_message,"(2A)") "Adding constant (0) forcing_field_type for ", &
                                 trim(this%marbl_varname)
       call document(subname, log_message)
       call forcing_constant_init(this%field_constant_info, c0)

    case('internal')
       if (.not.present(driver_varname)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          write(log_message, "(2A)") "Adding internal forcing_field_type for ",  &
                                    trim(this%marbl_varname)
          call document(subname, log_message)
          call this%field_driver_info%initialize(driver_varname)
       endif

    case('named_field')
       if (.not.present(named_field)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          write(log_message, "(2A)") "Adding named field forcing_field_type for ",  &
                                    trim(this%marbl_varname)
          call document(subname, log_message)
          call this%field_named_info%initialize(named_field)
       endif

    case('shr_stream','file_time_invariant') 
       if (.not.present(filename))     has_valid_inputs = .false.
       if (.not.present(file_varname)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          write(log_message,"(2A)") "Adding file forcing_field_type for ",     &
                                   trim(this%marbl_varname)
          call document(subname, log_message)
          call this%field_file_info%initialize(&
               filename, file_varname, &
               temporal=temporal, year_first=year_first,   &
               year_last=year_last, year_align=year_align, &
               date=date, time=time)
       endif

    case('POP monthly calendar') 
       if (.not.present(forcing_calendar_name)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          write(log_message,"(2A)") "Adding calendar forcing_field_type for ", &
                                   trim(this%marbl_varname)
          call document(subname, log_message)
          call this%field_monthly_calendar_info%initialize(forcing_calendar_name)
       endif

    end select

    if (.not.has_valid_inputs) then
      write(log_message,"(3A)") "Call to forcing_field%init does not have ",  &
                                "the correct optional arguments for ",        &
                                trim(field_source)
       call document(subname, log_message)
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

   end subroutine forcing_field_metadata_set

  !*****************************************************************************

  subroutine forcing_fields_add(this,       &
       field_source,                        &
       marbl_varname,                       &
       field_units,                         &
       id,                                  &
       unit_conv_factor,                    &
       temporal_interp,                     &
       field_constant,                      &
       driver_varname,                      &
       named_field,                         &
       filename, file_varname,              &
       temporal,                            &
       year_first, year_last, year_align,   &
       date, time,                          &
       forcing_calendar_name)

    class(forcing_fields_type)      , intent(inout) :: this
    character(len=*)                , intent(in)    :: field_source
    character(len=*)                , intent(in)    :: marbl_varname
    character(len=*)                , intent(in)    :: field_units
    integer(kind=int_kind)          , intent(in)    :: id
    real(kind=r8),          optional, intent(in)    :: unit_conv_factor
    character(len=*),       optional, intent(in)    :: temporal_interp
    real(kind=r8),          optional, intent(in)    :: field_constant
    character(len=*),       optional, intent(in)    :: driver_varname
    character(len=*),       optional, intent(in)    :: named_field
    character(len=*),       optional, intent(in)    :: filename
    character(len=*),       optional, intent(in)    :: file_varname
    character(len=*),       optional, intent(in)    :: temporal
    integer(kind=int_kind), optional, intent(in)    :: year_first
    integer(kind=int_kind), optional, intent(in)    :: year_last
    integer(kind=int_kind), optional, intent(in)    :: year_align
    integer(kind=int_kind), optional, intent(in)    :: date
    integer(kind=int_kind), optional, intent(in)    :: time
    type(forcing_monthly_every_ts), optional, target, intent(in) :: forcing_calendar_name

    integer (kind=int_kind) :: num_elem
    character(*), parameter :: subname = 'ecosys_forcing_mod:forcing_fields_add'
    character(len=char_len) :: log_message

    call this%metadata%set(                              &
         num_elem, field_source, marbl_varname,          &
         field_units, unit_conv_factor=unit_conv_factor, &
         temporal_interp=temporal_interp,                &
         field_constant=field_constant,                  &
         driver_varname=driver_varname,                  &
         named_field=named_field,                        &
         filename=filename, file_varname=file_varname,   &
         temporal=temporal, year_first=year_first,       &
         year_last=year_last, year_align=year_align,     &
         date=date, time=time,                           &
         forcing_calendar_name=forcing_calendar_name)

  end subroutine forcing_fields_add

  !*****************************************************************************

  subroutine init_monthly_surface_forcing_metadata(var)

    implicit none

    type(forcing_monthly_every_ts), intent(inout) :: var

    var%interp_type = 'linear'
    var%data_type   = 'monthly-calendar'
    var%interp_freq = 'every-timestep'
    var%filename    = 'not-used-for-monthly'
    var%data_label  = 'not-used-for-monthly'

  end subroutine init_monthly_surface_forcing_metadata

  !*****************************************************************************

  function ecosys_forcing_tracer_ref_val(ind)
    !
    ! !DESCRIPTION:
    !  return reference value for tracer with global tracer index ind
    !  this is used in virtual flux computations

    implicit none

    integer(int_kind) , intent(in) :: ind
    real(r8) :: ecosys_forcing_tracer_ref_val

    !  default value for reference value is 0

    ecosys_forcing_tracer_ref_val = c0
    if (vflux_flag(ind)) then
       ecosys_forcing_tracer_ref_val = surf_avg(ind)
    endif
       
  end function ecosys_forcing_tracer_ref_val

  !*****************************************************************************

end module ecosys_forcing_mod
